#include "opt_levelset.hpp"

using namespace std;
using namespace icarat;
using namespace icarat::opt_levelset;
using namespace Eigen;

void makeEleID(FEMLS &fem, vector<ElastLevelset> &element,
               vector<NodeLS> &node);

void topologicalDerivative(FEMLS &fem, vector<ElastLevelset> &element,
                           vector<NodeLS> &node, int optstep);

void updateLevelsetFunction(FEMLS &fem, vector<ElastLevelset> &element,
                            vector<NodeLS> &node, double dt);

void makeKMReact(int nel, vector<Triplet<double>> &tripletsK,
                 vector<Triplet<double>> &tripletsF, FEMLS &fem,
                 vector<ElastLevelset> &element, vector<NodeLS> &node,
                 double dt);

void output(FEMLS &fem, vector<ElastLevelset> &element, vector<NodeLS> &node,
            int optstep);

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <.toml file>\n";
    return 1;
  }

  auto t_start = chrono::system_clock::now(); // measure time
  int convergence = 1;                        // convergence judge
  int optstep = 0;

  // set input file
  const auto config = toml::parse(argv[1]);

  // set problem
  double V0 = 0.0; // initial volume

  FEMLS fem;
  fem.tau = toml::find<double>(config, "this_problem", "tau");
  // time increment
  double dt = toml::find<double>(config, "this_problem", "dt");
  // volume fraction
  double Vfrac = toml::find<double>(config, "this_problem", "Vfrac");
  // Rate of volume reduction
  double nvol = toml::find<double>(config, "this_problem", "nVol");
  // heaviside parameter
  double dd = toml::find<double>(config, "this_problem", "dd");
  double pp = toml::find<double>(config, "this_problem", "pp");

  MeshIcarat<ElastLevelset, NodeLS> mesh(fem, config);
  mesh.setParameter("structure");

  // set element
  vector<ElastLevelset> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      ElastLevelset actele(fem.voigt, mesh.ne(type), mesh.ipmax(type), numdof);
      actele.setParameter(config);
      element.push_back(actele);
    }
  }

  // set node
  vector<NodeLS> node(fem.numnp);

  // input mesh to element & node
  mesh.generate(element, node);

  // set element ID to each node
  makeEleID(fem, element, node);

  // boundary condition class
  BCIcarat<ElastLevelset, NodeLS> BC(fem, element, node);

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  Value dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  // initial levelset function field
  for(int i = 0; i < fem.numnp; i++)
    node[i].phi = 1.0;

  // analysis class
  LinearAnalysis<ElastLevelset, NodeLS> compt(fem, element, node, force, dis);

  // optimization start
  while(convergence == 1)
  {
    optstep++;

    std::cout << "============== optimization step ==============" << endl;
    std::cout << "step no. " << optstep << endl;

    compt.solve("Eigen_CG", 1);

    /// get  stress
    for(auto &e : element)
      e.makeMisesStress(fem, node);

    if(optstep == 1)
    {
      V0 = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        V0 += element[i].design_s() * element[i].volume;
      if(V0 == 0.0)
      {
        cerr << "V0 = 0. in optimization" << endl;
        exit(1);
      }
    }

    // get objective function
    Checker check;
    double object_f = dis.val.dot(force.fext);
    std::cout << "function value = " << object_f << endl;

    // output result
    addPlot(optstep, object_f, "function.csv", optstep);
    output(fem, element, node, optstep);

    convergence = check.isConvergence(object_f, optstep, 1.0e-4, 100);

    topologicalDerivative(fem, element, node, optstep);

    // compute volume & adjoint variable
    double V = 0.0;
    for(int i = 0; i < fem.nelx; i++)
      V += element[i].design_s() * element[i].volume;
    std::cout << "volume = " << V << endl;

    double ex = Vfrac * V0 +
                (V0 - Vfrac * V0) * std::max(0.0, 1.0 - (optstep + 1) / nvol);

    fem.lambda = fem.sum_tdn / fem.numnp * exp(pp * ((V - ex) / ex + dd));

    // update levelset function & characteristic function
    updateLevelsetFunction(fem, element, node, dt);

    // not convergence
    if(optstep == 1000)
    {
      cerr << "not convergence in this optimization" << endl;
      exit(1);
    }
  }

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  std::cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void makeEleID(FEMLS &fem, vector<ElastLevelset> &element, vector<NodeLS> &node)
{
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    // add element ID to each node
    for(int i = 0; i < element[nel].ne; i++)
      node[element[nel].nodeID[i]].eleID.push_back(element[nel].ID);
  }
}

void topologicalDerivative(FEMLS &fem, vector<ElastLevelset> &element,
                           vector<NodeLS> &node, int optstep)
{
  double kai_min = 1.0e-4;
  vector<double> TD(fem.nelx); ///< topological derivative in each element

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    double A1 =
        -1.5 * (1.0 - element[nel].poisson()) *
        (1.0 - 14.0 * element[nel].poisson() +
         15.0 * pow(element[nel].poisson(), 2.0)) *
        element[nel].young1() /
        ((1.0 + element[nel].poisson()) * (7.0 - 5.0 * element[nel].poisson()) *
         pow(1.0 - 2.0 * element[nel].poisson(), 2.0));

    double A2 =
        7.5 * (1.0 - element[nel].poisson()) * element[nel].young1() /
        ((1.0 + element[nel].poisson()) * (7.0 - 5.0 * element[nel].poisson()));

    // poisson ratio
    double cc = A1 / (A1 + 2.0 * A2);
    // young modulus
    double bb = (A1 + 2.0 * A2) * (1.0 - pow(cc, 2.0));

    // std::cout << "cc " << cc << "bb " << bb << endl;

    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * element[nel].ne);
    MatrixXd Be = MatrixXd::Zero(fem.voigt, element[nel].numdof);
    MatrixXd Ke = MatrixXd::Zero(element[nel].numdof, element[nel].numdof);
    MatrixXd X = MatrixXd::Zero(element[nel].ne, fem.ndim);
    VectorXd disp = VectorXd::Zero(element[nel].numdof);

    for(int i = 0; i < element[nel].ne; i++)
    {
      for(int j = 0; j < fem.ndim; j++)
        X.coeffRef(i, j) = node[element[nel].nodeID[i]].x[j];

      for(int j = 0; j < fem.dofnp; j++)
        disp.coeffRef(fem.dofnp * i + j) = node[element[nel].nodeID[i]].val[j];
    }

    MatrixXd De(fem.voigt, fem.voigt);
    element[nel].De(De, bb, cc);

    Bmatrix bmatrix;
    double jac = 0.0;
    for(int ip = 0; ip < element[nel].ipmax; ip++)
    {
      bmatrix.make(Be, Ne, jac, element[nel].eType, element[nel].ipmax, X, ip);
      Ke += Be.transpose() * De * Be * jac;
    }

    TD[nel] = (kai_min + element[nel].design_s() * (1.0 - kai_min)) *
              disp.transpose() * Ke * disp;
  }

  // element value to node value
  for(int i = 0; i < fem.numnp; i++)
  {
    node[i].tdn = 0.0;
    for(int j = 0; j < node[i].eleID.size(); j++)
      node[i].tdn += TD[node[i].eleID[j]] / (int)node[i].eleID.size();
  }

  // update parameter C
  double CC = 0.0;
  fem.sum_tdn = 0.0;
  for(int i = 0; i < fem.numnp; i++)
  {
    CC += fabs(node[i].tdn);
    fem.sum_tdn += node[i].tdn;
  }

  fem.C = fem.nelx / CC;

  // output sensivity graph
  if(optstep == 1)
    exportGraph(TD, "objective.csv");
}

void updateLevelsetFunction(FEMLS &fem, vector<ElastLevelset> &element,
                            vector<NodeLS> &node, double dt)
{
  // make element matrix
  vector<Triplet<double>> tripletsK, tripletsF;
  for(int nel = 0; nel < fem.nelx; nel++)
    makeKMReact(nel, tripletsK, tripletsF, fem, element, node, dt);

  // assembling
  SparseMatrix<double> globalK(fem.numnp, fem.numnp),
      fMat(fem.numnp, fem.numnp);

  globalK.setFromTriplets(tripletsK.begin(), tripletsK.end());
  globalK.makeCompressed();

  fMat.setFromTriplets(tripletsF.begin(), tripletsF.end());
  VectorXd F = fMat.diagonal();

  /// solve KU=F for design_s filter
  VectorXd val = VectorXd::Zero(fem.numnp);
  Solvers solv(globalK, F, val);
  solv.solve("Eigen_CG", 0);

  /// set phi for output
  for(int i = 0; i < fem.numnp; i++)
    node[i].phi = max(min(1.0, val.coeff(i)), -1.0);

  // element value to node value & update kai
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    double phie = 0.0;
    for(int j = 0; j < element[nel].ne; j++)
      phie += node[element[nel].nodeID[j]].phi / element[nel].ne;

    if(phie < 0.0)
    {
      double zero = 0.0;
      element[nel].setDesignS(zero);
    }
    else
    {
      double one = 1.0;
      element[nel].setDesignS(one);
    }
  }
}

void makeKMReact(int nel, vector<Triplet<double>> &tripletsK,
                 vector<Triplet<double>> &tripletsF, FEMLS &fem,
                 vector<ElastLevelset> &element, vector<NodeLS> &node,
                 double dt)
{
  MatrixXd Ne =
      MatrixXd::Zero(fem.ndim, fem.ndim * element[nel].ne); ///< shape function
  VectorXd NeArr = VectorXd::Zero(element[nel].ne); ///< shape function array
  MatrixXd Be = MatrixXd::Zero(fem.ndim, element[nel].ne); ///< B matrix
  VectorXd fe = VectorXd::Zero(element[nel].ne);           ///< inner force
  MatrixXd X = MatrixXd::Zero(element[nel].ne, fem.ndim);  ///< coordinate
  VectorXi idof = VectorXi::Zero(element[nel].ne); ///< dof in this element
  MatrixXd Ke = MatrixXd::Zero(element[nel].ne, element[nel].ne);
  MatrixXd Me = MatrixXd::Zero(element[nel].ne, element[nel].ne);
  VectorXd phie = VectorXd::Zero(element[nel].ne); ///< levelset function values
  VectorXd tde = VectorXd::Zero(element[nel].ne);  ///< topological derivative

  for(int i = 0; i < element[nel].ne; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
      X.coeffRef(i, j) = node[element[nel].nodeID[i]].x[j];

    idof.coeffRef(i) = element[nel].nodeID[i];
    phie.coeffRef(i) = node[element[nel].nodeID[i]].phi;
    tde.coeffRef(i) = node[element[nel].nodeID[i]].tdn;
  }

  double jac = 0.0;
  Bmatrix bmatrix;

  // gauss point loop start
  for(int ip = 0; ip < element[nel].ipmax; ip++)
  {
    bmatrix.make(Be, Ne, jac, element[nel].eType, element[nel].ipmax, X, ip);

    for(int i = 0; i < element[nel].ne; i++)
      NeArr.coeffRef(i) = Ne.coeff(0, fem.ndim * i);

    double uu = NeArr.dot(tde);

    fe += fem.C * (uu - fem.lambda) * NeArr * jac;
    Ke += fem.tau * fem.nelx * Be.transpose() * Be * jac;
    Me += NeArr * NeArr.transpose() * jac;
  } // gauss point loop end

  MatrixXd Te = Me / dt + Ke;
  VectorXd Ye = Me / dt * phie + fe;

  assembling(element[nel].ne, fem.numeq, idof, Te, Ye, tripletsK, tripletsF);
}

void output(FEMLS &fem, vector<ElastLevelset> &element, vector<NodeLS> &node,
            int optstep)
{
  ofstream fout("res" + to_string(optstep) + ".vtk");
  makeHeaderVTK(fout);

  // make grid
  setPointVTK(node, fout);
  setElementVTK(element, fout);
  setElementTypeVTK(fem.ndim, element, fout);

  // displacement
  addPointVectorVTK("displacement", node, fout, true,
                    [&](NodeLS &n) { return n.val; });

  // levelset function
  addPointScalarVTK("phi", node, fout, false, [&](NodeLS &n) { return n.phi; });

  // charastic function "kai"
  addElementScalarVTK("design_s", element, fout, true,
                      [&](ElastLevelset &e) { return e.design_s(); });

  // mises stress
  addElementScalarVTK("mises", element, fout, false,
                      [&](ElastLevelset &e) { return e.mises(); });

  fout.close();
}