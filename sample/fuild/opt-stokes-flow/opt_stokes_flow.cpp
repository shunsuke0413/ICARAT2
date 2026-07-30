#include "opt_stokes_flow.hpp"
using namespace icarat;
using namespace std;
using namespace Eigen;

void sensitivity(FEM &fem, vector<Stokes> &element, vector<Node> &node,
                 Optimize &optim, double V0);
void output(FEM &fem, vector<Stokes> &element, vector<Node> &node,
            Optimize &optim);

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <.toml file>\n";
    return 1;
  }

  // measure time
  auto t_start = chrono::system_clock::now();
  int convergence = 1; // convergence judge
  int numConst = 1;    // number of constraint
  double V0 = 0.0;     // first volume

  // set input file
  const auto config = toml::parse(argv[1]);
  double filterR = toml::find<double>(config, "this_problem", "filterR");

  // set problem
  FEM fem;
  MeshIcarat<Stokes, Node> mesh(fem, config);
  mesh.setParameter("fluid");

  vector<Stokes> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      Stokes actele(fem.voigt, mesh.ne(type), mesh.ipmax(type), numdof);
      actele.setParameter(config);
      element.push_back(actele);
    }
  }

  // set node
  vector<Node> node(fem.numnp);

  // input mesh to element & node
  mesh.generate(element, node);

  // boundary condition class
  BCIcarat<Stokes, Node> BC(fem, element, node);

  BC.input_dirich(config);

  // get number of effective dofs
  fem.numeq = BC.makeDOFFluid();

  // set force & displacement
  Force force(fem.numeq);
  Value dis(fem.numeq);
  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  // set optimization's classes
  Optimize optim(fem.nelx, numConst, 1.0e-4, 1.0);
  MMASolver mma(fem.nelx, numConst);
  FilterHelmholtz<Stokes, Node> filter(fem, element, node, filterR);

  // init
  optim.optstep = 0;
  for(int i = 0; i < fem.nelx; i++)
    optim.design_s[i] = element[i].design_s();

  LinearAnalysis<Stokes, Node> compt(fem, element, node, force, dis);

  output(fem, element, node, optim);

  // optimization start
  while(convergence == 1)
  {
    optim.optstep++;

    cout << "============== optimization step ==============" << endl;
    cout << "step no. " << optim.optstep << endl;

    // get filtered design_s
    auto design_filtered =
        filter.densityFilter(optim.design_s, "Eigen_CG", 1, false);

    /// set threshold function
    double beta = min(pow(2, (int)(optim.optstep / 40)), 16.0);
    auto design_threshold = filter.threshold(design_filtered, beta, 0.5);

    // set threshold design variable
    for(int i = 0; i < fem.nelx; i++)
      element[i].setDesign(design_threshold[i]);

    // structural analysis
    compt.solve("Eigen_CG", 1);

    // get objective function
    optim.object_f = dis.val.transpose() * compt.K() * dis.val;
    cout << "function value = " << optim.object_f << endl;

    // get const_func
    if(optim.optstep == 1)
    {
      V0 = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        V0 += element[i].design_s() * element[i].volume;

      optim.const_h[0] = 0.0;

      if(V0 == 0.0)
        cerr << "V0 = 0. in optimization" << endl;
    }
    else
    {
      optim.const_h[0] = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        optim.const_h[0] += element[i].design_s() * element[i].volume / V0;
      optim.const_h[0] -= 1.0;
    }
    cout << "Constraint value:" << optim.const_h[0] << endl;

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);

    output(fem, element, node, optim);

    exit(1);

    convergence = optim.isConvergence(1.0e-7, 50);

    // analytical
    vector<double> threDifferent =
        filter.thresholdDerivative(design_filtered, beta, 0.5);

    sensitivity(fem, element, node, optim, V0);

    for(int i = 0; i < fem.nelx; i++)
    {
      optim.dfds[i] *= threDifferent[i];
      optim.dgds[i] *= threDifferent[i];
    }

    optim.dfds = filter.densityFilter(optim.dfds, "Eigen_CG", 0, true);
    optim.dgds = filter.densityFilter(optim.dgds, "Eigen_CG", 0, true);

    // output sensivity graph
    if(optim.optstep == 1)
      exportGraph(optim.dfds, "objective.csv");

    // Get updated design variables with MMA
    mma.Update(optim.design_s.data(), optim.dfds.data(), optim.const_h.data(),
               optim.dgds.data(), optim.xmin.data(), optim.xmax.data());

    // not convergence
    if(optim.optstep == 1000)
    {
      cerr << "not convergence in this optimization" << endl;
      convergence = 0;
    }
  }

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void sensitivity(FEM &fem, vector<Stokes> &element, vector<Node> &node,
                 Optimize &optim, double V0)
{
  double jac = 1.0;

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * element[nel].ne);
    MatrixXd Be = MatrixXd::Zero(fem.voigt, element[nel].numdof);
    MatrixXd Ke = MatrixXd::Zero(element[nel].numdof, element[nel].numdof);
    VectorXd temp = VectorXd::Zero(element[nel].numdof);
    MatrixXd X = MatrixXd::Zero(element[nel].ne, fem.ndim);

    for(int i = 0; i < element[nel].ne; i++)
    {
      for(int j = 0; j < fem.dofnp; j++)
        temp[fem.dofnp * i + j] = node[element[nel].nodeID[i]].val[j];
      for(int j = 0; j < fem.ndim; j++)
        X.coeffRef(i, j) = node[element[nel].nodeID[i]].x[j];
    }

    Bmatrix bmatrix;
    for(int ip = 0; ip < element[nel].ipmax; ip++)
    {
      bmatrix.make(Be, Ne, jac, element[nel].eType, element[nel].ipmax, X, ip);
      //  Ke += Be.transpose() * dkds * Be * jac;
    }

    optim.dfds[nel] = -temp.transpose() * Ke * temp;
    optim.dgds[nel] = element[nel].volume / V0;
  }
}

void output(FEM &fem, vector<Stokes> &element, vector<Node> &node,
            Optimize &optim)
{
  string outFile = "res" + to_string(optim.optstep) + ".vtk";
  ofstream fout(outFile);
  makeHeaderVTK(fout);

  // make grid
  setPointVTK(node, fout);
  setElementVTK(element, fout);
  setElementTypeVTK(fem.ndim, element, fout);

  addPointVectorVTK("velocity", node, fout, true,
                    [&](Node &n)
                    {
                      if(fem.ndim == 3)
                        return n.val;
                      else
                      {
                        vector<double> three(3);
                        for(int i = 0; i < 2; i++)
                          three[i] = n.val[i];
                        return three.data();
                      }
                    });

  second_order_pressure = ;

  addPointScalarVTK("pressure", node, fout, false,
                    [&](Node &n) { return n.val[fem.ndim]; });

  addElementScalarVTK("design_s", element, fout, true,
                      [&](Stokes &e) { return e.design_s(); });

  fout.close();
}