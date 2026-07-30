///@sa T.Buhl et al., Stiffness design of geometrically nonlinear
/// structures using topology optimization, SMO, 2000.
///@sa F.Wang et al., Interpolation scheme for fictitious domain techniques and
/// topology optimization of finite strain elastic problems, CMAME, 2014.
#include "opt_finite_strain.hpp"
using namespace icarat;
using namespace icarat::opt_finite_strain;
using namespace std;
using namespace Eigen;

void sensitivity(FEM &fem, vector<FiniteNeoHookeTotal> &element,
                 vector<Node> &node, Force &force, Value &value,
                 const Eigen::SparseMatrix<double> &K, Optimize &optim,
                 double V0, vector<double> &fil_gamma);
void output(FEM &fem, vector<FiniteNeoHookeTotal> &element, vector<Node> &node,
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
  double xmin = 1.0e-4;
  double xmax = 1.0;
  double V0 = 0.0; // first volume
  double eta = 0.5;
  double beta1 = 50.0;
  double rho0 = 0.1;
  double tolNR = 1.0e-6;
  int iterNR = 100;

  // set input file
  const auto config = toml::parse(argv[1]);

  int numstep = toml::find<int>(config, "this_problem", "numstep");

  // set problem
  FEM fem;
  MeshIcarat<FiniteNeoHookeTotal, Node> mesh(fem, config);
  mesh.setParameter("structure");

  // set element
  vector<FiniteNeoHookeTotal> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      FiniteNeoHookeTotal actele(fem.voigt, mesh.ne(type), mesh.ipmax(type),
                                 numdof);
      actele.setParameter(config);
      element.push_back(actele);
    }
  }

  // set node
  vector<Node> node(fem.numnp);

  // input mesh to element & node
  mesh.generate(element, node);

  // boundary condition class
  BCIcarat<FiniteNeoHookeTotal, Node> BC(fem, element, node);

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  ValueNR dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  // set optimization's classes
  Optimize optim(fem.nelx, numConst, xmin, xmax);
  MMASolver mma(fem.nelx, numConst);
  double filterR = toml::find<double>(config, "this_problem", "filterR");
  Filter<FiniteNeoHookeTotal, Node> filter(fem, element, node, filterR);

  optim.optstep = 0;
  /// filtering procedure
  for(int i = 0; i < fem.nelx; i++)
    optim.design_s[i] = element[i].design_s();

  NewtonRaphson<FiniteNeoHookeTotal, Node> compt(fem, element, node, force,
                                                 dis);

  // optimization start
  while(convergence == 1)
  {
    optim.optstep++;

    std::cout << "============== optimization step ==============" << endl;
    std::cout << "step no. " << optim.optstep << endl;

    auto design_filtered = filter.densityFilter(optim.design_s, false);

    /// set threshold function
    double beta = min(1.0 + (double)(optim.optstep / 20), 16.0);
    std::cout << "beta" << beta << std::endl;

    auto design_threshold = filter.threshold(design_filtered, beta, eta);

    auto gamma = filter.threshold(design_threshold, beta1, rho0);

    // set threshold design variable
    for(int i = 0; i < fem.nelx; i++)
    {
      design_threshold[i] =
          std::max(std::min(design_threshold[i], xmax), xmin); // clamp
      element[i].setDesignS(design_threshold[i]);
      gamma[i] = std::max(std::min(gamma[i], xmax), xmin); // clamp
      element[i].setGamma(gamma[i]);
    }

    compt.initialization();

    for(int loadStep = 1; loadStep < numstep + 1; loadStep++)
    {
      force.fext = loadStep * force.fext_org / numstep;

      cout << "---------- load-step No. " << loadStep << " ----------" << endl;

      compt.solve(loadStep, numstep, tolNR, iterNR, "Eigen_LDLT", 0);
    }

    for(auto &e : element)
      e.makeMisesStress(fem, node);

    // get const_func
    if(optim.optstep == 1)
    {
      V0 = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        V0 += element[i].design_s() * element[i].volume;
      if(V0 == 0.0)
        cerr << "V0 = 0. in optimization" << endl;
    }
    else
    {
      optim.const_h[0] = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        optim.const_h[0] += element[i].design_s() * element[i].volume / V0;
      optim.const_h[0] -= 1;
      std::cout << "Constraint value:" << optim.const_h[0] << endl;
    }

    // get objective function
    optim.object_f = dis.val.dot(force.fext_org);
    std::cout << "function value = " << optim.object_f << endl;

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);
    output(fem, element, node, optim);

    convergence = optim.isConvergence(1.0e-5, 350);

    // analytical
    auto fil_gamma = filter.thresholdDerivative(design_threshold, beta1, rho0);
    sensitivity(fem, element, node, force, dis, compt.K(), optim, V0,
                fil_gamma);
    auto threDiff = filter.thresholdDerivative(design_filtered, beta, eta);
    for(int i = 0; i < fem.nelx; i++)
    {
      optim.dfds[i] *= threDiff[i];
      optim.dgds[i] *= threDiff[i];
    }
    optim.dfds = filter.densityFilter(optim.dfds, true);
    optim.dgds = filter.densityFilter(optim.dgds, true);

    // output sensivity graph
    if(optim.optstep == 1)
      exportGraph(optim.dfds, "objective.csv");

    // fdm
    if(toml::find<bool>(config, "this_problem", "fdm"))
      fdm(fem, element, node, force, dis, optim, filter, numstep);

    // Get updated design variables with MMA
    mma.Update(optim.design_s.data(), optim.dfds.data(), optim.const_h.data(),
               optim.dgds.data(), optim.xmin.data(), optim.xmax.data());

    // not convergence
    if(optim.optstep == 1000)
    {
      cerr << "not convergence in this optimization" << endl;
      convergence = 1;
    }

  } // optimization loop

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void sensitivity(FEM &fem, vector<FiniteNeoHookeTotal> &element,
                 vector<Node> &node, Force &force, Value &value,
                 const Eigen::SparseMatrix<double> &K, Optimize &optim,
                 double V0, vector<double> &fil_gamma)
{

  // get adjoint vector
  std::cout << "---- computing adjoint vector ----" << endl;
  Solvers solver(K, force.fext, value.val);
  solver.solve("Eigen_CG", 0);

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    double dens = element[nel].design_s();
    double young1 = element[nel].young1();
    double young2 = element[nel].young2();
    double p = element[nel].pp();

    assert(young2 >= young1);

    double young = (1.0 - pow(dens, p)) * young1 + pow(dens, p) * young2;
    double dyoung = p * pow(dens, p - 1.0) * (young2 - young1);

    double E = dyoung / young;

    /// Allocate final displacement and internal force to each element
    VectorXd valEle = VectorXd::Zero(element[nel].numdof);
    VectorXd disEle = VectorXd::Zero(element[nel].numdof);
    for(int i = 0; i < element[nel].ne; i++)
      for(int j = 0; j < fem.ndim; j++)
      {
        int dof = node[element[nel].nodeID[i]].dof[j];

        if(dof < fem.numeq)
          valEle.coeffRef(fem.ndim * i + j) = value.val[dof];
        // dirichlet is set.
        else
          valEle.coeffRef(fem.ndim * i + j) =
              node[element[nel].nodeID[i]].dirich->val[j];

        disEle.coeffRef(fem.ndim * i + j) = node[element[nel].nodeID[i]].val[j];
      }

    VectorXd fe_rho = element[nel].gamma() * element[nel].nfe() +
                      (1.0 - element[nel].gamma() * element[nel].gamma()) *
                          element[nel].lfe();

    VectorXd fe_gamma = element[nel].nfe() +
                        element[nel].Ke() * disEle * element[nel].gamma() -
                        2.0 * element[nel].gamma() * element[nel].lfe();

    optim.dfds[nel] =
        -valEle.transpose() * (E * fe_rho + fe_gamma * fil_gamma[nel]);
    optim.dgds[nel] = element[nel].volume / V0;
  }
}

void output(FEM &fem, vector<FiniteNeoHookeTotal> &element, vector<Node> &node,
            Optimize &optim)
{
  ofstream fout("res" + to_string(optim.optstep) + ".vtk");
  makeHeaderVTK(fout);

  // make grid
  setPointVTK(node, fout);
  setElementVTK(element, fout);
  setElementTypeVTK(fem.ndim, element, fout);

  // export vector data to node
  addPointVectorVTK("displacement", node, fout, true,
                    [&](Node &n) { return n.val; });

  addPointVectorVTK("traction", node, fout, false,
                    [&](Node &n) { return n.du; });

  // export scalar data to element
  addElementScalarVTK("density", element, fout, true,
                      [&](FiniteNeoHookeTotal &e) { return e.design_s(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](FiniteNeoHookeTotal &e) { return e.mises(); });

  fout.close();
}
