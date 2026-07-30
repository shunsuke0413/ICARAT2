#include "opt_stress_nonlinear.hpp"

using namespace std;
using namespace Eigen;
using namespace icarat;
using namespace icarat::stress_nonlinear;

void output(FEMStressConst &fem, std::vector<FiniteSaintVenantTotal> &element,
            std::vector<Node> &node, Optimize &optim);

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <.toml file>\n";
    return 1;
  }

  auto t_start = chrono::system_clock::now(); // measure time
  int convergence = 1;                        // convergence judge
  int numConst = 2;                           // number of constraint
  double V0 = 0.0;                            // first volume
  double tolNR = 1.0e-8;
  int iterNR = 100;
  double beta1 = 50.0;
  double rho0 = 0.1;
  double xmin = 1.0e-4;
  double xmax = 1.0;

  // set input file
  const auto config = toml::parse(argv[1]);
  int numstep = toml::find<int>(config, "this_problem", "numstep");

  // set problem
  FEMStressConst fem;
  fem.stressLim = toml::find<double>(config, "this_problem",
                                     "stress_limit"); // stress limit
  fem.pNorm =
      toml::find<double>(config, "this_problem", "p_norm"); // p-norm parameter
  fem.qp = toml::find<double>(config, "this_problem",
                              "qp"); // parameter for qp-relaxation

  MeshIcarat<FiniteSaintVenantTotal, Node> mesh(fem, config);
  mesh.setParameter("structure");

  vector<FiniteSaintVenantTotal> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      FiniteSaintVenantTotal actele(fem.voigt, mesh.ne(type), mesh.ipmax(type),
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
  BCIcarat<FiniteSaintVenantTotal, Node> BC(fem, element, node);

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  Value dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  Optimize optim(fem.nelx, numConst, xmin, xmax);

  MMASolver mma(fem.nelx, numConst);

  double filterR = 2.0 * mesh.length()[0] / mesh.div()[0];
  Filter<FiniteSaintVenantTotal, Node> filter(fem, element, node, filterR);

  // init
  optim.optstep = 0;
  for(int i = 0; i < fem.nelx; i++)
    optim.design_s[i] = element[i].design_s();

  NewtonRaphson<FiniteSaintVenantTotal, Node> compt(fem, element, node, force,
                                                    dis);

  // optimization start
  while(convergence == 1)
  {
    optim.optstep++;

    cout << "============== optimization step ==============" << endl;
    cout << "step no. " << optim.optstep << endl;

    // get filtered density
    auto design_filtered = filter.densityFilter(optim.design_s, false);
    double beta = min(pow(2.0, (int)(optim.optstep / 40)), 8.0);
    double T = 0.5;
    auto design_threshold = filter.threshold(design_filtered, beta, T);
    auto threDifferent = filter.thresholdDerivative(design_filtered, beta, T);

    auto gamma = filter.threshold(design_threshold, beta1, rho0);

    // set threshold design variable
    for(int i = 0; i < fem.nelx; i++)
    {
      design_threshold[i] =
          std::max(std::min(design_threshold[i], xmax), xmin); // clamp
      element[i].setDesignS(design_threshold[i]);
      gamma[i] = std::max(std::min(gamma[i], xmax), xmin); // clamp
      element[i].setGamma(1.0);
    }

    compt.initialization();

    for(int loadStep = 1; loadStep < numstep + 1; loadStep++)
    {
      force.fext = loadStep * force.fext_org / numstep;

      cout << "---------- load-step No. " << loadStep << " ----------" << endl;

      compt.solve(loadStep, numstep, tolNR, iterNR, "Eigen_LDLT", 0);
    }

    // get stress
    for(auto &e : element)
      e.makeMisesStress(fem, node);

    // get objective function
    optim.object_f = dis.val.dot(force.fext_org);
    std::cout << "function value: " << optim.object_f << endl;

    // get volume constraint function
    if(optim.optstep == 1)
    {
      V0 = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        V0 += element[i].design_s() * element[i].volume;

      cout << "first volume: " << V0 << endl;
      optim.const_h[0] = 0.0;
    }
    else
    {
      optim.const_h[0] = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        optim.const_h[0] += element[i].design_s() * element[i].volume;

      optim.const_h[0] = optim.const_h[0] / V0 - 1.0;
    }
    cout << "Constraint0 value: " << optim.const_h[0] << endl;

    // get stress constraint function
    optim.const_h[1] = funcStressConst(fem, element);
    cout << "Constraint1 value: " << optim.const_h[1] << endl;

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);
    output(fem, element, node, optim);

    convergence = optim.isConvergence(1.0e-6, 50);

    //  analytical sensitivity
    auto fil_gamma = filter.thresholdDerivative(design_threshold, beta1, rho0);
    auto dfds =
        sensCompliance(fem, element, node, force, dis, compt.K(), fil_gamma);
    auto dg2ds = sensStressConst(fem, element, node, compt.K());

    // add threshold derivative
    vector<double> dg1ds(fem.nelx);
    for(int i = 0; i < fem.nelx; i++)
    {
      optim.dfds[i] = dfds[i] * threDifferent[i];
      dg1ds[i] = element[i].volume * threDifferent[i] / V0;
      dg2ds[i] = dg2ds[i] * threDifferent[i];
    }

    // add filter's derivative
    optim.dfds = filter.densityFilter(optim.dfds, true);
    dg1ds = filter.densityFilter(dg1ds, true);
    dg2ds = filter.densityFilter(dg2ds, true);

    // output stress constraint's sensivity graph
    if(optim.optstep == 1)
      exportGraph(dg2ds, "objective.csv");

    // FDM
    if(toml::find<bool>(config, "this_problem", "fdm"))
      fdm(fem, element, node, force, dis, optim, filter, numstep, tolNR,
          iterNR);

    for(int i = 0; i < fem.nelx; i++)
    {
      optim.dgds[2 * i] = dg1ds[i];
      optim.dgds[2 * i + 1] = dg2ds[i];
    }

    // Get updated design variables with MMA
    mma.Update(optim.design_s.data(), optim.dfds.data(), optim.const_h.data(),
               optim.dgds.data(), optim.xmin.data(), optim.xmax.data());

    // not convergence
    if(optim.optstep == 1000)
    {
      cerr << "not convergence in this optimization" << endl;
      convergence = 1;
    }
  }

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void output(FEMStressConst &fem, vector<FiniteSaintVenantTotal> &element,
            vector<Node> &node, Optimize &optim)
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

  // export scalar data to element
  addElementScalarVTK("density", element, fout, true,
                      [&](FiniteSaintVenantTotal &e) { return e.design_s(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](FiniteSaintVenantTotal &e) { return e.mises(); });

  fout.close();
}
