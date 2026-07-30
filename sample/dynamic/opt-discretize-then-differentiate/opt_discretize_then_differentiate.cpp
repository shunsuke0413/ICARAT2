#include "opt_discretize_then_differentiate.hpp"
#include "fdm.hpp"
#include "sensitivity.hpp"

using namespace icarat;
using namespace std;
using namespace icarat::opt_discretize_then_differentiate;

void output(FEMDyna &fem, vector<DynamicElastic> &element, vector<Node> &node,
            int optstep, int time);

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <.toml file>\n";
    return 1;
  }

  // measure time
  auto t_start = chrono::system_clock::now();
  int convergence = 1;
  int numConst = 1;
  double V0 = 0.0;

  const auto config = toml::parse(argv[1]);

  FEMDyna fem;
  // parameters
  fem.gamma = 0.5; ///< newmark-beta parameter
  fem.beta = 0.25; ///< newmark-beta parameter
  fem.aa = 0;      ///< damping parameter
  fem.bb = 0;      ///< damping parameter
  fem.deltaT = toml::find<double>(config, "this_problem", "delta_t");
  fem.timestep = toml::find<int>(config, "this_problem", "timestep");
  // material density (kg/mm3)
  double density = toml::find<double>(config, "material", "density");
  bool fdm_is = toml::find<bool>(config, "this_problem", "fdm");
  double filterR = toml::find<double>(config, "this_problem", "filterR");

  // set problem
  MeshIcarat<DynamicElastic, Node> mesh(fem, config);
  mesh.setParameter("structure");

  // set element
  vector<DynamicElastic> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      DynamicElastic actele(fem.voigt, mesh.ne(type), mesh.ipmax(type), numdof,
                            density);
      actele.setParameter(config);
      element.push_back(actele);
    }
  }

  // set node
  vector<Node> node(fem.numnp);

  // input mesh to element & node
  mesh.generate(element, node);

  // boundary condition class
  BCIcarat<DynamicElastic, Node> BC(fem, element, node);
  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  ValueDyna dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  // set optimization class
  Optimize optim(fem.nelx, numConst, 1.0e-4, 1.0);
  MMASolver mma(fem.nelx, numConst);
  Filter<DynamicElastic, Node> filter(fem, element, node, filterR);

  // analysis class
  ImplicitDyna<DynamicElastic, Node> compt(fem, element, node, force, dis);

  objective_function::disp_squared_norm function(fem.numeq, fem.deltaT,
                                                 force.fext_org);
  Sensitivity<objective_function::disp_squared_norm> sensitivity(
      fem, element, node, force, dis, optim, function, compt.K(), compt.M());

  // optimization init
  optim.optstep = 0;
  for(int i = 0; i < fem.nelx; i++)
    optim.design_s[i] = element[i].design_s();

  // optimization start
  while(convergence == 1)
  {
    optim.optstep++;

    // initialize
    function.init();

    cout << "============== optimization step ==============" << endl;
    cout << "step no. " << optim.optstep << endl;

    auto design_filtered = filter.densityFilter(optim.design_s, false);

    /// set threshold function
    double beta = min(pow(2, (int)(optim.optstep / 40)), 16.0);
    cout << "beta" << beta << endl;
    auto design_threshold = filter.threshold(design_filtered, beta, 0.5);
    auto threDiff = filter.thresholdDerivative(design_filtered, beta, 0.5);

    // set threshold design variable
    for(int i = 0; i < fem.nelx; i++)
      element[i].setDesignS(design_threshold[i]);

    // compute dynamic analysis
    compt.initialization();
    for(int time = 1; time < fem.timestep + 1; time++)
    {
      force.fext = time * force.fext_org / fem.timestep;

      cout << "---------- timestep no. " << time << " ----------" << endl;

      compt.calRayleigh(time, fem.timestep, fem.aa, fem.bb, fem.deltaT,
                        fem.beta, fem.gamma, "Eigen_CG", 0);

      // update mises stress and store internal force
      for(auto &e : element)
      {
        e.makeMisesStress(fem, node);
        sensitivity.fint_h[time][e.ID] = e.fe();
      }
      // store global vectors
      sensitivity.disp_h[time] = dis.val;
      sensitivity.velo_h[time] = dis.velocity;
      sensitivity.acce_h[time] = dis.accel;

      // compute objective function
      function.update(dis.val);

      output(fem, element, node, optim.optstep, time);
    }
    optim.object_f = function.getValue();
    cout << "function value = " << optim.object_f << endl;

    // get const_func
    if(optim.optstep == 1)
    {
      V0 = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        V0 += element[i].design_s() * element[i].volume;
      if(V0 == 0.0)
        cerr << "V0 = 0. in optimization" << endl;

      optim.const_h[0] = 0.0;
    }
    else
    {
      optim.const_h[0] = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        optim.const_h[0] += element[i].design_s() * element[i].volume / V0;
      optim.const_h[0] -= 1.0;
    }
    cout << "const value: " << optim.const_h[0] << ",  "
         << "first volume: " << V0 << endl;

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);

    // convergence check
    convergence = optim.isConvergence(1.0e-5, 50);

    // analytical
    sensitivity.adjointAnalysis();

    // filtering sensitivity
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

    if(fdm_is)
      fdm(fem, element, node, force, dis, optim, function, filter);

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

void output(FEMDyna &fem, vector<DynamicElastic> &element, vector<Node> &node,
            int optstep, int time)
{
  ofstream fout("res" + to_string(optstep) + "_" + to_string(time) + ".vtk");
  makeHeaderVTK(fout);

  // make grid
  setPointVTK(node, fout);
  setElementVTK(element, fout);
  setElementTypeVTK(fem.ndim, element, fout);

  // export vector data to node
  addPointVectorVTK("displacement", node, fout, true,
                    [&](Node &n) { return n.val; });

  // export scalar data to element
  addElementScalarVTK("design_s", element, fout, true,
                      [&](DynamicElastic &e) { return e.design_s(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](DynamicElastic &e) { return e.mises(); });

  fout.close();
}
