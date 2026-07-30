#include "opt_stress_buckling.hpp"
using namespace std;
using namespace icarat;
using namespace Eigen;
using namespace icarat::opt_stress_buckling;

#define meshisICARAT false

void output(FEMBuckling &fem, vector<BucklingLinearElastic> &element,
            vector<Node> &node, ValueEigen &value, Optimize &optim);

int main(int argc, char *argv[])
{
  if(argc != 2 && argc != 3)
  {
    std::cout << "usage: " << argv[0] << " <.toml file>\n";
    std::cout << "usage: " << argv[0] << " <.toml file> <.vtk file>\n";
    return 1;
  }

  // measure time
  auto t_start = chrono::system_clock::now();
  int convergence = 1; // convergence judge
  double V0 = 0.0;     // first volume
  double smin = 1.0e-4;
  double smax = 1.0;

  const auto config = toml::parse(argv[1]);

  string pathVTK;
  if(argc == 3)
    pathVTK = argv[2];

  double filterR = 2.0 * toml::find<double>(config, "this_problem", "filterR");

  vector<bool> cnstflag(3);
  cnstflag[0] = true;
  cnstflag[1] = toml::find<bool>(config, "this_problem", "const_buckling");
  cnstflag[2] = toml::find<bool>(config, "this_problem", "const_stress");

  int numConst = 0;
  for(int i = 0; i < 3; i++)
    numConst += cnstflag[i];

  cout << "number of constraint: " << numConst << endl;

  // set problem
  FEMBuckling fem;
  fem.numeigen = toml::find<int>(config, "this_problem", "num_eigen");
  fem.pNorm = toml::find<double>(config, "this_problem", "p_norm");
  fem.limeigen = toml::find<double>(config, "this_problem", "lim_eigen") *
                 toml::find<double>(config, "this_problem", "alpha_eigen");
  fem.limstress = toml::find<double>(config, "this_problem", "lim_stress") *
                  toml::find<double>(config, "this_problem", "alpha_stress");
  fem.sigma = toml::find<double>(config, "this_problem", "shift");
  fem.paraconv =
      toml::find<int>(config, "this_problem", "shift_n") * fem.numeigen;

// ICARAT mesh
#if meshisICARAT
  MeshIcarat<BucklingLinearElastic, Node> mesh(fem, config);
#else
  // vtk
  MeshVTK<BucklingLinearElastic, Node> mesh(fem, pathVTK, config);
#endif
  mesh.setParameter("structure");

  // set element
  vector<BucklingLinearElastic> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      BucklingLinearElastic actele(fem.voigt, mesh.ne(type), mesh.ipmax(type),
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
  BCIcarat<BucklingLinearElastic, Node> BC(fem, element, node);

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  ValueEigen dis(fem.numeq, fem.numeigen);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  // set optimization's classes
  Optimize optim(fem.nelx, numConst, smin, smax);
  MMASolver mma(fem.nelx, numConst);
  mma.SetAsymptotes(0.1, 0.7, 1.2);

  FilterHelmholtz<BucklingLinearElastic, Node> filter(fem, element, node,
                                                      filterR);

  // analysis class
  LinearAnalysis<BucklingLinearElastic, Node> linear(fem, element, node, force,
                                                     dis);
  fem.paraconv =
      toml::find<int>(config, "this_problem", "shift_n") * fem.numeigen;
  GeneralEigen<BucklingLinearElastic, Node> eigen(
      fem, element, node, dis, fem.numeigen, fem.paraconv, fem.sigma);

  // initial setting of the optimization
  optim.optstep = 0;
  for(int i = 0; i < fem.nelx; i++)
    optim.design_s[i] = element[i].design_s();

  // optimization start
  while(convergence == 1)
  {
    optim.optstep++;

    std::cout << "============== optimization step ==============" << endl;
    std::cout << "step no. " << optim.optstep << endl;

    auto design_filtered =
        filter.densityFilter(optim.design_s, fem.solver, 0, false);

    // set threshold function
    double beta = min(pow(2, (int)(optim.optstep / 40)), 8.0);
    double T = 0.5;
    std::cout << "beta" << beta << std::endl;
    auto design_threshold = filter.threshold(design_filtered, beta, T);

    // set threshold design variable
    for(int i = 0; i < fem.nelx; i++)
      element[i].setDesignS(design_threshold[i]);

    // structural analysis
    linear.solve(fem.solverSt, 1);

    for(int i = 0; i < fem.nelx; i++)
      element[i].makeMisesStress(fem, node);

    // eigenvalue analysis
    if(cnstflag[1])
      eigen.solve();

    // get objective function
    optim.object_f = dis.val.dot(force.fext);
    std::cout << "function value: " << optim.object_f << endl;

    int counter = 0;
    // get const_func
    if(optim.optstep == 1)
    {
      V0 = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        V0 += element[i].design_s() * element[i].volume;
      if(V0 == 0.0)
        cerr << "V0 = 0. in optimization" << endl;

      optim.const_h[counter] = 0.0;
      counter++;
    }
    else
    {
      optim.const_h[counter] = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        optim.const_h[counter] +=
            element[i].design_s() * element[i].volume / V0;
      optim.const_h[counter] -= 1;
      std::cout << "volume const.: " << optim.const_h[0] << endl;
      counter++;
    }

    if(cnstflag[1])
    {
      optim.const_h[counter] = func_buckling(fem, dis);
      std::cout << "buckling const.: " << optim.const_h[counter] << endl;
      counter++;
    }

    if(cnstflag[2])
    {
      optim.const_h[counter] = func_stress(fem, element);
      std::cout << "stress const.: " << optim.const_h[counter] << endl;
      counter++;
    }

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);
    output(fem, element, node, dis, optim);

    convergence = optim.isConvergence(1.0e-6, 200);

    // analytical
    optim.dfds = sens_compliance(fem, element, node);

    vector<vector<double>> dgds(3);
    for(int i = 0; i < 3; i++)
      dgds[i].resize(fem.nelx);

    dgds[0] = sens_volume(fem, element, V0);

    if(cnstflag[1])
      dgds[1] = sens_buckling(fem, element, node, linear, dis);

    if(cnstflag[2])
      dgds[2] = sens_stress(fem, element, node, linear.K());

    // filter's sensitivity term
    auto threDiff = filter.thresholdDerivative(design_filtered, beta, T);
    for(int i = 0; i < fem.nelx; i++)
    {
      optim.dfds[i] *= threDiff[i];
      dgds[0][i] *= threDiff[i];
      dgds[1][i] *= threDiff[i];
      dgds[2][i] *= threDiff[i];
    }

    optim.dfds = filter.densityFilter(optim.dfds, fem.solver, 0, true);
    dgds[0] = filter.densityFilter(dgds[0], fem.solver, 0, true);

    if(cnstflag[1])
      dgds[1] = filter.densityFilter(dgds[1], fem.solver, 0, true);
    if(cnstflag[2])
      dgds[2] = filter.densityFilter(dgds[2], fem.solver, 0, true);

    // output sensivity graph
    if(optim.optstep == 1)
    {
      exportGraph(optim.dfds, "sens_comp.csv");
      exportGraph(dgds[0], "sens_volume.csv");
      exportGraph(dgds[1], "sens_buckling.csv");
      exportGraph(dgds[2], "sens_stress.csv");
    }
    // fdm
    if(toml::find<bool>(config, "this_problem", "fdm"))
      fdm(fem, element, node, force, dis, optim, filter);

    // make constraint sensitivity's array for MMA
    for(int i = 0; i < fem.nelx; i++)
    {
      int counter = 0;
      for(int j = 0; j < 3; j++)
      {
        if(cnstflag[j])
        {
          optim.dgds[numConst * i + counter] = dgds[j][i];
          counter++;
        }
      }
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

  if(!cnstflag[1])
  {
    eigen.solve();
    output(fem, element, node, dis, optim);
  }

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  std::cout << "Computation time: " << t_sec << " sec " << std::endl;

  return 0;
}

void output(FEMBuckling &fem, vector<BucklingLinearElastic> &element,
            vector<Node> &node, ValueEigen &value, Optimize &optim)
{
  ofstream fout("res" + to_string(optim.optstep) + ".vtk");
  makeHeaderVTK(fout);

  // make grid
  setPointVTK(node, fout);
  setElementVTK(element, fout);
  setElementTypeVTK(fem.ndim, element, fout);

  // export vector data to node
  addPointScalarVTK("ID", node, fout, true, [&](Node &n) { return n.ID; });
  addPointVectorVTK("displacement", node, fout, false,
                    [&](Node &n) { return n.val; });

  // make eigen mode's output
  for(int i = 0; i < fem.numeigen; i++)
    addPointVectorVTK("buckling_mode" + to_string(i), node, fout, false,
                      [&](Node &n)
                      {
                        vector<double> three(3);
                        for(int i = 0; i < fem.dofnp; i++)
                        {
                          // dirichlet is active.
                          if(n.dof[i] < fem.numeq)
                            three[i] = value.evectors[i](n.dof[i]);
                          else
                            three[i] = n.dirich->val[i];
                        }
                        return three;
                      });

  // export scalar data to element
  addElementScalarVTK("ID", element, fout, true,
                      [&](BucklingLinearElastic &e) { return e.ID; });

  addElementScalarVTK("design_s", element, fout, false,
                      [&](BucklingLinearElastic &e) { return e.design_s(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](BucklingLinearElastic &e) { return e.mises(); });

  fout.close();
}