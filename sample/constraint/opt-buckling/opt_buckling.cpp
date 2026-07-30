/// This implement is based on the below paper.
/// Mizutori's master thesis (2021)
/// @sa Gao X., Ma H. Topology optimization of continuum structures under
/// buckling constraints. Computers & Structures, Volume 157, 2015, pp.142-152
#include "opt_buckling.hpp"
using namespace std;
using namespace icarat;
using namespace Eigen;
using namespace icarat::opt_buckling;

void output(FEMBuckling &fem, vector<BucklingLinearElastic> &element,
            vector<Node> &node, ValueEigen &value, Optimize &optim);

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
  double V0 = 0.0;     // first volume
  double smin = 1.0e-4;
  double smax = 1.0;

  const auto config = toml::parse(argv[1]);

  bool isConst = toml::find<bool>(config, "this_problem", "constraint");

  int numConst; // number of constraint
  if(isConst)
    numConst = 2;
  else
    numConst = 1;

  // set problem
  FEMBuckling fem;
  fem.numeigen = toml::find<int>(config, "this_problem", "num_eigen");
  fem.pNorm = toml::find<double>(config, "this_problem", "p_norm");
  fem.limeigen = toml::find<double>(config, "this_problem", "lim_eigen") * toml::find<double>(config, "this_problem", "alpha");
  fem.sigma = toml::find<double>(config, "this_problem", "shift");
  fem.paraconv = toml::find<int>(config, "this_problem", "shift_n") * fem.numeigen;

  MeshIcarat<BucklingLinearElastic, Node> mesh(fem, config);
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

  double filterR = 2.0 * mesh.length()[0] / mesh.div()[0];
  FilterHelmholtz<BucklingLinearElastic, Node> filter(fem, element, node,
                                                      filterR);

  // analysis class
  LinearAnalysis<BucklingLinearElastic, Node> linear(fem, element, node, force,
                                                     dis);
  fem.paraconv = toml::find<int>(config, "this_problem", "shift_n") * fem.numeigen;
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
    double beta = min(pow(2, (int)(optim.optstep / 40)), 16.0);
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
    if(isConst)
      eigen.solve();

    // get objective function
    optim.object_f = dis.val.dot(force.fext);
    std::cout << "function value: " << optim.object_f << endl;

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
      optim.const_h[0] -= 1;
      std::cout << "constraint1: " << optim.const_h[0] << endl;
    }

    if(isConst)
    {
      optim.const_h[1] = func_buckling(fem, dis);
      std::cout << "constraint2: " << optim.const_h[1] << endl;
    }

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);
    output(fem, element, node, dis, optim);

    convergence = optim.isConvergence(1.0e-6, 250);

    // analytical
    vector<double> dg1ds(fem.nelx), dg2ds(fem.nelx);

    optim.dfds = sens_compliance(fem, element, node);

    dg1ds = sens_volume(fem, element, V0);

    if(isConst)
      dg2ds = sens_buckling(fem, element, node, linear, dis);

    // filter's sensitivity term
    auto threDiff = filter.thresholdDerivative(design_filtered, beta, T);
    for(int i = 0; i < fem.nelx; i++)
    {
      optim.dfds[i] *= threDiff[i];
      dg1ds[i] *= threDiff[i];
      dg2ds[i] *= threDiff[i];
    }

    optim.dfds = filter.densityFilter(optim.dfds, fem.solver, 0, true);
    dg1ds = filter.densityFilter(dg1ds, fem.solver, 0, true);
    if(isConst)
      dg2ds = filter.densityFilter(dg2ds, fem.solver, 0, true);

    // output sensivity graph
    if(optim.optstep == 1)
      exportGraph(dg2ds, "objective.csv");

    // fdm
    if(toml::find<bool>(config, "this_problem", "fdm"))
      fdm(fem, element, node, force, dis, optim, filter);

    // delete sensitivities in non-design area
    for(int i = 0; i < fem.nelx; i++)
    {
      if(element[i].isNonDesign == true)
      {
        optim.dfds[i] = 0.0;
        dg1ds[i] = 0.0;
        dg2ds[i] = 0.0;
      }
    }

    // make constraint sensitivity's array for MMA
    if(isConst)
      for(int i = 0; i < fem.nelx; i++)
      {
        optim.dgds[numConst * i] = dg1ds[i];
        optim.dgds[numConst * i + 1] = dg2ds[i];
      }
    else
      optim.dgds = dg1ds;

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

  // eigenvalue of final design
  if(!isConst)
    eigen.solve();

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