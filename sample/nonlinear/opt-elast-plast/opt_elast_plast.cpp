#include "opt_elast_plast.hpp"

using namespace std;
using namespace Eigen;
using namespace icarat;
using namespace icarat::opt_elast_plast;

void output(FEM &fem, vector<OptLinearElastPlast> &element, vector<Node> &node,
            Optimize &optim);

int main(int argc, char *argv[])
{
  if(argc == 1)
  {
    std::cout << "usage1: " << argv[0] << " <.toml file>\n";
    std::cout << "usage2: " << argv[0] << " <.toml file> <.mdpa file>\n";
    return 1;
  }

  // measure time
  auto t_start = chrono::system_clock::now();
  int convergence = 1; // convergence judge
  int numConst = 1;    // number of constraint
  double V0 = 0.0;     // first volume
  double T = 0.5;
  // set input file
  const auto config = toml::parse(argv[1]);
  int numStep = toml::find<int>(config, "this_problem", "numstep");
  double filterR = toml::find<double>(config, "this_problem", "filterR");

  FEM fem;
  MeshIcarat<OptLinearElastPlast, Node> mesh(fem, config);
  mesh.setParameter("structure");

  // set element(unstructured mesh version)
  vector<OptLinearElastPlast> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      OptLinearElastPlast actele(fem.voigt, mesh.ne(type), mesh.ipmax(type),
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
  BCIcarat<OptLinearElastPlast, Node> BC(fem, element, node);

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  ValueNR dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  Optimize optim(fem.nelx, numConst, 1.0e-4, 1.0);
  OC oc(fem.nelx);
  Filter<OptLinearElastPlast, Node> filter(fem, element, node, filterR);

  // initialization
  optim.optstep = 0;
  for(int i = 0; i < fem.nelx; i++)
    optim.design_s[i] = element[i].design_s();

  NewtonRaphson<OptLinearElastPlast, Node> compt(fem, element, node, force,
                                                 dis);

  // optimization start
  while(convergence == 1)
  {
    optim.optstep++;

    cout << "============== optimization step ==============" << endl;
    cout << "step no. " << optim.optstep << endl;

    auto design_filtered = filter.densityFilter(optim.design_s, false);

    /// set threshold function
    double beta = 1.0;
    cout << "beta: " << beta << endl;
    auto design_threshold = filter.threshold(design_filtered, beta, T);
    auto threDeriv = filter.thresholdDerivative(design_filtered, beta, T);

    // set threshold design variable
    for(int i = 0; i < fem.nelx; i++)
      element[i].setDesignS(design_threshold[i]);

    // Initialization of element internal variables
    compt.initialization();
    for(int nel = 0; nel < fem.nelx; nel++)
    {
      element[nel].initialization();
      element[nel].init_opt();
    }

    for(int loadStep = 1; loadStep < numStep + 1; loadStep++)
    {
      for(auto &e : element)
        e.iteration = 0;
      force.fext = loadStep * force.fext_org / numStep;

      cout << "---------- load-step No. " << loadStep << " ----------" << endl;

      compt.solve(loadStep, numStep, 1.0e-8, 10, "Eigen_CG", 0);

      vector<Triplet<double>> tripletsR;
      for(int nel = 0; nel < fem.nelx; nel++)
      {
        element[nel].stress_plastic(fem, node, tripletsR);
        element[nel].sensitivity(fem, node);
      }
    }

    // get objective function
    optim.object_f = 0.0;
    for(int i = 0; i < fem.nelx; i++)
      optim.object_f += element[i].func_new();

    cout << "object function: " << optim.object_f << endl;

    // get const_func
    if(optim.optstep == 1)
    {
      V0 = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        V0 += element[i].design_s() * element[i].volume;
      if(V0 == 0.0)
        cerr << "V0 = 0. in optimization" << endl;

      optim.const_h[0] = V0;
    }
    else
    {
      optim.const_h[0] = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        optim.const_h[0] += element[i].design_s() * element[i].volume;

      cout << "const value: " << optim.const_h[0] << ",  "
           << "first volume: " << V0 << endl;
    }

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);
    output(fem, element, node, optim);

    convergence = optim.isConvergence(1.0e-7, 50);

    for(int i = 0; i < fem.nelx; i++)
    {
      optim.dfds[i] = element[i].objective() * threDeriv[i];
      optim.dgds[i] = element[i].volume * threDeriv[i];
    }
    // filtering sensitivity
    optim.dfds = filter.densityFilter(optim.dfds, true);
    optim.dgds = filter.densityFilter(optim.dgds, true);

    exportGraph(optim.dfds, "objective.csv");

    if(toml::find<bool>(config, "this_problem", "fdm"))
      fdm(fem, element, node, force, dis, optim, numStep, filterR);

    // Get updated design variables with OC
    oc.UpdateVariables(optim.design_s, optim.dfds, V0, optim.dgds, optim.xmin,
                       optim.xmax,
                       [&](vector<double> &xkp1)
                       {
                         double g = 0.0;
                         auto in_fils = filter.densityFilter(xkp1, false);
                         auto in_thres = filter.threshold(in_fils, beta, T);
                         for(int i = 0; i < fem.nelx; i++)
                           g += in_thres[i] * element[i].volume;

                         return g;
                       });

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

void output(FEM &fem, vector<OptLinearElastPlast> &element, vector<Node> &node,
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

  // export scalar data to element
  addElementScalarVTK("density", element, fout, true,
                      [&](OptLinearElastPlast &e) { return e.design_s(); });

  addElementScalarVTK("epstn", element, fout, false,
                      [&](OptLinearElastPlast &e) { return e.eqstn(); });

  addElementScalarVTK("epsilon", element, fout, false,
                      [&](OptLinearElastPlast &e) { return e.epsilon(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](OptLinearElastPlast &e) { return e.mises(); });

  fout.close();
}
