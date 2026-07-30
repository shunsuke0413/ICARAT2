#include "opt_macro_oc.hpp"

using namespace std;
using namespace icarat;
using namespace Eigen;

void sensitivity(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
                 Optimize &optim);
void output(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
            Optimize &optim);

int main(int argc, char *argv[])
{
  if(argc != 3)
  {
    std::cout << "usage: " << argv[0] << " <.toml file> <.mdpa file>\n";
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
  vector<string> inputMDPA = read_file(argv[2]);

  // set problem
  FEM fem;
  MeshMDPA<LinearElastic, Node> mesh(fem, inputMDPA, config);
  mesh.setParameter("structure");

  // set element(unstructured mesh version)
  vector<LinearElastic> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      LinearElastic actele(fem.voigt, mesh.ne(type), mesh.ipmax(type), numdof);
      actele.setParameter(config);
      element.push_back(actele);
    }
  }

  // set node
  vector<Node> node(fem.numnp);

  // input mesh to element & node
  mesh.generate(element, node);

  // boundary condition class
  BCIcarat<LinearElastic, Node> BC(fem, element, node);

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  Value dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  // set optimization's class
  Optimize optim(fem.nelx, numConst, 1.0e-4, 1.0);
  OC oc(fem.nelx);
  double filterR = 1.0;
  FilterHelmholtz<LinearElastic, Node> filter(fem, element, node, filterR);

  // initialization
  optim.optstep = 0;
  for(int i = 0; i < fem.nelx; i++)
    optim.design_s[i] = element[i].design_s();
  LinearAnalysis<LinearElastic, Node> compt(fem, element, node, force, dis);

  // optimization start
  while(convergence == 1)
  {
    optim.optstep++;

    cout << "============== optimization step ==============" << endl;
    cout << "step no. " << optim.optstep << endl;

    auto design_filtered =
        filter.densityFilter(optim.design_s, "Eigen_CG", 0, false);

    /// set threshold function
    double beta = min(pow(2, (int)(optim.optstep / 40)), 16.0);
    cout << "beta: " << beta << endl;
    auto design_threshold = filter.threshold(design_filtered, beta, T);
    auto threDeriv = filter.thresholdDerivative(design_filtered, beta, T);

    // set threshold design variable
    for(int i = 0; i < fem.nelx; i++)
      element[i].setDesignS(design_threshold[i]);

    compt.solve("Eigen_CG", 1);

    /// get  stress
    for(auto &e : element)
      e.makeMisesStress(fem, node);

    // get objective function
    optim.object_f = dis.val.dot(force.fext);
    cout << "function value = " << optim.object_f << endl;

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

    convergence = optim.isConvergence(1.0e-5, 50);

    // analytical
    sensitivity(fem, element, node, optim);

    for(int i = 0; i < fem.nelx; i++)
    {
      optim.dfds[i] *= threDeriv[i];
      optim.dgds[i] *= threDeriv[i];
    }

    // filtering sensitivity
    optim.dfds = filter.densityFilter(optim.dfds, "Eigen_CG", 0, true);
    optim.dgds = filter.densityFilter(optim.dgds, "Eigen_CG", 0, true);

    // Get updated design variables with OC
    oc.UpdateVariables(
        optim.design_s, optim.dfds, V0, optim.dgds, optim.xmin, optim.xmax,
        [&](vector<double> &xkp1)
        {
          double g = 0.0;
          auto in_fils = filter.densityFilter(xkp1, "Eigen_CG", 0, false);
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

void sensitivity(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
                 Optimize &optim)
{
  assert(element[0].young2() > element[0].young1());

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    double dens = element[nel].design_s();
    double young1 = element[nel].young1();
    double young2 = element[nel].young2();
    double p = element[nel].pp();

    double young = (1.0 - pow(dens, p)) * young1 + pow(dens, p) * young2;
    double dyoung = p * pow(dens, p - 1.0) * (young2 - young1);

    double E = dyoung / young;

    Eigen::VectorXd disp = Eigen::VectorXd::Zero(element[nel].numdof);
    for(int i = 0; i < element[nel].ne; i++)
      for(int j = 0; j < fem.ndim; j++)
        disp.coeffRef(fem.ndim * i + j) = node[element[nel].nodeID[i]].val[j];

    optim.dfds[nel] = -E * disp.transpose() * element[nel].Ke() * disp;

    optim.dgds[nel] = element[nel].volume;
  }

  // output sensivity graph
  if(optim.optstep == 1)
    exportGraph(optim.dfds, "objective.csv");
}

void output(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
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
                      [&](LinearElastic &e) { return e.design_s(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](LinearElastic &e) { return e.mises(); });

  fout.close();
}