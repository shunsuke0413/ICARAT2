#include "meta_material.hpp"

using namespace std;
using namespace icarat;
using namespace icarat::meta_material;
using namespace Eigen;

void convertVoigtToTensor(FEM &fem, const VectorXd &voigt, MatrixXd &tensor);
void output(int optstep, FEM &fem, vector<mHomoElastic> &element,
            vector<HomoNode> &node, ValueHomo &value);

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <.toml file>\n";
    return 1;
  }

  auto t_start = chrono::system_clock::now(); // measure time

  // set input file
  const auto config = toml::parse(argv[1]);
  FEM fem;

  bool fdm_flag = toml::find<bool>(config, "this_problem", "fdm");
  int convergence = 1;
  int numConst = 1; // number of constraint
  double xmin = 1.0e-4;
  double xmax = 1.0;
  double V0 = 0.0;     // initial volume of solid
  double Vtotal = 0.0; // total volume of structure

  double tYoung = toml::find<double>(config, "this_problem", "target_young");
  double tPoisson =
      toml::find<double>(config, "this_problem", "target_poisson");

  // set problem & mesh
  MeshIcarat<mHomoElastic, HomoNode> mesh(fem, config);
  mesh.setParameter("structure");

  MatrixXd omega = MatrixXd::Zero(fem.voigt, fem.voigt);

  if(fem.ndim == 2)
    omega << 1.0, 1.0, 0.0, //
        1.0, 1.0, 0.0,      //
        0.0, 0.0, 0.0;
  else
    omega << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, //
        1.0, 1.0, 1.0, 0.0, 0.0, 0.0,      //
        1.0, 1.0, 1.0, 0.0, 0.0, 0.0,      //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,      //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,      //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

  // set element
  vector<mHomoElastic> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      mHomoElastic actele(fem.voigt, mesh.ne(type), mesh.ipmax(type), numdof);
      actele.setParameter(config);
      element.push_back(actele);
    }
  }

  // set node
  vector<HomoNode> node(fem.numnp);

  // input mesh to element & node
  mesh.generate(element, node);
  auto pairs = mesh.generateMasterSlavePairs(node);
  auto length = mesh.length();

  // boundary condition class
  BCIcarat<mHomoElastic, HomoNode> bc(fem, element, node);

  // set periodic dirichlet condition
  fem.numeq = bc.makeDOF_periodic(pairs);

  // set displacement
  ValueHomo dis(fem.numeq, fem.voigt);

  Homogenization<mHomoElastic, HomoNode> mCompt(fem, element, node, dis, pairs,
                                                length);

  Optimize optim(fem.nelx, numConst, xmin, xmax);
  double filterR = 2.0 * abs(node[0].x[0] - node[1].x[0]);
  FilterHelmholtz<mHomoElastic, HomoNode> filter(fem, element, node, filterR);

  MMASolver mma(fem.nelx, numConst);

  // set initial design variable
  optim.optstep = 0;
  optim.design_s =
      mesh.setDesignHomo(element, node, element[0].design_s(), 0.2, 0.15);

  // target elastic tensor
  MatrixXd tCH = MatrixXd::Zero(fem.voigt, fem.voigt);
  element[0].De(tCH, tYoung, tPoisson);

  cout << "target De:" << endl << tCH << endl;

  while(convergence == 1)
  {
    optim.optstep++;

    cout << "============== optimization step ==============" << endl;
    cout << "step no. " << optim.optstep << endl;

    // filtering procedure
    auto filtered_s =
        filter.densityFilter(optim.design_s, "Eigen_CG", 0, false);
    double beta = min(pow(2, (double)((optim.optstep) / 40)), 16.0);
    cout << "beta" << beta << endl;
    double T = 0.5;
    auto threshold_s = filter.threshold(filtered_s, beta, T);
    auto threDiff = filter.thresholdDerivative(filtered_s, beta, T);

    for(int i = 0; i < fem.nelx; i++)
      element[i].setDesignS(threshold_s[i]);

    // numerical material testing
    MatrixXd CH = mCompt.solveNMT(true);
    cout << "CH:" << endl << CH << endl;

    for(int i = 0; i < fem.nelx; i++)
      element[i].makeMicroValues(fem, node);

    optim.object_f = getObjectF(fem, tCH, CH, omega);
    cout << "function value = " << optim.object_f << endl;

    if(optim.optstep == 1)
    {
      V0 = 0.0;
      Vtotal = 0.0;
      for(int i = 0; i < fem.nelx; i++)
      {
        V0 += element[i].design_s() * element[i].volume;
        Vtotal += element[i].volume;
      }
      if(V0 == 0.0)
        cerr << "V0 = 0. in optimization" << endl;
      optim.const_h[0] = V0;
    }
    else
    {
      optim.const_h[0] = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        optim.const_h[0] += element[i].design_s() * element[i].volume / V0;
      optim.const_h[0] -= 1.0;
    }
    cout << "Volume error: " << optim.const_h[0] << endl;

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);
    output(optim.optstep, fem, element, node, dis);

    convergence = optim.isConvergence(1.0e-5, 200);

    sensitivity(fem, element, node, dis, optim, tCH, CH, omega, V0, Vtotal);

    // threshold's derivative
    for(int i = 0; i < fem.nelx; i++)
    {
      optim.dfds[i] *= threDiff[i];
      optim.dgds[i] *= threDiff[i];
    }
    // filter's derivative
    optim.dfds = filter.densityFilter(optim.dfds, "Eigen_CG", 0, true);
    optim.dgds = filter.densityFilter(optim.dgds, "Eigen_CG", 0, true);

    // output sensivity graph
    if(optim.optstep == 1)
      exportGraph(optim.dfds, "objective.csv");

    if(fdm_flag)
      fdm(config, fem, element, node, dis, optim, filter, omega, pairs, length);

    mma.Update(optim.design_s.data(), optim.dfds.data(), optim.const_h.data(),
               optim.dgds.data(), optim.xmin.data(), optim.xmax.data());

    // not convergence
    if(optim.optstep == 1000)
    {
      cerr << "not convergence in this optimization" << endl;
      exit(1);
    }
  }

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void output(int optstep, FEM &fem, vector<mHomoElastic> &element,
            vector<HomoNode> &node, ValueHomo &value)
{
  ofstream fout("res" + to_string(optstep) + ".vtk");
  makeHeaderVTK(fout);

  // make grid
  setPointVTK(node, fout);
  setElementVTK(element, fout);
  setElementTypeVTK(fem.ndim, element, fout);

  // export vector data to node
  addPointScalarVTK("ID", node, fout, true, [&](Node &n) { return n.ID; });

  // micro strain in each direction
  for(int d = 0; d < fem.voigt; d++)
    addPointVectorVTK("direction" + to_string(d), node, fout, false,
                      [&](HomoNode &n) { return n.nmtval[d]; });

  // export scalar data to element
  addElementScalarVTK("density", element, fout, true,
                      [&](mHomoElastic &e) { return e.design_s(); });

  addElementScalarVTK("young", element, fout, false,
                      [&](mHomoElastic &e) { return e.young(); });

  fout.close();
}
