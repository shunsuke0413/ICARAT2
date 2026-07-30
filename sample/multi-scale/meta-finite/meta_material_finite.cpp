#include "meta_material_finite.hpp"
#include <sys/stat.h>

using namespace std;
using namespace icarat;
using namespace icarat::meta_material_finite;
using namespace Eigen;

void convertVoigtToTensor(FEM &fem, const VectorXd &voigt, MatrixXd &tensor);
void output(int optstep, FEM &fem, vector<mHomoNeoHookeTotal> &element,
            vector<HomoNode> &node, ValueHomo &value);

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <.toml file>\n";
    return 1;
  }

  auto t_start = chrono::system_clock::now(); // measure time

  // make vtk folder
  string folder_name = ("./result_vtk");
  if(mkdir(folder_name.c_str(), 0755) == 0)
    std::cout << "making vtk folder done!" << endl;
  else
    std::cout << "making vtk folder fault" << endl;

  // set input file
  const auto config = toml::parse(argv[1]);
  FEM fem;

  bool fdm_flag = toml::find<bool>(config, "this_problem", "fdm");
  int convergence = 1;
  int numConst = 1; // number of constraint
  double xmin = 1.0e-4;
  double xmax = 1.0;
  double T = 0.5;      // threshold parameter
  double beta1 = 50.0; // threshold parameter for stabilize
  double rho0 = 0.1;   // threshold parameter for stabilize
  double V0 = 0.0;     // initial volume of solid
  double Vtotal = 0.0; // total volume of structure

  int numstep = toml::find<int>(config, "this_problem", "numstep");
  double MstrainFac =
      toml::find<double>(config, "this_problem", "Mstrain_factor");

  double tYoungIso =
      toml::find<double>(config, "this_problem", "target_young_iso");
  double tPoissonIso =
      toml::find<double>(config, "this_problem", "target_poisson_iso");

  auto tYoungOrtho = toml::find<std::array<double, 3>>(config, "this_problem",
                                                       "target_young_ortho");
  auto tPoissonOrtho = toml::find<std::array<double, 6>>(
      config, "this_problem", "target_poisson_ortho");

  string matdirection_ =
      toml::find<string>(config, "this_problem", "direction");

  // set problem & mesh
  MeshIcarat<mHomoNeoHookeTotal, HomoNode> mesh(fem, config);
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
  vector<mHomoNeoHookeTotal> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      mHomoNeoHookeTotal actele(fem.voigt, mesh.ne(type), mesh.ipmax(type),
                                numdof);
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
  BCIcarat<mHomoNeoHookeTotal, HomoNode> bc(fem, element, node);

  // set periodic dirichlet condition
  fem.numeq = bc.makeDOF_periodic(pairs);

  // set displacement
  ValueHomo dis(fem.numeq, fem.voigt);

  Homogenization<mHomoNeoHookeTotal, HomoNode> mCompt(fem, element, node, dis,
                                                      pairs, length);

  Optimize optim(fem.nelx, numConst, xmin, xmax);
  double filterR = 2.0 * abs(node[0].x[0] - node[1].x[0]);
  FilterHelmholtz<mHomoNeoHookeTotal, HomoNode> filter(fem, element, node,
                                                       filterR);

  MMASolver mma(fem.nelx, numConst);
  mma.SetAsymptotes(0.1, 0.7, 1.2);

  // set initial design variable
  optim.optstep = 0;
  optim.design_s =
      mesh.setDesignHomo(element, node, element[0].design_s(), 0.1, 0.15);
  // target elastic tensor
  MatrixXd tCH(fem.voigt, fem.voigt);
  if(matdirection_ == "isotropy")
  {
    tCH = isoDe(fem.voigt, tYoungIso, tPoissonIso, element[0].mattype());
    std::cout << "target De:" << endl << tCH << endl;
  }

  else if(matdirection_ == "orthotropy")
  {
    tCH = orthoDe(fem.voigt, tYoungOrtho, tPoissonOrtho, element[0].mattype());

    std::cout << "target De:" << endl << tCH << endl;
  }
  else
  {
    cerr << "Error in meta_material.cpp" << endl;
    exit(1);
  }

  while(convergence == 1)
  {
    optim.optstep++;

    std::cout << "============== optimization step ==============" << endl;
    std::cout << "step no. " << optim.optstep << endl;

    // filtering procedure
    auto filtered_s =
        filter.densityFilter(optim.design_s, "Eigen_CG", 0, false);
    double beta = min(pow(2, (double)((optim.optstep) / 40)), 8.0);
    std::cout << "beta" << beta << endl;

    auto threshold_s = filter.threshold(filtered_s, beta, T);
    auto gamma = filter.threshold(threshold_s, beta1, rho0);
    auto threDiff = filter.thresholdDerivative(filtered_s, beta, T);

    // set threshold design variable
    for(int i = 0; i < fem.nelx; i++)
    {
      threshold_s[i] = std::max(std::min(threshold_s[i], xmax), xmin); // clamp
      element[i].setDesignS(threshold_s[i]);
      gamma[i] = std::max(std::min(gamma[i], xmax), xmin); // clamp
      element[i].setGamma(gamma[i]);
    }

    // numerical material testing
    mCompt.initialization();

    // set macro strain
    VectorXd Mstrain(fem.voigt), one(fem.voigt);
    MatrixXd CH = tCH;
    CH.setZero();
    Mstrain.setZero();
    one.setOnes();

    // Setting incremental strain
    Mstrain = (1.0 / (double)numstep) * MstrainFac * one;
    // Incremental Analysis Loop
    for(int loadstep = 1; loadstep < numstep + 1; loadstep++)
    {
      cout << "loadstep " << to_string(loadstep) << endl;
      CH = mCompt.solveNLNMT(Mstrain, 1.0e-8, 100, "Eigen_CG", 0, true);
    }
    CH *= 1.0 / MstrainFac;
    std::cout << "CH: " << CH << std::endl;

    for(int i = 0; i < fem.nelx; i++)
      element[i].makeMicroValues(fem, node);

    optim.object_f = getObjectF(fem, tCH, CH, omega);
    std::cout << "function value = " << optim.object_f << endl;

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
    std::cout << "Volume error: " << optim.const_h[0] << endl;

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);
    output(optim.optstep, fem, element, node, dis);

    convergence = optim.isConvergence(1.0e-5, 200);

    sensitivity(fem, element, node, dis, optim, tCH, CH, omega, V0, Vtotal,
                MstrainFac);

    // threshold's derivative
    for(int i = 0; i < fem.nelx; i++)
    {
      optim.dfds[i] *= threDiff[i];
      optim.dgds[i] *= threDiff[i];
    }
    // filter's derivative
    optim.dfds = filter.densityFilter(optim.dfds, "Eigen_CG", 0, true);
    optim.dgds = filter.densityFilter(optim.dgds, "Eigen_CG", 0, true);

    // output CH to ch.csv
    if(optim.optstep == 1)
    {
      ofstream tfile;
      tfile.open("tCH.csv", std::ios::out);
      tfile << tCH << std::endl;
      tfile.close();
    }
    else
    {
      ofstream file;
      file.open("CH.csv", std::ios::out);
      file << CH << std::endl;
      file.close();
    }

    // output sensivity graph
    if(optim.optstep == 1)
      exportGraph(optim.dfds, "objective.csv");

    if(fdm_flag)
      fdm(config, fem, element, node, dis, optim, filter, omega, pairs, length,
          tCH, Mstrain, beta, MstrainFac);

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
  std::cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void output(int optstep, FEM &fem, vector<mHomoNeoHookeTotal> &element,
            vector<HomoNode> &node, ValueHomo &value)
{
  ofstream fout("result_vtk/res" + to_string(optstep) + ".vtk");
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
                      [&](mHomoNeoHookeTotal &e) { return e.design_s(); });

  addElementScalarVTK("young", element, fout, false,
                      [&](mHomoNeoHookeTotal &e) { return e.young(); });

  fout.close();
}
