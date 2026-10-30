///@sa T.Buhl et al., Stiffness design of geometrically nonlinear
/// structures using topology optimization, SMO, 2000.
///@sa F.Wang et al., Interpolation scheme for fictitious domain techniques and
/// topology optimization of finite strain elastic problems, CMAME, 2014.
#include "opt_finite_strain.hpp"
#include <filesystem>
using namespace icarat;
using namespace icarat::opt_finite_strain;
using namespace std;
using namespace Eigen;

void output(FEM &fem, vector<FiniteNeoHookedamage> &element,
            std::vector<Node> &node, Optimize &optim);

void outputReanalysis(FEM &fem, vector<FiniteNeoHookedamage> &element,
                      vector<Node> &node, int loadStep);

// main関数のループ前に宣言
std::vector<double> s_dx, s_fx;

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <.toml file>\n";
    return 1;
  }

  std::filesystem::create_directories("result"); // フォルダの作成

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
  double tolNR = 1.0e-4;
  int iterNR = 100;

  // set input file
  const auto config = toml::parse(argv[1]);

  int numstep = toml::find<int>(config, "this_problem", "numstep");

  // set problem
  FEM fem;
  MeshIcarat<FiniteNeoHookedamage, Node> mesh(fem, config);
  mesh.setParameter("structure");

  // set control node for load-displacement
  int direction;
  int control_ID;
  int reaction_axis;
  int reaction_direction;

  std::string dire =
      toml::find<std::string>(config, "load_displacement", "direction");

  std::array<double, 3> coordinate = toml::find<std::array<double, 3>>(
      config, "load_displacement", "observation_point");

  // 反力を集計する面の設定
  std::string reaction_axis_str =
      toml::find<std::string>(config, "load_displacement", "reaction_axis");

  double reaction_coordinate =
      toml::find<double>(config, "load_displacement", "reaction_coordinate");

  std::string reaction_direction_str = toml::find<std::string>(
      config, "load_displacement", "reaction_direction");

  // 変位方向を整数に変換
  if(dire == "X")
    direction = 0;
  else if(dire == "Y")
    direction = 1;
  else if(dire == "Z")
    direction = 2;
  else
    throw std::invalid_argument("Invalid load_displacement direction");

  // 反力を集計する座標軸を整数に変換
  if(reaction_axis_str == "X")
    reaction_axis = 0;
  else if(reaction_axis_str == "Y")
    reaction_axis = 1;
  else if(reaction_axis_str == "Z")
    reaction_axis = 2;
  else
    throw std::invalid_argument("Invalid reaction_axis");

  // 反力の方向を整数に変換
  if(reaction_direction_str == "X")
    reaction_direction = 0;
  else if(reaction_direction_str == "Y")
    reaction_direction = 1;
  else if(reaction_direction_str == "Z")
    reaction_direction = 2;
  else
    throw std::invalid_argument("Invalid reaction_direction");

  // set element
  vector<FiniteNeoHookedamage> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      FiniteNeoHookedamage actele(fem.voigt, mesh.ne(type), mesh.ipmax(type),
                                  numdof, numstep);
      actele.setParameter(config);
      element.push_back(actele);
    }
  }

  // set node
  vector<Node> node(fem.numnp, Node(1));

  // input mesh to element & node
  mesh.generate(element, node);

  // boundary condition class
  BCIcarat<FiniteNeoHookedamage, Node> BC(fem, element, node);

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq, 1);
  ValueNR dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  for(int i = 0; i < fem.numnp; i++)
  {
    if(node[i].x[0] == coordinate[0] && node[i].x[1] == coordinate[1] &&
       node[i].x[2] == coordinate[2])
    {
      control_ID = i;
      break;
    }
  }

  // set optimization's classes
  Optimize optim(fem.nelx, numConst, xmin, xmax);
  MMASolver mma(fem.nelx, numConst);
  double filterR = toml::find<double>(config, "this_problem", "filterR");
  double poisson = toml::find<double>(config, "material", "poisson");
  Filter<FiniteNeoHookedamage, Node> filter(fem, element, node, filterR);

  // set tau
  std::vector<double> tauFiltered;

  optim.optstep = 0;
  /// filtering procedure
  for(int i = 0; i < fem.nelx; i++)
    optim.design_s[i] = element[i].design_s();

  NewtonRaphson<FiniteNeoHookedamage, Node> compt(fem, element, node, force,
                                                  dis, numstep, tauFiltered);

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

#include <fstream>

    std::ofstream fout("load_displacement_damagefilter.csv");

    fout << "displacement force\n";
    fout.close();

    addPlot(0.0, 0.0, "load_displacement_damagefilter.csv", 0);

    std::vector<double> sum_f(numstep, 0);

    for(int loadStep = 1; loadStep < numstep + 1; loadStep++)
    {
      cout << "---------- load-step No. " << loadStep << " ----------" << endl;

      // compute tau
      auto element_tau_aves = compt.tau();

      auto result = filter.tauFilter(element_tau_aves);

      std::vector<double> tauFiltered = result.first;
      std::vector<double> dtauFiltered_dtau = result.second;

      compt.solve(loadStep, numstep, tolNR, iterNR, "Eigen_LDLT", 0,
                  tauFiltered, dtauFiltered_dtau);

      outputReanalysis(fem, element, node, loadStep);

      double totalReaction = 0.0;

      for(int i = 0; i < node.size(); i++)
      {
        if(std::abs(node[i].x[reaction_axis] - reaction_coordinate) < 1.0e-8)
        {
          totalReaction += node[i].reaction[reaction_direction];
        }
      }

      double displacement = node[control_ID].val[direction];

      addPlot(displacement, totalReaction, "load_displacement_damagefilter.csv",
              0);
    }

    // fout.close();

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

    optim.object_f = 0.0;
    optim.object_f = -compt.se().sum();

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);
    output(fem, element, node, optim);

    // 有限要素法のみ
    exit(1);

    convergence = optim.isConvergence(1.0e-5, 200);

    // analytical
    auto fil_gamma = filter.thresholdDerivative(design_threshold, beta1, rho0);
    sensitivity(fem, element, node, dis, compt.K(), optim, V0, fil_gamma,
                poisson);
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

void output(FEM &fem, vector<FiniteNeoHookedamage> &element, vector<Node> &node,
            Optimize &optim)
{

  // --- デバッグ追加 ---
  std::cout << "DEBUG: Exporting step " << optim.optstep
            << " Node[0] pos = " << node[0].x << endl;

  ofstream fout("res" + to_string(optim.optstep) + ".vtk");
  makeHeaderVTK(fout);

  // --- 修正箇所: ここで座標を明示的に書き出す ---
  fout << "\nPOINTS\t" << node.size() << "\tfloat\n";
  for(const auto &n : node)
  {
    // X, Y は必ず出力
    fout << n.x[0] << "\t" << n.x[1] << "\t";

    // Z成分は次元が3のときだけ計算値を出力し、それ以外は0.0にする
    if(fem.ndim == 3)
    {
      fout << n.x[2];
    }
    else
    {
      fout << "0.0";
    }
    fout << "\t" << std::endl;
  }
  // ------------------------------------------

  // for reaction stress
  /*
  fout << "VECTORS reaction float" << std::endl;//反力計算用一旦省略
  for(int i = 0; i < node.size(); i++) {
      // もし node[i] が固定境界であれば、その節点の反力を出力
      // 固定でなければ 0.0 0.0 0.0 を出力
      if(node[i].is_fixed()) { // このような判定関数がある場合
          fout << reaction[i] << " " << ... << endl;
      } else {
          fout << "0.0 0.0 0.0" << endl;
      }
  }
  */

  // make grid
  // setPointVTK(node, fout);
  setElementVTK(element, fout);
  setElementTypeVTK(fem.ndim, element, fout);

  // export vector data to node
  addPointVectorVTK("displacement", node, fout, true,
                    [&](Node &n) { return n.val; });

  addPointVectorVTK("traction", node, fout, false,
                    [&](Node &n) { return n.du; });

  // export scalar data to element
  addElementScalarVTK("density", element, fout, true,
                      [&](FiniteNeoHookedamage &e) { return e.design_s(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](FiniteNeoHookedamage &e) { return e.mises(); });

  addElementScalarVTK("damage", element, fout, false,
                      [&](FiniteNeoHookedamage &e)
                      { return e.getAverageDamage(); });

  fout.close();
}

void outputReanalysis(FEM &fem, vector<FiniteNeoHookedamage> &element,
                      vector<Node> &node, int loadStep)
{
  std::ostringstream filename;

  filename << "result/result_" << std::setw(4) << std::setfill('0') << loadStep
           << ".vtk";

  std::ofstream fout(filename.str());

  makeHeaderVTK(fout);

  // =========================
  // points
  // =========================
  fout << "\nPOINTS\t" << node.size() << "\tfloat\n";

  for(const auto &n : node)
  {
    fout << n.x[0] << "\t" << n.x[1] << "\t";

    if(fem.ndim == 3)
    {
      fout << n.x[2];
    }
    else
    {
      fout << "0.0";
    }

    fout << "\t" << std::endl;
  }

  setElementVTK(element, fout);
  setElementTypeVTK(fem.ndim, element, fout);

  addPointVectorVTK("displacement", node, fout, true,
                    [&](Node &n) { return n.val; });

  addElementScalarVTK("density", element, fout, true,
                      [&](FiniteNeoHookedamage &e) { return e.design_s(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](FiniteNeoHookedamage &e) { return e.mises(); });

  addElementScalarVTK("damage", element, fout, false,
                      [&](FiniteNeoHookedamage &e)
                      { return e.getAverageDamage(); });

  fout.close();
}