///@sa T.Buhl et al., Stiffness design of geometrically nonlinear
/// structures using topology optimization, SMO, 2000.
///@sa F.Wang et al., Interpolation scheme for fictitious domain techniques and
/// topology optimization of finite strain elastic problems, CMAME, 2014.
#include "opt_finite_strain.hpp"
#include <filesystem>
#include <limits>
#include <optional>
#include <stdexcept>
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
  std::filesystem::create_directories("optimization");

  // measure time
  auto t_start = chrono::system_clock::now();
  int numConst = 1;    // number of constraint
  double xmin = 1.0e-4;
  double xmax = 1.0;
  double V0 = 0.0; // first volume
  double eta = 0.5;
  double beta1 = 50.0;
  double rho0 = 0.1;
  double tolNR = 1.0e-8;
  int iterNR = 100;

  // set input file
  const auto config = toml::parse(argv[1]);

  const int numstep = toml::find<int>(config, "this_problem", "numstep");
  const double beta = toml::find_or<double>(config.at("this_problem"), "beta", 1.0);
  iterNR = toml::find_or<int>(config.at("this_problem"), "nr_max_iter", 100);
  if(numstep <= 0 || iterNR < 1 || !std::isfinite(beta) || !(beta > 0.0))
    throw std::runtime_error("Require positive fixed numstep, beta and nr_max_iter.");

  const int optimizationSteps = toml::find_or<int>(config.at("this_problem"), "optimization_steps", 50);
  const int designRetryLimit = toml::find_or<int>(config.at("this_problem"), "design_retry_limit", 6);
  const double moveShrink = toml::find_or<double>(config.at("this_problem"), "move_shrink", 0.5);
  const double minMoveLimit = toml::find_or<double>(config.at("this_problem"), "min_move_limit", 0.0005);
  if(optimizationSteps < 1 || designRetryLimit < 0 || !std::isfinite(moveShrink) ||
     !(moveShrink > 0.0 && moveShrink < 1.0) || !std::isfinite(minMoveLimit) || !(minMoveLimit > 0.0))
    throw std::runtime_error("Invalid optimization count or design retry settings.");

  // set problem
  FEM fem;
  MeshIcarat<FiniteNeoHookedamage, Node> mesh(fem, config);
  mesh.setParameter("structure");

  // set control node for load-displacement
  int direction;
  int control_ID = -1;
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
  if(fem.ndim == 2)
    for(auto &n : node)
    {
      n.x[2] = 0.0;
      n.val[2] = 0.0;
    }

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
    bool matches = true;
    for(int d = 0; d < fem.ndim; ++d)
      matches = matches && std::abs(node[i].x[d] - coordinate[d]) <= 1.e-8;
    if(matches)
    {
      control_ID = i;
      break;
    }
  }

  if(control_ID < 0 || numstep <= 0)
    throw std::runtime_error("Invalid observation node or load-step count.");

  if(direction >= fem.ndim || reaction_axis >= fem.ndim || reaction_direction >= fem.ndim)
    throw std::runtime_error("Observation and reaction directions must match mesh dimension.");
  const bool analysisOnly = toml::find_or<bool>(config.at("this_problem"), "analysis_only", false);

  // set optimization's classes
  Optimize optim(fem.nelx, numConst, xmin, xmax);
  std::optional<MMASolver> mma;
  mma.emplace(fem.nelx, numConst);
  const double moveLimit =
      toml::find_or<double>(config.at("this_problem"), "move_limit", 0.05);
  if(!(moveLimit > 0.0 && moveLimit <= 1.0))
    throw std::runtime_error("move_limit must be in (0,1].");
  if(minMoveLimit > moveLimit)
    throw std::runtime_error("min_move_limit must not exceed move_limit.");
  double currentMoveLimit = moveLimit;
  mma->SetMoveLimit(currentMoveLimit);
  std::optional<Optimize> acceptedOptim;
  std::optional<MMASolver> mmaBeforeUpdate;
  double filterR = toml::find<double>(config, "this_problem", "filterR");
  const bool checkFdm = toml::find<bool>(config, "this_problem", "fdm");
  std::cout << "Executable: " << argv[0] << '\n'
            << "Configuration: " << std::filesystem::absolute(argv[1]) << '\n'
            << "Mode: "
            << (checkFdm ? "FDM verification (one optimization step)"
                         : "MMA optimization")
            << std::endl;
  if(checkFdm)
    tolNR = 1.0e-10;
  Filter<FiniteNeoHookedamage, Node> filter(fem, element, node, filterR);

  const auto tauWeights = damageFilterWeights(fem, element, node, filterR);
  double previousObjective = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> previousDesign;
  std::ofstream diagnostics("optimization_diagnostics.csv");
  diagnostics << "step,beta,objective,objective_change,design_change,volume_"
                 "constraint,max_tau_input,max_damage,move_limit,design_retries\n";

  // set tau
  std::vector<double> tauFiltered;

  optim.optstep = 0;
  /// filtering procedure
  for(int i = 0; i < fem.nelx; i++)
    optim.design_s[i] = element[i].design_s();

  NewtonRaphson<FiniteNeoHookedamage, Node> compt(fem, element, node, force,
                                                  dis, numstep, tauFiltered);

  // optimization start
  for(int acceptedStep = 1; acceptedStep <= optimizationSteps; ++acceptedStep)
  {
    optim.optstep = acceptedStep;

    std::cout << "============== optimization step ==============" << endl;
    std::cout << "step no. " << optim.optstep << endl;

    std::vector<double> design_filtered, design_threshold, projected, fil_gamma;
    DamageHistory history;
    int designRetries = 0;
    for(;;)
    {
      design_filtered = filter.densityFilter(optim.design_s, false);

      // Use the same projection and displacement discretization at every step.
      std::cout << "Fixed beta=" << beta << ", fixed numstep=" << numstep << std::endl;

      design_threshold = filter.threshold(design_filtered, beta, eta);

      projected = design_threshold;
      for(int i = 0; i < fem.nelx; ++i)
      {
        design_threshold[i] = std::clamp(design_threshold[i], xmin, xmax);
        element[i].setDesignS(design_threshold[i]);
      }
      auto gamma = filter.threshold(design_threshold, beta1, rho0);
      fil_gamma = filter.thresholdDerivative(design_threshold, beta1, rho0);
      for(int i = 0; i < fem.nelx; ++i)
      {
        if(gamma[i] <= xmin || gamma[i] >= xmax)
          fil_gamma[i] = 0.0;
        gamma[i] = std::clamp(gamma[i], xmin, xmax);
        element[i].setGamma(gamma[i]);
      }

      try
      {
        compt.initialization();
        initializeDamagePath(fem, element, node);
        history.clear();
        history.push_back(captureDamageFrame(fem, element, node));

        std::ofstream fout("load_displacement_damagefilter.csv");

        fout << "displacement force\n";
        fout.close();

        addPlot(0.0, 0.0, "load_displacement_damagefilter.csv", 0);

        for(int loadStep = 1; loadStep < numstep + 1; loadStep++)
        {
          cout << "---------- optimization " << optim.optstep
               << ", load step " << loadStep << "/" << numstep
               << " ----------" << endl;

          // compute tau
          auto element_tau_aves = compt.tau();

          auto result = filter.tauFilter(element_tau_aves);

          std::vector<double> tauFiltered = result;

          force.fext = static_cast<double>(loadStep) * force.fext_org / numstep;
          const int iterations = compt.solve(loadStep, numstep, tolNR, iterNR,
                                             "Eigen_LDLT", 0, tauFiltered);
          if(iterations >= iterNR)
            throw std::runtime_error("NR load step " + std::to_string(loadStep) +
                "/" + std::to_string(numstep) + ": reached " +
                std::to_string(iterNR) + " iterations without convergence.");
          std::cout << "NR converged: load step " << loadStep << "/" << numstep
                    << ", iterations=" << iterations << std::endl;
          history.push_back(
              captureDamageFrame(fem, element, node, tauFiltered));
          outputReanalysis(fem, element, node, loadStep);

          double totalReaction = 0.0;

          for(int i = 0; i < node.size(); i++)
          {
            if(std::abs(node[i].x[reaction_axis] - reaction_coordinate) <
               1.0e-8)
            {
              totalReaction += node[i].reaction[reaction_direction];
            }
          }

          double displacement = node[control_ID].val[direction];

          addPlot(displacement, totalReaction,
                  "load_displacement_damagefilter.csv", 0);
        }

        break;
      }
      catch(const std::runtime_error &error)
      {
        const std::string reason = error.what();
        if(reason.find("NR load step ") != 0)
          throw;
        if(!acceptedOptim || !mmaBeforeUpdate)
          throw std::runtime_error("Initial design cannot be solved at the fixed load steps. " + reason);
        const double nextMove = std::max(minMoveLimit, currentMoveLimit * moveShrink);
        if(designRetries >= designRetryLimit || !(nextMove < currentMoveLimit))
        {
          optim = *acceptedOptim;
          mma.emplace(*mmaBeforeUpdate);
          throw std::runtime_error("Design retries exhausted at optimization step " +
              std::to_string(acceptedStep) + ". Last converged design is step " +
              std::to_string(optim.optstep) + ". Last failure: " + reason);
        }
        ++designRetries;
        currentMoveLimit = nextMove;
        // Restore the SAME accepted design, gradients and pre-update MMA
        // asymptote history, then generate a new proposal at the smaller move.
        optim = *acceptedOptim;
        mma.emplace(*mmaBeforeUpdate);
        mma->SetMoveLimit(currentMoveLimit);
        std::cout << "Rollback rejected design at optimization step " << acceptedStep
                  << "; retry " << designRetries << "/" << designRetryLimit
                  << ", move_limit=" << currentMoveLimit << "; reason: " << reason << std::endl;
        mma->Update(optim.design_s.data(), optim.dfds.data(), optim.const_h.data(),
                   optim.dgds.data(), optim.xmin.data(), optim.xmax.data());
        optim.optstep = acceptedStep;
        // The next attempt resets displacement and damage and rebuilds both
        // density filters from this proposal; beta and numstep remain fixed.
      }
    }

    // totalKirchhoff already computes the damaged Cauchy von Mises stress
    // exposed by neohookedamage::mises(). Do not call the linear Elastic
    // postprocessor here; it uses a different stress model and base storage.

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
    optim.object_f = -damageWork(fem, element, node, history);

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);
    output(fem, element, node, optim);

    if(analysisOnly)
    {
      std::cout << "Structural analysis completed (" << fem.ndim << "D)." << std::endl;
      return 0;
    }
    // The residual keeps the original lagged damage update. Its history
    // adjoint propagates nonlocal coupling to the preceding displacement.
    std::cout << "Starting sensitivity analysis for optimization step " << optim.optstep << std::endl;
    sensitivity(fem, element, node, optim, V0, fil_gamma, history, tauWeights);
    std::cout << "Sensitivity analysis completed." << std::endl;
    if(checkFdm)
      verifyDamageDerivatives(fem, element, node, history, tauWeights,
                              fil_gamma);
    auto threDiff = filter.thresholdDerivative(design_filtered, beta, eta);
    for(int i = 0; i < fem.nelx; i++)
    {
      if(projected[i] < xmin || projected[i] > xmax)
        threDiff[i] = 0.0;
      optim.dfds[i] *= threDiff[i];
      optim.dgds[i] *= threDiff[i];
    }
    optim.dfds = filter.densityFilter(optim.dfds, true);
    optim.dgds = filter.densityFilter(optim.dgds, true);

    // output sensivity graph
    if(optim.optstep == 1)
      exportGraph(optim.dfds, "objective.csv");

    // fdm
    if(checkFdm)
    {
      fdm(fem, element, node, force, dis, optim, filter, numstep, tauWeights,
          beta, V0);
      std::cout << "FDM verification completed; optimization was not run."
                << std::endl;
      return 0;
    }

    const double objectiveChange =
        std::isfinite(previousObjective)
            ? std::abs(optim.object_f - previousObjective) /
                  std::max(1.0e-30, std::abs(optim.object_f))
            : std::numeric_limits<double>::infinity();
    double designChange = std::numeric_limits<double>::infinity();
    if(!previousDesign.empty())
    {
      designChange = 0.0;
      for(int i = 0; i < fem.nelx; ++i)
        designChange = std::max(
            designChange, std::abs(optim.design_s[i] - previousDesign[i]));
    }
    double maxTau = 0.0, maxDamage = 0.0;
    for(const auto &frame : history)
    {
      maxTau = std::max(maxTau, frame.tauInput.maxCoeff());
      maxDamage = std::max(maxDamage, frame.damage.maxCoeff());
    }
    diagnostics << optim.optstep << ',' << beta << ',' << optim.object_f << ','
                << objectiveChange << ',' << designChange << ','
                << optim.const_h[0] << ',' << maxTau << ',' << maxDamage << ','
                << currentMoveLimit << ',' << designRetries << std::endl;
    std::cout << "Optimization: df=" << objectiveChange
              << " dx=" << designChange << " volume=" << optim.const_h[0]
              << " max filtered tau=" << maxTau << " max damage=" << maxDamage
              << std::endl;
    // A fixed number of successfully evaluated designs is requested.
    // Stopping at the specified count is not a convergence certificate.
    if(acceptedStep == optimizationSteps)
    {
      std::cout << "Completed requested optimization count: " << optimizationSteps
                << ". No unevaluated final MMA update is applied." << std::endl;
      break;
    }
    previousObjective = optim.object_f;
    previousDesign = optim.design_s;

    // Get updated design variables with MMA
    acceptedOptim = optim;
    mmaBeforeUpdate.emplace(*mma);
    mma->SetMoveLimit(currentMoveLimit);
    std::cout << "Updating design with MMA (move_limit=" << currentMoveLimit << ")." << std::endl;
    mma->Update(optim.design_s.data(), optim.dfds.data(), optim.const_h.data(),
               optim.dgds.data(), optim.xmin.data(), optim.xmax.data());
    std::cout << "Optimization step " << optim.optstep << " completed; continuing to next step."
              << std::endl;

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

  ofstream fout("optimization/res" + to_string(optim.optstep) + ".vtk");
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
