///  @file    fdm.cpp
///  @author	Daiki Watanabe
///  @date		May 27, 2023.

#include "opt_stress_nonlinear.hpp"
using namespace std;
using namespace Eigen;

namespace icarat
{
namespace stress_nonlinear
{
void fdm(FEMStressConst &fem, vector<FiniteSaintVenantTotal> &element,
         vector<Node> &node, Force &force, Value &dis, Optimize &optim,
         Filter<FiniteSaintVenantTotal, Node> &filter, int numstep,
         double tolNR, int iterNR)
{

  cout << "          Sensitivity analysis by FDM          " << endl;

  // init
  double delta = 1.0e-07;
  vector<double> object(fem.nelx);
  double func_old = optim.const_h[1]; // stress constraint function
  auto design_s_old = optim.design_s;

  /// cal df / ds
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    optim.design_s[nel] += delta;

    auto design_filtered = filter.densityFilter(optim.design_s, false);
    double beta = 1.0, T = 0.5;
    auto design_threshold = filter.threshold(design_filtered, beta, T);

    // set threshold design variable
    for(int i = 0; i < fem.nelx; i++)
      element[i].setDesignS(design_threshold[i]);

    /// cal f(s + delta)
    NewtonRaphson<FiniteSaintVenantTotal, Node> compt(fem, element, node, force,
                                                      dis);
    compt.initialization();
    for(int loadStep = 1; loadStep < numstep + 1; loadStep++)
    {
      force.fext = loadStep * force.fext_org / numstep;

      cout << "---------- load-step No. " << loadStep << " ----------" << endl;

      compt.solve(loadStep, numstep, tolNR, iterNR, "Eigen_LDLT", 0);
    }

    /// get  stress
    for(auto &e : element)
      e.makeMisesStress(fem, node);

    double func_new = funcStressConst(fem, element);

    /// dfds
    object[nel] = (func_new - func_old) / delta;

    /// init
    optim.design_s = design_s_old;

    cout << "Element ID : " << nel << " dfds macro = " << object[nel] << endl;
  }

  // output sensivity graph
  exportGraph(object, "objective_fdm.csv");

  cout << "========================================================" << endl;
  cout << "Finish computation succecfully at fdm.cpp               " << endl;
  cout << "========================================================" << endl;
  exit(0);
}
} // namespace stress_nonlinear
} // namespace icarat
