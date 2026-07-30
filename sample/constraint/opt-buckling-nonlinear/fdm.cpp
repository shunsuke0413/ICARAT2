///  @file    fdm.cpp
///  @author	Daiki Watanabe
///  @date		April 24, 2022.

#include "opt_buckling_nonlinear.hpp"
using namespace std;
using namespace Eigen;

namespace icarat
{
namespace opt_buckling_nonlinear
{
void fdm(FEMBuckling &fem, vector<BucklingFiniteSaintVenantTotal> &element,
         vector<Node> &node, Force &force, ValueEigen &dis, Optimize &optim,
         FilterHelmholtz<BucklingFiniteSaintVenantTotal, Node> &filter,
         double tolNR, int iterNR, int numstep)
{

  cout << "          Sensitivity analysis by FDM          " << endl;

  // init
  double delta = 1.0e-07;
  vector<double> object(fem.nelx);
  double func_old = func_buckling(fem, dis);
  auto design_s_old = optim.design_s;

  // analysis class
  NewtonRaphson<BucklingFiniteSaintVenantTotal, Node> newton(fem, element, node,
                                                             force, dis);
  int paraconv = 3 * fem.numeigen;
  GeneralEigen<BucklingFiniteSaintVenantTotal, Node> eigen(
      fem, element, node, dis, fem.numeigen, paraconv, fem.sigma);

  /// cal df / ds
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    optim.design_s[nel] += delta;

    auto design_filtered =
        filter.densityFilter(optim.design_s, "Eigen_LDLT", 0, false);

    // set threshold function
    double beta = 1.0;
    auto design_threshold = filter.threshold(design_filtered, beta, 0.5);

    // set threshold design variable
    for(int i = 0; i < fem.nelx; i++)
      element[i].setDesignS(design_threshold[i]);

    // structural analysis
    newton.initialization();

    for(int loadStep = 1; loadStep < numstep + 1; loadStep++)
    {
      force.fext = loadStep * force.fext_org / numstep;

      cout << "---------- load-step No. " << loadStep << " ----------" << endl;

      newton.solve(loadStep, numstep, tolNR, iterNR, "Eigen_LDLT", 0);
    }
    // eigenvalue analysis
    eigen.solve();

    double tmp = dis.evectors[0].transpose() * eigen.A() * dis.evectors[0];
    cout << "phiKLphi global: " << tmp << endl;

    double func_new = func_buckling(fem, dis);

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
} // namespace opt_buckling_nonlinear
} // namespace icarat
