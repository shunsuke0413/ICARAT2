///  @file    fdm.cpp
///  @author	Daiki Watanabe
///  @date		May 5, 2022.
#include "opt_finite_strain.hpp"
using namespace std;
using namespace Eigen;

namespace icarat
{
namespace opt_finite_se_damage_filter
{
void fdm(FEM &fem, vector<FiniteNeoHookedamage> &element, vector<Node> &node,
         Force &force, ValueNR &dis, Optimize &optim,
         Filter<FiniteNeoHookedamage, Node> &filter, int numstep,
         std::vector<double> tauFiltered, std::vector<double> dtauFiltered_dtau)
{
  double xmin = 1.0e-4;
  double xmax = 1.0;
  double delta = 1.0e-7;
  double beta1 = 50.0;
  double rho0 = 0.1;
  vector<double> object(fem.nelx);
  double func_old = optim.object_f;
  auto design_s_old = optim.design_s;
  NewtonRaphson<FiniteNeoHookedamage, Node> compt(fem, element, node, force,
                                                  dis, numstep, tauFiltered);

  cout << "------------------------------------------------" << endl;
  cout << "        Sensitivity analysis with FDM           " << endl;
  cout << "------------------------------------------------" << endl;

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    // set delta
    optim.design_s[nel] += delta;

    // set threshold design variable
    auto design_filtered = filter.densityFilter(optim.design_s, false);
    double beta = 1.0;
    auto design_threshold = filter.threshold(design_filtered, beta, 0.5);
    auto gamma = filter.threshold(design_threshold, beta1, rho0);

    int i;
    std::vector<Eigen::Triplet<double>> tripletC;

    for(int i = 0; i < fem.nelx; i++)
    {
      design_threshold[i] =
          std::max(std::min(design_threshold[i], xmax), xmin); // clamp
      element[i].setDesignS(design_threshold[i]);
      gamma[i] = std::max(std::min(gamma[i], xmax), xmin); // clamp
      element[i].setGamma(gamma[i]);
    }

    // cal f(s + delta)
    for(int loadStep = 1; loadStep < numstep + 1; loadStep++)
    {
      force.fext = loadStep * force.fext_org / numstep;

      cout << "---------- load-step No. " << loadStep << " ----------" << endl;

      compt.solve(loadStep, numstep, 1.0e-8, 10, "Eigen_CG", 0, tauFiltered,
                  dtauFiltered_dtau);
    }

    double func_new = -compt.se().sum();

    /// dfds
    object[nel] = (func_new - func_old) / delta;

    // reset design variables
    optim.design_s = design_s_old;

    cout << "Element ID : " << nel << " dfds macro = " << object[nel] << endl;
  }

  // output sensivity graph
  exportGraph(object, "objective_fdm.csv");

  cout << "------------------------------------------------" << endl;
  cout << "         Finish computation with FDM            " << endl;
  cout << "------------------------------------------------" << endl;
  exit(0);
}
} // namespace opt_finite_se_damage_filter
} // namespace icarat