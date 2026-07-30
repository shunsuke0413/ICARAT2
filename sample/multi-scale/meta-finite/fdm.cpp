///  @file    fdm.cpp
///  @author	Daiki Watanabe
///  @date		December 10, 2021.

#include "meta_material_finite.hpp"
using namespace std;
using namespace Eigen;

namespace icarat
{
namespace meta_material_finite
{
void fdm(const toml::value &config, FEM &fem,
         vector<mHomoNeoHookeTotal> &element, vector<HomoNode> &node,
         ValueHomo &dis, Optimize &optim,
         FilterHelmholtz<mHomoNeoHookeTotal, HomoNode> &filter,
         Eigen::MatrixXd &omega, vector<vector<pair<int, int>>> &pairs,
         array<double, 3> &length, Eigen::MatrixXd &tCH,
         Eigen::VectorXd &Mstrain, double beta, double fac)
{
  // init
  double delta = 1.0e-5;
  double xmin = 1.0e-4, xmax = 1.0;
  double beta1 = 50.0, rho0 = 0.1;
  vector<double> object(fem.nelx);
  auto design_s_old = optim.design_s;

  Homogenization<mHomoNeoHookeTotal, HomoNode> compt(fem, element, node, dis,
                                                     pairs, length);
  double func_old = optim.object_f;

  cout << "------------------------------------------------" << endl;
  cout << "        Sensitivity analysis with FDM           " << endl;
  cout << "------------------------------------------------" << endl;

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    // set delta
    optim.design_s[nel] += delta;

    // set threshold design variable
    auto filtered_s =
        filter.densityFilter(optim.design_s, "Eigen_CG", 0, false);

    auto threshold_s = filter.threshold(filtered_s, beta, 0.5);
    auto gamma = filter.threshold(threshold_s, beta1, rho0);

    for(int i = 0; i < fem.nelx; i++)
    {
      threshold_s[i] = std::max(std::min(threshold_s[i], xmax), xmin); // clamp
      element[i].setDesignS(threshold_s[i]);
      gamma[i] = std::max(std::min(gamma[i], xmax), xmin); // clamp
      element[i].setGamma(gamma[i]);
    }
    // cal f(s + delta)
    compt.initialization();
    MatrixXd CH = compt.solveNLNMT(Mstrain, 1.0e-8, 100, "Eigen_CG", 0, true);
    CH *= 1.0 / fac;

    double func_new = getObjectF(fem, tCH, CH, omega);

    // dfds
    object[nel] = (func_new - func_old) / delta;

    // reset design variables
    optim.design_s = design_s_old;

    cout << "Element ID : " << nel << " dfds macro = " << object[nel] << endl;
  }

  optim.dfds = object;
  // output sensivity graph
  exportGraph(object, "objective_fdm.csv");

  cout << "------------------------------------------------" << endl;
  cout << "         Finish computation with FDM            " << endl;
  cout << "------------------------------------------------" << endl;
  exit(0);
}
} // namespace meta_material_finite
} // namespace icarat