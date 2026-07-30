///  @file    fdm.cpp
///  @author	Daiki Watanabe
///  @date		December 10, 2021.

#include "meta_material.hpp"
using namespace std;
using namespace Eigen;

namespace icarat
{
namespace meta_material
{
void fdm(const toml::value &config, FEM &fem, vector<mHomoElastic> &element,
         vector<HomoNode> &node, ValueHomo &dis, Optimize &optim,
         FilterHelmholtz<mHomoElastic, HomoNode> &filter, MatrixXd &omega,
         vector<vector<pair<int, int>>> &pairs, array<double, 3> &length)
{
  // init
  double delta = 1.0e-5;
  vector<double> object(fem.nelx);
  auto design_s_old = optim.design_s;
  double tYoung = toml::find<double>(config, "this_problem", "target_young");
  double tPoisson =
      toml::find<double>(config, "this_problem", "target_poisson");
  MatrixXd tCH = MatrixXd::Zero(fem.voigt, fem.voigt);
  element[0].De(tCH, tYoung, tPoisson);
  Homogenization<mHomoElastic, HomoNode> compt(fem, element, node, dis, pairs,
                                               length);
  MatrixXd CH = compt.solveNMT(true);
  double func_old = optim.object_f;

  cout << "------------------------------------------------" << endl;
  cout << "        Sensitivity analysis with FDM           " << endl;
  cout << "------------------------------------------------" << endl;

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    // set delta
    optim.design_s[nel] += delta;

    // set threshold design variable
    auto design_filtered =
        filter.densityFilter(optim.design_s, "Eigen_CG", 0, false);
    double beta = 1.0;
    auto design_threshold = filter.threshold(design_filtered, beta, 0.5);

    for(int i = 0; i < fem.nelx; i++)
      element[i].setDesignS(design_threshold[i]);

    // cal f(s + delta)
    CH = compt.solveNMT(false);

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
} // namespace meta_material
} // namespace icarat