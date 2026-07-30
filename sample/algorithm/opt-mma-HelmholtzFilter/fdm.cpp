///  @file    fdm.cpp
///  @author	Daiki Watanabe
///  @date		July 17, 2021.
#include "opt_mma_Helmholtz.hpp"
using namespace std;
using namespace Eigen;

namespace icarat
{
namespace helmholtz
{
void fdm(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
         Force &force, Value &dis, Optimize &optim,
         FilterHelmholtz<LinearElastic, Node> &filter)
{
  double delta = 1.0e-7;
  vector<double> object(fem.nelx);
  double func_old = optim.object_f;
  auto design_s_old = optim.design_s;

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
    LinearAnalysis<LinearElastic, Node> compt(fem, element, node, force, dis);
    compt.solve("Eigen_CG", 0);

    double func_new = dis.val.dot(force.fext);

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
} // namespace helmholtz
} // namespace icarat