#include "opt_discretize_then_differentiate.hpp"
using namespace std;
using namespace Eigen;

namespace icarat {
namespace opt_discretize_then_differentiate {

template <class Obj>
void fdm(FEMDyna &fem, vector<DynamicElastic> &element, vector<Node> &node,
         Force &force, ValueDyna &dis, Optimize &optim, Obj &function,
         Filter<DynamicElastic, Node> &filter) {
  double delta = 1.0e-7;
  vector<double> object(fem.nelx);
  double func_old = optim.object_f;
  auto design_s_old = optim.design_s;
  ImplicitDyna<DynamicElastic, Node> compt(fem, element, node, force, dis);

  cout << "------------------------------------------------" << endl;
  cout << "        Sensitivity analysis with FDM           " << endl;
  cout << "------------------------------------------------" << endl;

  for (int nel = 0; nel < fem.nelx; nel++) {
    // init objective function
    function.init();

    // set delta
    optim.design_s[nel] += delta;

    // set threshold design variable
    auto design_filtered = filter.densityFilter(optim.design_s, false);
    double beta = 1.0;
    auto design_threshold = filter.threshold(design_filtered, beta, 0.5);

    for (int i = 0; i < fem.nelx; i++)
      element[i].setDesignS(design_threshold[i]);

    // cal f(s + delta)
    compt.initialization();
    for (int time = 1; time < fem.timestep + 1; time++) {
      force.fext = time * force.fext_org / fem.timestep;
      compt.calRayleigh(time, fem.timestep, fem.aa, fem.bb, fem.deltaT,
                        fem.beta, fem.gamma, "Eigen_CG", 0);

      function.update(dis.val);
    }
    double func_new = function.get();

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
} // namespace opt_discretize_then_differentiate
} // namespace icarat
