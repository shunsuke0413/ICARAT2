#include "opt_differentiate_then_discretize.hpp"
using namespace std;
using namespace Eigen;

namespace icarat {
namespace opt_differentiate_then_discretize {

template <class Obj>
void fdm(FEMDyna &fem, vector<DynamicElastic> &element, vector<Node> &node,
         Force &force, ValueDyna &dis, Optimize &optim,
         Filter<DynamicElastic, Node> &filter, std::vector<double> &dfds,
         history &history, Obj &function) {
  double delta = 1.0e-7;
  vector<double> object(fem.nelx);
  double func_old = optim.object_f;
  auto design_s_old = optim.design_s;
  ImplicitDyna<DynamicElastic, Node> compt(fem, element, node, force, dis);

  cout << "------------------------------------------------" << endl;
  cout << "        Sensitivity analysis with FDM           " << endl;
  cout << "------------------------------------------------" << endl;

  for (int nel = 0; nel < fem.nelx; nel++) {
    // set delta
    optim.design_s[nel] += delta;

    // set threshold design variable
    auto design_filtered = filter.densityFilter(optim.design_s, false);

    // set value
    for (int i = 0; i < fem.nelx; i++)
      element[i].setDesignS(optim.design_s[i]);

    // cal f(s + delta)
    compt.initialization();
    function.init();
    for (int time = 1; time <= fem.timestep; time++) {
      force.fext = sin(2.0 * M_PI * time / fem.timestep) * force.fext_org;
      compt.calRayleigh(time, fem.timestep, fem.aa, fem.bb, fem.deltaT,
                        fem.beta, fem.gamma, fem.solvName, 0);

      history.disp_h[time] = dis.val;
      history.velo_h[time] = dis.velocity;
      history.acce_h[time] = dis.accel;
      function.update(time);
    }
    double func_new = function.getValue();
    /// dfds
    object[nel] = (func_new - func_old) / delta;

    // reset design variables
    optim.design_s = design_s_old;

    cout << "Element ID : " << nel << " dfds macro = " << object[nel] << endl;
  }

  object = filter.densityFilter(object, true);
  std::ofstream ofs("sensitivity.dat");
  for (size_t i = 0; i < object.size(); i++) {
    ofs << i << " " << std::scientific << std::setprecision(14) << dfds[i]
        << " " << object[i] << std::endl;
  }

  // output sensivity graph
  exportGraph(object, "objective_fdm.csv");

  cout << "------------------------------------------------" << endl;
  cout << "         Finish computation with FDM            " << endl;
  cout << "------------------------------------------------" << endl;
  exit(0);
}
} // namespace opt_differentiate_then_discretize
} // namespace icarat
