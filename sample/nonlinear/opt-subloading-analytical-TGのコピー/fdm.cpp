///  @file    fdm.cpp
///  @author	Daiki Watanabe
///  @date		July 17, 2021.
/// You should run fdm using debug mode.
#include "subloading-analytical.hpp"
using namespace std;
using namespace Eigen;

namespace icarat
{
namespace opt_elast_plast
{
void fdm(FEM &fem, std::vector<OptLinearElastPlast> &element,
         std::vector<Node> &node, Force &force, ValueNR &dis, Optimize &optim,
         int numstep, double filterR)
{
  double delta = 1.0e-7;
  vector<double> object(fem.nelx);
  vector<double> error_of_object(fem.nelx);
  double func_old = optim.object_f;
  auto design_s_old = optim.design_s;

  cout << "------------------------------------------------" << endl;
  cout << "        Sensitivity analysis with FDM           " << endl;
  cout << "------------------------------------------------" << endl;

  NewtonRaphson<OptLinearElastPlast, Node> compt(fem, element, node, force,
                                                 dis);
  Filter<OptLinearElastPlast, Node> filter(fem, element, node, filterR);

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    // set delta
    optim.design_s[nel] += delta;

    auto design_filtered = filter.densityFilter(optim.design_s, false);
    // double beta = 1.0, T = 0.5;
    // //auto design_threshold = filter.threshold(design_filtered, beta, T);
    // auto design_threshold = design_filtered;

    for(int i = 0; i < fem.nelx; i++)
      element[i].setDesignS(design_filtered[i]);

    // cal f(s + delta)
    compt.initialization();

    /// Initialization of element internal variables
    for(int i = 0; i < fem.nelx; i++)
    {
      element[i].initialization();
      element[i].init_opt();
    }

    for(int step = 1; step < numstep + 1; step++)
    {
      force.fext = step * force.fext_org / numstep;

      compt.solve_M2(step, numstep, 1.0e-6, 10, "Eigen_LDLT", 0);

      vector<Triplet<double>> tripletsR;
      for(int i = 0; i < fem.nelx; i++)
      {
        element[i].stress_plastic(fem, node, tripletsR);
        element[i].sensitivity(fem, node);
      }
    }
    /// dfds
    double func_new = 0.0;
    for(int i = 0; i < fem.nelx; i++)
      func_new += element[i].func_new();

    object[nel] = (func_new - func_old) / delta;

    error_of_object[nel] = 100 * (optim.dfds[nel] - object[nel]) / object[nel];

    // reset design variables
    optim.design_s = design_s_old;

    std::cout << "Element ID : " << nel << " dfds macro = " << object[nel]
              << endl;
  }

  // output sensivity graph
  exportGraph(object, "objective_fdm.csv");

  exportGraph(error_of_object, "error_obj.csv");

  std::cout << "------------------------------------------------" << endl;
  std::cout << "         Finish computation with FDM            " << endl;
  std::cout << "------------------------------------------------" << endl;
  exit(0);
}
} // namespace opt_elast_plast
} // namespace icarat
