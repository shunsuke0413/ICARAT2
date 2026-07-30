///  @file    fdm.cpp
///  @author	Daiki Watanabe
///  @date		May 6, 2022.
#include "opt_multi_material.hpp"
using namespace std;
using namespace Eigen;

namespace icarat
{
namespace opt_multi_material
{
void fdm(FEM &fem, vector<LinearElastMultiMat> &element, vector<Node> &node,
         Force &force, Value &dis, Optimize &optim,
         Filter<LinearElastMultiMat, Node> &filter, int numMat)
{
  double delta = 1.0e-7;
  vector<double> object(fem.nelx);
  double func_old = optim.object_f;

  ij ss_old(numMat);
  for(int i = 0; i < numMat; i++)
  {
    ss_old[i].resize(fem.nelx);
    for(int j = 0; j < fem.nelx; j++)
      ss_old[i][j] = optim.design_s[i * fem.nelx + j];
  }

  auto ss = ss_old;

  LinearAnalysis<LinearElastMultiMat, Node> compt(fem, element, node, force,
                                                  dis);

  cout << "------------------------------------------------" << endl;
  cout << "        Sensitivity analysis with FDM           " << endl;
  cout << "------------------------------------------------" << endl;

  for(int matID = 0; matID < numMat; matID++)
  {
    for(int nel = 0; nel < fem.nelx; nel++)
    {
      // set delta
      ss[matID][nel] += delta;

      // set threshold design variable
      for(int matID2 = 0; matID2 < numMat; matID2++)
      {
        auto s_filtered = filter.densityFilter(ss[matID2], false);
        double beta = 1.0;
        auto design_threshold = filter.threshold(s_filtered, beta, 0.5);

        for(int i = 0; i < fem.nelx; i++)
          element[i].setDesignSS(matID2, design_threshold[i]);
      }

      // cal f(s + delta)
      compt.solve("Eigen_CG", 0);

      double func_new = dis.val.dot(force.fext);

      /// dfds
      object[nel] = (func_new - func_old) / delta;

      // reset design variables
      ss = ss_old;

      cout << "Element ID : " << nel << " dfds macro = " << object[nel] << endl;
    }

    // output sensivity graph
    exportGraph(object, "objective_fdm" + to_string(matID) + ".csv");
  }

  cout << "------------------------------------------------" << endl;
  cout << "         Finish computation with FDM            " << endl;
  cout << "------------------------------------------------" << endl;
  exit(0);
}
} // namespace opt_multi_material
} // namespace icarat