#include "opt_multiscale.hpp"

using namespace std;
using namespace Eigen;

namespace icarat
{
namespace opt_multiscale
{
void fdm(Macro &M, Micro &m, Force &Mforce, Value &Mdis, ValueHomo &mdis,
         Optimize &optim, Filter<mHomoElastic, HomoNode> &filter,
         vector<vector<pair<int, int>>> &pairs, std::array<double, 3> &length)
{
  // init
  double delta = 1.0e-5;
  vector<double> object(m.fem.nelx);
  auto design_s_old = optim.design_s;
  double func_old = optim.object_f;

  Homogenization<mHomoElastic, HomoNode> compt(m.fem, m.element, m.node, mdis,
                                               pairs, length);
  LinearAnalysis<MHomoElastic, Node> MCompt(M.fem, M.element, M.node, Mforce,
                                            Mdis);

  cout << "------------------------------------------------" << endl;
  cout << "        Sensitivity analysis with FDM           " << endl;
  cout << "------------------------------------------------" << endl;

  for(int nel = 0; nel < m.fem.nelx; nel++)
  {
    // set delta
    optim.design_s[nel] += delta;

    // set threshold design variable
    auto design_filtered = filter.densityFilter(optim.design_s, false);
    double beta = 1.0;
    auto design_threshold = filter.threshold(design_filtered, beta, 0.5);
    for(int i = 0; i < m.fem.nelx; i++)
      m.element[i].setDesignS(design_threshold[i]);

    // cal f(s + delta)
    MatrixXd CH = compt.solveNMT(false);

    // set CH to macro element
    for(int i = 0; i < M.fem.nelx; i++)
      M.element[i].setDe(CH);

    MCompt.solve("Eigen_LDLT", 0);
    for(int i = 0; i < m.fem.nelx; i++)
      m.element[i].makeMicroValues(m.fem, m.node);

    double func_new = Mdis.val.dot(Mforce.fext_org);

    // dfds
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
} // namespace opt_multiscale
} // namespace icarat