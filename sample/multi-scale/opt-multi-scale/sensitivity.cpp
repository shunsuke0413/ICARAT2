#include "opt_multiscale.hpp"

namespace icarat
{
namespace opt_multiscale
{
using namespace std;
using namespace Eigen;

void sensitivity(Macro &M, Force &Mforce, Value &Mdis, Micro &m,
                 Optimize &optim, double V0)
{
  assert(m.element[0].young2() > m.element[0].young1());

  LinearAnalysis<MHomoElastic, Node> MCompt(M.fem, M.element, M.node, Mforce,
                                            Mdis);

  for(int nel = 0; nel < m.fem.nelx; nel++)
  {
    // make derivative
    double dyoung = m.element[nel].pp() *
                    pow(m.element[nel].design_s(), m.element[nel].pp() - 1.0) *
                    (m.element[nel].young2() - m.element[nel].young1());

    MatrixXd dCHds = dyoung / m.element[nel].young() * m.element[nel].sEnergy();

    // set to macro & assembling
    for(int i = 0; i < M.fem.nelx; i++)
      M.element[i].setDe(dCHds);

    MCompt.assemblingK(0);

    optim.dfds[nel] = -Mdis.val.transpose() * MCompt.K() * Mdis.val;
    optim.dgds[nel] = m.element[nel].volume / V0;
  }
}

} // namespace opt_multiscale
} // namespace icarat