#include "meta_material.hpp"

namespace icarat
{
namespace meta_material
{
using namespace std;
using namespace Eigen;

double getObjectF(FEM &fem, MatrixXd &tCH, MatrixXd &CH, MatrixXd &omega)
{
  return 0.5 * (omega.array() * (CH - tCH).array().pow(2.0)).sum();
}

void sensitivity(FEM &fem, vector<mHomoElastic> &element,
                 vector<HomoNode> &node, ValueHomo &value, Optimize &optim,
                 MatrixXd &tCH, MatrixXd &CH, MatrixXd &omega, double V0,
                 double VTotal)
{
  assert(element[0].young2() > element[0].young1());

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    double dyoung = element[nel].pp() *
                    pow(element[nel].design_s(), element[nel].pp() - 1.0) *
                    (element[nel].young2() - element[nel].young1());

    MatrixXd dCHds = dyoung / element[nel].young() * element[nel].sEnergy();

    optim.dfds[nel] =
        1.0 / VTotal *
        (omega.array() * (CH - tCH).array() * dCHds.array()).sum();

    optim.dgds[nel] = element[nel].volume / V0;
  }
}

} // namespace meta_material
} // namespace icarat