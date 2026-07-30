#include "meta_material_finite.hpp"

namespace icarat
{
namespace meta_material_finite
{
using namespace std;
using namespace Eigen;

Eigen::MatrixXd isoDe(int voigt, double tYoung, double tPoisson, string mattype)
{
  MatrixXd De = MatrixXd::Zero(voigt, voigt);

  // 2D
  if(De.rows() == 3)
  {

    if(mattype == "plane_stress")
    {
      De.coeffRef(0, 0) = 1.0;
      De.coeffRef(0, 1) = tPoisson;
      De.coeffRef(1, 1) = 1.0;
      De.coeffRef(2, 2) = (0.5 * (1.0 - tPoisson));
      De.coeffRef(1, 0) = De(0, 1);
      De.coeffRef(2, 0) = De(0, 2);
      De.coeffRef(2, 1) = De(1, 2);
      De *= tYoung / (1.0 - tPoisson * tPoisson);
    }
    else if(mattype == "plane_strain")
    {
      De.coeffRef(0, 0) = 1.0;
      De.coeffRef(0, 1) = (tPoisson / (1.0 - tPoisson));
      De.coeffRef(1, 1) = 1.0;
      De.coeffRef(2, 2) = ((1.0 - 2.0 * tPoisson) / (2.0 * (1.0 - tPoisson)));
      De.coeffRef(1, 0) = De(0, 1);
      De.coeffRef(2, 0) = De(0, 2);
      De.coeffRef(2, 1) = De(1, 2);
      De *= tYoung *
            ((1.0 - tPoisson) / ((1.0 - 2.0 * tPoisson) * (1.0 + tPoisson)));
    }
    else
    {
      cerr << "Error in elasticity.cpp" << endl;
      exit(1);
    }
  }
  // 3D
  else if(De.rows() == 6)
  {
    // lame constant
    double lambda =
        tYoung * tPoisson / ((1.0 + tPoisson) * (1.0 - 2.0 * tPoisson));
    double mu = tYoung / (2 * (1.0 + tPoisson));

    De.coeffRef(0, 0) = lambda + 2 * mu;
    De.coeffRef(0, 1) = lambda;
    De.coeffRef(0, 2) = lambda;

    De.coeffRef(1, 0) = lambda;
    De.coeffRef(1, 1) = lambda + 2 * mu;
    De.coeffRef(1, 2) = lambda;

    De.coeffRef(2, 0) = lambda;
    De.coeffRef(2, 1) = lambda;
    De.coeffRef(2, 2) = lambda + 2 * mu;

    De.coeffRef(3, 3) = mu;
    De.coeffRef(4, 4) = mu;
    De.coeffRef(5, 5) = mu;
  }
  else
  {
    cerr << "Error in elasticity.cpp" << endl;
    exit(1);
  }

  return De;
}

Eigen::MatrixXd orthoDe(int voigt, array<double, 3> tYoung,
                        array<double, 6> tPoisson, string mattype)
{
  MatrixXd De = MatrixXd::Zero(voigt, voigt);

  // 2D
  if(De.rows() == 3)
  {

    if(mattype == "plane_stress")
    {
      De.coeffRef(0, 0) = tYoung[0] / (1.0 - tPoisson[0] * tPoisson[1]);
      De.coeffRef(0, 1) =
          tPoisson[1] * tYoung[0] / (1.0 - tPoisson[0] * tPoisson[1]);
      De.coeffRef(1, 1) = tYoung[1] / (1.0 - tPoisson[0] * tPoisson[1]);
      De.coeffRef(1, 0) =
          tPoisson[1] * tYoung[0] / (1.0 - tPoisson[0] * tPoisson[1]);
    }
    else if(mattype == "plane_strain")
    {
      double tdelta = (1.0 - tPoisson[0] * tPoisson[1] -
                       tPoisson[2] * tPoisson[3] - tPoisson[4] * tPoisson[5] -
                       2.0 * tPoisson[1] * tPoisson[3] * tPoisson[4]) /
                      (tYoung[0] * tYoung[1] * tYoung[2]);

      De.coeffRef(0, 0) =
          (1.0 - tPoisson[2] * tPoisson[3]) / (tYoung[1] * tYoung[2] * tdelta);
      De.coeffRef(0, 1) = (tPoisson[1] + tPoisson[5] * tPoisson[2]) /
                          (tYoung[1] * tYoung[2] * tdelta);
      De.coeffRef(1, 1) =
          (1.0 - tPoisson[5] * tPoisson[4]) / (tYoung[2] * tYoung[0] * tdelta);
      De.coeffRef(1, 0) = (tPoisson[0] + tPoisson[4] * tPoisson[3]) /
                          (tYoung[2] * tYoung[0] * tdelta);
    }
    else
    {
      cerr << "Error in elasticity.cpp" << endl;
      exit(1);
    }
  }
  // 3D
  else if(De.rows() == 6)
  {
    // lame constant
    double tdelta = (1.0 - tPoisson[0] * tPoisson[1] -
                     tPoisson[2] * tPoisson[3] - tPoisson[4] * tPoisson[5] -
                     2.0 * tPoisson[1] * tPoisson[3] * tPoisson[4]) /
                    (tYoung[0] * tYoung[1] * tYoung[2]);

    De.coeffRef(0, 0) =
        (1.0 - tPoisson[2] * tPoisson[3]) / (tYoung[1] * tYoung[2] * tdelta);
    De.coeffRef(0, 1) = (tPoisson[1] + tPoisson[2] * tPoisson[5]) /
                        (tYoung[1] * tYoung[2] * tdelta);
    De.coeffRef(0, 2) = (tPoisson[5] + tPoisson[1] * tPoisson[3]) /
                        (tYoung[1] * tYoung[2] * tdelta);

    De.coeffRef(1, 0) = De.coeffRef(0, 1);
    De.coeffRef(1, 1) =
        (1.0 - tPoisson[4] * tPoisson[5]) / (tYoung[0] * tYoung[2] * tdelta);
    De.coeffRef(1, 2) = (tPoisson[3] + tPoisson[0] * tPoisson[5]) /
                        (tYoung[0] * tYoung[2] * tdelta);

    De.coeffRef(2, 0) = De.coeffRef(0, 2);
    De.coeffRef(2, 1) = De.coeffRef(1, 2);
    De.coeffRef(2, 2) =
        (1.0 - tPoisson[0] * tPoisson[1]) / (tYoung[0] * tYoung[1] * tdelta);
  }
  else
  {
    cerr << "Error in elasticity.cpp" << endl;
    exit(1);
  }

  return De;
}

double getObjectF(FEM &fem, MatrixXd &tCH, MatrixXd &CH, MatrixXd &omega)
{
  return 0.5 * (omega.array() * (CH - tCH).array().pow(2.0)).sum();
}

void sensitivity(FEM &fem, vector<mHomoNeoHookeTotal> &element,
                 vector<HomoNode> &node, ValueHomo &value, Optimize &optim,
                 MatrixXd &tCH, MatrixXd &CH, MatrixXd &omega, double V0,
                 double VTotal, double fac)
{
  assert(element[0].young2() > element[0].young1());

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    double dyoung = element[nel].pp() *
                    pow(element[nel].design_s(), element[nel].pp() - 1.0) *
                    (element[nel].young2() - element[nel].young1());

    MatrixXd dCHds = dyoung / element[nel].young() * element[nel].sEnergy();

    optim.dfds[nel] =
        1.0 / fac / fac / VTotal *
        (omega.array() * (CH - tCH).array() * dCHds.array()).sum();

    optim.dgds[nel] = element[nel].volume / V0;
  }
}

} // namespace meta_material_finite
} // namespace icarat