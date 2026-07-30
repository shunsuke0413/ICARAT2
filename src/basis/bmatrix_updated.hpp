///  @file  bmatrix_updated.hpp
///  @author  Daiki Watanabe
///  @date  June 3, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php
#pragma once
#include "gauss_integration.hpp"
#include "shape_function.hpp"
#include <Eigen/Core>
#include <Eigen/LU>
#include <iostream>
#include <problem/base.hpp>
#include <vector>

namespace icarat
{
class BmatrixUpdated
{
public:
  /// main function of making B matrix and deformation gradient tensor F
  /// BL...constitutive low's term 'Bonet book (9.35)). Size is (voigt,numdof)
  /// BNL...initial stress term (Bonet book (9.44c)). Size is (ndim*ndim,numdof)
  void make(Eigen::MatrixXd &BL, Eigen::MatrixXd &BNL, Eigen::MatrixXd &Ne,
            double &jac, Eigen::MatrixXd &F, std::string &eType, int ipmax,
            Eigen::MatrixXd &X, Eigen::MatrixXd &U, int ip);

  /// Overload when using velocity gradient
  void make(Eigen::MatrixXd &BL, Eigen::MatrixXd &BNL, Eigen::MatrixXd &Ne,
            double &jac, Eigen::MatrixXd &F, Eigen::MatrixXd &L,
            Eigen::MatrixXd &du, std::string &eType, int ipmax,
            Eigen::MatrixXd &X, Eigen::MatrixXd &U, int ip);
};

template <class S, class G>
void bmatrix2dUpdated(Eigen::MatrixXd &BL, Eigen::MatrixXd &BNL,
                      Eigen::MatrixXd &LL, Eigen::MatrixXd &Ne, double &jac,
                      Eigen::MatrixXd &F, Eigen::MatrixXd &X,
                      Eigen::MatrixXd &U, Eigen::MatrixXd &du, int ip)
{
  S shape;
  G gauss;

  Eigen::MatrixXd x = X + U;
  Eigen::VectorXd NeArr = shape.N(gauss.points(ip));
  Eigen::MatrixXd dxdr = shape.dNdr(gauss.points(ip)) * x; ///<(d,d)
  Eigen::MatrixXd dNdx =
      dxdr.inverse() * shape.dNdr(gauss.points(ip)); ///<(d,n)

  ///@ F=(dN(X)/dr)^{-1} (dN(x)/dr)
  F = ((shape.dNdr(gauss.points(ip)) * X).inverse() *
       shape.dNdr(gauss.points(ip)) * x)
          .transpose();

  LL = dNdx * du; ///<(d,d)

  if(Ne.rows() == 1)
    for(int n = 0; n < shape.n; n++)
      Ne.coeffRef(n) = NeArr.coeff(n);
  else if(Ne.rows() == 2)
    for(int n = 0; n < shape.n; n++)
    {
      Ne.coeffRef(0, shape.d * n) = NeArr.coeff(n);
      Ne.coeffRef(1, shape.d * n + 1) = NeArr.coeff(n);
    }
  else
  {
    std::cerr << "Invalid row in bmatrix_updated.hpp" << std::endl;
    exit(1);
  }

  for(int n = 0; n < shape.n; n++)
  {
    BL.coeffRef(0, shape.d * n) = dNdx.coeff(0, n);
    BL.coeffRef(1, shape.d * n + 1) = dNdx.coeff(1, n);
    BL.coeffRef(2, shape.d * n) = dNdx.coeff(1, n);
    BL.coeffRef(2, shape.d * n + 1) = dNdx.coeff(0, n);

    BNL.coeffRef(0, shape.d * n) = dNdx.coeff(0, n);
    BNL.coeffRef(1, shape.d * n + 1) = dNdx.coeff(1, n);
    BNL.coeffRef(2, shape.d * n) = dNdx.coeff(0, n);
    BNL.coeffRef(3, shape.d * n + 1) = dNdx.coeff(1, n);
  }

  // jacobian * gauss weight
  jac = dxdr.determinant() * gauss.weights(ip);
}

template <class S, class G>
void bmatrix3dUpdated(Eigen::MatrixXd &BL, Eigen::MatrixXd &BNL,
                      Eigen::MatrixXd &LL, Eigen::MatrixXd &Ne, double &jac,
                      Eigen::MatrixXd &F, Eigen::MatrixXd &X,
                      Eigen::MatrixXd &U, Eigen::MatrixXd &du, int ip)
{
  S shape;
  G gauss;

  Eigen::MatrixXd x = X + U;
  Eigen::VectorXd NeArr = shape.N(gauss.points(ip));
  Eigen::MatrixXd dxdr = shape.dNdr(gauss.points(ip)) * x; ///<(d,d)
  Eigen::MatrixXd dNdx =
      dxdr.inverse() * shape.dNdr(gauss.points(ip)); ///<(d,n)

  ///@ F=(dN(X)/dr)^{-1} (dN(x)/dr)
  F = ((shape.dNdr(gauss.points(ip)) * X).inverse() *
       shape.dNdr(gauss.points(ip)) * x)
          .transpose();

  LL = dNdx * du; ///<(d,d)

  if(Ne.rows() == 1)
    for(int n = 0; n < shape.n; n++)
      Ne.coeffRef(n) = NeArr.coeff(n);
  else if(Ne.rows() == 3)
    for(int n = 0; n < shape.n; n++)
    {
      Ne.coeffRef(0, shape.d * n) = NeArr.coeff(n);
      Ne.coeffRef(1, shape.d * n + 1) = NeArr.coeff(n);
      Ne.coeffRef(2, shape.d * n + 2) = NeArr.coeff(n);
    }
  else
  {
    std::cerr << "Invalid row in bmatrix.hpp" << std::endl;
    exit(1);
  }

  for(int n = 0; n < shape.n; n++)
  {
    // xx,yy,zz
    BL.coeffRef(0, 3 * n) = dNdx.coeff(0, n);
    BL.coeffRef(1, 3 * n + 1) = dNdx.coeff(1, n);
    BL.coeffRef(2, 3 * n + 2) = dNdx.coeff(2, n);

    // xy
    BL.coeffRef(3, 3 * n) = dNdx.coeff(1, n);
    BL.coeffRef(3, 3 * n + 1) = dNdx.coeff(0, n);

    // yz
    BL.coeffRef(4, 3 * n + 1) = dNdx.coeff(2, n);
    BL.coeffRef(4, 3 * n + 2) = dNdx.coeff(1, n);

    // xz
    BL.coeffRef(5, 3 * n) = dNdx.coeff(2, n);
    BL.coeffRef(5, 3 * n + 2) = dNdx.coeff(0, n);

    // x?
    BNL.coeffRef(0, 3 * n) = dNdx.coeff(0, n);
    BNL.coeffRef(1, 3 * n + 1) = dNdx.coeff(1, n);
    BNL.coeffRef(2, 3 * n + 2) = dNdx.coeff(2, n);

    // y?
    BNL.coeffRef(3, 3 * n) = dNdx.coeff(0, n);
    BNL.coeffRef(4, 3 * n + 1) = dNdx.coeff(1, n);
    BNL.coeffRef(5, 3 * n + 2) = dNdx.coeff(2, n);

    // z?
    BNL.coeffRef(6, 3 * n) = dNdx.coeff(0, n);
    BNL.coeffRef(7, 3 * n + 1) = dNdx.coeff(1, n);
    BNL.coeffRef(8, 3 * n + 2) = dNdx.coeff(2, n);
  }

  // jacobian * gauss weight
  jac = dxdr.determinant() * gauss.weights(ip);
}

} // namespace icarat