///  @file  bmatrix_updated.cpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under tthickness MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "bmatrix_shell.hpp"
#include "gauss_integration.hpp"
#include "shape_function.hpp"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <iostream>
#include <misc/eigen_manipulation.hpp>

namespace icarat
{
void BmatrixShell::make(Eigen::MatrixXd &Be, Eigen::MatrixXd &Ne,
                        Eigen::MatrixXd &P, double &jac, int ndim, int ne,
                        int ipmax, Eigen::MatrixXd &X, Eigen::Vector3d &length,
                        double thickness, int ip, int ipz)
{
  Be.setZero();
  Ne.setZero();

  if(ne == 4 && ipmax == 4)
    square4gauss4(Be, Ne, P, jac, X, length, thickness, ip, ipz);
  else
  {
    std::cerr
        << "Tthickness combination of ne and ndim is invalid in bmatrix shell. "
        << std::endl;
    exit(1);
  }
}

inline void BmatrixShell::square4gauss4(Eigen::MatrixXd &Be,
                                        Eigen::MatrixXd &Ne, Eigen::MatrixXd &P,
                                        double &jac, Eigen::MatrixXd &X,
                                        Eigen::Vector3d &length,
                                        double thickness, int ip, int ipz)
{
  // Gauss1Square gauss;
  Gauss4Square gauss;
  Gauss2Line gaussLine;
  ShapeFunction4Square shape;

  Eigen::VectorXd NeArr = shape.N(gauss.points(ip));

  Eigen::MatrixXd v1 = Eigen::MatrixXd::Zero(shape.n, 3);
  Eigen::MatrixXd v2 = Eigen::MatrixXd::Zero(shape.n, 3);
  Eigen::MatrixXd v3 = Eigen::MatrixXd::Zero(shape.n, 3);

  for(int nn = 0; nn < shape.n; nn++)
  {
    Eigen::MatrixXd dXdr = shape.dNdr(gauss.points(ip)) * X; ///<(2,3)
    Eigen::Vector3d vecTmp1 = (dXdr.transpose()).block(0, 0, 3, 1);
    Eigen::Vector3d vecTmp2 = (dXdr.transpose()).block(0, 1, 3, 1);

    Eigen::Vector3d v3i = (vecTmp1.cross(vecTmp2)).normalized();
    for(int i = 0; i < 3; i++)
      v3(nn, i) = v3i(i);

    Eigen::Vector3d v1i = (length.cross(v3i)).normalized();
    for(int i = 0; i < 3; i++)
      v1(nn, i) = v1i(i);

    Eigen::Vector3d v2i = (v3i.cross(v1i)).normalized();
    for(int i = 0; i < 3; i++)
      v2(nn, i) = v2i(i);
  }

  Eigen::VectorXd zeta = gaussLine.points(ipz);
  Eigen::MatrixXd hx1 = shape.dNdr(gauss.points(ip)) *
                        (X + 0.5 * thickness * zeta(0) * v3);     ///<(2,3)
  Eigen::MatrixXd hx2 = 0.5 * thickness * NeArr.transpose() * v3; ///<(1,3)

  Eigen::MatrixXd J = vstack(hx1, hx2); ///<(3,3)
  Eigen::MatrixXd invJ = J.inverse();

  Eigen::MatrixXd dNdx = invJ.block(0, 0, 3, 2) * shape.dNdr(gauss.points(ip));
  Eigen::MatrixXd dNzdx =
      dNdx * zeta(0) + invJ.block(0, 2, 3, 1) * NeArr.transpose();

  Eigen::Vector3d J0 = J.block(0, 0, 1, 3).transpose();
  Eigen::Vector3d J1 = J.block(1, 0, 1, 3).transpose();
  Eigen::Vector3d P2 = (J0.cross(J1)).normalized();
  Eigen::Vector3d P0 = (length.cross(P2)).normalized();
  Eigen::Vector3d P1 = (P2.cross(P0)).normalized();

  for(int i = 0; i < shape.n; i++)
  {
    if(Ne.rows() == 1)
      for(int n = 0; n < shape.n; n++)
        Ne.coeffRef(n) = NeArr.coeff(n);
    else if(Ne.rows() == 2)
      for(int n = 0; n < shape.n; n++)
      {
        Ne.coeffRef(0, shape.d * n) = NeArr.coeff(n);
        Ne.coeffRef(1, shape.d * n + 1) = NeArr.coeff(n);
      }

    // rotation tensor
    // xx,yy,zz,xy,yz,xz
    P(0, 0) = P0(0) * P0(0);
    P(0, 1) = P0(1) * P0(1);
    P(0, 2) = P0(2) * P0(2);
    P(0, 3) = P0(0) * P0(1);
    P(0, 4) = P0(1) * P0(2);
    P(0, 5) = P0(0) * P0(2);

    P(1, 0) = P1(0) * P1(0);
    P(1, 1) = P1(1) * P1(1);
    P(1, 2) = P1(2) * P1(2);
    P(1, 3) = P1(0) * P1(1);
    P(1, 4) = P1(1) * P1(2);
    P(1, 5) = P1(0) * P1(2);

    P(2, 0) = 2.0 * P0(0) * P1(0);
    P(2, 1) = 2.0 * P0(1) * P1(1);
    P(2, 2) = 2.0 * P0(2) * P1(2);
    P(2, 3) = P0(0) * P1(1) + P0(1) * P1(0);
    P(2, 4) = P0(1) * P1(2) + P0(2) * P1(1);
    P(2, 5) = P0(0) * P1(2) + P0(2) * P1(0);

    P(3, 0) = 2.0 * P0(0) * P2(0);
    P(3, 1) = 2.0 * P0(1) * P2(1);
    P(3, 2) = 2.0 * P0(2) * P2(2);
    P(3, 3) = P0(0) * P2(1) + P0(1) * P2(0);
    P(3, 4) = P0(1) * P2(2) + P0(2) * P2(1);
    P(3, 5) = P0(0) * P2(2) + P0(2) * P2(0);

    P(4, 0) = 2.0 * P1(0) * P2(0);
    P(4, 1) = 2.0 * P1(1) * P2(1);
    P(4, 2) = 2.0 * P1(2) * P2(2);
    P(4, 3) = P1(0) * P2(1) + P1(1) * P2(0);
    P(4, 4) = P1(1) * P2(2) + P1(2) * P2(1);
    P(4, 5) = P1(0) * P2(2) + P1(2) * P2(0);

    /// xx
    Be(0, 6 * i) = dNdx(0, i);
    Be(0, 6 * i + 1) = 0.0;
    Be(0, 6 * i + 2) = 0.0;
    Be(0, 6 * i + 3) = -0.5 * thickness * dNzdx(0, i) * v2(i, 0) * v1(i, 0) +
                       0.5 * thickness * dNzdx(0, i) * v1(i, 0) * v2(i, 0);
    Be(0, 6 * i + 4) = -0.5 * thickness * dNzdx(0, i) * v2(i, 0) * v1(i, 1) +
                       0.5 * thickness * dNzdx(0, i) * v1(i, 0) * v2(i, 1);
    Be(0, 6 * i + 5) = -0.5 * thickness * dNzdx(0, i) * v2(i, 0) * v1(i, 2) +
                       0.5 * thickness * dNzdx(0, i) * v1(i, 0) * v2(i, 2);
    /// yy
    Be(1, 6 * i) = 0.0;
    Be(1, 6 * i + 1) = dNdx(1, i);
    Be(1, 6 * i + 2) = 0.0;
    Be(1, 6 * i + 3) = -0.5 * thickness * dNzdx(1, i) * v2(i, 1) * v1(i, 0) +
                       0.5 * thickness * dNzdx(1, i) * v1(i, 1) * v2(i, 0);
    Be(1, 6 * i + 4) = -0.5 * thickness * dNzdx(1, i) * v2(i, 1) * v1(i, 1) +
                       0.5 * thickness * dNzdx(1, i) * v1(i, 1) * v2(i, 1);
    Be(1, 6 * i + 5) = -0.5 * thickness * dNzdx(1, i) * v2(i, 1) * v1(i, 2) +
                       0.5 * thickness * dNzdx(1, i) * v1(i, 1) * v2(i, 2);
    /// zz
    Be(2, 6 * i) = 0.0;
    Be(2, 6 * i + 1) = 0.0;
    Be(2, 6 * i + 2) = dNdx(2, i);
    Be(2, 6 * i + 3) = -0.5 * thickness * dNzdx(2, i) * v2(i, 2) * v1(i, 0) +
                       0.5 * thickness * dNzdx(2, i) * v1(i, 2) * v2(i, 0);
    Be(2, 6 * i + 4) = -0.5 * thickness * dNzdx(2, i) * v2(i, 2) * v1(i, 1) +
                       0.5 * thickness * dNzdx(2, i) * v1(i, 2) * v2(i, 1);
    Be(2, 6 * i + 5) = -0.5 * thickness * dNzdx(2, i) * v2(i, 2) * v1(i, 2) +
                       0.5 * thickness * dNzdx(2, i) * v1(i, 2) * v2(i, 2);
    /// xy
    Be(3, 6 * i) = dNdx(1, i);
    Be(3, 6 * i + 1) = dNdx(0, i);
    Be(3, 6 * i + 2) = 0.0;
    Be(3, 6 * i + 3) =
        -0.5 * thickness * (dNzdx(1, i) * v2(i, 0) + dNzdx(0, i) * v2(i, 1)) *
            v1(i, 0) +
        0.5 * thickness * (dNzdx(1, i) * v1(i, 0) + dNzdx(0, i) * v1(i, 1)) *
            v2(i, 0);
    Be(3, 6 * i + 4) =
        -0.5 * thickness * (dNzdx(1, i) * v2(i, 0) + dNzdx(0, i) * v2(i, 1)) *
            v1(i, 1) +
        0.5 * thickness * (dNzdx(1, i) * v1(i, 0) + dNzdx(0, i) * v1(i, 1)) *
            v2(i, 1);
    Be(3, 6 * i + 5) =
        -0.5 * thickness * (dNzdx(1, i) * v2(i, 0) + dNzdx(0, i) * v2(i, 1)) *
            v1(i, 2) +
        0.5 * thickness * (dNzdx(1, i) * v1(i, 0) + dNzdx(0, i) * v1(i, 1)) *
            v2(i, 2);
    /// yz
    Be(4, 6 * i) = 0.0;
    Be(4, 6 * i + 1) = dNdx(2, i);
    Be(4, 6 * i + 2) = dNdx(1, i);
    Be(4, 6 * i + 3) =
        -0.5 * thickness * (dNzdx(2, i) * v2(i, 1) + dNzdx(1, i) * v2(i, 2)) *
            v1(i, 0) +
        0.5 * thickness * (dNzdx(2, i) * v1(i, 1) + dNzdx(1, i) * v1(i, 2)) *
            v2(i, 0);
    Be(4, 6 * i + 4) =
        -0.5 * thickness * (dNzdx(2, i) * v2(i, 1) + dNzdx(1, i) * v2(i, 2)) *
            v1(i, 1) +
        0.5 * thickness * (dNzdx(2, i) * v1(i, 1) + dNzdx(1, i) * v1(i, 2)) *
            v2(i, 1);
    Be(4, 6 * i + 5) =
        -0.5 * thickness * (dNzdx(2, i) * v2(i, 1) + dNzdx(1, i) * v2(i, 2)) *
            v1(i, 2) +
        0.5 * thickness * (dNzdx(2, i) * v1(i, 1) + dNzdx(1, i) * v1(i, 2)) *
            v2(i, 2);
    /// xz
    Be(5, 6 * i) = dNdx(2, i);
    Be(5, 6 * i + 1) = 0.0;
    Be(5, 6 * i + 2) = dNdx(0, i);
    Be(5, 6 * i + 3) =
        -0.5 * thickness * (dNzdx(0, i) * v2(i, 2) + dNzdx(2, i) * v2(i, 0)) *
            v1(i, 0) +
        0.5 * thickness * (dNzdx(0, i) * v1(i, 2) + dNzdx(2, i) * v1(i, 0)) *
            v2(i, 0);
    Be(5, 6 * i + 4) =
        -0.5 * thickness * (dNzdx(0, i) * v2(i, 2) + dNzdx(2, i) * v2(i, 0)) *
            v1(i, 1) +
        0.5 * thickness * (dNzdx(0, i) * v1(i, 2) + dNzdx(2, i) * v1(i, 0)) *
            v2(i, 1);
    Be(5, 6 * i + 5) =
        -0.5 * thickness * (dNzdx(0, i) * v2(i, 2) + dNzdx(2, i) * v2(i, 0)) *
            v1(i, 2) +
        0.5 * thickness * (dNzdx(0, i) * v1(i, 2) + dNzdx(2, i) * v1(i, 0)) *
            v2(i, 2);
  }

  jac = J.determinant();
  if(jac <= 0.0)
  {
    std::cerr << " jacobian <= 0.0 in bmatrix_shell.cpp" << std::endl;
    exit(1);
  }
  jac *= gauss.weights(ip) * gaussLine.weights(ipz);
}

} // namespace icarat
