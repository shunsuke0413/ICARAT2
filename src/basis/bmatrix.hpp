///  @file  bmatrix.hpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
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
class Bmatrix
{
public:
  /// making N & B matrix
  void make(Eigen::MatrixXd &Be, Eigen::MatrixXd &Ne, double &jac,
            std::string eType, int ipmax, Eigen::MatrixXd &X, int ip);
};

/// make B matrix for 2D
template <class S, class G>
void bmatrix2d(Eigen::MatrixXd &Be, Eigen::MatrixXd &Ne, double &jac,
               Eigen::MatrixXd &X, int ip)
{
  S shape;
  G gauss;

  Eigen::MatrixXd dXdr = shape.dNdr(gauss.points(ip)) * X;
  Eigen::MatrixXd dNdX = dXdr.inverse() * shape.dNdr(gauss.points(ip));
  Eigen::VectorXd NeArr = shape.N(gauss.points(ip));

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
    std::cerr << "Invalid Ne row in bmatrix.hpp" << std::endl;
    exit(1);
  }

  // for heat Bmatrix
  if(Be.rows() == 2)
    for(int n = 0; n < shape.n; n++)
    {
      Be.coeffRef(0, n) = dNdX.coeff(0, n);
      Be.coeffRef(1, n) = dNdX.coeff(1, n);
    }
  // for structure Bmatrix
  else if(Be.rows() == 3)
    for(int n = 0; n < shape.n; n++)
    {
      Be.coeffRef(0, shape.d * n) = dNdX.coeff(0, n);
      Be.coeffRef(1, shape.d * n + 1) = dNdX.coeff(1, n);
      Be.coeffRef(2, shape.d * n) = dNdX.coeff(1, n);
      Be.coeffRef(2, shape.d * n + 1) = dNdX.coeff(0, n);
    }
  // for geometoric Bmatrix
  else if(Be.rows() == 4)
    for(int n = 0; n < shape.n; n++)
    {
      Be.coeffRef(0, shape.d * n) = dNdX.coeff(0, n);
      Be.coeffRef(1, shape.d * n) = dNdX.coeff(1, n);
      Be.coeffRef(2, shape.d * n + 1) = dNdX.coeff(0, n);
      Be.coeffRef(3, shape.d * n + 1) = dNdX.coeff(1, n);
    }
  else
  {
    std::cerr << "Invalid Be row in bmatrix.hpp" << std::endl;
    exit(1);
  }

  // jacobian * gauss weight
  jac = dXdr.determinant() * gauss.weights(ip);
}

/// make B matrix for 3D
template <class S, class G>
void bmatrix3d(Eigen::MatrixXd &Be, Eigen::MatrixXd &Ne, double &jac,
               Eigen::MatrixXd &X, int ip)
{
  S shape;
  G gauss;

  Eigen::MatrixXd dXdr = shape.dNdr(gauss.points(ip)) * X;
  Eigen::MatrixXd dNdX = dXdr.inverse() * shape.dNdr(gauss.points(ip));
  Eigen::VectorXd NeArr = shape.N(gauss.points(ip));

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
    std::cerr << "Invalid Ne row in bmatrix.hpp" << std::endl;
    exit(1);
  }
  // for heat Bmatrix
  if(Be.rows() == 3)
    for(int n = 0; n < shape.n; n++)
    {
      Be.coeffRef(0, n) = dNdX.coeff(0, n);
      Be.coeffRef(1, n) = dNdX.coeff(1, n);
      Be.coeffRef(2, n) = dNdX.coeff(2, n);
    }
  // for structure Bmatrix
  else if(Be.rows() == 6)
    for(int n = 0; n < shape.n; n++)
    {
      Be.coeffRef(0, shape.d * n) = dNdX.coeff(0, n);
      Be.coeffRef(1, shape.d * n + 1) = dNdX.coeff(1, n);
      Be.coeffRef(2, shape.d * n + 2) = dNdX.coeff(2, n);

      Be.coeffRef(3, shape.d * n) = dNdX.coeff(1, n);
      Be.coeffRef(3, shape.d * n + 1) = dNdX.coeff(0, n);

      Be.coeffRef(4, shape.d * n + 1) = dNdX.coeff(2, n);
      Be.coeffRef(4, shape.d * n + 2) = dNdX.coeff(1, n);

      Be.coeffRef(5, shape.d * n) = dNdX.coeff(2, n);
      Be.coeffRef(5, shape.d * n + 2) = dNdX.coeff(0, n);
    }
  // for geometric Bmatrix
  else if(Be.rows() == 9)
    for(int n = 0; n < shape.n; n++)
    {
      Be.coeffRef(0, shape.d * n) = dNdX.coeff(0, n);
      Be.coeffRef(1, shape.d * n) = dNdX.coeff(1, n);
      Be.coeffRef(2, shape.d * n) = dNdX.coeff(2, n);
      Be.coeffRef(3, shape.d * n + 1) = dNdX.coeff(0, n);
      Be.coeffRef(4, shape.d * n + 1) = dNdX.coeff(1, n);
      Be.coeffRef(5, shape.d * n + 1) = dNdX.coeff(2, n);
      Be.coeffRef(6, shape.d * n + 2) = dNdX.coeff(0, n);
      Be.coeffRef(7, shape.d * n + 2) = dNdX.coeff(1, n);
      Be.coeffRef(8, shape.d * n + 2) = dNdX.coeff(2, n);
    }
  else
  {
    std::cerr << "Invalid Be row in bmatrix.hpp" << std::endl;
    exit(1);
  }
  // jacobian * gauss weight
  jac = dXdr.determinant() * gauss.weights(ip);
}

} // namespace icarat
