///  @file  bmatrix_updated.cpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "bmatrix_updated.hpp"
#include "gauss_integration.hpp"
#include "shape_function.hpp"
#include <Eigen/LU>
#include <iostream>

namespace icarat
{
void BmatrixUpdated::make(Eigen::MatrixXd &BL, Eigen::MatrixXd &BNL,
                          Eigen::MatrixXd &Ne, double &jac, Eigen::MatrixXd &F,
                          std::string &eType, int ipmax, Eigen::MatrixXd &X,
                          Eigen::MatrixXd &U, int ip)
{
  // dummy variables
  Eigen::MatrixXd du = X;
  Eigen::MatrixXd L(X.cols(), X.cols());

  make(BL, BNL, Ne, jac, F, L, du, eType, ipmax, X, U, ip);
}

void BmatrixUpdated::make(Eigen::MatrixXd &BL, Eigen::MatrixXd &BNL,
                          Eigen::MatrixXd &Ne, double &jac, Eigen::MatrixXd &F,
                          Eigen::MatrixXd &L, Eigen::MatrixXd &du,
                          std::string &eType, int ipmax, Eigen::MatrixXd &X,
                          Eigen::MatrixXd &U, int ip)
{
  // init
  BL.setZero();
  BNL.setZero();
  Ne.setZero();
  L.setZero();
  jac = 0.0;
  F.setZero();

  if(eType == "tria3")
    bmatrix2dUpdated<ShapeFunction3Triangle, Gauss1Triangle>(
        BL, BNL, L, Ne, jac, F, X, U, du, ip);
  else if(eType == "tria6")
    bmatrix2dUpdated<ShapeFunction6Triangle, Gauss3Triangle>(
        BL, BNL, L, Ne, jac, F, X, U, du, ip);
  else if(eType == "quad4")
  {
    if(ipmax == 1)
      bmatrix2dUpdated<ShapeFunction4Square, Gauss1Square>(BL, BNL, L, Ne, jac,
                                                           F, X, U, du, ip);
    else if(ipmax == 4)
      bmatrix2dUpdated<ShapeFunction4Square, Gauss4Square>(BL, BNL, L, Ne, jac,
                                                           F, X, U, du, ip);
    else
    {
      std::cerr << "Ipmax is invalid in bmatrix. " << std::endl;
      exit(1);
    }
  }
  else if(eType == "quad8")
    bmatrix2dUpdated<ShapeFunction8Square, Gauss9Square>(BL, BNL, L, Ne, jac, F,
                                                         X, U, du, ip);
  else if(eType == "tetra4")
    bmatrix3dUpdated<ShapeFunction4Tetrahedron, Gauss1Tetrahedron>(
        BL, BNL, L, Ne, jac, F, X, U, du, ip);
  else if(eType == "hexa8")
    bmatrix3dUpdated<ShapeFunction8Cubic, Gauss8Cubic>(BL, BNL, L, Ne, jac, F,
                                                       X, U, du, ip);
  else if(eType == "hexa20")
    bmatrix3dUpdated<ShapeFunction20Cubic, Gauss27Cubic>(BL, BNL, L, Ne, jac, F,
                                                         X, U, du, ip);
  else
  {
    std::cerr << "The combination of ne and ndim is invalid in bmatrix. "
              << std::endl;
    exit(1);
  }
}
} // namespace icarat
