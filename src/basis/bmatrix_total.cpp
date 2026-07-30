///  @file  bmatrix_total.cpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "bmatrix_total.hpp"
#include "gauss_integration.hpp"
#include "shape_function.hpp"
#include <Eigen/LU>
#include <iostream>
#include <problem/base.hpp>

namespace icarat
{
void BmatrixTotal::make(Eigen::MatrixXd &BL, Eigen::MatrixXd &BNL,
                        Eigen::MatrixXd &Ne, double &jac, Eigen::MatrixXd &F,
                        std::string eType, int ipmax, Eigen::MatrixXd &X,
                        Eigen::MatrixXd &U, int ip)
{
  // init
  BL.setZero();
  BNL.setZero();
  Ne.setZero();
  jac = 0.0;
  F.setZero();

  if(eType == "tria3")
    bmatrix2dTotal<ShapeFunction3Triangle, Gauss1Triangle>(BL, BNL, Ne, jac, F,
                                                           X, U, ip);
  else if(eType == "tria6")
    bmatrix2dTotal<ShapeFunction6Triangle, Gauss3Triangle>(BL, BNL, Ne, jac, F,
                                                           X, U, ip);
  else if(eType == "quad4")
  {
    if(ipmax == 1)
      bmatrix2dTotal<ShapeFunction4Square, Gauss1Square>(BL, BNL, Ne, jac, F, X,
                                                         U, ip);
    else if(ipmax == 4)
      bmatrix2dTotal<ShapeFunction4Square, Gauss4Square>(BL, BNL, Ne, jac, F, X,
                                                         U, ip);
  }
  else if(eType == "quad8")
    bmatrix2dTotal<ShapeFunction8Square, Gauss9Square>(BL, BNL, Ne, jac, F, X,
                                                       U, ip);
  else if(eType == "tetra4")
    bmatrix3dTotal<ShapeFunction4Tetrahedron, Gauss1Tetrahedron>(
        BL, BNL, Ne, jac, F, X, U, ip);
  else if(eType == "prism6")
    bmatrix3dTotal<ShapeFunction6Prism, Gauss6Prism>(BL, BNL, Ne, jac, F, X, U,
                                                     ip);
  else if(eType == "hexa8")
    bmatrix3dTotal<ShapeFunction8Cubic, Gauss8Cubic>(BL, BNL, Ne, jac, F, X, U,
                                                     ip);
  else if(eType == "tetra10")
    bmatrix3dTotal<ShapeFunction10Tetrahedron, Gauss4Tetrahedron>(
        BL, BNL, Ne, jac, F, X, U, ip);
  else if(eType == "prism15")
    bmatrix3dTotal<ShapeFunction15Prism, Gauss21Prism>(BL, BNL, Ne, jac, F, X,
                                                       U, ip);
  else if(eType == "hexa20")
    bmatrix3dTotal<ShapeFunction20Cubic, Gauss27Cubic>(BL, BNL, Ne, jac, F, X,
                                                       U, ip);
  else
  {
    std::cerr << "eType is invalid in bmatrix. " << std::endl;
    exit(1);
  }
}

void BmatrixTotal::makeDifferential(Eigen::MatrixXd &dBLdu,
                                    Eigen::MatrixXd &dFdu, double &dJdu,
                                    std::string eType, int ipmax,
                                    Eigen::MatrixXd &X, Eigen::MatrixXd &U,
                                    int ip, int dofid)
{
  // init
  dBLdu.setZero();
  dFdu.setZero();

  if(eType == "tria3")
    dBdu2dTotal<ShapeFunction3Triangle, Gauss1Triangle>(dBLdu, dFdu, dJdu, X, U,
                                                        ip, dofid);
  else if(eType == "tria6")
    dBdu2dTotal<ShapeFunction6Triangle, Gauss3Triangle>(dBLdu, dFdu, dJdu, X, U,
                                                        ip, dofid);
  if(eType == "quad4")
  {
    if(ipmax == 1)
      dBdu2dTotal<ShapeFunction4Square, Gauss1Square>(dBLdu, dFdu, dJdu, X, U,
                                                      ip, dofid);
    else if(ipmax == 4)
      dBdu2dTotal<ShapeFunction4Square, Gauss4Square>(dBLdu, dFdu, dJdu, X, U,
                                                      ip, dofid);
  }
  else if(eType == "quad8")
    dBdu2dTotal<ShapeFunction8Square, Gauss9Square>(dBLdu, dFdu, dJdu, X, U, ip,
                                                    dofid);
  else if(eType == "tetra4")
    dBdu3dTotal<ShapeFunction4Tetrahedron, Gauss1Tetrahedron>(dBLdu, dFdu, dJdu,
                                                              X, U, ip, dofid);
  else if(eType == "prism6")
    dBdu3dTotal<ShapeFunction6Prism, Gauss6Prism>(dBLdu, dFdu, dJdu, X, U, ip,
                                                  dofid);
  else if(eType == "hexa8")
    dBdu3dTotal<ShapeFunction8Cubic, Gauss8Cubic>(dBLdu, dFdu, dJdu, X, U, ip,
                                                  dofid);
  else if(eType == "tetra10")
    dBdu3dTotal<ShapeFunction10Tetrahedron, Gauss4Tetrahedron>(
        dBLdu, dFdu, dJdu, X, U, ip, dofid);
  else if(eType == "prism15")
    dBdu3dTotal<ShapeFunction15Prism, Gauss21Prism>(dBLdu, dFdu, dJdu, X, U, ip,
                                                    dofid);
  else if(eType == "hexa20")
    dBdu3dTotal<ShapeFunction20Cubic, Gauss27Cubic>(dBLdu, dFdu, dJdu, X, U, ip,
                                                    dofid);
  else
  {
    std::cerr << "The combination of ne and ndim is invalid in bmatrix. "
              << std::endl;
    exit(1);
  }
}

} // namespace icarat
