///  @file  bmatrix.cpp
///  @author  Daiki Watanabe
///  @date  June 11, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "nmatrix_surface.hpp"
#include "gauss_integration.hpp"
#include "shape_function.hpp"
#include <Eigen/LU>
#include <iostream>
#include <problem/base.hpp>
namespace icarat
{

void SurfaceNmatrix::make(Eigen::MatrixXd &Ne, double &jac, std::string eType,
                          int ipmax, Eigen::MatrixXd &X, int ip)
{
  Ne.setZero();

  if(eType == "line2")
    snmatrix1d<ShapeFunction2Line, Gauss2Line>(Ne, jac, X, ip);
  else if(eType == "line3")
    snmatrix1d<ShapeFunction3Line, Gauss3Line>(Ne, jac, X, ip);
  else if(eType == "tria3")
    snmatrix2d<ShapeFunction3Triangle, Gauss1Triangle>(Ne, jac, X, ip);
  else if(eType == "tria6")
    snmatrix2d<ShapeFunction6Triangle, Gauss3Triangle>(Ne, jac, X, ip);
  else if(eType == "quad4")
    snmatrix2d<ShapeFunction4Square, Gauss4Square>(Ne, jac, X, ip);
  else if(eType == "quad8")
    snmatrix2d<ShapeFunction8Square, Gauss9Square>(Ne, jac, X, ip);
  else
  {
    std::cerr << "Element type or ipmax is invalid in nmatrix_surface.cpp. "
              << std::endl;
    exit(1);
  }
}

} // namespace icarat
