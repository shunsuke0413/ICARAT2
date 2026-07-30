///  @file  bmatrix.cpp
///  @author  Daiki Watanabe
///  @date  June 11, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "bmatrix.hpp"
#include "gauss_integration.hpp"
#include "shape_function.hpp"
#include <Eigen/LU>
#include <iostream>
#include <problem/base.hpp>
namespace icarat
{

void Bmatrix::make(Eigen::MatrixXd &Be, Eigen::MatrixXd &Ne, double &jac,
                   std::string eType, int ipmax, Eigen::MatrixXd &X, int ip)
{
  Be.setZero();
  Ne.setZero();

  if(eType == "tria3")
    bmatrix2d<ShapeFunction3Triangle, Gauss1Triangle>(Be, Ne, jac, X, ip);
  else if(eType == "tria6")
    bmatrix2d<ShapeFunction6Triangle, Gauss3Triangle>(Be, Ne, jac, X, ip);
  else if(eType == "quad4")
  {
    if(ipmax == 1)
      bmatrix2d<ShapeFunction4Square, Gauss1Square>(Be, Ne, jac, X, ip);
    else if(ipmax == 4)
      bmatrix2d<ShapeFunction4Square, Gauss4Square>(Be, Ne, jac, X, ip);
    else if(ipmax == 9)
      bmatrix2d<ShapeFunction4Square, Gauss9Square>(Be, Ne, jac, X, ip);
  }
  else if(eType == "quad8")
    bmatrix2d<ShapeFunction8Square, Gauss9Square>(Be, Ne, jac, X, ip);
  else if(eType == "tetra4")
    bmatrix3d<ShapeFunction4Tetrahedron, Gauss1Tetrahedron>(Be, Ne, jac, X, ip);
  else if(eType == "prism6")
    bmatrix3d<ShapeFunction6Prism, Gauss6Prism>(Be, Ne, jac, X, ip);
  else if(eType == "hexa8")
    bmatrix3d<ShapeFunction8Cubic, Gauss8Cubic>(Be, Ne, jac, X, ip);
  else if(eType == "tetra10")
    bmatrix3d<ShapeFunction10Tetrahedron, Gauss4Tetrahedron>(Be, Ne, jac, X,
                                                             ip);
  else if(eType == "prism15")
    bmatrix3d<ShapeFunction15Prism, Gauss21Prism>(Be, Ne, jac, X, ip);
  else if(eType == "hexa20")
    bmatrix3d<ShapeFunction20Cubic, Gauss27Cubic>(Be, Ne, jac, X, ip);
  else
  {
    std::cerr << "Element type or ipmax is invalid in bmatrix.cpp. "
              << std::endl;
    exit(1);
  }
}

} // namespace icarat
