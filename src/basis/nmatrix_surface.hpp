#pragma once
#include "gauss_integration.hpp"
#include "shape_function.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <iostream>
#include <problem/base.hpp>
#include <vector>

namespace icarat
{
class SurfaceNmatrix
{
public:
  /// making N & B matrix
  void make(Eigen::MatrixXd &Ne, double &jac, std::string eType, int ipmax,
            Eigen::MatrixXd &X, int ip);
};

/// make surface shape function (line in 2D)
template <class S, class G>
void snmatrix1d(Eigen::MatrixXd &Ne, double &jac, Eigen::MatrixXd &X, int ip)
{
  S shape;
  G gauss;

  Eigen::MatrixXd dXdr = shape.dNdr(gauss.points(ip)) * X;
  Eigen::VectorXd NeArr = shape.N(gauss.points(ip));

  if(Ne.rows() == 1)
    for(int n = 0; n < shape.n; n++)
      Ne.coeffRef(n) = NeArr.coeff(n);
  else if(Ne.rows() == 2)
    for(int n = 0; n < shape.n; n++)
    {
      Ne.coeffRef(0, 2 * n) = NeArr.coeff(n);
      Ne.coeffRef(1, 2 * n + 1) = NeArr.coeff(n);
    }
  else if(Ne.rows() == 3)
    for(int n = 0; n < shape.n; n++)
    {
      Ne.coeffRef(0, 3 * n) = NeArr.coeff(n);
      Ne.coeffRef(1, 3 * n + 1) = NeArr.coeff(n);
      Ne.coeffRef(2, 3 * n + 2) = NeArr.coeff(n);
    }
  else
  {
    std::cerr << "Invalid Ne row in nmatrix_surface.hpp" << std::endl;
    exit(1);
  }

  // jacobian * gauss weight
  jac = sqrt(pow(dXdr.coeff(0, 0), 2.0) + pow(dXdr.coeff(0, 1), 2.0)) *
        gauss.weights(ip);
}

/// make surface shape function (surface in 3D)
template <class S, class G>
void snmatrix2d(Eigen::MatrixXd &Ne, double &jac, Eigen::MatrixXd &X, int ip)
{
  S shape;
  G gauss;
  Eigen::MatrixXd dXdr = shape.dNdr(gauss.points(ip)) * X;
  Eigen::VectorXd NeArr = shape.N(gauss.points(ip));

  if(Ne.rows() == 1)
    for(int n = 0; n < shape.n; n++)
      Ne.coeffRef(n) = NeArr.coeff(n);
  else if(Ne.rows() == 3)
    for(int n = 0; n < shape.n; n++)
    {
      Ne.coeffRef(0, 3 * n) = NeArr.coeff(n);
      Ne.coeffRef(1, 3 * n + 1) = NeArr.coeff(n);
      Ne.coeffRef(2, 3 * n + 2) = NeArr.coeff(n);
    }
  else
  {
    std::cerr << "Invalid row in nmatrix_surface.hpp" << std::endl;
    exit(1);
  }

  // jacobian * gauss weight
  Eigen::Vector3d dXdr1, dXdr2;
  dXdr1 << dXdr(0, 0), dXdr(0, 1), dXdr(0, 2);
  dXdr2 << dXdr(1, 0), dXdr(1, 1), dXdr(1, 2);

  jac = dXdr1.cross(dXdr2).norm() * gauss.weights(ip);
}

} // namespace icarat