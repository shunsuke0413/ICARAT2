///  @file  eigen_manipulation.cpp
///  @author  Daiki Watanabe
///  @date  November 20, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once

#include <Eigen/Core>
#include <iostream>
#include <problem/base.hpp>

namespace icarat
{
/// This is a reference to numpy.
inline Eigen::MatrixXd hstack(Eigen::MatrixXd const &a,
                              Eigen::MatrixXd const &b)
{
  Eigen::MatrixXd c(a.rows(), a.cols() + b.cols());
  c << a, b;
  return c;
}

/// This is a reference to numpy.
inline Eigen::MatrixXd vstack(Eigen::MatrixXd const &a,
                              Eigen::MatrixXd const &b)
{
  Eigen::MatrixXd c(a.rows() + b.rows(), a.cols());
  c << a, b;
  return c;
}

/// Conversion from Voigt to Tensor Notation
///@sa https://en.wikipedia.org/wiki/Voigt_notation
/// isStress ... true = stress conversion, false = strain conversion
inline Eigen::MatrixXd Tensor(Eigen::VectorXd const &voi, bool isStress)
{
  Eigen::MatrixXd ten;
  if(isStress)
  {
    if(voi.size() == 3)
    {
      ten.resize(2, 2);
      ten << voi[0], voi[2], // xx, xy
          voi[2], voi[1];    // xy, yy
    }
    else if(voi.size() == 6)
    {
      ten.resize(3, 3);
      ten << voi[0], voi[3], voi[5], // xx, xy, xz
          voi[3], voi[1], voi[4],    // xy, yy, yz
          voi[5], voi[4], voi[2];    // xz, yz, zz
    }
    else
    {
      std::cerr << "error in eigen_manipulation.hpp" << std::endl;
      exit(1);
    }
  }
  else // strain
  {
    if(voi.size() == 3)
    {
      ten.resize(2, 2);
      ten << voi[0], 0.5 * voi[2], // xx, xy
          0.5 * voi[2], voi[1];    // xy, yy
    }
    else if(voi.size() == 6)
    {
      ten.resize(3, 3);
      ten << voi[0], 0.5 * voi[3], 0.5 * voi[5], // xx, xy, xz
          0.5 * voi[3], voi[1], 0.5 * voi[4],    // xy, yy, yz
          0.5 * voi[5], 0.5 * voi[4], voi[2];    // xz, yz, zz
    }
    else
    {
      std::cerr << "error in eigen_manipulation.hpp" << std::endl;
      exit(1);
    }
  }
  return ten;
}

/// Conversion from Tensor to Voigt Notation
///@sa https://en.wikipedia.org/wiki/Voigt_notation
/// isStress ... true = stress conversion, false = strain conversion
inline Eigen::VectorXd Voigt(Eigen::MatrixXd const &ten, bool isStress)
{
  Eigen::VectorXd voi;
  if(isStress)
  {
    if(ten.rows() == 2)
    {
      voi.resize(3);
      voi << ten.coeff(0, 0), ten.coeff(1, 1), ten.coeff(0, 1); // xx, yy, xy
    }
    else if(ten.rows() == 3)
    {
      voi.resize(6);
      voi << voi.coeff(0, 0), voi.coeff(1, 1), voi.coeff(2, 2), voi.coeff(0, 1),
          voi.coeff(1, 2), voi.coeff(0, 2); // xx, yy, zz, xy, yz, xz
    }
    else
    {
      std::cerr << "error in eigen_manipulation.hpp" << std::endl;
      exit(1);
    }
  }
  else
  {
    if(ten.rows() == 2)
    {
      voi.resize(3);
      voi << ten.coeff(0, 0), 2.0 * ten.coeff(1, 1),
          2.0 * ten.coeff(0, 1); // xx, yy, xy
    }
    else if(ten.rows() == 3)
    {
      voi.resize(6);
      voi << voi.coeff(0, 0), voi.coeff(1, 1), voi.coeff(2, 2),
          2.0 * voi.coeff(0, 1), 2.0 * voi.coeff(1, 2),
          2.0 * voi.coeff(0, 2); // xx, yy, zz, xy, yz, xz
    }
    else
    {
      std::cerr << "error in eigen_manipulation.hpp" << std::endl;
      exit(1);
    }
  }
  return voi;
}

} // namespace icarat