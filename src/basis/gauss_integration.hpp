///  @file  gauss_integration.hpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php
/// For the list of shape functions, following references were refferred.
/// @sa M. Maździarz, Unified Isoparametric 3D LagrangeFinite Elements, 2010.
#pragma once
#include <Eigen/Core>
#include <cmath>
#include <vector>

namespace icarat
{
//********************GaussIntegrationConstant1D2PointsForLine********************
class Gauss2Line
{
public:
  Gauss2Line();
  static const int N = 2; // Number of integration point
  const Eigen::VectorXd &points(int ip) const
  {
    return points_[ip];
  } ///< Cordinate values of integration point
  const double &weights(int ip) const
  {
    return weights_[ip];
  } ///< Cordinate values of integration point
private:
  std::vector<Eigen::VectorXd> points_;
  std::vector<double> weights_;
};

//********************GaussIntegrationConstant1D3PointsForLine********************
class Gauss3Line
{
public:
  Gauss3Line();
  static const int N = 3; // Number of integration point
  const Eigen::VectorXd &points(int ip) const
  {
    return points_[ip];
  } ///< Cordinate values of integration point
  const double &weights(int ip) const
  {
    return weights_[ip];
  } ///< Cordinate values of integration point
private:
  std::vector<Eigen::VectorXd> points_;
  std::vector<double> weights_;
};

///******************GaussIntegrationConstant2D1PointsForTriangle********************
class Gauss1Triangle
{
public:
  Gauss1Triangle();
  const int N = 1; ///< Number of integration point
  const Eigen::VectorXd &points(int ip) const
  {
    return points_[ip];
  } ///< Cordinate values of integration point
  const double &weights(int ip) const
  {
    return weights_[ip];
  } ///< Cordinate values of integration point
private:
  std::vector<Eigen::VectorXd> points_;
  std::vector<double> weights_;
};

///******************GaussIntegrationConstant2D3PointsForTriangle********************
class Gauss3Triangle
{
public:
  Gauss3Triangle();
  const int N = 3; ///< Number of integration point
  const Eigen::VectorXd &points(int ip) const
  {
    return points_[ip];
  } ///< Cordinate values of integration point
  const double &weights(int ip) const
  {
    return weights_[ip];
  } ///< Cordinate values of integration point
private:
  std::vector<Eigen::VectorXd> points_;
  std::vector<double> weights_;
};

///******************GaussIntegrationConstant2D1Points********************
class Gauss1Square
{
public:
  Gauss1Square();
  const int N = 1; ///< Number of integration point
  const Eigen::VectorXd &points(int ip) const
  {
    return points_[ip];
  } ///< Cordinate values of integration point
  const double &weights(int ip) const
  {
    return weights_[ip];
  } ///< Cordinate values of integration point
private:
  std::vector<Eigen::VectorXd> points_;
  std::vector<double> weights_;
};

///******************GaussIntegrationConstant2D4Points********************
class Gauss4Square
{
public:
  Gauss4Square();
  const int N = 4; ///< Number of integration point
  const Eigen::VectorXd &points(int ip) const
  {
    return points_[ip];
  } ///< Cordinate values of integration point
  const double &weights(int ip) const
  {
    return weights_[ip];
  } ///< Cordinate values of integration point
private:
  std::vector<Eigen::VectorXd> points_;
  std::vector<double> weights_;
};

///******************GaussIntegrationConstant2D9Points********************
class Gauss9Square
{
public:
  Gauss9Square();
  const int N = 9; ///< Number of integration point
  const Eigen::VectorXd &points(int ip) const
  {
    return points_[ip];
  } ///< Cordinate values of integration point
  const double &weights(int ip) const
  {
    return weights_[ip];
  } ///< Cordinate values of integration point
private:
  std::vector<Eigen::VectorXd> points_;
  std::vector<double> weights_;
};

///******************GaussIntegrationConstant3D1PointsForTetrahedron********************
class Gauss1Tetrahedron
{
public:
  Gauss1Tetrahedron();
  const int N = 1; ///< Number of integration point
  const Eigen::VectorXd &points(int ip) const
  {
    return points_[ip];
  } ///< Cordinate values of integration point
  const double &weights(int ip) const
  {
    return weights_[ip];
  } ///< Cordinate values of integration point
private:
  std::vector<Eigen::VectorXd> points_;
  std::vector<double> weights_;
};

///******************GaussIntegrationConstant3D1PointsForTetrahedron********************
class Gauss4Tetrahedron
{
public:
  Gauss4Tetrahedron();
  const int N = 4; ///< Number of integration point
  const Eigen::VectorXd &points(int ip) const
  {
    return points_[ip];
  } ///< Cordinate values of integration point
  const double &weights(int ip) const
  {
    return weights_[ip];
  } ///< Cordinate values of integration point
private:
  std::vector<Eigen::VectorXd> points_;
  std::vector<double> weights_;
};

///******************GaussIntegrationConstant3D2PointsForPrism********************
class Gauss6Prism
{
public:
  Gauss6Prism();
  const int N = 6; ///< Number of integration point
  const Eigen::VectorXd &points(int ip) const
  {
    return points_[ip];
  } ///< Cordinate values of integration point
  const double &weights(int ip) const
  {
    return weights_[ip];
  } ///< Cordinate values of integration point
private:
  std::vector<Eigen::VectorXd> points_;
  std::vector<double> weights_;
};

///******************GaussIntegrationConstant3D9PointsForPrism********************
class Gauss21Prism
{
public:
  Gauss21Prism();
  const int N = 21; ///< Number of integration point
  const Eigen::VectorXd &points(int ip) const
  {
    return points_[ip];
  } ///< Cordinate values of integration point
  const double &weights(int ip) const
  {
    return weights_[ip];
  } ///< Cordinate values of integration point
private:
  std::vector<Eigen::VectorXd> points_;
  std::vector<double> weights_;
};

///******************GaussIntegrationConstant3D8Points********************

class Gauss8Cubic
{
public:
  Gauss8Cubic();
  const int N = 8; ///< Number of integration point
  const Eigen::VectorXd &points(int ip) const
  {
    return points_[ip];
  } ///< Cordinate values of integration point
  const double &weights(int ip) const
  {
    return weights_[ip];
  } ///< Cordinate values of integration point
private:
  std::vector<Eigen::VectorXd> points_;
  std::vector<double> weights_;
};

///******************GaussIntegrationConstant27D8Points********************

class Gauss27Cubic
{
public:
  Gauss27Cubic();
  const int N = 27; ///< Number of integration point
  const Eigen::VectorXd &points(int ip) const
  {
    return points_[ip];
  } ///< Cordinate values of integration point
  const double &weights(int ip) const
  {
    return weights_[ip];
  } ///< Cordinate values of integration point
private:
  std::vector<Eigen::VectorXd> points_;
  std::vector<double> weights_;
};

} // namespace icarat