///  @file  shape_function.hpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// For the list of shape functions, following references were refferred.
/// @sa K.J. Bathe, Finite Element Procedures, 1995.
/// @sa M. Maździarz, Unified Isoparametric 3D LagrangeFinite Elements, 2010.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <Eigen/Core>
#include <vector>
namespace icarat
{
//********************2NodesLine********************
class ShapeFunction2Line
{
public:
  const int d = 1;
  const int n = 2;
  Eigen::VectorXd N(const Eigen::VectorXd &r);
  Eigen::MatrixXd dNdr(const Eigen::VectorXd &r);
};

//********************3NodesLine********************

class ShapeFunction3Line
{
public:
  const int d = 1;
  const int n = 3;
  Eigen::VectorXd N(const Eigen::VectorXd &r);
  Eigen::MatrixXd dNdr(const Eigen::VectorXd &r);
};
//********************3NodesTriangle********************

class ShapeFunction3Triangle
{
public:
  const int d = 2;
  const int n = 3;
  Eigen::VectorXd N(const Eigen::VectorXd &r);
  Eigen::MatrixXd dNdr(const Eigen::VectorXd &r);
};

//********************6NodesTriangle********************

class ShapeFunction6Triangle
{
public:
  const int d = 2;
  const int n = 6;
  Eigen::VectorXd N(const Eigen::VectorXd &r);
  Eigen::MatrixXd dNdr(const Eigen::VectorXd &r);
};

//********************4NodesSquare********************

class ShapeFunction4Square
{
public:
  const int d = 2;
  const int n = 4;
  Eigen::VectorXd N(const Eigen::VectorXd &r);
  Eigen::MatrixXd dNdr(const Eigen::VectorXd &r);
};

//********************8NodesSquare********************

class ShapeFunction8Square
{
public:
  const int d = 2;
  const int n = 8;
  Eigen::VectorXd N(const Eigen::VectorXd &r);
  Eigen::MatrixXd dNdr(const Eigen::VectorXd &r);
};

//********************4NodesTetrahedron********************

class ShapeFunction4Tetrahedron
{
public:
  const int d = 3;
  const int n = 4;
  Eigen::VectorXd N(const Eigen::VectorXd &r);
  Eigen::MatrixXd dNdr(const Eigen::VectorXd &r);
};

//********************10NodesTetrahedron********************

class ShapeFunction10Tetrahedron
{
public:
  const int d = 3;
  const int n = 10;
  Eigen::VectorXd N(const Eigen::VectorXd &r);
  Eigen::MatrixXd dNdr(const Eigen::VectorXd &r);
};

//********************6NodesPrism********************

class ShapeFunction6Prism
{
public:
  const int d = 3;
  const int n = 6;
  Eigen::VectorXd N(const Eigen::VectorXd &r);
  Eigen::MatrixXd dNdr(const Eigen::VectorXd &r);
};

//********************15NodesPrism********************

class ShapeFunction15Prism
{
public:
  const int d = 3;
  const int n = 15;
  Eigen::VectorXd N(const Eigen::VectorXd &r);
  Eigen::MatrixXd dNdr(const Eigen::VectorXd &r);
};

//********************8NodesCubic********************

class ShapeFunction8Cubic
{
public:
  const int d = 3;
  const int n = 8;
  Eigen::VectorXd N(const Eigen::VectorXd &r);
  Eigen::MatrixXd dNdr(const Eigen::VectorXd &r);
};

//********************20NodesIsoParametricElement********************

class ShapeFunction20Cubic
{
public:
  const int d = 3;
  const int n = 20;
  Eigen::VectorXd N(const Eigen::VectorXd &r);
  Eigen::MatrixXd dNdr(const Eigen::VectorXd &r);
};
} // namespace icarat