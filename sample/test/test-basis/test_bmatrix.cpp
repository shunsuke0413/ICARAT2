#include "make_element.hpp"
#include <basis/bmatrix.hpp>
#include <cmath>
#include <gtest/gtest.h>
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace icarat;
using namespace icarat::test;

TEST(BmatrixTest, quad4)
{
  FEM fem;
  vector<Element> element;
  vector<Node> node;
  make1quad4(fem, element, node);

  MatrixXd Be = MatrixXd::Zero(fem.voigt, element[0].numdof);
  MatrixXd X = MatrixXd::Zero(element[0].ne, fem.ndim);
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, element[0].numdof);
  for(int i = 0; i < element[0].ne; i++)
    for(int j = 0; j < fem.ndim; j++)
      X(i, j) = node[element[0].nodeID[i]].x[j];

  Bmatrix bmatrix;
  double volume = 0.0;
  for(int ip = 0; ip < element[0].ipmax; ip++)
  {
    double jac = 0.0;
    bmatrix.make(Be, Ne, jac, element[0].eType, element[0].ipmax, X, ip);
    volume += jac;
  }
  EXPECT_LT(abs(volume - 1.0), 1.0e-10);

  // reduced integral
  int ipmax = 1;
  volume = 0.0;
  for(int ip = 0; ip < ipmax; ip++)
  {
    double jac = 0.0;
    bmatrix.make(Be, Ne, jac, element[0].eType, ipmax, X, ip);
    volume += jac;
  }
  EXPECT_LT(abs(volume - 1.0), 1.0e-10);
}

TEST(BmatrixTest, quad8)
{
  FEM fem;
  vector<Element> element;
  vector<Node> node;
  make1quad8(fem, element, node);

  MatrixXd Be = MatrixXd::Zero(fem.voigt, element[0].numdof);
  MatrixXd X = MatrixXd::Zero(element[0].ne, fem.ndim);
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, element[0].numdof);
  for(int i = 0; i < element[0].ne; i++)
    for(int j = 0; j < fem.ndim; j++)
      X(i, j) = node[element[0].nodeID[i]].x[j];

  Bmatrix bmatrix;
  double volume = 0.0;
  for(int ip = 0; ip < element[0].ipmax; ip++)
  {
    double jac = 0.0;
    bmatrix.make(Be, Ne, jac, element[0].eType, element[0].ipmax, X, ip);
    volume += jac;
  }
  EXPECT_LT(abs(volume - 1.0), 1.0e-10);
}

TEST(BmatrixTest, tria3)
{
  FEM fem;
  vector<Element> element;
  vector<Node> node;
  make1tria3(fem, element, node);

  MatrixXd Be = MatrixXd::Zero(fem.voigt, element[0].numdof);
  MatrixXd X = MatrixXd::Zero(element[0].ne, fem.ndim);
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, element[0].numdof);
  for(int i = 0; i < element[0].ne; i++)
    for(int j = 0; j < fem.ndim; j++)
      X(i, j) = node[element[0].nodeID[i]].x[j];

  Bmatrix bmatrix;
  double volume = 0.0;
  for(int ip = 0; ip < element[0].ipmax; ip++)
  {
    double jac = 0.0;
    bmatrix.make(Be, Ne, jac, element[0].eType, element[0].ipmax, X, ip);
    volume += jac;
  }
  EXPECT_LT(abs(volume - 0.5), 1.0e-10);
}

TEST(BmatrixTest, tria6)
{
  FEM fem;
  vector<Element> element;
  vector<Node> node;
  make1tria6(fem, element, node);

  MatrixXd Be = MatrixXd::Zero(fem.voigt, element[0].numdof);
  MatrixXd X = MatrixXd::Zero(element[0].ne, fem.ndim);
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, element[0].numdof);
  for(int i = 0; i < element[0].ne; i++)
    for(int j = 0; j < fem.ndim; j++)
      X(i, j) = node[element[0].nodeID[i]].x[j];

  Bmatrix bmatrix;
  double volume = 0.0;
  for(int ip = 0; ip < element[0].ipmax; ip++)
  {
    double jac = 0.0;
    bmatrix.make(Be, Ne, jac, element[0].eType, element[0].ipmax, X, ip);
    volume += jac;
  }
  EXPECT_LT(abs(volume - 0.5), 1.0e-10);
}

TEST(BmatrixTest, cubic8)
{
  FEM fem;
  vector<Element> element;
  vector<Node> node;
  make1cubic8(fem, element, node);

  MatrixXd Be = MatrixXd::Zero(fem.voigt, element[0].numdof);
  MatrixXd X = MatrixXd::Zero(element[0].ne, fem.ndim);
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, element[0].numdof);
  for(int i = 0; i < element[0].ne; i++)
    for(int j = 0; j < fem.ndim; j++)
      X(i, j) = node[element[0].nodeID[i]].x[j];

  Bmatrix bmatrix;
  double volume = 0.0;
  for(int ip = 0; ip < element[0].ipmax; ip++)
  {
    double jac = 0.0;
    bmatrix.make(Be, Ne, jac, element[0].eType, element[0].ipmax, X, ip);
    volume += jac;
  }
  EXPECT_LT(abs(volume - 1.0), 1.0e-10);
}

TEST(BmatrixTest, cubic20)
{
  FEM fem;
  vector<Element> element;
  vector<Node> node;
  make1cubic20(fem, element, node);

  MatrixXd Be = MatrixXd::Zero(fem.voigt, element[0].numdof);
  MatrixXd X = MatrixXd::Zero(element[0].ne, fem.ndim);
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, element[0].numdof);
  for(int i = 0; i < element[0].ne; i++)
    for(int j = 0; j < fem.ndim; j++)
      X(i, j) = node[element[0].nodeID[i]].x[j];

  Bmatrix bmatrix;
  double volume = 0.0;
  for(int ip = 0; ip < element[0].ipmax; ip++)
  {
    double jac = 0.0;
    bmatrix.make(Be, Ne, jac, element[0].eType, element[0].ipmax, X, ip);
    volume += jac;
  }
  EXPECT_LT(abs(volume - 1.0), 1.0e-10);
}

TEST(BmatrixTest, tetra4)
{
  FEM fem;
  vector<Element> element;
  vector<Node> node;
  make1tetra4(fem, element, node);

  MatrixXd Be = MatrixXd::Zero(fem.voigt, element[0].numdof);
  MatrixXd X = MatrixXd::Zero(element[0].ne, fem.ndim);
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, element[0].numdof);
  for(int i = 0; i < element[0].ne; i++)
    for(int j = 0; j < fem.ndim; j++)
      X(i, j) = node[element[0].nodeID[i]].x[j];

  Bmatrix bmatrix;
  double volume = 0.0;
  for(int ip = 0; ip < element[0].ipmax; ip++)
  {
    double jac = 0.0;
    bmatrix.make(Be, Ne, jac, element[0].eType, element[0].ipmax, X, ip);
    volume += jac;
  }

  EXPECT_LT(abs(volume - (1.0 / 6.0)), 1.0e-10);
}

TEST(BmatrixTest, tetra10)
{
  FEM fem;
  vector<Element> element;
  vector<Node> node;
  make1tetra10(fem, element, node);

  MatrixXd Be = MatrixXd::Zero(fem.voigt, element[0].numdof);
  MatrixXd X = MatrixXd::Zero(element[0].ne, fem.ndim);
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, element[0].numdof);
  for(int i = 0; i < element[0].ne; i++)
    for(int j = 0; j < fem.ndim; j++)
      X(i, j) = node[element[0].nodeID[i]].x[j];

  Bmatrix bmatrix;
  double volume = 0.0;
  for(int ip = 0; ip < element[0].ipmax; ip++)
  {
    double jac = 0.0;
    bmatrix.make(Be, Ne, jac, element[0].eType, element[0].ipmax, X, ip);
    volume += jac;
  }

  EXPECT_LT(abs(volume - (1.0 / 6.0)), 1.0e-8);
}

TEST(BmatrixTest, prism6)
{
  FEM fem;
  vector<Element> element;
  vector<Node> node;
  make1prism6(fem, element, node);

  MatrixXd Be = MatrixXd::Zero(fem.voigt, element[0].numdof);
  MatrixXd X = MatrixXd::Zero(element[0].ne, fem.ndim);
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, element[0].numdof);
  for(int i = 0; i < element[0].ne; i++)
    for(int j = 0; j < fem.ndim; j++)
      X(i, j) = node[element[0].nodeID[i]].x[j];

  Bmatrix bmatrix;
  double volume = 0.0;
  for(int ip = 0; ip < element[0].ipmax; ip++)
  {
    double jac = 0.0;
    bmatrix.make(Be, Ne, jac, element[0].eType, element[0].ipmax, X, ip);
    volume += jac;
  }
  EXPECT_LT(abs(volume - 0.5), 1.0e-9);
}

TEST(BmatrixTest, prism15)
{
  FEM fem;
  vector<Element> element;
  vector<Node> node;
  make1prism15(fem, element, node);

  MatrixXd Be = MatrixXd::Zero(fem.voigt, element[0].numdof);
  MatrixXd X = MatrixXd::Zero(element[0].ne, fem.ndim);
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, element[0].numdof);
  for(int i = 0; i < element[0].ne; i++)
    for(int j = 0; j < fem.ndim; j++)
      X(i, j) = node[element[0].nodeID[i]].x[j];

  Bmatrix bmatrix;
  double volume = 0.0;
  for(int ip = 0; ip < element[0].ipmax; ip++)
  {
    double jac = 0.0;
    bmatrix.make(Be, Ne, jac, element[0].eType, element[0].ipmax, X, ip);
    volume += jac;
  }
  EXPECT_LT(abs(volume - 0.5), 1.0e-8);
}
