#include "make_element.hpp"
#include <basis/nmatrix_surface.hpp>
#include <cmath>
#include <gtest/gtest.h>
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace icarat;
using namespace icarat::test;

TEST(SNmatrixTest, line2)
{
  int ndim = 2;
  int ne = 2;
  int ipmax = 2;
  string Etype = "line2";
  MatrixXd X = MatrixXd::Zero(ne, ndim);
  MatrixXd Ne = MatrixXd::Zero(ndim, ndim * ne);
  X << 0.0, 0.0, //
      1.0, 1.0;

  SurfaceNmatrix snmatrix;
  double volume = 0.0;
  for(int ip = 0; ip < ipmax; ip++)
  {
    double jac = 0.0;
    snmatrix.make(Ne, jac, Etype, ipmax, X, ip);
    volume += jac;
  }
  EXPECT_LT(abs(volume - sqrt(2.0)), 1.0e-10);
}

TEST(SNmatrixTest, line3)
{
  int ndim = 2;
  int ne = 3;
  int ipmax = 3;
  string Etype = "line3";
  MatrixXd X = MatrixXd::Zero(ne, ndim);
  MatrixXd Ne = MatrixXd::Zero(ndim, ndim * ne);
  MatrixXd Be = MatrixXd::Zero(1, ne);
  X << 0.0, 0.0, //
      1.0, 1.0,  //
      0.5, 0.5;

  SurfaceNmatrix snmatrix;
  double volume = 0.0;
  for(int ip = 0; ip < ipmax; ip++)
  {
    double jac = 0.0;
    snmatrix.make(Ne, jac, Etype, ipmax, X, ip);
    volume += jac;
  }
  EXPECT_LT(abs(volume - sqrt(2.0)), 1.0e-10);
}