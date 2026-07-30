///  @file  shape_function.hpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php
#include "shape_function.hpp"
#include <vector>

namespace icarat
{
Eigen::VectorXd ShapeFunction2Line::N(const Eigen::VectorXd &r)
{
  Eigen::VectorXd N = Eigen::VectorXd(n);
  N(0) = 0.5 * (1.0 - r(0));
  N(1) = 0.5 * (1.0 + r(0));
  return N;
}

Eigen::MatrixXd ShapeFunction2Line::dNdr(const Eigen::VectorXd &r)
{
  Eigen::MatrixXd dNdr = Eigen::MatrixXd(d, n);
  dNdr(0, 0) = -0.5;
  dNdr(0, 1) = 0.5;
  return dNdr;
}

Eigen::VectorXd ShapeFunction3Line::N(const Eigen::VectorXd &r)
{
  Eigen::VectorXd N = Eigen::VectorXd(n);
  N(0) = -0.5 * (1.0 - r(0)) * r(0);
  N(1) = 0.5 * r(0) * (1.0 + r(0));
  N(2) = (1.0 - r(0)) * (1.0 + r(0));
  return N;
}

Eigen::MatrixXd ShapeFunction3Line::dNdr(const Eigen::VectorXd &r)
{
  Eigen::MatrixXd dNdr = Eigen::MatrixXd(d, n);
  dNdr(0, 0) = -0.5 * (1.0 - 2.0 * r(0));
  dNdr(0, 1) = 0.5 * (1.0 + 2.0 * r(0));
  dNdr(0, 2) = -2.0 * r(0);
  return dNdr;
}

Eigen::VectorXd ShapeFunction3Triangle::N(const Eigen::VectorXd &r)
{
  Eigen::VectorXd N = Eigen::VectorXd(n);
  N(0) = r(0);
  N(1) = r(1);
  N(2) = 1.0 - r(0) - r(1);
  return N;
}

Eigen::MatrixXd ShapeFunction3Triangle::dNdr(const Eigen::VectorXd &r)
{
  Eigen::MatrixXd dNdr = Eigen::MatrixXd(d, n);
  dNdr(0, 0) = 1.0;
  dNdr(0, 1) = 0.0;
  dNdr(0, 2) = -1.0;
  dNdr(1, 0) = 0.0;
  dNdr(1, 1) = 1.0;
  dNdr(1, 2) = -1.0;
  return dNdr;
}

Eigen::VectorXd ShapeFunction6Triangle::N(const Eigen::VectorXd &r)
{
  Eigen::VectorXd N = Eigen::VectorXd(n);
  N(0) = r(0) * (2.0 * r(0) - 1.0);
  N(1) = r(1) * (2.0 * r(1) - 1.0);
  N(2) = (1.0 - r(0) - r(1)) * (1.0 - 2.0 * r(0) - 2.0 * r(1));
  N(3) = 4.0 * r(0) * r(1);
  N(4) = 4.0 * r(1) * (1.0 - r(0) - r(1));
  N(5) = 4.0 * (1.0 - r(0) - r(1)) * r(0);
  return N;
}

Eigen::MatrixXd ShapeFunction6Triangle::dNdr(const Eigen::VectorXd &r)
{
  Eigen::MatrixXd dNdr = Eigen::MatrixXd(d, n);
  dNdr(0, 0) = 4.0 * r(0) - 1.0;
  dNdr(0, 1) = 0.0;
  dNdr(0, 2) = -3.0 + 4.0 * r(0) + 4.0 * r(1);
  dNdr(0, 3) = 4.0 * r(1);
  dNdr(0, 4) = -4.0 * r(1);
  dNdr(0, 5) = 4.0 * (1.0 - 2.0 * r(0) - r(1));
  dNdr(1, 0) = 0.0;
  dNdr(1, 1) = 4.0 * r(1) - 1.0;
  dNdr(1, 2) = -3.0 + 4.0 * r(0) + 4.0 * r(1);
  dNdr(1, 3) = 4.0 * r(0);
  dNdr(1, 4) = 4.0 * (1.0 - r(0) - 2.0 * r(1));
  dNdr(1, 5) = -4.0 * r(0);
  return dNdr;
}

Eigen::VectorXd ShapeFunction4Square::N(const Eigen::VectorXd &r)
{
  Eigen::VectorXd N = Eigen::VectorXd(n);
  N(0) = 0.25 * (1.0 - r(0)) * (1.0 - r(1));
  N(1) = 0.25 * (1.0 + r(0)) * (1.0 - r(1));
  N(2) = 0.25 * (1.0 + r(0)) * (1.0 + r(1));
  N(3) = 0.25 * (1.0 - r(0)) * (1.0 + r(1));
  return N;
}

Eigen::MatrixXd ShapeFunction4Square::dNdr(const Eigen::VectorXd &r)
{
  Eigen::MatrixXd dNdr = Eigen::MatrixXd(d, n);
  dNdr(0, 0) = -0.25 * (1.0 - r(1));
  dNdr(0, 1) = 0.25 * (1.0 - r(1));
  dNdr(0, 2) = 0.25 * (1.0 + r(1));
  dNdr(0, 3) = -0.25 * (1.0 + r(1));
  dNdr(1, 0) = -0.25 * (1.0 - r(0));
  dNdr(1, 1) = -0.25 * (1.0 + r(0));
  dNdr(1, 2) = 0.25 * (1.0 + r(0));
  dNdr(1, 3) = 0.25 * (1.0 - r(0));
  return dNdr;
}

Eigen::VectorXd ShapeFunction8Square::N(const Eigen::VectorXd &r)
{
  Eigen::VectorXd N = Eigen::VectorXd(n);
  N(0) = 0.25 * (1.0 - r(0)) * (1.0 - r(1)) * (-r(0) - r(1) - 1.0);
  N(1) = 0.25 * (1.0 + r(0)) * (1.0 - r(1)) * (r(0) - r(1) - 1.0);
  N(2) = 0.25 * (1.0 + r(0)) * (1.0 + r(1)) * (r(0) + r(1) - 1.0);
  N(3) = 0.25 * (1.0 - r(0)) * (1.0 + r(1)) * (-r(0) + r(1) - 1.0);
  N(4) = 0.5 * (1.0 - r(0)) * (1.0 + r(0)) * (1.0 - r(1));
  N(5) = 0.5 * (1.0 + r(0)) * (1.0 + r(1)) * (1.0 - r(1));
  N(6) = 0.5 * (1.0 + r(0)) * (1.0 - r(0)) * (1.0 + r(1));
  N(7) = 0.5 * (1.0 - r(0)) * (1.0 + r(1)) * (1.0 - r(1));
  return N;
}

Eigen::MatrixXd ShapeFunction8Square::dNdr(const Eigen::VectorXd &r)
{
  Eigen::MatrixXd dNdr = Eigen::MatrixXd(d, n);
  dNdr(0, 0) = 0.25 * (-(1.0 - r(1)) * (-r(0) - r(1) - 1.0) -
                       (1.0 - r(0)) * (1.0 - r(1)));
  dNdr(0, 1) =
      0.25 * ((1.0 - r(1)) * (r(0) - r(1) - 1.0) + (1.0 + r(0)) * (1.0 - r(1)));
  dNdr(0, 2) =
      0.25 * ((1.0 + r(1)) * (r(0) + r(1) - 1.0) + (1.0 + r(0)) * (1.0 + r(1)));
  dNdr(0, 3) = 0.25 * (-(1.0 + r(1)) * (-r(0) + r(1) - 1.0) -
                       (1.0 - r(0)) * (1.0 + r(1)));
  dNdr(0, 4) = -r(0) * (1.0 - r(1));
  dNdr(0, 5) = 0.5 * (1.0 + r(1)) * (1.0 - r(1));
  dNdr(0, 6) = -r(0) * (1.0 + r(1));
  dNdr(0, 7) = -0.5 * (1.0 + r(1)) * (1.0 - r(1));

  dNdr(1, 0) = 0.25 * (-(1.0 - r(0)) * (-r(0) - r(1) - 1.0) -
                       (1.0 - r(0)) * (1.0 - r(1)));
  dNdr(1, 1) = 0.25 * (-(1.0 + r(0)) * (r(0) - r(1) - 1.0) -
                       (1.0 + r(0)) * (1.0 - r(1)));
  dNdr(1, 2) =
      0.25 * ((1.0 + r(0)) * (r(0) + r(1) - 1.0) + (1.0 + r(0)) * (1.0 + r(1)));
  dNdr(1, 3) = 0.25 * ((1.0 - r(0)) * (-r(0) + r(1) - 1.0) +
                       (1.0 - r(0)) * (1.0 + r(1)));
  dNdr(1, 4) = -0.5 * (1.0 + r(0)) * (1.0 - r(0));
  dNdr(1, 5) = -r(1) * (1.0 + r(0));
  dNdr(1, 6) = 0.5 * (1.0 + r(0)) * (1.0 - r(0));
  dNdr(1, 7) = -r(1) * (1.0 - r(0));
  return dNdr;
}

Eigen::VectorXd ShapeFunction4Tetrahedron::N(const Eigen::VectorXd &r)
{
  Eigen::VectorXd N = Eigen::VectorXd(n);
  N(0) = r(0);
  N(1) = r(2);
  N(2) = r(1);
  N(3) = 1.0 - r(0) - r(1) - r(2);
  return N;
}

Eigen::MatrixXd ShapeFunction4Tetrahedron::dNdr(const Eigen::VectorXd &r)
{
  Eigen::MatrixXd dNdr = Eigen::MatrixXd(d, n);
  dNdr(0, 0) = 1.0;
  dNdr(0, 1) = 0.0;
  dNdr(0, 2) = 0.0;
  dNdr(0, 3) = -1.0;
  dNdr(1, 0) = 0.0;
  dNdr(1, 1) = 0.0;
  dNdr(1, 2) = 1.0;
  dNdr(1, 3) = -1.0;
  dNdr(2, 0) = 0.0;
  dNdr(2, 1) = 1.0;
  dNdr(2, 2) = 0.0;
  dNdr(2, 3) = -1.0;
  return dNdr;
}

Eigen::VectorXd ShapeFunction10Tetrahedron::N(const Eigen::VectorXd &r)
{
  Eigen::VectorXd N = Eigen::VectorXd(n);
  N(4) = 4.0 * r(0) * (1.0 - r(0) - r(1) - r(2));
  N(5) = 4.0 * r(0) * r(1);
  N(6) = 4.0 * r(1) * (1.0 - r(0) - r(1) - r(2));
  N(7) = 4.0 * r(0) * r(2);
  N(8) = 4.0 * r(1) * r(2);
  N(6) = 4.0 * r(2) * (1.0 - r(0) - r(1) - r(2));
  N(0) = 1.0 - r(0) - r(1) - r(2) - 0.5 * N(4) - 0.5 * N(6) - 0.5 * N(9);
  N(1) = r(1) - 0.5 * N(5) - 0.5 * N(6) - 0.5 * N(8);
  N(2) = r(0) - 0.5 * N(4) - 0.5 * N(5) - 0.5 * N(7);
  N(3) = r(2) - 0.5 * N(7) - 0.5 * N(8) - 0.5 * N(9);

  return N;
}

Eigen::MatrixXd ShapeFunction10Tetrahedron::dNdr(const Eigen::VectorXd &r)
{
  Eigen::MatrixXd dNdr = Eigen::MatrixXd(d, n); // 3,10
  dNdr(0, 4) = 4.0 * (1.0 - r(0) - r(1) - r(2)) - 4.0 * r(0);
  dNdr(0, 5) = 4.0 * r(1);
  dNdr(0, 6) = -4.0 * r(1);
  dNdr(0, 7) = 4.0 * r(2);
  dNdr(0, 8) = 0.0;
  dNdr(0, 9) = -4.0 * r(2);
  dNdr(0, 0) = -1.0 - 0.5 * dNdr(0, 4) - 0.5 * dNdr(0, 6) - 0.5 * dNdr(0, 9);
  dNdr(0, 1) = 1.0 - 0.5 * dNdr(0, 4) - 0.5 * dNdr(0, 5) - 0.5 * dNdr(0, 7);
  dNdr(0, 2) = -0.5 * dNdr(0, 5) - 0.5 * dNdr(0, 6) - 0.5 * dNdr(0, 8);
  dNdr(0, 3) = -0.5 * dNdr(0, 7) - 0.5 * dNdr(0, 8) - 0.5 * dNdr(0, 9);

  dNdr(1, 4) = -4.0 * r(0);
  dNdr(1, 5) = 4.0 * r(0);
  dNdr(1, 6) = 4.0 * (1.0 - r(0) - r(1) - r(2)) - 4.0 * r(1);
  dNdr(1, 7) = 0.0;
  dNdr(1, 8) = 4.0 * r(2);
  dNdr(1, 9) = -4.0 * r(2);
  dNdr(1, 0) = -1.0 - 0.5 * dNdr(1, 4) - 0.5 * dNdr(1, 6) - 0.5 * dNdr(1, 9);
  dNdr(1, 1) = -0.5 * dNdr(1, 4) - 0.5 * dNdr(1, 5) - 0.5 * dNdr(1, 7);
  dNdr(1, 2) = 1.0 - 0.5 * dNdr(1, 5) - 0.5 * dNdr(1, 6) - 0.5 * dNdr(1, 8);
  dNdr(1, 3) = -0.5 * dNdr(1, 7) - 0.5 * dNdr(1, 8) - 0.5 * dNdr(1, 9);

  dNdr(2, 4) = -4.0 * r(0);
  dNdr(2, 5) = 0.0;
  dNdr(2, 6) = -4.0 * r(1);
  dNdr(2, 7) = 4.0 * r(0);
  dNdr(2, 8) = 4.0 * r(1);
  dNdr(2, 9) = 4.0 * (1.0 - r(0) - r(1) - r(2)) - 4.0 * r(2);
  dNdr(2, 0) = -1.0 - 0.5 * dNdr(2, 4) - 0.5 * dNdr(2, 6) - 0.5 * dNdr(2, 9);
  dNdr(2, 1) = -0.5 * dNdr(2, 4) - 0.5 * dNdr(2, 5) - 0.5 * dNdr(2, 7);
  dNdr(2, 2) = -0.5 * dNdr(2, 5) - 0.5 * dNdr(2, 6) - 0.5 * dNdr(2, 8);
  dNdr(2, 3) = 1.0 - 0.5 * dNdr(2, 7) - 0.5 * dNdr(2, 8) - 0.5 * dNdr(2, 9);

  return dNdr;
}

Eigen::VectorXd ShapeFunction6Prism::N(const Eigen::VectorXd &r)
{
  Eigen::VectorXd N = Eigen::VectorXd(n);
  N(0) = r(0) * (1.0 - r(2)) / 2.0;
  N(1) = r(1) * (1.0 - r(2)) / 2.0;
  N(2) = (1.0 - r(0) - r(1)) * (1.0 - r(2)) / 2.0;

  N(3) = r(0) * (1.0 + r(2)) / 2.0;
  N(4) = r(1) * (1.0 + r(2)) / 2.0;
  N(5) = (1.0 - r(0) - r(1)) * (1.0 + r(2)) / 2.0;
  return N;
}

Eigen::MatrixXd ShapeFunction6Prism::dNdr(const Eigen::VectorXd &r)
{
  Eigen::MatrixXd dNdr = Eigen::MatrixXd(d, n); // 3,6
  dNdr(0, 0) = (1.0 - r(2)) / 2.0;
  dNdr(0, 1) = 0.0;
  dNdr(0, 2) = -(1.0 - r(2)) / 2.0;
  dNdr(0, 3) = (1.0 + r(2)) / 2.0;
  dNdr(0, 4) = 0.0;
  dNdr(0, 5) = -(1.0 + r(2)) / 2.0;

  dNdr(1, 0) = 0.0;
  dNdr(1, 1) = (1.0 - r(2)) / 2.0;
  dNdr(1, 2) = -(1.0 - r(2)) / 2.0;
  dNdr(1, 3) = 0.0;
  dNdr(1, 4) = (1.0 + r(2)) / 2.0;
  dNdr(1, 5) = -(1.0 + r(2)) / 2.0;

  dNdr(2, 0) = -r(0) / 2.0;
  dNdr(2, 1) = -r(1) / 2.0;
  dNdr(2, 2) = -(1.0 - r(0) - r(1)) / 2.0;
  dNdr(2, 3) = r(0) / 2.0;
  dNdr(2, 4) = r(1) / 2.0;
  dNdr(2, 5) = (1.0 - r(0) - r(1)) / 2.0;

  return dNdr;
}

Eigen::VectorXd ShapeFunction15Prism::N(const Eigen::VectorXd &r)
{
  Eigen::VectorXd N = Eigen::VectorXd(n);
  double L1 = r(0);
  double L2 = r(1);
  double rr = r(2);

  double L3 = 1.0 - r(0) - r(1);

  N(0) = 0.5 * L1 * ((2.0 * L1 - 1.0) * (1.0 - rr) - (1.0 - rr * rr));
  N(1) = 0.5 * L2 * ((2.0 * L2 - 1.0) * (1.0 - rr) - (1.0 - rr * rr));
  N(2) = 0.5 * L3 * ((2.0 * L3 - 1.0) * (1.0 - rr) - (1.0 - rr * rr));

  N(3) = 0.5 * L1 * ((2.0 * L1 - 1.0) * (1.0 + rr) - (1.0 - rr * rr));
  N(4) = 0.5 * L2 * ((2.0 * L2 - 1.0) * (1.0 + rr) - (1.0 - rr * rr));
  N(5) = 0.5 * L3 * ((2.0 * L3 - 1.0) * (1.0 + rr) - (1.0 - rr * rr));

  N(6) = 2.0 * L1 * L2 * (1.0 - rr);
  N(7) = 2.0 * L2 * L3 * (1.0 - rr);
  N(8) = 2.0 * L1 * L3 * (1.0 - rr);

  N(9) = 2.0 * L1 * L2 * (1.0 + rr);
  N(10) = 2.0 * L2 * L3 * (1.0 + rr);
  N(11) = 2.0 * L1 * L3 * (1.0 + rr);

  N(12) = L1 * (1.0 - rr * rr);
  N(13) = L2 * (1.0 - rr * rr);
  N(14) = L3 * (1.0 - rr * rr);

  return N;
}

Eigen::MatrixXd ShapeFunction15Prism::dNdr(const Eigen::VectorXd &r)
{
  Eigen::MatrixXd dNdr = Eigen::MatrixXd(d, n); // 3,6
  double L1 = r(0);
  double L2 = r(1);
  double L3 = 1.0 - r(0) - r(1);
  double rr = r(2);
  // u1
  dNdr(0, 0) =
      0.5 * ((2 * L1 - 1) * (1 - rr) - (1 - rr * rr) + 2 * L1 * (1 - rr));
  dNdr(0, 1) = 0.0;
  dNdr(0, 2) =
      0.5 * (-(2 * L3 - 1) * (1 - rr) + (1 - rr * rr) - 2 * L3 * (1 - rr));

  dNdr(0, 3) =
      0.5 * ((2 * L1 - 1) * (1 + rr) - (1 - rr * rr) + 2 * L1 * (1 + rr));
  dNdr(0, 4) = 0.0;
  dNdr(0, 5) =
      0.5 * (-(2 * L3 - 1) * (1 + rr) + (1 - rr * rr) - 2 * L3 * (1 + rr));

  dNdr(0, 6) = 2 * L2 * (1 - rr);
  dNdr(0, 7) = -2 * L2 * (1 - rr);
  dNdr(0, 8) = 2 * L3 * (1 - rr) - 2 * L1 * (1 - rr);

  dNdr(0, 9) = 2 * L2 * (1 + rr);
  dNdr(0, 10) = -2 * L2 * (1 + rr);
  dNdr(0, 11) = 2 * L3 * (1 + rr) - 2 * L1 * (1 + rr);

  dNdr(0, 12) = 1 - rr * rr;
  dNdr(0, 13) = 0.0;
  dNdr(0, 14) = -1 + rr * rr;

  // u2
  dNdr(1, 0) = 0.0;
  dNdr(1, 1) =
      0.5 * ((2 * L2 - 1) * (1 - rr) - (1 - rr * rr) + 2 * L2 * (1 - rr));
  dNdr(1, 2) =
      0.5 * (-(2 * L3 - 1) * (1 - rr) + (1 - rr * rr) - 2 * L3 * (1 - rr));

  dNdr(1, 3) = 0.0;
  dNdr(1, 4) =
      0.5 * ((2 * L2 - 1) * (1 + rr) - (1 - rr * rr) + 2 * L2 * (1 + rr));
  dNdr(1, 5) =
      0.5 * (-(2 * L3 - 1) * (1 + rr) + (1 - rr * rr) - 2 * L3 * (1 + rr));

  dNdr(1, 6) = 2 * L1 * (1 - rr);
  dNdr(1, 7) = 2 * L3 * (1 - rr) - 2 * L2 * (1 - rr);
  dNdr(1, 8) = -2 * L1 * (1 - rr);

  dNdr(1, 9) = 2 * L1 * (1 + rr);
  dNdr(1, 10) = 2 * L3 * (1 + rr) - 2 * L2 * (1 + rr);
  dNdr(1, 11) = -2 * L1 * (1 + rr);

  dNdr(1, 12) = 0.0;
  dNdr(1, 13) = 1 - rr * rr;
  dNdr(1, 14) = -1 + rr * rr;

  // u3
  dNdr(2, 0) = 0.5 * L1 * (-2 * L1 + 1 + 2 * rr);
  dNdr(2, 1) = 0.5 * L2 * (-2 * L2 + 1 + 2 * rr);
  dNdr(2, 2) = 0.5 * L3 * (-2 * L3 + 1 + 2 * rr);

  dNdr(2, 3) = 0.5 * L1 * (2 * L1 - 1 + 2 * rr);
  dNdr(2, 4) = 0.5 * L2 * (2 * L2 - 1 + 2 * rr);
  dNdr(2, 5) = 0.5 * L3 * (2 * L3 - 1 + 2 * rr);

  dNdr(2, 6) = -2 * L2 * L1;
  dNdr(2, 7) = -2 * L2 * L3;
  dNdr(2, 8) = -2 * L3 * L1;

  dNdr(2, 9) = 2 * L1 * L2;
  dNdr(2, 10) = 2 * L2 * L3;
  dNdr(2, 11) = 2 * L1 * L3;

  dNdr(2, 12) = -2 * L1 * rr;
  dNdr(2, 13) = -2 * L2 * rr;
  dNdr(2, 14) = -2 * L3 * rr;

  return dNdr;
}

Eigen::VectorXd ShapeFunction8Cubic::N(const Eigen::VectorXd &r)
{
  Eigen::VectorXd N = Eigen::VectorXd(n);
  N(0) = 0.125 * (1.0 - r(0)) * (1.0 - r(1)) * (1.0 - r(2));
  N(1) = 0.125 * (1.0 + r(0)) * (1.0 - r(1)) * (1.0 - r(2));
  N(2) = 0.125 * (1.0 + r(0)) * (1.0 + r(1)) * (1.0 - r(2));
  N(3) = 0.125 * (1.0 - r(0)) * (1.0 + r(1)) * (1.0 - r(2));
  N(4) = 0.125 * (1.0 - r(0)) * (1.0 - r(1)) * (1.0 + r(2));
  N(5) = 0.125 * (1.0 + r(0)) * (1.0 - r(1)) * (1.0 + r(2));
  N(6) = 0.125 * (1.0 + r(0)) * (1.0 + r(1)) * (1.0 + r(2));
  N(7) = 0.125 * (1.0 - r(0)) * (1.0 + r(1)) * (1.0 + r(2));
  return N;
}

Eigen::MatrixXd ShapeFunction8Cubic::dNdr(const Eigen::VectorXd &r)
{
  Eigen::MatrixXd dNdr = Eigen::MatrixXd(d, n);
  dNdr(0, 0) = -0.125 * (1.0 - r(1)) * (1.0 - r(2));
  dNdr(0, 1) = 0.125 * (1.0 - r(1)) * (1.0 - r(2));
  dNdr(0, 2) = 0.125 * (1.0 + r(1)) * (1.0 - r(2));
  dNdr(0, 3) = -0.125 * (1.0 + r(1)) * (1.0 - r(2));
  dNdr(0, 4) = -0.125 * (1.0 - r(1)) * (1.0 + r(2));
  dNdr(0, 5) = 0.125 * (1.0 - r(1)) * (1.0 + r(2));
  dNdr(0, 6) = 0.125 * (1.0 + r(1)) * (1.0 + r(2));
  dNdr(0, 7) = -0.125 * (1.0 + r(1)) * (1.0 + r(2));

  dNdr(1, 0) = -0.125 * (1.0 - r(2)) * (1.0 - r(0));
  dNdr(1, 1) = -0.125 * (1.0 - r(2)) * (1.0 + r(0));
  dNdr(1, 2) = 0.125 * (1.0 - r(2)) * (1.0 + r(0));
  dNdr(1, 3) = 0.125 * (1.0 - r(2)) * (1.0 - r(0));
  dNdr(1, 4) = -0.125 * (1.0 + r(2)) * (1.0 - r(0));
  dNdr(1, 5) = -0.125 * (1.0 + r(2)) * (1.0 + r(0));
  dNdr(1, 6) = 0.125 * (1.0 + r(2)) * (1.0 + r(0));
  dNdr(1, 7) = 0.125 * (1.0 + r(2)) * (1.0 - r(0));

  dNdr(2, 0) = -0.125 * (1.0 - r(0)) * (1.0 - r(1));
  dNdr(2, 1) = -0.125 * (1.0 + r(0)) * (1.0 - r(1));
  dNdr(2, 2) = -0.125 * (1.0 + r(0)) * (1.0 + r(1));
  dNdr(2, 3) = -0.125 * (1.0 - r(0)) * (1.0 + r(1));
  dNdr(2, 4) = 0.125 * (1.0 - r(0)) * (1.0 - r(1));
  dNdr(2, 5) = 0.125 * (1.0 + r(0)) * (1.0 - r(1));
  dNdr(2, 6) = 0.125 * (1.0 + r(0)) * (1.0 + r(1));
  dNdr(2, 7) = 0.125 * (1.0 - r(0)) * (1.0 + r(1));
  return dNdr;
}

Eigen::VectorXd ShapeFunction20Cubic::N(const Eigen::VectorXd &r)
{
  Eigen::VectorXd N = Eigen::VectorXd(n);
  N(0) = -0.125 * (1.0 - r(0)) * (1.0 - r(1)) * (1.0 - r(2)) *
         (2.0 + r(0) + r(1) + r(2));
  N(1) = -0.125 * (1.0 + r(0)) * (1.0 - r(1)) * (1.0 - r(2)) *
         (2.0 - r(0) + r(1) + r(2));
  N(2) = -0.125 * (1.0 + r(0)) * (1.0 + r(1)) * (1.0 - r(2)) *
         (2.0 - r(0) - r(1) + r(2));
  N(3) = -0.125 * (1.0 - r(0)) * (1.0 + r(1)) * (1.0 - r(2)) *
         (2.0 + r(0) - r(1) + r(2));
  N(4) = -0.125 * (1.0 - r(0)) * (1.0 - r(1)) * (1.0 + r(2)) *
         (2.0 + r(0) + r(1) - r(2));
  N(5) = -0.125 * (1.0 + r(0)) * (1.0 - r(1)) * (1.0 + r(2)) *
         (2.0 - r(0) + r(1) - r(2));
  N(6) = -0.125 * (1.0 + r(0)) * (1.0 + r(1)) * (1.0 + r(2)) *
         (2.0 - r(0) - r(1) - r(2));
  N(7) = -0.125 * (1.0 - r(0)) * (1.0 + r(1)) * (1.0 + r(2)) *
         (2.0 + r(0) - r(1) - r(2));
  N(8) = 0.25 * (1.0 + r(0)) * (1.0 - r(0)) * (1.0 - r(1)) * (1.0 - r(2));
  N(9) = 0.25 * (1.0 + r(0)) * (1.0 + r(1)) * (1.0 - r(1)) * (1.0 - r(2));
  N(10) = 0.25 * (1.0 + r(0)) * (1.0 - r(0)) * (1.0 + r(1)) * (1.0 - r(2));
  N(11) = 0.25 * (1.0 - r(0)) * (1.0 + r(1)) * (1.0 - r(1)) * (1.0 - r(2));
  N(12) = 0.25 * (1.0 + r(0)) * (1.0 - r(0)) * (1.0 - r(1)) * (1.0 + r(2));
  N(13) = 0.25 * (1.0 + r(0)) * (1.0 + r(1)) * (1.0 - r(1)) * (1.0 + r(2));
  N(14) = 0.25 * (1.0 + r(0)) * (1.0 - r(0)) * (1.0 + r(1)) * (1.0 + r(2));
  N(15) = 0.25 * (1.0 - r(0)) * (1.0 + r(1)) * (1.0 - r(1)) * (1.0 + r(2));
  N(16) = 0.25 * (1.0 - r(0)) * (1.0 - r(1)) * (1.0 + r(2)) * (1.0 - r(2));
  N(17) = 0.25 * (1.0 + r(0)) * (1.0 - r(1)) * (1.0 + r(2)) * (1.0 - r(2));
  N(18) = 0.25 * (1.0 + r(0)) * (1.0 + r(1)) * (1.0 + r(2)) * (1.0 - r(2));
  N(19) = 0.25 * (1.0 - r(0)) * (1.0 + r(1)) * (1.0 + r(2)) * (1.0 - r(2));
  return N;
}

Eigen::MatrixXd ShapeFunction20Cubic::dNdr(const Eigen::VectorXd &r)
{
  Eigen::MatrixXd dNdr = Eigen::MatrixXd(d, n);
  dNdr(0, 0) =
      0.125 * (1.0 - r(1)) * (1.0 - r(2)) * (1.0 + 2.0 * r(0) + r(1) + r(2));
  dNdr(0, 1) =
      -0.125 * (1.0 - r(1)) * (1.0 - r(2)) * (1.0 - 2.0 * r(0) + r(1) + r(2));
  dNdr(0, 2) =
      -0.125 * (1.0 + r(1)) * (1.0 - r(2)) * (1.0 - 2.0 * r(0) - r(1) + r(2));
  dNdr(0, 3) =
      0.125 * (1.0 + r(1)) * (1.0 - r(2)) * (1.0 + 2.0 * r(0) - r(1) + r(2));
  dNdr(0, 4) =
      0.125 * (1.0 - r(1)) * (1.0 + r(2)) * (1.0 + 2.0 * r(0) + r(1) - r(2));
  dNdr(0, 5) =
      -0.125 * (1.0 - r(1)) * (1.0 + r(2)) * (1.0 - 2.0 * r(0) + r(1) - r(2));
  dNdr(0, 6) =
      -0.125 * (1.0 + r(1)) * (1.0 + r(2)) * (1.0 - 2.0 * r(0) - r(1) - r(2));
  dNdr(0, 7) =
      0.125 * (1.0 + r(1)) * (1.0 + r(2)) * (1.0 + 2.0 * r(0) - r(1) - r(2));
  dNdr(0, 8) = -0.5 * r(0) * (1.0 - r(1)) * (1.0 - r(2));
  dNdr(0, 9) = 0.25 * (1.0 - r(1) * r(1)) * (1.0 - r(2));
  dNdr(0, 10) = -0.5 * r(0) * (1.0 + r(1)) * (1.0 - r(2));
  dNdr(0, 11) = -0.25 * (1.0 - r(1) * r(1)) * (1.0 - r(2));
  dNdr(0, 12) = -0.5 * r(0) * (1.0 - r(1)) * (1.0 + r(2));
  dNdr(0, 13) = 0.25 * (1.0 - r(1) * r(1)) * (1.0 + r(2));
  dNdr(0, 14) = -0.5 * r(0) * (1.0 + r(1)) * (1.0 + r(2));
  dNdr(0, 15) = -0.25 * (1.0 - r(1) * r(1)) * (1.0 + r(2));
  dNdr(0, 16) = -0.25 * (1.0 - r(1)) * (1.0 - r(2) * r(2));
  dNdr(0, 17) = 0.25 * (1.0 - r(1)) * (1.0 - r(2) * r(2));
  dNdr(0, 18) = 0.25 * (1.0 + r(1)) * (1.0 - r(2) * r(2));
  dNdr(0, 19) = -0.25 * (1.0 + r(1)) * (1.0 - r(2) * r(2));

  dNdr(1, 0) =
      0.125 * (1.0 - r(0)) * (1.0 - r(2)) * (1.0 + r(0) + 2.0 * r(1) + r(2));
  dNdr(1, 1) =
      0.125 * (1.0 + r(0)) * (1.0 - r(2)) * (1.0 - r(0) + 2.0 * r(1) + r(2));
  dNdr(1, 2) =
      -0.125 * (1.0 + r(0)) * (1.0 - r(2)) * (1.0 - r(0) - 2.0 * r(1) + r(2));
  dNdr(1, 3) =
      -0.125 * (1.0 - r(0)) * (1.0 - r(2)) * (1.0 + r(0) - 2.0 * r(1) + r(2));
  dNdr(1, 4) =
      0.125 * (1.0 - r(0)) * (1.0 + r(2)) * (1.0 + r(0) + 2.0 * r(1) - r(2));
  dNdr(1, 5) =
      0.125 * (1.0 + r(0)) * (1.0 + r(2)) * (1.0 - r(0) + 2.0 * r(1) - r(2));
  dNdr(1, 6) =
      -0.125 * (1.0 + r(0)) * (1.0 + r(2)) * (1.0 - r(0) - 2.0 * r(1) - r(2));
  dNdr(1, 7) =
      -0.125 * (1.0 - r(0)) * (1.0 + r(2)) * (1.0 + r(0) - 2.0 * r(1) - r(2));
  dNdr(1, 8) = -0.25 * (1.0 - r(0) * r(0)) * (1.0 - r(2));
  dNdr(1, 9) = -0.5 * (1.0 + r(0)) * r(1) * (1.0 - r(2));
  dNdr(1, 10) = 0.25 * (1.0 - r(0) * r(0)) * (1.0 - r(2));
  dNdr(1, 11) = -0.5 * (1.0 - r(0)) * r(1) * (1.0 - r(2));
  dNdr(1, 12) = -0.25 * (1.0 - r(0) * r(0)) * (1.0 + r(2));
  dNdr(1, 13) = -0.5 * (1.0 + r(0)) * r(1) * (1.0 + r(2));
  dNdr(1, 14) = 0.25 * (1.0 - r(0) * r(0)) * (1.0 + r(2));
  dNdr(1, 15) = -0.5 * (1.0 - r(0)) * r(1) * (1.0 + r(2));
  dNdr(1, 16) = -0.25 * (1.0 - r(0)) * (1.0 - r(2) * r(2));
  dNdr(1, 17) = -0.25 * (1.0 + r(0)) * (1.0 - r(2) * r(2));
  dNdr(1, 18) = 0.25 * (1.0 + r(0)) * (1.0 - r(2) * r(2));
  dNdr(1, 19) = 0.25 * (1.0 - r(0)) * (1.0 - r(2) * r(2));

  dNdr(2, 0) =
      0.125 * (1.0 - r(0)) * (1.0 - r(1)) * (1.0 + r(0) + r(1) + 2.0 * r(2));
  dNdr(2, 1) =
      0.125 * (1.0 + r(0)) * (1.0 - r(1)) * (1.0 - r(0) + r(1) + 2.0 * r(2));
  dNdr(2, 2) =
      0.125 * (1.0 + r(0)) * (1.0 + r(1)) * (1.0 - r(0) - r(1) + 2.0 * r(2));
  dNdr(2, 3) =
      0.125 * (1.0 - r(0)) * (1.0 + r(1)) * (1.0 + r(0) - r(1) + 2.0 * r(2));
  dNdr(2, 4) =
      -0.125 * (1.0 - r(0)) * (1.0 - r(1)) * (1.0 + r(0) + r(1) - 2.0 * r(2));
  dNdr(2, 5) =
      -0.125 * (1.0 + r(0)) * (1.0 - r(1)) * (1.0 - r(0) + r(1) - 2.0 * r(2));
  dNdr(2, 6) =
      -0.125 * (1.0 + r(0)) * (1.0 + r(1)) * (1.0 - r(0) - r(1) - 2.0 * r(2));
  dNdr(2, 7) =
      -0.125 * (1.0 - r(0)) * (1.0 + r(1)) * (1.0 + r(0) - r(1) - 2.0 * r(2));
  dNdr(2, 8) = -0.25 * (1.0 - r(0) * r(0)) * (1.0 - r(1));
  dNdr(2, 9) = -0.25 * (1.0 + r(0)) * (1.0 - r(1) * r(1));
  dNdr(2, 10) = -0.25 * (1.0 - r(0) * r(0)) * (1.0 + r(1));
  dNdr(2, 11) = -0.25 * (1.0 - r(0)) * (1.0 - r(1) * r(1));
  dNdr(2, 12) = 0.25 * (1.0 - r(0) * r(0)) * (1.0 - r(1));
  dNdr(2, 13) = 0.25 * (1.0 + r(0)) * (1.0 - r(1) * r(1));
  dNdr(2, 14) = 0.25 * (1.0 - r(0) * r(0)) * (1.0 + r(1));
  dNdr(2, 15) = 0.25 * (1.0 - r(0)) * (1.0 - r(1) * r(1));
  dNdr(2, 16) = -0.5 * (1.0 - r(0)) * (1.0 - r(1)) * r(2);
  dNdr(2, 17) = -0.5 * (1.0 + r(0)) * (1.0 - r(1)) * r(2);
  dNdr(2, 18) = -0.5 * (1.0 + r(0)) * (1.0 + r(1)) * r(2);
  dNdr(2, 19) = -0.5 * (1.0 - r(0)) * (1.0 + r(1)) * r(2);
  return dNdr;
}
} // namespace icarat
