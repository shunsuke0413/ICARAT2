///  @file gauss_integration.cpp
///  @author	Daiki Watanabe
///  @date		May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "gauss_integration.hpp"
namespace icarat
{

Gauss2Line::Gauss2Line()
{
  Eigen::VectorXd point(1);
  point << -1.0 / sqrt(3.0);
  points_.push_back(point);
  point << 1.0 / sqrt(3.0);
  points_.push_back(point);

  weights_.push_back(1.0);
  weights_.push_back(1.0);
}

Gauss3Line::Gauss3Line()
{
  Eigen::VectorXd point(1);
  point << -sqrt(0.6);
  points_.push_back(point);
  point << 0.0;
  points_.push_back(point);
  point << sqrt(0.6);
  points_.push_back(point);

  weights_.push_back(5.0 / 9.0);
  weights_.push_back(8.0 / 9.0);
  weights_.push_back(5.0 / 9.0);
}

Gauss1Triangle::Gauss1Triangle()
{
  Eigen::VectorXd point(2);
  point << 1.0 / 3.0, 1.0 / 3.0;
  points_.push_back(point);

  weights_.push_back(1.0 / 2.0);
}

Gauss3Triangle::Gauss3Triangle()
{
  Eigen::VectorXd point(2);
  point << 1.0 / 6.0, 1.0 / 6.0;
  points_.push_back(point);
  point << 2.0 / 3.0, 1.0 / 6.0;
  points_.push_back(point);
  point << 1.0 / 6.0, 2.0 / 3.0;
  points_.push_back(point);

  weights_.push_back(1.0 / 6.0);
  weights_.push_back(1.0 / 6.0);
  weights_.push_back(1.0 / 6.0);
}

Gauss1Square::Gauss1Square()
{
  Eigen::VectorXd point(2);
  point << 0.0, 0.0;
  points_.push_back(point);

  weights_.push_back(2.0 * 2.0);
}

Gauss4Square::Gauss4Square()
{
  Eigen::VectorXd point(8);
  point << -1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0),
      -1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0),
      1.0 / sqrt(3.0);
  for(int i = 0; i < 4; i++)
  {
    Eigen::VectorXd tmp(2);
    for(int j = 0; j < 2; j++)
      tmp[j] = point(2 * i + j);
    points_.push_back(tmp);
  }

  for(int i = 0; i < 4; i++)
    weights_.push_back(1.0 * 1.0);
}

Gauss9Square::Gauss9Square()
{
  Eigen::VectorXd point(18);
  point << -sqrt(0.6), -sqrt(0.6), //
      0.0, -sqrt(0.6),             //
      sqrt(0.6), -sqrt(0.6),       //
      -sqrt(0.6), 0.0,             //
      0.0, 0.0,                    //
      sqrt(0.6), 0.0,              //
      -sqrt(0.6), sqrt(0.6),       //
      0.0, sqrt(0.6),              //
      sqrt(0.6), sqrt(0.6);
  for(int i = 0; i < 9; i++)
  {
    Eigen::VectorXd tmp(2);
    for(int j = 0; j < 2; j++)
      tmp[j] = point(2 * i + j);
    points_.push_back(tmp);
  }

  Eigen::VectorXd weight(18);
  weight << 5.0 / 9.0, 5.0 / 9.0, //
      8.0 / 9.0, 5.0 / 9.0,       //
      5.0 / 9.0, 5.0 / 9.0,       //
      5.0 / 9.0, 8.0 / 9.0,       //
      8.0 / 9.0, 8.0 / 9.0,       //
      5.0 / 9.0, 8.0 / 9.0,       //
      5.0 / 9.0, 5.0 / 9.0,       //
      8.0 / 9.0, 5.0 / 9.0,       //
      5.0 / 9.0, 5.0 / 9.0;
  for(int i = 0; i < 9; i++)
  {
    double tmp = 1.0;
    for(int j = 0; j < 2; j++)
      tmp *= weight(2 * i + j);
    weights_.push_back(tmp);
  }
}

Gauss1Tetrahedron::Gauss1Tetrahedron()
{
  Eigen::VectorXd point(3);
  point << 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0;
  points_.push_back(point);

  weights_.push_back(1.0 / 6.0);
}

Gauss4Tetrahedron::Gauss4Tetrahedron()
{
  Eigen::VectorXd point(3);
  point << 0.58541020, 0.13819660, 0.13819660;
  points_.push_back(point);
  point << 0.13819660, 0.58541020, 0.13819660;
  points_.push_back(point);
  point << 0.13819660, 0.13819660, 0.58541020;
  points_.push_back(point);
  point << 0.13819660, 0.13819660, 0.13819660;
  points_.push_back(point);

  for(int i = 0; i < 4; i++)
    weights_.push_back(0.041666667);
}

Gauss6Prism::Gauss6Prism()
{
  Eigen::VectorXd point(3);
  point << 1.0 / 2.0, 0.0, -0.5773502692;
  points_.push_back(point);

  point << 0.0, 1.0 / 2.0, -0.5773502692;
  points_.push_back(point);

  point << 1.0 / 2.0, 1.0 / 2.0, -0.5773502692;
  points_.push_back(point);

  point << 1.0 / 2.0, 0.0, 0.5773502692;
  points_.push_back(point);

  point << 0.0, 1.0 / 2.0, 0.5773502692;
  points_.push_back(point);

  point << 1.0 / 2.0, 1.0 / 2.0, 0.5773502692;
  points_.push_back(point);

  for(int i = 0; i < 6; i++)
    weights_.push_back(1.0 / 6.0);
}

Gauss21Prism::Gauss21Prism()
{
  Eigen::VectorXd point(3);
  // 1
  point << 0.3333333333, 0.3333333333, -0.7745696692;
  points_.push_back(point);
  // 2
  point << 0.0597158718, 0.4701420641, -0.7745696692;
  points_.push_back(point);
  // 3
  point << 0.4701420641, 0.0597158718, -0.7745696692;
  points_.push_back(point);
  // 4
  point << 0.4701420641, 0.4701420641, -0.7745696692;
  points_.push_back(point);
  // 5
  point << 0.7974269854, 0.1012865073, -0.7745696692;
  points_.push_back(point);
  // 6
  point << 0.1012865073, 0.7974269854, -0.7745696692;
  points_.push_back(point);
  // 7
  point << 0.1012865073, 0.1012865073, -0.7745696692;
  points_.push_back(point);
  // 8
  point << 0.3333333333, 0.3333333333, 0.0000000000;
  points_.push_back(point);
  // 9
  point << 0.0597158718, 0.4701420641, 0.0000000000;
  points_.push_back(point);
  // 10
  point << 0.4701420641, 0.0597158718, 0.0000000000;
  points_.push_back(point);
  // 11
  point << 0.4701420641, 0.4701420641, 0.0000000000;
  points_.push_back(point);
  // 12
  point << 0.7974269854, 0.1012865073, 0.0000000000;
  points_.push_back(point);
  // 13
  point << 0.1012865073, 0.7974269854, 0.0000000000;
  points_.push_back(point);
  // 14
  point << 0.1012865073, 0.1012865073, 0.0000000000;
  points_.push_back(point);
  // 15
  point << 0.3333333333, 0.3333333333, 0.7745696692;
  points_.push_back(point);
  // 16
  point << 0.0597158718, 0.4701420641, 0.7745696692;
  points_.push_back(point);
  // 17
  point << 0.4701420641, 0.0597158718, 0.7745696692;
  points_.push_back(point);
  // 18
  point << 0.4701420641, 0.4701420641, 0.7745696692;
  points_.push_back(point);
  // 19
  point << 0.7974269854, 0.1012865073, 0.7745696692;
  points_.push_back(point);
  // 20
  point << 0.1012865073, 0.7974269854, 0.7745696692;
  points_.push_back(point);
  // 21
  point << 0.1012865073, 0.1012865073, 0.7745696692;
  points_.push_back(point);

  // 0~6
  weights_.push_back(0.0625000000);
  for(int i = 0; i < 3; i++)
    weights_.push_back(0.0367761536);
  for(int i = 0; i < 3; i++)
    weights_.push_back(0.0349831057);

  // 7~13
  weights_.push_back(0.1);
  for(int i = 0; i < 3; i++)
    weights_.push_back(0.0588418457);
  for(int i = 0; i < 3; i++)
    weights_.push_back(0.0559729691);

  // 14~20
  weights_.push_back(0.0625000000);
  for(int i = 0; i < 3; i++)
    weights_.push_back(0.0367761536);
  for(int i = 0; i < 3; i++)
    weights_.push_back(0.0349831057);
}

Gauss8Cubic::Gauss8Cubic()
{
  Eigen::VectorXd point(24);
  point << -1.0 / sqrt(3.0), -1.0 / sqrt(3.0), -1.0 / sqrt(3.0), //
      1.0 / sqrt(3.0), -1.0 / sqrt(3.0), -1.0 / sqrt(3.0),       //
      1.0 / sqrt(3.0), 1.0 / sqrt(3.0), -1.0 / sqrt(3.0),        //
      -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), -1.0 / sqrt(3.0),       //
      -1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0),       //
      1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0),        //
      1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0),         //
      -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0);
  for(int i = 0; i < 8; i++)
  {
    Eigen::VectorXd tmp(3);
    for(int j = 0; j < 3; j++)
      tmp[j] = point(3 * i + j);
    points_.push_back(tmp);
  }

  for(int i = 0; i < 8; i++)
    weights_.push_back(1.0 * 1.0 * 1.0);
}

Gauss27Cubic::Gauss27Cubic()
{
  Eigen::VectorXd point(81);
  point << -sqrt(0.6), -sqrt(0.6), -sqrt(0.6), //
      0.0, -sqrt(0.6), -sqrt(0.6),             //
      sqrt(0.6), -sqrt(0.6), -sqrt(0.6),       //
      -sqrt(0.6), 0.0, -sqrt(0.6),             //
      0.0, 0.0, -sqrt(0.6),                    //
      sqrt(0.6), 0.0, -sqrt(0.6),              //
      -sqrt(0.6), sqrt(0.6), -sqrt(0.6),       //
      0.0, sqrt(0.6), -sqrt(0.6),              //
      sqrt(0.6), sqrt(0.6), -sqrt(0.6),        //
      -sqrt(0.6), -sqrt(0.6), 0.0,             //
      0.0, -sqrt(0.6), 0.0,                    //
      sqrt(0.6), -sqrt(0.6), 0.0,              //
      -sqrt(0.6), 0.0, 0.0,                    //
      0.0, 0.0, 0.0,                           //
      sqrt(0.6), 0.0, 0.0,                     //
      -sqrt(0.6), sqrt(0.6), 0.0,              //
      0.0, sqrt(0.6), 0.0,                     //
      sqrt(0.6), sqrt(0.6), 0.0,               //
      -sqrt(0.6), -sqrt(0.6), sqrt(0.6),       //
      0.0, -sqrt(0.6), sqrt(0.6),              //
      sqrt(0.6), -sqrt(0.6), sqrt(0.6),        //
      -sqrt(0.6), 0.0, sqrt(0.6),              //
      0.0, 0.0, sqrt(0.6),                     //
      sqrt(0.6), 0.0, sqrt(0.6),               //
      -sqrt(0.6), sqrt(0.6), sqrt(0.6),        //
      0.0, sqrt(0.6), sqrt(0.6),               //
      sqrt(0.6), sqrt(0.6), sqrt(0.6);
  for(int i = 0; i < 27; i++)
  {
    Eigen::VectorXd tmp(3);
    for(int j = 0; j < 3; j++)
      tmp[j] = point(3 * i + j);
    points_.push_back(tmp);
  }

  Eigen::VectorXd weight(81);
  weight << 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0,
      5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0,
      8.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0,
      5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0,
      5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0,
      8.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0,
      5.0 / 9.0, 8.0 / 9.0, 8.0 / 9.0, 8.0 / 9.0, 8.0 / 9.0, 8.0 / 9.0,
      5.0 / 9.0, 8.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0,
      8.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0,
      5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0,
      5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0,
      8.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0,
      5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0,
      5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0;
  for(int i = 0; i < 27; i++)
  {
    double tmp = 1.0;
    for(int j = 0; j < 3; j++)
      tmp *= weight(3 * i + j);
    weights_.push_back(tmp);
  }
}

} // namespace icarat
