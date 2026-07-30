///  @file  elast_levelset.cpp
///  @author  Daiki Watanabe
///  @date  October 13, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "elast_levelset.hpp"
#include <basis/bmatrix.hpp>
#include <basis/bmatrix_shell.hpp>
#include <iostream>
#include <misc/assembling.hpp>

namespace icarat
{
namespace opt_levelset
{
using namespace std;
using namespace Eigen;

void ElastLevelset::makeKeFint(vector<Triplet<double>> &tripletsK,
                               vector<Triplet<double>> &tripletsF,
                               vector<Triplet<double>> &tripletsR, FEM &fem,
                               vector<NodeLS> &node)
{
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
  MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  VectorXd disp = VectorXd::Zero(numdof);
  VectorXi idof = VectorXi::Zero(numdof);

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.dofnp; j++)
    {
      disp.coeffRef(fem.dofnp * i + j) = node[nodeID[i]].val[j];
      idof.coeffRef(fem.dofnp * i + j) = node[nodeID[i]].dof[j];
    }
    for(int j = 0; j < fem.ndim; j++)
      X.coeffRef(i, j) = node[nodeID[i]].x[j];
  }

  young_ = (1.0 - design_s_) * young1_ + design_s_ * young2_;

  /// make elastic coefficient tensor
  De(De_, young_, poisson_);

  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;
  double jac = 0.0;
  Bmatrix bmatrix;

  for(int ip = 0; ip < ipmax; ip++)
  {
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);

    /// calculation
    fe_ += Be.transpose() * De_ * Be * disp * jac;
    Ke_ += Be.transpose() * De_ * Be * jac;
    volume += jac;
  }

  assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
}

void ElastLevelset::makeMisesStress(FEMLS &fem, vector<NodeLS> &node)
{
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
  MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  VectorXd disp = VectorXd::Zero(numdof);
  MatrixXd V = MatrixXd::Zero(fem.voigt, fem.voigt);

  if(fem.ndim == 2)
    V << 1, -0.5, 0, -0.5, 1, 0, 0, 0, 3;
  else if(fem.ndim == 3)
    V << 1, -0.5, -0.5, 0, 0, 0, -0.5, 1, -0.5, 0, 0, 0, -0.5, -0.5, 1, 0, 0, 0,
        0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 3;

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.dofnp; j++)
    {
      disp.coeffRef(fem.dofnp * i + j) = node[nodeID[i]].val[j];
      X.coeffRef(i, j) = node[nodeID[i]].x[j];
    }
  }

  strain_.setZero();
  stress_.setZero();
  /// gauss point loop start
  double fac = (double)(1.0 / ipmax);
  double dummy;
  Bmatrix bmatrix;
  for(int ip = 0; ip < ipmax; ip++)
  {
    bmatrix.make(Be, Ne, dummy, eType, ipmax, X, ip);
    strain_ += Be * disp * fac;
    stress_ += De_ * Be * disp * fac;
  }

  mises_ = sqrt(stress_.transpose() * V * stress_);
}
} // namespace opt_levelset
} // namespace icarat