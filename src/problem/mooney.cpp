///  @file  mooney.cpp
///  @author  Daiki Watanabe
///  @date  October 3, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "mooney.hpp"
#include <Eigen/LU>
#include <basis/bmatrix_total.hpp>
#include <iostream>
#include <misc/assembling.hpp>
#include <unsupported/Eigen/KroneckerProduct>

namespace icarat
{
using namespace std;
using namespace Eigen;

void MooneyTotal::setParameter(const toml::value &config)
{
  const auto &mat = toml::find(config, "material");
  assert(toml::find<string>(mat, "type") == "Mooney");

  design_s_ = toml::find<double>(mat, "designs_0");
  pp_ = toml::find<double>(mat, "penalty");

  C10s_ = toml::find<vector<double>>(mat, "C10");
  C01s_ = toml::find<vector<double>>(mat, "C01");
  D1s_ = toml::find<vector<double>>(mat, "D01");

  mattype_ = toml::find<string>(mat, "state");
  if(mattype_ == "plane_stress" || mattype_ == "plane_strain")
    cerr << "This state is invalid in mooney.cpp" << endl;
  else
    mattype_ = "3D";
  if(toml::find<int>(config, "mesh", "dimension") == 3)
    mattype_ = "3D";
}

void MooneyTotal::totalKirchhoff(FEM &fem, VectorXd &sigma, MatrixXd &F,
                                 MatrixXd &sigmaMatrix)
{
  // set SIMP coefficient
  SIMPMooney();

  // right Cauchy Greeen tensor
  MatrixXd C = F.transpose() * F;

  double C1 = C(0, 0);
  double C2 = C(1, 1);
  double C3 = C(2, 2);
  double C4 = C(0, 1);
  double C5 = C(1, 2);
  double C6 = C(0, 2);
  double I1 = C1 + C2 + C3;
  double I2 = C1 * C2 + C1 * C3 + C2 * C3 - C4 * C4 - C5 * C5 - C6 * C6;
  double I3 = C.determinant();
  if(I3 < 0.0)
  {
    cout << "jacobian is negative: " << I3 << endl;
    exit(1);
  }
  double J3 = sqrt(I3);
  double J3M1 = J3 - 1.0;

  VectorXd I1E(6), I2E(6), I3E(6);
  I1E << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;
  I1E *= 2.0;
  I2E << C2 + C3, C3 + C1, C1 + C2, -C4, -C5, -C6;
  I2E *= 2.0;
  I3E << C2 * C3 - C5 * C5, C3 * C1 - C6 * C6, C1 * C2 - C4 * C4,
      C5 * C6 - C3 * C4, C6 * C4 - C1 * C5, C4 * C5 - C2 * C6;
  I3E *= 2.0;

  double X12 = 1.0 / 2.0;
  double X13 = 1.0 / 3.0;
  double X23 = 2.0 / 3.0;
  double X43 = 4.0 / 3.0;
  double X53 = 5.0 / 3.0;
  double X89 = 8.0 / 9.0;

  double W1 = pow(I3, -X13);
  double W2 = X13 * I1 * pow(I3, -X43);
  double W3 = pow(I3, -X23);
  double W4 = X23 * I2 * pow(I3, -X53);
  double W5 = X12 * pow(I3, -X12);

  VectorXd J1E = W1 * I1E - W2 * I3E;
  VectorXd J2E = W3 * I2E - W4 * I3E;
  VectorXd J3E = W5 * I3E;

  // voigt expression
  sigma = C10_ * J1E + C01_ * J2E + D1_ * J3M1 * J3E;

  // matrix expression
  MatrixXd eye = MatrixXd::Identity(fem.ndim, fem.ndim);
  MatrixXd sig(3, 3);
  sig << sigma(0), sigma(3), sigma(5), //
      sigma(3), sigma(1), sigma(4),    //
      sigma(5), sigma(4), sigma(2);

  sigmaMatrix = kroneckerProduct(eye, sig);

  //////////////////////////////
  /// make De matrix
  //////////////////////////////
  MatrixXd I2EE(fem.voigt, fem.voigt), I3EE(fem.voigt, fem.voigt);
  I2EE << 0.0, 4.0, 4.0, 0.0, 0.0, 0.0, //
      4.0, 0.0, 4.0, 0.0, 0.0, 0.0,     //
      4.0, 4.0, 0.0, 0.0, 0.0, 0.0,     //
      0.0, 0.0, 0.0, -2.0, 0.0, 0.0,    //
      0.0, 0.0, 0.0, 0.0, -2.0, 0.0,    //
      0.0, 0.0, 0.0, 0.0, 0.0, -2.0;

  I3EE << 0.0, 4.0 * C3, 4.0 * C2, 0.0, -4.0 * C5, 0.0,   //
      4.0 * C3, 0.0, 4.0 * C1, 0.0, 0.0, -4.0 * C6,       //
      4.0 * C2, 4.0 * C1, 0.0, -4.0 * C4, 0.0, 0.0,       //
      0.0, 0.0, -4.0 * C4, -2.0 * C3, 2.0 * C6, 2.0 * C5, //
      -4.0 * C5, 0.0, 0.0, 2.0 * C6, -2.0 * C1, 2.0 * C4, //
      0.0, -4.0 * C6, 0.0, 2.0 * C5, 2.0 * C4, -2.0 * C2;

  W1 = X23 * pow(I3, -X12);
  W2 = X89 * I1 * pow(I3, -X43);
  W3 = X13 * I1 * pow(I3, -X43);
  W4 = X43 * pow(I3, -X12);
  W5 = X89 * I2 * pow(I3, -X53);
  double W6 = pow(I3, -X23);
  double W7 = X23 * I2 * pow(I3, -X53);
  double W8 = pow(I3, -X12);
  double W9 = X12 * pow(I3, -X12);

  MatrixXd J1EE = -W1 * (J1E * J3E.transpose() + J3E * J1E.transpose()) +
                  W2 * (J3E * J3E.transpose()) - W3 * I3EE;
  MatrixXd J2EE = -W4 * (J2E * J3E.transpose() + J3E * J2E.transpose()) +
                  W5 * (J3E * J3E.transpose()) + W6 * I2EE - W7 * I3EE;
  MatrixXd J3EE = -W8 * (J3E * J3E.transpose()) + W9 * I3EE;

  De_.setZero();
  De_ = C10_ * J1EE + C01_ * J2EE + D1_ * (J3E * J3E.transpose()) +
        D1_ * J3M1 * J3EE;
}

void MooneyTotal::SIMPMooney()
{
  if(C10s_[1] > C10s_[0] && C01s_[1] > C01s_[0] && D1s_[1] > D1s_[0])
  {
    C10_ = pow(design_s_, pp_) * C10s_[1] + (1.0 - design_s_) * C10s_[0];
    C01_ = pow(design_s_, pp_) * C01s_[1] + (1.0 - design_s_) * C01s_[0];
    D1_ = pow(design_s_, pp_) * D1s_[1] + (1.0 - design_s_) * D1s_[0];
  }
  else
  {
    cerr << "All parameters of material 2 must be greater than 1." << endl;
    exit(1);
  }
}

void FiniteMooneyTotal::makeKeFint(vector<Triplet<double>> &tripletsK,
                                   vector<Triplet<double>> &tripletsF,
                                   vector<Triplet<double>> &tripletsR, FEM &fem,
                                   vector<Node> &nodes)
{
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd U = MatrixXd::Zero(ne, fem.ndim);
  VectorXi idof = VectorXi::Zero(numdof);

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
    {
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
      U.coeffRef(i, j) = nodes[nodeID[i]].val[j];
      idof.coeffRef(fem.ndim * i + j) = nodes[nodeID[i]].dof[j];
    }
  }

  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;
  BmatrixTotal bmatrix;

  // gauss point loop start
  for(int ip = 0; ip < ipmax; ip++)
  {
    double jac = 1.0;
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
    MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
    // second piola-Kirchhoff stress tensor
    MatrixXd Sbar = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    VectorXd sigma = VectorXd::Zero(fem.voigt);

    // make bmatrix for finite strain
    bmatrix.make(BeL, BeNL, Ne, jac, F, eType, ipmax, X, U, ip);

    totalKirchhoff(fem, sigma, F, Sbar);

    // inner stress
    fe_ += BeL.transpose() * sigma * jac;

    // refer bonet's book formula (8.10)
    Ke_ += BeL.transpose() * De_ * BeL * jac;
    Ke_ += BeNL.transpose() * Sbar * BeNL * jac;

    volume += jac;
  } // gauss point loop end

  assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
  assemblingReact(numdof, idof, fe_, tripletsR);
}
} // namespace icarat