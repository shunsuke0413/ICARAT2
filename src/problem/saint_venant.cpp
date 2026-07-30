///  @file  saint_venant.cpp
///  @author  Daiki Watanabe
///  @date  June 4, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "saint_venant.hpp"
#include <Eigen/LU>
#include <basis/bmatrix.hpp>
#include <basis/bmatrix_total.hpp>
#include <basis/bmatrix_updated.hpp>
#include <misc/assembling.hpp>
#include <misc/eigen_manipulation.hpp>
#include <unsupported/Eigen/KroneckerProduct>

namespace icarat
{
using namespace std;
using namespace Eigen;

void SaintVenantTotal::totalKirchhoff(FEM &fem, VectorXd &sigmaVector,
                                      MatrixXd &F, MatrixXd &sigmaMatrix)
{
  MatrixXd I = MatrixXd::Identity(fem.ndim, fem.ndim);

  // Cauchy Greeen tensor E
  MatrixXd E = 0.5 * (F.transpose() * F - I);
  double jac = F.determinant();

  if(fem.ndim == 2)
  {
    // convert to voigt expression
    VectorXd Ev(fem.voigt);
    Ev << E(0, 0), E(1, 1), E(0, 1) + E(1, 0);

    // Create elasticity tensor.
    double mu0 = 0.5 * young_ / (1 + poisson_);
    double lambda0 = 2.0 * mu0 * poisson_ / (1.0 - 2.0 * poisson_);

    lameDe(De_, mu0, lambda0);

    // Saint-Venent Kirchhoff
    sigmaVector = De_ * Ev;

    // make 2nd PK stress (bonet book(5.43))
    MatrixXd PK2(fem.ndim, fem.ndim);
    PK2 << sigmaVector(0), sigmaVector(2), //
        sigmaVector(2), sigmaVector(1);

    // matrix expression
    sigmaMatrix = kroneckerProduct(PK2, I);
  }
  else
  {
    // convert to voigt expression
    VectorXd Ev(fem.voigt);
    Ev << E(0, 0), E(1, 1), E(2, 2), E(0, 1) + E(1, 0), E(1, 2) + E(2, 1),
        E(2, 0) + E(0, 2);

    // Create elasticity tensor.
    young_ = SIMP(young1_, young2_, design_s_, pp_);
    double mu0 = 0.5 * young_ / (1 + poisson_);
    double lambda0 = 2.0 * mu0 * poisson_ / (1.0 - 2.0 * poisson_);

    lameDe(De_, mu0, lambda0);

    // Saint-Venent Kirchhoff
    sigmaVector = De_ * Ev;

    // make 9*9 matrix (refer bonet book(5.43))
    MatrixXd PK2(fem.ndim, fem.ndim);

    PK2 << sigmaVector(0), sigmaVector(3), sigmaVector(5), //
        sigmaVector(3), sigmaVector(1), sigmaVector(4),    //
        sigmaVector(5), sigmaVector(4), sigmaVector(2);

    // matrix expression
    sigmaMatrix = kroneckerProduct(PK2, I);
  }
}

void SaintVenantTotal::makeMisesStress(FEM &fem, vector<Node> &node)
{
  MatrixXd V = MatrixXd::Zero(fem.voigt, fem.voigt);

  if(fem.ndim == 2)
    V << 1.0, -0.5, 0.0, //
        -0.5, 1.0, 0.0,  //
        0.0, 0.0, 3.0;
  else if(fem.ndim == 3)
    V << 1.0, -0.5, -0.5, 0.0, 0.0, 0.0, //
        -0.5, 1.0, -0.5, 0.0, 0.0, 0.0,  //
        -0.5, -0.5, 1.0, 0.0, 0.0, 0.0,  //
        0.0, 0.0, 0.0, 3.0, 0.0, 0.0,    //
        0.0, 0.0, 0.0, 0.0, 3.0, 0.0,    //
        0.0, 0.0, 0.0, 0.0, 0.0, 3.0;

  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd U = MatrixXd::Zero(ne, fem.ndim);
  VectorXd Uarr = VectorXd::Zero(ne * fem.ndim);

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
    {
      X.coeffRef(i, j) = node[nodeID[i]].x[j];
      U.coeffRef(i, j) = node[nodeID[i]].val[j];
      Uarr.coeffRef(fem.ndim * i + j) = node[nodeID[i]].val[j];
    }
  }

  // init
  mises_ = 0.0;
  strain_.setZero();
  stress_.setZero();
  BmatrixTotal bTotal;
  // gauss point loop start
  for(int ip = 0; ip < ipmax; ip++)
  {
    // nonlinear part
    double jac = 1.0;
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
    MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
    // second piola-Kirchhoff stress tensor
    MatrixXd Sbar = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);

    // make bmatrix for finite strain
    bTotal.make(BeL, BeNL, Ne, jac, F, eType, ipmax, X, U, ip);
    assert(jac > 0.0);

    totalKirchhoff(fem, sigmaVector, F, Sbar);
    VectorXd PK2 = De_ * BeL * Uarr;
    MatrixXd cauchy = F * Tensor(PK2, 1) * F.transpose() / F.determinant();
    strain_ += BeL * Uarr / ipmax;
    stress_ += Voigt(cauchy, 1) / ipmax;
    mises_ += sqrt(Voigt(cauchy, 1).transpose() * V * Voigt(cauchy, 1)) / ipmax;
    // mises_ += sqrt(PK2.transpose() * V * PK2) / ipmax;
  }
}

void FiniteSaintVenantTotal::makeKeFint(vector<Triplet<double>> &tripletsK,
                                        vector<Triplet<double>> &tripletsF,
                                        vector<Triplet<double>> &tripletsR,
                                        FEM &fem, vector<Node> &nodes)
{
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd U = MatrixXd::Zero(ne, fem.ndim);
  VectorXd disp = VectorXd::Zero(numdof);
  VectorXi idof = VectorXi::Zero(numdof);

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
    {
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
      U.coeffRef(i, j) = gamma_ * nodes[nodeID[i]].val[j];
      disp.coeffRef(fem.ndim * i + j) = nodes[nodeID[i]].val[j];
      idof.coeffRef(fem.ndim * i + j) = nodes[nodeID[i]].dof[j];
    }
  }

  volume = 0.0;
  Ke_.setZero();
  fe_.setZero();

  // linear init
  l_.init();
  Bmatrix bLinear;

  // nonlinear init
  MatrixXd nKe = Ke_;
  VectorXd nfe = fe_;
  BmatrixTotal bTotal;

  // gauss point loop start
  for(int ip = 0; ip < ipmax; ip++)
  {
    young_ = SIMP(young1_, young2_, design_s_, pp_);

    // linear part
    double jac = 1.0;
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
    MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
    De(l_.De, young_, poisson_);
    bLinear.make(Be, Ne, jac, eType, ipmax, X, ip);
    l_.strain = Be * disp;
    l_.stress = l_.De * l_.strain;
    l_.fe += Be.transpose() * l_.stress * jac;
    l_.Ke += Be.transpose() * l_.De * Be * jac;
    volume += jac;

    // nonlinear part
    jac = 1.0;
    MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
    // second piola-Kirchhoff stress tensor
    MatrixXd Sbar = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);

    // make bmatrix for finite strain
    bTotal.make(BeL, BeNL, Ne, jac, F, eType, ipmax, X, U, ip);
    assert(jac > 0.0);

    totalKirchhoff(fem, sigmaVector, F, Sbar);

    // inner stress(Bonet's book (9.15b))
    nfe += BeL.transpose() * sigmaVector * jac;

    // element stiffness matrix(Bonet's book (9.31),(9.35),(9.44c))
    nKe += BeL.transpose() * De_ * BeL * jac +
           BeNL.transpose() * Sbar * BeNL * jac;
  }

  // element stiffness matrix(Bonet's book (9.31),(9.35),(9.44c))
  Ke_ = (1.0 - gamma_ * gamma_) * l_.Ke + gamma_ * gamma_ * nKe;

  // inner stress(Bonet's book (9.15b))
  fe_ = (1.0 - gamma_ * gamma_) * l_.fe + gamma_ * nfe;

  lfe_ = l_.fe;
  nfe_ = nfe;

  assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
  assemblingReact(numdof, idof, fe_, tripletsR);
}

void BucklingFiniteSaintVenantTotal::makeEigenCoeffs(
    std::vector<Eigen::Triplet<double>> &tripletsA,
    std::vector<Eigen::Triplet<double>> &tripletsB, FEM &fem,
    std::vector<Node> &nodes)
{
  // get stress
  makeMisesStress(fem, nodes);

  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
  MatrixXd Be = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
  MatrixXd Se = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  VectorXd disp = VectorXd::Zero(numdof);
  VectorXi idof = VectorXi::Zero(numdof);

  if(fem.ndim == 2)
    for(int i = 0; i < 2; i++)
    {
      Se.coeffRef(2 * i, 2 * i) = stress_.coeff(0);
      Se.coeffRef(2 * i, 2 * i + 1) = stress_.coeff(2);
      Se.coeffRef(2 * i + 1, 2 * i) = stress_.coeff(2);
      Se.coeffRef(2 * i + 1, 2 * i + 1) = stress_.coeff(1);
    }
  else
  {
    for(int i = 0; i < 3; i++)
    {
      Se.coeffRef(3 * i, 3 * i) = stress_.coeff(0);
      Se.coeffRef(3 * i, 3 * i + 1) = stress_.coeff(3);
      Se.coeffRef(3 * i, 3 * i + 2) = stress_.coeff(5);
      Se.coeffRef(3 * i + 1, 3 * i) = stress_.coeff(3);
      Se.coeffRef(3 * i + 1, 3 * i + 1) = stress_.coeff(1);
      Se.coeffRef(3 * i + 1, 3 * i + 2) = stress_.coeff(4);
      Se.coeffRef(3 * i + 2, 3 * i) = stress_.coeff(5);
      Se.coeffRef(3 * i + 2, 3 * i + 1) = stress_.coeff(4);
      Se.coeffRef(3 * i + 2, 3 * i + 2) = stress_.coeff(2);
    }
  }

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.dofnp; j++)
    {
      disp.coeffRef(fem.dofnp * i + j) = nodes[nodeID[i]].val[j];
      idof.coeffRef(fem.dofnp * i + j) = nodes[nodeID[i]].dof[j];
    }
    for(int j = 0; j < fem.ndim; j++)
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
  }

  Bmatrix bmatrix;
  double jac = 0.0;
  KGe_.setZero();
  for(int ip = 0; ip < ipmax; ip++)
  {
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);
    KGe_ += Be.transpose() * Se * Be * jac;
  }

  assembling(numdof, fem.numeq, idof, Ke_, tripletsA);
  assembling(numdof, fem.numeq, idof, KGe_, tripletsB);
}

void mHomoSaintVenantTotal::makeKeHomo(vector<Triplet<double>> &tripletsK,
                                       vector<Triplet<double>> &tripletsF,
                                       FEM &fem, vector<HomoNode> &nodes,
                                       int dir)
{
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd U = MatrixXd::Zero(ne, fem.ndim);
  VectorXi idof = VectorXi::Zero(numdof);

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.dofnp; j++)
    {
      idof.coeffRef(fem.dofnp * i + j) = nodes[nodeID[i]].dof[j];
      U.coeffRef(i, j) = nodes[nodeID[i]].nmtval[dir][j];
    }

    for(int j = 0; j < fem.ndim; j++)
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
  }

  young_ = SIMP(young1_, young2_, design_s_, pp_);
  De(De_, young_, poisson_);

  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;
  BmatrixTotal bmatrix;
  for(int ip = 0; ip < ipmax; ip++)
  {
    double jac = 1.0;
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
    MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
    // second piola-Kirchhoff stress tensor
    MatrixXd Sbar = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    // stess for internal force
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);

    // make bmatrix for finite strain
    bmatrix.make(BeL, BeNL, Ne, jac, F, eType, ipmax, X, U, ip);

    totalKirchhoff(fem, sigmaVector, F, Sbar);

    // inner stress(Bonet's book (9.15b))
    fe_ += BeL.transpose() * sigmaVector * jac;

    // element stiffness matrix(Bonet's book (9.31),(9.35),(9.44c))
    Ke_ += BeL.transpose() * De_ * BeL * jac +
           BeNL.transpose() * Sbar * BeNL * jac;

    volume += jac;
  }

  assemblingHomo(ne, numdof, idof, Ke_, fe_, tripletsK, tripletsF, fem, nodeID,
                 nodes);
}

void mHomoSaintVenantTotal::makeMicroValues(FEM &fem, vector<HomoNode> &nodes)
{
  VectorXi idof = VectorXi::Zero(numdof);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd U = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd disps = MatrixXd::Zero(numdof, fem.voigt); ///< for nmt

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.dofnp; j++)
    {
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
      idof.coeffRef(i * fem.dofnp + j) = nodes[nodeID[i]].dof[j];

      for(int k = 0; k < fem.voigt; k++)
        disps(fem.dofnp * i + j, k) = nodes[nodeID[i]].nmtval[k][j];
    }
  }

  sEnergy_.setZero();
  strain_.setZero();
  stress_.setZero();
  volume = 0.0;
  double fac = (double)(1.0 / ipmax);
  BmatrixTotal bmatrix;
  for(int ip = 0; ip < ipmax; ip++)
  {
    double jac = 1.0;
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
    MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
    // second piola-Kirchhoff stress tensor
    MatrixXd Sbar = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    // stess for internal force
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);

    bmatrix.make(BeL, BeNL, Ne, jac, F, eType, ipmax, X, U, ip);

    totalKirchhoff(fem, sigmaVector, F, Sbar);

    MatrixXd mstrain = BeL * disps;
    volume += jac;
    // strain_ += BeL * disp * fac;
    //  stress_ += De_ * Be * disp * fac;

    sEnergy_ += mstrain.transpose() * De_ * mstrain * jac;
  }

  // make mises stress
  // MatrixXd V = MatrixXd::Zero(fem.voigt, fem.voigt);
  // if(fem.ndim == 2)
  //   V << 1.0, -0.5, 0.0, //
  //       -0.5, 1.0, 0.0,  //
  //       0.0, 0.0, 3.0;
  // else if(fem.ndim == 3)
  //   V << 1.0, -0.5, -0.5, 0.0, 0.0, 0.0, //
  //       -0.5, 1.0, -0.5, 0.0, 0.0, 0.0,  //
  //       -0.5, -0.5, 1.0, 0.0, 0.0, 0.0,  //
  //       0.0, 0.0, 0.0, 3.0, 0.0, 0.0,    //
  //       0.0, 0.0, 0.0, 0.0, 3.0, 0.0,    //
  //       0.0, 0.0, 0.0, 0.0, 0.0, 3.0;

  // mises_ = sqrt(stress_.transpose() * V * stress_);
}

} // namespace icarat