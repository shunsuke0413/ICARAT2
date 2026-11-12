///  @file  neo_hooke.cpp
///  @author  Daiki Watanabe
///  @date  June 4, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "neo_hooke_damage.hpp"
#include <Eigen/LU>
#include <basis/bmatrix.hpp>
#include <basis/bmatrix_total.hpp>
#include <basis/bmatrix_updated.hpp>
#include <iostream>
#include <misc/assembling.hpp>
#include <unsupported/Eigen/KroneckerProduct>
namespace icarat
{
using namespace std;
using namespace Eigen;

double FiniteNeoHookedamage::computeTau(FEM &fem, std::vector<Node> &nodes)
{
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd U = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd Uorg = MatrixXd::Zero(ne, fem.ndim);
  VectorXd disp = VectorXd::Zero(numdof);
  VectorXi idof = VectorXi::Zero(numdof);
  Gauss9Square gauss;

  young_ = SIMP(young1_, young2_, design_s_, pp_);
  double mu = 0.5 * young_ / (1 + poisson_);
  double lambda =
      (poisson_ * young_) / ((1.0 + poisson_) * (1.0 - 2.0 * poisson_)) +
      2.0 / 3.0 * mu; // 本当はκ、とりあえず名前変えていないだけ

  mu_ = mu;
  lambda_ = lambda;

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
    {
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j]; // place of element's node point
      U.coeffRef(i, j) = gamma_ * nodes[nodeID[i]].val[j];
      disp.coeffRef(fem.ndim * i + j) = nodes[nodeID[i]].val[j];
      idof.coeffRef(fem.ndim * i + j) =
          nodes[nodeID[i]].dof[j]; // global dof number
    }
  }

  VectorXd nfe = fe_;
  BmatrixTotal bTotal;

  double jac = 1.0;

  // linear part
  double J = 1.0;
  double J_m23 = 1.0;
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne); // shape function
  MatrixXd Be =
      MatrixXd::Zero(fem.voigt, numdof); // strain-displacement matrix ε = B u

  // nonlinear part
  MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
  MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
  MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);

  MatrixXd C = MatrixXd::Zero(fem.ndim, fem.ndim); // 右コーシー

  double tau_ave = 0.0;

  // gauss point loop start
  for(int ip = 0; ip < ipmax; ip++)
  {
    young_ = SIMP(young1_, young2_, design_s_, pp_);
    double w = gauss.weights(ip);

    // make bmatrix for finite strain
    bTotal.make(BeL, BeNL, Ne, jac, F, eType, ipmax, X, U, ip);
    assert(jac > 0.0);

    J = F.determinant();
    C = F.transpose() * F;
    // Plane strain embeds the 2D deformation in 3D with Fzz = Czz = 1.
    traC_ = C.trace() + (fem.ndim == 2 ? 1.0 : 0.0);
    J_m23 = std::pow(J, -2.0 / 3.0);
    double Psi_iso = 0.5 * mu * (traC_ * J_m23 - 3.0);
    double Psi_vol = 0.5 * lambda * (J - 1.0) * (J - 1.0);
    double Psi = Psi_iso + Psi_vol;

    double tau = sqrt(2.0 * std::abs(Psi));

    tau_ave += tau / ipmax;

    tau_all_[ip] = tau;
  }
  return tau_ave;
}

double neohookedamage::calculateTensileDamage(double tau)
{
  return evaluateDamageLaw(tau).trial;
}

double neohookedamage::differentiateD(double tau)
{
  return evaluateDamageLaw(tau).derivative;
}

//----------------------- for S:E integration ------------------//
void neohookedamage::totalKirchhoff(FEM &fem, MatrixXd &PK2,
                                    VectorXd &sigmaVector, MatrixXd &F,
                                    MatrixXd &sigmaMatrix, double &product,
                                    Eigen::MatrixXd &Ctan, int ip,
                                    double tauFiltered)
//--------------------------------------------------------------//
{
  // lame constant for initial & current position　ラメ定数p.194
  double mu = 0.5 * young_ / (1 + poisson_);
  double lambda =
      (poisson_ * young_) / ((1.0 + poisson_) * (1.0 - 2.0 * poisson_)) +
      2.0 / 3.0 * mu; // 本当はκ、とりあえず名前変えていないだけ

  mu_ = mu;
  lambda_ = lambda;

  MatrixXd I = MatrixXd::Identity(fem.ndim, fem.ndim); // 単位行列
  MatrixXd C = MatrixXd::Zero(fem.ndim, fem.ndim);     // 右コーシー
  MatrixXd E = MatrixXd::Zero(fem.ndim, fem.ndim);     // グリーン・ラグランジュ

  // right C-G tensor
  C = F.transpose() * F;
  traC_ = C.trace();
  // Green strain　ラグランジュ・グリーンひずみテンソル
  E = 0.5 * (C - I);

  double J = F.determinant();
  double J_m23 = std::pow(J, -2.0 / 3.0);

  // ガウス点ごとのリストを更新
  double D_trial = calculateTensileDamage(tauFiltered);

  D_all_[ip] = std::max(D_all_[ip], D_trial);

  // 今後の計算で使う現在の損傷度
  double current_D = D_all_[ip];
  // double current_D = 0.0;

  if(fem.ndim == 2)
  {
    if(mattype_ != "plane_strain")
      throw std::runtime_error("2D Neo-Hookean damage currently supports plane_strain only.");
    if(fem.voigt != 3 || F.rows() != 2 || F.cols() != 2 ||
       !F.allFinite() || !std::isfinite(J) || !(J > 0.0))
      throw std::runtime_error("Invalid 2D deformation in Neo-Hookean damage.");

    // Restrict the same 3D isochoric/volumetric law to Fzz = 1.
    // Retain J^(-2/3) and the full 3D trace; using trace(C2D) alone
    // would produce nonzero stress and energy in the undeformed state.
    const double trace3D = C.trace() + 1.0;
    const MatrixXd C_inv = C.inverse();
    const double retained = 1.0 - current_D;
    const double volumetric = lambda * J * (J - 1.0); // lambda denotes bulk modulus.
    PK2 = retained * (mu * J_m23 * (I - trace3D / 3.0 * C_inv) +
                      volumetric * C_inv);
    sigmaVector.resize(3);
    sigmaVector << PK2(0, 0), PK2(1, 1), PK2(0, 1);
    sigmaMatrix = kroneckerProduct(PK2, I);

    // Ctan maps [dExx, dEyy, 2*dExy] to [dSxx, dSyy, dSxy].
    // For update B, filtered tau and damage are fixed during this NR solve.
    // The derivative of the damage law belongs to the history coupling,
    // not to this partial derivative at fixed damage.
    const int vi[3] = {0, 1, 0};
    const int vj[3] = {0, 1, 1};
    Ctan = MatrixXd::Zero(3, 3);
    for(int row = 0; row < 3; ++row)
      for(int col = 0; col < 3; ++col)
      {
        const int i = vi[row], j = vj[row], k = vi[col], l = vj[col];
        const double inverseSym = C_inv(i, k) * C_inv(j, l) +
                                  C_inv(i, l) * C_inv(j, k);
        const double vol = lambda * J * (2.0 * J - 1.0) * C_inv(i, j) * C_inv(k, l)
                           - volumetric * inverseSym;
        const double iso1 = -C_inv(k, l) / 3.0 *
                               (I(i, j) - trace3D / 3.0 * C_inv(i, j));
        const double iso2 = -I(k, l) / 3.0 * C_inv(i, j) +
                           trace3D / 6.0 * inverseSym;
        Ctan(row, col) = retained * (vol + 2.0 * mu * J_m23 * (iso1 + iso2));
      }

    product += PK2(0, 0) * E(0, 0) + PK2(1, 1) * E(1, 1) +
               2.0 * PK2(0, 1) * E(0, 1);

    // Plane strain has nonzero out-of-plane stress; include it in 3D von Mises.
    const double Szz = retained * (mu * J_m23 * (1.0 - trace3D / 3.0) + volumetric);
    const MatrixXd cauchy = F * PK2 * F.transpose() / J;
    const double sxx = cauchy(0, 0), syy = cauchy(1, 1), szz = Szz / J;
    mises_ = std::sqrt(0.5 * ((sxx - syy) * (sxx - syy) +
                              (syy - szz) * (syy - szz) +
                              (szz - sxx) * (szz - sxx)) +
                       3.0 * cauchy(0, 1) * cauchy(0, 1));
  }
  else // 3D
  {

    // tensol -> vector
    VectorXd strainvec = VectorXd::Zero(fem.voigt);
    strainvec[0] = E(0, 0);
    strainvec[1] = E(1, 1);
    strainvec[2] = E(2, 2);
    strainvec[3] = 2.0 * E(0, 1);
    strainvec[4] = 2.0 * E(1, 2);
    strainvec[5] = 2.0 * E(2, 0);
    // cal stress S
    double jac = F.determinant();
    double jacC = pow(jac, 2.0);
    MatrixXd C_inv = MatrixXd::Zero(fem.ndim, fem.ndim);
    C_inv(0, 0) = (C(1, 1) * C(2, 2) - C(1, 2) * C(2, 1)) / jacC;
    C_inv(0, 1) = (C(0, 2) * C(2, 1) - C(0, 1) * C(2, 2)) / jacC;
    C_inv(0, 2) = (C(0, 1) * C(1, 2) - C(0, 2) * C(1, 1)) / jacC;

    C_inv(1, 0) = (C(1, 2) * C(2, 0) - C(1, 0) * C(2, 2)) / jacC;
    C_inv(1, 1) = (C(0, 0) * C(2, 2) - C(0, 2) * C(2, 0)) / jacC;
    C_inv(1, 2) = (C(0, 2) * C(1, 0) - C(0, 0) * C(1, 2)) / jacC;

    C_inv(2, 0) = (C(1, 0) * C(2, 1) - C(1, 1) * C(2, 0)) / jacC;
    C_inv(2, 1) = (C(0, 1) * C(2, 0) - C(0, 0) * C(2, 1)) / jacC;
    C_inv(2, 2) = (C(0, 0) * C(1, 1) - C(0, 1) * C(1, 0)) / jacC;

    MatrixXd PK2 = MatrixXd::Zero(fem.ndim, fem.ndim);
    MatrixXd PK2_original = MatrixXd::Zero(fem.ndim, fem.ndim);
    MatrixXd PK2_iso = MatrixXd::Zero(fem.ndim, fem.ndim);
    MatrixXd PK2_vol = MatrixXd::Zero(fem.ndim, fem.ndim);

    for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
      {
        PK2_iso(i, j) =
            mu * J_m23 * (I(i, j) - (1.0 / 3.0) * traC_ * C_inv(i, j));
        PK2_vol(i, j) = lambda * jac * (jac - 1.0) * C_inv(i, j);

        PK2_original(i, j) = PK2_iso(i, j) + PK2_vol(i, j);
        PK2(i, j) = PK2_original(i, j) * (1.0 - current_D);
      }
    }

    sigmaVector(0) = PK2(0, 0);
    sigmaVector(1) = PK2(1, 1);
    sigmaVector(2) = PK2(2, 2);
    sigmaVector(3) = PK2(0, 1);
    sigmaVector(4) = PK2(1, 2);
    sigmaVector(5) = PK2(2, 0);

    sigmaMatrix = kroneckerProduct(PK2, I);

    // De matrix
    // make Voigt-Tensor index data
    MatrixXi IJdata = MatrixXi::Zero(fem.ndim, fem.ndim);
    IJdata << 0, 3, 5, // xx, xy, xz
        3, 1, 4,       // xy, yy, yz
        5, 4, 2;       // xz, yz, zz

    // tauFiltered comes from the previous converged load step and is fixed
    // during this NR solve. Thus dD/du_current = 0. The history adjoint
    // differentiates D(tauFiltered) with the shared sigmoid derivative.

    double p = lambda * (1.0 - J);

    double A = lambda * J * (2.0 * J - 1.0);        // 体積弾性
    double B = 2.0 * (mu - lambda * J * (J - 1.0)); // せん断弾性

    Ctan = MatrixXd::Zero(6, 6);

    for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
      {
        for(int k = 0; k < 3; k++)
        {
          for(int l = 0; l < 3; l++)
          {


            double term_vol1 =
                lambda * (2.0 * J - 1.0) * J * C_inv(i, j) * C_inv(k, l);
            double term_vol2 =
                -lambda * J * (J - 1.0) *
                (C_inv(i, k) * C_inv(j, l) + C_inv(i, l) * C_inv(j, k));
            double Ctan_vol_original = (term_vol1 + term_vol2);

            double term_iso1 = -1.0 / 3.0 * C_inv(k, l) *
                               (I(i, j) - 1.0 / 3.0 * traC_ * C_inv(i, j));
            double term_iso2 =
                -1.0 / 3.0 * I(k, l) * C_inv(i, j) +
                1.0 / 6.0 * traC_ *
                    (C_inv(i, k) * C_inv(j, l) + C_inv(i, l) * C_inv(j, k));

            double Ctan_iso_original =
                2.0 * mu * J_m23 * (term_iso1 + term_iso2);

            // Scale by retained damage also during unloading.
            Ctan(IJdata(i, j), IJdata(k, l)) =
                (1.0 - current_D) * (Ctan_vol_original + Ctan_iso_original);

            /*
            De_(IJdata(i, j), IJdata(k, l)) =
                A * C_inv(i, j) * C_inv(k, l) +
                (1.0 - current_D) * B * 0.5 *
                    (C_inv(i, k) * C_inv(j, l) + C_inv(i, l) * C_inv(j, k));
                    */
          }
        }
      }
    }

    //----------------------- for S:E integration ------------------//
    product += PK2(0, 0) * strainvec(0);
    product += PK2(1, 1) * strainvec(1);
    product += PK2(2, 2) * strainvec(2);

    product += 0.5 * (PK2(0, 1) + PK2(1, 0)) * strainvec(3);
    product += 0.5 * (PK2(1, 2) + PK2(2, 1)) * strainvec(4);
    product += 0.5 * (PK2(0, 2) + PK2(2, 0)) * strainvec(5);
    //--------------------------------------------------------------//

    // calculate Cauchy von Mises
    MatrixXd sigmaCauchy = (1.0 / J) * F * PK2 * F.transpose();

    // von Mises stress
    double sigma_xx = sigmaCauchy(0, 0);
    double sigma_yy = sigmaCauchy(1, 1);
    double sigma_zz = sigmaCauchy(2, 2);

    double sigma_xy = sigmaCauchy(0, 1);
    double sigma_yz = sigmaCauchy(1, 2);
    double sigma_zx = sigmaCauchy(2, 0);

    mises_ = std::sqrt(0.5 * ((sigma_xx - sigma_yy) * (sigma_xx - sigma_yy) +
                              (sigma_yy - sigma_zz) * (sigma_yy - sigma_zz) +
                              (sigma_zz - sigma_xx) * (sigma_zz - sigma_xx)) +
                       3.0 * (sigma_xy * sigma_xy + sigma_yz * sigma_yz +
                              sigma_zx * sigma_zx));
  }
}

//----------------------- for S:E integration ------------------//
void FiniteNeoHookedamage::makeKeFint(
    std::vector<Eigen::Triplet<double, int>> &tripletsK,
    std::vector<Eigen::Triplet<double, int>> &tripletsF,
    std::vector<Eigen::Triplet<double, int>> &tripletsR, FEM &fem,
    std::vector<Node> &nodes, int thisTime, double &product,
    double tauFilteredElement)
{
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd U = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd Uorg = MatrixXd::Zero(ne, fem.ndim);
  VectorXd disp = VectorXd::Zero(numdof);
  VectorXi idof = VectorXi::Zero(numdof);
  Gauss9Square gauss;

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
    {
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j]; // place of element's node point
      U.coeffRef(i, j) = gamma_ * nodes[nodeID[i]].val[j];
      disp.coeffRef(fem.ndim * i + j) = nodes[nodeID[i]].val[j];
      idof.coeffRef(fem.ndim * i + j) =
          nodes[nodeID[i]].dof[j]; // global dof number
    }
  }

  volume = 0.0;
  Ke_.setZero(); // elemental stiffness matrix
  fe_.setZero(); // elemental internal force matrix

  // linear init
  l_.init();       // object for linear analysis (l_)
  Bmatrix bLinear; // object for linear analysis(strain-displacement matrix)

  // nonlinear init
  MatrixXd nKe = Ke_;
  MatrixXd dnKe = Ke_;
  MatrixXd SKe = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
  MatrixXd dSnKe = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
  VectorXd nfe = fe_;
  BmatrixTotal bTotal;
  double Product_SE = 0.0;
  product = 0.0;

  double jac = 1.0;
  product_tempv.resize(ipmax);
  product_tempv.setZero();

  // gauss point loop start
  for(int ip = 0; ip < ipmax; ip++)
  {
    young_ = SIMP(young1_, young2_, design_s_, pp_);
    double w = gauss.weights(ip);
    double product_temp = 0.0;

    // linear part
    double J = 1.0;
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne); // shape function
    MatrixXd Be =
        MatrixXd::Zero(fem.voigt, numdof); // strain-displacement matrix ε = B u
    De(l_.De, young_, poisson_); // elastic matrix(stress-strain) σ = D ε
    bLinear.make(Be, Ne, jac, eType, ipmax, X, ip);
    l_.strain = Be * disp;
    l_.stress = l_.De * l_.strain;
    l_.fe += Be.transpose() * l_.stress * jac;
    l_.Ke += Be.transpose() * l_.De * Be * jac;

    // nonlinear part
    MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);

    // make bmatrix for finite strain
    bTotal.make(BeL, BeNL, Ne, jac, F, eType, ipmax, X, U, ip);
    assert(jac > 0.0);

    //  second piola-Kirchhoff stress tensor
    MatrixXd sigmaMatrix =
        MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd PK2 = MatrixXd::Zero(fem.ndim, fem.ndim);
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);
    MatrixXd Ctan = MatrixXd::Zero(fem.voigt, fem.voigt);

    double product = 0.0;

    totalKirchhoff(fem, PK2, sigmaVector, F, sigmaMatrix, product, Ctan, ip,
                   tauFilteredElement);

    Product_SE += product * w * jac;
    product_tempv[ip] = product;

    // inner stress(Bonet's book (9.15b))
    nfe += BeL.transpose() * sigmaVector * jac;

    // element stiffness matrix(Bonet's book (9.31),(9.35),(9.44c))
    nKe += BeL.transpose() * Ctan * BeL * jac +
           BeNL.transpose() * sigmaMatrix * BeNL * jac;

    volume += jac;
  }
  ///@sa Wang F et al., Interpolation scheme for fictitious domain techniques
  /// and topology optimization of finite strain elastic problems, 2014.
  Ke_ = (1.0 - gamma_ * gamma_) * l_.Ke + gamma_ * gamma_ * nKe;

  fe_ = (1.0 - gamma_ * gamma_) * l_.fe + gamma_ * nfe;
  nKe_ = nKe;
  nKel_ = nKe_;
  lKe_ = l_.Ke;

  lfe_ = l_.fe;
  nfe_ = nfe;

  product += Product_SE;

  assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
  assemblingReact(numdof, idof, fe_, tripletsR);
}

//--------------------------------------------------------------//

void DynamicNeoHookedamage::makeKeMeFint(
    std::vector<Eigen::Triplet<double>> &tripletsK,
    std::vector<Eigen::Triplet<double>> &tripletsM,
    std::vector<Eigen::Triplet<double>> &tripletsF,
    std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
    std::vector<Node> &nodes, std::string ExpOrImp, double &product)
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

  Ke_.setZero();
  Me_.setZero();
  fe_.setZero();
  volume = 0.0;

  // linear part
  l_.init();
  Bmatrix bLinear;

  // nonlinear part
  MatrixXd nKe = Ke_;
  VectorXd nfe = fe_;
  BmatrixTotal bTotal;

  Gauss9Square gauss;

  for(int ip = 0; ip < ipmax; ip++)
  {
    young_ = RAMP(young1_, young2_, design_s_, pp_);
    double w = gauss.weights(ip);

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
    Me_ += design_s_ * density_ * Ne.transpose() * Ne * jac;

    // nonlinear part
    jac = 1.0;
    MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
    // second piola-Kirchhoff stress tensor
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    MatrixXd PK2 = MatrixXd::Zero(fem.ndim, fem.ndim);
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);
    MatrixXd sigmaMatrix =
        MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd Ctan = MatrixXd::Zero(fem.voigt, fem.voigt);
    double tauFiltered;

    // make bmatrix for finite strain
    bTotal.make(BeL, BeNL, Ne, jac, F, eType, ipmax, X, U, ip);
    assert(jac > 0.0);

    double product = 0.0; // 積分点ごとのエネルギー
    totalKirchhoff(fem, PK2, sigmaVector, F, sigmaMatrix, product, Ctan, ip,
                   tauFiltered);
    product += product * w * jac;

    // inner stress(Bonet's book (9.15b))
    nfe += BeL.transpose() * sigmaVector * jac;

    // element stiffness matrix(Bonet's book (9.31),(9.35),(9.44c))
    nKe += BeL.transpose() * Ctan * BeL * jac +
           BeNL.transpose() * sigmaMatrix * BeNL * jac;

    volume += jac;
  }

  ///@sa Wang F et al., Interpolation scheme for fictitious domain techniques
  /// and topology optimization of finite strain elastic problems, 2014.
  Ke_ = (1.0 - gamma_ * gamma_) * l_.Ke + gamma_ * gamma_ * nKe;
  fe_ = (1.0 - gamma_ * gamma_) * l_.fe + gamma_ * nfe;

  lfe_ = l_.fe;
  nfe_ = nfe;
  nke_ = nKe;

  if(ExpOrImp == "explicit")
  {
    // concentrated mass matrix
    VectorXd invM = VectorXd::Constant(numdof, 1.0 / Me_.sum());
    assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
    assembling(numdof, fem.numeq, idof, invM, tripletsM);
  }
  else if(ExpOrImp == "implicit")
  {
    // consistent mass matrix
    assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
    assembling(numdof, fem.numeq, idof, Me_, tripletsM);
  }
  else
  {
    cerr << "error in neo_hooke_damage.cpp" << endl;
    exit(1);
  }

  assemblingReact(numdof, idof, fe_, tripletsR);
}

// 均質化法　ミクロな構造（細かい網目や細胞のような構造）を持つ材料を、シミュレーション上であたかも「均一な一つの材料」として扱うための計算手法//
void mHomoNeoHookedamage::makeKeHomo(vector<Triplet<double>> &tripletsK,
                                     vector<Triplet<double>> &tripletsF,
                                     vector<Triplet<double>> &tripletsR,
                                     FEM &fem, vector<HomoNode> &nodes, int dir)
{
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd U = MatrixXd::Zero(ne, fem.ndim);
  VectorXd disp = VectorXd::Zero(numdof);
  VectorXi idof = VectorXi::Zero(numdof);

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.dofnp; j++)
    {
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
      U.coeffRef(i, j) = gamma_ * nodes[nodeID[i]].nmtval[dir][j];
      disp.coeffRef(fem.ndim * i + j) = nodes[nodeID[i]].nmtval[dir][j];
      idof.coeffRef(fem.dofnp * i + j) = nodes[nodeID[i]].dof[j];
    }
  }

  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;

  // linear init
  l_.init();
  Bmatrix bLinear;

  // nonlinear init
  MatrixXd nKe = Ke_;
  VectorXd nfe = fe_;
  BmatrixTotal bTotal;
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
    MatrixXd sigmaMatrix =
        MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd Ctan = MatrixXd::Zero(fem.voigt, fem.voigt);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    MatrixXd PK2 = MatrixXd::Zero(fem.ndim, fem.ndim);
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);
    double tauFiltered;

    // make bmatrix
    bTotal.make(BeL, BeNL, Ne, jac, F, eType, ipmax, X, U, ip);
    assert(jac > 0.0);

    double product = 0.0;

    totalKirchhoff(fem, PK2, sigmaVector, F, sigmaMatrix, product, Ctan, ip,
                   tauFiltered);

    // inner stress(Bonet's book (9.15b))
    nfe += BeL.transpose() * sigmaVector * jac;

    // element stiffness matrix(Bonet's book (9.31),(9.35),(9.44c))
    nKe += BeL.transpose() * Ctan * BeL * jac +
           BeNL.transpose() * sigmaMatrix * BeNL * jac;
  }

  // element stiffness matrix(Bonet's book (9.31),(9.35),(9.44c))
  Ke_ = (1.0 - gamma_ * gamma_) * l_.Ke + gamma_ * gamma_ * nKe;

  // inner stress(Bonet's book (9.15b))
  fe_ = (1.0 - gamma_ * gamma_) * l_.fe + gamma_ * nfe;

  assemblingHomo(ne, numdof, idof, Ke_, fe_, tripletsK, tripletsF, fem, nodeID,
                 nodes);
  assemblingReactHomo(ne, numdof, idof, fe_, tripletsR, fem, nodeID, nodes);
}

void mHomoNeoHookedamage::makeMicroValues(FEM &fem, vector<HomoNode> &nodes,
                                          double tauFiltered)
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

  sEnergy_.setZero(); // ひずみエネ初期化
  strain_.setZero();  // ひずみ初期化
  stress_.setZero();  // 応力初期化
  volume = 0.0;
  BmatrixTotal bmatrix; // BmatrixTotalは設計図　bmatrixは実施に計算に利用
  for(int ip = 0; ip < ipmax; ip++)
  {
    double jac = 1.0;
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne); // 初期化
    MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
    // second piola-Kirchhoff stress tensor
    MatrixXd sigmaMatrix =
        MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd Ctan = MatrixXd::Zero(fem.voigt, fem.voigt);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    // stess for internal force
    MatrixXd PK2 = MatrixXd::Zero(fem.ndim, fem.ndim);
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);
    double tauFiltered;

    bmatrix.make(BeL, BeNL, Ne, jac, F, eType, ipmax, X, U, ip);
    // input: x,u,ip積分点番号,ipmax,etype要素の種類,
    // output:bel線形B,benl非線形B,ne形状関数,jacヤコビアン,f

    double product = 0.0;

    totalKirchhoff(fem, PK2, sigmaVector, F, sigmaMatrix, product, Ctan, ip,
                   tauFiltered);

    MatrixXd mstrain = BeL * disps;
    volume += jac;
    // strain_ += BeL * disp * fac;
    //  stress_ += De_ * Be * disp * fac;

    sEnergy_ += mstrain.transpose() * Ctan * mstrain * jac;
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