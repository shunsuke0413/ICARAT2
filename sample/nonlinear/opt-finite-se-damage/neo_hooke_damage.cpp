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

double neohookedamage::calculateTensileDamage(double tau)
{
  const double tau_0 = 0.42;

  double D = (tau >= tau_0)
                 ? 1.0 - ((-1.083 * tau + 51.438) / (54.703 * tau + 21.477))
                 : 0.0;
  D = std::max(0.0, std::min(0.999, D));

  return D;
}

double neohookedamage::differentiateD(double tau) {
  const double tau_0 = 0.42;
  // D = 1.0 - ((A*tau + B) / (C*tau + D)) の微分
  // D' = - ( (A*(C*tau+D) - (A*tau+B)*C) / (C*tau+D)^2 )

  double A = -1.083; double B = 51.438;
  double C = 54.703; double D = 21.477;
  
  double denom = C * tau + D;
  double d_f_d_tau = (A * D - B * C) / (denom * denom);
  
  return -d_f_d_tau;
}


//----------------------- for S:E integration ------------------//
void neohookedamage::totalKirchhoff(FEM &fem, VectorXd &sigmaVector,
                                    MatrixXd &F, MatrixXd &sigmaMatrix,
                                    double &product, int ip)
//--------------------------------------------------------------//
{
  // lame constant for initial & current position　ラメ定数p.194
  double mu = 0.5 * young_ / (1 + poisson_);
  double lambda =
      (poisson_ * young_) / ((1.0 + poisson_) * (1.0 - 2.0 * poisson_));

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

  double jac = F.determinant();
  double jac_m23 = pow(jac, -2.0 / 3.0);
  // 1. せん断ひずみエネルギー密度 (Psi_iso) の算出
  // 偏差成分の計算に必要な指標です
  double Psi_iso = 0.5 * mu * (traC_ * jac_m23 - 3.0);

  // 2. Energy norm
  double tau = sqrt(2.0 * std::abs(Psi_iso));

  // ガウス点ごとのリストを更新
  D_all_[ip] = calculateTensileDamage(tau);

  // 今後の計算で使う現在の損傷度
  double current_D = D_all_[ip];

  if(fem.ndim == 2)
  {
    // tensol -> vector */
    VectorXd strainvec = VectorXd::Zero(fem.voigt);
    strainvec[0] = E(0, 0);
    strainvec[1] = E(1, 1);
    strainvec[2] = 2.0 * E(0, 1);

    // calculate the inverse of matrix c

    double jacC = pow(jac, 2.0);
    MatrixXd C_inv = MatrixXd::Zero(fem.ndim, fem.ndim);
    C_inv(0, 0) = C(1, 1) / jacC;
    C_inv(0, 1) = -C(0, 1) / jacC;
    C_inv(1, 0) = -C(1, 0) / jacC;
    C_inv(1, 1) = C(0, 0) / jacC;

    cinv_ = C_inv;

    // cal stress　第二ピオラ・キルヒホッフ応力テンソル？
    sigmaVector(0) = (1.0 - current_D) * 0.5 * mu * jac_m23 * (1.0 - C_inv(0, 0)) +
                     lambda * jac * (jac - 1.0) * C_inv(0, 0);
    sigmaVector(1) = (1.0 - current_D) * 0.5 * mu * jac_m23 * (1.0 - C_inv(1, 1)) +
                     lambda * jac * (jac - 1.0) * C_inv(1, 1);
    sigmaVector(2) = (1.0 - current_D) * 0.5 * mu * jac_m23 * (-C_inv(0, 1)) +
                     lambda * jac * (jac - 1.0) * C_inv(0, 1);

    // cal 2nd PK stress S
    MatrixXd PK2 = MatrixXd::Zero(fem.ndim, fem.ndim);
    PK2(0, 0) = sigmaVector(0);
    PK2(0, 1) = sigmaVector(2);
    PK2(1, 0) = sigmaVector(2);
    PK2(1, 1) = sigmaVector(1);
    sigmaMatrix = kroneckerProduct(PK2, I);

    // De matrix
    // make IJKL data
    MatrixXi IJdata = MatrixXi::Zero(fem.ndim, fem.ndim);
    IJdata(0, 0) = 0;
    IJdata(0, 1) = 2;
    IJdata(1, 0) = 2;
    IJdata(1, 1) = 1;

    double p = lambda * (1.0 - jac);

    double A = lambda * jac * (2.0 * jac - 1.0);        // 体積弾性
    double B = 2.0 * (mu - lambda * jac * (jac - 1.0)); // せん断弾性

    // Ctan : 物質弾性係数テンソル

    for(int i = 0; i < 2; i++)
    {
      for(int j = 0; j < 2; j++)
      {
        for(int k = 0; k < 2; k++)
        {
          for(int l = 0; l < 2; l++)
          {
            double term1_vol = p * jac * (C_inv(i, j) * C_inv(k, l) - (C_inv(i, k) * C_inv(j, l) + C_inv(i, l) * C_inv(j, k)));
            double term2_vol = -lambda * (jac * jac) * C_inv(i, j) * C_inv(k, l);

            double Ctan_vol_ijkl = term_vol1 + term_vol2;



            /*
            // 損傷の進展による微分値を計算
            double dD_dtau = differentiateD(tau);

            // 既存のDeの計算値（弾性成分）
            double elastic_part = A * C_inv(i, j) * C_inv(k, l) +
                                  (1.0 - current_D) * B * 0.5 *
                                      (C_inv(i, k) * C_inv(j, l) + C_inv(i, l) * C_inv(j, k));

            // 【重要】損傷の進展項（補正項）を引く
            // 理論上：補正項 = (dD/dtau) * (dTau/dSigma) * (Sigma ⊗ Sigma)
            // ここでは簡略化のため、この軟化特性を反映させる項を適用します
            double damping = 0.1;
            double correction = damping * dD_dtau * (PK2(i, j) * PK2(k, l));

            De_(IJdata(i, j), IJdata(k, l)) = elastic_part - correction;
            */
          }
        }
      }
    }

    //----------------------- for S:E integration ------------------//
    product += PK2(0, 0) * strainvec(0);
    product += PK2(1, 1) * strainvec(1);
    product += 0.5 * (PK2(0, 1) + PK2(1, 0)) * strainvec(2);
    //--------------------------------------------------------------//
  }
  else // 3D
  {
    // tensol -> vector */
    VectorXd strainvec = VectorXd::Zero(fem.voigt);
    strainvec[0] = E(0, 0);
    strainvec[1] = E(1, 1);
    strainvec[2] = E(2, 2);
    strainvec[3] = 2.0 * E(0, 1);
    strainvec[4] = 2.0 * E(1, 2);
    strainvec[5] = 2.0 * E(2, 0);
    // cal stress S */
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

    sigmaVector(0) = (1.0 - current_D) * mu * (1.0 - C_inv(0, 0)) +
                     lambda * jac * (jac - 1.0) * C_inv(0, 0);
    sigmaVector(1) = (1.0 - current_D) * mu * (1.0 - C_inv(1, 1)) +
                     lambda * jac * (jac - 1.0) * C_inv(1, 1);
    sigmaVector(2) = (1.0 - current_D) * mu * (1.0 - C_inv(2, 2)) +
                     lambda * jac * (jac - 1.0) * C_inv(2, 2);
    sigmaVector(3) = (1.0 - current_D) * mu * (-C_inv(0, 1)) +
                     lambda * jac * (jac - 1.0) * C_inv(0, 1);
    sigmaVector(4) = (1.0 - current_D) * mu * (-C_inv(1, 2)) +
                     lambda * jac * (jac - 1.0) * C_inv(1, 2);
    sigmaVector(5) = (1.0 - current_D) * mu * (-C_inv(2, 0)) +
                     lambda * jac * (jac - 1.0) * C_inv(2, 0);

    MatrixXd PK2 = MatrixXd::Zero(fem.ndim, fem.ndim);
    PK2 << sigmaVector[0], sigmaVector[3], sigmaVector[5], // xx, xy, xz
        sigmaVector[3], sigmaVector[1], sigmaVector[4],    // xy, yy, yz
        sigmaVector[5], sigmaVector[4], sigmaVector[2];    // xz, yz, zz

    sigmaMatrix = kroneckerProduct(PK2, I);

    // De matrix
    // make Voigt-Tensor index data
    MatrixXi IJdata = MatrixXi::Zero(fem.ndim, fem.ndim);
    IJdata << 0, 3, 5, // xx, xy, xz
        3, 1, 4,       // xy, yy, yz
        5, 4, 2;       // xz, yz, zz

    double A = lambda * jac * (2.0 * jac - 1.0);
    double B = 2.0 * (mu - lambda * jac * (jac - 1.0));

    for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
      {
        for(int k = 0; k < 3; k++)
        {
          for(int l = 0; l < 3; l++)
          {
            De_(IJdata(i, j), IJdata(k, l)) =
                A * C_inv(i, j) * C_inv(k, l) +
                (1.0 - current_D) * B * 0.5 *
                    (C_inv(i, k) * C_inv(j, l) + C_inv(i, l) * C_inv(j, k));
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
  }
}

//----------------------- for S:E integration ------------------//
void FiniteNeoHookedamage::makeKeFint(
    std::vector<Eigen::Triplet<double, int>> &tripletsK,
    std::vector<Eigen::Triplet<double, int>> &tripletsF,
    std::vector<Eigen::Triplet<double, int>> &tripletsR, FEM &fem,
    std::vector<Node> &nodes, int thisTime, double &product)
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
    double jac = 1.0;
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
    //  second piola-Kirchhoff stress tensor
    MatrixXd Sbar = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);

    // make bmatrix for finite strain
    bTotal.make(BeL, BeNL, Ne, jac, F, eType, ipmax, X, U, ip);
    assert(jac > 0.0);

    double product = 0.0;

    totalKirchhoff(fem, sigmaVector, F, Sbar, product, ip);

    Product_SE += product_temp * w * jac;
    product_tempv[ip] = product_temp;

    // inner stress(Bonet's book (9.15b))
    nfe += BeL.transpose() * sigmaVector * jac;

    // element stiffness matrix(Bonet's book (9.31),(9.35),(9.44c))
    nKe += BeL.transpose() * De_ * BeL * jac +
           BeNL.transpose() * Sbar * BeNL * jac;

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
    MatrixXd Sbar = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);

    // make bmatrix for finite strain
    bTotal.make(BeL, BeNL, Ne, jac, F, eType, ipmax, X, U, ip);
    assert(jac > 0.0);

    double product = 0.0; // 積分点ごとのエネルギー
    totalKirchhoff(fem, sigmaVector, F, Sbar, product, ip);
    product += product * w * jac;

    // inner stress(Bonet's book (9.15b))
    nfe += BeL.transpose() * sigmaVector * jac;

    // element stiffness matrix(Bonet's book (9.31),(9.35),(9.44c))
    nKe += BeL.transpose() * De_ * BeL * jac +
           BeNL.transpose() * Sbar * BeNL * jac;

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
    MatrixXd Sbar = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);

    // make bmatrix
    bTotal.make(BeL, BeNL, Ne, jac, F, eType, ipmax, X, U, ip);
    assert(jac > 0.0);

    double product = 0.0;

    totalKirchhoff(fem, sigmaVector, F, Sbar, product, ip);

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

  assemblingHomo(ne, numdof, idof, Ke_, fe_, tripletsK, tripletsF, fem, nodeID,
                 nodes);
  assemblingReactHomo(ne, numdof, idof, fe_, tripletsR, fem, nodeID, nodes);
}

void mHomoNeoHookedamage::makeMicroValues(FEM &fem, vector<HomoNode> &nodes)
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
    MatrixXd Sbar = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    // stess for internal force
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);

    bmatrix.make(BeL, BeNL, Ne, jac, F, eType, ipmax, X, U, ip);
    // input: x,u,ip積分点番号,ipmax,etype要素の種類,
    // output:bel線形B,benl非線形B,ne形状関数,jacヤコビアン,f

    double product = 0.0;

    totalKirchhoff(fem, sigmaVector, F, Sbar, product, ip);

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