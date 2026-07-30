#include "multipli_elasto_plast_iso.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <basis/bmatrix_updated.hpp>
#include <misc/assembling.hpp>
#include <unsupported/Eigen/KroneckerProduct>

namespace icarat
{
using namespace std;
using namespace Eigen;

void MultipliEPUpdated::setParameter(const toml::value &config)
{
  const auto &mat = toml::find(config, "material");
  assert(toml::find<string>(mat, "type") == "MultipliEP");

  design_s_ = toml::find<double>(mat, "designs_0");
  pp_ = toml::find<double>(mat, "penalty");
  poisson_ = toml::find<double>(mat, "poisson");
  beta_ = toml::find<double>(mat, "beta");

  using pair = std::array<double, 2>;
  pair input;
  // young
  input = toml::find<pair>(mat, "young");
  young1_ = get<0>(input);
  young2_ = get<1>(input);

  // Hs
  Hs_ = toml::find<vector<double>>(mat, "Hs");

  // Ys
  Y0s_ = toml::find<vector<double>>(mat, "y0");
  lambdas_[0] =
      young1_ * poisson_ / ((1.0 + poisson_) * (1.0 - 2.0 * poisson_));
  mus_[0] = young1_ / (2 * (1.0 + poisson_));
  lambdas_[1] =
      young2_ * poisson_ / ((1.0 + poisson_) * (1.0 - 2.0 * poisson_));
  mus_[1] = young2_ / (2 * (1.0 + poisson_));

  mattype_ = toml::find<string>(mat, "state");
  if(toml::find<int>(config, "mesh", "dimension") == 3)
    mattype_ = "3D";
}

void MultipliEPUpdated::initialization()
{
  for(int ip = 0; ip < ipmax; ip++)
  {
    pstress_[ip].setZero();
    for(int i = 0; i < 3; i++)
    {
      pstress_[ip].coeffRef(i) = 1.0;
    }

    alpha_[ip].setZero();
    ep_[ip] = 0.0;
  }
}

void MultipliEPUpdated::updateVariables(FEM &fem, vector<Node> &nodes)
{
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd U = MatrixXd::Zero(ne, fem.ndim); ///<(24)
  VectorXi idof = VectorXi::Zero(numdof);
  // displacement in this increment
  MatrixXd du = MatrixXd::Zero(ne, fem.ndim); ///<(24)

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
    {
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
      U.coeffRef(i, j) = nodes[nodeID[i]].val[j];
      du.coeffRef(i, j) = nodes[nodeID[i]].du[j];
    }
  }

  // SIMP lame constant
  SIMPEP();

  // initialize
  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;
  BmatrixUpdated bmatrix;

  // gauss point loop start
  for(int ip = 0; ip < ipmax; ip++)
  {
    double jac = 0.0;
    lameDe(De_, mu_, lambda_);
    // shape function
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
    MatrixXd dNdx = MatrixXd::Zero(fem.ndim, ne);
    MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim); // Deformation gradient
    // second piola-Kirchhoff stress tensor
    MatrixXd sigmaMatrix =
        MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    // Directional Differentiation of Green Lagrangian Strain(DE[u])
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);
    MatrixXd LL = MatrixXd::Zero(fem.ndim, fem.ndim); // velocity gradient

    bmatrix.make(BeL, BeNL, Ne, jac, F, LL, du, eType, ipmax, X, U, ip);

    // Prepare for return mapping & tangent stiffness
    MatrixXd DDR = De_.block(0, 0, 3, 3);
    MatrixXd LLR = LL;
    VectorXd bbR = pstress_[ip];
    VectorXd alphaR = alpha_[ip];
    double epR = ep_[ip];

    // return mapping
    returnMapping(fem, sigmaVector, sigmaMatrix, DDR, LLR, bbR, alphaR, epR,
                  true);

    // stock values
    pstress_[ip] = bbR;
    alpha_[ip] = alphaR;
    ep_[ip] = epR;

  } // gauss loop end
}

void MultipliEPUpdated::returnMapping(FEM &fem, VectorXd &sigma,
                                      MatrixXd &sigmaMatrix, MatrixXd &DD,
                                      MatrixXd &LL, VectorXd &bb,
                                      VectorXd &alpha, double &ep, bool flag)
{
  double EPS = 1.0e-12;
  VectorXd Iden = VectorXd::Ones(3);
  MatrixXd eye = MatrixXd::Identity(3, 3);
  double two3 = 2.0 / 3.0;
  double stwo3 = sqrt(two3);
  double ftol = Y0_ * 1.0e-6;        ///< tolerance for yield
  MatrixXd R = (eye - LL).inverse(); ///< inc. deformation gradient

  MatrixXd bm(3, 3);
  bm << bb(0), bb(3), bb(5), //
      bb(3), bb(1), bb(4),   //
      bb(5), bb(4), bb(2);

  bm = R * bm * R.transpose(); ///< trial elastic left Cauchy-Green tensor

  bb << bm(0, 0), bm(1, 1), bm(2, 2), bm(0, 1), bm(1, 2), bm(0, 2);

  // cal eigenvalues & sort descending-order
  VectorXd P = bm.eigenvalues().real();
  vector<double> eigen((int)P.size());
  VectorXd::Map(&eigen[0], P.size()) = P;
  std::sort(eigen.begin(), eigen.end(), std::greater<double>());

  // stock original eigenvalues(principal stretch)
  auto eigen2 = eigen;

  // Duplicated eigenvalues (are shifted.)
  int TMP = -1;
  for(int i = 0; i < 2; i++)
  {
    if(abs(eigen[0] - eigen[2]) < EPS)
    {
      eigen[i] = eigen[i] + TMP * EPS;
      TMP = -TMP;
    }
  }
  if(abs(eigen[0] - eigen[1]) < EPS)
    eigen[1] = eigen[1] + EPS;
  if(abs(eigen[1] - eigen[2]) < EPS)
    eigen[1] = eigen[1] + EPS;

  // eigenvalue matrixM [6,3]
  MatrixXd M = MatrixXd::Zero(6, 3);
  for(int K = 0; K < 3; K++)
  {
    int KB = 1 + (K + 1) % 3;
    int KC = 1 + KB % 3;
    KB -= 1;
    KC -= 1;
    double EA = eigen[K];
    double EB = eigen[KB];
    double EC = eigen[KC];
    double D1 = EB - EA;
    double D2 = EC - EA;

    double DA = 1.0 / (D1 * D2);

    M(0, K) =
        ((bb(0) - EB) * (bb(0) - EC) + bb(3) * bb(3) + bb(5) * bb(5)) * DA;
    M(1, K) =
        ((bb(1) - EB) * (bb(1) - EC) + bb(3) * bb(3) + bb(4) * bb(4)) * DA;
    M(2, K) =
        ((bb(2) - EB) * (bb(2) - EC) + bb(4) * bb(4) + bb(5) * bb(5)) * DA;
    M(3, K) = (bb(3) * (bb(0) - EB + bb(1) - EC) + bb(4) * bb(5)) * DA;
    M(4, K) = (bb(4) * (bb(1) - EB + bb(2) - EC) + bb(3) * bb(5)) * DA;
    M(5, K) = (bb(5) * (bb(2) - EB + bb(0) - EC) + bb(3) * bb(4)) * DA;
  }

  //////////////////////////////////////////////////////////////
  /////////////////////return mapping///////////////////////////
  //////////////////////////////////////////////////////////////
  if(flag == true)
  {

    VectorXd deps(3); ///< logarithmic
    for(int i = 0; i < 3; i++)
      deps[i] = 0.5 * log(eigen2[i]);

    VectorXd sigtr = DD * deps; ///< trial proncipal stress
    VectorXd eta = sigtr - alpha - sigtr.sum() * Iden / 3.0; ///< shifted stress
    double etat = eta.norm();                                ///< norm of eta
    double fyld =
        etat - stwo3 * (Y0_ + (1.0 - beta_) * H_ * ep); ///< trial ield function

    if(fyld < ftol) ///< yield test is ok
    {
      VectorXd sig = sigtr; ///< trial stress are final
      sigma = M * sig;      ///< stress(6*1)
    }
    else ///< not ok
    {
      double gamma =
          fyld / (2 * mu_ + two3 * H_); ///< plastic consistency param
      ep = ep + gamma * stwo3;          ///< updated eff. plastic strain
      VectorXd N = eta / etat;          ///< unit vector normal to f
      deps = deps - gamma * N;          ///< updated elastic strain
      VectorXd sig = sigtr - 2.0 * mu_ * gamma * N;  ///< updated stress
      alpha = alpha + two3 * beta_ * H_ * gamma * N; ///< updated back stress
      VectorXd sigmaTmp = M * sig;                   ///< stress(6*1)
      sigma << sigmaTmp(0), sigmaTmp(1), sigmaTmp(2), sigmaTmp(3), sigmaTmp(5),
          sigmaTmp(4);
      VectorXd tmp(3);
      tmp << exp(2 * deps[0]), exp(2 * deps[1]), exp(2 * deps[2]);
      bb = M * tmp; ///< updated elastic left Cauchy-Green tensor
    }

    MatrixXd sig(3, 3);
    sig << sigma(0), sigma(3), sigma(5), //
        sigma(3), sigma(1), sigma(4),    //
        sigma(5), sigma(4), sigma(2);

    sigmaMatrix = kroneckerProduct(eye, sig);
  }
  //////////////////////////////////////////////////////////////
  /////////////////////tangent stiffness D//////////////////////
  //////////////////////////////////////////////////////////////
  else
  {
    VectorXd deps(3); ///< logarithmic
    for(int i = 0; i < 3; i++)
      deps[i] = 0.5 * log(eigen[i]);

    VectorXd sigtr = DD * deps; ///< trial proncipal stress
    VectorXd eta = sigtr - alpha - sigtr.sum() * Iden / 3.0; ///< shifted stress
    double etat = eta.norm();                                ///< norm of eta
    double fyld =
        etat - stwo3 * (Y0_ + (1.0 - beta_) * H_ * ep); ///< trial ield function

    VectorXd sig = sigtr;

    if(fyld >= ftol) ///< yield test is not ok
    {
      double gamma = fyld / (2.0 * mu_ + two3 * H_);
      VectorXd NN = eta / etat;             ///< unit vector normal to f
      sig = sigtr - 2.0 * mu_ * gamma * NN; ///< updated stress
      double var1 = 4.0 * pow(mu_, 2.0) / (2 * mu_ + two3 * H_); ///< coeff
      double var2 = 4.0 * pow(mu_, 2.0) * gamma / etat;          ///< coeff
      DD = DD - (var1 - var2) * (NN * NN.transpose()) +
           var2 * (Iden * Iden.transpose()) / 3.0; ///< tangent stiffness(3,3)
      DD(0, 0) = DD(0, 0) - var2;                  ///< contr. from 4th-order I
      DD(1, 1) = DD(1, 1) - var2;
      DD(2, 2) = DD(2, 2) - var2;
    }

    double J1 = 0.0;
    for(int i = 0; i < 3; i++)
      J1 += eigen[i];

    double J3 = eigen[0] * eigen[1] * eigen[2];

    VectorXd I2(6);
    I2 << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;

    // VectorXd I3 = VectorXd::Ones(6);

    MatrixXd I4 = MatrixXd::Identity(6, 6);
    I4(3, 3) = 0.5;
    I4(4, 4) = 0.5;
    I4(5, 5) = 0.5;

    MatrixXd Ibb(6, 6);
    Ibb << 0.0, pow(bb(3), 2.0) - bb(0) * bb(1),
        pow(bb(5), 2.0) - bb(0) * bb(2), 0.0, bb(3) * bb(5) - bb(0) * bb(4),
        0.0, ///< row 1

        bb(3) * bb(3) - bb(0) * bb(1), 0.0, pow(bb(4), 2.0) - bb(1) * bb(2),
        0.0, 0.0, bb(3) * bb(4) - bb(1) * bb(5), ///< row 2

        pow(bb(5), 2.0) - bb(0) * bb(2), pow(bb(4), 2.0) - bb(1) * bb(2), 0.0,
        bb(4) * bb(5) - bb(2) * bb(3), 0.0, 0.0, ///< row 3

        0.0, 0.0, bb(4) * bb(5) - bb(2) * bb(3),
        (bb(0) * bb(1) - pow(bb(3), 2.0)) / 2.0,
        (bb(1) * bb(5) - bb(3) * bb(4)) / 2.0,
        (bb(0) * bb(4) - bb(3) * bb(5)) / 2.0, ///< row 4

        bb(3) * bb(5) - bb(0) * bb(4), 0.0, 0.0,
        (bb(1) * bb(5) - bb(3) * bb(4)) / 2.0,
        (bb(1) * bb(2) - pow(bb(4), 2.0)) / 2.0,
        (bb(2) * bb(3) - bb(4) * bb(5)) / 2.0, ///< row 5

        0.0, bb(3) * bb(4) - bb(1) * bb(5), 0.0,
        (bb(0) * bb(4) - bb(3) * bb(5)) / 2.0,
        (bb(2) * bb(3) - bb(4) * bb(5)) / 2.0,
        (bb(0) * bb(2) - pow(bb(5), 2.0)) / 2.0; ///< row 6

    vector<double> dd(3), t1(3), t2(3), t3(3);
    dd[0] = 1.0 / ((eigen[1] - eigen[0]) * (eigen[2] - eigen[0]));
    dd[1] = 1.0 / ((eigen[2] - eigen[1]) * (eigen[0] - eigen[1]));
    dd[2] = 1.0 / ((eigen[0] - eigen[2]) * (eigen[1] - eigen[2]));
    t1[0] = -J3 * dd[0] / eigen[0];
    t1[1] = -J3 * dd[1] / eigen[1];
    t1[2] = -J3 * dd[2] / eigen[2];
    t2[0] = dd[0] * eigen[0];
    t2[1] = dd[1] * eigen[1];
    t2[2] = dd[2] * eigen[2];
    t3[0] = t2[0] * (J1 - 4.0 * eigen[0]);
    t3[1] = t2[1] * (J1 - 4.0 * eigen[1]);
    t3[2] = t2[2] * (J1 - 4.0 * eigen[2]);

    vector<MatrixXd> CTs;
    for(int i = 0; i < 3; i++)
    {
      MatrixXd CT =
          dd[i] * Ibb +
          t1[i] * (I4 - (I2 - M.col(i)) * (I2 - M.col(i)).transpose()) +
          t2[i] * (bb * M.col(i).transpose() + M.col(i) * bb.transpose()) +
          t3[i] * M.col(i) * M.col(i).transpose(); ///<(6,6)

      CTs.push_back(CT);
    }

    De_ = M * DD * M.transpose() +
          2.0 * (CTs[0] * sig(0) + CTs[1] * sig(1) + CTs[2] * sig(2)); ///<(6,6)
  }
}

void MultipliEPUpdated::SIMPEP()
{
  assert(young2_ >= young1_);
  young_ =
      pow(design_s_, pp_) * young2_ + (1.0 - pow(design_s_, pp_)) * young1_;
  lambda_ = pow(design_s_, pp_) * lambdas_[1] + (1.0 - design_s_) * lambdas_[0];
  mu_ = pow(design_s_, pp_) * mus_[1] + (1.0 - design_s_) * mus_[0];
  H_ = pow(design_s_, pp_) * Hs_[1] + (1.0 - design_s_) * Hs_[0];
  Y0_ = pow(design_s_, pp_) * Y0s_[1] + (1.0 - design_s_) * Hs_[0];
}

void FiniteMultipliEPUpdated::makeKeFint(vector<Triplet<double>> &tripletsK,
                                         vector<Triplet<double>> &tripletsF,
                                         vector<Triplet<double>> &tripletsR,
                                         FEM &fem, vector<Node> &nodes)
{
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd U = MatrixXd::Zero(ne, fem.ndim); ///<(24)
  VectorXi idof = VectorXi::Zero(numdof);
  // displacement in this increment
  MatrixXd du = MatrixXd::Zero(ne, fem.ndim); ///<(24)

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
    {
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
      U.coeffRef(i, j) = nodes[nodeID[i]].val[j];
      idof.coeffRef(fem.ndim * i + j) = nodes[nodeID[i]].dof[j];

      du.coeffRef(i, j) = nodes[nodeID[i]].du[j];
    }
  }

  // SIMP lame constant
  SIMPEP();

  // initialize
  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;
  BmatrixUpdated bmatrix;

  // gauss point loop start
  for(int ip = 0; ip < ipmax; ip++)
  {
    double jac = 0.0;
    lameDe(De_, mu_, lambda_);
    // shape function
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
    MatrixXd dNdx = MatrixXd::Zero(fem.ndim, ne);
    MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim); // Deformation gradient
    // second piola-Kirchhoff stress tensor
    MatrixXd sigmaMatrix =
        MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
    // Directional Differentiation of Green Lagrangian Strain(DE[u])
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, numdof);
    VectorXd sigmaVector = VectorXd::Zero(fem.voigt);
    MatrixXd LL = MatrixXd::Zero(fem.ndim, fem.ndim); // velocity gradient

    bmatrix.make(BeL, BeNL, Ne, jac, F, LL, du, eType, ipmax, X, U, ip);

    // Prepare for return mapping & tangent stiffness
    MatrixXd DDR = De_.block(0, 0, 3, 3);
    MatrixXd LLR = LL;
    VectorXd bbR = pstress_[ip];
    VectorXd alphaR = alpha_[ip];
    double epR = ep_[ip];

    MatrixXd DDT = De_.block(0, 0, 3, 3);
    MatrixXd LLT = LL;
    VectorXd bbT = pstress_[ip];
    VectorXd alphaT = alpha_[ip];
    double epT = ep_[ip];

    // return mapping
    returnMapping(fem, sigmaVector, sigmaMatrix, DDR, LLR, bbR, alphaR, epR,
                  true);

    fe_ += BeL.transpose() * sigmaVector * jac;
    Ke_ += BeNL.transpose() * sigmaMatrix * BeNL * jac;

    // update De matrix
    returnMapping(fem, sigmaVector, sigmaMatrix, DDT, LLT, bbT, alphaT, epT,
                  false);

    Ke_ += BeL.transpose() * De_ * BeL * jac;

    volume += jac;
  } // gauss point loop end

  assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
  assemblingReact(numdof, idof, fe_, tripletsR);
}

} // namespace icarat