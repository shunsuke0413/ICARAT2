///  @file  elasticity.cpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "elasticity.hpp"
#include <basis/bmatrix.hpp>
#include <basis/bmatrix_shell.hpp>
#include <iostream>
#include <misc/assembling.hpp>

namespace icarat
{
using namespace std;
using namespace Eigen;

void Elastic::setParameter(const toml::value &config)
{
  const auto &mat = toml::find(config, "material");
  assert(toml::find<string>(mat, "type") == "Elastic");

  design_s_ = toml::find<double>(mat, "designs_0");
  pp_ = toml::find<double>(mat, "penalty");
  auto young = toml::find<std::array<double, 2>>(mat, "young");
  young1_ = get<0>(young);
  young2_ = get<1>(young);
  poisson_ = toml::find<double>(mat, "poisson");
  mattype_ = toml::find<string>(mat, "state");
  if(toml::find<int>(config, "mesh", "dimension") == 3)
    mattype_ = "3D";
}

void Elastic::De(MatrixXd &De, double E, double nu)
{
  De.setZero();

  // 2D
  if(De.rows() == 3)
  {

    if(mattype_ == "plane_stress")
    {
      De.coeffRef(0, 0) = 1.0;
      De.coeffRef(0, 1) = nu;
      De.coeffRef(1, 1) = 1.0;
      De.coeffRef(2, 2) = (0.5 * (1.0 - nu));
      De.coeffRef(1, 0) = De(0, 1);
      De.coeffRef(2, 0) = De(0, 2);
      De.coeffRef(2, 1) = De(1, 2);
      De *= E / (1.0 - nu * nu);
    }
    else if(mattype_ == "plane_strain")
    {
      De.coeffRef(0, 0) = 1.0;
      De.coeffRef(0, 1) = (nu / (1.0 - nu));
      De.coeffRef(1, 1) = 1.0;
      De.coeffRef(2, 2) = ((1.0 - 2.0 * nu) / (2.0 * (1.0 - nu)));
      De.coeffRef(1, 0) = De(0, 1);
      De.coeffRef(2, 0) = De(0, 2);
      De.coeffRef(2, 1) = De(1, 2);
      De *= E * ((1.0 - nu) / ((1.0 - 2.0 * nu) * (1.0 + nu)));
    }
    else
    {
      cerr << "Error in elasticity.cpp" << endl;
      exit(1);
    }
  }
  // 3D
  else if(De.rows() == 6)
  {
    // lame constant
    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double mu = E / (2 * (1.0 + nu));

    De.coeffRef(0, 0) = lambda + 2 * mu;
    De.coeffRef(0, 1) = lambda;
    De.coeffRef(0, 2) = lambda;

    De.coeffRef(1, 0) = lambda;
    De.coeffRef(1, 1) = lambda + 2 * mu;
    De.coeffRef(1, 2) = lambda;

    De.coeffRef(2, 0) = lambda;
    De.coeffRef(2, 1) = lambda;
    De.coeffRef(2, 2) = lambda + 2 * mu;

    De.coeffRef(3, 3) = mu;
    De.coeffRef(4, 4) = mu;
    De.coeffRef(5, 5) = mu;
  }
  else
  {
    cerr << "Error in elasticity.cpp" << endl;
    exit(1);
  }
}

void Elastic::lameDe(MatrixXd &De, double mu, double lambda)
{
  De.setZero();

  // 2D
  if(De.rows() == 3)
  {
    // convert to E and nu
    double nu = lambda / (2.0 * (lambda + mu));
    double E = 2.0 * mu * (1.0 + nu);

    if(mattype_ == "plane_stress")
    {
      De.coeffRef(0, 0) = 1.0;
      De.coeffRef(0, 1) = nu;
      De.coeffRef(1, 1) = 1.0;
      De.coeffRef(2, 2) = (0.5 * (1.0 - nu));
      De.coeffRef(1, 0) = De(0, 1);
      De.coeffRef(2, 0) = De(0, 2);
      De.coeffRef(2, 1) = De(1, 2);
      De *= E / (1.0 - nu * nu);
    }
    else if(mattype_ == "plane_strain")
    {
      De.coeffRef(0, 0) = 1.0;
      De.coeffRef(0, 1) = (nu / (1.0 - nu));
      De.coeffRef(1, 1) = 1.0;
      De.coeffRef(2, 2) = ((1.0 - 2.0 * nu) / (2.0 * (1.0 - nu)));
      De.coeffRef(1, 0) = De(0, 1);
      De.coeffRef(2, 0) = De(0, 2);
      De.coeffRef(2, 1) = De(1, 2);
      De *= E * ((1.0 - nu) / ((1.0 - 2.0 * nu) * (1.0 + nu)));
    }
    else
    {
      cerr << "Error in elasticity.cpp" << endl;
      exit(1);
    }
  }
  // 3D
  else if(De.rows() == 6)
  {
    De.coeffRef(0, 0) = lambda + 2 * mu;
    De.coeffRef(0, 1) = lambda;
    De.coeffRef(0, 2) = lambda;

    De.coeffRef(1, 0) = lambda;
    De.coeffRef(1, 1) = lambda + 2 * mu;
    De.coeffRef(1, 2) = lambda;

    De.coeffRef(2, 0) = lambda;
    De.coeffRef(2, 1) = lambda;
    De.coeffRef(2, 2) = lambda + 2 * mu;

    De.coeffRef(3, 3) = mu;
    De.coeffRef(4, 4) = mu;
    De.coeffRef(5, 5) = mu;
  }
  else
  {
    cerr << "Error in elasticity.cpp" << endl;
    exit(1);
  }
}

void Elastic::makeMisesStress(FEM &fem, vector<Node> &node)
{
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
  MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  VectorXd disp = VectorXd::Zero(numdof);
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

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
    {
      disp.coeffRef(fem.ndim * i + j) = node[nodeID[i]].val[j];
      X.coeffRef(i, j) = node[nodeID[i]].x[j];
    }
  }

  strain_.setZero();
  stress_.setZero();

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

void LinearElastic::makeKeFint(vector<Triplet<double>> &tripletsK,
                               vector<Triplet<double>> &tripletsF,
                               vector<Triplet<double>> &tripletsR, FEM &fem,
                               vector<Node> &nodes)
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
      disp.coeffRef(fem.dofnp * i + j) = nodes[nodeID[i]].val[j];
      idof.coeffRef(fem.dofnp * i + j) = nodes[nodeID[i]].dof[j];
    }
    for(int j = 0; j < fem.ndim; j++)
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
  }

  young_ = SIMP(young1_, young2_, design_s_, pp_);

  De(De_, young_, poisson_);

  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;
  double jac = 0.0;
  Bmatrix bmatrix;

  // gauss loop
  for(int ip = 0; ip < ipmax; ip++)
  {
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);

    fe_ += Be.transpose() * De_ * Be * disp * jac;
    Ke_ += Be.transpose() * De_ * Be * jac;
    volume += jac;
  }

  assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
  assemblingReact(numdof, idof, fe_, tripletsR);
}

void BucklingLinearElastic::makeEigenCoeffs(
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

void ModalLinearElastic::makeEigenCoeffs(vector<Triplet<double>> &tripletsA,
                                         vector<Triplet<double>> &tripletsB,
                                         FEM &fem, vector<Node> &nodes)
{
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
  MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  VectorXd disp = VectorXd::Zero(numdof);
  VectorXi idof = VectorXi::Zero(numdof);

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
    {
      disp.coeffRef(fem.ndim * i + j) = nodes[nodeID[i]].val[j];
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
      idof.coeffRef(fem.ndim * i + j) = nodes[nodeID[i]].dof[j];
    }
  }

  young_ = SIMP(young1_, young2_, design_s_, pp_);

  De(De_, young_, poisson_);

  Ke_.setZero();
  Me_.setZero();
  volume = 0.0;
  double jac = 0.0;
  Bmatrix bmatrix;

  for(int ip = 0; ip < ipmax; ip++)
  {
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);

    Ke_ += Be.transpose() * De_ * Be * jac;
    Me_ += design_s_ * density_ * Ne.transpose() * Ne * jac;

    volume += jac;
  }

  assembling(numdof, fem.numeq, idof, Ke_, tripletsA);
  assembling(numdof, fem.numeq, idof, Me_, tripletsB);
}

void DynamicElastic::makeKeMeFint(vector<Triplet<double>> &tripletsK,
                                  vector<Triplet<double>> &tripletsM,
                                  vector<Triplet<double>> &tripletsF,
                                  vector<Triplet<double>> &tripletsR, FEM &fem,
                                  vector<Node> &nodes, string ExpOrImp)
{
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
  MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  VectorXd disp = VectorXd::Zero(numdof);
  VectorXi idof = VectorXi::Zero(numdof);

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
    {
      disp.coeffRef(fem.ndim * i + j) = nodes[nodeID[i]].val[j];
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
      idof.coeffRef(fem.ndim * i + j) = nodes[nodeID[i]].dof[j];
    }
  }

  young_ = SIMP(young1_, young2_, design_s_, pp_);

  De(De_, young_, poisson_);

  Ke_.setZero();
  Me_.setZero();
  fe_.setZero();
  volume = 0.0;
  double jac = 0.0;
  Bmatrix bmatrix;

  for(int ip = 0; ip < ipmax; ip++)
  {
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);

    fe_ += Be.transpose() * De_ * Be * disp * jac;
    Ke_ += Be.transpose() * De_ * Be * jac;
    Me_ += design_s_ * density_ * Ne.transpose() * Ne * jac;

    volume += jac;
  }

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
    cerr << "error in elasticity.cpp" << endl;
    exit(1);
  }

  assemblingReact(numdof, idof, fe_, tripletsR);
}

void mHomoElastic::makeKeHomo(vector<Triplet<double>> &tripletsK,
                              vector<Triplet<double>> &tripletsF,
                              vector<Triplet<double>> &tripletsR, FEM &fem,
                              vector<HomoNode> &nodes, int dir)
{
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
  MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  VectorXi idof = VectorXi::Zero(numdof);
  VectorXd disp = VectorXd::Zero(numdof);

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.dofnp; j++)
    {
      idof.coeffRef(fem.dofnp * i + j) = nodes[nodeID[i]].dof[j];
      disp.coeffRef(fem.dofnp * i + j) = nodes[nodeID[i]].nmtval[dir][j];
    }

    for(int j = 0; j < fem.ndim; j++)
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
  }

  young_ = SIMP(young1_, young2_, design_s_, pp_);
  De(De_, young_, poisson_);

  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;
  double jac = 0.0;
  Bmatrix bmatrix;
  for(int ip = 0; ip < ipmax; ip++)
  {
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);
    Ke_ += Be.transpose() * De_ * Be * jac;
    fe_ += Be.transpose() * De_ * Be * disp * jac;
    volume += jac;
  }

  assemblingHomo(ne, numdof, idof, Ke_, fe_, tripletsK, tripletsF, fem, nodeID,
                 nodes);
  assemblingReactHomo(ne, numdof, idof, fe_, tripletsR, fem, nodeID, nodes);
}

void mHomoElastic::makeMicroValues(FEM &fem, vector<HomoNode> &nodes)
{
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
  MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
  VectorXi idof = VectorXi::Zero(numdof);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  VectorXd disp = VectorXd::Zero(numdof);             ///< for localization
  MatrixXd disps = MatrixXd::Zero(numdof, fem.voigt); ///< for nmt

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.dofnp; j++)
    {
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
      idof.coeffRef(i * fem.dofnp + j) = nodes[nodeID[i]].dof[j];
      disp.coeffRef(i * fem.dofnp + j) = nodes[nodeID[i]].val[j];

      for(int k = 0; k < fem.voigt; k++)
        disps(fem.dofnp * i + j, k) = nodes[nodeID[i]].nmtval[k][j];
    }
  }

  sEnergy_.setZero();
  strain_.setZero();
  stress_.setZero();
  volume = 0.0;
  double jac = 0.0;
  double fac = (double)(1.0 / ipmax);
  Bmatrix bmatrix;
  for(int ip = 0; ip < ipmax; ip++)
  {
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);
    volume += jac;
    strain_ += Be * disp * fac;
    stress_ += De_ * Be * disp * fac;

    MatrixXd mstrain = Be * disps;
    sEnergy_ += mstrain.transpose() * De_ * mstrain * jac;
  }

  // make mises stress
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

  mises_ = sqrt(stress_.transpose() * V * stress_);
}

void MHomoElastic::makeKeFint(vector<Triplet<double>> &tripletsK,
                              vector<Triplet<double>> &tripletsF,
                              vector<Triplet<double>> &tripletsR, FEM &fem,
                              vector<Node> &nodes)
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
      disp.coeffRef(fem.dofnp * i + j) = nodes[nodeID[i]].val[j];
      idof.coeffRef(fem.dofnp * i + j) = nodes[nodeID[i]].dof[j];
    }
    for(int j = 0; j < fem.ndim; j++)
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
  }

  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;
  double jac = 0.0;
  Bmatrix bmatrix;
  for(int ip = 0; ip < ipmax; ip++)
  {
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);

    fe_ += Be.transpose() * De_ * Be * disp * jac;
    Ke_ += Be.transpose() * De_ * Be * jac;
    volume += jac;
  }

  assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
  assemblingReact(numdof, idof, fe_, tripletsR);
}

double ThermalStressElastic::thermalStress(const double &designS,
                                           const double &young1,
                                           const double &young2,
                                           const double &alpha1,
                                           const double &alpha2)
{
  double beta1 = young1 * alpha1;
  double beta2 = young2 * alpha2;
  double beta;
  if(beta1 >= beta2)
    beta = (beta1 - beta2) * designS + beta2;
  else
    beta = (beta2 - beta1) * designS + beta1;
  return beta;
}

VectorXd ThermalStressElastic::generateThermalLoad(FEM &fem, vector<Node> &node,
                                                   double deltaT)
{
  VectorXd phi(fem.voigt);
  if(fem.ndim == 2)
    phi << 1.0, 1.0, 0.0;
  else
    phi << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;

  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
  MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  VectorXi idof = VectorXi::Zero(numdof);

  for(int j = 0; j < ne; j++)
  {
    for(int k = 0; k < fem.ndim; k++)
    {
      idof.coeffRef(fem.dofnp * j + k) = node[nodeID[j]].dof[k];
      X.coeffRef(j, k) = node[nodeID[j]].x[k];
    }
  }

  VectorXd fthi = VectorXd::Zero(numdof);
  De(De_, 1.0, poisson_);
  Bmatrix bmatrix;
  for(int ip = 0; ip < ipmax; ip++)
  {
    double jac;
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);
    fthi += Be.transpose() * De_ * phi * jac;
  }

  beta_ = thermalStress(design_s_, young1_, young2_, alpha_[0], alpha_[1]);
  fthi *= deltaT * beta_;
  fthi_ = fthi;

  // assembling
  VectorXd fth = VectorXd::Zero(fem.numeq);
  for(int i = 0; i < numdof; i++)
  {
    if(idof[i] < fem.numeq)
      fth[idof[i]] = fthi[i];
    else
      continue;
  }

  return fth;
}

void ThermalStressElastic::makeKeFint(vector<Triplet<double>> &tripletsK,
                                      vector<Triplet<double>> &tripletsF,
                                      vector<Triplet<double>> &tripletsR,
                                      FEM &fem, vector<Node> &nodes)
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
      disp.coeffRef(fem.dofnp * i + j) = nodes[nodeID[i]].val[j];
      idof.coeffRef(fem.dofnp * i + j) = nodes[nodeID[i]].dof[j];
    }
    for(int j = 0; j < fem.ndim; j++)
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
  }

  young_ = RAMP(young1_, young2_, design_s_, pp_);

  De(De_, young_, poisson_);

  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;
  double jac = 0.0;
  Bmatrix bmatrix;

  // gauss loop
  for(int ip = 0; ip < ipmax; ip++)
  {
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);

    fe_ += Be.transpose() * De_ * Be * disp * jac;
    Ke_ += Be.transpose() * De_ * Be * jac;
    volume += jac;
  }

  assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
  assemblingReact(numdof, idof, fe_, tripletsR);
}
} // namespace icarat