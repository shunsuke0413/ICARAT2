///  @file  elasto_plastic_mises.cpp
///  @author  Daiki Watanabe
///  @date  October 3, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "elasto_plast_iso.hpp"
#include <basis/bmatrix.hpp>
#include <iostream>
#include <misc/assembling.hpp>

namespace icarat
{
using namespace std;
using namespace Eigen;

void ElastPlast::setParameter(const toml::value &config)
{
  const auto &mat = toml::find(config, "material");
  assert(toml::find<string>(mat, "type") == "ElastPlast");

  design_s_ = toml::find<double>(mat, "designs_0");
  pp_ = toml::find<double>(mat, "penalty");
  poisson_ = toml::find<double>(mat, "poisson");
  betaH_ = toml::find<double>(mat, "beta");

  using pair = std::array<double, 2>;
  pair input;
  input = toml::find<pair>(mat, "young");
  young1_ = get<0>(input);
  young2_ = get<1>(input);
  Ehs_ = toml::find<vector<double>>(mat, "h");
  sigma_ys_ = toml::find<vector<double>>(mat, "y0");

  // state
  mattype_ = toml::find<string>(mat, "state");
  if(toml::find<int>(config, "mesh", "dimension") == 4)
    mattype_ = "3D";
}

void ElastPlast::initialization()
{
  epsilon_ = 0.0;
  eqstn_ = 0.0;
  mises_ = 0.0;
  //  tr_stress_ = 0.0;
  for(int ip = 0; ip < ipmax; ip++)
  {
    strain_p_[ip].setZero();
    beta_[ip].setZero();
    strain_p_bar_[ip] = 0.0;
    xi_tri_norm_[ip] = 0.0;
  }
}

void ElastPlast::plasto_simo(const VectorXd &strain, VectorXd &stress, int ip)
{
  if(mattype_ == "plane_strain")
    plasto_simo_plane_strain(strain, stress, ip);
  else if(mattype_ == "3D")
    plasto_simo_3d(strain, stress, ip);
  else
  {
    cerr << "mattype is incorrect in elasto_plast_iso.cpp" << endl;
    exit(1);
  }
}

void ElastPlast::stress_plastic(FEM &fem, std::vector<Node> &nodes,
                                std::vector<Eigen::Triplet<double>> &tripletsR)
{
  if(mattype_ == "plane_strain")
    stress_plastic_plane_strain(fem, nodes, tripletsR);
  else if(mattype_ == "3D")
    stress_plastic_3d(fem, nodes, tripletsR);
  else
  {
    cerr << "mattype is incorrect in elasto_plast_iso.cpp" << endl;
    exit(1);
  }
}

void ElastPlast::SIMPEP()
{
  assert(young2_ >= young1_);
  young_ =
      pow(design_s_, pp_) * young2_ + (1.0 - pow(design_s_, pp_)) * young1_;
  Eh_ = pow(design_s_, pp_) * Ehs_[1] + (1.0 - pow(design_s_, pp_)) * Ehs_[0];
  sigma_y_ = pow(design_s_, pp_) * sigma_ys_[1] +
             (1.0 - pow(design_s_, pp_)) * sigma_ys_[0];
}

/////////////////////////////////////////////////////////////
/////////////////////////// private /////////////////////////
/////////////////////////////////////////////////////////////

void ElastPlast::plasto_simo_plane_strain(const VectorXd &strain,
                                          VectorXd &stress, int ip)
{

  SIMPEP();

  double hards = (young_ * Eh_) / (young_ - Eh_);
  double Ehi = betaH_ * hards;
  double G = (0.5 * young_) / (1.0 + poisson_);
  double K = young_ / (3.0 * (1.0 - 2.0 * poisson_));

  // A matrix
  VectorXd I(4);
  I << 1.0, 1.0, 0.0, 1.0;

  // I matrix
  MatrixXd I_s = MatrixXd::Identity(4, 4);
  I_s(2, 2) = 0.5;

  // cal nomal strain
  VectorXd strain_n(4);
  strain_n[0] = strain[0];
  strain_n[1] = strain[1];
  strain_n[2] = 0.5 * strain[2];
  strain_n[3] = 0.0;
  // trace of strain
  double tr_strain = strain_n[0] + strain_n[1] + strain_n[3];

  // deviation strain
  VectorXd strain_n1(4);
  strain_n1[0] = strain_n[0] - (1.0 / 3.0) * tr_strain;
  strain_n1[1] = strain_n[1] - (1.0 / 3.0) * tr_strain;
  strain_n1[2] = strain_n[2];
  strain_n1[3] = strain_n[3] - (1.0 / 3.0) * tr_strain;

  // deviation trial stress
  VectorXd stress_tri1(4);
  stress_tri1[0] = 2.0 * G * (strain_n1[0] - strain_p_[ip][0]);
  stress_tri1[1] = 2.0 * G * (strain_n1[1] - strain_p_[ip][1]);
  stress_tri1[2] = 2.0 * G * (strain_n1[2] - strain_p_[ip][2]);
  stress_tri1[3] = -(stress_tri1[0] + stress_tri1[1]);

  // tryal xi
  VectorXd xi_tri(4);
  for(int i = 0; i < 4; i++)
    xi_tri[i] = stress_tri1[i] - beta_[ip][i];

  //  equivalent stress */
  double xi_tri_norm = 0.0;
  xi_tri_norm = xi_tri[0] * xi_tri[0] + xi_tri[1] * xi_tri[1] +
                xi_tri[3] * xi_tri[3] + 2.0 * xi_tri[2] * xi_tri[2];
  xi_tri_norm = sqrt(xi_tri_norm);

  // cal yield condition
  double f =
      xi_tri_norm - sqrt(2.0 / 3.0) * (sigma_y_ + Ehi * strain_p_bar_[ip]);

  // Not yield
  if(f <= 0.0)
  {
    // final stress */
    stress[0] = K * tr_strain + xi_tri[0];
    stress[1] = K * tr_strain + xi_tri[1];
    stress[2] = xi_tri[2];

    // make Ce */
    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++)
        De_(i, j) =
            K * I[i] * I[j] + 2.0 * G * (I_s(i, j) - (1.0 / 3.0) * I[i] * I[j]);
  }
  // Yield Plastic
  else
  {
    // delta gamma
    double gamma = f / (2.0 * G + (2.0 / 3.0) * hards);

    if(gamma < 0.0)
    {
      cerr << "gamma is negative in mat_plasto_simo.cpp" << endl;
      exit(1);
    }

    VectorXd n(3);
    for(int i = 0; i < 3; i++)
      n[i] = xi_tri[i] / xi_tri_norm;

    // final stress
    stress[0] = K * tr_strain + xi_tri[0] - 2.0 * G * gamma * n[0];
    stress[1] = K * tr_strain + xi_tri[1] - 2.0 * G * gamma * n[1];
    stress[2] = xi_tri[2] - 2.0 * G * gamma * n[2];

    // theta
    double theta = 1.0 - (2.0 * G * gamma) / xi_tri_norm;

    // theta bar
    double theta1 = (1.0 / (1.0 + (hards / (3.0 * G)))) - (1.0 - theta);

    // make Cep
    for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
      {
        De_(i, j) = K * I[i] * I[j] +
                    2.0 * G * theta * (I_s(i, j) - (1.0 / 3.0) * I[i] * I[j]) -
                    2.0 * G * theta1 * n[i] * n[j];
      }
    }
  }

  return;
}

void ElastPlast::plasto_simo_3d(const VectorXd &strain, VectorXd &stress,
                                int ip)
{

  SIMPEP();

  double hards = (young_ * Eh_) / (young_ - Eh_);
  double Ehi = betaH_ * hards;
  double G = (0.5 * young_) / (1.0 + poisson_);
  double K = young_ / (3.0 * (1.0 - 2.0 * poisson_));

  // A matrix
  VectorXd I(6);
  I << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;

  // I matrix
  MatrixXd I_s = MatrixXd::Identity(6, 6);
  I_s(3, 3) = 0.5;
  I_s(4, 4) = 0.5;
  I_s(5, 5) = 0.5;

  // cal nomal strain
  VectorXd strain_n(6);
  strain_n[0] = strain[0];
  strain_n[1] = strain[1];
  strain_n[2] = strain[2];
  strain_n[3] = 0.5 * strain[3];
  strain_n[4] = 0.5 * strain[4];
  strain_n[5] = 0.5 * strain[5];

  // trace of strain
  double tr_strain = strain_n[0] + strain_n[1] + strain_n[2];

  // deviation strain
  VectorXd strain_n1(6);
  strain_n1[0] = strain_n[0] - (1.0 / 3.0) * tr_strain;
  strain_n1[1] = strain_n[1] - (1.0 / 3.0) * tr_strain;
  strain_n1[2] = strain_n[2] - (1.0 / 3.0) * tr_strain;
  strain_n1[3] = strain_n[3];
  strain_n1[4] = strain_n[4];
  strain_n1[5] = strain_n[5];

  // deviation trial stress(size6)
  VectorXd stress_tri1 = 2.0 * G * (strain_n1 - strain_p_[ip]);

  // tryal xi(size6)
  VectorXd xi_tri = stress_tri1 - beta_[ip];

  //  equivalent stress */
  double xi_tri_norm = 0.0;
  xi_tri_norm = xi_tri[0] * xi_tri[0] + xi_tri[1] * xi_tri[1] +
                xi_tri[2] * xi_tri[2] +
                2.0 * (xi_tri[3] * xi_tri[3] + xi_tri[4] * xi_tri[4] +
                       xi_tri[5] * xi_tri[5]);
  xi_tri_norm = sqrt(xi_tri_norm);

  // cal yield condition
  double f =
      xi_tri_norm - sqrt(2.0 / 3.0) * (sigma_y_ + Ehi * strain_p_bar_[ip]);

  // Not yield
  if(f <= 0.0)
  {
    // final stress */
    stress[0] = K * tr_strain + stress_tri1[0];
    stress[1] = K * tr_strain + stress_tri1[1];
    stress[2] = K * tr_strain + stress_tri1[2];
    stress[3] = stress_tri1[3];
    stress[4] = stress_tri1[4];
    stress[5] = stress_tri1[5];

    // make Ce */
    for(int i = 0; i < 6; i++)
      for(int j = 0; j < 6; j++)
        De_(i, j) =
            K * I[i] * I[j] + 2.0 * G * (I_s(i, j) - (1.0 / 3.0) * I[i] * I[j]);
  }
  // Yield Plastic
  else
  {
    // delta gamma
    double gamma = f / (2.0 * G + (2.0 / 3.0) * hards);

    if(gamma < 0.0)
    {
      cerr << "gamma is negative in mat_plasto_simo.cpp" << endl;
      exit(1);
    }

    VectorXd n = xi_tri / xi_tri_norm;

    // final stress
    stress[0] = K * tr_strain + stress_tri1[0] - 2.0 * G * gamma * n[0];
    stress[1] = K * tr_strain + stress_tri1[1] - 2.0 * G * gamma * n[1];
    stress[2] = K * tr_strain + stress_tri1[2] - 2.0 * G * gamma * n[2];
    stress[3] = stress_tri1[3] - 2.0 * G * gamma * n[3];
    stress[4] = stress_tri1[4] - 2.0 * G * gamma * n[4];
    stress[5] = stress_tri1[5] - 2.0 * G * gamma * n[5];

    // theta
    double theta = 1.0 - (2.0 * G * gamma) / xi_tri_norm;

    // theta bar
    double theta1 = (1.0 / (1.0 + (hards / (3.0 * G)))) - (1.0 - theta);

    // make Cep
    for(int i = 0; i < 6; i++)
    {
      for(int j = 0; j < 6; j++)
      {
        De_(i, j) = K * I[i] * I[j] +
                    2.0 * G * theta * (I_s(i, j) - (1.0 / 3.0) * I[i] * I[j]) -
                    2.0 * G * theta1 * n[i] * n[j];
      }
    }
  }

  return;
}

void ElastPlast::stress_plastic_plane_strain(
    FEM &fem, std::vector<Node> &nodes,
    std::vector<Eigen::Triplet<double>> &tripletsR)
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

  SIMPEP();

  double hards = (young_ * Eh_) / (young_ - Eh_);
  double Ehi = betaH_ * hards;
  double Ehk = (1.0 - betaH_) * hards;
  double G = (0.5 * young_) / (1.0 + poisson_);
  double K = young_ / (3.0 * (1.0 - 2.0 * poisson_));

  // init
  Bmatrix bmatrix;
  double jac = 0.0;
  double fac = 1.0 / ipmax; //<factor for element-wise value
  eqstn_ = 0.0;
  mises_ = 0.0;
  tr_stress_ = 0.0;
  fe_.setZero();
  stress_.setZero();
  strain_.setZero();

  // gauss loop start
  for(int ip = 0; ip < ipmax; ip++)
  {
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);

    VectorXd strain = Be * disp;

    // normal strain
    VectorXd strain_n(4);
    strain_n[0] = strain[0];
    strain_n[1] = strain[1];
    strain_n[2] = 0.5 * strain[2];
    strain_n[3] = 0.0;

    // trace of strain
    double tr_strain = strain_n[0] + strain_n[1] + strain_n[3];

    // deviation strain
    VectorXd strain_n1(4);
    strain_n1[0] = strain_n[0] - (1.0 / 3.0) * tr_strain;
    strain_n1[1] = strain_n[1] - (1.0 / 3.0) * tr_strain;
    strain_n1[2] = strain_n[2];
    strain_n1[3] = strain_n[3] - (1.0 / 3.0) * tr_strain;

    // cal deviation trial stress
    VectorXd stress_tri1(4);
    stress_tri1[0] = 2.0 * G * (strain_n1[0] - strain_p_[ip][0]);
    stress_tri1[1] = 2.0 * G * (strain_n1[1] - strain_p_[ip][1]);
    stress_tri1[2] = 2.0 * G * (strain_n1[2] - strain_p_[ip][2]);
    stress_tri1[3] = -(stress_tri1[0] + stress_tri1[1]);

    // trial xi
    xi_tri_[ip].setZero();
    for(int i = 0; i < 4; i++)
      xi_tri_[ip][i] = stress_tri1[i] - beta_[ip][i];

    // equivalent stress
    xi_tri_norm_[ip] =
        xi_tri_[ip][0] * xi_tri_[ip][0] + xi_tri_[ip][1] * xi_tri_[ip][1] +
        xi_tri_[ip][3] * xi_tri_[ip][3] + 2.0 * xi_tri_[ip][2] * xi_tri_[ip][2];
    xi_tri_norm_[ip] = sqrt(xi_tri_norm_[ip]);

    // cal yield condition
    yield_[ip] = xi_tri_norm_[ip] -
                 sqrt(2.0 / 3.0) * (sigma_y_ + Ehi * strain_p_bar_[ip]);

    VectorXd stress(fem.voigt);
    stress.setZero();
    // start of yield criterion
    // yield Elastic
    if(yield_[ip] <= 0.0)
    {
      // final stress
      stress[0] = K * tr_strain + xi_tri_[ip][0];
      stress[1] = K * tr_strain + xi_tri_[ip][1];
      stress[2] = xi_tri_[ip][2];
      sig3_[ip] = poisson_ * (stress[0] + stress[1]);
      gamma_[ip] = 0.0;
    }
    // Yield Plastic
    else
    {
      // delta gamma
      gamma_[ip] = yield_[ip] / (2.0 * G + (2.0 / 3.0) * hards);

      VectorXd n(4);
      for(int i = 0; i < 4; i++)
        n[i] = xi_tri_[ip][i] / xi_tri_norm_[ip];

      // equivalent plastic strain
      strain_p_bar_[ip] += sqrt(2.0 / 3.0) * gamma_[ip];
      epsilon_ += sqrt(2.0 / 3.0) * gamma_[ip] * fac;

      // update back stress,plastic strain and stress
      // beta
      for(int i = 0; i < 4; i++)
        beta_[ip][i] += (2.0 / 3.0) * Ehk * gamma_[ip] * n[i];

      // plastic strain
      for(int i = 0; i < 4; i++)
        strain_p_[ip][i] += gamma_[ip] * n[i];
      //  final stress
      stress[0] = K * tr_strain + xi_tri_[ip][0] - 2.0 * G * gamma_[ip] * n[0];
      stress[1] = K * tr_strain + xi_tri_[ip][1] - 2.0 * G * gamma_[ip] * n[1];
      stress[2] = xi_tri_[ip][2] - 2.0 * G * gamma_[ip] * n[2];
      sig3_[ip] = poisson_ * (stress[0] + stress[1]);
    }
    // end of yield criterion

    // deviation strain & stress */
    double tr_stress = stress[0] + stress[1] + sig3_[ip];

    VectorXd strain1(3);
    strain1[0] = strain[0] - (1.0 / 3.0) * tr_strain;
    strain1[1] = strain[1] - (1.0 / 3.0) * tr_strain;
    strain1[2] = strain[2];

    VectorXd stress1(4);
    stress1[0] = stress[0] - (1.0 / 3.0) * tr_stress;
    stress1[1] = stress[1] - (1.0 / 3.0) * tr_stress;
    stress1[2] = stress[2];
    stress1[3] = -(stress1[0] + stress1[1]);

    // cal stress
    tr_stress_ += (1.0 / 3.0) * tr_stress * fac;

    // cal Mises stress and equivalent strain
    double mises = stress1[0] * stress1[0] + stress1[1] * stress1[1] +
                   stress1[3] * stress1[3] + 2.0 * stress1[2] * stress1[2];

    double eqstn = strain1[0] * strain1[0] + strain[1] * strain[1] +
                   2.0 * strain[2] * strain[2];

    stress_fin_bar_[ip] = sqrt((3.0 / 2.0) * mises);
    mises_ += stress_fin_bar_[ip] * fac;
    eqstn_ += sqrt((3.0 / 2.0) * eqstn) * fac;

    // stock strain & stress
    strain_gp_[ip] = strain;
    stress_gp_[ip] = stress;
    for(int i = 0; i < 3; i++)
    {
      stress_[i] += stress[i] * fac;
      strain_[i] += strain[i] * fac;
    }

    fe_ += Be.transpose() * stress * jac;

  } // gauss loop end

  assemblingReact(numdof, idof, fe_, tripletsR);

  return;
}

void ElastPlast::stress_plastic_3d(
    FEM &fem, std::vector<Node> &nodes,
    std::vector<Eigen::Triplet<double>> &tripletsR)
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

  SIMPEP();

  double hards = (young_ * Eh_) / (young_ - Eh_);
  double Ehi = betaH_ * hards;
  double Ehk = (1.0 - betaH_) * hards;
  double G = (0.5 * young_) / (1.0 + poisson_);
  double K = young_ / (3.0 * (1.0 - 2.0 * poisson_));

  // init
  Bmatrix bmatrix;
  double jac = 0.0;
  double fac = 1.0 / ipmax; //<factor for element-wise value
  eqstn_ = 0.0;
  mises_ = 0.0;
  tr_stress_ = 0.0;
  fe_.setZero();
  stress_.setZero();
  strain_.setZero();

  // gauss loop start
  for(int ip = 0; ip < ipmax; ip++)
  {
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);

    VectorXd strain = Be * disp;

    // normal strain
    VectorXd strain_n(6);
    strain_n[0] = strain[0];
    strain_n[1] = strain[1];
    strain_n[2] = strain[2];
    strain_n[3] = 0.5 * strain[3];
    strain_n[4] = 0.5 * strain[4];
    strain_n[5] = 0.5 * strain[5];

    // trace of strain
    double tr_strain = strain_n[0] + strain_n[1] + strain_n[2];

    // deviation strain
    VectorXd strain_n1(6);
    strain_n1[0] = strain_n[0] - (1.0 / 3.0) * tr_strain;
    strain_n1[1] = strain_n[1] - (1.0 / 3.0) * tr_strain;
    strain_n1[2] = strain_n[2] - (1.0 / 3.0) * tr_strain;
    strain_n1[3] = strain_n[3];
    strain_n1[4] = strain_n[4];
    strain_n1[5] = strain_n[5];

    // cal deviation trial stress(size6)
    VectorXd stress_tri1 = 2.0 * G * (strain_n1 - strain_p_[ip]);

    // trial xi
    xi_tri_[ip] = stress_tri1 - beta_[ip];

    // equivalent stress
    for(int i = 0; i < 3; i++)
      xi_tri_norm_[ip] += xi_tri_[ip][i] * xi_tri_[ip][i];
    for(int i = 3; i < 6; i++)
      xi_tri_norm_[ip] += 2.0 * xi_tri_[ip][i] * xi_tri_[ip][i];

    xi_tri_norm_[ip] = sqrt(xi_tri_norm_[ip]);

    // cal yield condition
    yield_[ip] = xi_tri_norm_[ip] -
                 sqrt(2.0 / 3.0) * (sigma_y_ + Ehi * strain_p_bar_[ip]);

    VectorXd stress(fem.voigt);
    stress.setZero();
    // start of yield criterion
    // yield Elastic
    if(yield_[ip] <= 0.0)
    {
      // final stress
      stress[0] = K * tr_strain + stress_tri1[0];
      stress[1] = K * tr_strain + stress_tri1[1];
      stress[2] = K * tr_strain + stress_tri1[2];
      stress[3] = stress_tri1[3];
      stress[4] = stress_tri1[4];
      stress[5] = stress_tri1[5];
      gamma_[ip] = 0.0;
    }
    // Yield Plastic
    else
    {
      // delta gamma
      gamma_[ip] = yield_[ip] / (2.0 * G + (2.0 / 3.0) * hards);

      VectorXd n = xi_tri_[ip] / xi_tri_norm_[ip];

      // equivalent plastic strain
      strain_p_bar_[ip] += sqrt(2.0 / 3.0) * gamma_[ip];
      epsilon_ += sqrt(2.0 / 3.0) * gamma_[ip] * fac;

      // update back stress,plastic strain and stress
      // beta
      beta_[ip] += (2.0 / 3.0) * Ehk * gamma_[ip] * n;

      // plastic strain
      strain_p_[ip] += gamma_[ip] * n;

      //  final stress
      stress[0] = K * tr_strain + stress_tri1[0] - 2.0 * G * gamma_[ip] * n[0];
      stress[1] = K * tr_strain + stress_tri1[1] - 2.0 * G * gamma_[ip] * n[1];
      stress[2] = K * tr_strain + stress_tri1[2] - 2.0 * G * gamma_[ip] * n[2];
      stress[3] = stress_tri1[3] - 2.0 * G * gamma_[ip] * n[3];
      stress[4] = stress_tri1[4] - 2.0 * G * gamma_[ip] * n[4];
      stress[5] = stress_tri1[5] - 2.0 * G * gamma_[ip] * n[5];
    }
    // end of yield criterion

    tr_strain = 0.0;
    double tr_stress = 0.0;
    VectorXd strain1(6), stress1(6);
    for(int i = 0; i < 3; i++)
    {
      tr_strain += strain[i];
      tr_stress += stress[i];
    }
    for(int i = 0; i < 3; i++)
    {
      strain1[i] = strain[i] - (1.0 / 3.0) * tr_strain;
      stress1[i] = stress[i] - (1.0 / 3.0) * tr_stress;
    }
    for(int i = 3; i < 6; i++)
    {
      strain1[i] = strain[i];
      stress1[i] = stress[i];
    }

    tr_stress_ += (1.0 / 3.0) * tr_stress * fac;

    // cal Mises stress and equivalent strain
    double mises = 0.0, eqstn = 0.0;
    for(int i = 0; i < 3; i++)
    {
      mises += stress1[i] * stress1[i];
      eqstn += strain1[i] * strain1[i];
    }
    for(int i = 3; i < 6; i++)
    {
      mises += 2.0 * stress1[i] * stress1[i];
      eqstn += 2.0 * strain1[i] * strain1[i];
    }

    stress_fin_bar_[ip] = sqrt((3.0 / 2.0) * mises);
    mises_ += stress_fin_bar_[ip] * fac;
    eqstn_ += sqrt((3.0 / 2.0) * eqstn) * fac;

    // stock strain & stress
    strain_gp_[ip] = strain;
    stress_gp_[ip] = stress;

    stress_ += stress * fac;
    strain_ += strain * fac;

    fe_ += Be.transpose() * stress * jac;

  } // gauss loop end

  assemblingReact(numdof, idof, fe_, tripletsR);

  return;
}

void ElastPlast::elasto(const VectorXd &strain, VectorXd &stress, int ip)
{
  De(De_, young_, poisson_);

  stress += De_ * strain;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

void LinearElastPlast::makeKeFint(vector<Triplet<double>> &tripletsK,
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

  SIMPEP();

  // init
  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;
  Bmatrix bmatrix;

  iteration++;

  // gauss loop
  for(int ip = 0; ip < ipmax; ip++)
  {
    double jac = 0.0;
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);

    VectorXd strain = Be * disp;
    VectorXd stress = VectorXd::Zero(fem.voigt);

    // make De_ & stress
    if(iteration == 1) // The first step of Newton-Raphson method is computed by
                       // elastic constitutive law
    {
      elasto(strain, stress, ip);
    }
    else
    {
      plasto_simo(strain, stress, ip);
    }

    // inner force
    fe_ += Be.transpose() * stress * jac;
    Ke_ += Be.transpose() * De_ * Be * jac;
    volume += jac;

  } // gauss point loop end

  assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
  assemblingReact(numdof, idof, fe_, tripletsR);
}

} // namespace icarat