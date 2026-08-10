///  @file  elasto_plastic_mises.cpp
///  @author  Daiki Watanabe
///  @date  October 3, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "subloading-analytical_iso.hpp"
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
  Re = toml::find<double>(mat, "Re");
  alpha_inter_ = toml::find<double>(mat, "alpha");
  beta_inter_ = toml::find<double>(mat, "betainter");
  us_ = toml::find<double>(mat, "u");

  using pair = std::array<double, 2>;
  pair input;
  input = toml::find<pair>(mat, "young");
  young1_ = get<0>(input);
  young2_ = get<1>(input);
  Ehs_ = toml::find<vector<double>>(mat, "h");
  sigma_ys_ = toml::find<vector<double>>(mat, "y0");
  // us_ = toml::find<vector<double>>(mat, "u");

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
  mises_gauss_ = 0.0;
  eqstn_gauss_ = 0.0;
  for(int ip = 0; ip < ipmax; ip++)
  {
    strain_p_[ip].setZero();
    beta_[ip].setZero();
    strain_p_bar_[ip] = 0.0;
    R_[ip] = Re;
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
  if(Ehs_[1] >= Ehs_[0])
  Eh_ = pow(design_s_, pp_) * Ehs_[1] + (1.0 - pow(design_s_, pp_)) * Ehs_[0];
  else
  Eh_ = (1.0 - pow(1.0 - design_s_, pp_)) * Ehs_[1] + pow(1.0 - design_s_, pp_) * Ehs_[0];
  if(sigma_ys_[1] >= sigma_ys_[0])
  sigma_y_ = pow(design_s_, pp_) * sigma_ys_[1] + (1.0 - pow(design_s_, pp_)) * sigma_ys_[0];
  else
  sigma_y_ = (1.0 - pow(1.0 - design_s_, pp_)) * sigma_ys_[1] + pow(1.0 - design_s_, pp_) * sigma_ys_[0];
  //u_ = pow(design_s_, pp_) * us_[1] + (1.0 - pow(design_s_, pp_)) * us_[0]; ///< SIMP
  u_ = us_ - alpha_inter_ * exp(- (design_s_ - 0.5) * (design_s_ - 0.5) / beta_inter_); ///< new interporation of u
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

  // I matrix for strain
  MatrixXd I_s2 = MatrixXd::Identity(4, 4);
  I_s2(2, 2) = 2.0;

  //deviation tensor
  MatrixXd I_d = MatrixXd::Identity(4, 4);
  for(int i = 0; i < 4; i++)
  {
    for(int j = 0; j < 4; j++)
    {
      I_d(i, j) = I_s2(i, j) - (1.0 / 3.0) * I[i] * I[j];
    }
  }

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
  stress_tri1[3] = 2.0 * G * (strain_n1[3] - strain_p_[ip][3]);

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
      xi_tri_norm - sqrt(2.0 / 3.0) * R_[ip] * (sigma_y_ + Ehi * strain_p_bar_[ip]);

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
    double gamma = -1.0;
    double R = R_[ip];

    double A1 = (2.0 / 3.0) * Ehi;
    double A2 = 2.0 * G - (2.0 / 3.0) * Ehi * Re - (2.0 / 3.0) * Ehi * R - sqrt(2.0 / 3.0) * u_ * (sigma_y_ + Ehi * strain_p_bar_[ip]);
    double A3 = - 2.0 * G * (Re + R) + (2.0 / 3.0) * R * Re * Ehi + u_ * xi_tri_norm + sqrt(2.0 / 3.0) * u_ * (sigma_y_ + Ehi * strain_p_bar_[ip]);
    double A4 = 2.0 * G * R * Re - u_ * xi_tri_norm;

    std::complex<double> solR[3] = {};
    VectorXd gamma_temp(3);
    VectorXd gamma_den(3);
    VectorXd gamma_num(3);

    gamma_num = VectorXd::Zero(3);
    gamma_den = VectorXd::Zero(3);
    gamma_temp = VectorXd::Zero(3);

    solve_cubic(solR, A1, A2, A3, A4);
    for(int i = 0; i < 4; i++)
    {
      int count_i = 0;
      if(i == 3)
      {
        cerr << "No satisgying R in mat_plasto_simo.cpp" << endl;
        cerr << count_i << endl;
        cerr << solR[0] << solR[1] << solR[2] << endl;
        cerr << xi_tri_norm << endl;
        cerr << strain_p_bar_[ip] << endl;
        cerr << sigma_y_ << endl;
        cerr << "Ehi:" << Ehi << endl;
        cerr << "G:" << G << endl;
        cerr << "Re:" << Re << endl;
        cerr << gamma_temp[0] << " " << gamma_temp[1] << " " << gamma_temp[2] << endl;
        cerr << gamma_den[0] << " " << gamma_den[1] << " " << gamma_den[2] << endl;
        cerr << gamma_num[0] << " " << gamma_num[1] << " " << gamma_num[2] << endl;
        exit(1);
      }
      double realR = solR[i].real();
      double imagR = solR[i].imag();
      if (realR > 1.0 && realR < 1.0 + 1e-4) 
      {
        realR = 1.0;
      }
      if (realR < Re && realR > Re - 1e-4)
      {
        realR = Re;
      }
      if(std::abs(imagR) < 1e-10 && (Re <= realR) && (realR <= 1.0))
      {
        count_i += 1;
        R = realR;
        gamma_den[i] = xi_tri_norm - sqrt(2.0 / 3.0) * R * sigma_y_ - sqrt(2.0 / 3.0) * R * Ehi * strain_p_bar_[ip];
        gamma_num[i] = 2.0 * G + (2.0 / 3.0) * R * Ehi;
        gamma_temp[i] = gamma_den[i] / gamma_num[i];
        if(gamma_temp[i] >= 0.0)
        {
          gamma = gamma_temp[i];
          break;
        }
        if(gamma_temp[i] < 0.0 && gamma_temp[i] > -1e-8)
        {
          gamma = 0.0;
          break;
        }
      }
    }

    if(gamma < 0.0)
    {
      cerr << "gamma is negative in mat_plasto_simo.cpp" << endl;
      cerr << solR[0] << solR[1] << solR[2] << endl;
      cerr << xi_tri_norm << endl;
      cerr << strain_p_bar_[ip] << endl;
      cerr << sigma_y_ << endl;
      cerr << "Ehi:" << Ehi << endl;
      cerr << "G:" << G << endl;
      cerr << gamma_temp[0] << " " << gamma_temp[1] << " " << gamma_temp[2] << endl;
      cerr << gamma_den[0] << " " << gamma_den[1] << " " << gamma_den[2] << endl;
      cerr << gamma_num[0] << " " << gamma_num[1] << " " << gamma_num[2] << endl;
      cerr << gamma << endl;
      cerr << xi_tri_norm - sqrt(2.0 / 3.0) * R * (sigma_y_ + Ehi * strain_p_bar_[ip]) << endl;
      exit(1);
    }

    // equivalent plastic strain
    double strain_p_bar = strain_p_bar_[ip] + sqrt(2.0 / 3.0) * gamma;

    VectorXd n = xi_tri / xi_tri_norm;

    // final stress
    stress[0] = K * tr_strain + stress_tri1[0] - 2.0 * G * gamma * n[0];
    stress[1] = K * tr_strain + stress_tri1[1] - 2.0 * G * gamma * n[1];
    stress[2] = stress_tri1[2] - 2.0 * G * gamma * n[2];

    //Consistent tangent modulus
    n[2] = 2.0 * xi_tri[2] / xi_tri_norm;

    MatrixXd De_inv = MatrixXd::Zero(3, 3);
    De(De_, young_, poisson_);
      double U_R = - u_ * (R - 1.0) / (R - Re);
      double dU_dR = - u_ * (1.0 - Re) / ((R - Re) * (R - Re));

      for(int i = 0; i < 3; i++)
      {
        for(int j = 0; j < 3; j++)
        {
          De_inv(i, j) = De_.inverse()(i, j) + gamma * I_d(i, j) / xi_tri_norm + (1.0 / ((2.0 / 3.0) * R * Ehi + sqrt(2.0 / 3.0) * (sigma_y_ + Ehi * strain_p_bar) * U_R / (1.0 - gamma * dU_dR)) - gamma / xi_tri_norm) * (n[i] * n[j]);
        }
      }
    De_ = De_inv.inverse();

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

  // I matrix for strain
  MatrixXd I_s2 = MatrixXd::Identity(6, 6);
  I_s2(3, 3) = 2.0;
  I_s2(4, 4) = 2.0;
  I_s2(5, 5) = 2.0;
  
  //deviation tensor
  MatrixXd I_d = MatrixXd::Identity(6, 6);
  for(int i = 0; i < 6; i++)
  {
    for(int j = 0; j < 6; j++)
    {
      I_d(i, j) = I_s2(i, j) - (1.0 / 3.0) * I[i] * I[j];
    }
  }

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
      xi_tri_norm - sqrt(2.0 / 3.0) * R_[ip] * (sigma_y_ + Ehi * strain_p_bar_[ip]);
  
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
    double gamma = -1.0;
    double R = R_[ip];

    double A1 = (2.0 / 3.0) * Ehi;
    double A2 = 2.0 * G - (2.0 / 3.0) * Ehi * Re - (2.0 / 3.0) * Ehi * R - sqrt(2.0 / 3.0) * u_ * (sigma_y_ + Ehi * strain_p_bar_[ip]);
    double A3 = - 2.0 * G * (Re + R) + (2.0 / 3.0) * R * Re * Ehi + u_ * xi_tri_norm + sqrt(2.0 / 3.0) * u_ * (sigma_y_ + Ehi * strain_p_bar_[ip]);
    double A4 = 2.0 * G * R * Re - u_ * xi_tri_norm;

    std::complex<double> solR[3] = {};
    VectorXd gamma_temp(3);
    VectorXd gamma_den(3);
    VectorXd gamma_num(3);

    gamma_num = VectorXd::Zero(3);
    gamma_den = VectorXd::Zero(3);
    gamma_temp = VectorXd::Zero(3);

    solve_cubic(solR, A1, A2, A3, A4);
    for(int i = 0; i < 4; i++)
    {
      int count_i = 0;
      if(i == 3)
      {
        cerr << "No satisgying R in mat_plasto_simo.cpp" << endl;
        cerr << count_i << endl;
        cerr << solR[0] << solR[1] << solR[2] << endl;
        cerr << xi_tri_norm << endl;
        cerr << strain_p_bar_[ip] << endl;
        cerr << sigma_y_ << endl;
        cerr << "Ehi:" << Ehi << endl;
        cerr << "G:" << G << endl;
        cerr << "Re:" << Re << endl;
        cerr << gamma_temp[0] << " " << gamma_temp[1] << " " << gamma_temp[2] << endl;
        cerr << gamma_den[0] << " " << gamma_den[1] << " " << gamma_den[2] << endl;
        cerr << gamma_num[0] << " " << gamma_num[1] << " " << gamma_num[2] << endl;
        exit(1);
      }
      double realR = solR[i].real();
      double imagR = solR[i].imag();
      if (realR > 1.0 && realR < 1.0 + 1e-4) 
      {
        realR = 1.0;
      }
      if (realR < Re && realR > Re - 1e-4)
      {
        realR = Re;
      }
      if(std::abs(imagR) < 1e-10 && (Re <= realR) && (realR <= 1.0))
      {
        count_i += 1;
        R = realR;
        gamma_den[i] = xi_tri_norm - sqrt(2.0 / 3.0) * R * sigma_y_ - sqrt(2.0 / 3.0) * R * Ehi * strain_p_bar_[ip];
        gamma_num[i] = 2.0 * G + (2.0 / 3.0) * R * Ehi;
        gamma_temp[i] = gamma_den[i] / gamma_num[i];
        if(gamma_temp[i] >= 0.0)
        {
          gamma = gamma_temp[i];
          break;
        }
        if(gamma_temp[i] < 0.0 && gamma_temp[i] > -1e-8)
        {
          gamma = 0.0;
          break;
        }
      }
    }

    if(gamma < 0.0)
    {
      cerr << "gamma is negative in mat_plasto_simo.cpp" << endl;
      cerr << solR[0] << solR[1] << solR[2] << endl;
      cerr << xi_tri_norm << endl;
      cerr << strain_p_bar_[ip] << endl;
      cerr << sigma_y_ << endl;
      cerr << "Ehi:" << Ehi << endl;
      cerr << "G:" << G << endl;
      cerr << gamma_temp[0] << " " << gamma_temp[1] << " " << gamma_temp[2] << endl;
      cerr << gamma_den[0] << " " << gamma_den[1] << " " << gamma_den[2] << endl;
      cerr << gamma_num[0] << " " << gamma_num[1] << " " << gamma_num[2] << endl;
      cerr << gamma << endl;
      cerr << xi_tri_norm - sqrt(2.0 / 3.0) * R * (sigma_y_ + Ehi * strain_p_bar_[ip]) << endl;
      exit(1);
    }

    // equivalent plastic strain
    double strain_p_bar = strain_p_bar_[ip] + sqrt(2.0 / 3.0) * gamma;

    VectorXd n = xi_tri / xi_tri_norm;

    // final stress
    stress[0] = K * tr_strain + stress_tri1[0] - 2.0 * G * gamma * n[0];
    stress[1] = K * tr_strain + stress_tri1[1] - 2.0 * G * gamma * n[1];
    stress[2] = K * tr_strain + stress_tri1[2] - 2.0 * G * gamma * n[2];
    stress[3] = stress_tri1[3] - 2.0 * G * gamma * n[3];
    stress[4] = stress_tri1[4] - 2.0 * G * gamma * n[4];
    stress[5] = stress_tri1[5] - 2.0 * G * gamma * n[5];

    //Consistent tangent modulus
    for(int i = 0; i < 3; i++)
      n[i] = xi_tri[i] / xi_tri_norm;
    for(int i = 3; i < 6; i++)
      n[i] = 2.0 * xi_tri[i] / xi_tri_norm;

    MatrixXd De_inv = MatrixXd::Zero(6, 6);
    De(De_, young_, poisson_);
      double U_R = - u_ * (R - 1.0) / (R - Re);
      double dU_dR = - u_ * (1.0 - Re) / ((R - Re) * (R - Re));

      for(int i = 0; i < 6; i++)
      {
        for(int j = 0; j < 6; j++)
        {
          De_inv(i, j) = De_.inverse()(i, j) + gamma * I_d(i, j) / xi_tri_norm + (1.0 / ((2.0 / 3.0) * R * Ehi + sqrt(2.0 / 3.0) * (sigma_y_ + Ehi * strain_p_bar) * U_R / (1.0 - gamma * dU_dR)) - gamma / xi_tri_norm) * (n[i] * n[j]);
        }
      }
    De_ = De_inv.inverse();
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
    stress_tri1[3] = 2.0 * G * (strain_n1[3] - strain_p_[ip][3]);

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
                 sqrt(2.0 / 3.0) * R_[ip] * (sigma_y_ + Ehi * strain_p_bar_[ip]);

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
      sig3_[ip] = K * tr_strain + xi_tri_[ip][3];
      gamma_[ip] = 0.0;
    }
    // Yield Plastic
    else
    {
      gamma_[ip] = -1.0;

      double A1 = (2.0 / 3.0) * Ehi;
      double A2 = 2.0 * G - (2.0 / 3.0) * Ehi * Re - (2.0 / 3.0) * Ehi * R_[ip] - sqrt(2.0 / 3.0) * u_ * (sigma_y_ + Ehi * strain_p_bar_[ip]);
      double A3 = - 2.0 * G * (Re + R_[ip]) + (2.0 / 3.0) * R_[ip] * Re * Ehi + u_ * xi_tri_norm_[ip] + sqrt(2.0 / 3.0) * u_ * (sigma_y_ + Ehi * strain_p_bar_[ip]);
      double A4 = 2.0 * G * R_[ip] * Re - u_ * xi_tri_norm_[ip];

      std::complex<double> solR[3] = {};

      solve_cubic(solR, A1, A2, A3, A4);
      for(int i = 0; i < 3; i++)
      {
        double realR = solR[i].real();
        double imagR = solR[i].imag();
        if (realR > 1.0 && realR < 1.0 + 1e-4) 
        {
          realR = 1.0;
        }
        if (realR < Re && realR > Re - 1e-4)
        {
          realR = Re;
        }
        if(std::abs(imagR) < 1e-10 && (Re <= realR) && (realR <= 1.0))
        {
          R_[ip] = realR;
          double gamma_den = xi_tri_norm_[ip] - sqrt(2.0 / 3.0) * R_[ip] * sigma_y_ - sqrt(2.0 / 3.0) * R_[ip] * Ehi * strain_p_bar_[ip];
          double gamma_num = 2.0 * G + (2.0 / 3.0) * R_[ip] * Ehi;
          double gamma_temp = gamma_den / gamma_num;
          if(gamma_temp >= 0.0)
          {
            gamma_[ip] = gamma_temp;
            break;
          }
          if(gamma_temp < 0.0 && gamma_temp > -1e-8)
          {
            gamma_[ip] = 0.0;
            break;
          }
        }
      }
      
      if(gamma_[ip] < 0.0)
      {
        cerr << "gamma_[ip] is negative in mat_plasto_simo.cpp" << endl;
        exit(1);
      }

      VectorXd n = VectorXd::Zero(3);
      n = xi_tri_[ip] / xi_tri_norm_[ip];

      // equivalent plastic strain
      strain_p_bar_[ip] += sqrt(2.0 / 3.0) * gamma_[ip];
      epsilon_ += strain_p_bar_[ip] * fac;

      // update back stress,plastic strain and stress
      // beta
      for(int i = 0; i < 4; i++)
        beta_[ip][i] += (2.0 / 3.0) * Ehk * gamma_[ip] * n[i];

      // plastic strain
      for(int i = 0; i < 4; i++)
        strain_p_[ip][i] += gamma_[ip] * n[i];
      
      //  final stress
      stress[0] = K * tr_strain + stress_tri1[0] - 2.0 * G * gamma_[ip] * n[0];
      stress[1] = K * tr_strain + stress_tri1[1] - 2.0 * G * gamma_[ip] * n[1];
      stress[2] = stress_tri1[2] - 2.0 * G * gamma_[ip] * n[2];
      sig3_[ip] = K * tr_strain + stress_tri1[3] - 2.0 * G * gamma_[ip] * n[3];
    }
    // end of yield criterion

    // deviation strain & stress */
    double tr_stress = stress[0] + stress[1] + sig3_[ip];

    VectorXd strain1(4);
    strain1[0] = strain[0] - (1.0 / 3.0) * tr_strain;
    strain1[1] = strain[1] - (1.0 / 3.0) * tr_strain;
    strain1[2] = 0.5 * strain[2];
    strain1[3] = - (1.0 / 3.0) * tr_strain;

    VectorXd stress1(4);
    stress1[0] = stress[0] - (1.0 / 3.0) * tr_stress;
    stress1[1] = stress[1] - (1.0 / 3.0) * tr_stress;
    stress1[2] = stress[2];
    stress1[3] = sig3_[ip] - (1.0 / 3.0) * tr_stress;

    // cal stress
    tr_stress_ += (1.0 / 3.0) * tr_stress * fac;

    // cal Mises stress and equivalent strain
    double mises = stress1[0] * stress1[0] + stress1[1] * stress1[1] +
                   stress1[3] * stress1[3] + 2.0 * stress1[2] * stress1[2];

    double eqstn = strain1[0] * strain1[0] + strain1[1] * strain1[1] +
                   2.0 * strain1[2] * strain1[2] + strain1[3] * strain1[3];

    stress_fin_bar_[ip] = sqrt((3.0 / 2.0) * mises);
    if(ip == 0)
    {
      mises_gauss_ = stress_fin_bar_[ip];
      eqstn_gauss_ = sqrt((3.0 / 2.0) * eqstn);
    }
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

    double temp = 0.0;
    // equivalent stress
    for(int i = 0; i < 3; i++)
      temp += xi_tri_[ip][i] * xi_tri_[ip][i];
    for(int i = 3; i < 6; i++)
      temp += 2.0 * xi_tri_[ip][i] * xi_tri_[ip][i];

    xi_tri_norm_[ip] = sqrt(temp);

    // cal yield condition
    yield_[ip] = xi_tri_norm_[ip] -
                 sqrt(2.0 / 3.0) * R_[ip] * (sigma_y_ + Ehi * strain_p_bar_[ip]); 

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
      //R_[ip] = xi_tri_norm_[ip] / (sqrt(2.0 / 3.0) * (sigma_y_ + Ehi * strain_p_bar_[ip]));
    }
    // Yield Plastic
    else
    {
      gamma_[ip] = -1.0;

      double A1 = (2.0 / 3.0) * Ehi;
      double A2 = 2.0 * G - (2.0 / 3.0) * Ehi * Re - (2.0 / 3.0) * Ehi * R_[ip] - sqrt(2.0 / 3.0) * u_ * (sigma_y_ + Ehi * strain_p_bar_[ip]);
      double A3 = - 2.0 * G * (Re + R_[ip]) + (2.0 / 3.0) * R_[ip] * Re * Ehi + u_ * xi_tri_norm_[ip] + sqrt(2.0 / 3.0) * u_ * (sigma_y_ + Ehi * strain_p_bar_[ip]);
      double A4 = 2.0 * G * R_[ip] * Re - u_ * xi_tri_norm_[ip];

      std::complex<double> solR[3] = {};

      solve_cubic(solR, A1, A2, A3, A4);
      for(int i = 0; i < 3; i++)
      {
        double realR = solR[i].real();
        double imagR = solR[i].imag();
        if (realR > 1.0 && realR < 1.0 + 1e-4) 
        {
          realR = 1.0;
        }
        if (realR < Re && realR > Re - 1e-4)
        {
          realR = Re;
        }
        if(std::abs(imagR) < 1e-10 && (Re <= realR) && (realR <= 1.0))
        {
          R_[ip] = realR;
          double gamma_den = xi_tri_norm_[ip] - sqrt(2.0 / 3.0) * R_[ip] * sigma_y_ - sqrt(2.0 / 3.0) * R_[ip] * Ehi * strain_p_bar_[ip];
          double gamma_num = 2.0 * G + (2.0 / 3.0) * R_[ip] * Ehi;
          double gamma_temp = gamma_den / gamma_num;
          if(gamma_temp >= 0.0)
          {
            gamma_[ip] = gamma_temp;
            break;
          }
          if(gamma_temp < 0.0 && gamma_temp > -1e-8)
          {
            gamma_[ip] = 0.0;
            break;
          }
        }
      }
      
      if(gamma_[ip] < 0.0)
      {
        cerr << "gamma_[ip] is negative in mat_plasto_simo.cpp" << endl;
        exit(1);
      }

      VectorXd n = VectorXd::Zero(6);
      n = xi_tri_[ip] / xi_tri_norm_[ip];

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
      strain1[i] = 0.5 * strain[i];
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
    if(ip == 0)
    {
      // mises_gauss_ = stress_fin_bar_[ip];
      // eqstn_gauss_ = sqrt((3.0 / 2.0) * eqstn);
      mises_gauss_ = stress[1];
      eqstn_gauss_ = strain[1];
    }
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

  stress = De_ * strain;
}

void ElastPlast::solve_cubic(std::complex<double> x[3], double A1, double A2, double A3, double A4)
{
  double B2, B3, B4, p, q, Judge;
  B2 = A2 / A1;
  B3 = A3 / A1;
  B4 = A4 / A1;
  p = B3 - B2 * B2 / 3.0;
  q = 2.0 * B2 * B2 * B2 / 27.0 - B2 * B3 / 3.0 + B4;
  Judge = q * q / 4.0 + p * p * p / 27.0;
  if(Judge < 0.0)
  {
    double theta = atan2(sqrt(- Judge),- q * 0.5);
    x[0] = 2.0 * sqrt(- p / 3.0) * cos(theta / 3.0) - B2 / 3.0;
		x[1] = 2.0 * sqrt(- p / 3.0) * cos((theta + 2.0 * M_PI) / 3.0) - B2 / 3.0;
		x[2] = 2.0 * sqrt(- p / 3.0) * cos((theta + 4.0 * M_PI) / 3.0) - B2 / 3.0;
  }
  else
  {
    double u1 = std::cbrt(- q * 0.5 + sqrt(Judge));
    double v1 = std::cbrt(- q * 0.5 - sqrt(Judge));
		x[0] = u1 + v1 - B2 / 3.0;
		x[1] = - 0.5 * (u1 + v1) + sqrt(3.0) * 0.5 * 1i * (u1 - v1) - B2 / 3.0;
		x[2] = - 0.5 * (u1 + v1) - sqrt(3.0) * 0.5 * 1i * (u1 - v1) - B2 / 3.0;
  }
}

double ElastPlast::Cuberoot(double rootx)
{
  if(rootx > 0.0)
		return pow(rootx, 1.0 / 3.0);
  else
		return -pow(- rootx, 1.0 / 3.0);
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
