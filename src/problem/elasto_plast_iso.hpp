///  @file  elasto_plastic_mises.hpp
///  @author  Daiki Watanabe
///  @date  October 3, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include "elasticity.hpp"

namespace icarat
{
/// infinitesimal plasticity with linear combined hardening material class
/// This class was implemented by referring to the following book.
/// @cite J.C. Simo, Computational Inelasticity.
class ElastPlast : public Elastic
{
public:
  ElastPlast(int voigt_, int ne_, int ipmax_, int numdof_)
      : Elastic(voigt_, ne_, ipmax_, numdof_), Ehs_(2), sigma_ys_(2),
        strain_p_(ipmax_), beta_(ipmax_), strain_p_bar_(ipmax_),
        strain_gp_(ipmax_), stress_gp_(ipmax_), xi_tri_norm_(ipmax_),
        stress_fin_bar_(ipmax_), xi_tri_(ipmax_), gamma_(ipmax_),
        yield_(ipmax_), sig3_(ipmax_)
  {
    for(int i = 0; i < ipmax_; i++)
    {
      strain_p_[i].resize(6); ///< in plane strain
      strain_p_[i].setZero();

      beta_[i].resize(6); ///< in plane strain
      beta_[i].setZero();

      strain_gp_[i].resize(voigt_);
      strain_gp_[i].setZero();

      strain_gp_[i].resize(voigt_);
      strain_gp_[i].setZero();

      xi_tri_[i].resize(6); ///< in plane strain
      xi_tri_[i].setZero();
    }
  }

  /// set material parameter from .dat file
  void setParameter(const toml::value &config);

  /// Initialize private variables
  void initialization();

  void plasto_simo(const Eigen::VectorXd &strain, Eigen::VectorXd &stress,
                   int ip);

  void stress_plastic(FEM &fem, std::vector<Node> &nodes,
                      std::vector<Eigen::Triplet<double>> &tripletsR);

  void elasto(const Eigen::VectorXd &strain, Eigen::VectorXd &stress, int ip);

  /// SIMP model for elasto-plastic model
  void SIMPEP();

  /// getters
  const std::vector<double> &Ehs() const { return Ehs_; }
  const std::vector<double> &sigma_ys() const { return sigma_ys_; }
  const double &Eh() const { return Eh_; }
  const double &sigma_y() const { return sigma_y_; }
  const double &betaH() const { return betaH_; }
  const double &tr_stress() const { return tr_stress_; }
  const double &epsilon() const { return epsilon_; }
  const double &eqstn() const { return eqstn_; }
  const Eigen::VectorXd &strain_gp(int ip) const { return strain_gp_[ip]; }
  const Eigen::VectorXd &stress_gp(int ip) const { return stress_gp_[ip]; }
  const double &xi_tri_norm(int ip) const { return xi_tri_norm_[ip]; }
  const double &strain_p_bar(int ip) const { return strain_p_bar_[ip]; }
  const double &gamma(int ip) const { return gamma_[ip]; }
  const double &yield(int ip) const { return yield_[ip]; }
  const double &sig3(int ip) const { return sig3_[ip]; }
  const double &stress_fin_bar(int ip) const { return stress_fin_bar_[ip]; }
  const Eigen::VectorXd &xi_tri(int ip) const { return xi_tri_[ip]; }

private:
  void plasto_simo_plane_strain(const Eigen::VectorXd &strain,
                                Eigen::VectorXd &stress, int ip);

  void plasto_simo_3d(const Eigen::VectorXd &strain, Eigen::VectorXd &stress,
                      int ip);

  void
  stress_plastic_plane_strain(FEM &fem, std::vector<Node> &nodes,
                              std::vector<Eigen::Triplet<double>> &tripletsR);

  void stress_plastic_3d(FEM &fem, std::vector<Node> &nodes,
                         std::vector<Eigen::Triplet<double>> &tripletsR);

protected:
  std::vector<double> Ehs_;      ///<
  std::vector<double> sigma_ys_; ///< Yield Stress
  double Eh_;                    ///< hardening parameter
  double sigma_y_;               ///< Yield Stress
  double betaH_; ///< parameter which combines kinematic and isotropic hardening

  double epsilon_;
  double tr_stress_;
  double eqstn_;

  /// stock values for analysis
  std::vector<Eigen::VectorXd> strain_p_; ///< plastic strain[ipmax][voigt]
  std::vector<Eigen::VectorXd> beta_;     ///<[ipmax][voigt]
  std::vector<double> strain_p_bar_;      ///<[ipmax]

  /// stock values for TO
  std::vector<Eigen::VectorXd> strain_gp_; ///< strain used for TO[ipmax][voigt]
  std::vector<Eigen::VectorXd> stress_gp_; ///< stress used for TO[ipmax][voigt]
  std::vector<double> xi_tri_norm_;
  std::vector<double> gamma_;
  std::vector<double> yield_;
  std::vector<double> sig3_;
  std::vector<double> stress_fin_bar_;
  std::vector<Eigen::VectorXd> xi_tri_;
};

class LinearElastPlast : public ElastPlast
{
public:
  LinearElastPlast(int voigt_, int ne_, int ipmax_, int numdof_)
      : ElastPlast(voigt_, ne_, ipmax_, numdof_)
  {
  }

  /// genarate element stiffness Ke, inner force Fe
  void makeKeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<Node> &nodes);
  int iteration = 0;
};

} // namespace icarat