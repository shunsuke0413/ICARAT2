///  @file  elast_heat.hpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include "base.hpp"
#include <toml/toml.hpp>
namespace icarat
{
/// heat material model
class Heat : public Element
{
public:
  /// voigt...size of voigt expression. 2 in 2D, 3 in 3D problems
  /// numdof...number of DOFs in this element
  Heat(int voigt, int ne_, int ipmax_, int numdof_)
      : De_(voigt, voigt), De1_(voigt, voigt), De2_(voigt, voigt),
        Ke_(numdof_, numdof_), fe_(numdof_), strain_(voigt), stress_(voigt)
  {
    ne = ne_;
    ipmax = ipmax_;
    numdof = numdof_;
  }

  /// set De parameters to this element form input file
  void setParameter(const toml::value &config);

  /// make elastic modulus tensor (De_) with SIMP model
  void makeSIMPDe();

  /// setters
  void setDesign(double d) { this->design_s_ = d; }

  /// getters
  const Eigen::MatrixXd &De() const { return De_; }
  const Eigen::MatrixXd &De1() const { return De1_; }
  const Eigen::MatrixXd &De2() const { return De2_; }
  const Eigen::MatrixXd &Ke() const { return Ke_; }
  const Eigen::VectorXd &fe() const { return fe_; }
  const double &design_s() const { return design_s_; }
  const double &pp() const { return pp_; }
  const Eigen::VectorXd &strain() const { return strain_; }
  const Eigen::VectorXd &stress() const { return stress_; }
  const double &conductivity(int i) const { return specificHeat_[i]; }

protected:
  Eigen::MatrixXd De_;               ///< composite flux matrix
  Eigen::MatrixXd De1_;              ///< Flux matrix of material1
  Eigen::MatrixXd De2_;              ///< Flux matrix of material1
  Eigen::MatrixXd Ke_;               ///< stiffness matrix
  Eigen::VectorXd fe_;               ///< inner force
  double design_s_;                  ///< design_s
  double pp_;                        ///< penalty parameter for SIMP method
  Eigen::VectorXd strain_;           ///< strain
  Eigen::VectorXd stress_;           ///< stress
  std::vector<double> specificHeat_; ///< specific heat
};

/// heat matrial model for linear & static problem
class LinearHeat : public Heat
{
public:
  LinearHeat(int voigt_, int ne_, int ipmax_, int numdof_)
      : Heat(voigt_, ne_, ipmax_, numdof_)
  {
  }

  /// make element stiffness matrix & internal force(Ke_ and fe_)
  void makeKeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<Node> &nodes);
};
} // namespace icarat