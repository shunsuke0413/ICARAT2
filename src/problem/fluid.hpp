///  @file  fluid.hpp
///  @author  Daiki Watanabe
///  @date  October 28, 2023.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include "base.hpp"
#include <toml/toml.hpp>
namespace icarat
{
/// heat material model
class Fluid : public Element
{
public:
  /// voigt...size of voigt expression. 2 in 2D, 3 in 3D problems
  /// numdof...number of DOFs in this element
  Fluid(int voigt, int ne_, int ipmax_, int numdof_)
      : De_(voigt, voigt), Ke_(numdof_, numdof_), fe_(numdof_), strain_(voigt),
        stress_(voigt)
  {
    ne = ne_;
    ipmax = ipmax_;
    numdof = numdof_;
  }

  /// set De parameters to this element form input file
  void setParameter(const toml::value &config);

  /// setters
  void setDesign(double d) { this->design_s_ = d; }

  /// getters
  const Eigen::MatrixXd &De() const { return De_; }
  const Eigen::MatrixXd &Ke() const { return Ke_; }
  const Eigen::VectorXd &fe() const { return fe_; }
  const double &design_s() const { return design_s_; }
  const double &pp() const { return pp_; }
  const Eigen::VectorXd &strain() const { return strain_; }
  const Eigen::VectorXd &stress() const { return stress_; }
  const double &eta() const { return eta_; }
  const double &rho() const { return rho_; }

protected:
  Eigen::MatrixXd De_;     ///< kinematic viscosity matrix
  Eigen::MatrixXd Ke_;     ///< element stiffness matrix
  Eigen::VectorXd fe_;     ///< inner force
  double design_s_;        ///< design_s
  double pp_;              ///< penalty parameter for SIMP method
  Eigen::VectorXd strain_; ///< strain
  Eigen::VectorXd stress_; ///< stress
  double eta_;             ///< kinematic viscosity coefficient
  double rho_;             ///< density
  std::string eTypep_;     ///< eType for pressure
  int nep_;                ///< ne for pressure
};

/// compute element stiffness & inner force of Stokes flow
class Stokes : public Fluid
{
public:
  Stokes(int voigt_, int ne_, int ipmax_, int numdof_)
      : Fluid(voigt_, ne_, ipmax_, numdof_)
  {
  }

  /// make element stiffness matrix & internal force(Ke_ and fe_)
  void makeKeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<Node> &nodes);
};
} // namespace icarat