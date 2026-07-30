///  @file  mooney.hpp
///  @author  Daiki Watanabe
///  @date   October 3, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include "saint_venant.hpp"

namespace icarat
{
/// This class is Mooney-Rivlin model
class MooneyTotal : public Elastic
{
public:
  MooneyTotal(int voigt_, int ne_, int ipmax_, int numdof_)
      : Elastic(voigt_, ne_, ipmax_, numdof_),
        C10s_(2), /// 2-phase material in SIMP
        C01s_(2), D1s_(2)
  {
  }

  /// set material parameter from .dat file
  void setParameter(const toml::value &config);

  /// compute 2nd PK stress (sigmaVector, sigmaMatrix) and Material
  /// stiffness(De_)
  void totalKirchhoff(FEM &fem, Eigen::VectorXd &sigmaVector,
                      Eigen::MatrixXd &F, Eigen::MatrixXd &sigmaMatrix);

  /// SIMP model in Mooney
  void SIMPMooney();

  /// getters
  const double &C10s(int i) const { return C10s_[i]; }
  const double &C01s(int i) const { return C01s_[i]; }
  const double &D1s_inv(int i) const { return D1s_[i]; }
  const double &C10() const { return C10_; }
  const double &C01() const { return C01_; }
  const double &D1_inv() const { return D1_; }

protected:
  std::vector<double> C10s_; ///< each material parameter
  std::vector<double> C01s_; ///< each material parameter
  std::vector<double> D1s_;  ///< each inverse of incompressibility parameter
  double C10_;               ///< material parameter
  double C01_;               ///< material parameter
  double D1_;                ///< incompressibility parameter
};

/// Mooney-Rivlin model for finite deformation
class FiniteMooneyTotal : public MooneyTotal
{
public:
  FiniteMooneyTotal(int voigt_, int ne_, int ipmax_, int numdof_)
      : MooneyTotal(voigt_, ne_, ipmax_, numdof_)
  {
  }

  void makeKeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<Node> &nodes);
};
} // namespace icarat