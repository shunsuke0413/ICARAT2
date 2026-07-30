///  @file  m _plasto_.hpp
///  @author  Daiki Watanabe
///  @date  October 3, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include "elasticity.hpp"

namespace icarat
{
/// Elasto-plastic material model class based on multiplicative decomposition
/// This class was implemented by referring to the following book.
/// @cite N. H. Kim, introduction to Nonlinear Finite Element Analysis.
class MultipliEPUpdated : public Elastic
{
public:
  MultipliEPUpdated(int voigt_, int ne_, int ipmax_, int numdof_)
      : Elastic(voigt_, ne_, ipmax_, numdof_), mus_(2), lambdas_(2), Hs_(2),
        Y0s_(2), pstress_(ipmax_), alpha_(ipmax_), ep_(ipmax_)
  {
    for(int i = 0; i < ipmax_; i++)
    {
      pstress_[i].resize(6);
      pstress_[i].setZero();

      alpha_[i].resize(3);
      alpha_[i].setZero();
    }
  }

  /// set material parameter from .toml file
  void setParameter(const toml::value &config);

  void initialization();

  void updateVariables(FEM &fem, std::vector<Node> &nodes);

  /// return mapping algorithm
  ///@param [in] DD = elasticity matrix b/w prin stress & log prin stretch (3x3)
  ///@param [in] LL = [dui / dxj] velocity gradient
  ///@param [in,out] b = elastic left C-G deformation vector(6x1)
  /// @param [in,out] alpha = principal back stress(3x1)
  ///@param [in,out] ep = effective plastic strain
  ///@param [in,out] flag = return mapping(true) or tangent stiffness De_(false)
  void returnMapping(FEM &fem, Eigen::VectorXd &sigma,
                     Eigen::MatrixXd &sigmaMatrix, Eigen::MatrixXd &DD,
                     Eigen::MatrixXd &LL, Eigen::VectorXd &bb,
                     Eigen::VectorXd &alpha, double &ep, bool flag);

  /// SIMP model for elasto-plastic model
  void SIMPEP();

  /// getters
  const std::vector<double> &ep() const { return ep_; }

protected:
  std::vector<double> mus_;     ///< lame constant
  std::vector<double> lambdas_; ///< lame constant
  std::vector<double> Hs_;      ///<
  std::vector<double> Y0s_;     ///< Yield Stress

  double mu_;     ///< lame constant
  double lambda_; ///< lame constant
  double beta_;   ///<
  double H_;      ///<
  double Y0_;     ///< Yield Stress

  /// stock values
  std::vector<Eigen::VectorXd> pstress_; ///< previous stress
  std::vector<Eigen::VectorXd> alpha_;   ///<
  std::vector<double> ep_;               ///<
};

class FiniteMultipliEPUpdated : public MultipliEPUpdated
{
public:
  FiniteMultipliEPUpdated(int voigt_, int ne_, int ipmax_, int numdof_)
      : MultipliEPUpdated(voigt_, ne_, ipmax_, numdof_)
  {
  }

  /// genarate element stiffness Ke, inner force Fe
  void makeKeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<Node> &nodes);
};
} // namespace icarat