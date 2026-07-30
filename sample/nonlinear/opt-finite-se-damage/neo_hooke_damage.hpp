///  @file  neo_hooke.hpp
///  @author  Daiki Watanabe
///  @date  June 4, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <problem/elasticity.hpp>
#include <unsupported/Eigen/CXX11/Tensor>

namespace icarat
{
/// compressible Neo-Hookean model of total lagrange
class neohookedamage : public Elastic
{
public:
  neohookedamage(int voigt_, int ne_, int ipmax_, int numdof_)
      : Elastic(voigt_, ne_, ipmax_, numdof_)
  {
    D_all_.assign(ipmax_, 0.0);
  }

  //----------------------- for S:E integration ------------------//
  void totalKirchhoff(FEM &fem, Eigen::VectorXd &sigmaVector,
                      Eigen::MatrixXd &F, Eigen::MatrixXd &sigmaMatrix,
                      double &product, int ip);
  //--------------------------------------------------------------//

  //----------------------- for damage model-----------------------//
  double calculateTensileDamage(double tau);

  double differentiateD(double tau);

  // 損傷変数の管理用
  std::vector<double> D_all_;

  // 今現在のガウス点インデックスを保持する変数（計算中のポイントを追跡）
  int current_ip_ = 0;

  double getAverageDamage_() const
  {
    if(D_all_.empty())
      return 0.0;
    double sum = 0.0;
    for(const double &d : D_all_)
      sum += d;
    return sum / static_cast<double>(D_all_.size());
  }
  //---------------------------------------------------------------//
};

/// stabilizer class for numerical instability in the finite deformation's TO
///  @sa Wang G. et al.: Interpolation scheme for fictitious domain techniques
///  and topology optimization of finite strain elastic problems (2014).
struct LinearStabilizer
{
  LinearStabilizer(int voigt_, int numdof_)
      : strain(voigt_), stress(voigt_), fe(numdof_), dufe(numdof_),
        De(voigt_, voigt_), Ke(numdof_, numdof_)
  {
  }

  void init()
  {
    strain.setZero();
    stress.setZero();
    fe.setZero();
    dufe.setZero();
    De.setZero(); // linear elastic matrix
    Ke.setZero(); // linear stiffness matrix
  }

  /// linear parts
  Eigen::VectorXd strain;
  Eigen::VectorXd stress;
  Eigen::VectorXd fe;
  Eigen::VectorXd dufe;
  Eigen::MatrixXd De;
  Eigen::MatrixXd Ke;
};

/// compressible Neo-Hookean model for finite deformation
class FiniteNeoHookedamage : public neohookedamage
{
public:
  FiniteNeoHookedamage(int voigt_, int ne_, int ipmax_, int numdof_,
                       int totalTime)
      : NeoHookedamage(voigt_, ne_, ipmax_, numdof_), gamma_(0.0),
        l_(voigt_, numdof_), lfe_(numdof_), Kegamma_(numdof_, numdof_)
  {
  }

  // --- この中継関数を追加してください ---
  double getAverageDamage() const
  {
    return NeoHookedamage::getAverageDamage_();
  }

  void setGamma(double gamma) { gamma_ = gamma; }
  //----------------------- for S:E integration ------------------//
  void makeKeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<Node> &nodes, int thisTime, double &product);
  //--------------------------------------------------------------//

  Eigen::VectorXd product_tempv;
  const double &gamma() const { return gamma_; }
  const Eigen::VectorXd &lfe() const { return lfe_; }
  const Eigen::VectorXd &nfe() const { return nfe_; }
  const Eigen::MatrixXd &kegamma() const { return Kegamma_; }

  double getAverageDamage2() const
  {
    if(D_all_.empty())
      return 0.0;
    double sum = 0.0;
    for(const double &d : D_all_)
      sum += d;
    return sum / static_cast<double>(D_all_.size());
  }

protected:
  double gamma_;            ///< threshold parameter
  LinearStabilizer l_;      ///< fictitious linear class
  Eigen::VectorXd lfe_;     ///< linear inner force
  Eigen::VectorXd nfe_;     ///< nonlinear inner force
  Eigen::MatrixXd Kegamma_; ///< nonlinear nonlinear tangential stiffness
  Eigen::MatrixXd lKe_;
  Eigen::MatrixXd nKe_;
};

class DynamicNeoHookedamage : public NeoHookedamage
{
public:
  DynamicNeoHookedamage(int voigt_, int ne_, int ipmax_, int numdof_,
                        double density)
      : NeoHookedamage(voigt_, ne_, ipmax_, numdof_), gamma_(0.0),
        l_(voigt_, numdof_), lfe_(numdof_), nfe_(numdof_),
        Me_(numdof_, numdof_), density_(density)
  {
  }

  void setGamma(double gamma) { gamma_ = gamma; }

  /// make element stiffness "Ke_", mass matrix "Me_" and inner force "fe_"
  /// @param  ExpOrImp = "implicit" or "explicit"
  /// implicit...calc the consistent mass matrix for "tripletsM"
  /// explicit...calc inverse of the concentrated mass matrix for "tripletsM"
  void makeKeMeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                    std::vector<Eigen::Triplet<double>> &tripletsM,
                    std::vector<Eigen::Triplet<double>> &tripletsF,
                    std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                    std::vector<Node> &nodes, std::string ExpOrImp,
                    double &product);

  /// getters
  const Eigen::MatrixXd &Me() const { return Me_; }
  const double &density() const { return density_; }
  const double &gamma() const { return gamma_; }
  const Eigen::VectorXd &lfe() const { return lfe_; }
  const Eigen::VectorXd &nfe() const { return nfe_; }
  const Eigen::MatrixXd &nke() const { return nke_; }

protected:
  double gamma_;        ///< threshold parameter
  LinearStabilizer l_;  ///< fictitious linear class
  Eigen::VectorXd lfe_; ///< linear inner force
  Eigen::VectorXd nfe_; ///< nonlinear inner force
  Eigen::MatrixXd nke_; ///< nonlinear stiffness matrix
  Eigen::MatrixXd Me_;  ///< elastic modulus tensor
  double density_;      ///< material density (not design variable!)
};

/// micro scale elastic model class based on homogenization theory
///@sa  寺田賢二郎ら，「数値材料試験」，丸善出版，2021.
class mHomoNeoHookedamage : public NeoHookedamage
{
public:
  mHomoNeoHookedamage(int voigt_, int ne_, int ipmax_, int numdof_)
      : NeoHookedamage(voigt_, ne_, ipmax_, numdof_), gamma_(0.0),
        l_(voigt_, numdof_), sEnergy_(voigt_, voigt_)
  {
  }

  void setGamma(double gamma) { gamma_ = gamma; }

  /// used for homogenization in micro analysis
  void makeKeHomo(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<HomoNode> &nodes, int dir);

  /// generate the homogenized elasticity tensor in micro analysis
  void makeMicroValues(FEM &fem, std::vector<HomoNode> &nodes);

  /// getter
  const double &gamma() const { return gamma_; }
  const Eigen::MatrixXd &sEnergy() const { return sEnergy_; }

protected:
  double gamma_;            ///< threshold parameter
  LinearStabilizer l_;      ///< fictitious linear class
  Eigen::MatrixXd sEnergy_; ///< strain energy
};

} // namespace icarat