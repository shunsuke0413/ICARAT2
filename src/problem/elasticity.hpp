///  @file  elasticity.hpp
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
/// SIMP for TO (is default model in ICARAT)
inline double SIMP(double young1, double young2, double design_s, double pp)
{
  assert(young2 - young1 >= 0.0);
  return (1.0 - pow(design_s, pp)) * young1 + pow(design_s, pp) * young2;
}
/// RAMP for TO
inline double RAMP(double young1, double young2, double design_s, double pp)
{
  assert(young2 - young1 >= 0.0);
  return young1 +
         (design_s / (1.0 + pp * (1.0 - design_s))) * (young2 - young1);
}
///  isotropic elastic material model class
class Elastic : public Element
{
public:
  /// @param voigt_ ... 3 in 2D, 6 in 3D problems
  /// @param ne_ ... node in this element
  /// @param ipmax_ ... number of gaussian integration points
  /// @param numdof_ ... number of DOFs in this element
  Elastic(int voigt_, int ne_, int ipmax_, int numdof_)
      : De_(voigt_, voigt_), De_dU_(voigt_, voigt_), Ke_(numdof_, numdof_),
        fe_(numdof_), strain_(voigt_), stress_(voigt_)
  {
    ne = ne_;
    ipmax = ipmax_;
    numdof = numdof_;
  }
  /// set De parameters to this element form input file
  void setParameter(const toml::value &config);
  /// isotropic elastic modulus tensor (De_)
  void De(Eigen::MatrixXd &De, double E, double nu);
  /// isotropic elastic modulus tensor (De_) with using Lame's constants
  void lameDe(Eigen::MatrixXd &De, double mu, double lambda);
  /// make mises stress in this elmenet (mises_)
  /// The stresses calculated for each Gaussian point are averaged.
  void makeMisesStress(FEM &fem, std::vector<Node> &node);
  /// setters
  void setDesignS(double &d) { this->design_s_ = d; }
  /// getters
  const Eigen::MatrixXd &De() const { return De_; }
  const Eigen::MatrixXd &Ke() const { return Ke_; }
  const Eigen::MatrixXd &nKe() const { return nKel_; }
  const Eigen::VectorXd &fe() const { return fe_; }
  const std::string &mattype() const { return mattype_; }
  const double &design_s() const { return design_s_; }
  const double &young() const { return young_; }
  const double &young1() const { return young1_; }
  const double &young2() const { return young2_; }
  const double &poisson() const { return poisson_; }
  const double &pp() const { return pp_; }
  const double &mu() const { return mu_; }
  const double &lambda() const { return lambda_; }
  const double &trac_inv() const { return traC_inv; }
  const double &trac() const { return traC_; }
  const Eigen::VectorXd &strain() const { return strain_; }
  const Eigen::VectorXd &stress() const { return stress_; }
  const double &mises() const { return mises_; }

protected:
  Eigen::MatrixXd De_; ///< elastic modulus tensor
  Eigen::MatrixXd De_dU_;
  Eigen::MatrixXd dNdX_;
  Eigen::MatrixXd Ke_; ///< stiffness matrix
  Eigen::MatrixXd cinv_;
  Eigen::MatrixXd cres_;
  Eigen::MatrixXd est_;
  Eigen::VectorXd fe_; ///< inner force
  Eigen::VectorXd dufe_;
  Eigen::VectorXd cresv_;
  Eigen::MatrixXd nKel_;
  std::string mattype_; ///< 1:PLANE_STRESS 2:PLANE_STRAIN 3:3D
  double design_s_;     ///< design variable
  double young_;        ///< young modulus of SIMP
  double young1_;       ///< young modulus of material1
  double young2_;       ///< young modulus of material2
  double poisson_;      ///< poisson rate
  double pp_;           ///< penalty parameter for SIMP method
  double mu_;
  double lambda_;
  double laJ_;
  double traC_inv;
  double traC_;
  Eigen::VectorXd strain_; ///< element-mean strain
  Eigen::VectorXd stress_; ///< element-mean stress
  double mises_;           ///< element-mean Mises stress
};
/// elastic material model class used in linear analysis
class LinearElastic : public Elastic
{
public:
  LinearElastic(int voigt_, int ne_, int ipmax_, int numdof_)
      : Elastic(voigt_, ne_, ipmax_, numdof_)
  {
  }
  /// make element stiffness matrix & internal force(Ke_ and fe_)
  void makeKeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<Node> &nodes);
};
/// material model for modal analysis class
class ModalLinearElastic : public LinearElastic
{
public:
  ModalLinearElastic(int voigt_, int ne_, int ipmax_, int numdof_,
                     double density)
      : LinearElastic(voigt_, ne_, ipmax_, numdof_), Me_(numdof_, numdof_),
        density_(density)
  {
  }
  /// make geometric stiffness (KGe_)
  void makeEigenCoeffs(std::vector<Eigen::Triplet<double>> &tripletsA,
                       std::vector<Eigen::Triplet<double>> &tripletsB, FEM &fem,
                       std::vector<Node> &nodes);
  /// getters
  const Eigen::MatrixXd &Me() const { return Me_; }
  const double &density() const { return density_; }

protected:
  Eigen::MatrixXd Me_; ///< elastic modulus tensor
  double density_;     ///< material density (not design variable!)
};
/// linear buckling material model class
class BucklingLinearElastic : public LinearElastic
{
public:
  BucklingLinearElastic(int voigt_, int ne_, int ipmax_, int numdof_)
      : LinearElastic(voigt_, ne_, ipmax_, numdof_), KGe_(numdof_, numdof_)
  {
  }
  /// make geometric stiffness (KGe_)
  void makeEigenCoeffs(std::vector<Eigen::Triplet<double>> &tripletsA,
                       std::vector<Eigen::Triplet<double>> &tripletsB, FEM &fem,
                       std::vector<Node> &nodes);
  /// Non-design area flag
  bool isNonDesign = false;
  // getters
  const Eigen::MatrixXd &KGe() const { return KGe_; }

protected:
  Eigen::MatrixXd KGe_; ///< element geometric stiffness matrix
};
class DynamicElastic : public Elastic
{
public:
  DynamicElastic(int voigt_, int ne_, int ipmax_, int numdof_, double density)
      : Elastic(voigt_, ne_, ipmax_, numdof_), Me_(numdof_, numdof_),
        density_(density)
  {
  }
  /// make element stiffness "Ke_", mass matrix "Me_" and inner force "fe_"
  /// @param  ExpOrImp = "implicit" or "explicit"
  /// implicit...calc the consistent mass matrix for "tripletsM"
  /// explicit...calc inverse of the concentrated mass matrix for "tripletsM"
  void makeKeMeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                    std::vector<Eigen::Triplet<double>> &tripletsM,
                    std::vector<Eigen::Triplet<double>> &tripletsF,
                    std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                    std::vector<Node> &nodes, std::string ExpOrImp);
  /// getters
  const Eigen::MatrixXd &Me() const { return Me_; }
  const double &density() const { return density_; }

protected:
  Eigen::MatrixXd Me_; ///< elastic modulus tensor
  double density_;     ///< material density (not design variable!)
};
/// node variable list
struct HomoNode : public Node
{
  HomoNode() : Node(1) {}
  /// solved values obtained by NMT
  double nmtval[6][3];      ///< size: [voigt][dofnp]
  std::shared_ptr<MPC> mpc; ///< MPC structure
};
/// micro scale elastic model class based on homogenization theory
///@sa  寺田賢二郎ら，「数値材料試験」，丸善出版，2021.
class mHomoElastic : public Elastic
{
public:
  mHomoElastic(int voigt_, int ne_, int ipmax_, int numdof_)
      : Elastic(voigt_, ne_, ipmax_, numdof_), sEnergy_(voigt_, voigt_)
  {
  }
  /// used for homogenization in micro analysis
  void makeKeHomo(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<HomoNode> &nodes, int dir);
  /// generate the homogenized elasticity tensor in micro analysis
  void makeMicroValues(FEM &fem, std::vector<HomoNode> &nodes);
  /// getter
  const Eigen::MatrixXd &sEnergy() const { return sEnergy_; }

protected:
  Eigen::MatrixXd sEnergy_; ///< strain energy
};
/// macro scale elastic model class used in multi-scale
class MHomoElastic : public Elastic
{
public:
  MHomoElastic(int voigt_, int ne_, int ipmax_, int numdof_)
      : Elastic(voigt_, ne_, ipmax_, numdof_)
  {
  }
  /// used in macro analysis
  void makeKeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<Node> &nodes); ///< element stiffness
  /// setters
  /// De used for Macro analysis
  void setDe(Eigen::MatrixXd &De) { this->De_ = De; }

protected:
};
/// Thermo-structural coupled material model considering thermal stress
class ThermalStressElastic : public Elastic
{
public:
  ThermalStressElastic(int voigt_, int ne_, int ipmax_, int numdof_,
                       std::vector<double> alpha)
      : Elastic(voigt_, ne_, ipmax_, numdof_), alpha_(alpha), fthi_(numdof_)
  {
  }
  /// make thermal stress coefficient beta
  double thermalStress(const double &designS, const double &young1,
                       const double &young2, const double &alpha1,
                       const double &alpha2);
  /// generate thermal load
  Eigen::VectorXd generateThermalLoad(FEM &fem, std::vector<Node> &node,
                                      double deltaT);
  /// make element stiffness matrix & internal force(Ke_ and fe_)
  void makeKeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<Node> &nodes);
  /// getters
  const double &alpha(int i) const { return alpha_[i]; }
  const double &beta() const { return beta_; }
  const Eigen::VectorXd &fthi() const { return fthi_; }

private:
  std::vector<double> alpha_; ///< thermal expansion of materials
  double beta_;               ///< thermal stress coefficient
  Eigen::VectorXd fthi_;      /// thermal load
};
} // namespace icarat