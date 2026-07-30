#pragma once
#include <problem/elasticity.hpp>

namespace icarat
{

/// elastic multi-material model class with using DMO interpolation scheme
///@sa Gao T. et al. Topology optimization of thermo-elastic structures with
/// multiple materials under mass constraint, 2016.
class ThermalStressElastMultiMat : public Elastic
{
public:
  ThermalStressElastMultiMat(int voigt_, int ne_, int ipmax_, int numdof_,
                             int numMat)
      : Elastic(voigt_, ne_, ipmax_, numdof_), numMat_(numMat),
        design_s_(numMat), eachYoung_(numMat), alpha_(numMat), fthi_(numdof_)
  {
  }

  /// set parameters to this element form config file
  void setParameter(const toml::value &config);

  Eigen::VectorXd generateThermalLoad(FEM &fem, std::vector<Node> &node,
                                      double deltaT);

  /// make DMO young Modulus
  std::vector<double> makeCoeff(double pp);

  /// @return "matID1"'s derivative of DMO by material number "matID2"
  double dwds(int matID1, int matID2, double pp);

  /// setter
  void setDesignS(int solidID, double d) { this->design_s_[solidID] = d; }

  /// getters
  const int &numMat() const { return numMat_; }
  const double &design_s(int matID) const { return design_s_(matID); }
  const double &eachYoung(int matID) const { return eachYoung_(matID); }
  const double &alpha(int matID) const { return alpha_(matID); }
  const double &beta() const { return beta_; }
  const Eigen::VectorXd &fthi() const { return fthi_; }

protected:
  /// young1_...Young's modulus of void
  /// young2_...Young's modulus of solids using DMO function
  /// young_...Young's modulus of all materials & void
  int numMat_;                ///< number of solid material
  Eigen::VectorXd design_s_;  ///< design variable of each material
  Eigen::VectorXd eachYoung_; ///< Young modulus of solid material
  Eigen::VectorXd alpha_;     ///< thermal expansion of solid material
  double beta_;               ///< thermal stress coefficient
  Eigen::VectorXd fthi_;      /// thermal load
};

/// elastic multi-material model class with using DMO interpolation scheme
class LinearThermalStressElastMultiMat : public ThermalStressElastMultiMat
{
public:
  LinearThermalStressElastMultiMat(int voigt_, int ne_, int ipmax_, int numdof_,
                                   int numMat)
      : ThermalStressElastMultiMat(voigt_, ne_, ipmax_, numdof_, numMat)
  {
  }

  /// make element stiffness matrix & internal force(Ke_ and fe_)
  void makeKeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<Node> &nodes);
};

} // namespace icarat