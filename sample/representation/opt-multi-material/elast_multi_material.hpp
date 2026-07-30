#pragma once
#include <problem/elasticity.hpp>

namespace icarat
{
namespace opt_multi_material
{
/// elastic multi-material model class with using DMO interpolation scheme
class ElastMultiMat : public Elastic
{
public:
  ElastMultiMat(int voigt_, int ne_, int ipmax_, int numdof_, int numMat,
                std::string type)
      : Elastic(voigt_, ne_, ipmax_, numdof_), numMat_(numMat),
        design_ss_(numMat), eachYoung_(numMat), type_(type)
  {
    // Only the two "type" below are available.
    assert(type_ == "DMO" || type_ == "extendedSIMP");
  }

  /// set parameters to this element form config file
  void setParameter(const toml::value &config);

  /// make DMO young Modulus
  std::vector<double> makeCoeff(double pp);

  /// @return "matID1"'s derivative of DMO by material number "matID2"
  double dwds(int matID1, int matID2, double pp);
  double dYoung(int matID);

  /// setter
  void setDesignSS(int solidID, double d) { this->design_ss_[solidID] = d; }

  /// getters
  const std::string &type() const { return type_; }
  const int &numMat() const { return numMat_; }
  const double &design_ss(int matID) const { return design_ss_(matID); }
  const double &eachYoung(int matID) const { return eachYoung_(matID); }

protected:
  std::string type_; ///<"DMO" or "extendedSIMP"
  /// young1_...Young's modulus of void
  /// young2_...Young's modulus of solids using DMO function
  /// young_...Young's modulus of all materials & void
  int numMat_;                ///< number of solid material
  Eigen::VectorXd design_ss_; ///< design variable of each material
  Eigen::VectorXd eachYoung_; ///< Young modulus of solid material
};

/// elastic multi-material model class with using DMO interpolation scheme
class LinearElastMultiMat : public ElastMultiMat
{
public:
  LinearElastMultiMat(int voigt_, int ne_, int ipmax_, int numdof_, int numMat,
                      std::string type)
      : ElastMultiMat(voigt_, ne_, ipmax_, numdof_, numMat, type)
  {
  }

  /// make element stiffness matrix & internal force(Ke_ and fe_)
  void makeKeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<Node> &nodes);
};

} // namespace opt_multi_material
} // namespace icarat