#include "thermal_stress.hpp"
#include <basis/bmatrix.hpp>
#include <misc/assembling.hpp>

namespace icarat
{
using namespace std;
using namespace Eigen;

namespace
{
inline double RAMPcoeff(double &s, double &p)
{
  return s / (1.0 + p * (1.0 - s));
}
inline double dRAMPcoeff(double &s, double &p)
{
  return (1.0 + p) / pow(1.0 + p * (1.0 - s), 2.0);
}
} // namespace

void ThermalStressElastMultiMat::setParameter(const toml::value &config)
{
  assert(toml::find<string>(config, "material", "type") ==
         "ThermalStressElastMultiMat");

  poisson_ = toml::find<double>(config, "material", "poisson");
  pp_ = toml::find<double>(config, "material", "penalty");
  young1_ = toml::find<double>(config, "material", "young_void");

  auto designs_tmp = toml::find<vector<double>>(config, "material", "designs");
  auto eaYoung_tmp = toml::find<vector<double>>(config, "material", "young");
  auto alpha_tmp = toml::find<vector<double>>(config, "material", "alpha");

  for(int matID = 0; matID < numMat_; matID++)
  {
    design_s_(matID) = designs_tmp[matID];
    eachYoung_(matID) = eaYoung_tmp[matID];
    alpha_(matID) = alpha_tmp[matID];
  }

  mattype_ = toml::find<string>(config, "material", "state");
  if(toml::find<int>(config, "mesh", "dimension") == 3)
    mattype_ = "3D";
}

VectorXd ThermalStressElastMultiMat::generateThermalLoad(FEM &fem,
                                                         vector<Node> &node,
                                                         double deltaT)
{
  if(isDesignable)
  {
    beta_ = 0.0;
    auto ww = makeCoeff(0.0);
    for(int i = 0; i < numMat_; i++)
      beta_ += ww[i] * alpha_[i] * eachYoung_[i];
  }
  else
  {
    poisson_ = 0.2;        // glass
    beta_ = 1.0e-4 * 72.5; // glass
  }

  VectorXd phi(fem.voigt);
  if(fem.ndim == 2)
    phi << 1.0, 1.0, 0.0;
  else
    phi << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;

  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
  MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  VectorXi idof = VectorXi::Zero(numdof);

  for(int j = 0; j < ne; j++)
  {
    for(int k = 0; k < fem.ndim; k++)
    {
      idof.coeffRef(fem.dofnp * j + k) = node[nodeID[j]].dof[k];
      X.coeffRef(j, k) = node[nodeID[j]].x[k];
    }
  }

  VectorXd fthi = VectorXd::Zero(numdof);
  De(De_, 1.0, poisson_);
  Bmatrix bmatrix;
  for(int ip = 0; ip < ipmax; ip++)
  {
    double jac;
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);
    fthi += Be.transpose() * De_ * phi * jac;
  }

  fthi *= deltaT * beta_;
  fthi_ = fthi;

  // assembling
  VectorXd fth = VectorXd::Zero(fem.numeq);
  for(int i = 0; i < numdof; i++)
  {
    if(idof[i] < fem.numeq)
      fth[idof[i]] = fthi[i];
    else
      continue;
  }

  return fth;
}

vector<double> ThermalStressElastMultiMat::makeCoeff(double pp)
{
  vector<double> dens;
  // solid material
  for(int i = 0; i < numMat_; i++)
  {
    // double tmpDensity = pow(design_s_[i], pp);
    double tmpDensity = RAMPcoeff(design_s_[i], pp);
    for(int j = 0; j < numMat_; j++)
    {
      if(i == j)
        continue;
      else
        tmpDensity *= (1.0 - RAMPcoeff(design_s_[j], pp));
    }
    dens.push_back(tmpDensity);
  }
  return dens;
}

double ThermalStressElastMultiMat::dwds(int matID1, int matID2, double pp)
{
  vector<double> dwds;
  for(int i = 0; i < numMat_; i++)
  {
    double tmp = 0.0;
    if(matID2 == i)
    {
      // tmp = pp * pow(design_s_[matID2], pp - 1.0);
      tmp = dRAMPcoeff(design_s_[matID2], pp);
      for(int j = 0; j < numMat_; j++)
      {
        if(i == j)
          continue;
        else
          // tmp *= (1.0 - pow(design_s_[j], pp));
          tmp *= (1.0 - RAMPcoeff(design_s_[j], pp));
      }
    }
    else
    {
      // tmp = -pp * pow(design_s_[matID2], pp - 1.0) * pow(design_s_[i], pp);
      tmp = -dRAMPcoeff(design_s_[matID2], pp) * RAMPcoeff(design_s_[i], pp);
      for(int j = 0; j < numMat_; j++)
      {
        if(i == j)
          continue;
        else if(matID2 == j)
          continue;
        else
          tmp *= (1.0 - RAMPcoeff(design_s_[j], pp));
        // tmp *= (1.0 - pow(design_s_[j], pp));
      }
    }
    dwds.push_back(tmp);
  }

  return dwds[matID1];
}

void LinearThermalStressElastMultiMat::makeKeFint(
    std::vector<Eigen::Triplet<double>> &tripletsK,
    std::vector<Eigen::Triplet<double>> &tripletsF,
    std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
    std::vector<Node> &nodes)
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

  if(isDesignable)
  {
    // get multi-material Young's modulus
    auto coeff = makeCoeff(pp_);
    young2_ = 0.0;
    for(int i = 0; i < numMat_; i++)
      young2_ += coeff[i] * eachYoung_(i);

    young_ = young1_ + young2_;
  }
  else
  {
    poisson_ = 0.2; // glass
    young_ = 72.5;  // glass
  }

  De(De_, young_, poisson_);

  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;
  double jac = 0.0;
  Bmatrix bmatrix;

  // gauss loop
  for(int ip = 0; ip < ipmax; ip++)
  {
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);

    fe_ += Be.transpose() * De_ * Be * disp * jac;
    Ke_ += Be.transpose() * De_ * Be * jac;
    volume += jac;
  }

  assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
  assemblingReact(numdof, idof, fe_, tripletsR);
}

} // namespace icarat