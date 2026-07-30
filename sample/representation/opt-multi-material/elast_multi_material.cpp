#include "elast_multi_material.hpp"
#include <basis/bmatrix.hpp>
#include <misc/assembling.hpp>

namespace icarat
{
namespace opt_multi_material
{
using namespace std;
using namespace Eigen;

void ElastMultiMat::setParameter(const toml::value &config)
{
  assert(toml::find<string>(config, "material", "type") == "ElastMultiMat");

  poisson_ = toml::find<double>(config, "material", "poisson");
  pp_ = toml::find<double>(config, "material", "penalty");
  young1_ = toml::find<double>(config, "material", "young_void");

  auto designs_tmp = toml::find<vector<double>>(config, "material", "designs");
  auto eaYoung_tmp = toml::find<vector<double>>(config, "material", "young");

  for(int matID = 0; matID < numMat_; matID++)
  {
    design_ss_(matID) = designs_tmp[matID];
    eachYoung_(matID) = eaYoung_tmp[matID];
  }

  // cout << "nummat: " << numMat_ << endl;
  // cout << "design_s: " << design_ss_ << endl;
  // cout << "eachY: " << eachYoung_ << endl;

  mattype_ = toml::find<string>(config, "material", "state");
  if(toml::find<int>(config, "mesh", "dimension") == 3)
    mattype_ = "3D";
}

vector<double> ElastMultiMat::makeCoeff(double pp)
{
  vector<double> dens;
  // DMO
  if(type_ == "DMO")
  {
    // void material
    dens.push_back(1.0);

    // solid material
    for(int i = 0; i < numMat_; i++)
    {
      double tmpDensity = pow(design_ss_[i], pp);
      for(int j = 0; j < numMat_; j++)
      {
        if(i == j)
          continue;
        else
          tmpDensity *= (1.0 - pow(design_ss_[j], pp));
      }
      dens.push_back(tmpDensity);
    }
  }
  // extended SIMP
  else
  {
    // i=0 ... void
    for(int i = 0; i < numMat_ + 1; i++)
    {
      double tmpDensity;
      if(i == numMat_)
        tmpDensity = 1.0;
      else
        tmpDensity = 1.0 - pow(design_ss_[i], pp);

      if(i == 0)
        tmpDensity *= 1.0;
      else
        for(int j = 0; j < i; j++)
          tmpDensity *= pow(design_ss_[j], pp);

      dens.push_back(tmpDensity);
    }
  }
  return dens;
}

double ElastMultiMat::dwds(int matID1, int matID2, double pp)
{
  vector<double> dwds;
  // DMO
  if(type_ == "DMO")
  {
    // void material
    dwds.push_back(0.0);

    for(int i = 0; i < numMat_; i++)
    {
      double tmp = 0.0;
      if(matID2 == i)
      {
        tmp = pp * pow(design_ss_[matID2], pp - 1.0);
        for(int j = 0; j < numMat_; j++)
        {
          if(i == j)
            continue;
          else
            tmp *= (1.0 - pow(design_ss_[j], pp));
        }
      }
      else
      {
        tmp = -pp * pow(design_ss_[matID2], pp - 1.0) * pow(design_ss_[i], pp);
        for(int j = 0; j < numMat_; j++)
        {
          if(i == j)
            continue;
          else if(matID2 == j)
            continue;
          else
            tmp *= (1.0 - pow(design_ss_[j], pp));
        }
      }
      dwds.push_back(tmp);
    }
  }
  // extended SIMP
  else
  {
    for(int i = 0; i < numMat_ + 1; i++)
    {
      double tmp;
      if(i == matID2)
      {
        if(i == numMat_)
          tmp = 0.0;
        else
          tmp = -pp * pow(design_ss_[i], pp - 1.0);
      }
      else
      {
        if(i == numMat_)
          tmp = 1.0;
        else
          tmp = 1.0 - pow(design_ss_[i], pp);
      }

      if(i == 0)
        tmp *= 1.0;
      else
      {
        for(int j = 0; j < i; j++)
        {
          if(j == matID2)
            tmp *= pp * pow(design_ss_[j], pp - 1.0);
          else
            tmp *= pow(design_ss_[j], pp);
        }
      }

      dwds.push_back(tmp);
    }
  }
  return dwds[matID1];
}

double ElastMultiMat::dYoung(int matID2)
{
  double dYoung = 0.0;

  // void
  dYoung += dwds(0, matID2, pp_) * young1_;

  // solid
  for(int matID1 = 0; matID1 < numMat_; matID1++)
  {
    dYoung += dwds(matID1 + 1, matID2, pp_) * eachYoung_(matID1);
  }

  return dYoung;
}

void LinearElastMultiMat::makeKeFint(
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

  // get multi-material Young's modulus
  // void material
  auto coeff = makeCoeff(pp_);
  young_ = coeff[0] * young1_;

  // solid material
  young2_ = 0.0;
  for(int i = 0; i < numMat_; i++)
    young2_ += coeff[i + 1] * eachYoung_(i);

  // void & solid composite
  young_ += young2_;

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

} // namespace opt_multi_material
} // namespace icarat