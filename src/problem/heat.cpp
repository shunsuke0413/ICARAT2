///  @file  elast_heat.cpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "heat.hpp"
#include <basis/bmatrix.hpp>
#include <iostream>
#include <misc/assembling.hpp>

namespace icarat
{
using namespace std;
using namespace Eigen;

void Heat::setParameter(const toml::value &config)
{
  const auto &mat = toml::find(config, "material");
  assert(toml::find<string>(mat, "type") == "Heat");

  design_s_ = toml::find<double>(mat, "designs_0");
  pp_ = toml::find<double>(mat, "penalty");
  using pair = std::tuple<std::vector<double>, std::vector<double>>;
  auto conduct = toml::find<pair>(mat, "conductivity");
  auto conduct1 = get<0>(conduct);
  auto conduct2 = get<1>(conduct);
  auto specific = toml::find<std::array<double, 2>>(mat, "specific_heat");

  specificHeat_.push_back(get<0>(specific));
  specificHeat_.push_back(get<1>(specific));

  // 2D
  if(De1_.rows() == 2)
  {
    De1_ << conduct1[0], conduct1[3], conduct1[3], conduct1[1];
    De2_ << conduct2[0], conduct2[3], conduct2[3], conduct2[1];
  }
  // 3D
  else if(De1_.rows() == 3)
  {
    De1_ << conduct1[0], conduct1[3], conduct1[5], conduct1[3], conduct1[1],
        conduct1[4], conduct1[5], conduct1[4], conduct1[2];

    De1_ << conduct2[0], conduct2[3], conduct2[5], conduct2[3], conduct2[1],
        conduct2[4], conduct2[5], conduct2[4], conduct2[2];
  }
  else
  {
    cerr << "error in heat.cpp" << endl;
    exit(1);
  }
}

void Heat::makeSIMPDe()
{
  De_.setZero();
  if(De1_(0, 0) < De2_(0, 0))
  {
    cerr << "The thermal conductivity of material 1 needs to be greater than 2."
         << endl;
    exit(1);
  }

  De_ = pow(design_s_, pp_) * De1_ + (1.0 - pow(design_s_, pp_)) * De2_;
}

void LinearHeat::makeKeFint(vector<Triplet<double>> &tripletsK,
                            vector<Triplet<double>> &tripletsF,
                            vector<Triplet<double>> &tripletsR, FEM &fem,
                            vector<Node> &nodes)
{
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
  MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
  VectorXd temp = VectorXd::Zero(numdof);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  VectorXi idof = VectorXi::Zero(numdof);

  for(int i = 0; i < ne; i++)
  {
    temp.coeffRef(i) = nodes[nodeID[i]].val[0];
    idof.coeffRef(i) = nodes[nodeID[i]].dof[0];

    for(int j = 0; j < fem.ndim; j++)
      X.coeffRef(i, j) = nodes[nodeID[i]].x[j];
  }

  // make De matrix
  makeSIMPDe();

  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;
  double jac = 0.0;
  Bmatrix bmatrix;

  // gauss point loop start
  for(int ip = 0; ip < ipmax; ip++)
  {
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);

    fe_ += Be.transpose() * De_ * Be * temp * jac;
    Ke_ += Be.transpose() * De_ * Be * jac;
    volume += jac;
  } // gauss point loop end

  assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
  assemblingReact(numdof, idof, fe_, tripletsR);
}

} // namespace icarat