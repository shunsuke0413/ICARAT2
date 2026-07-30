///  @file  fluid.cpp
///  @author  Daiki Watanabe
///  @date  October 28, 2023.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "fluid.hpp"
#include <basis/bmatrix.hpp>
#include <io/boundary_condition_base.hpp>
#include <iostream>
#include <misc/assembling.hpp>

namespace icarat
{
using namespace std;
using namespace Eigen;

void Fluid::setParameter(const toml::value &config)
{
  const auto &mat = toml::find(config, "material");
  assert(toml::find<string>(mat, "type") == "Fluid");

  design_s_ = toml::find<double>(mat, "designs_0");
  pp_ = toml::find<double>(mat, "penalty");
  eta_ = toml::find<double>(mat, "kine_visco");
  rho_ = toml::find<double>(mat, "density");
}

void Stokes::makeKeFint(vector<Triplet<double>> &tripletsK,
                        vector<Triplet<double>> &tripletsF,
                        vector<Triplet<double>> &tripletsR, FEM &fem,
                        vector<Node> &nodes)
{
  tie(eTypep_, nep_) = reducedOrderEtypeNe(eType, ne);

  int numdofp = nep_;             // numdof of pressure
  int numdofu = numdof - numdofp; // total numdof
  VectorXi idof = VectorXi::Zero(numdof);

  // velocity variables
  MatrixXd Xu = MatrixXd::Zero(ne, fem.ndim);
  MatrixXd Neu = MatrixXd::Zero(fem.ndim, numdofu);
  MatrixXd Beu = MatrixXd::Zero(fem.voigt, numdofu);
  VectorXd u = VectorXd::Zero(numdofu);

  // pressure variables
  MatrixXd Xp = MatrixXd::Zero(nep_, fem.ndim);
  MatrixXd Nep = MatrixXd::Zero(1, numdofp);
  MatrixXd Bep = MatrixXd::Zero(fem.voigt, fem.ndim * numdofp); // dummy
  VectorXd p = VectorXd::Zero(numdofp);

  // velocity part
  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
    {
      Xu.coeffRef(i, j) = nodes[nodeID[i]].x[j];
      u.coeffRef(fem.ndim * i + j) = nodes[nodeID[i]].val[j];
      idof.coeffRef(fem.ndim * i + j) = nodes[nodeID[i]].dof[j];
    }
  }
  // pressure part
  for(int i = 0; i < nep_; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
      Xp.coeffRef(i, j) = nodes[nodeID[i]].x[j];

    p.coeffRef(i) = nodes[nodeID[i]].val[fem.ndim];
    idof.coeffRef(numdofu + i) = nodes[nodeID[i]].dof[fem.ndim];
  }

  Ke_.setZero();
  fe_.setZero();
  volume = 0.0;
  double jac = 0.0, dummy = 0.0;
  Bmatrix bmatrix;

  De_ = MatrixXd::Identity(fem.voigt, fem.voigt);
  De_ *= 2.0 * eta_;

  VectorXd m(fem.voigt);
  if(fem.ndim == 2)
    m << 1.0, 1.0, 0.0;
  else // 3D
    m << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;

  // gauss point loop start
  for(int ip = 0; ip < ipmax; ip++)
  {
    MatrixXd Ketmp = MatrixXd::Zero(numdof, numdof);
    VectorXd fetmp = VectorXd::Zero(numdof);

    bmatrix.make(Beu, Neu, jac, eType, ipmax, Xu, ip);
    MatrixXd Keu = Beu.transpose() * De_ * Beu * jac;
    VectorXd feu = Beu.transpose() * De_ * Beu * u * jac;
    Ketmp.block(0, 0, numdofu, numdofu) = Keu;
    fetmp.block(0, 0, numdofu, 1) = feu;

    bmatrix.make(Bep, Nep, dummy, eTypep_, ipmax, Xp, ip);
    MatrixXd Kep = -Beu.transpose() * m * Nep * jac;
    VectorXd fep = -Beu.transpose() * m * Nep * u * jac;
    Ketmp.block(0, numdofu, numdofu, numdofp) = Kep;
    Ketmp.block(numdofu, 0, numdofp, numdofu) = Kep.transpose();
    fetmp.block(numdofu, 0, numdofp, 1) = fep;

    Ke_ += Ketmp;
    fe_ += fetmp;
    volume += jac;
  } // gauss point loop end

  assembling(numdof, fem.numeq, idof, Ke_, fe_, tripletsK, tripletsF);
  assemblingReact(numdof, idof, fe_, tripletsR);
}

} // namespace icarat