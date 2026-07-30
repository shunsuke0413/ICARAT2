///  @file  filter_helm.hpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include "opt_base.hpp"
#include <analysis/solver/solvers.hpp>
#include <basis/bmatrix.hpp>
#include <iostream>
#include <misc/assembling.hpp>
#include <problem/base.hpp>
#include <vector>

namespace icarat
{
/// Helmholtz filters in TO
/// This method is preferred for complex geometries or large-scale models.
/// @sa B.S. Lazarov, O. Sigmund, Filters in topology optimization based on
/// Helmholtz-type differential equations, 2010.
template <class E, class N> class FilterHelmholtz
{
public:
  /// This constructor needs
  /// @param [in] radious_: filter radious_ which effects surrounding elements
  FilterHelmholtz(FEM &fem, std::vector<E> &element, std::vector<N> &node,
                  double radious_);

  /// apply density filter
  /// @returns filtered design_s and sensitivity derivative term of design_s
  /// @param  [in] sensFlag true...sensitivity derivative term,
  /// false...density filter procedure
  std::vector<double> densityFilter(const std::vector<double> &design_s,
                                    std::string solvName, int log,
                                    bool sensFlag);

  /// update design variable using threshold function
  std::vector<double> threshold(const std::vector<double> &design_s,
                                double beta, double T);

  /// differential of threshold function
  std::vector<double>
  thresholdDerivative(const std::vector<double> &design_filtered, double beta,
                      double T);

  /// setter(This is used for shell structure)
  void setNdim(int n) { ndim_ = n; }

private:
  void makeKeFint(int nel, double design_s,
                  std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  bool sensFlag);

  ///@return filtered design variable
  double divideToElement(int nel, Eigen::VectorXd &val, bool sensFlag);

  FEM &fem;
  std::vector<E> &element;
  std::vector<N> &node;
  double radious_; ///< filter radious
  double ndim_;    ///< dimension
};

template <class E, class N>
FilterHelmholtz<E, N>::FilterHelmholtz(FEM &fem, std::vector<E> &element,
                                       std::vector<N> &node, double radious_)
    : fem(fem), element(element), node(node), radious_(radious_)
{
  ndim_ = fem.ndim;
}

template <class E, class N>
std::vector<double>
FilterHelmholtz<E, N>::densityFilter(const std::vector<double> &design_s,
                                     std::string solvName, int log,
                                     bool sensFlag)
{
  if(log == 1)
    std::cout << "Helmholtz filtering..." << std::endl;

  Eigen::SparseMatrix<double> globalK(fem.numnp, fem.numnp),
      globalF(fem.numnp, fem.numnp);
  std::vector<Eigen::Triplet<double>> tripletsK, tripletsF;
  Eigen::VectorXd val = Eigen::VectorXd::Zero(fem.numnp);

#define EIGEN_DONT_PARALLELIZE

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    /// private variables
    std::vector<Eigen::Triplet<double>> tripK_pri;
    std::vector<Eigen::Triplet<double>> tripF_pri;

#ifdef _OPENMP
#pragma omp for nowait
#endif
    /// make element matrix & force
    for(int nel = 0; nel < fem.nelx; nel++)
      makeKeFint(nel, design_s[nel], tripK_pri, tripF_pri, sensFlag);

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      /// connect
      tripletsK.insert(tripletsK.end(), tripK_pri.begin(), tripK_pri.end());
      tripletsF.insert(tripletsF.end(), tripF_pri.begin(), tripF_pri.end());
    }
  }

#undef EIGEN_DONT_PARALLELIZE

  globalK.setFromTriplets(tripletsK.begin(), tripletsK.end());
  globalK.makeCompressed();

  globalF.setFromTriplets(tripletsF.begin(), tripletsF.end());
  Eigen::VectorXd fint = globalF.diagonal();

  /// solve KU=F for design_s filter
  Solvers solv(globalK, fint, val);
  solv.solve(solvName, log);

  std::vector<double> filteredDensity((int)design_s.size());

/// Allocate to each element
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int nel = 0; nel < fem.nelx; nel++)
    filteredDensity[nel] = divideToElement(nel, val, sensFlag);

  return filteredDensity;
}

template <class E, class N>
std::vector<double>
FilterHelmholtz<E, N>::threshold(const std::vector<double> &design_filtered,
                                 double beta, double T)
{
  std::vector<double> hat(fem.nelx);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int nel = 0; nel < fem.nelx; nel++)
    hat[nel] =
        (0.5 * (tanh(T * beta) + tanh(beta * (design_filtered[nel] - T))) /
         tanh(T * beta));

  return hat;
}

template <class E, class N>
std::vector<double> FilterHelmholtz<E, N>::thresholdDerivative(
    const std::vector<double> &design_filtered, double beta, double T)
{
  std::vector<double> dshat(fem.nelx);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int nel = 0; nel < fem.nelx; nel++)
    dshat[nel] =
        (0.5 *
         (beta * (1.0 - pow(tanh(beta * (design_filtered[nel] - T)), 2.0))) /
         tanh(T * beta));

  return dshat;
}

/////////////////////////////////////////////////////////
//////////////////////////private////////////////////////
/////////////////////////////////////////////////////////
template <class E, class N>
void FilterHelmholtz<E, N>::makeKeFint(
    int nel, double design_s, std::vector<Eigen::Triplet<double>> &tripletsK,
    std::vector<Eigen::Triplet<double>> &tripletsF, bool sensFlag)
{
  Eigen::MatrixXd Ne = Eigen::MatrixXd::Zero(ndim_, ndim_ * element[nel].ne);
  Eigen::VectorXd NeArr = Eigen::VectorXd::Zero(element[nel].ne);
  Eigen::MatrixXd Be = Eigen::MatrixXd::Zero(ndim_, element[nel].ne);
  /// force for Helmholtz design_s filter
  Eigen::VectorXd fe = Eigen::VectorXd::Zero(element[nel].ne);
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(element[nel].ne, ndim_);
  Eigen::VectorXi idof = Eigen::VectorXi::Zero(element[nel].ne);
  Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(element[nel].ne, element[nel].ne);
  Eigen::MatrixXd IKe = Eigen::MatrixXd::Identity(ndim_, ndim_);
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(ndim_);

  for(int i = 0; i < element[nel].ne; i++)
  {
    for(int j = 0; j < ndim_; j++)
      X.coeffRef(i, j) = node[element[nel].nodeID[i]].x[j];

    idof.coeffRef(i) = element[nel].nodeID[i];
  }

  double jac = 0.0;
  Bmatrix bmatrix;
  IKe *= pow(radious_ / (2 * sqrt(3)), 2);

  // gauss point loop start
  for(int ip = 0; ip < element[nel].ipmax; ip++)
  {
    bmatrix.make(Be, Ne, jac, element[nel].eType, element[nel].ipmax, X, ip);

    for(int n = 0; n < element[nel].ne; n++)
      NeArr.coeffRef(n) = Ne.coeff(0, ndim_ * n);

    Ke += (NeArr * NeArr.transpose() + Be.transpose() * IKe * Be) * jac;
    fe += design_s * NeArr * jac;
  } // gauss point loop end

  if(sensFlag == true)
    fe = fe / element[nel].ipmax;

  assembling(element[nel].ne, fem.numnp, idof, Ke, fe, tripletsK, tripletsF);
}

template <class E, class N>
double FilterHelmholtz<E, N>::divideToElement(int nel, Eigen::VectorXd &val,
                                              bool sensFlag)
{
  Eigen::MatrixXd NeArr = Eigen::MatrixXd::Zero(1, element[nel].ne);
  Eigen::MatrixXd Be = Eigen::MatrixXd::Zero(ndim_, element[nel].ne);
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(element[nel].ne, ndim_);
  Eigen::VectorXd valEle = Eigen::VectorXd::Zero(element[nel].ne);

  for(int i = 0; i < element[nel].ne; i++)
  {
    for(int j = 0; j < ndim_; j++)
      X.coeffRef(i, j) = node[element[nel].nodeID[i]].x[j];

    valEle.coeffRef(i) = val.coeff(element[nel].nodeID[i]);
  }

  Bmatrix bmatrix;
  double dummy, filteredDens = 0.0;
  // gauss point loop start
  for(int ip = 0; ip < element[nel].ipmax; ip++)
  {
    bmatrix.make(Be, NeArr, dummy, element[nel].eType, element[nel].ipmax, X,
                 ip);

    filteredDens += (NeArr * valEle).coeff(0, 0);
  } // gauss point loop end

  if(sensFlag == false)
    filteredDens = filteredDens / element[nel].ipmax;

  return filteredDens;
}

} // namespace icarat