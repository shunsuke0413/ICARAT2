///  @file  solvers.hpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
///   You can choose solver from Eigen's solver or icarat's solver
///  The lists of Eigen's solver and how to use is in the below URL.
///  @sa http://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <iostream>
#include <problem/base.hpp>
namespace icarat
{
class Solvers
{
public:
  /// constructor
  Solvers(Eigen::SparseMatrix<double> const &K, Eigen::VectorXd &rhs,
          Eigen::VectorXd &sol);

  Solvers(Eigen::SparseMatrix<double> const &K, Eigen::VectorXd &rhs,
          Eigen::VectorXd &sol, double tol, int maxIter);

  /// do computing solver
  /// For a list of "solvName" arguments, see the comment out below.
  ///"log" argument = 1 ... display log
  ///"log" argument = 0 (≠1) ... not display log
  void solve(std::string solvName, int log);

  /// classic CG method (solvName="CG")
  /// slow & only symmetric matrix
  void CG(int log);

  /// diagonal scaling CG method (solvName="Scaling_CG")
  /// fast in topology optimization & only symmetric matrix
  void ScalingCG(int log);

  /// ILU0CG method (solvName="ICCG")
  /// slow & only symmetric matrix
  void ICCG(int log);

  /// Eigen's CG method (solvName="Eigen_CG")
  /// Eigen's default pre-conditioner is diagonal
  /// fast & only symmetric matrix & openmp is available
  void EigenCG(int log);

  /// Eigen's CG method (solvName="Eigen_ICCG")
  /// Eigen's default pre-conditioner is diagonal
  /// fast & only symmetric matrix & openmp is available
  void EigenICCG(int log);

  /// Eigen's A bi conjugate gradient stabilized solver for sparse square
  /// problems. (solvName="Eigen_BiCGSTAB")
  /// slow & non-symmetric matrix is also ok & openmp is available
  void EigenBiCGSTAB(int log);

  /// Eigen's LU method (solvName="Eigen_SparseLU")
  /// slow & any matrix is OK.
  void EigenSparseLU(int log);

  /// Eigen's LDLT method (solvName="Eigen_LDLT")
  /// slow & any matrix is OK.
  void EigenLDLT(int log);

  /// Cholesky decomposition (solvName="Pardiso_LDLT")
  /// so fast with intel MKL, and it must be real positive definite matrix
#ifdef __INTEL_LLVM_COMPILER
  void PardisoLDLT(int log);
#endif

  /// Cholesky decomposition (solvName="Pardiso_LLT")
  /// fast with intel MKL, and it must be positive definite matrix
#ifdef __INTEL_LLVM_COMPILER
  void PardisoLLT(int log);
#endif

  /// A sparse direct LU factorization and solver based on the PARDISO
  /// library. (solvName="Pardiso_LU")
#ifdef __INTEL_LLVM_COMPILER
  void PardisoLU(int log);
#endif

private:
  Eigen::SparseMatrix<double> const &A;
  Eigen::VectorXd &rhs;  ///< residual force
  Eigen::VectorXd &sol;  ///< solution value
  double tol_ = 1.0e-10; ///< torelance in iterative method
  int maxIter_ = 10000;  ///< max iteration in iterative method
};
} // namespace icarat