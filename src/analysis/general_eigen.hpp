///  @file  general_eigen.hpp
///  @author  Daiki Watanabe
///  @date  January 16, 2022.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php
#pragma once
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SymShiftInvert.h>
#include <Spectra/SymGEigsShiftSolver.h>
#include <problem/base.hpp>

namespace icarat
{
using namespace Spectra;

/// FEM analysis class for general eigenvalue problem A x =\lambda B x
/// with using Shifted-Invert Lanczos method.
/// result will be stored to "ValueEigen" class
template <class E, class N> class GeneralEigen
{
public:
  GeneralEigen(FEM &fem, std::vector<E> &element, std::vector<N> &node,
               ValueEigen &value, int numeigen, int paraconv, double sigma);

  /// main method for eigen problem
  void solve();

  ///   initializaiton
  void initialization();

  /// assembling A and B matrix
  void assemblingAB();

  /// getters
  const Eigen::SparseMatrix<double> &A() const { return globalA; }
  const Eigen::SparseMatrix<double> &B() const { return globalB; }
  const int &numeigen() const { return numeigen_; }

private:
  FEM &fem;
  std::vector<E> &element;
  std::vector<N> &node;
  ValueEigen &valueE;
  Eigen::SparseMatrix<double> globalA; ///< global stiffness matrix
  Eigen::SparseMatrix<double> globalB; ///< global geometory stiffness matrix
  int numeigen_;                       ///< number of evaluated eigenvalues
  int paraconv_; ///< convergence parameter (paraconv_≥2⋅numeigen_)
  double sigma_; ///< start value in  Lanczos algorithm
};

template <class E, class N>
GeneralEigen<E, N>::GeneralEigen(FEM &fem, std::vector<E> &element,
                                 std::vector<N> &node, ValueEigen &value,
                                 int numeigen, int paraconv, double sigma)
    : fem(fem), element(element), node(node), valueE(value),
      globalA(fem.numeq, fem.numeq), globalB(fem.numeq, fem.numeq),
      numeigen_(numeigen), paraconv_(paraconv), sigma_(sigma)
{
}

template <class E, class N> void GeneralEigen<E, N>::solve()
{
  std::cout << "--------- general eigenvalue analysis ---------" << std::endl;

  // init
  initialization();

  // set dirichlet boundary condition
  for(auto &n : node)
  {
    if(n.dirich == nullptr)
      continue;
    for(int j = 0; j < 3; j++)
      n.val[j] = n.dirich->val[j];
  }

  assemblingAB();

  using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
  using BOpType = SparseSymMatProd<double>;
  OpType op(globalA, globalB);
  BOpType Bop(globalA);

  // Construct generalized eigen solver object, seeking three generalized
  // eigenvalues that are closest to and larger than 1.0. This is equivalent to
  // specifying a shift sigma = 1.0 combined with the SortRule::LargestAlge
  // selection rule
  // op...The matrix operation object that computes 𝑦=(𝐾−𝜎𝐾𝐺)−1𝑣 for any vector
  // 𝑣.
  // Users could either create the object from the wrapper class SymShiftInvert,
  // or define their own that implements all the public members as in
  // SymShiftInvert.
  // Bop...The 𝐾 matrix operation object that implements the matrix-vector
  // multiplication 𝐾𝑣.
  // Users could either create the object from the wrapper classes such as
  // DenseSymMatProd and SparseSymMatProd, or define their own that implements
  // all the public member functions as in DenseSymMatProd. 𝐾 needs to be
  // positive definite.
  // nev...Number of eigenvalues requested. This should satisfy 1≤𝑛𝑒𝑣≤𝑛−1, where
  // 𝑛 is the size of matrix. ncv...Parameter that controls the convergence
  // speed of the algorithm. Typically a larger ncv means faster convergence,
  // but it may also result in greater memory use and more matrix operations in
  // each iteration. This parameter must satisfy 𝑛𝑒𝑣<𝑛𝑐𝑣≤𝑛, and is advised to
  // take 𝑛𝑐𝑣≥2⋅𝑛𝑒𝑣. sigma...The value of the shift.
  int nev = numeigen_;
  int ncv = paraconv_;
  SymGEigsShiftSolver<OpType, BOpType, GEigsMode::Buckling> geigs(op, Bop, nev,
                                                                  ncv, sigma_);

  // Initialize and compute
  geigs.init();
  int nconv = geigs.compute(SortRule::LargestAlge);

  // Retrieve results to ValueEigen class
  if(geigs.info() == CompInfo::Successful)
  {
    for(int i = 0; i < numeigen_; i++)
    {
      valueE.evalues[numeigen_ - i - 1] = geigs.eigenvalues()(i);
      valueE.evectors[numeigen_ - i - 1] = geigs.eigenvectors().col(i);

      // for(int j = 0; j < fem.numeq; j++)
      // {
      //   // too small value is omitted for visualization
      //   if(abs(valueE.evectors[i].coeff(j)) < 1.0e-20)
      //     valueE.evectors[i].coeffRef(j) = 0.0;
      // }
    }
  }
  else
  {
    std::cerr << "Eigen value computation failed in linear_buckling.hpp."
              << std::endl;
    exit(1);
  }

  std::cout << "eigen values: " << geigs.eigenvalues().transpose() << std::endl;

  return;
}

template <class E, class N> void GeneralEigen<E, N>::initialization()
{
  // ValueEigen class
  for(int i = 0; i < (int)valueE.evalues.size(); i++)
  {
    valueE.evalues[i] = 0.0;
    valueE.evectors[i].setZero();
  }
}

template <class E, class N> void GeneralEigen<E, N>::assemblingAB()
{
  std::cout << "assembling...";
  fflush(stdout);

  globalA.setZero();
  globalB.setZero();
  std::vector<Eigen::Triplet<double>> tripletsA, tripletsB;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    /// private variables
    std::vector<Eigen::Triplet<double>> tripA_pri, tripB_pri;

#ifdef _OPENMP
#pragma omp for nowait
#endif
    for(int nel = 0; nel < fem.nelx; nel++)
      element[nel].makeEigenCoeffs(tripA_pri, tripB_pri, fem, node);

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      /// connect
      tripletsA.insert(tripletsA.end(), tripA_pri.begin(), tripA_pri.end());
      tripletsB.insert(tripletsB.end(), tripB_pri.begin(), tripB_pri.end());
    }
  }

  globalA.setFromTriplets(tripletsA.begin(), tripletsA.end());
  globalB.setFromTriplets(tripletsB.begin(), tripletsB.end());

  globalA.makeCompressed();
  globalB.makeCompressed();

  std::cout << "done" << std::endl;
  fflush(stdout);
}

} // namespace icarat