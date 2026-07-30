///  @file  solvers.cpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php
#include "solvers.hpp"
#include <Eigen/Sparse>
#ifdef __INTEL_LLVM_COMPILER
#include <Eigen/PardisoSupport>
#endif
namespace icarat
{

Solvers::Solvers(Eigen::SparseMatrix<double> const &A_, Eigen::VectorXd &rhs,
                 Eigen::VectorXd &sol)
    : A(A_), rhs(rhs), sol(sol)
{
}

Solvers::Solvers(Eigen::SparseMatrix<double> const &A_, Eigen::VectorXd &rhs,
                 Eigen::VectorXd &sol, double tol, int maxIter)
    : A(A_), rhs(rhs), sol(sol), tol_(tol), maxIter_(maxIter)
{
}

void Solvers::solve(std::string solvName, int log)
{
  if(solvName == "CG")
    CG(log);
  else if(solvName == "Scaling_CG")
    ScalingCG(log);
  else if(solvName == "ICCG")
    ICCG(log);
  /// Eigen's solvers
  else if(solvName == "Eigen_CG")
    EigenCG(log);
  else if(solvName == "Eigen_ICCG")
    EigenICCG(log);
  else if(solvName == "Eigen_BiCGSTAB")
    EigenBiCGSTAB(log);
  else if(solvName == "Eigen_SparseLU")
    EigenSparseLU(log);
  else if(solvName == "Eigen_LDLT")
    EigenLDLT(log);
#ifdef __INTEL_LLVM_COMPILER
  else if(solvName == "Pardiso_LDLT")
    PardisoLDLT(log);

  else if(solvName == "Pardiso_LLT")
    PardisoLLT(log);

  else if(solvName == "Pardiso_LU")
    PardisoLU(log);
#endif
  else
  {
    std::cerr << "Error ! Solver type name is wrong." << std::endl;
    exit(1);
  }
}

void Solvers::EigenCG(int log)
{
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper,
                           Eigen::DiagonalPreconditioner<double>>
      solver;
  solver.setTolerance(tol_);
  solver.setMaxIterations(maxIter_);
  solver.compute(A);

  if(log == 1)
  {
    if(solver.info() != Eigen::Success)
      std::cerr << "decomposition failed" << std::endl;
    else
      std::cerr << "decomposition succeeded" << std::endl;
  }
  sol = solver.solve(rhs);
  if(log == 1)
  {
    if(solver.info() != Eigen::Success)
      std::cerr << "solving failed" << std::endl;
    else
      std::cerr << "solving succeeded" << std::endl;

    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error() << std::endl;
  }
}

void Solvers::EigenICCG(int log)
{
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper,
                           Eigen::IncompleteCholesky<double>>
      solver;
  solver.setTolerance(tol_);
  solver.setMaxIterations(maxIter_);
  solver.compute(A);

  if(log == 1)
  {
    if(solver.info() != Eigen::Success)
      std::cerr << "decomposition failed" << std::endl;
    else
      std::cerr << "decomposition succeeded" << std::endl;
  }
  sol = solver.solve(rhs);
  if(log == 1)
  {
    if(solver.info() != Eigen::Success)
      std::cerr << "solving failed" << std::endl;
    else
      std::cerr << "solving succeeded" << std::endl;

    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error() << std::endl;
  }
}

void Solvers::EigenBiCGSTAB(int log)
{
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
  solver.setTolerance(tol_);
  solver.setMaxIterations(maxIter_);
  solver.compute(A);

  if(log == 1)
  {
    if(solver.info() != Eigen::Success)
      std::cerr << "decomposition failed" << std::endl;
    else
      std::cerr << "decomposition succeeded" << std::endl;
  }

  sol = solver.solve(rhs);

  if(log == 1)
  {
    if(solver.info() != Eigen::Success)
      std::cerr << "solving failed" << std::endl;
    else
      std::cerr << "solving succeeded" << std::endl;

    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error() << std::endl;
  }
}

void Solvers::EigenSparseLU(int log)
{
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;
  solver.compute(A);

  if(log == 1)
  {
    if(solver.info() != Eigen::Success)
      std::cerr << "decomposition failed" << std::endl;
    else
      std::cerr << "decomposition succeeded" << std::endl;
  }

  sol = solver.solve(rhs);

  if(log == 1)
  {
    if(solver.info() != Eigen::Success)
      std::cerr << "solving failed" << std::endl;
    else
      std::cerr << "solving succeeded" << std::endl;
  }
}

void Solvers::EigenLDLT(int log)
{
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);

  if(log == 1)
  {
    if(solver.info() != Eigen::Success)
      std::cerr << "decomposition failed" << std::endl;
    else
      std::cerr << "decomposition succeeded" << std::endl;
  }

  sol = solver.solve(rhs);

  if(log == 1)
  {
    if(solver.info() != Eigen::Success)
      std::cerr << "solving failed" << std::endl;
    else
      std::cerr << "solving succeeded" << std::endl;
  }
}

#ifdef __INTEL_LLVM_COMPILER
void Solvers::PardisoLDLT(int log)
{
  Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);

  if(log == 1)
  {
    if(solver.info() != Eigen::Success)
      std::cerr << "decomposition failed" << std::endl;
  }
  sol = solver.solve(rhs);

  if(log == 1)
  {
    if(solver.info() != Eigen::Success)
      std::cerr << "solving failed" << std::endl;
  }
}

void Solvers::PardisoLLT(int log)
{
  Eigen::PardisoLLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);

  if(log == 1)
  {
    if(solver.info() != Eigen::Success)
      std::cerr << "decomposition failed" << std::endl;
  }

  sol = solver.solve(rhs);

  if(log == 1)
  {
    if(solver.info() != Eigen::Success)
      std::cerr << "solving failed" << std::endl;
  }
}

void Solvers::PardisoLU(int log)
{
  Eigen::PardisoLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);

  if(log == 1)
  {
    if(solver.info() != Eigen::Success)
      std::cerr << "decomposition failed" << std::endl;
  }

  sol = solver.solve(rhs);

  if(log == 1)
  {
    if(solver.info() != Eigen::Success)
      std::cerr << "solving failed" << std::endl;
  }
}
#endif
} // namespace icarat
