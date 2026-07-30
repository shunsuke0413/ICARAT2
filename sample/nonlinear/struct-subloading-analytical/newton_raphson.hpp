///  @file  newton_raphson.hpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php
#pragma once
#include <analysis/analysis_utility.hpp>
#include <analysis/solver/solvers.hpp>
#include <problem/base.hpp>
#include <Eigen/SparseCholesky>

namespace icarat
{
/// FEM procedure for linear elasticity
template <class E, class N> class NewtonRaphson
{
public:
  NewtonRaphson(FEM &fem, std::vector<E> &element, std::vector<N> &node,
                Force &force, Value &value);

  ///  main function (Incremental load or displacement method)
  ///@return Newton-Raphson iteration
  int solve(int thisTime, int totalTime, double tolNR, int maxIter,
            std::string solvName, int log);

  int solve_M(int thisTime, int totalTime, double tolNR, int maxIter,
              std::string solvName, int log);

  int solve_M2(int thisTime, int totalTime, double tolNR, int maxIter,
              std::string solvName, int log);

  int solve_cyclic(int thisTime, int totalTime, double tolNR, int maxIter,
                   std::string solvName, int log, double sign);

  /// main function (arc length method)
  void solveArcLength(int thisStep, int totalStep, std::string solvName,
                      int log);

  /// initializaiton
  void initialization();

  double calculateTotalReaction(int direction, double target_x, double target_y, double target_z);

  /// getters
  const Eigen::SparseMatrix<double> &K() const { return globalK; }

private:
  void assemblingK(int log);

  void divideValIntoNode(int dofnp, Eigen::VectorXd &val, Eigen::VectorXd &du,
                         int thisStep, int totalStep);

  void integrateValFromNode(int dofnp, Eigen::VectorXd &val,
                            Eigen::VectorXd &du);

  void divideValIntoNodeCyclic(int dofnp, Eigen::VectorXd &val,
                               Eigen::VectorXd &du, int thisStep, int totalStep,
                               double sign);

  FEM &fem;
  std::vector<E> &element;
  std::vector<N> &node;
  Force &force;
  Value &value;
  Eigen::SparseMatrix<double> globalK; ///< global stiffness matrix
  Eigen::VectorXd reaction;            ///< reaction force

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_;
  bool is_solver_analyzed_ = false;

  std::vector<Eigen::Triplet<double>> tripletsK_;
  std::vector<Eigen::Triplet<double>> tripletsF_;
  std::vector<Eigen::Triplet<double>> tripletsR_;
};

template <class E, class N>
NewtonRaphson<E, N>::NewtonRaphson(FEM &fem, std::vector<E> &element,
                                   std::vector<N> &node, Force &force,
                                   Value &value)
    : fem(fem), element(element), node(node), force(force), value(value),
      globalK(fem.numeq, fem.numeq), reaction(fem.neq)
{
  std::cout << "Element class name is " << typeid(E).name() << std::endl;
  std::cout << "Node class name is " << typeid(N).name() << std::endl;

  is_solver_analyzed_ = false;

  // Reserve memory for triplets to avoid reallocations
  // Estimated size: element count * (dof per element)^2
  // Roughly estimating dof per element as 24 (3D hex) -> 24*24 = 576
  // Allocating a safe upper bound to prevent reallocations during push_back
  long long estimated_nnz = (long long)fem.nelx * 576; 
  if (estimated_nnz > 0) {
      tripletsK_.reserve(estimated_nnz);
      tripletsF_.reserve(fem.nelx * 24); 
      tripletsR_.reserve(fem.nelx * 24);
  }
}

template <class E, class N>
int NewtonRaphson<E, N>::solve(int thisTime, int totalTime, double tolNR,
                               int maxIter, std::string solvName, int log)
{

  int iterstep = 0;
  double firstNorm, residual = 1;
  Eigen::VectorXd sol = Eigen::VectorXd::Zero(fem.numeq);
  Eigen::VectorXd du = Eigen::VectorXd::Zero(fem.numeq);

  // set dirichlet boundary condition
  for(auto &n : node)
  {
    if(n.dirich != nullptr)
      for(int j = 0; j < 3; j++)
      {
        if(n.dirich->flag[j] == 1)
        {
          n.val[j] = thisTime * n.dirich->val[j] / totalTime;
          n.du[j] = n.dirich->val[j] / totalTime;
        }
        else
          n.du[j] = 0.0;
      }
    else
      for(int j = 0; j < 3; j++)
      {
        n.du[j] = 0.0;
      }
  }

  /////////////////////////////////////////////////////////////////
  //////////////////////////start NR loop//////////////////////////
  /////////////////////////////////////////////////////////////////
  std::cout << "------------ NR analysis ------------" << std::endl;
  std::cout << " iter "
            << " first norm "
            << "      res " << std::endl;

  while(residual > tolNR)
  {
    iterstep++;

    // assembling
    assemblingK(log);

    // cal residual force
    force.fres = force.fext - force.fint;

    /// convergence check
    double norm = sqrt(force.fres.dot(force.fres));

    if(iterstep == 1)
      firstNorm = norm;

    residual = norm / firstNorm;

    std::cout << "   " << iterstep << "      " << firstNorm << "        "
              << residual << std::endl;

    if(iterstep == maxIter)
    {
      std::cerr << "itaration step reach " << maxIter
                << " in newton_raphson.hpp" << std::endl;
      return iterstep;
    }

    sol.setZero();
    Solvers solv(fem, globalK, force.fres, sol);
    solv.solve(solvName, log);

    // add corrector
    value.val += sol;
    du += sol;

    // make result to node
    divideValIntoNode(fem.dofnp, value.val, du, thisTime, totalTime);
    divideReactionIntoNode(fem.dofnp, reaction, node);
  } ///< NR loop end

  // make result to node
  divideValIntoNode(fem.dofnp, value.val, du, thisTime, totalTime);
  divideReactionIntoNode(fem.dofnp, reaction, node);

  return iterstep;
}

// Using previous displacement increment as initial value of Newton-Raphson
// written by Fujiwara
template <class E, class N>
int NewtonRaphson<E, N>::solve_M(int thisTime, int totalTime, double tolNR,
                                 int maxIter, std::string solvName, int log)
{

  int iterstep = 0;
  double firstNorm, residual = 1;
  Eigen::VectorXd sol = Eigen::VectorXd::Zero(fem.numeq);
  Eigen::VectorXd du = Eigen::VectorXd::Zero(fem.numeq);

  // set dirichlet boundary condition for firstNorm
  for(auto &n : node)
  {
    for(int j = 0; j < 3; j++)
    {
      n.du_old[j] = n.du[j];
    }
    if(n.dirich != nullptr)
      for(int j = 0; j < 3; j++)
      {
        if(n.dirich->flag[j] == 1)
        {
          n.val[j] = thisTime * n.dirich->val[j] / totalTime;
          n.du[j] = n.dirich->val[j] / totalTime;
        }
        else
          n.du[j] = 0.0;
      }
    else
      for(int j = 0; j < 3; j++)
      {
        n.du[j] = 0.0;
      }
  }

  // assembling
  assemblingK(log);

  // cal residual force
  force.fres = force.fext - force.fint;

  // cal first norm
  firstNorm = sqrt(force.fres.dot(force.fres));

  /////////////////////////////////////////////////////////////////
  //////////////////////////start NR loop//////////////////////////
  /////////////////////////////////////////////////////////////////
  std::cout << "------------ NR analysis ------------" << std::endl;
  std::cout << " iter "
            << " first norm "
            << "      res " << std::endl;


  if(thisTime == 1)
  {
        std::cout << "   " << iterstep << "      " << firstNorm << "        "
              << residual << std::endl;
    sol.setZero();
    Solvers solv(fem, globalK, force.fres, sol);
    solv.solve(solvName, log);

    // add corrector
    value.val += sol;
    du += sol;

    // make result to node
    divideValIntoNode(fem.dofnp, value.val, du, thisTime, totalTime);
    divideReactionIntoNode(fem.dofnp, reaction, node);
  }
  else if(thisTime > 1)
  {
    // set dirichlet boundary condition
    for(auto &n : node)
    {
      if(n.dirich != nullptr)
        for(int j = 0; j < 3; j++)
        {
          if(n.dirich->flag[j] == 1)
          {
            n.val[j] = thisTime * n.dirich->val[j] / totalTime;
            n.du[j] = n.dirich->val[j] / totalTime;
          }
          else
          {
            n.du[j] = n.du_old[j];
            n.val[j] += n.du[j];
          }
        }
      else
        for(int j = 0; j < 3; j++)
        {
          n.du[j] = n.du_old[j];
          n.val[j] += n.du[j];
        }
    }
    integrateValFromNode(fem.dofnp, value.val, du);
  }


  while(residual > tolNR)
  {
    iterstep++;

    // assembling
    assemblingK(log);

    // cal residual force
    force.fres = force.fext - force.fint;

    /// convergence check
    double norm = sqrt(force.fres.dot(force.fres));

    residual = norm / firstNorm;

    std::cout << "   " << iterstep << "      " << firstNorm << "        "
              << residual << std::endl;

    if(iterstep == maxIter)
    {
      std::cerr << "itaration step reach " << maxIter
                << " in newton_raphson.hpp" << std::endl;
      return iterstep;
    }

    sol.setZero();
    Solvers solv(fem, globalK, force.fres, sol);
    solv.solve(solvName, log);

    // add corrector
    value.val += sol;
    du += sol;

    // make result to node
    divideValIntoNode(fem.dofnp, value.val, du, thisTime, totalTime);
    divideReactionIntoNode(fem.dofnp, reaction, node);
  } ///< NR loop end

  // make result to node
  divideValIntoNode(fem.dofnp, value.val, du, thisTime, totalTime);
  divideReactionIntoNode(fem.dofnp, reaction, node);

  return iterstep;
}

template <class E, class N>
int NewtonRaphson<E, N>::solve_M2(int thisTime, int totalTime, double tolNR,
                                 int maxIter, std::string solvName, int log)
{
  int iterstep = 0;
  double firstNorm, residual = 1;
  double norm = 0.0;
  Eigen::VectorXd sol = Eigen::VectorXd::Zero(fem.numeq);
  Eigen::VectorXd du = Eigen::VectorXd::Zero(fem.numeq);

  // set dirichlet boundary condition for firstNorm
  for(auto &n : node)
  {
    for(int j = 0; j < 3; j++)
    {
      n.du_old[j] = n.du[j];
    }
    if(n.dirich != nullptr)
      for(int j = 0; j < 3; j++)
      {
        if(n.dirich->flag[j] == 1)
        {
          n.val[j] = thisTime * n.dirich->val[j] / totalTime;
          n.du[j] = n.dirich->val[j] / totalTime;
        }
        else
          n.du[j] = 0.0;
      }
    else
      for(int j = 0; j < 3; j++)
      {
        n.du[j] = 0.0;
      }
  }

  // assembling
  assemblingK(log);

  if (!is_solver_analyzed_) {
      solver_.analyzePattern(globalK);
      is_solver_analyzed_ = true;
  }

  // cal residual force
  force.fres = force.fext - force.fint;

  // cal first norm
  firstNorm = sqrt(force.fres.dot(force.fres));

  /////////////////////////////////////////////////////////////////
  //////////////////////////start NR loop//////////////////////////
  /////////////////////////////////////////////////////////////////
  std::cout << "------------ NR analysis ------------" << std::endl;
  std::cout << " iter "
            << " first norm "
            << "      res " << std::endl;


  if(thisTime == 1)
  {
        std::cout << "   " << iterstep << "      " << firstNorm << "        "
              << residual << std::endl;
    sol.setZero();
    // Solvers solv(fem, globalK, force.fres, sol);
    // solv.solve(solvName, log);

    solver_.factorize(globalK); // 数値分解のみ（超高速・メモリ確保なし）
    if(solver_.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        exit(1);
    }
    sol = solver_.solve(force.fres);

    // add corrector
    value.val += sol;
    du += sol;

    // make result to node
    divideValIntoNode(fem.dofnp, value.val, du, thisTime, totalTime);
    divideReactionIntoNode(fem.dofnp, reaction, node);

    iterstep++;

    // assembling
    assemblingK(log);

    // cal residual force
    force.fres = force.fext - force.fint;

    /// convergence check
    norm = sqrt(force.fres.dot(force.fres));

    residual = norm / firstNorm;

    std::cout << "   " << iterstep << "      " << firstNorm << "        "
              << residual << std::endl;
  }
  else if(thisTime > 1)
  {
    // set dirichlet boundary condition
    for(auto &n : node)
    {
      if(n.dirich != nullptr)
        for(int j = 0; j < 3; j++)
        {
          if(n.dirich->flag[j] == 1)
          {
            n.val[j] = thisTime * n.dirich->val[j] / totalTime;
            n.du[j] = n.dirich->val[j] / totalTime;
          }
          else
          {
            n.du[j] = n.du_old[j];
            n.val[j] += n.du[j];
          }
        }
      else
        for(int j = 0; j < 3; j++)
        {
          n.du[j] = n.du_old[j];
          n.val[j] += n.du[j];
        }
    }
    integrateValFromNode(fem.dofnp, value.val, du);

    iterstep++;

    // assembling
    assemblingK(log);

    // cal residual force
    force.fres = force.fext - force.fint;

    /// convergence check
    norm = sqrt(force.fres.dot(force.fres));

    residual = norm / firstNorm;

    std::cout << "   " << iterstep << "      " << firstNorm << "        "
              << residual << std::endl;
  }

  double min_residual = residual;
  Eigen::VectorXd best_val = value.val;
  Eigen::VectorXd best_du = du;
  const double emergency_tol = 1e-3;
  const double divergence_limit = 1.0;

  while(residual > tolNR)
  {
    iterstep++;

    sol.setZero();
    // Solvers solv(fem, globalK, force.fres, sol);
    // solv.solve(solvName, log);

    solver_.factorize(globalK); 
    
    if(solver_.info() != Eigen::Success) { break; }
    sol = solver_.solve(force.fres);
    
    // add corrector
    value.val += sol;
    du += sol;

    // make result to node
    divideValIntoNode(fem.dofnp, value.val, du, thisTime, totalTime);
    divideReactionIntoNode(fem.dofnp, reaction, node);

    // assembling
    assemblingK(log);

    // cal residual force
    force.fres = force.fext - force.fint;

    /// convergence check
    norm = sqrt(force.fres.dot(force.fres));

    residual = norm / firstNorm;

    std::cout << "   " << iterstep << "      " << firstNorm << "        "
              << residual << std::endl;

    if (residual < min_residual)
    {
        min_residual = residual;
        best_val = value.val;
        best_du = du;
        // ※この時点ではまだ node への書き戻しは確定しない
    }

    if(iterstep == maxIter || residual > divergence_limit)
    {
      if (min_residual < emergency_tol)
      {
          // A. 許容範囲内なら、ベストな状態を復元して進む
          std::cout << "----------------------------------------------------" << std::endl;
          std::cout << " WARNING: NR failed to converge fully." << std::endl;
          std::cout << "          Restoring best intermediate state." << std::endl;
          std::cout << "          Best Residual = " << min_residual << " (Tolerance: " << emergency_tol << ")" << std::endl;
          std::cout << "----------------------------------------------------" << std::endl;

          // ベストな状態を変数に戻す
          value.val = best_val;
          du = best_du;

          // 重要: 復元した値で node の状態も整合させる
          divideValIntoNode(fem.dofnp, value.val, du, thisTime, totalTime);
          
          // ループを強制的に抜ける
          break;
      }
      else
      {
          // B. 許容範囲を超えていたらエラー終了
          std::cerr << "ERROR: Iteration reached maxIter (" << maxIter << ")." << std::endl;
          std::cerr << "       Best residual found (" << min_residual << ") is still > " << emergency_tol << std::endl;
          exit(1); // または return -1;
      }
      // std::cerr << "itaration step reach " << maxIter
      //           << " in newton_raphson.hpp" << std::endl;
      // return iterstep;
    }
  } ///< NR loop end

  // make result to node
  divideValIntoNode(fem.dofnp, value.val, du, thisTime, totalTime);
  divideReactionIntoNode(fem.dofnp, reaction, node);

  return iterstep;
}

template <class E, class N>
int NewtonRaphson<E, N>::solve_cyclic(int thisTime, int totalTime, double tolNR,
                                      int maxIter, std::string solvName,
                                      int log, double sign)
{

  int iterstep = 0;
  double firstNorm, residual = 1;
  Eigen::VectorXd sol = Eigen::VectorXd::Zero(fem.numeq);
  Eigen::VectorXd du = Eigen::VectorXd::Zero(fem.numeq);

  // set dirichlet boundary condition
  for(auto &n : node)
  {
    if(n.dirich != nullptr)
      for(int j = 0; j < 3; j++)
      {
        if(n.dirich->flag[j] == 1)
        {
          n.du[j] = sign * n.dirich->val[j] / totalTime;
          n.val[j] += n.du[j];
        }
        else
          n.du[j] = 0.0;
      }
    else
      for(int j = 0; j < 3; j++)
      {
        n.du[j] = 0.0;
      }
  }

  /////////////////////////////////////////////////////////////////
  //////////////////////////start NR loop//////////////////////////
  /////////////////////////////////////////////////////////////////
  std::cout << "------------ NR analysis ------------" << std::endl;
  std::cout << " iter "
            << " first norm "
            << "      res " << std::endl;

  while(residual > tolNR)
  {
    iterstep++;

    // assembling
    assemblingK(log);

    // cal residual force
    force.fres = force.fext - force.fint;

    /// convergence check
    double norm = sqrt(force.fres.dot(force.fres));

    if(iterstep == 1)
      firstNorm = norm;

    residual = norm / firstNorm;

    std::cout << "   " << iterstep << "      " << firstNorm << "        "
              << residual << std::endl;

    if(iterstep == maxIter)
    {
      std::cerr << "itaration step reach " << maxIter
                << " in newton_raphson.hpp" << std::endl;
      return iterstep;
    }

    sol.setZero();
    Solvers solv(fem, globalK, force.fres, sol);
    solv.solve(solvName, log);

    // add corrector
    value.val += sol;
    du += sol;

    // make result to node
    divideValIntoNodeCyclic(fem.dofnp, value.val, du, thisTime, totalTime,
                            sign);
    divideReactionIntoNode(fem.dofnp, reaction, node);
  } ///< NR loop end

  // make result to node
  divideValIntoNodeCyclic(fem.dofnp, value.val, du, thisTime, totalTime, sign);
  divideReactionIntoNode(fem.dofnp, reaction, node);

  return iterstep;
}

template <class E, class N> void NewtonRaphson<E, N>::initialization()
{
  /// each nodal values
  for(auto &n : node)
    for(int j = 0; j < 3; j++)
    {
      n.val[j] = 0.0;
      n.du[j] = 0.0;
      n.du_old[j] = 0.0;
      n.reaction[j] = 0.0;
    }

  /// ValueNR class
  value.val.setZero();
  is_solver_analyzed_ = false;
}

template <class E, class N>
double NewtonRaphson<E, N>::calculateTotalReaction(int direction, double target_x, double target_y, double target_z)
{
  double total_reaction = 0.0;
  double tol = 1e-6; // 浮動小数点比較の許容誤差

  for(int i = 0; i < fem.numnp; i++)
  {
    // 1. X座標のチェック（0以上が指定されている場合のみ、一致しなければスルー）
    if(target_x >= 0.0 && std::abs(node[i].x[0] - target_x) > tol) continue;

    // 2. Y座標のチェック（0以上が指定されている場合のみ、一致しなければスルー）
    if(target_y >= 0.0 && std::abs(node[i].x[1] - target_y) > tol) continue;

    // 3. Z座標のチェック（0以上が指定されている場合のみ、一致しなければスルー）
    if(target_z >= 0.0 && std::abs(node[i].x[2] - target_z) > tol) continue;

    // すべての条件をクリアした（指定された支点・固定端にある）節点の反力だけを足す
    total_reaction += node[i].reaction[direction];
  }
  return total_reaction;
}

template <class E, class N> void NewtonRaphson<E, N>::assemblingK(int log)
{
  if(log == 1)
    std::cout << "assembling...";
  fflush(stdout);

  reaction.setZero();
  force.fint.setZero();
  globalK.setZero();
  Eigen::SparseMatrix<double> fintMat(fem.numeq, fem.numeq);
  Eigen::SparseMatrix<double> reactMat(fem.neq, fem.neq);
  tripletsK_.clear();
  tripletsF_.clear();
  tripletsR_.clear();
  // std::vector<Eigen::Triplet<double>> tripletsK, tripletsF, tripletsR;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    /// private variables
    std::vector<Eigen::Triplet<double>> tripK_pri, tripF_pri, tripR_pri;

    int estimated_per_thread = (int)(tripletsK_.capacity() / omp_get_max_threads());
    if (estimated_per_thread > 0) tripK_pri.reserve(estimated_per_thread);

#ifdef _OPENMP
#pragma omp for nowait
#endif
    for(int nel = 0; nel < fem.nelx; nel++)
      element[nel].makeKeFint(tripK_pri, tripF_pri, tripR_pri, fem, node);

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      /// connect
      tripletsK_.insert(tripletsK_.end(), tripK_pri.begin(), tripK_pri.end());
      tripletsF_.insert(tripletsF_.end(), tripF_pri.begin(), tripF_pri.end());
      tripletsR_.insert(tripletsR_.end(), tripR_pri.begin(), tripR_pri.end());
    }
  }

  globalK.setFromTriplets(tripletsK_.begin(), tripletsK_.end());
  fintMat.setFromTriplets(tripletsF_.begin(), tripletsF_.end());
  reactMat.setFromTriplets(tripletsR_.begin(), tripletsR_.end());

  reaction = reactMat.diagonal();
  force.fint = fintMat.diagonal();
  globalK.makeCompressed();

  if(log == 1)
    std::cout << "done" << std::endl;
  fflush(stdout);
}

template <class E, class N>
void NewtonRaphson<E, N>::divideValIntoNode(int dofnp, Eigen::VectorXd &val,
                                            Eigen::VectorXd &du, int thisStep,
                                            int totalStep)
{
  int dofnptmp = std::min(3, dofnp);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < (int)node.size(); i++)
  {
    for(int j = 0; j < dofnptmp; j++)
    {
      if(node[i].dof[j] >= fem.numeq)
      {
        node[i].val[j] = thisStep * node[i].dirich->val[j] / totalStep;
        node[i].du[j] = node[i].dirich->val[j] / totalStep;
      }
      else
      {
        node[i].val[j] = val.coeff(node[i].dof[j]);
        node[i].du[j] = du.coeff(node[i].dof[j]);
      }
    }
  }
}

template <class E, class N>
void NewtonRaphson<E, N>::integrateValFromNode(int dofnp, Eigen::VectorXd &val,
                                               Eigen::VectorXd &du)
{
  int dofnptmp = std::min(3, dofnp);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < (int)node.size(); i++)
  {
    for(int j = 0; j < dofnptmp; j++)
    {
      if(node[i].dof[j] < fem.numeq)
      {
        val.coeffRef(node[i].dof[j]) = node[i].val[j];
        du.coeffRef(node[i].dof[j]) = node[i].du[j];
      }
    }
  }
}

template <class E, class N>
void NewtonRaphson<E, N>::divideValIntoNodeCyclic(int dofnp,
                                                  Eigen::VectorXd &val,
                                                  Eigen::VectorXd &du,
                                                  int thisStep, int totalStep,
                                                  double sign)
{
  int dofnptmp = std::min(3, dofnp);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < (int)node.size(); i++)
  {
    for(int j = 0; j < dofnptmp; j++)
    {
      if(node[i].dof[j] < fem.numeq)
      {
        node[i].val[j] = val.coeff(node[i].dof[j]);
        node[i].du[j] = du.coeff(node[i].dof[j]);
      }
    }
  }
}

} // namespace icarat
