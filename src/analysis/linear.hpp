///  @file  linear.hpp
///  @author  Daiki Watanabe
///  @date  January 22, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php
#pragma once
#include "analysis_utility.hpp"
#include "solver/solvers.hpp"
#include <problem/base.hpp>

namespace icarat
{
///   FEM procedure for linear elasticity
template <class E, class N> class LinearAnalysis
{
public:
  LinearAnalysis(FEM &fem, std::vector<E> &element, std::vector<N> &node,
                 Force &force, Value &value);

  ///   main function
  void solve(std::string solvName, bool log);

  ///   initializaiton
  void initialization();

  /// assembling stiffness matrix & inner force
  void assemblingK(bool log);

  /// getters
  const Eigen::SparseMatrix<double> &K() const { return globalK; }

private:
  FEM &fem;
  std::vector<E> &element;
  std::vector<N> &node;
  Force &force;
  Value &value;
  Eigen::SparseMatrix<double> globalK; ///< global stiffness matrix
  Eigen::VectorXd reaction;            ///< reaction force
};

template <class E, class N>
LinearAnalysis<E, N>::LinearAnalysis(FEM &fem, std::vector<E> &element,
                                     std::vector<N> &node, Force &force,
                                     Value &value)
    : fem(fem), element(element), node(node), force(force), value(value),
      globalK(fem.numeq, fem.numeq), reaction(fem.neq)
{
  std::cout << "Element class name is " << typeid(E).name() << std::endl;
  std::cout << "Node class name is " << typeid(N).name() << std::endl;
}

template <class E, class N>
void LinearAnalysis<E, N>::solve(std::string solvName, bool log)
{
  // init
  initialization();

  std::cout << "--------------- linear analysis ---------------" << std::endl;

  // set dirichlet boundary condition
  for(auto &n : node)
  {
    if(n.dirich != nullptr)
      for(int j = 0; j < 6; j++)
      {
        if(n.dirich->flag[j] == 1)
          n.val[j] = n.dirich->val[j];
      }
  }

  // cal external force at 1 step
  force.fext = force.fext_org;

  // assembling
  assemblingK(log);

  // cal residual force
  force.fres = force.fext - force.fint;
  double fnorm = sqrt(force.fres.dot(force.fres));

  // solve KU=F with
  Solvers solv(globalK, force.fres, value.val);
  solv.solve(solvName, log);

  divideValIntoNode(fem.dofnp, fem.numeq, value.val, node);

  // compute reaction force
  assemblingK(log);
  divideReactionIntoNode(fem.dofnp, reaction, node);

  // check convergence of residual force
  if(log)
  {
    force.fres = force.fext - force.fint;
    double norm = sqrt(force.fres.dot(force.fres));
    double residual = norm / fnorm;

    std::cout << "first norm "
              << "      res " << std::endl;
    std::cout << fnorm << "       " << residual << std::endl;
  }
}

template <class E, class N> void LinearAnalysis<E, N>::initialization()
{
  /// each nodal values
  for(auto &n : node)
    for(int j = 0; j < 6; j++)
      n.val[j] = 0.0;

  /// Value class
  value.val.setZero();
}

template <class E, class N> void LinearAnalysis<E, N>::assemblingK(bool log)
{
  if(log == 1)
    std::cout << "assembling...";
  fflush(stdout);

  force.fint.setZero();
  globalK.setZero();
  reaction.setZero();

  Eigen::SparseMatrix<double> fintMat(fem.numeq, fem.numeq);
  Eigen::SparseMatrix<double> globalR(fem.neq, fem.neq);
  std::vector<Eigen::Triplet<double>> tripletsK, tripletsF, tripletsR;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    /// private variables
    std::vector<Eigen::Triplet<double>> tripK_pri, tripF_pri, tripR_pri;

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
      tripletsK.insert(tripletsK.end(), tripK_pri.begin(), tripK_pri.end());
      tripletsF.insert(tripletsF.end(), tripF_pri.begin(), tripF_pri.end());
      tripletsR.insert(tripletsR.end(), tripR_pri.begin(), tripR_pri.end());
    }
  }

  globalK.setFromTriplets(tripletsK.begin(), tripletsK.end());
  fintMat.setFromTriplets(tripletsF.begin(), tripletsF.end());
  globalR.setFromTriplets(tripletsR.begin(), tripletsR.end());

  force.fint = fintMat.diagonal();
  reaction = globalR.diagonal();
  globalK.makeCompressed();

  if(log)
    std::cout << "done" << std::endl;
  fflush(stdout);
}

} // namespace icarat