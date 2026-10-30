///  @file  newton_raphson.hpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php
#pragma once
#include <analysis/analysis_utility.hpp>
#include <analysis/solver/solvers.hpp>
#include <cmath>
#include <problem/base.hpp>

namespace icarat
{
/// FEM procedure for linear elasticity
template <class E, class N> class NewtonRaphson
{
public:
  NewtonRaphson(FEM &fem, std::vector<E> &element, std::vector<N> &node,
                Force &force, Value &value, int numstep,
                std::vector<double> tauFiltered);

  std::vector<double> tau();
  ///  main function (Incremental load or displacement method)
  ///@return Newton-Raphson iteration
  int solve(int thisTime, int totalTime, double tolNR, int maxIter,
            std::string solvName, int log, std::vector<double> tauFiltered,
            std::vector<double> dtauFiltered_dtau);

  /// main function (arc length method)
  void solveArcLength(int thisTime, int totalTime, std::string solvName,
                      int log);

  /// initializaiton
  void initialization();

  /// getters
  const Eigen::SparseMatrix<double> &K() const { return globalK; }

  //----------------------- for reaction -------------------------//
  /*const Eigen::VectorXd &reaction_() const { return reaction; }
  const Eigen::SparseMatrix<double> &KFull() const { return globalKFull; }*/
  //--------------------------------------------------------------//
  //----------------------- for S:E integration ------------------//
  const Eigen::VectorXd &se() const { return productse_; }
  //--------------------------------------------------------------//

private:
  //----------------------- for reaction -------------------------//
  /*void assemblingK(int log, int thisTime);*/
  //--------------------------------------------------------------//
  //----------------------- for S:E integration ------------------//
  void assemblingK(int log, int thisTime, int totalTime,
                   const std::vector<double> &tauFiltered,
                   const std::vector<double> &dtauFiltered_dtau);
  //--------------------------------------------------------------//

  void divideValIntoNode(int dofnp, Eigen::VectorXd &val, Eigen::VectorXd &du,
                         int thisTime, int totalTime);

  void integrateValFromNode(int dofnp, Eigen::VectorXd &val,
                            Eigen::VectorXd &du);

  void divideValIntoNodeCyclic(int dofnp, Eigen::VectorXd &val,
                               Eigen::VectorXd &du, int thisStep, int totalStep,
                               double sign);

  FEM &fem;
  std::vector<E> &element;
  std::vector<N> &node;
  std::vector<double> tau_all;
  Force &force;
  Value &value;
  Eigen::SparseMatrix<double> globalK; ///< global stiffness matrix
  Eigen::VectorXd reaction;            ///< reaction force

  //----------------------- for reaction -------------------------//
  /*Eigen::SparseMatrix<double> globalKFull;
  std::vector<Eigen::Triplet<double>> tripletsKfull;*/
  //--------------------------------------------------------------//

  //----------------------- for S:E integration ------------------//
  Eigen::VectorXd productse_;
  //--------------------------------------------------------------//
};

template <class E, class N>
NewtonRaphson<E, N>::NewtonRaphson(FEM &fem, std::vector<E> &element,
                                   std::vector<N> &node, Force &force,
                                   Value &value, int numstep,
                                   std::vector<double> tauFiltered)
    : fem(fem), element(element), node(node), force(force), value(value),
      globalK(fem.numeq, fem.numeq), reaction(fem.neq),
      //----------------------- for reaction -------------------------//
      // globalKFull(fem.neq, fem.neq)
      //--------------------------------------------------------------//
      //----------------------- for S:E integration ------------------//
      productse_(fem.nelx)
//--------------------------------------------------------------//
{
  std::cout << "Element class name is " << typeid(E).name() << std::endl;
  std::cout << "Node class name is " << typeid(N).name() << std::endl;
}

template <class E, class N> std::vector<double> NewtonRaphson<E, N>::tau()
{
  std::vector<double> tau_all(fem.nelx, 0.0);

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    double test = element[nel].computeTau(fem, node);
    tau_all[nel] = test;
  }

  return tau_all;
}

template <class E, class N>
int NewtonRaphson<E, N>::solve(int thisTime, int totalTime, double tolNR,
                               int maxIter, std::string solvName, int log,
                               std::vector<double> tauFiltered,
                               std::vector<double> dtauFiltered_dtau)
{
  int iterstep = 0;
  double firstNorm, residual = 1;
  Eigen::VectorXd sol = Eigen::VectorXd::Zero(fem.numeq);
  Eigen::VectorXd du = Eigen::VectorXd::Zero(fem.numeq);

  // set dirichlet boundary condition
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
  //----------------------- for reaction -------------------------//
  // assemblingK(log, thisTime);
  //--------------------------------------------------------------//
  //----------------------- for S:E integration ------------------//
  assemblingK(log, thisTime, totalTime, tauFiltered, dtauFiltered_dtau);
  //--------------------------------------------------------------//

  // cal residual force
  force.fres = force.fext - force.fint;

  // cal first norm
  firstNorm = sqrt(force.fres.dot(force.fres));

  /////////////////////////////////////////////////////////////////
  //////////////////////////start NR loop//////////////////////////
  /////////////////////////////////////////////////////////////////
  std::cout << "------------ NR analysis ------------" << std::endl;
  std::cout << " iter " << " first norm " << "      res " << std::endl;

  if(thisTime == 1)
  {
    std::cout << "   " << iterstep << "      " << firstNorm << "        "
              << residual << std::endl;
    sol.setZero();
    Solvers solv(globalK, force.fres, sol);
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
    //----------------------- for reaction -------------------------//
    // assemblingK(log, thisTime);
    //--------------------------------------------------------------//
    //----------------------- for S:E integration ------------------//
    assemblingK(log, thisTime, totalTime, tauFiltered, dtauFiltered_dtau);
    //--------------------------------------------------------------//

    // cal residual force
    force.fres = force.fext - force.fint;

    /// convergence check
    double norm = sqrt(force.fres.dot(force.fres));
    if(iterstep == 1)
    {
      firstNorm = norm;
    }

    residual = norm / firstNorm;

    std::cout << "   " << iterstep << "      " << firstNorm << "        "
              << residual << std::endl;

    if(iterstep == maxIter)
    {
      std::cerr << "itaration step reach " << maxIter
                << " in newton_raphson.hpp" << std::endl;
      return iterstep;
    }

    /// solving the linear equation
    sol.setZero();
    Solvers solv(globalK, force.fres, sol);
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

/*template <class E, class N>
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
  std::cout << " iter " << " first norm " << "      res " << std::endl;

  while(residual > tolNR)
  {
    iterstep++;

    // assembling
    assemblingK(log, thisTime);

    // cal residual force
    // force.fres = force.fext - force.fint;
    force.fres =
        force.fexternal.col(thisTime - 1) - force.fintarray.col(thisTime - 1);

    /// convergence check
    double norm = sqrt(force.fres.dot(force.fres));
    if(iterstep == 1)
    {
      firstNorm = norm;
    }

    residual = norm / firstNorm;

    std::cout << "   " << iterstep << "      " << firstNorm << "        "
              << residual << std::endl;

    if(iterstep == maxIter)
    {
      std::cerr << "itaration step reach " << maxIter
                << " in newton_raphson.hpp" << std::endl;
      return iterstep;
    }

    /// solving the linear equation
    sol.setZero();
    Solvers solv(globalK, force.fres, sol);
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
}*/

template <class E, class N> void NewtonRaphson<E, N>::initialization()
{
  /// each nodal values
  for(auto &n : node)
    for(int j = 0; j < 3; j++)
    {
      n.val[j] = 0.0;
      n.du[j] = 0.0;
      n.reaction[j] = 0.0;
      n.reaeach.setZero();
      n.valeach.setZero();
    }

  /// ValueNR class
  value.val.setZero();
}

//----------------------- for reaction -------------------------//
/*template <class E, class N>
void NewtonRaphson<E, N>::assemblingK(int log, int thisTime)
{
  if(log == 1)
    std::cout << "assembling...";
  fflush(stdout);

  Eigen::SparseMatrix<double> fintMat(fem.numeq, fem.numeq);
  Eigen::SparseMatrix<double> reactMat(fem.neq, fem.neq);
  std::vector<Eigen::Triplet<double>> tripletsK, tripletsF, tripletsR;
  reaction.setZero();
  force.fint.setZero();
  globalK.setZero();
  globalKFull.setZero();
  fintMat.setZero();
  reactMat.setZero();
  tripletsK.clear();
  tripletsKfull.clear();
  tripletsF.clear();
  tripletsR.clear();

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    /// private variables
    std::vector<Eigen::Triplet<double>> tripK_pri, tripK_prifull, tripF_pri,
        tripR_pri;

#ifdef _OPENMP
#pragma omp for nowait
#endif
    for(int nel = 0; nel < fem.nelx; nel++)
    {
      element[nel].makeKeFint(tripK_pri, tripK_prifull, tripF_pri, tripR_pri,
                              fem, node, thisTime);
    }

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      /// connect
      tripletsK.insert(tripletsK.end(), tripK_pri.begin(), tripK_pri.end());
      tripletsKfull.insert(tripletsKfull.end(), tripK_prifull.begin(),
                           tripK_prifull.end());
      tripletsF.insert(tripletsF.end(), tripF_pri.begin(), tripF_pri.end());
      tripletsR.insert(tripletsR.end(), tripR_pri.begin(), tripR_pri.end());
    }
  }

  globalK.setFromTriplets(tripletsK.begin(), tripletsK.end());
  globalKFull.setFromTriplets(tripletsKfull.begin(), tripletsKfull.end());
  fintMat.setFromTriplets(tripletsF.begin(), tripletsF.end());
  reactMat.setFromTriplets(tripletsR.begin(), tripletsR.end());
  reaction = reactMat.diagonal();
  force.fint = fintMat.diagonal();

  globalK.makeCompressed();
  globalKFull.makeCompressed();

  if(log == 1)
    std::cout << "done" << std::endl;
  fflush(stdout);
}*/
//--------------------------------------------------------------//

//----------------------- for S:E integration ------------------//
template <class E, class N>
void NewtonRaphson<E, N>::assemblingK(
    int log, int thisTime, int totalTime,
    const std::vector<double> &tauFiltered,
    const std::vector<double> &dtauFiltered_dtau)
{
  if(log == 1)
    std::cout << "assembling...";
  fflush(stdout);

  reaction.setZero();
  force.fint.setZero();
  globalK.setZero();
  Eigen::SparseMatrix<double> fintMat(fem.numeq, fem.numeq);
  Eigen::SparseMatrix<double> reactMat(fem.neq, fem.neq);
  std::vector<Eigen::Triplet<double>> tripletsK, tripletsF, tripletsR;
  double product_all = 0.0;

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
    {
      double product_energy = 0.0;
      double tauFilteredElement = tauFiltered[nel];
      double dtauFiltered_dtauElement = dtauFiltered_dtau[nel];
      element[nel].makeKeFint(tripK_pri, tripF_pri, tripR_pri, fem, node,
                              thisTime, product_energy, tauFilteredElement,
                              dtauFiltered_dtauElement);
      productse_[nel] = product_energy;
    }

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
  reactMat.setFromTriplets(tripletsR.begin(), tripletsR.end());
  reaction = reactMat.diagonal();
  force.fint = fintMat.diagonal();

  globalK.makeCompressed();

  if(log == 1)
    std::cout << "done" << std::endl;
  fflush(stdout);
}
//--------------------------------------------------------------//

template <class E, class N>
void NewtonRaphson<E, N>::divideValIntoNode(int dofnp, Eigen::VectorXd &val,
                                            Eigen::VectorXd &du, int thisTime,
                                            int totalTime)
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
        node[i].val[j] = thisTime * node[i].dirich->val[j] / totalTime;
        node[i].du[j] = node[i].dirich->val[j] / totalTime;
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

} // namespace icarat
