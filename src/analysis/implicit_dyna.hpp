///  @file  implicit_dyna.hpp
///  @author  Daiki Watanabe
///  @date  May 14, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php
#pragma once
#include "analysis/analysis_utility.hpp"
#include "solver/solvers.hpp"
#include <problem/base.hpp>

namespace icarat
{
/// FEM analysis class with newmark beta method
template <class E, class N> class ImplicitDyna
{
public:
  ImplicitDyna(FEM &fem, std::vector<E> &element, std::vector<N> &node,
               Force &force, ValueDyna &value);

  void initialization();

  /// main function for Dynamic problem Ma+Kx=F
  void solve(int thisStep, int totalStep, double deltaT, double beta,
             double gamma, std::string solvName, bool log);

  /// main function for Rayleigh damping
  void calRayleigh(int thisStep, int totalStep, double aa, double bb,
                   double deltaT, double beta, double gamma,
                   std::string solvName, bool log);

  /// getters
  const Eigen::SparseMatrix<double> &K() const { return globalK; }
  const Eigen::SparseMatrix<double> &M() const { return globalM; }

  ///  assembling stiffness matrix & mass matrix
  void assemblingMK(std::string solName, bool log);

  void divideValIntoNode(int dofnp, int numeq, Eigen::VectorXd &value,
                         std::vector<N> &node);

protected:
  FEM &fem;
  std::vector<E> &element;
  std::vector<N> &node;
  Force &force;
  ValueDyna &value;
  Eigen::SparseMatrix<double> globalK; ///< global stiffness matrix
  Eigen::SparseMatrix<double> globalM; ///< global mass matrix
  Eigen::VectorXd reaction;            ///< reaction force
};

template <class E, class N>
ImplicitDyna<E, N>::ImplicitDyna(FEM &fem, std::vector<E> &element,
                                 std::vector<N> &node, Force &force,
                                 ValueDyna &value)
    : fem(fem), element(element), node(node), force(force), value(value),
      globalK(fem.numeq, fem.numeq), globalM(fem.numeq, fem.numeq),
      reaction(fem.neq)
{
  std::cout << "Element class name is " << typeid(E).name() << std::endl;
  std::cout << "Node class name is " << typeid(N).name() << std::endl;
}

template <class E, class N>
void ImplicitDyna<E, N>::solve(int thisStep, int totalStep, double deltaT,
                               double beta, double gamma, std::string solvName,
                               bool log)
{
  calRayleigh(thisStep, totalStep, 0.0, 0.0, deltaT, beta, gamma, solvName,
              log);
}

template <class E, class N>
void ImplicitDyna<E, N>::calRayleigh(int thisStep, int totalStep, double aa,
                                     double bb, double deltaT, double beta,
                                     double gamma, std::string solvName,
                                     bool log)
{

  int iterstep = 0;
  double firstNorm, residual = 1;

  value.du_incre.setZero();

  // set coeff for new mark beta
  double A1 = 1.0 / (beta * pow(deltaT, 2.0));
  double A2 = gamma / (beta * deltaT);
  double A3 = (gamma / beta) - 1.0;
  double A4 = 1.0 / (beta * deltaT);
  double A5 = ((0.5 * gamma) / beta - 1.0) * deltaT;
  double A6 = (0.5 / beta) - 1.0;

  // val at time t
  Eigen::VectorXd val_old = value.val;

  /////////////////////////////////////////////////////////////////
  //////////////////////////start NR loop//////////////////////////
  /////////////////////////////////////////////////////////////////
  std::cout << "------------ NR analysis ------------" << std::endl;
  std::cout << " iter "
            << " first res "
            << "      res " << std::endl;

  Eigen::VectorXd du = Eigen::VectorXd::Zero(fem.numeq);

  while(residual > 1.0e-7)
  {
    iterstep++;

    // assembling globalK, globalM and Fint
    assemblingMK("implicit", log);

    // make rayleigh matrix
    Eigen::SparseMatrix<double> globalMK(fem.numeq, fem.numeq);
    globalMK = A1 * globalM + globalK + A2 * (aa * globalK + bb * globalM);
    globalMK.makeCompressed();

    divideReactionIntoNode(fem.dofnp, reaction, node);

    Eigen::VectorXd j1 = A6 * value.accel + A4 * value.velocity;
    Eigen::VectorXd j2 = A5 * value.accel + A3 * value.velocity;
    Eigen::VectorXd Fm = globalM * j1;
    Eigen::VectorXd Fd = (aa * globalK + bb * globalM) * j2;
    Eigen::VectorXd F1 = (globalMK - globalK) * val_old;
    Eigen::VectorXd F2 = (globalMK - globalK) * value.val;

    // cal residual force
    force.fres = force.fext - force.fint + Fm + Fd + F1 - F2;

    // convergence check
    double norm = sqrt(force.fres.dot(force.fres));
    if(iterstep == 1)
      firstNorm = norm;
    residual = norm / firstNorm;

    std::cout << "   " << iterstep << "   " << firstNorm << "   " << residual
              << std::endl;

    Solvers solv(globalMK, force.fres, du);
    solv.solve(solvName, log);

    value.val += du;
    value.du_incre += du;
    this->divideValIntoNode(fem.dofnp, fem.numeq, value.val, node);

    if(iterstep == 100)
    {
      std::cerr << "itaration step > max itaration step" << std::endl;
      exit(1);
    }
  }

  value.accel = A1 * value.du_incre - A4 * value.velocity - A6 * value.accel;
  value.velocity = A2 * value.du_incre - A3 * value.velocity - A5 * value.accel;

  divideValIntoNode(fem.dofnp, fem.numeq, value.val, node);
}

template <class E, class N> void ImplicitDyna<E, N>::initialization()
{
  for(auto &n : node)
    for(int j = 0; j < 3; j++)
      n.val[j] = 0.0;

  value.accel.setZero();    // solved value
  value.velocity.setZero(); // solved value
  value.val.setZero();
  value.du_incre.setZero();

  force.fint.setZero();
  force.fres.setZero();
}

template <class E, class N>
void ImplicitDyna<E, N>::assemblingMK(std::string solName, bool log)
{
  if(log == 1)
    std::cout << "assembling...";
  fflush(stdout);

  force.fint.setZero();
  globalK.setZero();
  globalM.setZero();
  Eigen::SparseMatrix<double> fintMat(fem.numeq, fem.numeq);
  Eigen::SparseMatrix<double> globalR(fem.neq, fem.neq);
  std::vector<Eigen::Triplet<double>> tripletsK, tripletsM, tripletsF,
      tripletsR;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // private variables
    std::vector<Eigen::Triplet<double>> tripK_pri, tripM_pri, tripF_pri,
        tripR_pri;

#ifdef _OPENMP
#pragma omp for nowait
#endif
    for(int nel = 0; nel < fem.nelx; nel++)
      element[nel].makeKeMeFint(tripK_pri, tripM_pri, tripF_pri, tripR_pri, fem,
                                node, solName);

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      // connect
      tripletsK.insert(tripletsK.end(), tripK_pri.begin(), tripK_pri.end());
      tripletsM.insert(tripletsM.end(), tripM_pri.begin(), tripM_pri.end());
      tripletsF.insert(tripletsF.end(), tripF_pri.begin(), tripF_pri.end());
      tripletsR.insert(tripletsR.end(), tripR_pri.begin(), tripR_pri.end());
    }
  }

  globalK.setFromTriplets(tripletsK.begin(), tripletsK.end());
  globalK.makeCompressed();
  globalM.setFromTriplets(tripletsM.begin(), tripletsM.end());
  globalM.makeCompressed();
  fintMat.setFromTriplets(tripletsF.begin(), tripletsF.end());
  force.fint = fintMat.diagonal();
  globalR.setFromTriplets(tripletsR.begin(), tripletsR.end());
  reaction = globalR.diagonal();

  if(log)
    std::cout << "done" << std::endl;
  fflush(stdout);
}

template <class E, class N>
void ImplicitDyna<E, N>::divideValIntoNode(int dofnp, int numeq,
                                           Eigen::VectorXd &value,
                                           std::vector<N> &node)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < (int)node.size(); i++)
  {
    for(int j = 0; j < 6; j++)
    {
      if(node[i].dof[j] < numeq)
        node[i].val[j] = value.coeff(node[i].dof[j]);
      else
        continue;

      // Too small value is deleted for paraview
      if(abs(node[i].val[j]) < 1.0e-15)
        node[i].val[j] = 0.0;
    }
  }
}

} // namespace icarat
