///  @file  arc_length.hpp
///  @author  Daiki Watanabe
///  @date  July 2, 2024.
/// Copyright (c) 2024 SAODLab-Nagoya. All Rights reserved.

#pragma once
#include "analysis_utility.hpp"
#include "solver/solvers.hpp"
#include <problem/base.hpp>

namespace icarat
{
/// FEM procedure for linear elasticity
template <class E, class N> class ArcLength
{
public:
  ArcLength(FEM &fem, std::vector<E> &element, std::vector<N> &node,
            Force &force, Value &value);

  ///  main function (Incremental load or displacement method)
  int solve(double tolNR, int iterNR, std::string solvName, int isFirst,
            int log);

  /// initializaiton
  void initialization();

  /// getters
  const Eigen::SparseMatrix<double> &K() const { return globalK; }
  const Eigen::VectorXd &R() const { return reaction; }

  double lambda_;
  double lambda1_;
  double lambda2_;

protected:
  void divideValIntoNode(int dofnp, Eigen::VectorXd &val, Eigen::VectorXd &du);

  /// assembling stiffness matrix
  void assemblingK(bool log);

  double getArcLength(Eigen::VectorXd &incre, Eigen::VectorXd &incre_lin);

  FEM &fem;
  std::vector<E> &element;
  std::vector<N> &node;
  Force &force;
  Value &value;
  Eigen::SparseMatrix<double> globalK; ///< global stiffness matrix
  Eigen::VectorXd reaction;            ///< reaction force
  Eigen::VectorXd disp1;               ///< disp increment in a previous step
  Eigen::VectorXd disp2;               ///< disp increment in a previous step
  int nitar;                           ///< iteration in previous step
  double ll_;                          ///< length in previous step
  double l_max;                        ///< maximum arc length
  double length_;                      ///< current arc length
  double factor = 8.0, pp = 5.0, dlambda0 = 1.0e-3; ///< arc length parameters
};

template <class E, class N>
ArcLength<E, N>::ArcLength(FEM &fem, std::vector<E> &element,
                           std::vector<N> &node, Force &force, Value &value)
    : fem(fem), element(element), node(node), force(force), value(value),
      globalK(fem.numeq, fem.numeq), reaction(fem.neq), disp1(fem.numeq),
      disp2(fem.numeq), nitar(0), ll_(0.0), l_max(0.0), length_(0.0),
      lambda_(0.0), lambda1_(0.0), lambda2_(0.0)
{
}

template <class E, class N>
int ArcLength<E, N>::solve(double tolNR, int iterNR, std::string solvName,
                           int isFirst, int log)
{
  lambda1_ = 0.0;
  Eigen::VectorXd incre_lin = Eigen::VectorXd::Zero(fem.numeq);
  Eigen::VectorXd incre = Eigen::VectorXd::Zero(fem.numeq);
  disp1.setZero();
  disp2.setZero();
  Solvers solver0(globalK, force.fext, incre_lin);
  Solvers solver1(globalK, force.fres, incre);

  force.fext = force.fext_org;

  assemblingK(log);

  // cal incre_lin
  solver0.solve(solvName, log);

  double aa = incre_lin.dot(incre_lin);

  if(isFirst == 1)
  {
    ll_ = sqrt(pow(dlambda0, 2.0) * aa);
    l_max = pp * ll_;
  }
  else
  {
    ll_ = length_ * factor / (double)nitar;
    if(ll_ > l_max)
      ll_ = l_max;
  }

  double sign = incre.dot(incre_lin);
  if(sign >= 0.0)
    sign = 1.0;
  else
    sign = -1.0;

  if(aa <= 0.0)
  {
    std::cerr << "error in arc_length.hpp" << std::endl;
    exit(1);
  }

  double dlambda = sign * fabs(ll_ / sqrt(aa));
  lambda1_ += dlambda;
  lambda2_ = lambda_ + lambda1_;

  disp1 = dlambda * incre_lin;
  disp2 = value.val + disp1;

  // make result to node
  divideValIntoNode(fem.dofnp, disp2, disp1);

  double fnorm = fabs(lambda2_) * sqrt(force.fext.dot(force.fext));

  /////////////////////////////////////////////////////////////////
  //////////////////////////start NR loop//////////////////////////
  /////////////////////////////////////////////////////////////////
  std::cout << "------------ NR analysis ------------" << std::endl;
  std::cout << " iter " << " first res " << "      res " << std::endl;

  int iterstep = 0;
  double residual = 1.0;
  while(residual > tolNR)
  {
    iterstep++;

    // assembling
    assemblingK(log);

    // cal residual force
    force.fres = lambda2_ * force.fext - force.fint;

    solver0.solve(solvName, log); // incle_lin
    solver1.solve(solvName, log); // incre

    dlambda = getArcLength(incre, incre_lin);

    // update length parameters
    lambda1_ = lambda1_ + dlambda;
    lambda2_ = lambda_ + lambda1_;

    incre += dlambda * incre_lin; // fix incre

    disp1 += incre;
    disp2 = value.val + disp1;

    // make result to node
    divideValIntoNode(fem.dofnp, disp2, disp1);
    divideReactionIntoNode(fem.dofnp, reaction, node);

    /// convergence check
    double resnorm = sqrt(force.fres.dot(force.fres));
    residual = resnorm / fnorm;

    std::cout << "   " << iterstep << "     " << fnorm << "       " << residual
              << std::endl;

    if(iterstep == iterNR)
    {
      std::cerr << "itaration step > " << iterstep << " in arc_length.hpp"
                << std::endl;
      return iterstep;
    }
  } ///< NR loop end

  // store value
  length_ = ll_;
  nitar = iterstep;
  lambda_ = lambda2_;

  value.val = disp2;
  // make result to node
  divideValIntoNode(fem.dofnp, disp2, disp1);
  divideReactionIntoNode(fem.dofnp, reaction, node);

  return iterstep;
}

template <class E, class N> void ArcLength<E, N>::initialization()
{
  /// each nodal values
  for(auto &n : node)
    for(int j = 0; j < 3; j++)
    {
      n.val[j] = 0.0;
      n.du[j] = 0.0;
    }

  /// ValueNR class
  value.val.setZero();
}

template <class E, class N>
void ArcLength<E, N>::divideValIntoNode(int dofnp, Eigen::VectorXd &val,
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
      if(node[i].dof[j] >= fem.numeq)
      {
        node[i].val[j] = fabs(lambda2_) * node[i].dirich->val[j];
      }
      else
      {
        node[i].val[j] = val.coeff(node[i].dof[j]);
      }
    }
  }
}

template <class E, class N> void ArcLength<E, N>::assemblingK(bool log)
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
  reactMat.setFromTriplets(tripletsR.begin(), tripletsR.end());

  reaction = reactMat.diagonal();
  force.fint = fintMat.diagonal();
  globalK.makeCompressed();

  if(log == 1)
    std::cout << "done" << std::endl;
  fflush(stdout);
}

template <class E, class N>
double ArcLength<E, N>::getArcLength(Eigen::VectorXd &incre,
                                     Eigen::VectorXd &incre_lin)
{
  double dlambda;
  double tol = 1.0e-8;

  double aa = 0.0, bb = 0.0, cc = 0.0, dd = 0.0;
  for(int i = 0; i < this->fem.numeq; i++)
  {
    aa += incre_lin.coeff(i) * incre_lin.coeff(i);
    bb += 2.0 * (disp1.coeff(i) + incre.coeff(i)) * incre_lin.coeff(i);
    cc += std::pow(disp1(i) + incre(i), 2.0) - std::pow(ll_, 2.0);
  }

  double val0 = pow(bb, 2.0) - 4.0 * aa * cc;

  /// cal dlambda
  std::vector<double> dl(2);
  if(val0 > tol)
  {
    if(bb >= 0.0)
    {
      dl[0] = (-bb - sqrt(val0)) / (2.0 * aa);
      dl[1] = (cc / aa) / dl[0];
    }
    else
    {
      dl[0] = (-bb + sqrt(val0)) / (2.0 * aa);
      dl[1] = (cc / aa) / dl[0];
    }

    double val1 = (disp1 + incre + dl[0] * incre_lin).dot(disp1);
    double val2 = (disp1 + incre + dl[1] * incre_lin).dot(disp1);

    if(val1 > val2)
      dlambda = dl[0];
    else
      dlambda = dl[1];
  }
  else if(fabs(val0) <= tol)
  {
    dlambda = -bb / (2.0 * aa);
  }
  else // if(val0 < tol)
  {
    std::cerr << "arc length parameter is negative in newton_raphson.hpp"
              << std::endl;
    exit(1);
  }

  return dlambda;
}
} // namespace icarat
