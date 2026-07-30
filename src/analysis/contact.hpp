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
///   fem_m procedure for linear elasticity
template <class E, class N> class ContactAnalysis
{
public:
  ContactAnalysis(FEM &fem_m, std::vector<E> &element_m, std::vector<N> &node_m,
                  Force &force_m, Value &value_m, FEM &fem_s,
                  std::vector<E> &element_s, std::vector<N> &node_s,
                  Force &force_s, Value &value_s)
      : fem_m(fem_m), element_m(element_m), node_m(node_m), force_m(force_m),
        value_m(value_m), globalK(fem_m.numeq, fem_m.numeq),
        reaction_m(fem_m.neq)
  {
    std::cout << "Element class name is " << typeid(E).name() << std::endl;
    std::cout << "Node class name is " << typeid(N).name() << std::endl;
  }

  ///  main function (Incremental load or displacement method)
  ///@return Newton-Raphson iteration
  int solve(int thisTime, int totalTime, double tolNR, int maxIter,
            std::string solvName, int log);

  ///   initializaiton
  void initialization();

  /// assembling stiffness matrix & inner force
  void assemblingK(bool log);

  /// @brief ペナルティ法で接触力と接触分の剛性行列を作成する
  /// @param omegan 垂直応力のペナルティ係数
  /// @param omegat 摩擦力のペナルティ係数
  ///@param elxy masterとslaveの現在位置．1列目がslave節点，2,3列目がmaster
  ///@param elxyp 前ステップの~~
  /// @param cfri 摩擦係数.cfri=0.0なら摩擦を考慮しない．
  void cntelm2d(double omegan, double omegat, Eigen::MatrixXd &elxy,
                Eigen::MatrixXd &elxyp, double cfri);

  /// getters
  const Eigen::SparseMatrix<double> &K() const { return globalK; }

private:
  FEM &fem_m, &fem_s;
  std::vector<E> &element_m, &element_s;
  std::vector<N> &node_m, node_s;
  Force &force_m, &force_s;
  Value &value_m, &value_s;
  Eigen::SparseMatrix<double> globalK;    ///< global stiffness matrix
  Eigen::VectorXd reaction_m, reaction_s; ///< reaction force
};

template <class E, class N>
int ContactAnalysis<E, N>::solve(int thisTime, int totalTime, double tolNR,
                                 int maxIter, std::string solvName, int log)
{

  int iterstep = 0;
  double firstNorm, residual = 1;
  Eigen::VectorXd sol = Eigen::VectorXd::Zero(fem_m.numeq);
  Eigen::VectorXd du = Eigen::VectorXd::Zero(fem_m.numeq);

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
    force_m.fres = force_m.fext - force_m.fint;

    /// convergence check
    double norm = sqrt(force_m.fres.dot(force_m.fres));
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
    Solvers solv(globalK, force_m.fres, sol);
    solv.solve(solvName, log);

    // add corrector
    value.val += sol;
    du += sol;

    // make result to node
    divideValIntoNode(fem_m.dofnp, value.val, du, thisTime, totalTime);
    divideReactionIntoNode(fem_m.dofnp, reaction, node);
  } ///< NR loop end

  // make result to node
  divideValIntoNode(fem_m.dofnp, value.val, du, thisTime, totalTime);
  divideReactionIntoNode(fem_m.dofnp, reaction, node);

  return iterstep;
}

template <class E, class N> void ContactAnalysis<E, N>::initialization()
{
  /// each nodal values
  for(auto &n : node)
    for(int j = 0; j < 3; j++)
      n.val[j] = 0.0;

  /// Value class
  value.val.setZero();
}

template <class E, class N> void ContactAnalysis<E, N>::assemblingK(bool log)
{
  if(log == 1)
    std::cout << "assembling...";
  fflush(stdout);

  force.fint.setZero();
  globalK.setZero();
  reaction.setZero();

  Eigen::SparseMatrix<double> fintMat(fem_m.numeq, fem_m.numeq);
  Eigen::SparseMatrix<double> globalR(fem_m.neq, fem_m.neq);
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
    for(int nel = 0; nel < fem_m.nelx; nel++)
      element[nel].makeKeFint(tripK_pri, tripF_pri, tripR_pri, fem_m, node);

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

template <class E, class N>
void ContactAnalysis<E, N>::cntelm2d(double omegan, double omegat,
                                     Eigen::MatrixXd &elxy,
                                     Eigen::MatrixXd &elxyp, double cfri)
{
  double para_seg = 0.05; // NTS法のセグメントに関するパラメータ
  double eps = 1.0e-6;

  Eigen::VectorXd xt = elxy.col(2) - elxy.col(1); // 接線ベクトル
  Eigen::VectorXd xtp = elxyp.col(2) - elxyp.col(1); // 前ステップの接線ベクトル
  double xlen = xt.norm(), xlenp = xtp.norm();

  if(xlen < eps) // masterの2節点間が小さすぎる場合は中断
    return;

  xt = xt / xlen; // 単位接線ベクトル
  Eigen::VectorXd xn(2);
  xn << -xt(1), xt(0); // 単位法線ベクトル

  // 法線隙間関数
  Eigen::VectorXd gap = elxy.col(0) - elxy.col(1);
  Eigen::VectorXd gapp = elxyp.col(0) - elxyp.col(1);
  double gapn = gap.dot(xn);

  // 接触していない or 貫入しすぎている場合は計算しない
  if(gapn = > 0.0 || gapn <= -xlen)
    return;

  // 接触点での自然座標
  double alpha = gap.dot(xt) / xlen;
  double alpha0 = (gapp.dot(xtp)) / xlenp;

  // 自然座標の情報から，slaveがmasterのセグメント外なら計算しない
  if(alpha > 1.0 + para_seg || alpha < -para_seg)
    return;

  // 接触力(ペナルティ)の大きさ
  double xlambn = -omegan * gapn;

  int lfric = 1; //=1なら摩擦考慮する
  if(cfri == 0.0)
    lfric = 0;

  int lslide = 0; //=1なら摩擦力発生
  double gapt = 0.0;
  double xlambt = 0.0; // 摩擦力の大きさ
  if(lfric)
  {
    gapt = (alpha - alpha0) * xlenp;
    xlambt = -omegat * gapt;
    double frtol = xlambn * cfri; // 摩擦の滑りクライテリア
    if(abs(xlambt) > frtol)       // 滑る場合
    {
      lslide = 1;
      xlambt = -frtol * sign(one, gapt);
    }
  }

  Eigen::MatrixXd nn(2, 3), tt(2, 3), pp(2, 3), qq(2, 3);
  nn.col(0) = xn;
  nn.col(1) = -(1.0 - alpha) * xn;
  nn.col(2) = -alpha * xn;

  tt.col(0) = xt;
  tt.col(1) = -(one - alpha) * xt;
  tt.col(2) = -alpha * xt;

  pp.setZero();
  pp.col(1) = -xn;
  pp.col(2) = xn;

  qq.setZero();
  qq.col(1) = -xt;
  qq.col(2) = xt;

  Eigen::MatrixXd cn = nn - gapn * qq / xlen;
  Eigen::MatrixXd ct = tt + gapn * pp / xlen;

  // 接触力の要素節点値ベクトル
  Eigen::VectorXd force = xlambn * cn + xlambt * ct;

  // 垂直応力分の剛性行列
  Eigen::MatrixXd stiff = omegan * cn * cn.transpose();

  // 摩擦分の剛性行列
  if(lfric)
  {
    double tmp1 = -cfri * omegan * sign(1.0, gapt);
    double tmp2 = -xlambt / xlen;
    if(lslide) // 摩擦力発生
      stiff + = tmp1 * (ct * cn.transpose()) +
                tmp2 * (cn * pp.transpose() + pp * cn.transpose() -
                        ct * qq.transpose() - qq * ct.transpose());
    else // 摩擦力発生しない
      stiff += omegat * (ct * ct.transpose()) +
               tmp2 * (cn * pp.transpose() + pp * cn.transpose() -
                       ct * qq.transpose() - qq * ct.transpose());
  }
}

} // namespace icarat