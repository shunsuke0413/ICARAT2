///  @file  homogenization.hpp
///  @author  Daiki Watanabe
///  @date  December 5, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php
#pragma once
#include "solver/solvers.hpp"
#include <Eigen/Sparse>
#include <problem/elasticity.hpp>
#include <vector>
#ifdef __INTEL_LLVM_COMPILER
#include <Eigen/PardisoSupport>
#endif
namespace icarat
{
///   A class of numerical material tests in homogenization theory
/// @cite 寺田賢二郎ら, 「数値材料試験」,丸善出版(2021).
template <class E, class N> class Homogenization
{
public:
  Homogenization(FEM &fem, std::vector<E> &element, std::vector<N> &node,
                 ValueHomo &value,
                 std::vector<std::vector<std::pair<int, int>>> &pairs,
                 std::array<double, 3> &length);

  void initialization();

  /// solve Numerical material testing: Dirichlet problem
  Eigen::MatrixXd solveNMT(int log);

  /// solve localization analysis
  /// @param [in] barE macro strain
  /// @param [in] barS macro stress
  void solveLocal(Eigen::VectorXd &barE, Eigen::VectorXd &barS, int log);

  /// solve non-linear Numerical material testing
  /// @param [in] MstrainInce macro strain increment
  /// @param [in] forcedConti whether to force to continue(true) or not(false)
  Eigen::MatrixXd solveNLNMT(const Eigen::VectorXd &MstrainIncre, double tolNR,
                             int maxIter, std::string solvName, int log,
                             bool forcedConti);

  /// solve non-linear localization analysis
  void solveNLLocal(Eigen::VectorXd &MstrainIncre, Eigen::VectorXd &Mstress,
                    double tolNR, int maxIter, int log);

  /// getters
  const Eigen::SparseMatrix<double> &K() const { return reducedK; }

private:
  /// assign MPC object & dof to slave nodes
  void setControlPointDOF();

  /// set macrostrain as control point value & voigt expression
  Eigen::VectorXd addControlPointValue(const Eigen::VectorXd &Mstrain,
                                       int direc, bool to_nmtval);

  /// get macro stresses from the nodal load obtained.
  Eigen::VectorXd getMstress(Eigen::VectorXd &react);

  void divideValIntoNodeHomo(Eigen::VectorXd &value, int dir);

  /// assembling K matrix & internal force
  void assemblingKHomo(int log, int dir);

  FEM &fem;
  std::vector<E> &element;
  std::vector<N> &node;
  ValueHomo &value;
  std::vector<std::vector<std::pair<int, int>>> &pairs;
  std::array<double, 3> length;                  ///< length of unit cell
  Eigen::SparseMatrix<double> globalK, reducedK; ///< global stiffness matrix
  Eigen::VectorXd fint;                          ///< internal force
  Eigen::VectorXd reaction;                      ///< reaction force
};

template <class E, class N>
Homogenization<E, N>::Homogenization(
    FEM &fem, std::vector<E> &element, std::vector<N> &node, ValueHomo &value,
    std::vector<std::vector<std::pair<int, int>>> &pairs,
    std::array<double, 3> &length)
    : fem(fem), element(element), node(node), value(value), pairs(pairs),
      length(length), globalK(fem.numeq + fem.dofnp * fem.ndim,
                              fem.numeq + fem.dofnp * fem.ndim),
      reducedK(fem.numeq, fem.numeq), fint(fem.numeq + fem.dofnp * fem.ndim),
      reaction(fem.numeq + fem.dofnp * fem.ndim)
{
  std::cout << "Element class name is " << typeid(E).name() << std::endl;
  std::cout << "Node class name is " << typeid(N).name() << std::endl;

  assert(length.size() == 3);

  setControlPointDOF();
}

template <class E, class N>
Eigen::MatrixXd Homogenization<E, N>::solveNMT(int log)
{
  std::cout << "-- numerical material testing --" << std::endl;

  // init
  initialization();

// solver
#ifdef __INTEL_LLVM_COMPILER
  Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> solv;
  std::cout << "pardiso solver" << std::endl;
#else
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solv;
  std::cout << "Eigen solver" << std::endl;
#endif

  // assembling & reducing K matrix
  assemblingKHomo(log, 0);
  Eigen::SparseMatrix<double> reducedK = globalK;
  reducedK.conservativeResize(fem.numeq, fem.numeq);

  // compute inverse matrix
  solv.compute(reducedK);
  if(log)
  {
    if(solv.info() != Eigen::Success)
      std::cerr << "decomposition failed" << std::endl;
    else
      std::cout << "decomposition succeeded" << std::endl;
  }

  Eigen::MatrixXd CH(fem.voigt, fem.voigt);
  for(int d = 0; d < fem.voigt; d++) // direction loop
  {
    Eigen::VectorXd one(fem.voigt);
    one.setOnes();

    // set macro strain
    Eigen::VectorXd Mdisp = addControlPointValue(one, d, true);

    // make fint
    assemblingKHomo(log, d);
    Eigen::VectorXd reducedfint = fint;
    reducedfint.conservativeResize(fem.numeq);

    double fnorm = sqrt(reducedfint.dot(reducedfint));

    // get reaction displacement
    reducedfint = -reducedfint;
    value.values[d] = solv.solve(reducedfint);

    // divide displacement & add the forced displacement to each node
    divideValIntoNodeHomo(value.values[d], d);

    // get reaction force
    assemblingKHomo(log, d);

    // check convergence of residual force
    if(log)
    {
      reducedfint = fint;
      reducedfint.conservativeResize(fem.numeq);
      double norm = sqrt(reducedfint.dot(reducedfint));
      double residual = norm / fnorm;

      std::cout << "first norm " << "      res " << std::endl;
      std::cout << fnorm << "       " << residual << std::endl;
    }

    // get macro stress
    CH.col(d) = getMstress(reaction);
  }

  return CH;
}

template <class E, class N>
void Homogenization<E, N>::solveLocal(Eigen::VectorXd &barE,
                                      Eigen::VectorXd &barS, int log)
{
  std::cout << "----- localization analysis -----" << std::endl;

  // init
  initialization();

// solver
#ifdef __INTEL_LLVM_COMPILER
  Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> solv;
#else
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solv;
#endif

  // compute micro boundary values from macro strain
  Eigen::VectorXd Mstrain(fem.dofnp * fem.ndim),
      sum_fint(fem.numeq + fem.dofnp * fem.ndim);
  Mstrain.setZero();
  sum_fint.setZero();
  for(int d = 0; d < fem.voigt; d++)
  {
    Mstrain += addControlPointValue(barE, d, true);
    // make K & fint
    assemblingKHomo(log, d);

    sum_fint += fint;
  }

  // compute inverse matrix
  Eigen::SparseMatrix<double> reducedK = globalK;
  reducedK.conservativeResize(fem.numeq, fem.numeq);
  solv.compute(reducedK);
  if(log)
  {
    if(solv.info() != Eigen::Success)
      std::cerr << "decomposition failed" << std::endl;
    else
      std::cout << "decomposition succeeded" << std::endl;
  }

  sum_fint = -sum_fint;
  Eigen::VectorXd reducedfint = sum_fint;
  reducedfint.conservativeResize(fem.numeq);

  // get reaction displacement
  value.val = solv.solve(reducedfint);

  // get reaction force
  Eigen::VectorXd allVal(fem.numeq + fem.dofnp * fem.ndim);
  allVal << value.val, Mstrain;
  Eigen::VectorXd react = globalK * allVal;

  // divide displacement & add the forced displacement to each node
  divideValIntoNode(fem.dofnp, fem.numeq, value.val, node);

  for(int d = 0; d < fem.voigt; d++)
    Eigen::VectorXd dummy = addControlPointValue(barE, d, false);

  return;
}

template <class E, class N>
Eigen::MatrixXd Homogenization<E, N>::solveNLNMT(
    const Eigen::VectorXd &MstrainIncre, double tolNR, int maxIter,
    std::string solvName, int log, bool forcedConti)
{
  std::cout << "-- numerical material testing --" << std::endl;

  Eigen::MatrixXd CH(fem.voigt, fem.voigt);
  for(int d = 0; d < fem.voigt; d++) // direction loop
  {
    // initialization
    int iterstep = 0;
    double firstNorm, residual = 1, residualold;
    Eigen::VectorXd du = Eigen::VectorXd::Zero(fem.numeq);
    Eigen::VectorXd sol = Eigen::VectorXd::Zero(fem.numeq);

    // set macro strain
    Eigen::VectorXd Mdisp = addControlPointValue(MstrainIncre, d, true);

    /////////////////////////////////////////////////////////////////
    //////////////////////////start NR loop//////////////////////////
    /////////////////////////////////////////////////////////////////
    std::cout << "------------ NR analysis ------------" << std::endl;
    std::cout << " iter " << " first norm " << "      res " << std::endl;

    while(residual > tolNR)
    {
      iterstep++;

      // assembling K & fint
      assemblingKHomo(log, d);

      // contraction
      Eigen::VectorXd reducedfint = fint;
      reducedfint.conservativeResize(fem.numeq);
      Eigen::SparseMatrix<double> reducedKtmp = globalK;
      reducedKtmp.conservativeResize(fem.numeq, fem.numeq);
      reducedK = reducedKtmp;

      /// convergence check
      double norm = sqrt(reducedfint.dot(reducedfint));
      if(iterstep == 1)
        firstNorm = norm;

      residual = norm / firstNorm;
      std::cout << "   " << iterstep << "      " << firstNorm << "        "
                << residual << std::endl;

      if(residual < tolNR)
        break;

      if(iterstep == maxIter)
      {
        std::cerr << "itaration step reach " << maxIter
                  << " in homogenization.hpp" << std::endl;
        if(forcedConti)
          return CH;
        else
          exit(1);
      }
      if(iterstep > 3 && residual > residualold)
      {
        std::cerr << "NR is not converge in homogenization.hpp" << std::endl;
        if(forcedConti)
          return CH;
        else
          exit(1);
      }
      residualold = residual;

      // get reaction displacement
      reducedfint = -reducedfint;

      // solve
      sol.setZero();
      Solvers solv(reducedK, reducedfint, sol);
      solv.solve(solvName, log);

      // update displacement
      value.values[d] += sol;
      du += sol;

      // divide displacement & add the forced displacement to each node
      divideValIntoNodeHomo(sol, d);
    }

    // get macro stress
    CH.col(d) = getMstress(reaction);
  }

  return CH;
}

template <class E, class N>
void Homogenization<E, N>::solveNLLocal(Eigen::VectorXd &barE,
                                        Eigen::VectorXd &barS, double tolNR,
                                        int maxIter, int log)
{
}

template <class E, class N> void Homogenization<E, N>::initialization()
{
  // each nodal values
  for(int i = 0; i < fem.numnp; i++)
    for(int j = 0; j < 3; j++)
    {
      node[i].val[j] = 0.0;
      for(int k = 0; k < 6; k++)
      {
        node[i].nmtval[k][j] = 0.0;
      }
    }
  // Value class
  for(int i = 0; i < (int)value.values.size(); i++)
    value.values[i].setZero();
}

////////////////////////////////////////////////////////////////////////
/////////////////////////////////private////////////////////////////////
////////////////////////////////////////////////////////////////////////

template <class E, class N> void Homogenization<E, N>::setControlPointDOF()
{
  // reset previous MPC class
  for(int i = 0; i < fem.numnp; i++)
    if(node[i].mpc != nullptr)
      node[i].mpc.reset();

  for(int dim = 0; dim < (int)pairs.size(); dim++)
  {
    for(int i = 0; i < (int)pairs[dim].size(); i++)
    {
      int slave = pairs[dim][i].second;
      if(node[slave].mpc == nullptr)
      {
        std::shared_ptr<MPC> mpc(new MPC);
        node[slave].mpc = mpc;
        node[slave].mpc->cdof.resize(fem.dofnp);
      }
      for(int j = 0; j < fem.dofnp; j++)
      {
        int dof = fem.numeq + fem.dofnp * dim + j;
        node[slave].mpc->cdof[j].push_back(dof);
      }
    }
  }
}

template <class E, class N>
Eigen::VectorXd
Homogenization<E, N>::addControlPointValue(const Eigen::VectorXd &Mstrain,
                                           int dir, bool to_nmtval)
{
  Eigen::VectorXd val(fem.dofnp * fem.ndim);
  // 2D
  if(fem.ndim == 2)
  {
    // heat problem
    if(fem.voigt == 2)
    {
      if(dir == 0) // x
        val << length[0], 0.0;
      else if(dir == 1) // y
        val << 0.0, length[1];
      else
      {
        std::cerr << "error in homogenization.hpp" << std::endl;
        exit(1);
      }
    }
    // structural problem
    else if(fem.voigt == 3)
    {
      if(dir == 0) // x
        val << length[0], 0.0, 0.0, 0.0;
      else if(dir == 1) // y
        val << 0.0, 0.0, 0.0, length[1];
      else if(dir == 2) // xy
        val << 0.0, length[0] / 2.0, length[1] / 2.0, 0.0;
      else
      {
        std::cerr << "error in homogenization.hpp" << std::endl;
        exit(1);
      }
    }
    else
    {
      std::cerr << "error in homogenization.hpp" << std::endl;
      exit(1);
    }
  }
  // 3D
  else if(fem.ndim == 3)
  {
    // heat problem
    if(fem.voigt == 3)
    {
      if(dir == 0) // x
        val << length[0], 0.0, 0.0;
      else if(dir == 1) // y
        val << 0.0, length[1], 0.0;
      else if(dir == 2) // z
        val << 0.0, 0.0, length[2];
      else
      {
        std::cerr << "error in homogenization.hpp" << std::endl;
        exit(1);
      }
    }
    // structure problem
    else if(fem.voigt == 6)
    {
      if(dir == 0)                        // xx
        val << length[0], 0.0, 0.0,       // point 1
            0.0, 0.0, 0.0,                // point 2
            0.0, 0.0, 0.0;                // point 3
      else if(dir == 1)                   // yy
        val << 0.0, 0.0, 0.0,             // point 1
            0.0, length[1], 0.0,          // point 2
            0.0, 0.0, 0.0;                // point 3
      else if(dir == 2)                   // zz
        val << 0.0, 0.0, 0.0,             // point 1
            0.0, 0.0, 0.0,                // point 2
            0.0, 0.0, length[2];          // point 3
      else if(dir == 3)                   // xy
        val << 0.0, 0.5 * length[0], 0.0, // point 1
            0.5 * length[1], 0.0, 0.0,    // point 2
            0.0, 0.0, 0.0;                // point 3
      else if(dir == 4)                   // yz
        val << 0.0, 0.0, 0.0,             // point 1
            0.0, 0.0, 0.5 * length[1],    // point 2
            0.0, 0.5 * length[2], 0.0;    // point 3
      else if(dir == 5)                   // xz
        val << 0.0, 0.0, 0.5 * length[0], // point 1
            0.0, 0.0, 0.0,                // point 2
            0.5 * length[2], 0.0, 0.0;    // point 3
      else
      {
        std::cerr << "error in homogenization.hpp" << std::endl;
        exit(1);
      }
    }
    else
    {
      std::cerr << "error in homogenization.hpp" << std::endl;
      exit(1);
    }
  }

  // add E bar
  val = Mstrain[dir] * val;

  // add control point's values
  for(int i = 0; i < fem.numnp; i++)
  {
    if(node[i].mpc != nullptr)
    {
      for(int j = 0; j < fem.dofnp; j++)
      {
        for(int k = 0; k < (int)node[i].mpc->cdof[j].size(); k++)
        {
          if(to_nmtval)
            node[i].nmtval[dir][j] += val(node[i].mpc->cdof[j][k] - fem.numeq);
          else
            node[i].val[j] += val(node[i].mpc->cdof[j][k] - fem.numeq);
        }
      }
    }
  }

  return val;
}

template <class E, class N>
Eigen::VectorXd Homogenization<E, N>::getMstress(Eigen::VectorXd &react)
{
  Eigen::VectorXd F = react.block(fem.numeq, 0, fem.ndim * fem.dofnp, 1);

  Eigen::VectorXd Mstress(fem.voigt);

  if(fem.ndim == 2)
  {
    if(fem.voigt == 3) // structure: xx yy xy
      Mstress << F[0] / length[1], F[3] / length[0], F[1] / length[1];
    else if(fem.voigt == 2) // heat: xx yy
      Mstress << F[0] / length[1], F[1] / length[0];
    else
    {
      std::cerr << "error in homogenization.hpp" << std::endl;
      exit(1);
    }
  }
  else
  {
    if(fem.voigt == 6) // structure: xx yy zz xy yz xz
      Mstress << F[0] / (length[1] * length[2]), F[4] / (length[2] * length[0]),
          F[8] / (length[0] * length[1]), F[1] / (length[1] * length[2]),
          F[5] / (length[0] * length[2]), F[2] / (length[1] * length[2]);
    else if(fem.voigt == 3) // heat
      Mstress << F[0] / (length[1] * length[2]), F[1] / (length[0] * length[2]),
          F[2] / (length[0] * length[1]);
    else
    {
      std::cerr << "error in homogenization.hpp" << std::endl;
      exit(1);
    }
  }

  return Mstress;
}

template <class E, class N>
void Homogenization<E, N>::divideValIntoNodeHomo(Eigen::VectorXd &value,
                                                 int dir)
{
  int dofnptmp = std::min(3, fem.dofnp);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < (int)node.size(); i++)
    for(int j = 0; j < dofnptmp; j++)
      node[i].nmtval[dir][j] += value.coeff(node[i].dof[j]);
}

template <class E, class N>
void Homogenization<E, N>::assemblingKHomo(int log, int dir)
{
  if(log)
    std::cout << "assembling...";
  fflush(stdout);

  fint.setZero();
  globalK.setZero();
  reaction.setZero();

  Eigen::SparseMatrix<double> fintMat(fem.numeq + fem.dofnp * fem.ndim,
                                      fem.numeq + fem.dofnp * fem.ndim),
      globalR(fem.numeq + fem.dofnp * fem.ndim,
              fem.numeq + fem.dofnp * fem.ndim);
  std::vector<Eigen::Triplet<double>> tripletsK, tripletsF, tripletsR;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // private variables
    std::vector<Eigen::Triplet<double>> tripK_pri, tripF_pri, tripR_pri;

#ifdef _OPENMP
#pragma omp for nowait
#endif
    for(int nel = 0; nel < fem.nelx; nel++)
      element[nel].makeKeHomo(tripK_pri, tripF_pri, tripR_pri, fem, node, dir);

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      // connect
      tripletsK.insert(tripletsK.end(), tripK_pri.begin(), tripK_pri.end());
      tripletsF.insert(tripletsF.end(), tripF_pri.begin(), tripF_pri.end());
      tripletsR.insert(tripletsR.end(), tripR_pri.begin(), tripR_pri.end());
    }
  }

  globalK.setFromTriplets(tripletsK.begin(), tripletsK.end());
  fintMat.setFromTriplets(tripletsF.begin(), tripletsF.end());
  globalR.setFromTriplets(tripletsR.begin(), tripletsR.end());

  fint = fintMat.diagonal();
  reaction = globalR.diagonal();
  globalK.makeCompressed();

  if(log)
    std::cout << "done" << std::endl;
  fflush(stdout);
}

} // namespace icarat