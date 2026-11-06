#include "neo_hooke_damage.hpp"
#include <Eigen/SparseLU>
#include <basis/bmatrix.hpp>
#include <basis/bmatrix_total.hpp>
#include <cmath>
#include <stdexcept>

namespace icarat
{
namespace opt_finite_se_damage_filter
{
using namespace Eigen;
using namespace std;

void calCFJ(FEM &, MatrixXd &C, MatrixXd &C_inv, MatrixXd &F, MatrixXd &F_inv,
            double &J)
{
  C = F.transpose() * F;
  C_inv = C.inverse();
  F_inv = F.inverse();
  J = F.determinant();
}

namespace
{
VectorXd stressVoigt(const MatrixXd &S)
{
  VectorXd v(6);
  v << S(0, 0), S(1, 1), S(2, 2), S(0, 1), S(1, 2), S(2, 0);
  return v;
}
} // namespace

// Conditional sensitivity of -integral(S:E), with D_all_ held fixed.
// tau is filtered BEFORE each load step and held fixed during the NR solve.
// Differentiating the entire loading path requires a separate history adjoint;
// the final K alone cannot supply derivatives of those previous states.
void sensitivity(FEM &fem, vector<FiniteNeoHookedamage> &element,
                 vector<Node> &node, Value &value,
                 const SparseMatrix<double> &K, Optimize &optim, double V0,
                 vector<double> &fil_gamma, double poisson)
{
  (void)value; // Do not overwrite the converged primal displacement.
  (void)K; // Its damage tangent is not the derivative at fixed filtered tau.
  if(fem.ndim != 3 || fem.voigt != 6)
    throw runtime_error(
        "Damage sensitivity requires the implemented 3D stress model.");
  if(!(V0 > 0.0) || fil_gamma.size() != element.size())
    throw runtime_error("Invalid reference volume or gamma derivative size.");

  vector<VectorXd> forceRho(fem.nelx);
  VectorXd objectiveRho = VectorXd::Zero(fem.nelx);
  VectorXd objectiveU = VectorXd::Zero(fem.numeq);
  vector<Triplet<double>> triplets;
  const MatrixXd I = MatrixXd::Identity(3, 3);
  Gauss9Square gauss; // Match the objective weights in makeKeFint exactly.

  for(int nel = 0; nel < fem.nelx; ++nel)
  {
    auto &e = element[nel];
    // The forward implementation uses Gauss9Square weights even in 3D.
    if(e.ipmax > 9 || e.D_all_.size() != static_cast<size_t>(e.ipmax))
      throw runtime_error("Unsupported quadrature or missing damage state.");
    const double rho = e.design_s();
    const double young =
        e.young1() + pow(rho, e.pp()) * (e.young2() - e.young1());
    if(!(young > 0.0))
      throw runtime_error("Nonpositive interpolated Young's modulus.");
    const double ratio =
        e.pp() * pow(rho, e.pp() - 1.0) * (e.young2() - e.young1()) / young;
    const double mu = young / (2.0 * (1.0 + poisson));
    const double kappa =
        poisson * young / ((1.0 + poisson) * (1.0 - 2.0 * poisson)) +
        2.0 * mu / 3.0;
    const double gamma = e.gamma();
    const double dgamma = fil_gamma[nel];
    MatrixXd X(e.ne, 3), U(e.ne, 3);
    VectorXd u(e.numdof);
    vector<int> dofs(e.numdof);
    for(int i = 0; i < e.ne; ++i)
      for(int j = 0; j < 3; ++j)
      {
        const int s = 3 * i + j;
        X(i, j) = node[e.nodeID[i]].x[j];
        u(s) = node[e.nodeID[i]].val[j];
        U(i, j) = gamma * u(s);
        dofs[s] = node[e.nodeID[i]].dof[j];
      }

    MatrixXd nonlinearK = MatrixXd::Zero(e.numdof, e.numdof);
    MatrixXd linearK = nonlinearK;
    VectorXd nonlinearForce = VectorXd::Zero(e.numdof);
    VectorXd productU = nonlinearForce;
    double product = 0.0;
    BmatrixTotal total;
    Bmatrix linear;
    for(int ip = 0; ip < e.ipmax; ++ip)
    {
      double jac = 0.0;
      MatrixXd N(3, e.numdof), BL(6, e.numdof), BNL(9, e.numdof), F(3, 3);
      total.make(BL, BNL, N, jac, F, e.eType, e.ipmax, X, U, ip);
      const double J = F.determinant();
      if(!(jac > 0.0) || !(J > 0.0))
        throw runtime_error("Invalid deformation in damage sensitivity.");
      const MatrixXd C = F.transpose() * F;
      const MatrixXd Cinv = C.inverse();
      const MatrixXd Finv = F.inverse();
      const MatrixXd strain = 0.5 * (C - I);
      const double a = pow(J, -2.0 / 3.0);
      const double damageFactor = 1.0 - e.D_all_[ip];
      const MatrixXd S = damageFactor * (mu * a * (I - C.trace() / 3.0 * Cinv) +
                                         kappa * J * (J - 1.0) * Cinv);
      const VectorXd sv = stressVoigt(S);
      const double weight = gauss.weights(ip);
      product += (S.array() * strain.array()).sum() * weight * jac;
      // makeKeFint uses jac alone for forces and stiffness, but w*jac for S:E.
      nonlinearForce += BL.transpose() * sv * jac;
      for(int s = 0; s < e.numdof; ++s)
      {
        // The library's 3D makeDifferential does not populate dF and uses
        // 2D indices for dBL. Recover reference gradients from BNL instead.
        MatrixXd dBL = MatrixXd::Zero(6, e.numdof);
        MatrixXd dF = MatrixXd::Zero(3, 3);
        for(int k = 0; k < 3; ++k)
          dF(s % 3, k) = BNL(3 * k, 3 * (s / 3));
        const int vi[6] = {0, 1, 2, 0, 1, 2};
        const int vj[6] = {0, 1, 2, 1, 2, 0};
        for(int t = 0; t < e.numdof; ++t)
          for(int v = 0; v < 6; ++v)
          {
            const int i = vi[v], j = vj[v];
            dBL(v, t) = dF(t % 3, i) * BNL(3 * j, 3 * (t / 3));
            if(i != j)
              dBL(v, t) += dF(t % 3, j) * BNL(3 * i, 3 * (t / 3));
          }
        const MatrixXd dC = dF.transpose() * F + F.transpose() * dF;
        const MatrixXd dCinv = -Cinv * dC * Cinv;
        const double dJ = J * (Finv * dF).trace();
        const MatrixXd dS =
            damageFactor *
            (mu * a *
                 (-2.0 / 3.0 * dJ / J * (I - C.trace() / 3.0 * Cinv) -
                  dC.trace() / 3.0 * Cinv - C.trace() / 3.0 * dCinv) +
             kappa * ((2.0 * J - 1.0) * dJ * Cinv + J * (J - 1.0) * dCinv));
        productU(s) += ((dS.array() * strain.array()).sum() +
                        0.5 * (S.array() * dC.array()).sum()) *
                       weight * jac;
        nonlinearK.col(s) +=
            (dBL.transpose() * sv + BL.transpose() * stressVoigt(dS)) * jac;
      }
      MatrixXd B(6, e.numdof), D(6, 6);
      linear.make(B, N, jac, e.eType, e.ipmax, X, ip);
      e.De(D, young, poisson);
      linearK += B.transpose() * D * B * jac;
    }

    const VectorXd linearForce = linearK * u;
    const VectorXd force =
        gamma * nonlinearForce + (1.0 - gamma * gamma) * linearForce;
    forceRho[nel] =
        ratio * force + dgamma * (nonlinearForce + gamma * nonlinearK * u -
                                  2.0 * gamma * linearForce);
    objectiveRho(nel) = ratio * product + dgamma * productU.dot(u);
    const MatrixXd tangent =
        gamma * gamma * nonlinearK + (1.0 - gamma * gamma) * linearK;
    for(int s = 0; s < e.numdof; ++s)
      if(dofs[s] >= 0 && dofs[s] < fem.numeq)
      {
        objectiveU(dofs[s]) += gamma * productU(s);
        for(int t = 0; t < e.numdof; ++t)
          if(dofs[t] >= 0 && dofs[t] < fem.numeq)
            triplets.emplace_back(dofs[s], dofs[t], tangent(s, t));
      }
  }

  SparseMatrix<double> tangent(fem.numeq, fem.numeq);
  tangent.setFromTriplets(triplets.begin(), triplets.end());
  SparseMatrix<double> transpose = tangent.transpose();
  SparseLU<SparseMatrix<double>> solver;
  solver.compute(transpose);
  if(solver.info() != Success)
    throw runtime_error("Damage sensitivity adjoint factorization failed: " +
                        solver.lastErrorMessage());
  const VectorXd adjoint = solver.solve(objectiveU);
  if(solver.info() != Success || !adjoint.allFinite())
    throw runtime_error("Damage sensitivity adjoint solve failed.");

  for(int nel = 0; nel < fem.nelx; ++nel)
  {
    const auto &e = element[nel];
    double indirect = 0.0;
    for(int i = 0; i < e.ne; ++i)
      for(int j = 0; j < 3; ++j)
      {
        const int dof = node[e.nodeID[i]].dof[j];
        // Prescribed displacements have zero adjoint, including nonzero BCs.
        if(dof >= 0 && dof < fem.numeq)
          indirect += adjoint(dof) * forceRho[nel](3 * i + j);
      }
    optim.dfds[nel] = -objectiveRho(nel) + indirect;
    optim.dgds[nel] = e.volume / V0;
  }
}
} // namespace opt_finite_se_damage_filter
} // namespace icarat
