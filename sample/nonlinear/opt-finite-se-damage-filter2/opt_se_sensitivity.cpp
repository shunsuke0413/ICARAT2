#include "opt_finite_strain.hpp"
#include <Eigen/SparseLU>
#include <Eigen/LU>
#include <basis/bmatrix.hpp>
#include <basis/bmatrix_total.hpp>
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace icarat
{
namespace opt_finite_strain
{
using namespace Eigen;
using namespace std;
namespace
{
DamageLawValue damageLaw(double tau)
{
  return evaluateDamageLaw(tau);
}
VectorXd stressVoigt(const MatrixXd &S)
{
  VectorXd v(S.rows() == 2 ? 3 : 6);
  if(S.rows() == 2)
    v << S(0, 0), S(1, 1), S(0, 1);
  else
    v << S(0, 0), S(1, 1), S(2, 2), S(0, 1), S(1, 2), S(2, 0);
  return v;
}
VectorXd strainVoigt(const MatrixXd &E)
{
  VectorXd v = stressVoigt(E);
  v.tail(E.rows() == 2 ? 1 : 3) *= 2.0;
  return v;
}
struct Point
{
  VectorXd strain, stress0;
  MatrixXd strainU, stress0U; // derivatives w.r.t. effective displacement gamma*u
  double weight;
};
struct ElementState
{
  vector<Point> points;
  vector<int> dofs;
  VectorXd u, force, forceD, forceRho, tauU;
  MatrixXd tangent;
  double tau = 0.0, tauRho = 0.0, ratio = 0.0, gamma = 0.0, dgamma = 0.0;
};
ElementState evaluate(FEM &fem, FiniteNeoHookedamage &e,
                      const vector<Node> &node, const DamageFrame &frame,
                      int nel, double dgamma, bool derivatives)
{
  const int ndim = fem.ndim;
  if((ndim != 2 && ndim != 3) || fem.voigt != (ndim == 2 ? 3 : 6) ||
     e.numdof != ndim * e.ne || e.ipmax <= 0 || e.ipmax > 9)
    throw runtime_error("Unsupported dimension, DOFs or quadrature in damage adjoint.");
  if(ndim == 2 && e.mattype() != "plane_strain")
    throw runtime_error("2D damage adjoint supports plane_strain only.");
  const double rho = e.design_s();
  const double young = SIMP(e.young1(), e.young2(), rho, e.pp());
  if(!(young > 0.0))
    throw runtime_error("Nonpositive interpolated Young's modulus.");
  ElementState result;
  result.ratio = e.pp() * pow(rho, e.pp() - 1.0) *
                 (e.young2() - e.young1()) / young;
  result.gamma = e.gamma();
  result.dgamma = dgamma;
  result.u.resize(e.numdof);
  result.dofs.resize(e.numdof);
  result.tauU = VectorXd::Zero(e.numdof);
  const double mu = young / (2.0 * (1.0 + e.poisson()));
  const double kappa = young * e.poisson() /
                           ((1.0 + e.poisson()) * (1.0 - 2.0 * e.poisson())) +
                       2.0 * mu / 3.0;
  MatrixXd X(e.ne, ndim), U(e.ne, ndim);
  for(int i = 0; i < e.ne; ++i)
    for(int j = 0; j < fem.ndim; ++j)
    {
      const int s = ndim * i + j;
      X(i, j) = node[e.nodeID[i]].x[j];
      result.u(s) = frame.displacement(e.nodeID[i], j);
      U(i, j) = result.gamma * result.u(s);
      result.dofs[s] = node[e.nodeID[i]].dof[j];
    }
  const MatrixXd I = MatrixXd::Identity(ndim, ndim);
  MatrixXd nK = MatrixXd::Zero(e.numdof, e.numdof), lK = nK;
  VectorXd nForce0 = VectorXd::Zero(e.numdof);
  BmatrixTotal total;
  Bmatrix linear;
  // BmatrixTotal returns det(dX/dr) times the actual element Gauss weight.
  for(int ip = 0; ip < e.ipmax; ++ip)
  {
    double jac = 0.0;
    MatrixXd N = MatrixXd::Zero(ndim, e.numdof);
    MatrixXd BL = MatrixXd::Zero(fem.voigt, e.numdof);
    MatrixXd BNL = MatrixXd::Zero(ndim * ndim, e.numdof), F(ndim, ndim);
    total.make(BL, BNL, N, jac, F, e.eType, e.ipmax, X, U, ip);
    const double J = F.determinant();
    if(!(jac > 0.0) || !(J > 0.0))
      throw runtime_error("Invalid reference or deformation Jacobian in damage adjoint.");
    const MatrixXd C = F.transpose() * F, Cinv = C.inverse(), Finv = F.inverse();
    const double a = pow(J, -2.0 / 3.0);
    const double trace3D = C.trace() + (ndim == 2 ? 1.0 : 0.0);
    const MatrixXd S0 = mu * a * (I - trace3D / 3.0 * Cinv) +
                        kappa * J * (J - 1.0) * Cinv;
    const double psi = 0.5 * mu * (trace3D * a - 3.0) +
                       0.5 * kappa * (J - 1.0) * (J - 1.0);
    const double tau = sqrt(2.0 * abs(psi));
    result.tau += tau / e.ipmax;
    Point p;
    p.strain = strainVoigt(0.5 * (C - I));
    p.stress0 = stressVoigt(S0);
    p.weight = jac; // Integrate the new work objective with each weight once.
    nForce0 += BL.transpose() * p.stress0 * jac;
    if(derivatives)
    {
      p.strainU = BL;
      p.stress0U.resize(fem.voigt, e.numdof);
      // abs(Psi) requires its sign. At exactly zero use the zero subgradient;
      // no epsilon is inserted into the forward constitutive law.
      if(tau > 0.0)
        result.tauU += (psi > 0.0 ? 1.0 : -1.0) *
                       result.gamma * BL.transpose() * p.stress0 / (tau * e.ipmax);
      for(int s = 0; s < e.numdof; ++s)
      {
        MatrixXd dF = MatrixXd::Zero(ndim, ndim);
        for(int k = 0; k < ndim; ++k)
          dF(s % ndim, k) = BNL(ndim * k, ndim * (s / ndim));
        const MatrixXd dC = dF.transpose() * F + F.transpose() * dF;
        const MatrixXd dCinv = -Cinv * dC * Cinv;
        const double dJ = J * (Finv * dF).trace();
        const MatrixXd dS0 =
            mu * a * (-2.0 / 3.0 * dJ / J * (I - trace3D / 3.0 * Cinv) -
                      dC.trace() / 3.0 * Cinv - trace3D / 3.0 * dCinv) +
            kappa * ((2.0 * J - 1.0) * dJ * Cinv + J * (J - 1.0) * dCinv);
        p.stress0U.col(s) = stressVoigt(dS0);
        MatrixXd dBL = MatrixXd::Zero(fem.voigt, e.numdof);
        const int vi2[3] = {0, 1, 0}, vj2[3] = {0, 1, 1};
        const int vi3[6] = {0, 1, 2, 0, 1, 2}, vj3[6] = {0, 1, 2, 1, 2, 0};
        const int *vi = ndim == 2 ? vi2 : vi3;
        const int *vj = ndim == 2 ? vj2 : vj3;
        for(int t = 0; t < e.numdof; ++t)
          for(int v = 0; v < fem.voigt; ++v)
          {
            const int i = vi[v], j = vj[v];
            dBL(v, t) = dF(t % ndim, i) * BNL(ndim * j, ndim * (t / ndim));
            if(i != j)
              dBL(v, t) += dF(t % ndim, j) * BNL(ndim * i, ndim * (t / ndim));
          }
        nK.col(s) += (1.0 - frame.damage(nel)) *
                    (dBL.transpose() * p.stress0 +
                     BL.transpose() * p.stress0U.col(s)) * jac;
      }
      result.tauRho += result.ratio * tau / (2.0 * e.ipmax);
      if(tau > 0.0)
        result.tauRho += (psi > 0.0 ? 1.0 : -1.0) * dgamma *
                        p.stress0.dot(BL * result.u) / (tau * e.ipmax);
    }
    result.points.push_back(std::move(p));
    MatrixXd B = MatrixXd::Zero(fem.voigt, e.numdof), D(fem.voigt, fem.voigt);
    linear.make(B, N, jac, e.eType, e.ipmax, X, ip);
    e.De(D, young, e.poisson());
    lK += B.transpose() * D * B * jac;
  }
  const VectorXd nForce = (1.0 - frame.damage(nel)) * nForce0;
  const VectorXd lForce = lK * result.u;
  const double gamma = result.gamma;
  result.force = gamma * nForce + (1.0 - gamma * gamma) * lForce;
  result.forceD = -gamma * nForce0;
  if(derivatives)
  {
    result.tangent = gamma * gamma * nK + (1.0 - gamma * gamma) * lK;
    result.forceRho = result.ratio * result.force +
                      dgamma * (nForce + gamma * nK * result.u - 2.0 * gamma * lForce);
  }
  return result;
}
void scatter(VectorXd &global, const VectorXd &local, const vector<int> &dofs)
{
  for(int s = 0; s < local.size(); ++s)
    if(dofs[s] >= 0 && dofs[s] < global.size())
      global(dofs[s]) += local(s);
}
VectorXd gather(const VectorXd &global, const vector<int> &dofs)
{
  VectorXd local = VectorXd::Zero(dofs.size());
  for(int s = 0; s < local.size(); ++s)
    if(dofs[s] >= 0 && dofs[s] < global.size())
      local(s) = global(dofs[s]);
  return local;
}
struct Partials
{
  double work = 0.0;
  vector<VectorXd> u, damage;
  VectorXd rho;
};
Partials workPartials(FEM &fem, vector<FiniteNeoHookedamage> &element,
                      const vector<Node> &node, const DamageHistory &history,
                      const vector<double> &dgamma, bool derivatives)
{
  if(history.size() < 2)
    throw runtime_error("Damage work requires initial and converged load states.");
  Partials f;
  f.rho = VectorXd::Zero(fem.nelx);
  f.u.assign(history.size(), VectorXd::Zero(fem.numeq));
  f.damage.assign(history.size(), VectorXd::Zero(fem.nelx));
  for(int nel = 0; nel < fem.nelx; ++nel)
  {
    ElementState a = evaluate(fem, element[nel], node, history[0], nel,
                               dgamma[nel], derivatives);
    for(size_t n = 1; n < history.size(); ++n)
    {
      ElementState b = evaluate(fem, element[nel], node, history[n], nel,
                                 dgamma[nel], derivatives);
      VectorXd fa = VectorXd::Zero(a.u.size()), fb = fa;
      const double da = 1.0 - history[n - 1].damage(nel);
      const double db = 1.0 - history[n].damage(nel);
      for(size_t ip = 0; ip < b.points.size(); ++ip)
      {
        const auto &pa = a.points[ip];
        const auto &pb = b.points[ip];
        const VectorXd delta = pb.strain - pa.strain;
        const VectorXd average = 0.5 * (da * pa.stress0 + db * pb.stress0);
        const double w = pb.weight;
        f.work += average.dot(delta) * w;
        if(derivatives)
        {
          fa += a.gamma * (0.5 * da * pa.stress0U.transpose() * delta -
                           pa.strainU.transpose() * average) * w;
          fb += b.gamma * (0.5 * db * pb.stress0U.transpose() * delta +
                           pb.strainU.transpose() * average) * w;
          f.damage[n - 1](nel) -= 0.5 * pa.stress0.dot(delta) * w;
          f.damage[n](nel) -= 0.5 * pb.stress0.dot(delta) * w;
          const VectorXd ds = 0.5 * dgamma[nel] *
              (da * pa.stress0U * a.u + db * pb.stress0U * b.u);
          const VectorXd de = dgamma[nel] *
              (pb.strainU * b.u - pa.strainU * a.u);
          f.rho(nel) += (a.ratio * average.dot(delta) +
                         ds.dot(delta) + average.dot(de)) * w;
        }
      }
      if(derivatives)
      {
        scatter(f.u[n - 1], fa, a.dofs);
        scatter(f.u[n], fb, b.dofs);
      }
      a = std::move(b);
    }
  }
  return f;
}
void branches(const DamageHistory &history, size_t n, VectorXd &active,
              VectorXd &q)
{
  const int nelx = history[n].damage.size();
  active = VectorXd::Zero(nelx);
  q = active;
  for(int e = 0; e < nelx; ++e)
  {
    const auto law = damageLaw(history[n].tauInput(e));
    active(e) = law.trial > history[n - 1].damage(e) ? 1.0 : 0.0;
    q(e) = active(e) * law.derivative;
    if(abs(history[n].damage(e) - max(history[n - 1].damage(e), law.trial)) > 1.e-10)
      throw runtime_error("Recorded damage does not satisfy the lagged history law.");
  }
}
} // namespace

SparseMatrix<double> damageFilterWeights(
    FEM &fem, const vector<FiniteNeoHookedamage> &element,
    const vector<Node> &node, double radius)
{
  if(!std::isfinite(radius) || radius < 0.0)
    throw runtime_error("Invalid damage filter radius.");
  MatrixXd centers = MatrixXd::Zero(fem.nelx, fem.ndim);
  for(int e = 0; e < fem.nelx; ++e)
    for(int d = 0; d < fem.ndim; ++d)
      for(int i = 0; i < element[e].ne; ++i)
        centers(e, d) += node[element[e].nodeID[i]].x[d] / element[e].ne;
  vector<Triplet<double>> entries;
  for(int e = 0; e < fem.nelx; ++e)
  {
    vector<pair<int, double>> row;
    double sum = 0.0;
    for(int j = 0; j < fem.nelx; ++j)
    {
      const double distance = (centers.row(e) - centers.row(j)).norm();
      if(distance <= radius)
      {
        const double weight = radius - distance;
        row.emplace_back(j, weight);
        sum += weight;
      }
    }
    if(sum > 0.0)
    {
      for(const auto &entry : row)
        if(entry.second > 0.0)
          entries.emplace_back(e, entry.first, entry.second / sum);
    }
    else
      entries.emplace_back(e, e, 1.0); // tauFilter's actual fallback is identity.
  }
  SparseMatrix<double> W(fem.nelx, fem.nelx);
  W.setFromTriplets(entries.begin(), entries.end());
  return W;
}

void initializeDamagePath(FEM &fem, vector<FiniteNeoHookedamage> &element,
                          vector<Node> &node)
{
  for(auto &e : element)
  {
    fill(e.D_all_.begin(), e.D_all_.end(), 0.0);
    fill(e.tau_all_.begin(), e.tau_all_.end(), 0.0);
    // computeTau reads cached young_ before assigning its SIMP value. Prime
    // that cache at zero displacement using the unchanged element routine.
    vector<Triplet<double>> k, f, r;
    double product = 0.0;
    e.makeKeFint(k, f, r, fem, node, 0, product, 0.0);
  }
}

DamageFrame captureDamageFrame(FEM &fem,
    const vector<FiniteNeoHookedamage> &element, const vector<Node> &node,
    const vector<double> &tauInput)
{
  DamageFrame frame;
  frame.displacement.resize(node.size(), fem.ndim);
  frame.damage.resize(fem.nelx);
  frame.tauInput = VectorXd::Zero(fem.nelx);
  if(!tauInput.empty() && tauInput.size() != element.size())
    throw runtime_error("Invalid filtered tau history size.");
  for(size_t i = 0; i < node.size(); ++i)
    for(int j = 0; j < fem.ndim; ++j)
      frame.displacement(i, j) = node[i].val[j];
  for(int e = 0; e < fem.nelx; ++e)
  {
    if(element[e].D_all_.empty())
      throw runtime_error("Missing Gauss-point damage state.");
    frame.damage(e) = element[e].D_all_[0];
    for(double d : element[e].D_all_)
      if(!std::isfinite(d) || abs(d - frame.damage(e)) > 1.e-12)
        throw runtime_error("Element damage must be equal at all Gauss points.");
    if(!tauInput.empty())
      frame.tauInput(e) = tauInput[e];
  }
  if(!frame.displacement.allFinite() || !frame.damage.allFinite())
    throw runtime_error("Nonfinite converged damage state.");
  return frame;
}

double damageWork(FEM &fem, vector<FiniteNeoHookedamage> &element,
                  const vector<Node> &node, const DamageHistory &history)
{
  return workPartials(fem, element, node, history,
                      vector<double>(element.size(), 0.0), false).work;
}

void sensitivity(FEM &fem, vector<FiniteNeoHookedamage> &element,
                 const vector<Node> &node, Optimize &optim, double V0,
                 const vector<double> &fil_gamma, const DamageHistory &history,
                 const SparseMatrix<double> &weights)
{
  if(!(V0 > 0.0) || fil_gamma.size() != element.size() ||
     weights.rows() != fem.nelx || weights.cols() != fem.nelx)
    throw runtime_error("Invalid volume, gamma derivative or damage filter.");
  cout << "Sensitivity: computing objective partial derivatives for "
       << history.size() - 1 << " load steps." << endl;
  auto f = workPartials(fem, element, node, history, fil_gamma, true);
  // Minimize negative internal work. L = -work + lambda^T R + eta^T H.
  VectorXd gradient = -f.rho;
  VectorXd etaNext = VectorXd::Zero(fem.nelx);
  VectorXd activeNext = VectorXd::Zero(fem.nelx), qNext = activeNext;
  for(int n = static_cast<int>(history.size()) - 1; n >= 1; --n)
  {
    cout << "Sensitivity: backward load step " << n << "/" << history.size() - 1 << endl;
    VectorXd active, q;
    branches(history, n, active, q);
    const VectorXd filteredFuture = weights.transpose() *
                                     (qNext.array() * etaNext.array()).matrix();
    VectorXd rhs = f.u[n]; // -f_u for objective -work
    vector<Triplet<double>> entries;
    vector<ElementState> state;
    for(int e = 0; e < fem.nelx; ++e)
    {
      state.push_back(evaluate(fem, element[e], node, history[n], e,
                               fil_gamma[e], true));
      const auto &s = state.back();
      // H[n+1] depends on u[n], whereas H[n] does not. Therefore the
      // same-time reduced matrix is K^(D), with future coupling on the RHS.
      scatter(rhs, filteredFuture(e) * s.tauU, s.dofs);
      for(int i = 0; i < s.tangent.rows(); ++i)
        if(s.dofs[i] >= 0 && s.dofs[i] < fem.numeq)
          for(int j = 0; j < s.tangent.cols(); ++j)
            if(s.dofs[j] >= 0 && s.dofs[j] < fem.numeq)
              entries.emplace_back(s.dofs[i], s.dofs[j], s.tangent(i, j));
    }
    SparseMatrix<double> tangent(fem.numeq, fem.numeq);
    tangent.setFromTriplets(entries.begin(), entries.end());
    VectorXd lambda = VectorXd::Zero(fem.numeq);
    if(fem.numeq > 0)
    {
      SparseMatrix<double> transpose = tangent.transpose();
      SparseLU<SparseMatrix<double>> solver;
      solver.compute(transpose);
      if(solver.info() != Success)
        throw runtime_error("Damage adjoint factorization failed: " + solver.lastErrorMessage());
      lambda = solver.solve(rhs);
      if(solver.info() != Success || !lambda.allFinite() ||
         (transpose * lambda - rhs).norm() > 1.e-8 * max(1.0, rhs.norm()))
        throw runtime_error("Damage adjoint linear residual check failed.");
    }
    VectorXd eta = f.damage[n] +
                   ((VectorXd::Ones(fem.nelx) - activeNext).array() *
                    etaNext.array()).matrix();
    for(int e = 0; e < fem.nelx; ++e)
    {
      const VectorXd localLambda = gather(lambda, state[e].dofs);
      eta(e) -= state[e].forceD.dot(localLambda);
      gradient(e) += state[e].forceRho.dot(localLambda);
    }
    const VectorXd filteredEta = weights.transpose() *
                                (q.array() * eta.array()).matrix();
    VectorXd previousTau(fem.nelx);
    for(int e = 0; e < fem.nelx; ++e)
    {
      const auto previous = evaluate(fem, element[e], node, history[n - 1], e,
                                      fil_gamma[e], true);
      previousTau(e) = previous.tau;
      gradient(e) -= filteredEta(e) * previous.tauRho;
    }
    const VectorXd predicted = weights * previousTau;
    if((predicted - history[n].tauInput).norm() >
       1.e-9 * max(1.0, history[n].tauInput.norm()))
      throw runtime_error("Stored damage input differs from W*tau(u_previous).");
    etaNext = std::move(eta);
    activeNext = std::move(active);
    qNext = std::move(q);
  }
  if(!gradient.allFinite())
    throw runtime_error("Nonfinite full-history design sensitivity.");
  optim.object_f = -f.work;
  for(int e = 0; e < fem.nelx; ++e)
  {
    optim.dfds[e] = gradient(e);
    optim.dgds[e] = element[e].volume / V0;
  }
}

void verifyDamageDerivatives(FEM &fem, vector<FiniteNeoHookedamage> &element,
                            const vector<Node> &node,
                            const DamageHistory &history,
                            const SparseMatrix<double> &weights,
                            const vector<double> &fil_gamma)
{
  auto near = [](double error, double scale, const string &name) {
    const double normalized = error / max(1.0, scale);
    cout << "Derivative check " << name << ": " << normalized << '\n';
    if(!std::isfinite(normalized) || normalized > 2.e-5)
      throw runtime_error("Derivative check failed: " + name);
  };
  for(double tau : {0.2, 0.36, 0.55, 0.8, 2.0, 1000.0})
  {
    const double h = 1.e-6;
    const double fd = (damageLaw(tau + h).trial - damageLaw(tau - h).trial) / (2.0 * h);
    near(abs(fd - damageLaw(tau).derivative), abs(fd), "damage/sigmoid");
    near(abs(damageLaw(tau).trial - element[0].calculateTensileDamage(tau)),
         1.0, "forward damage law");
  }
  const VectorXd ones = VectorXd::Ones(fem.nelx);
  near((weights * ones - ones).norm(), ones.norm(), "filter normalization");
  for(int nel = 0; nel < min(fem.nelx, 3); ++nel)
  {
    const int n = history.size() - 1;
    const auto s = evaluate(fem, element[nel], node, history[n], nel, fil_gamma[nel], true);
    VectorXd direction(s.u.size());
    for(int j = 0; j < direction.size(); ++j)
      direction(j) = sin(0.7 * (j + 1));
    direction.normalize();
    const double h = 1.e-6;
    auto plus = history[n], minus = plus;
    for(int i = 0; i < element[nel].ne; ++i)
      for(int d = 0; d < fem.ndim; ++d)
      {
        plus.displacement(element[nel].nodeID[i], d) += h * direction(fem.ndim * i + d);
        minus.displacement(element[nel].nodeID[i], d) -= h * direction(fem.ndim * i + d);
      }
    const auto sp = evaluate(fem, element[nel], node, plus, nel, 0.0, false);
    const auto sm = evaluate(fem, element[nel], node, minus, nel, 0.0, false);
    const VectorXd fd = (sp.force - sm.force) / (2.0 * h);
    near((fd - s.tangent * direction).norm(), fd.norm(), "fixed-damage residual tangent");
    if(s.tau > 1.e-8)
      near(abs((sp.tau - sm.tau) / (2.0 * h) - s.tauU.dot(direction)),
           abs(s.tauU.dot(direction)), "energy-norm displacement");
    plus = history[n]; minus = plus;
    plus.damage(nel) += h; minus.damage(nel) -= h;
    const auto dp = evaluate(fem, element[nel], node, plus, nel, 0.0, false);
    const auto dm = evaluate(fem, element[nel], node, minus, nel, 0.0, false);
    near(((dp.force - dm.force) / (2.0 * h) - s.forceD).norm(),
         s.forceD.norm(), "residual damage");
    auto ep = element[nel], em = element[nel];
    double rp = ep.design_s() + h, rm = em.design_s() - h;
    ep.setDesignS(rp); em.setDesignS(rm);
    ep.setGamma(ep.gamma() + h * fil_gamma[nel]);
    em.setGamma(em.gamma() - h * fil_gamma[nel]);
    const auto srp = evaluate(fem, ep, node, history[n], nel, 0.0, false);
    const auto srm = evaluate(fem, em, node, history[n], nel, 0.0, false);
    near(((srp.force - srm.force) / (2.0 * h) - s.forceRho).norm(),
         s.forceRho.norm(), "residual design/gamma");
    if(s.tau > 1.e-8)
      near(abs((srp.tau - srm.tau) / (2.0 * h) - s.tauRho),
           abs(s.tauRho), "energy-norm design/gamma");
    // Compare with the unchanged primal element, freezing damage via tau=0.
    auto primalElement = element[nel];
    fill(primalElement.D_all_.begin(), primalElement.D_all_.end(), history[n].damage(nel));
    auto primalNode = node;
    for(size_t i = 0; i < node.size(); ++i)
      for(int d = 0; d < fem.ndim; ++d)
        primalNode[i].val[d] = history[n].displacement(i, d);
    vector<Triplet<double>> k, f, r;
    double product = 0.0;
    primalElement.makeKeFint(k, f, r, fem, primalNode, n, product, 0.0);
    near((primalElement.fe() - s.force).norm(), s.force.norm(), "unchanged primal residual");
    near((primalElement.Ke() - s.tangent).norm(), s.tangent.norm(),
         "primal fixed-input tangent");
    // Exercise retained damage even on a mesh whose actual path is elastic.
    auto damagedFrame = history[n];
    damagedFrame.damage(nel) = 0.4;
    const auto damaged = evaluate(fem, element[nel], node, damagedFrame, nel, 0.0, true);
    fill(primalElement.D_all_.begin(), primalElement.D_all_.end(), 0.4);
    k.clear(); f.clear(); r.clear();
    primalElement.makeKeFint(k, f, r, fem, primalNode, n, product, 0.0);
    near((primalElement.Ke() - damaged.tangent).norm(), damaged.tangent.norm(),
         "primal retained-damage unloading tangent");
  }
  if(history.size() * element.size() <= 2000)
  {
    auto f = workPartials(fem, element, node, history, fil_gamma, true);
    const double h = 1.e-6;
    const int n = max(1, static_cast<int>(history.size() / 2));
    auto plus = history, minus = history;
    VectorXd direction(fem.numeq);
    for(int j = 0; j < direction.size(); ++j)
      direction(j) = cos(0.3 * (j + 1)) / max(1, fem.numeq);
    for(size_t i = 0; i < node.size(); ++i)
      for(int d = 0; d < fem.ndim; ++d)
        if(node[i].dof[d] >= 0 && node[i].dof[d] < fem.numeq)
        {
          plus[n].displacement(i, d) += h * direction(node[i].dof[d]);
          minus[n].displacement(i, d) -= h * direction(node[i].dof[d]);
        }
    const double fdU = (damageWork(fem, element, node, plus) -
                        damageWork(fem, element, node, minus)) / (2.0 * h);
    near(abs(fdU - f.u[n].dot(direction)), abs(fdU), "trapezoidal work displacement");
    plus = history; minus = history;
    plus[n].damage(0) += h; minus[n].damage(0) -= h;
    const double fdD = (damageWork(fem, element, node, plus) -
                        damageWork(fem, element, node, minus)) / (2.0 * h);
    near(abs(fdD - f.damage[n](0)), abs(fdD), "trapezoidal work damage");
    auto ep = element, em = element;
    double rp = ep[0].design_s() + h, rm = em[0].design_s() - h;
    ep[0].setDesignS(rp); em[0].setDesignS(rm);
    ep[0].setGamma(ep[0].gamma() + h * fil_gamma[0]);
    em[0].setGamma(em[0].gamma() - h * fil_gamma[0]);
    const double fdRho = (damageWork(fem, ep, node, history) -
                          damageWork(fem, em, node, history)) / (2.0 * h);
    near(abs(fdRho - f.rho(0)), abs(fdRho), "trapezoidal work design/gamma");
  }
  else
    cout << "Global objective partial checks require a smaller validation mesh.\n";

  // Method A: assemble the complete time-coupled (u,D) Jacobian. Compare its
  // adjoint gradient with elimination (method B) on a small validation mesh.
  const int block = fem.numeq + fem.nelx;
  const int count = history.size() - 1;
  if(block * count > 400)
  {
    cout << "Coupled/reduced comparison requires <=400 total history unknowns.\n";
    return;
  }
  auto f = workPartials(fem, element, node, history, fil_gamma, true);
  MatrixXd A = MatrixXd::Zero(block * count, block * count);
  MatrixXd Arho = MatrixXd::Zero(block * count, fem.nelx);
  VectorXd rhs = VectorXd::Zero(block * count);
  for(int n = 1; n <= count; ++n)
  {
    const int offset = (n - 1) * block;
    rhs.segment(offset, fem.numeq) = f.u[n];
    rhs.segment(offset + fem.numeq, fem.nelx) = f.damage[n];
    VectorXd active, q;
    branches(history, n, active, q);
    MatrixXd tauU = MatrixXd::Zero(fem.nelx, fem.numeq);
    VectorXd tauRho(fem.nelx);
    for(int e = 0; e < fem.nelx; ++e)
    {
      const auto s = evaluate(fem, element[e], node, history[n], e, fil_gamma[e], true);
      const auto previous = evaluate(fem, element[e], node, history[n - 1], e,
                                      fil_gamma[e], true);
      tauRho(e) = previous.tauRho;
      for(int i = 0; i < s.u.size(); ++i)
        if(s.dofs[i] >= 0 && s.dofs[i] < fem.numeq)
        {
          const int row = offset + s.dofs[i];
          A(row, offset + fem.numeq + e) += s.forceD(i);
          Arho(row, e) += s.forceRho(i);
          tauU(e, s.dofs[i]) += previous.tauU(i);
          for(int j = 0; j < s.u.size(); ++j)
            if(s.dofs[j] >= 0 && s.dofs[j] < fem.numeq)
              A(row, offset + s.dofs[j]) += s.tangent(i, j);
        }
      A(offset + fem.numeq + e, offset + fem.numeq + e) = 1.0;
      if(n > 1)
        A(offset + fem.numeq + e, offset - block + fem.numeq + e) = -(1.0 - active(e));
    }
    const MatrixXd denseWeights = MatrixXd(weights);
    Arho.block(offset + fem.numeq, 0, fem.nelx, fem.nelx) =
        (-q).asDiagonal() * denseWeights * tauRho.asDiagonal();
    if(n > 1)
      A.block(offset + fem.numeq, offset - block, fem.nelx, fem.numeq) =
          (-q).asDiagonal() * denseWeights * tauU;
  }
  FullPivLU<MatrixXd> solver(A.transpose());
  if(!solver.isInvertible())
    throw runtime_error("Small coupled history Jacobian is singular.");
  const VectorXd adjoint = solver.solve(rhs);
  near((A.transpose() * adjoint - rhs).norm(), rhs.norm(), "coupled adjoint residual");
  const VectorXd coupledGradient = -f.rho + Arho.transpose() * adjoint;
  Optimize reduced(fem.nelx, 1, 1.e-4, 1.0);
  sensitivity(fem, element, node, reduced, 1.0, fil_gamma, history, weights);
  const VectorXd reducedGradient = Map<const VectorXd>(reduced.dfds.data(), fem.nelx);
  near((coupledGradient - reducedGradient).norm(), coupledGradient.norm(),
       "coupled versus reduced full-history gradient");
}
} // namespace opt_finite_strain
} // namespace icarat
