#include "solvers.hpp"
#include <Eigen/Sparse>

namespace icarat
{
using namespace std;
using namespace Eigen;

//////////////////////////////////////////////////
///////////////////////public/////////////////////
//////////////////////////////////////////////////

void Solvers::CG(int log)
{
  sol.setZero();

  int iter = 0;
  double alpha;
  double up, low;
  double beta;
  Eigen::VectorXd p(A.cols());
  Eigen::VectorXd r(A.cols());
  Eigen::VectorXd rb(A.cols());
  p.setZero();
  r.setZero();
  rb.setZero();

  r = rhs - A * sol;
  p = r;
  while(true)
  {
    iter++;
    up = r.dot(p);
    low = p.dot(A * p);
    alpha = up / low;
    sol += alpha * p;
    rb = r - alpha * A * p;
    if(rb.norm() < tol_)
      break;
    if(iter > maxIter_)
      break;
    up = rb.dot(rb);
    low = r.dot(r);
    beta = up / low;
    p = rb + beta * p;
    r = rb;
  }
  if(log == 1)
  {
    std::cout << "iteration: " << iter << std::endl;
    std::cout << "residual: " << rb.norm() << std::endl;
  }
  return;
}

void Solvers::ScalingCG(int log)
{
  // diagonal scaling
  DiagonalPreconditioner<double> precon(A);

  //----------Initialize----------
  sol.setZero();
  int iter = 0;
  Eigen::VectorXd D = A.diagonal(); // Scaling A matrix
  Eigen::VectorXd rk = rhs - A * sol;
  Eigen::VectorXd pk = precon.solve(rk); // Scaling rk
  Eigen::VectorXd Mrk = pk;
  double bnorm = sqrt(rhs.dot(rhs));
  double Mrkrk = Mrk.dot(rk);

  double rnorm, alpha, beta, Mrkp1rkp1;
  Eigen::VectorXd Apk(A.rows());
  Apk.setZero();
  while(true)
  {
    iter++;
    Apk = A * pk;
    alpha = Mrkrk / pk.dot(Apk);
    sol += alpha * pk;
    rk -= alpha * Apk;
    Mrk = precon.solve(rk); // Scaling rk
    Mrkp1rkp1 = Mrk.dot(rk);
    beta = Mrkp1rkp1 / Mrkrk;
    pk = beta * pk + Mrk;
    Mrkrk = Mrkp1rkp1;
    rnorm = sqrt(rk.dot(rk));

    if(rnorm < tol_ * bnorm)
      break;
    if(iter > maxIter_)
      break;
  }
  if(log == 1)
  {
    std::cout << "iteration: " << iter << std::endl;
    std::cout << "residual: " << rnorm << std::endl;
  }

  return;
}

void Solvers::ICCG(int log)
{
  // Cholesky decomposition
  IncompleteCholesky<double> precon(A);

  //----------Initialize----------
  sol.setZero();
  int iter = 0;
  VectorXd rk = rhs - A * sol;
  VectorXd pk((int)rk.size());
  pk.setZero();
  pk = precon.solve(rk);
  Eigen::VectorXd Mrk = pk;
  double bnorm = sqrt(rhs.dot(rhs));
  double Mrkrk = Mrk.dot(rk);

  double rnorm, alpha, beta, Mrkp1rkp1;
  Eigen::VectorXd Apk(A.rows());
  Apk.setZero();
  while(true)
  {
    iter++;
    Apk = A * pk;
    alpha = Mrkrk / pk.dot(Apk);
    sol += alpha * pk;
    rk -= alpha * Apk;
    Mrk = precon.solve(rk);
    Mrkp1rkp1 = Mrk.dot(rk);
    beta = Mrkp1rkp1 / Mrkrk;
    pk = beta * pk + Mrk;
    Mrkrk = Mrkp1rkp1;
    rnorm = sqrt(rk.dot(rk));

    if(rnorm < tol_ * bnorm)
      break;
    if(iter > maxIter_)
      break;
  }
  if(log == 1)
  {
    std::cout << "iteration: " << iter << std::endl;
    std::cout << "residual: " << rnorm << std::endl;
  }

  return;
}

} // namespace icarat