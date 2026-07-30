#include "opt_buckling_nonlinear.hpp"
#include <Eigen/Sparse>
#include <basis/bmatrix.hpp>
#include <basis/bmatrix_total.hpp>
#include <misc/eigen_manipulation.hpp>
#ifdef __INTEL_LLVM_COMPILER
#include <Eigen/PardisoSupport>
#endif
using namespace std;
using namespace Eigen;
namespace icarat
{
namespace opt_buckling_nonlinear
{
namespace
{
double MM(FEMBuckling &fem, ValueEigen &value)
{
  double sumEigen = 0.0;

  for(int i = 0; i < fem.numeigen; i++)
    sumEigen += pow(1.0 / value.evalues[i], fem.pNorm);

  return pow(sumEigen, 1.0 / fem.pNorm);
}

double dSIMP(BucklingFiniteSaintVenantTotal &actele)
{
  assert(actele.young2() > actele.young1());

  double dens = actele.design_s();
  double young1 = actele.young1();
  double young2 = actele.young2();
  double p = actele.pp();

  double young = (1.0 - pow(dens, p)) * young1 + pow(dens, p) * young2;
  double dyoung = p * pow(dens, p - 1.0) * (young2 - young1);

  double E = dyoung / young;

  return E;
}

MatrixXd dkgdu(FEM &fem, BucklingFiniteSaintVenantTotal &actele,
               std::vector<Node> &node, int ID)
{
  /// 1. get stress in this element
  MatrixXd X = MatrixXd::Zero(actele.ne, fem.ndim);
  MatrixXd U = MatrixXd::Zero(actele.ne, fem.ndim);
  VectorXd Uarr = VectorXd::Zero(actele.ne * fem.ndim);
  VectorXd Uone = VectorXd::Zero(actele.numdof);
  Uone[ID] = 1.0;

  for(int i = 0; i < actele.ne; i++)
  {
    for(int j = 0; j < fem.ndim; j++)
    {
      U.coeffRef(i, j) = node[actele.nodeID[i]].val[j];
      Uarr.coeffRef(fem.ndim * i + j) = node[actele.nodeID[i]].val[j];
      X.coeffRef(i, j) = node[actele.nodeID[i]].x[j];
    }
  }

  double fac = (double)(1.0 / actele.ipmax);
  BmatrixTotal bTotal;
  VectorXd stress = VectorXd::Zero(fem.voigt);
  // gauss point loop start
  for(int ip = 0; ip < actele.ipmax; ip++)
  {
    // nonlinear part
    double jac = 1.0;
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * actele.ne);
    MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
    MatrixXd dFdu = MatrixXd::Zero(fem.ndim, fem.ndim);
    MatrixXd BeL = MatrixXd::Zero(fem.voigt, actele.numdof);
    MatrixXd dBLdu = MatrixXd::Zero(fem.voigt, actele.numdof);
    MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, actele.numdof);
    double dJdu = 0.0;

    // make bmatrix and its derivative for finite strain
    bTotal.make(BeL, BeNL, Ne, jac, F, actele.eType, actele.ipmax, X, U, ip);
    bTotal.makeDifferential(dBLdu, dFdu, dJdu, actele.eType, actele.ipmax, X, U,
                            ip, ID);
    assert(jac > 0.0);

    VectorXd dSdu = actele.De() * (BeL * Uone + dBLdu * Uarr);
    VectorXd S = actele.De() * BeL * Uarr;
    MatrixXd tmp =
        dFdu * Tensor(S, 1) * F.transpose() / F.determinant() +
        F * Tensor(dSdu, 1) * F.transpose() / F.determinant() +
        F * Tensor(S, 1) * dFdu.transpose() / F.determinant() -
        F * Tensor(S, 1) * F.transpose() * dJdu / pow(F.determinant(), 2.0);
    stress += dSdu / actele.ipmax;
  }

  // 2. get stress matrix
  MatrixXd Se = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
  if(fem.ndim == 2)
    for(int i = 0; i < 2; i++)
    {
      Se.coeffRef(2 * i, 2 * i) = stress.coeff(0);
      Se.coeffRef(2 * i, 2 * i + 1) = stress.coeff(2);
      Se.coeffRef(2 * i + 1, 2 * i) = stress.coeff(2);
      Se.coeffRef(2 * i + 1, 2 * i + 1) = stress.coeff(1);
    }
  else
  {
    for(int i = 0; i < 3; i++)
    {
      Se.coeffRef(3 * i, 3 * i) = stress.coeff(0);
      Se.coeffRef(3 * i, 3 * i + 1) = stress.coeff(3);
      Se.coeffRef(3 * i, 3 * i + 2) = stress.coeff(5);
      Se.coeffRef(3 * i + 1, 3 * i) = stress.coeff(3);
      Se.coeffRef(3 * i + 1, 3 * i + 1) = stress.coeff(1);
      Se.coeffRef(3 * i + 1, 3 * i + 2) = stress.coeff(4);
      Se.coeffRef(3 * i + 2, 3 * i) = stress.coeff(5);
      Se.coeffRef(3 * i + 2, 3 * i + 1) = stress.coeff(4);
      Se.coeffRef(3 * i + 2, 3 * i + 2) = stress.coeff(2);
    }
  }

  /// 3. get geometory stiffness matrix dKeG/du
  MatrixXd Beg = MatrixXd::Zero(fem.ndim * fem.ndim, actele.numdof);
  MatrixXd dkegdu = MatrixXd::Zero(actele.numdof, actele.numdof);
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * actele.ne);
  Bmatrix bmatrix;
  double jac = 0.0;
  for(int ip = 0; ip < actele.ipmax; ip++)
  {
    bmatrix.make(Beg, Ne, jac, actele.eType, actele.ipmax, X, ip);
    dkegdu -= Beg.transpose() * Se * Beg * jac;
  }

  return dkegdu;
}
} // namespace

vector<double> sens_compliance(FEM &fem,
                               vector<BucklingFiniteSaintVenantTotal> &element,
                               vector<Node> &node, Force &force, Value &value,
                               const Eigen::SparseMatrix<double> &K,
                               vector<double> &fil_gamma)
{
  vector<double> dfds(fem.nelx);

  // get adjoint vector
  std::cout << "---- computing adjoint vector ----" << endl;
  Solvers solver(K, force.fext, value.val);
  solver.solve("Eigen_CG", 0);

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    double dens = element[nel].design_s();
    double young1 = element[nel].young1();
    double young2 = element[nel].young2();
    double p = element[nel].pp();

    assert(young2 >= young1);

    double young = (1.0 - pow(dens, p)) * young1 + pow(dens, p) * young2;
    double dyoung = p * pow(dens, p - 1.0) * (young2 - young1);

    double E = dyoung / young;

    /// Allocate final displacement and internal force to each element
    VectorXd valEle = VectorXd::Zero(element[nel].numdof);
    for(int i = 0; i < element[nel].ne; i++)
      for(int j = 0; j < fem.ndim; j++)
      {
        int dof = node[element[nel].nodeID[i]].dof[j];

        if(dof < fem.numeq)
          valEle.coeffRef(fem.ndim * i + j) = value.val[dof];
        // dirichlet is set.
        else
          valEle.coeffRef(fem.ndim * i + j) =
              node[element[nel].nodeID[i]].dirich->val[j];
      }

    VectorXd fe_rho = element[nel].gamma() * element[nel].nfe() +
                      (1.0 - element[nel].gamma() * element[nel].gamma()) *
                          element[nel].lfe();

    VectorXd fe_gamma =
        element[nel].nfe() - 2.0 * element[nel].gamma() * element[nel].lfe();

    dfds[nel] = -valEle.transpose() * (E * fe_rho + fe_gamma * fil_gamma[nel]);
  }

  return dfds;
}

vector<double> sens_volume(FEMBuckling &fem,
                           std::vector<BucklingFiniteSaintVenantTotal> &element,
                           double V0)
{
  vector<double> dfds(fem.nelx);
  for(int nel = 0; nel < fem.nelx; nel++)
    dfds[nel] = element[nel].volume / V0;

  return dfds;
}

double func_buckling(FEMBuckling &fem, ValueEigen &value)
{
  double sumEigen = 0.0;

  for(int i = 0; i < fem.numeigen; i++)
    sumEigen += pow(1.0 / value.evalues[i], fem.pNorm);

  double constValue = pow(sumEigen, 1.0 / fem.pNorm) * (fem.limeigen) - 1.0;

  return constValue;
}

vector<double>
sens_buckling(FEMBuckling &fem, vector<BucklingFiniteSaintVenantTotal> &element,
              vector<Node> &node,
              NewtonRaphson<BucklingFiniteSaintVenantTotal, Node> &newton,
              ValueEigen &value)
{
  auto t_start = chrono::system_clock::now();

  std::cout << "----- computing inverse K -----" << endl;
#ifdef __INTEL_LLVM_COMPILER
  PardisoLDLT<SparseMatrix<double>> solv;
  std::cout << "pardiso solver" << std::endl;
#else
  SimplicialLDLT<SparseMatrix<double>> solv;
  std::cout << "Eigen solver" << std::endl;
#endif

  solv.compute(newton.K());

  VectorXd dkds = VectorXd::Zero(fem.nelx);

  // eigenvalue loop start
  for(int eig = 0; eig < fem.numeigen; eig++)
  {
    double kappa = -1.0 / value.evalues[eig];
    ////////////////////////////////////
    /// 1. make force for adjoint problem
    ////////////////////////////////////
    vector<VectorXd> eleDualForce(fem.nelx);
    for(int nel = 0; nel < fem.nelx; nel++)
    {
      // eigenvector in an element
      VectorXd evecEle = VectorXd::Zero(element[nel].numdof);
      for(int i = 0; i < element[nel].ne; i++)
      {
        for(int j = 0; j < fem.ndim; j++)
        {
          int dof = node[element[nel].nodeID[i]].dof[j];
          // not constrainted value
          if(dof < fem.numeq)
            evecEle.coeffRef(fem.ndim * i + j) = value.evectors[eig].coeff(dof);
          // constrainted value
          else
            evecEle.coeffRef(fem.ndim * i + j) = 0.0;
        }
      }

      eleDualForce[nel].resize(element[nel].numdof);
      for(int i = 0; i < element[nel].numdof; i++)
      {
        // partial difference term dKG/du
        MatrixXd dkgdutmp = dkgdu(fem, element[nel], node, i);
        eleDualForce[nel].coeffRef(i) =
            evecEle.transpose() * dkgdutmp * evecEle;
      }
    }

    // assembling force of dual problem &solve
    VectorXd dualForce = elementToNumeq<BucklingFiniteSaintVenantTotal, Node>(
        fem, element, node, eleDualForce);

    std::cout << "----- computing adjoint vector" + to_string(eig) + " -----"
              << endl;
    VectorXd vv = solv.solve(dualForce);

    // send adjoint vector to each element
    vector<VectorXd> vvEle =
        numeqToElement<BucklingFiniteSaintVenantTotal, Node>(fem, element, node,
                                                             vv);

    //////////////////////////////////
    /// 2. make dkdsj term in sensitivity
    //////////////////////////////////
    VectorXd dkdsj = VectorXd::Zero(fem.nelx);
    for(int nel = 0; nel < fem.nelx; nel++)
    {
      double dEds = dSIMP(element[nel]);

      // eigenvector in an element
      VectorXd evecEle = VectorXd::Zero(element[nel].numdof);
      VectorXd disp = VectorXd::Zero(element[nel].numdof);
      for(int i = 0; i < element[nel].ne; i++)
      {
        for(int j = 0; j < fem.ndim; j++)
        {
          int dof = node[element[nel].nodeID[i]].dof[j];
          // not constrainted value
          if(dof < fem.numeq)
            evecEle.coeffRef(fem.ndim * i + j) = value.evectors[eig].coeff(dof);
          // constrainted value
          else
            evecEle.coeffRef(fem.ndim * i + j) = 0.0;

          disp.coeffRef(fem.ndim * i + j) = node[element[nel].nodeID[i]].val[j];
        }
      }

      double dkds_former = dEds * evecEle.transpose() *
                           (element[nel].KGe() + kappa * element[nel].Ke()) *
                           evecEle;
      double dkds_latter = dEds * vvEle[nel].dot(element[nel].fe());

      dkdsj.coeffRef(nel) = dkds_former - dkds_latter;
    }
    dkdsj *= pow(-kappa, fem.pNorm - 1.0);
    dkds += dkdsj;
  } // eigenvalue loop end

  //////////////////////////////////
  // 3. finally get buckling sensitivity dgds
  // & get dfds(sensitivity of objective function)
  //////////////////////////////////
  vector<double> dfds(fem.nelx);
  for(int nel = 0; nel < fem.nelx; nel++)
    dfds[nel] = fem.limeigen * pow(MM(fem, value), 1.0 - fem.pNorm) * dkds[nel];

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  std::cout << "Computation time: " << t_sec << " sec " << std::endl;

  return dfds;
}

} // namespace opt_buckling_nonlinear
} // namespace icarat
