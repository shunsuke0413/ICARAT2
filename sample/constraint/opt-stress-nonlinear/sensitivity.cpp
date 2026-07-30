#include "opt_stress_nonlinear.hpp"
#include <basis/bmatrix_total.hpp>
#include <misc/eigen_manipulation.hpp>

using namespace std;
using namespace Eigen;
namespace icarat
{
namespace stress_nonlinear
{
double funcStressConst(FEMStressConst &fem,
                       vector<FiniteSaintVenantTotal> &element)
{
  double misesPsum = 0.0;
  for(int i = 0; i < fem.nelx; i++)
  {
    double misesQP = pow(element[i].design_s(), fem.qp) * element[i].young2() *
                     element[i].mises() / element[i].young();
    misesPsum += pow(misesQP / fem.stressLim, fem.pNorm);
  }
  return pow(misesPsum, 1.0 / fem.pNorm) - 1.0;
}

vector<double> sensCompliance(FEMStressConst &fem,
                              vector<FiniteSaintVenantTotal> &element,
                              vector<Node> &node, Force &force, Value &value,
                              const Eigen::SparseMatrix<double> &K,
                              vector<double> &fil_gamma)
{

  /// get adjoint vector
  std::cout << "---- computing adjoint vector ----" << endl;
  Solvers solver(K, force.fext, value.val);
  solver.solve("Eigen_CG", 0);

  vector<double> dfds(fem.nelx);

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

vector<double> sensStressConst(FEMStressConst &fem,
                               vector<FiniteSaintVenantTotal> &element,
                               vector<Node> &node,
                               const SparseMatrix<double> &K)
{
  double misesPsum = 0.0;
  vector<double> misesQP(fem.nelx);
  for(int i = 0; i < fem.nelx; i++)
  {
    misesQP[i] = pow(element[i].design_s(), fem.qp) * element[i].young2() *
                 element[i].mises() / element[i].young();
    misesPsum += pow(misesQP[i] / fem.stressLim, fem.pNorm);
  }

  double dpndvm = pow(misesPsum, (1.0 / fem.pNorm) - 1.0);

  /// make factor1 for stress constraint's sensitivity
  vector<VectorXd> factor1, factor2;
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    VectorXd facEle1 = VectorXd::Zero(fem.voigt);
    vector<double> dvmdsig;
    if(fem.ndim == 2 && element[nel].mattype() == "plane_stress")
    {
      dvmdsig.push_back(
          (2.0 * element[nel].stress()[0] - element[nel].stress()[1]) /
          (2.0 * element[nel].mises()));
      dvmdsig.push_back(
          (2.0 * element[nel].stress()[1] - element[nel].stress()[0]) /
          (2.0 * element[nel].mises()));
      dvmdsig.push_back(3.0 * element[nel].stress()[2] / element[nel].mises());
    }
    else if(fem.ndim == 3)
    {
      dvmdsig.push_back((2.0 * element[nel].stress()[0] -
                         element[nel].stress()[1] - element[nel].stress()[2]) /
                        (2.0 * element[nel].mises()));
      dvmdsig.push_back((2.0 * element[nel].stress()[1] -
                         element[nel].stress()[0] - element[nel].stress()[2]) /
                        (2.0 * element[nel].mises()));
      dvmdsig.push_back((2.0 * element[nel].stress()[2] -
                         element[nel].stress()[0] - element[nel].stress()[1]) /
                        (2.0 * element[nel].mises()));
      dvmdsig.push_back(3.0 * element[nel].stress()[3] / element[nel].mises());
      dvmdsig.push_back(3.0 * element[nel].stress()[4] / element[nel].mises());
      dvmdsig.push_back(3.0 * element[nel].stress()[5] / element[nel].mises());
    }
    else
    {
      cerr << "error in sensitivity.cpp" << endl;
      exit(1);
    }
    for(int i = 0; i < fem.voigt; i++)
      facEle1[i] = dpndvm * pow(misesQP[nel] / fem.stressLim, fem.pNorm - 1) *
                   dvmdsig[i] / fem.stressLim;

    factor1.push_back(facEle1);

    /// make factor2 for stress constraint's sensitivity
    VectorXd facEle2 = VectorXd::Zero(fem.voigt);

    double dvds;
    if(element[nel].young1() < element[nel].young2())
      dvds = fem.qp * pow(element[nel].design_s(), fem.qp - 1.0) *
             element[nel].young2();
    else
    {
      cerr << "young2 need larger parameter than young1" << endl;
      exit(1);
    }

    facEle2 = dvds * element[nel].stress() / element[nel].young();

    factor2.push_back(facEle2);
  }

  /// make factor12
  VectorXd factor12 = VectorXd::Zero(fem.nelx);
  for(int i = 0; i < fem.nelx; i++)
    factor12[i] = factor1[i].dot(factor2[i]);

  /// make factor 3
  vector<VectorXd> dualFes;
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    MatrixXd Ne = MatrixXd::Zero(fem.ndim, element[nel].numdof);
    MatrixXd X = MatrixXd::Zero(element[nel].ne, fem.ndim);
    MatrixXd U = MatrixXd::Zero(element[nel].ne, fem.ndim);
    VectorXd Uarr = VectorXd::Zero(element[nel].ne * fem.ndim);
    double jac = 1.0;
    VectorXd Fe = VectorXd::Zero(element[nel].numdof);
    BmatrixTotal bTotal;

    for(int i = 0; i < element[nel].ne; i++)
      for(int j = 0; j < fem.ndim; j++)
      {
        X.coeffRef(i, j) = node[element[nel].nodeID[i]].x[j];
        U.coeffRef(i, j) = node[element[nel].nodeID[i]].val[j];
        Uarr.coeffRef(fem.ndim * i + j) = node[element[nel].nodeID[i]].val[j];
      }

    for(int ip = 0; ip < element[nel].ipmax; ip++)
    {
      MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
      MatrixXd dFdu = MatrixXd::Zero(fem.ndim, fem.ndim);
      MatrixXd BeL = MatrixXd::Zero(fem.voigt, element[nel].numdof);
      MatrixXd dBLdu = MatrixXd::Zero(fem.voigt, element[nel].numdof);
      MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, element[nel].numdof);
      MatrixXd DdBdUU = MatrixXd::Zero(fem.voigt, element[nel].numdof);
      double dJdu;
      for(int k = 0; k < element[nel].ne * fem.ndim; k++)
      {
        VectorXd Uone = VectorXd::Zero(element[nel].ne * fem.ndim);
        Uone[k] = 1.0;
        bTotal.make(BeL, BeNL, Ne, jac, F, element[nel].eType,
                    element[nel].ipmax, X, U, ip);
        bTotal.makeDifferential(dBLdu, dFdu, dJdu, element[nel].eType,
                                element[nel].ipmax, X, U, ip, k);

        VectorXd dSdu = element[nel].De() * (BeL * Uone + dBLdu * Uarr);
        VectorXd S = element[nel].De() * BeL * Uarr;
        MatrixXd tmp =
            dFdu * Tensor(S, 1) * F.transpose() / F.determinant() +
            F * Tensor(dSdu, 1) * F.transpose() / F.determinant() +
            F * Tensor(S, 1) * dFdu.transpose() / F.determinant() -
            F * Tensor(S, 1) * F.transpose() * dJdu / pow(F.determinant(), 2.0);

        DdBdUU.col(k) = Voigt(tmp, 1);
      }

      Fe.transpose() += (pow(element[nel].design_s(), fem.qp) *
                         element[nel].young2() / element[nel].young()) *
                        factor1[nel].transpose() * DdBdUU / element[nel].ipmax;
    }
    dualFes.push_back(Fe);
  }

  VectorXd dualProbLoad =
      elementToNumeq<FiniteSaintVenantTotal, Node>(fem, element, node, dualFes);

  std::cout << "---- computing adjoint vector ----" << endl;
  VectorXd lambda = VectorXd::Zero(fem.numeq);
  Solvers solver(K, dualProbLoad, lambda);
  solver.solve("Eigen_CG", 0);
  std::cout << "------- finish to compute -------" << endl;

  vector<VectorXd> lambdaEle =
      numeqToElement<FiniteSaintVenantTotal, Node>(fem, element, node, lambda);

  // make factor 4 and combine with factor3
  vector<double> dgds(fem.nelx);
  VectorXd factor34 = VectorXd::Zero(fem.nelx);
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    int counter = 0;
    VectorXd dis = VectorXd::Zero(fem.ndim * element[nel].ne);
    for(int i = 0; i < element[nel].ne; i++)
    {
      for(int j = 0; j < fem.ndim; j++)
      {
        dis[counter] = node[element[nel].nodeID[i]].val[j];
        counter++;
      }
    }
    double dvpds = element[nel].pp() *
                   pow(element[nel].design_s(), element[nel].pp() - 1.0) *
                   (element[nel].young2() - element[nel].young1());

    VectorXd factor4 = (dvpds / element[nel].young()) * element[nel].fe();
    factor34[nel] = factor4.dot(lambdaEle[nel]);

    dgds[nel] = factor12[nel] - factor34[nel];
  } /*element loop end*/

  return dgds;
}
} // namespace stress_nonlinear
} // namespace icarat