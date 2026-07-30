#include "opt_stress_constraint.hpp"
using namespace std;
using namespace Eigen;
namespace icarat
{
namespace stress
{
double calConstraintValue(FEMStressConst &fem, vector<LinearElastic> &element)
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
                              vector<LinearElastic> &element,
                              vector<Node> &node)
{
  assert(element[0].young2() > element[0].young1());

  vector<double> dfds(fem.nelx);
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    double dens = element[nel].design_s();
    double young1 = element[nel].young1();
    double young2 = element[nel].young2();
    double p = element[nel].pp();

    double young = (1.0 - pow(dens, p)) * young1 + pow(dens, p) * young2;
    double dyoung = p * pow(dens, p - 1.0) * (young2 - young1);

    double E = dyoung / young;

    Eigen::VectorXd disp = Eigen::VectorXd::Zero(element[nel].numdof);
    for(int i = 0; i < element[nel].ne; i++)
      for(int j = 0; j < fem.ndim; j++)
        disp.coeffRef(fem.ndim * i + j) = node[element[nel].nodeID[i]].val[j];

    dfds[nel] = -E * disp.transpose() * element[nel].Ke() * disp;
  }

  return dfds;
}

vector<double> sensStressConst(FEMStressConst &fem,
                               vector<LinearElastic> &element,
                               vector<Node> &node,
                               const SparseMatrix<double> &K)
{
  vector<double> dgds(fem.nelx);

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
  vector<VectorXd> factor1;
  vector<VectorXd> factor2;
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    VectorXd facEle1 = VectorXd::Zero(fem.voigt);
    if(fem.ndim == 2)
    {
      vector<double> dvmdsig;
      dvmdsig.push_back(
          (2.0 * element[nel].stress()[0] - element[nel].stress()[1]) /
          (2.0 * element[nel].mises()));
      dvmdsig.push_back(
          (2.0 * element[nel].stress()[1] - element[nel].stress()[0]) /
          (2.0 * element[nel].mises()));
      dvmdsig.push_back(3.0 * element[nel].stress()[2] / element[nel].mises());

      for(int i = 0; i < fem.voigt; i++)
        facEle1[i] = dpndvm * pow(misesQP[nel] / fem.stressLim, fem.pNorm - 1) *
                     dvmdsig[i] / fem.stressLim;

      factor1.push_back(facEle1);
    }
    else
    {
      cerr << "stressconstrainted 3D is not applied." << endl;
      exit(1);
    }

    /// make factor2 for stress constraint's sensitivity
    VectorXd facEle2 = VectorXd::Zero(fem.voigt);
    if(fem.ndim == 2)
    {
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
    else
    {
      cerr << "stressconstrainted 3D is not applied." << endl;
      exit(1);
    }
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
    MatrixXd Be = MatrixXd::Zero(fem.voigt, element[nel].numdof);
    MatrixXd X = MatrixXd::Zero(element[nel].ne, fem.ndim);
    double jac = 1.0;
    VectorXd Fe = VectorXd::Zero(element[nel].numdof);
    Bmatrix bmatrix;

    for(int i = 0; i < element[nel].ne; i++)
      for(int j = 0; j < fem.ndim; j++)
        X.coeffRef(i, j) = node[element[nel].nodeID[i]].x[j];

    for(int ip = 0; ip < element[nel].ipmax; ip++)
    {
      bmatrix.make(Be, Ne, jac, element[nel].eType, element[nel].ipmax, X, ip);
      Fe.transpose() += (pow(element[nel].design_s(), fem.qp) *
                         element[nel].young2() / element[nel].young()) *
                        factor1[nel].transpose() * element[nel].De() * Be /
                        element[nel].ipmax;
    }
    dualFes.push_back(Fe);
  }

  VectorXd dualProbLoad =
      elementToNumeq<LinearElastic, Node>(fem, element, node, dualFes);

  std::cout << "---- computing adjoint vector ----" << endl;
  VectorXd lambda = VectorXd::Zero(fem.numeq);
  Solvers solver(K, dualProbLoad, lambda);
  solver.solve("Eigen_CG", 0);
  std::cout << "------- finish to compute -------" << endl;

  vector<VectorXd> lambdaEle =
      numeqToElement<LinearElastic, Node>(fem, element, node, lambda);

  /// make factor 4 and combine with factor3
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

    VectorXd factor4 = (dvpds / element[nel].young()) * element[nel].Ke() * dis;
    factor34[nel] = factor4.dot(lambdaEle[nel]);

    dgds[nel] = factor12[nel] - factor34[nel];
  } /*element loop end*/

  return dgds;
}
} // namespace stress
} // namespace icarat