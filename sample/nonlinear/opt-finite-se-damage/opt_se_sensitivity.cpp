#include "neo_hooke_damage.hpp"
#include <Eigen/Sparse>
#include <analysis/newton_raphson.hpp>
#include <basis/bmatrix.hpp>
#include <basis/bmatrix_total.hpp>
#include <basis/bmatrix_updated.hpp>
#include <misc/assembling.hpp>
#include <optimization/opt_base.hpp>
#include <problem/elasticity.hpp>

namespace icarat
{
namespace opt_finite_strain
{
using namespace std;
using namespace Eigen;

void calCFJ(FEM &fem, MatrixXd &C, MatrixXd &C_inv, MatrixXd &F,
            MatrixXd &F_inv, double &J)
{
  // right C-G tensor
  C = F.transpose() * F;
  J = F.determinant();

  // cal inverse of C
  if(fem.ndim == 2)
  {
    C_inv = C.inverse();
    F_inv = F.inverse();
  }
  else // 3D
  {
    C_inv = C.inverse();
    F_inv = F.inverse();
  }
}

void sensitivity(FEM &fem, vector<FiniteNeoHookedamage> &element,
                 vector<Node> &node, Value &value,
                 const Eigen::SparseMatrix<double> &K, Optimize &optim,
                 double V0, std::vector<double> &fil_gamma, double poisson)
{
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// 1. cal dSigmadu
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  vector<VectorXd> eledSigmadu(fem.nelx);

  Gauss9Square gauss;
  if(fem.ndim == 3)
  {
    Gauss8Cubic gauss;
  }

  std::vector<double> dPdgamma_vec(fem.nelx, 0.0);
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    double dens = element[nel].design_s();
    double young1 = element[nel].young1();
    double young2 = element[nel].young2();
    double p = element[nel].pp();

    // SIMP
    double young = (1.0 - pow(dens, p)) * young1 + pow(dens, p) * young2;

    double mu = 0.5 * young / (1 + poisson);
    double lambda =
        (poisson * young) / ((1.0 + poisson) * (1.0 - 2.0 * poisson));

    eledSigmadu[nel].resize(element[nel].numdof);

    MatrixXd X = MatrixXd::Zero(element[nel].ne, fem.ndim);
    MatrixXd U = MatrixXd::Zero(element[nel].ne, fem.ndim);
    VectorXd disEle = VectorXd::Zero(element[nel].numdof);

    for(int i = 0; i < element[nel].ne; i++)
    {
      for(int j = 0; j < fem.ndim; j++)
      {
        X.coeffRef(i, j) = node[element[nel].nodeID[i]].x[j];
        U.coeffRef(i, j) =
            element[nel].gamma() * node[element[nel].nodeID[i]].val[j];
        disEle.coeffRef(fem.ndim * i + j) = node[element[nel].nodeID[i]].val[j];
      }
    }

    VectorXd dPdU = VectorXd::Zero(element[nel].numdof);

    for(int s = 0; s < element[nel].numdof; s++)
    {
      BmatrixTotal bTotal;
      double dSigmadu = 0.0;

      for(int ip = 0; ip < element[nel].ipmax; ip++)
      {
        VectorXd r = gauss.points(ip);
        double w = gauss.weights(ip);
        double jac = 1.0;
        MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * element[nel].ne);
        MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
        MatrixXd F_inv = MatrixXd::Zero(fem.ndim, fem.ndim);
        MatrixXd dFdu = MatrixXd::Zero(fem.ndim, fem.ndim);
        MatrixXd BeL = MatrixXd::Zero(fem.voigt, element[nel].numdof);
        MatrixXd dBLdu = MatrixXd::Zero(fem.voigt, element[nel].numdof);
        MatrixXd BeNL =
            MatrixXd::Zero(fem.ndim * fem.ndim, element[nel].numdof);
        double dJdu = 0.0;

        bTotal.make(BeL, BeNL, Ne, jac, F, element[nel].eType,
                    element[nel].ipmax, X, U, ip);
        bTotal.makeDifferential(dBLdu, dFdu, dJdu, element[nel].eType,
                                element[nel].ipmax, X, U, ip, s);
        assert(jac > 0.0);

        MatrixXd I = MatrixXd::Identity(fem.ndim, fem.ndim);
        MatrixXd C = MatrixXd::Zero(fem.ndim, fem.ndim);
        MatrixXd C_inv = MatrixXd::Zero(fem.ndim, fem.ndim);
        MatrixXd D = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
        MatrixXd E = MatrixXd::Zero(fem.ndim * fem.ndim, fem.ndim * fem.ndim);
        double J = 0.0;
        calCFJ(fem, C, C_inv, F, F_inv, J);

        MatrixXd E_mat = 0.5 * (C - I);
        MatrixXd S_mat = mu * (I - C_inv) + lambda * J * (J - 1.0) * C_inv;

        double dSdutmp = 0.0;
        for(int m = 0; m < fem.ndim; m++)
        {
          for(int n = 0; n < fem.ndim; n++)
          {
            double dSE_dFmn = 0.0;
            for(int k = 0; k < fem.ndim; k++)
            {
              for(int l = 0; l < fem.ndim; l++)
              {
                double dE_dF =
                    0.5 * ((l == n ? F(m, k) : 0.0) + (k == n ? F(m, l) : 0.0));

                double dCinv_dF = -C_inv(k, n) * F_inv.transpose()(m, l) -
                                  F_inv.transpose()(m, k) * C_inv(n, l);

                double dS_dF = (lambda * J * (J - 1.0) - mu) * dCinv_dF +
                               lambda * (2.0 * J - 1.0) * J *
                                   F_inv.transpose()(m, n) * C_inv(k, l);

                dSE_dFmn += dS_dF * E_mat(k, l) + S_mat(k, l) * dE_dF;
              }
            }
            dSdutmp += dSE_dFmn * dFdu(m, n) * w * jac;
          }
        }
        dSigmadu += dSdutmp;

        /*for(int k = 0; k < fem.ndim; k++)
        {
          for(int l = 0; l < fem.ndim; l++)
          {
            for(int m = 0; m < fem.ndim; m++)
            {
              for(int n = 0; n < fem.ndim; n++)
              {
                int idx_S = k * fem.ndim + l;
                int idx_F = m * fem.ndim + n;

                double first_term = 0.0;
                double dEdu = 0.0;
                for(int o = 0; o < fem.ndim; o++)
                {
                  for(int p = 0; p < fem.ndim; p++)
                  {
                    first_term +=
                        -C_inv(k, o) * (F(o, p) + F(o, p)) * C_inv(p, l);
                    dEdu += (o == m ? F(n, p) : 0) + (n == p ? F(m, o) : 0);
                  }
                }
                first_term *= (-mu + lambda * J * (J - 1));
                dEdu *= 0.5;

                double second_term = lambda * (2 * J - 1) * J *
                                     F_inv.transpose()(m, n) * C_inv(k, l);

                D(idx_S, idx_F) = first_term / element[nel].ipmax +
                                  second_term / element[nel].ipmax;

                E(idx_S, idx_F) = dEdu / element[nel].ipmax;
              }
            }
          }
        }

        double dSdutmp = 0.0;
        for(int k = 0; k < fem.ndim; k++)
        {
          for(int l = 0; l < fem.ndim; l++)
          {
            int idx_S = k * fem.ndim + l;
            for(int m = 0; m < fem.ndim; m++)
            {
              for(int n = 0; n < fem.ndim; n++)
              {
                int idx_F = m * fem.ndim + n;
                dSdutmp += D(idx_S, idx_F) * dFdu(m, n) * 0.5 *
                               (C(m, n) - I(m, n)) * w * jac +
                           E(idx_S, idx_F) * dFdu(m, n) *
                               (mu * (I(m, n) - C_inv(m, n)) +
                                lambda * J * (J - 1) * C_inv(m, n)) *
                               w * jac;
              }
            }
          }
        }

        dSigmadu += dSdutmp;*/
      }

      dPdU(s) = dSigmadu;
      eledSigmadu[nel].coeffRef(s) = element[nel].gamma() * dSigmadu;
      // eledSigmadu[nel].coeffRef(s) = dSigmadu;
    }
    dPdgamma_vec[nel] = dPdU.dot(disEle);
  }

  VectorXd dSigmadu_nodetmp = elementToNumeq<FiniteNeoHookedamage, Node>(
      fem, element, node, eledSigmadu);
  VectorXd dSigmadu_node = dSigmadu_nodetmp;

  cout << "---- computing adjoint vector ----" << endl;
  Solvers solver(K, dSigmadu_node, value.val);
  solver.solve("Eigen_CG", 0);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// 2. cal dSigmadrho
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  VectorXd dSigmadrho(fem.nelx);
  for(int nelj = 0; nelj < fem.nelx; nelj++)
  {
    double dens = element[nelj].design_s();
    double young1 = element[nelj].young1();
    double young2 = element[nelj].young2();
    double p = element[nelj].pp();
    double young = (1.0 - pow(dens, p)) * young1 + pow(dens, p) * young2;
    double dyoung = p * pow(dens, p - 1.0) * (young2 - young1);
    double E = dyoung / young;

    MatrixXd X = MatrixXd::Zero(element[nelj].ne, fem.ndim);
    MatrixXd U = MatrixXd::Zero(element[nelj].ne, fem.ndim);
    for(int i = 0; i < element[nelj].ne; i++)
    {
      for(int j = 0; j < fem.ndim; j++)
      {
        X(i, j) = node[element[nelj].nodeID[i]].x[j];
        U(i, j) = element[nelj].gamma() * node[element[nelj].nodeID[i]].val[j];
      }
    }

    double S_E = 0.0;
    for(int ip = 0; ip < element[nelj].ipmax; ip++)
    {
      VectorXd r = gauss.points(ip);
      double w = gauss.weights(ip);
      double jac = 1.0;
      MatrixXd Ne = MatrixXd::Zero(fem.ndim, element[nelj].numdof);
      MatrixXd F = MatrixXd::Zero(fem.ndim, fem.ndim);
      MatrixXd F_inv = MatrixXd::Zero(fem.ndim, fem.ndim);
      MatrixXd BeL = MatrixXd::Zero(fem.voigt, element[nelj].numdof);
      MatrixXd BeNL = MatrixXd::Zero(fem.ndim * fem.ndim, element[nelj].numdof);

      BmatrixTotal bTotal;
      bTotal.make(BeL, BeNL, Ne, jac, F, element[nelj].eType,
                  element[nelj].ipmax, X, U, ip);
      assert(jac > 0.0);

      MatrixXd I = MatrixXd::Identity(fem.ndim, fem.ndim);
      MatrixXd C = MatrixXd::Zero(fem.ndim, fem.ndim);
      MatrixXd C_inv = MatrixXd::Zero(fem.ndim, fem.ndim);
      double J = 0.0;
      calCFJ(fem, C, C_inv, F, F_inv, J);

      MatrixXd E_mat = 0.5 * (C - I);
      S_E += element[nelj].product_tempv[ip] * w * jac;
    }
    // dSigmadrho(nelj) = E * S_E;
    dSigmadrho(nelj) = E * S_E + dPdgamma_vec[nelj] * fil_gamma[nelj];
  }

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    double dens = element[nel].design_s();
    double young1 = element[nel].young1();
    double young2 = element[nel].young2();
    double p = element[nel].pp();
    double young = (1.0 - pow(dens, p)) * young1 + pow(dens, p) * young2;
    double dyoung = p * pow(dens, p - 1.0) * (young2 - young1);
    double E = dyoung / young;

    VectorXd valEle = VectorXd::Zero(element[nel].numdof);
    VectorXd disEle = VectorXd::Zero(element[nel].numdof);
    for(int i = 0; i < element[nel].ne; i++)
    {
      for(int j = 0; j < fem.ndim; j++)
      {
        int dof = node[element[nel].nodeID[i]].dof[j];
        if(dof < fem.numeq)
        {
          valEle(fem.ndim * i + j) = value.val[dof];
        }
        else
        {
          valEle(fem.ndim * i + j) =
              node[element[nel].nodeID[i]].dirich->val[j];
          if(node[element[nel].nodeID[i]].dirich->flag[j] == 1)
          {
            valEle(fem.ndim * i + j) = 0.0;
          }
          disEle(fem.ndim * i + j) = node[element[nel].nodeID[i]].val[j];
        }
      }
    }

    VectorXd fe_rho = element[nel].gamma() * element[nel].nfe() +
                      (1.0 - element[nel].gamma() * element[nel].gamma()) *
                          element[nel].lfe();
    VectorXd fe_gamma = element[nel].nfe() +
                        element[nel].nKe() * disEle * element[nel].gamma() -
                        2.0 * element[nel].gamma() * element[nel].lfe();

    double indirect =
        valEle.transpose() * (E * fe_rho + fe_gamma * fil_gamma[nel]);

    optim.dfds[nel] = -dSigmadrho(nel) + indirect;
    optim.dgds[nel] = element[nel].volume / V0;
  }
}
} // namespace opt_finite_strain
} // namespace icarat