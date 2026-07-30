#include "opt_discretize_then_differentiate_finite.hpp"

namespace icarat {
namespace opt_discretize_then_differentiate_finite {
using namespace std;
using namespace Eigen;

template <class Obj> class Sensitivity {
public:
  Sensitivity(FEMDyna &fem, std::vector<DynamicNeoHookeTotal> &element,
              std::vector<Node> &node, Force &force, ValueDyna &dis,
              Optimize &optim, Obj &function,
              const Eigen::SparseMatrix<double> &globalM)
      : fem_(fem), element_(element), node_(node), force_(force), dis_(dis),
        optim_(optim), function_(function), Ru(fem.numeq, fem.numeq),
        R1u(fem.numeq, fem.numeq), R1v(fem.numeq, fem.numeq),
        R1a(fem.numeq, fem.numeq), globalM_(globalM),
        globalK_h(fem.timestep + 1,
                  Eigen::SparseMatrix<double>(fem.numeq, fem.numeq)),
        disp_h(fem.timestep + 1, Eigen::VectorXd(fem.numeq)),
        velo_h(fem.timestep + 1, Eigen::VectorXd(fem.numeq)),
        acce_h(fem.timestep + 1, Eigen::VectorXd(fem.numeq)),
        fint_L_h(fem.timestep + 1,
                 std::vector<Eigen::VectorXd>(
                     fem.nelx, Eigen::VectorXd(element_[0].numdof))),
        fint_NL_h(fem.timestep + 1,
                  std::vector<Eigen::VectorXd>(
                      fem.nelx, Eigen::VectorXd(element_[0].numdof))),
        ke_h(fem.timestep + 1,
             std::vector<Eigen::MatrixXd>(
                 fem.nelx,
                 Eigen::MatrixXd(element_[0].numdof, element_[0].numdof))),
        adj_lambda(fem.timestep + 1, Eigen::VectorXd(fem.numeq)),
        adj_mu(fem.timestep + 1, Eigen::VectorXd(fem.numeq)),
        adj_eta(fem.timestep + 1, Eigen::VectorXd(fem.numeq)) {

    // newmark params
    A1 = 1.0 / (fem_.beta * fem_.deltaT * fem_.deltaT);
    A2 = fem_.gamma / (fem_.beta * fem_.deltaT);
    A3 = (fem_.gamma / fem_.beta) - 1.0;
    A4 = 1.0 / (fem_.beta * fem_.deltaT);
    A5 = (((0.5 * fem_.gamma) / fem_.beta) - 1.0) * fem_.deltaT;
    A6 = (0.5 / fem_.beta) - 1.0;

    Hu = A1;
    H1u = -1.0 * A1;
    H1v = -1.0 * A4;
    H1a = -1.0 * A6;

    Lu = A2;
    L1u = -1.0 * A2;
    L1v = -1.0 * A3;
    L1a = -1.0 * A5;
  }

  // adjoint problem
  void adjointAnalysis(std::vector<double> &fil_gamma);

  FEMDyna &fem_;
  std::vector<DynamicNeoHookeTotal> &element_;
  std::vector<Node> &node_;
  Force &force_;
  ValueDyna &dis_;
  Optimize &optim_;
  Obj &function_;

  // global const
  Eigen::SparseMatrix<double> Ru, R1u, R1v, R1a;
  const Eigen::SparseMatrix<double> &globalM_;

  // global wise time dependent
  std::vector<Eigen::SparseMatrix<double>> globalK_h;
  std::vector<Eigen::VectorXd> disp_h;
  std::vector<Eigen::VectorXd> velo_h;
  std::vector<Eigen::VectorXd> acce_h;

  // element wise time dependent
  std::vector<std::vector<Eigen::VectorXd>> fint_L_h;
  std::vector<std::vector<Eigen::VectorXd>> fint_NL_h;
  std::vector<std::vector<Eigen::MatrixXd>> ke_h;

  // adjoint vectors
  std::vector<Eigen::VectorXd> adj_lambda;
  std::vector<Eigen::VectorXd> adj_mu;
  std::vector<Eigen::VectorXd> adj_eta;

  // const variables
  double A1, A2, A3, A4, A5, A6;
  double Hu, H1u, H1v, H1a;
  double Lu, L1u, L1v, L1a;
};

template <class Obj>
void Sensitivity<Obj>::adjointAnalysis(std::vector<double> &fil_gamma) {
  std::cout << "Path Dependent Sensitivity Analysis" << std::endl;

  // const value
  R1u = -1.0 * A1 * globalM_;
  R1v = -1.0 * A4 * globalM_;
  R1a = -1.0 * A6 * globalM_;

  // terminal adjoint equation
  std::printf("\rSTEP No. %3d", fem_.timestep);
  std::fflush(stdout);
  Ru = A1 * globalM_ + globalK_h[fem_.timestep];
  Eigen::VectorXd rhs = -function_.pfpu(disp_h[fem_.timestep]);
  adj_eta[fem_.timestep].setZero();
  adj_mu[fem_.timestep].setZero();
  Solvers solv0(Ru, rhs, adj_lambda[fem_.timestep]);
  solv0.solve("Eigen_CG", 1);


  // inverse time integration
  for (int t = 1; t <= fem_.timestep; t++) {
    int k = fem_.timestep - t;
    std::printf("\rSTEP No. %3d", k);
    std::fflush(stdout);

    Ru = A1 * globalM_ + globalK_h[k];

    adj_mu[k] =
        R1a * adj_lambda[k + 1] + H1a * adj_mu[k + 1] + L1a * adj_eta[k + 1];
    adj_eta[k] =
        R1v * adj_lambda[k + 1] + H1v * adj_mu[k + 1] + L1v * adj_eta[k + 1];
    rhs = -function_.pfpu(disp_h[k]) - R1u * adj_lambda[k + 1] -
          Hu * adj_mu[k] - H1u * adj_mu[k + 1] - Lu * adj_eta[k] -
          L1u * adj_eta[k + 1];
    Solvers solv(Ru, rhs, adj_lambda[k]);
    solv.solve("Eigen_CG", 0);
  } // time step loop end
#ifdef _OPENMP
#pragma omp parallel for default(none), shared(optim, node)
#endif
  for (auto &e : element_) {
    int id = e.ID;
    double design = e.design_s();

    double young1 = e.young1();
    double young2 = e.young2();
    double p = e.pp();
    double young =
        young1 + (design / (1.0 + p * (1.0 - design))) * (young2 - young1);
    double dyoung = (1 + p) * (young2 - young1) / (1 + p * (1 - design)) /
                    (1 + p * (1 - design));
    double E = dyoung / young;

    double density2 = e.density();
    double density = design * density2;
    double ddensity = density2;
    double D = ddensity / density;

    // gathered values
    Eigen::VectorXd disp_e(e.numdof), accel_e(e.numdof), lambda_e(e.numdof);

    optim_.dfds[id] = 0.0;

    // time integration start
    for (int t = 0; t <= fem_.timestep; t++) {
      // extract value from global equation
      int counter = 0;
      for (int np = 0; np < e.ne; np++) {
        for (int j = 0; j < fem_.ndim; j++) {
          if (node_[e.nodeID[np]].dof[j] < fem_.numeq) {
            disp_e[counter] = disp_h[t][node_[e.nodeID[np]].dof[j]];
            accel_e[counter] = acce_h[t][node_[e.nodeID[np]].dof[j]];
            lambda_e[counter] = adj_lambda[t][node_[e.nodeID[np]].dof[j]];
          } else {
            disp_e[counter] = 0.0;
            accel_e[counter] = 0.0;
            lambda_e[counter] = 0.0;
          }
          counter++;
        }
      }
      // extract element information
      auto fint_rho = (e.gamma() * fint_NL_h[t][id] +
                       (1.0 - e.gamma() * e.gamma()) * fint_L_h[t][id]);
      auto fint_gamma = fint_NL_h[t][id] + ke_h[t][id] * disp_e * e.gamma() -
                        2.0 * fint_L_h[t][id] * e.gamma();
      // compute dfdx
      optim_.dfds[id] += lambda_e.dot(e.Me() * D * accel_e + fint_rho * E +
                                      fint_gamma * fil_gamma[id]);
    } // time loop end
    optim_.dgds[id] = e.volume;
  std::cout << optim_.dfds[id] << std::endl;
  }
}

} // namespace opt_discretize_then_differentiate_finite
} // namespace icarat
