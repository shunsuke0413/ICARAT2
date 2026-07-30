#include "opt_discretize_then_differentiate.hpp"

namespace icarat {
namespace opt_discretize_then_differentiate {
using namespace std;
using namespace Eigen;

template <class Obj> class Sensitivity {
public:
  Sensitivity(FEMDyna &fem, std::vector<DynamicElastic> &element,
              std::vector<Node> &node, Force &force, ValueDyna &dis,
              Optimize &optim, Obj &function,
              const Eigen::SparseMatrix<double> &globalK,
              const Eigen::SparseMatrix<double> &globalM)
      : fem_(fem), element_(element), node_(node), force_(force), dis_(dis),
        optim_(optim), function_(function), globalK_(globalK),
        globalM_(globalM), Ru(fem.numeq, fem.numeq), R1u(fem.numeq, fem.numeq),
        R1v(fem.numeq, fem.numeq), R1a(fem.numeq, fem.numeq),
        disp_h(fem.timestep + 1, Eigen::VectorXd(fem.numeq)),
        velo_h(fem.timestep + 1, Eigen::VectorXd(fem.numeq)),
        acce_h(fem.timestep + 1, Eigen::VectorXd(fem.numeq)),
        fint_h(fem.timestep + 1,
               std::vector<Eigen::VectorXd>(
                   fem.nelx, Eigen::VectorXd(element_[0].numdof))),
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
  void adjointAnalysis();

  FEMDyna &fem_;
  std::vector<DynamicElastic> &element_;
  std::vector<Node> &node_;
  Force &force_;
  ValueDyna &dis_;
  Optimize &optim_;
  Obj &function_;

  // global const
  const Eigen::SparseMatrix<double> &globalK_, &globalM_;

  // global wise time dependent
  Eigen::SparseMatrix<double> Ru, R1u, R1v, R1a;
  std::vector<Eigen::VectorXd> disp_h;
  std::vector<Eigen::VectorXd> velo_h;
  std::vector<Eigen::VectorXd> acce_h;
  // element wise time dependent
  std::vector<std::vector<Eigen::VectorXd>> fint_h;

  // adjoint vectors
  std::vector<Eigen::VectorXd> adj_lambda;
  std::vector<Eigen::VectorXd> adj_mu;
  std::vector<Eigen::VectorXd> adj_eta;

  // const variables
  double A1, A2, A3, A4, A5, A6;
  double Hu, H1u, H1v, H1a;
  double Lu, L1u, L1v, L1a;
};

template <class Obj> void Sensitivity<Obj>::adjointAnalysis() {
  std::cout << "Path Dependent Sensitivity Analysis" << std::endl;

  // const value
  R1u = -1.0 * A1 * globalM_;
  R1v = -1.0 * A4 * globalM_;
  R1a = -1.0 * A6 * globalM_;

  // terminal adjoint equation
  std::printf("\rSTEP No. %3d", fem_.timestep);
  std::fflush(stdout);
  Ru = A1 * globalM_ + globalK_;
  Eigen::VectorXd rhs = -function_.pfpu(disp_h[fem_.timestep]);
  adj_eta[fem_.timestep].setZero();
  adj_mu[fem_.timestep].setZero();
  Solvers solv0(globalK_, rhs, adj_lambda[fem_.timestep]);
  solv0.solve("Eigen_CG", 1);

  // inverse time integration
  for (int t = 1; t <= fem_.timestep; t++) {
    int k = fem_.timestep - t;
    std::printf("\rSTEP No. %3d", k);
    std::fflush(stdout);

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
    double young = (1.0 - pow(design, p)) * young1 + pow(design, p) * young2;
    double dyoung = p * pow(design, p - 1.0) * (young2 - young1);
    double E = dyoung / young;

    double density1 = 0.0; // void density
    double density2 = e.density();
    double density = design * density2;
    double ddensity = density2;
    double D = ddensity / density;

    // gathered values
    Eigen::VectorXd accel_e(e.numdof), lambda_e(e.numdof);
    Eigen::VectorXd fint_e(e.numdof);

    optim_.dfds[id] = 0.0;

    // time integration start
    for (int t = 0; t <= fem_.timestep; t++) {
      // extract element information
      fint_e = fint_h[t][id];
      // extract value from global equation
      int counter = 0;
      for (int np = 0; np < e.ne; np++) {
        for (int j = 0; j < fem_.ndim; j++) {
          if (node_[e.nodeID[np]].dof[j] < fem_.numeq) {
            accel_e[counter] = acce_h[t][node_[e.nodeID[np]].dof[j]];
            lambda_e[counter] = adj_lambda[t][node_[e.nodeID[np]].dof[j]];
          } else {
            accel_e[counter] = 0.0;
            lambda_e[counter] = 0.0;
          }
          counter++;
        }
      }
      // compute dfdx
      optim_.dfds[id] += lambda_e.dot(e.Me() * D * accel_e + fint_e * E);
    } // time loop end
    optim_.dgds[id] = e.volume;
  }
}

} // namespace opt_discretize_then_differentiate
} // namespace icarat
