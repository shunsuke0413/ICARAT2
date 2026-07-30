#pragma once

#include "opt_differentiate_then_discretize.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;
namespace icarat {
namespace opt_differentiate_then_discretize {

template <class Obj> class sensitivity {
public:
  sensitivity(FEMDyna &fem, std::vector<DynamicElastic> &element,
              std::vector<Node> &node, Force &force, ValueDyna &dis,
              Optimize &optim, const Eigen::SparseMatrix<double> &globalK,
              const Eigen::SparseMatrix<double> &globalM,
              const history &history, Obj &function)
      : fem_(fem), element_(element), node_(node), force_(force), dis_(dis),
        optim_(optim), history_(history), function_(function),
        globalK_(globalK), globalM_(globalM),
        globalC_(-fem_.aa * globalK + -fem.bb * globalM),
        adj_disp(fem_.timestep + 1, Eigen::VectorXd(globalK.cols())),
        adj_velo(fem_.timestep + 1, Eigen::VectorXd(globalK.cols())),
        adj_acce(fem_.timestep + 1, Eigen::VectorXd(globalK.cols())) {

    A1 = 1.0 / (fem_.beta * pow(fem_.deltaT, 2.0));
    A2 = fem_.gamma / (fem_.beta * fem_.deltaT);
    A3 = (fem_.gamma / fem_.beta) - 1.0;
    A4 = 1.0 / (fem_.beta * fem_.deltaT);
    A5 = ((0.5 * fem_.gamma) / fem_.beta - 1.0) * fem_.deltaT;
    A6 = (0.5 / fem_.beta) - 1.0;
  };

  void adjointAnalysis();

private:
  // common params
  FEMDyna &fem_;
  std::vector<DynamicElastic> &element_;
  std::vector<Node> &node_;
  Force &force_;
  ValueDyna &dis_;
  Optimize &optim_;
  const history &history_;
  Obj &function_;

  // coeff for new mark beta
  double A1, A2, A3, A4, A5, A6;

  // for adjoint analysis
  const Eigen::SparseMatrix<double> &globalK_, &globalM_, globalC_;
  std::vector<Eigen::VectorXd> adj_disp, adj_velo, adj_acce;
};

template <class Obj> void sensitivity<Obj>::adjointAnalysis() {
  std::cout << "---------- computing adjoint vectors ----------" << endl;

  // terminal adjoint equation
  int k = fem_.timestep;
  adj_disp[k].setZero();
  adj_velo[k] = function_.getTerminalAdjointDot();
  adj_acce[k].setZero();
  adj_velo[k - 1] = -A3 * adj_velo[k];
  adj_acce[k - 1] = -A4 * adj_velo[k];
  
  // middle steps adjoint equation
  for (int t = 1; t < fem_.timestep + 1; t++) {
    k = fem_.timestep - t;

    Eigen::SparseMatrix<double> globalMK =
        A1 * globalM_ + globalK_ + A2 * globalC_;
    Eigen::VectorXd rhs =
        function_.getFext(k) +
        globalM_ * (A1 * adj_disp[k + 1] + A4 * adj_velo[k + 1] +
                    A6 * adj_acce[k + 1]) +
        globalC_ * (A2 * adj_disp[k + 1] + A3 * adj_velo[k + 1] +
                    A5 * adj_acce[k + 1]);
    Solvers solve(globalMK, rhs, adj_disp[k]);
    solve.solve("Eigen_LDLT", 0);
    auto incre = adj_disp[k] - adj_disp[k + 1];
    adj_velo[k] = A2 * incre - A3 * adj_velo[k + 1] - A5 * adj_acce[k + 1];
    adj_acce[k] = A1 * incre - A4 * adj_velo[k + 1] - A6 * adj_acce[k + 1];
  }

  // element loop start
  for (auto &e : element_) {
    double dfds = 0.0;
    int id = e.ID;
    double design = e.design_s();
    double young1 = e.young1();
    double young2 = e.young2();
    double p = e.pp();
    double young = (1.0 - pow(design, p)) * young1 + pow(design, p) * young2;
    double dyoung = p * pow(design, p - 1.0) * (young2 - young1);
    double E = dyoung / young;
    double D = 1.0 / e.density();

    // time integration start
    for (auto t = 0; t <= fem_.timestep; t++) {

      Eigen::VectorXd disp_e(e.numdof), velo_e(e.numdof), accel_e(e.numdof);
      Eigen::VectorXd lambda_e(e.numdof);

      // extract history values
      int counter = 0;
      for (int np = 0; np < e.ne; np++) {
        for (int j = 0; j < fem_.ndim; j++) {
          if (node_[e.nodeID[np]].dof[j] < fem_.numeq) {
            disp_e[counter] = history_.disp_h[t][node_[e.nodeID[np]].dof[j]];
            velo_e[counter] = history_.velo_h[t][node_[e.nodeID[np]].dof[j]];
            accel_e[counter] = history_.acce_h[t][node_[e.nodeID[np]].dof[j]];
            lambda_e[counter] = adj_disp[t][node_[e.nodeID[np]].dof[j]];
          } else {
            disp_e[counter] = 0.0;
            velo_e[counter] = 0.0;
            accel_e[counter] = 0.0;
            lambda_e[counter] = 0.0;
          }
          counter++;
        }
      }

      // sensitivity computation
      auto inner = D * e.Me() * accel_e +
                   (D * fem_.bb * e.Me() + E * fem_.aa * e.Ke()) * velo_e +
                   E * e.Ke() * disp_e;
      dfds += function_.pfps(e, E, disp_e, velo_e, accel_e) * fem_.deltaT +
              lambda_e.dot(inner) * fem_.deltaT;
    } // time loop end
    optim_.dfds[id] = dfds;
  } // element loop end
}
} // namespace opt_differentiate_then_discretize
} // namespace icarat
