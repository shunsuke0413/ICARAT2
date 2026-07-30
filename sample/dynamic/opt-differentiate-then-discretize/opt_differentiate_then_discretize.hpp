#pragma once
#include <cmath>
#include <iostream>
#include <vector>

#include <analysis/implicit_dyna.hpp>
#include <io/boundary_condition_icarat.hpp>
#include <io/export_graph.hpp>
#include <io/export_vtk.hpp>
#include <io/file_reader.hpp>
#include <io/mesh_icarat.hpp>
#include <optimization/filter.hpp>
#include <optimization/opt-solver/MMASolver.hpp>
#include <optimization/opt_base.hpp>
#include <problem/elasticity.hpp>

namespace icarat {
namespace opt_differentiate_then_discretize {
struct FEMDyna : FEM {
  std::string solvName;
  int timestep;
  double deltaT;
  double gamma = 0.5; // newmark-beta parameter
  double beta = 0.25; // newmark-beta parameter
  double aa = 0.0;    // damping parameter
  double bb = 0.0;    // damping parameter
  Eigen::VectorXd ll; // coefficient vector used in the objective function
};

struct history {
  history(int timestep, int numeq)
      : disp_h(timestep + 1, Eigen::VectorXd(numeq)),
        velo_h(timestep + 1, Eigen::VectorXd(numeq)),
        acce_h(timestep + 1, Eigen::VectorXd(numeq)) {}
  // history variables
  std::vector<Eigen::VectorXd> disp_h;
  std::vector<Eigen::VectorXd> velo_h;
  std::vector<Eigen::VectorXd> acce_h;
};
 
// \int_0^T d \cdot Kd \mathrm{d}t
namespace objective_function {
class strain_energy {
public:
  strain_energy(double deltaT, const Eigen::SparseMatrix<double> &globalK,
                const history &history)
      : deltaT_(deltaT), history_(history), globalK_(globalK){};
  auto init() { val = 0.0; };
  auto update(int t) {
    val +=
        0.5 * history_.disp_h[t].dot(globalK_ * history_.disp_h[t]) * deltaT_;
  }
  auto getValue() { return val; }
  Eigen::VectorXd getFext(int t) { return -globalK_ * history_.disp_h[t]; }
  auto getTerminalAdjointDot() {
    return Eigen::VectorXd::Zero(globalK_.cols());
  }
  auto pfps([[maybe_unused]] DynamicElastic &e, [[maybe_unused]] double E,
            [[maybe_unused]] Eigen::VectorXd &disp,
            [[maybe_unused]] Eigen::VectorXd &velo,
            [[maybe_unused]] Eigen::VectorXd &accel) {
    return 0.5 * E * disp.dot(e.Ke() * disp);
  }

private:
  double val;
  double deltaT_;
  const history &history_;
  const Eigen::SparseMatrix<double> &globalK_;
};

// \int_0^T d_\mathrm{tar}^2 \mathrm{d}t
class disp_squared_norm {
public:
  disp_squared_norm(int numeq_, double deltaT_, Eigen::VectorXd &fext,
                    const history &history)
      : ID(numeq_), history_(history), deltaT(deltaT_) {
    ID.setZero();
    for (int i = 0; i < numeq_; i++) {
      if (std::fabs(fext[i]) > 1.0e-7)
        ID[i] = 1.0;
    }
  }
  auto init() { value = 0.0; }
  auto update(int t) {
    value += (ID.array() * history_.disp_h[t].array()).matrix().squaredNorm() *
             deltaT;
  }
  double getValue() { return value; }

  auto getFext(int t) {
    return -(2.0 * history_.disp_h[t].array() * ID.array()).matrix();
  }
  auto getTerminalAdjointDot() { return Eigen::VectorXd::Zero(ID.size()); }
  auto pfps([[maybe_unused]] DynamicElastic &e, [[maybe_unused]] double E,
            [[maybe_unused]] Eigen::VectorXd &disp,
            [[maybe_unused]] Eigen::VectorXd &velo,
            [[maybe_unused]] Eigen::VectorXd &accel) {
    return 0.0;
  }

private:
  Eigen::VectorXd ID; // has value of 0 or 1 to identify target dof
  const history &history_;
  double deltaT;
  double value;
};

// \int_0^T f^\mathrm{ext} \cdot d \mathrm{d}t
class compliance_norm {
public:
  compliance_norm(double deltaT, Eigen::VectorXd &fext_org,
                  const history &history, int timestep)
      : history_(history), fext_org_(fext_org), deltaT_(deltaT),
        timestep(timestep) {}
  auto init() { value = 0.0; }
  auto update(int t) {
    value += sin(2.0 * M_PI * t / timestep) *
             fext_org_.dot(history_.disp_h[t]) * deltaT_;
  }
  double getValue() { return value; }

  auto getFext([[maybe_unused]] int t) {
    return -fext_org_ * sin(2.0 * M_PI * t / timestep);
  }
  auto getTerminalAdjointDot() {
    return Eigen::VectorXd::Zero(fext_org_.size());
  }
  auto pfps([[maybe_unused]] DynamicElastic &e, [[maybe_unused]] double E,
            [[maybe_unused]] Eigen::VectorXd &disp,
            [[maybe_unused]] Eigen::VectorXd &velo,
            [[maybe_unused]] Eigen::VectorXd &accel) {
    return 0.0;
  }

private:
  const history &history_;
  Eigen::VectorXd fext_org_; // has value of 0 or 1 to identify target dof
  double deltaT_;
  double value;
  int timestep;
};

} // namespace objective_function
} // namespace opt_differentiate_then_discretize
} // namespace icarat
