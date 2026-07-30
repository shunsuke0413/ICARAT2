#pragma once
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include <analysis/newton_raphson.hpp>
#include <analysis/solver/solvers.hpp>
#include <basis/bmatrix.hpp>
#include <io/boundary_condition_icarat.hpp>
#include <io/export_graph.hpp>
#include <io/export_vtk.hpp>
#include <io/file_reader.hpp>
#include <io/mesh_icarat.hpp>
#include <misc/assembling.hpp>
#include <optimization/filter.hpp>
#include <optimization/opt-solver/MMASolver.hpp>
#include <optimization/opt_base.hpp>
#include <problem/saint_venant.hpp>

namespace icarat
{
namespace stress_nonlinear
{

struct FEMStressConst : FEM
{
  double stressLim; ///< stress limit
  double pNorm;     /// p-norm parameter (4.....32 are often used)
  double qp;        /// parameter for qp-relaxation(cf.) Holmberg et al.
};

std::vector<double> sensCompliance(FEMStressConst &fem,
                                   std::vector<FiniteSaintVenantTotal> &element,
                                   std::vector<Node> &node, Force &force,
                                   Value &value,
                                   const Eigen::SparseMatrix<double> &K,
                                   std::vector<double> &fil_gamma);

std::vector<double>
sensStressConst(FEMStressConst &fem,
                std::vector<FiniteSaintVenantTotal> &element,
                std::vector<Node> &node, const Eigen::SparseMatrix<double> &K);

double funcStressConst(FEMStressConst &fem,
                       std::vector<FiniteSaintVenantTotal> &element);

void fdm(FEMStressConst &fem, std::vector<FiniteSaintVenantTotal> &element,
         std::vector<Node> &node, Force &force, Value &value, Optimize &optim,
         Filter<FiniteSaintVenantTotal, Node> &filter, int numstep,
         double tolNR, int iterNR);

} // namespace stress_nonlinear
} // namespace icarat
