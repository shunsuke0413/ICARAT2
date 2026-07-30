#pragma once
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include <analysis/linear.hpp>
#include <analysis/solver/solvers.hpp>
#include <basis/bmatrix.hpp>
#include <io/boundary_condition_icarat.hpp>
#include <io/export_graph.hpp>
#include <io/export_vtk.hpp>
#include <io/file_reader.hpp>
#include <io/mesh_icarat.hpp>
#include <misc/assembling.hpp>
#include <optimization/filter.hpp>
#include <optimization/opt-solver/GCMMASolver.hpp>
#include <optimization/opt_base.hpp>
#include <problem/elasticity.hpp>

namespace icarat
{
namespace stress
{

struct FEMStressConst : FEM
{
  double stressLim; ///< stress limit
  double pNorm;     /// p-norm parameter (4.....32 are often used)
  double qp;        /// parameter for qp-relaxation(cf.) Holmberg et al.
};

std::vector<double> sensCompliance(FEMStressConst &fem,
                                   std::vector<LinearElastic> &element,
                                   std::vector<Node> &node);

std::vector<double> sensStressConst(FEMStressConst &fem,
                                    std::vector<LinearElastic> &element,
                                    std::vector<Node> &node,
                                    const Eigen::SparseMatrix<double> &K);

double calConstraintValue(FEMStressConst &fem,
                          std::vector<LinearElastic> &element);

void fdm(FEMStressConst &fem, std::vector<LinearElastic> &element,
         std::vector<Node> &node, Force &force, Value &value, Optimize &optim,
         Filter<LinearElastic, Node> &filter);

} // namespace stress
} // namespace icarat
