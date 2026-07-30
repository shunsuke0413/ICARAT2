#pragma once
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include <analysis/general_eigen.hpp>
#include <analysis/linear.hpp>
#include <io/boundary_condition_icarat.hpp>
#include <io/export_graph.hpp>
#include <io/export_vtk.hpp>
#include <io/file_reader.hpp>
#include <io/mesh_icarat.hpp>
#include <io/mesh_mdpa.hpp>
#include <misc/assembling.hpp>
#include <optimization/filter_helm.hpp>
#include <optimization/opt-solver/MMASolver.hpp>
#include <optimization/opt_base.hpp>
#include <problem/elasticity.hpp>

namespace icarat
{
namespace opt_buckling
{
struct FEMBuckling : FEM
{
  int numeigen;                        ///< number of eigenvalue
  double limeigen;                     ///< Lower limit of buckling load
  double pNorm;                        ///< p-norm paramteter(odd value)
  double sigma;                        ///< parameter for shift-invert solver
  double paraconv;                     ///< parameter for shift-invert solver
  std::string solverSt = "Eigen_LDLT"; ///< solver for structural analysis
  std::string solver = "Eigen_CG";     ///< solver for filtering
};

std::vector<double> sens_compliance(FEMBuckling &fem,
                                    std::vector<BucklingLinearElastic> &element,
                                    std::vector<Node> &node);

std::vector<double> sens_volume(FEMBuckling &fem,
                                std::vector<BucklingLinearElastic> &element,
                                double V0);

double func_buckling(FEMBuckling &fem, ValueEigen &value);

std::vector<double>
sens_buckling(FEMBuckling &fem, std::vector<BucklingLinearElastic> &element,
              std::vector<Node> &node,
              LinearAnalysis<BucklingLinearElastic, Node> &linear,
              ValueEigen &value);

void fdm(FEMBuckling &fem, std::vector<BucklingLinearElastic> &element,
         std::vector<Node> &node, Force &force, ValueEigen &dis,
         Optimize &optim, FilterHelmholtz<BucklingLinearElastic, Node> &filter);

} // namespace opt_buckling
} // namespace icarat
