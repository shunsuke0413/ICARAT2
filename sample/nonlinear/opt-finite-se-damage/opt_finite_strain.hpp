#pragma once
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include "neo_hooke_damage.hpp"
#include <analysis/newton_raphson.hpp>
#include <io/boundary_condition_icarat.hpp>
#include <io/export_graph.hpp>
#include <io/export_vtk.hpp>
#include <io/file_reader.hpp>
#include <io/mesh_icarat.hpp>
#include <optimization/filter.hpp>
#include <optimization/opt-solver/MMASolver.hpp>
#include <optimization/opt_base.hpp>

namespace icarat
{
namespace opt_finite_strain
{
void fdm(FEM &fem, std::vector<FiniteNeoHookedamage> &element,
         std::vector<Node> &node, Force &force, ValueNR &dis, Optimize &optim,
         Filter<FiniteNeoHookedamage, Node> &filter, int numstep);
extern void calCFJ(FEM &fem, Eigen::MatrixXd &C, Eigen::MatrixXd &C_inv,
                   Eigen::MatrixXd &F, Eigen::MatrixXd &F_inv, double &J);
extern double delta(int i, int j);
void sensitivity(FEM &fem, std::vector<FiniteNeoHookedamage> &element,
                 std::vector<Node> &node, Value &value,
                 const Eigen::SparseMatrix<double> &K, Optimize &optim,
                 double V0, std::vector<double> &fil_gamma, double poisson);
} // namespace opt_finite_strain
} // namespace icarat