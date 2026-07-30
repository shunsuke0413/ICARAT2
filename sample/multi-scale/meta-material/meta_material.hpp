#pragma once
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include <analysis/homogenization.hpp>
#include <analysis/linear.hpp>
#include <io/boundary_condition_icarat.hpp>
#include <io/export_graph.hpp>
#include <io/export_vtk.hpp>
#include <io/file_reader.hpp>
#include <io/mesh_icarat.hpp>
#include <io/mesh_mdpa.hpp>
#include <optimization/filter_helm.hpp>
#include <optimization/opt-solver/MMASolver.hpp>
#include <optimization/opt_base.hpp>
#include <problem/elasticity.hpp>

namespace icarat
{
namespace meta_material
{
double getObjectF(FEM &fem, Eigen::MatrixXd &tCH, Eigen::MatrixXd &CH,
                  Eigen::MatrixXd &omega);

void sensitivity(FEM &fem, std::vector<mHomoElastic> &element,
                 std::vector<HomoNode> &node, ValueHomo &value, Optimize &optim,
                 Eigen::MatrixXd &tCH, Eigen::MatrixXd &CH,
                 Eigen::MatrixXd &omega, double V0, double VTotal);

void fdm(const toml::value &config, FEM &fem,
         std::vector<mHomoElastic> &element, std::vector<HomoNode> &node,
         ValueHomo &dis, Optimize &optim,
         FilterHelmholtz<mHomoElastic, HomoNode> &filter,
         Eigen::MatrixXd &omega,
         std::vector<std::vector<std::pair<int, int>>> &pairs,
         std::array<double, 3> &length);

} // namespace meta_material
} // namespace icarat