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
#include <problem/neo_hooke.hpp>

namespace icarat
{
namespace meta_material_finite
{
double getObjectF(FEM &fem, Eigen::MatrixXd &tCH, Eigen::MatrixXd &CH,
                  Eigen::MatrixXd &omega);

void sensitivity(FEM &fem, std::vector<mHomoNeoHookeTotal> &element,
                 std::vector<HomoNode> &node, ValueHomo &value, Optimize &optim,
                 Eigen::MatrixXd &tCH, Eigen::MatrixXd &CH,
                 Eigen::MatrixXd &omega, double V0, double VTotal, double fac);

void fdm(const toml::value &config, FEM &fem,
         std::vector<mHomoNeoHookeTotal> &element, std::vector<HomoNode> &node,
         ValueHomo &dis, Optimize &optim,
         FilterHelmholtz<mHomoNeoHookeTotal, HomoNode> &filter,
         Eigen::MatrixXd &omega,
         std::vector<std::vector<std::pair<int, int>>> &pairs,
         std::array<double, 3> &length, Eigen::MatrixXd &tCH,
         Eigen::VectorXd &Mstrain, double beta, double fac);

Eigen::MatrixXd isoDe(int voigt, double tYoung, double tPoisson,
                      std::string mattype);

Eigen::MatrixXd orthoDe(int voigt, std::array<double, 3> tYoung,
                        std::array<double, 6> tPoisson, std::string mattype);
} // namespace meta_material_finite
} // namespace icarat