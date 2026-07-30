#pragma once
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include <analysis/newton_raphson.hpp>
#include <io/boundary_condition_icarat.hpp>
#include <io/export_graph.hpp>
#include <io/export_vtk.hpp>
#include <io/file_reader.hpp>
#include <io/mesh_icarat.hpp>
#include <optimization/filter_helm.hpp>
#include <optimization/opt-solver/MMASolver.hpp>
#include <optimization/opt_base.hpp>
#include <problem/thermal_stress.hpp>

namespace icarat
{
namespace opt_multi_material
{
typedef std::vector<std::vector<std::vector<double>>> ijk;
typedef std::vector<std::vector<double>> ij;

void fdm(FEM &fem, std::vector<LinearThermalStressElastMultiMat> &element,
         std::vector<Node> &node, Force &force, Value &dis, Optimize &optim,
         FilterHelmholtz<LinearThermalStressElastMultiMat, Node> &filter,
         int numMat, std::vector<int> &nodeID, int direction, double deltaT);
} // namespace opt_multi_material
} // namespace icarat