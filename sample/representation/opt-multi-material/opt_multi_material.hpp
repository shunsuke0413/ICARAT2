#pragma once
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include "elast_multi_material.hpp"
#include <analysis/linear.hpp>
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
namespace opt_multi_material
{
typedef std::vector<std::vector<std::vector<double>>> ijk;
typedef std::vector<std::vector<double>> ij;

void fdm(FEM &fem, std::vector<LinearElastMultiMat> &element,
         std::vector<Node> &node, Force &force, Value &dis, Optimize &optim,
         Filter<LinearElastMultiMat, Node> &filter, int numMat);
} // namespace opt_multi_material
} // namespace icarat