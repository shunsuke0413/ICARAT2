#pragma once
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include <analysis/linear.hpp>
#include <io/boundary_condition_icarat.hpp>
#include <io/export_graph.hpp>
#include <io/export_vtk.hpp>
#include <io/file_reader.hpp>
#include <io/mesh_icarat.hpp>
#include <optimization/filter_helm.hpp>
#include <optimization/opt-solver/MMASolver.hpp>
#include <optimization/opt_base.hpp>
#include <problem/elasticity.hpp>

namespace icarat
{
namespace helmholtz
{
void fdm(FEM &fem, std::vector<LinearElastic> &element, std::vector<Node> &node,
         Force &force, Value &dis, Optimize &optim,
         FilterHelmholtz<LinearElastic, Node> &filter);
}
} // namespace icarat