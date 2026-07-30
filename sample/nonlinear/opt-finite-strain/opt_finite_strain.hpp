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
#include <optimization/filter.hpp>
#include <optimization/opt-solver/MMASolver.hpp>
#include <optimization/opt_base.hpp>
#include <problem/neo_hooke.hpp>

namespace icarat
{
namespace opt_finite_strain
{
void fdm(FEM &fem, std::vector<FiniteNeoHookeTotal> &element,
         std::vector<Node> &node, Force &force, ValueNR &dis, Optimize &optim,
         Filter<FiniteNeoHookeTotal, Node> &filter, int numstep);
}
} // namespace icarat