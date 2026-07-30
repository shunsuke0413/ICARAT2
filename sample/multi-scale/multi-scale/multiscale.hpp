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
#include <problem/elasticity.hpp>

namespace icarat
{
namespace multiscale
{
struct Micro
{
  FEM fem;
  std::vector<mHomoElastic> element;
  std::vector<HomoNode> node;
};

struct Macro
{
  FEM fem;
  std::vector<MHomoElastic> element;
  std::vector<Node> node;
};

} // namespace multiscale
} // namespace icarat