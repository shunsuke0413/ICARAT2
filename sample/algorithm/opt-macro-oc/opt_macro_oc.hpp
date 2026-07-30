#pragma once
#include <chrono>
#include <cmath>
#include <iostream>
#include <toml/toml.hpp>
#include <vector>

#include <analysis/linear.hpp>
#include <io/boundary_condition_icarat.hpp>
#include <io/export_graph.hpp>
#include <io/export_vtk.hpp>
#include <io/mesh_icarat.hpp>
#include <io/mesh_mdpa.hpp>
#include <optimization/filter_helm.hpp>
#include <optimization/opt-solver/OCSolver.hpp>
#include <optimization/opt_base.hpp>
#include <problem/elasticity.hpp>