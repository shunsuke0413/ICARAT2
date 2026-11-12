#pragma once
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include "filter.hpp"
#include "neo_hooke_damage.hpp"
#include "newton_raphson.hpp"
#include <io/boundary_condition_icarat.hpp>
#include <io/export_graph.hpp>
#include <io/export_vtk.hpp>
#include <io/file_reader.hpp>
#include <io/mesh_icarat.hpp>
#include <optimization/opt-solver/MMASolver.hpp>
#include <optimization/opt_base.hpp>

namespace icarat
{
namespace opt_finite_strain
{
// Converged states, including the undeformed state at index zero.
// tauInput[n] drives D[n] and is evaluated from displacement[n-1].
struct DamageFrame
{
  Eigen::MatrixXd displacement;
  Eigen::VectorXd damage;
  Eigen::VectorXd tauInput;
};
using DamageHistory = std::vector<DamageFrame>;

Eigen::SparseMatrix<double> damageFilterWeights(
    FEM &fem, const std::vector<FiniteNeoHookedamage> &element,
    const std::vector<Node> &node, double radius);
void initializeDamagePath(FEM &fem, std::vector<FiniteNeoHookedamage> &element,
                          std::vector<Node> &node);
DamageFrame captureDamageFrame(
    FEM &fem, const std::vector<FiniteNeoHookedamage> &element,
    const std::vector<Node> &node, const std::vector<double> &tauInput = {});
double damageWork(FEM &fem, std::vector<FiniteNeoHookedamage> &element,
                  const std::vector<Node> &node, const DamageHistory &history);
void sensitivity(FEM &fem, std::vector<FiniteNeoHookedamage> &element,
                 const std::vector<Node> &node, Optimize &optim, double V0,
                 const std::vector<double> &fil_gamma,
                 const DamageHistory &history,
                 const Eigen::SparseMatrix<double> &weights);
void fdm(FEM &fem, std::vector<FiniteNeoHookedamage> &element,
         std::vector<Node> &node, Force &force, ValueNR &dis, Optimize &optim,
         Filter<FiniteNeoHookedamage, Node> &filter, int numstep,
         const Eigen::SparseMatrix<double> &weights, double beta, double V0);
// Local derivative, filter, and coupled/reduced-system checks; no files written.
void verifyDamageDerivatives(FEM &fem,
                            std::vector<FiniteNeoHookedamage> &element,
                            const std::vector<Node> &node,
                            const DamageHistory &history,
                            const Eigen::SparseMatrix<double> &weights,
                            const std::vector<double> &fil_gamma);
} // namespace opt_finite_strain
} // namespace icarat