#pragma once
#include <cmath>
#include <iostream>
#include <vector>

#include <analysis/implicit_dyna.hpp>
#include <io/boundary_condition_icarat.hpp>
#include <io/export_graph.hpp>
#include <io/export_vtk.hpp>
#include <io/file_reader.hpp>
#include <io/mesh_icarat.hpp>
#include <optimization/filter.hpp>
#include <optimization/opt-solver/MMASolver.hpp>
#include <optimization/opt_base.hpp>
#include <problem/elasticity.hpp>

namespace icarat
{
namespace opt_discretize_then_differentiate
{
struct FEMDyna : FEM
{
  int timestep;
  double deltaT;
  double gamma;
  double beta;
  double aa;
  double bb;
};

namespace objective_function
{
// objective function class for minimizing time integration of squared norm of
// displacement vector
class disp_squared_norm
{
public:
  disp_squared_norm(int numeq_, double deltaT_, Eigen::VectorXd &fext)
      : ID(numeq_), deltaT(deltaT_)
  {
    ID.setZero();
    for(int i = 0; i < numeq_; i++)
    {
      if(std::fabs(fext[i]) > 1.0e-7)
        ID[i] = 1.0;
    }
  }
  void init() { value = 0.0; }
  void update(Eigen::VectorXd &disp)
  {
    value += (ID.array() * disp.array()).matrix().squaredNorm() * deltaT;
  }
  double getValue() { return value; }

  auto pfpu(Eigen::VectorXd &disp)
  {
    return (2.0 * disp.array() * ID.array() * deltaT).matrix();
  }

private:
  Eigen::VectorXd ID; // has value of 0 or 1 to identify target dof
  double deltaT;
  double value;
};

// objective function class for minimizing dynamic compliance which integrate
// external work in time
class dynamic_compliance
{
};
} // namespace objective_function
void fdm(FEMDyna &fem, std::vector<DynamicElastic> &element,
         std::vector<Node> &node, Force &force, ValueDyna &dis, Optimize &optim,
         Filter<DynamicElastic, Node> &filter);

} // namespace opt_discretize_then_differentiate
} // namespace icarat
