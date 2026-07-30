#pragma once
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include "subloading-analytical_iso.hpp"
#include "newton_raphson.hpp"
#include <basis/bmatrix.hpp>
#include <io/boundary_condition_icarat.hpp>
#include <io/export_graph.hpp>
#include <io/export_vtk.hpp>
#include <io/file_reader.hpp>
#include <io/mesh_icarat.hpp>
#include <io/mesh_vtk.hpp>
#include <optimization/filter.hpp>
#include <optimization/opt-solver/OCSolver.hpp>
#include <optimization/opt-solver/MMASolver.hpp>
#include <optimization/opt_base.hpp>

namespace icarat
{
namespace opt_elast_plast
{
class OptLinearElastPlast : public LinearElastPlast
{
public:
  OptLinearElastPlast(int voigt_, int ne_, int ipmax_, int numdof_)
      : LinearElastPlast(voigt_, ne_, ipmax_, numdof_), stress1_obj_(ipmax_),
        stress_obj_(ipmax_), strain_p_bar_obj_(ipmax_), R_obj_(ipmax_), R_old_(ipmax_), strain_old_(ipmax_),
        stress_old_(ipmax_)
  {
    for(int i = 0; i < ipmax_; i++)
    {
      stress1_obj_[i].resize(6);
      stress1_obj_[i].setZero();

      stress_obj_[i].resize(6);
      stress_obj_[i].setZero();

      strain_old_[i].resize(voigt_);
      strain_old_[i].setZero();

      stress_old_[i].resize(6);
      stress_old_[i].setZero();
    }
  }

  void init_opt()
  {
    func_new_ = 0.0;
    objective_ = 0.0;

    for(int i = 0; i < ipmax; i++)
    {
      stress1_obj_[i].setZero();
      stress_obj_[i].setZero();
      strain_p_bar_obj_[i] = 0.0;
      R_obj_[i] = 0.0;
      R_old_[i] = 0.0;
      strain_old_[i].setZero();
      stress_old_[i].setZero();
    }
  }

  void sensitivity(FEM &fem, std::vector<Node> &node);

  const double &func_new() const { return func_new_; }
  const double &objective() const { return objective_; }

private:
  void sensitivity_plane_strain(FEM &fem, std::vector<Node> &node);
  void sensitivity_3d(FEM &fem, std::vector<Node> &node);

protected:
  double func_new_;
  double objective_;
  std::vector<Eigen::VectorXd> stress1_obj_;
  std::vector<Eigen::VectorXd> stress_obj_;
  std::vector<double> strain_p_bar_obj_;
  std::vector<double> R_obj_;
  std::vector<double> R_old_;
  std::vector<Eigen::VectorXd> strain_old_; ///< previous strain[ipmax][voigt]
  std::vector<Eigen::VectorXd> stress_old_; ///< previous stress[ipmax][voigt]
};

void fdm(FEM &fem, std::vector<OptLinearElastPlast> &element,
         std::vector<Node> &node, Force &force, ValueNR &dis, Optimize &optim,
         int numstep, double filterR);
// void taylor_test(FEM &fem, std::vector<OptLinearElastPlast> &element,
//                  std::vector<Node> &node, Force &force, ValueNR &dis, Optimize &optim,
//                  int numstep, double filterR);

} // namespace opt_elast_plast
} // namespace icarat
