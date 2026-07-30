///  @file  elast_levelset.hpp
///  @author  Daiki Watanabe
///  @date  October 13, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <iostream>
#include <problem/elasticity.hpp>
namespace icarat
{
namespace opt_levelset
{

///   FEM with levelset TO variable list
struct FEMLS : FEM
{
  double C;       ///< parameter
  double lambda;  ///< adjoint variable
  double sum_tdn; ///< sum total of tdn
  double tau;     ///< regularization parameter
};

///   node variable list in levelset Topopt
struct NodeLS : Node
{
  std::vector<int> eleID; ///< element ID this node has
  double phi;             ///< levelset function
  double tdn;             ///< topological derivative at node
};

/// convergence check class
class Checker
{
public:
  int isConvergence(double object_f, int optstep, double epsilon, int minstep)
  {
    // convergence judge
    if(abs(object_old - object_f) / object_f < epsilon && optstep > minstep)
    {
      std::cout << std::endl
                << "--------------------Optimized--------------------"
                << std::endl;
      return 0;
    }
    else
    {
      object_old = object_f;
      return 1;
    }
  }

private:
  double object_old;
};

///   element variable list in levelset Topopt
class ElastLevelset : public Elastic
{
public:
  ElastLevelset(int voigt_, int ne_, int ipmax_, int numdof_)
      : Elastic(voigt_, ne_, ipmax_, numdof_)
  {
  }

  /// make element stiffness matrix & internal force(Ke_ and fe_)
  void makeKeFint(std::vector<Eigen::Triplet<double>> &tripletsK,
                  std::vector<Eigen::Triplet<double>> &tripletsF,
                  std::vector<Eigen::Triplet<double>> &tripletsR, FEM &fem,
                  std::vector<NodeLS> &nodes);

  /// make mises stress in this elmenet (mises_)
  void makeMisesStress(FEMLS &fem, std::vector<NodeLS> &node);
};

} // namespace opt_levelset
} // namespace icarat