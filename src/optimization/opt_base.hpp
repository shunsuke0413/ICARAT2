///  @file  opt_base.hpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <cmath>
#include <iostream>
#include <vector>
namespace icarat
{
class Optimize
{
public:
  /// this constructor needs the number of some properties
  /// The arguments are "design: number of design variable
  /// constraint: number of constraint function
  /// xmin: minimum value of design variable
  /// xmax: maximum value of design variable
  Optimize(int design, int constraint, double xmin, double xmax)
      : design_s(design), const_h(constraint), dfds(design),
        dgds(constraint * design), xmin(design, xmin), xmax(design, xmax)
  {
  }

  int isConvergence(double epsilon, double minstep)
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
  };

  // for MPI interface
  int isConvergence(int procid, double epsilon, double minstep)
  {
    // convergence judge
    if(abs(object_old - object_f) / object_f < epsilon && optstep > minstep)
    {
      if(procid == 0)
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
  };

  std::vector<double> design_s; ///< design variable
  std::vector<double> const_h;  ///< constraint value
  std::vector<double> dfds;     ///< sensitivity of objective function
  std::vector<double> dgds;     ///< sensitivity of constraint function
  std::vector<double> xmin;     ///< minimum value of design variable
  std::vector<double> xmax;     ///< maximum value of design variable
  int optstep;                  ///< optimization step
  double object_f;              ///< new objective function

private:
  double object_old; ///< old objective function
};
} // namespace icarat