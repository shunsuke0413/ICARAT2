///  @file  OCSolver.hpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// doublehis software is released under the MIdouble License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <vector>

namespace icarat
{
class OC
{
public:
  OC(int _n) : n(_n) {}
  ~OC() {}

  /// F is lambda function for computing the volume (Please see sample.)
  template <class F>
  void UpdateVariables(std::vector<double> &_xk, std::vector<double> &_dfdx,
                       double _g, std::vector<double> &_dgdx,
                       std::vector<double> &xmin, std::vector<double> &xmax,
                       F _gkp1)
  {
    //----------Get updated design variables with OC method----------
    double lambda0 = this->lambdamin, lambda1 = this->lambdamax, lambda;
    std::vector<double> xkp1 = std::vector<double>(this->n);
    int counter = 0;
    while((lambda1 - lambda0) / (lambda1 + lambda0) > this->lambdaeps &&
          this->maxiter > counter)
    {
      lambda = 0.5 * (lambda1 + lambda0);

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int i = 0; i < this->n; i++)
      {
        xkp1[i] =
            pow(std::max(0.0, -_dfdx[i] / (_dgdx[i] * lambda)), this->iota) *
            _xk[i];
        if(xkp1[i] < std::max(xmin[i], (1.0 - this->movelimit) * _xk[i]))
        {
          xkp1[i] = std::max(xmin[i], (1.0 - this->movelimit) * _xk[i]);
        }
        else if(xkp1[i] > std::min(xmax[i], (1.0 + this->movelimit) * _xk[i]))
        {
          xkp1[i] = std::min(xmax[i], (1.0 + this->movelimit) * _xk[i]);
        }
      }
      if(_gkp1(xkp1) - _g > 0.0)
      {
        lambda0 = lambda;
      }
      else
      {
        lambda1 = lambda;
      }
      counter++;
      //   std::cout << "counter: " << counter << std::endl;
    }

    if(counter >= this->maxiter)
      std::cerr << "not converge in OCSolver.hpp" << std::endl;
    // std::cout << "counter: " << counter << std::endl;
    _xk = xkp1;
  }

  void setMoveLimit(double m) { movelimit = m; }

private:
  int n; //  Number of design variables
  double iota = 0.5;
  double lambdamin = 0.0;
  double lambdamax = 1.0e4;
  double lambdaeps = 1.0e-5;
  double movelimit = 0.15;
  int maxiter = 100;
};

} // namespace icarat
