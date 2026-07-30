///  @file  explicit_dyna.hpp
///  @author  Daiki Watanabe
///  @date  May 14, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php
#pragma once
#include "implicit_dyna.hpp"
namespace icarat
{
/// dynamic analysis class in explicit method
template <class E, class N> class ExplicitDyna : public ImplicitDyna<E, N>
{
public:
  ExplicitDyna(FEM &fem, std::vector<E> &element, std::vector<N> &node,
               Force &force, ValueDyna &value);

  ///   main function for structural analysis
  void solve(int thisStep, int totalStep, double deltaT);
};

template <class E, class N>
ExplicitDyna<E, N>::ExplicitDyna(FEM &fem, std::vector<E> &element,
                                 std::vector<N> &node, Force &force,
                                 ValueDyna &value)
    : ImplicitDyna<E, N>(fem, element, node, force, value)
{
}

template <class E, class N>
void ExplicitDyna<E, N>::solve(int thisStep, int totalStep, double deltaT)
{
  this->assemblingMK("explicit", log);

  // update fint
  this->force.fint = this->force.fext - this->globalK * this->value.val;
  // update accel
  this->value.accel = this->globalM * this->force.fint;
  // update velocity
  this->value.velocity = this->value.velocity + deltaT * this->value.accel;
  // update val
  this->value.val = this->value.val + deltaT * this->value.velocity;
  // value.du_incre = deltaT * value.velocity;

  divideValIntoNode(this->fem.dofnp, this->fem.numeq, this->value.val,
                    this->node);
}

} // namespace icarat
