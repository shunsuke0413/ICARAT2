#include <cmath>
#include <gtest/gtest.h>
#include <iostream>
#include <optimization/opt-solver/GCMMASolver.hpp>
#include <optimization/opt-solver/MMASolver.hpp>
#include <optimization/opt-solver/OCSolver.hpp>

using namespace std;
using namespace icarat;

//*********************************************************
//  Minimize    x0^2 + x1^2 + x2^2
//  Subject to  (x0 - 5)^2 + (x1 - 2)^2 + (x2 - 1)^2 <= 9
//              (x0 - 3)^2 + (x1 - 4)^2 + (x2 - 3)^2 <= 9
//              0 <= xj <= 5
//*********************************************************
TEST(MMASolverTest, toy)
{
  // number of design variable &constraint
  int n = 3, m = 2;
  MMASolver mma(n, m);
  mma.SetMoveLimit(0.5);
  mma.SetOscillation(0.1);
  mma.SetAsymptotes(0.5, 0.7, 1.2);
  vector<double> s({4.0, 3.0, 2.0});
  vector<double> smin = vector<double>(n, 0.0);
  vector<double> smax = vector<double>(n, 5.0);
  double objective;                   // Function value of compliance
  vector<double> dobjective(n);       // Sensitivities of compliance
  vector<double> constraints(m);      // Function values of weight
  vector<double> dconstraints(m * n); // Sensitivities of weight

  //----------Optimize loop----------
  for(int k = 0; k < 10; k++)
  {
    // cout << endl << "k = " << k << "\t";

    //  Objective
    objective = pow(s[0], 2.0) + pow(s[1], 2.0) + pow(s[2], 2.0);
    dobjective[0] = 2.0 * s[0];
    dobjective[1] = 2.0 * s[1];
    dobjective[2] = 2.0 * s[2];

    //  Constraint1
    constraints[0] = pow(s[0] - 5.0, 2.0) + pow(s[1] - 2.0, 2.0) +
                     pow(s[2] - 1.0, 2.0) - 9.0;
    dconstraints[0] = 2.0 * (s[0] - 5.0);
    dconstraints[2] = 2.0 * (s[1] - 2.0);
    dconstraints[4] = 2.0 * (s[2] - 1.0);

    //  Constraint2
    constraints[1] = pow(s[0] - 3.0, 2.0) + pow(s[1] - 4.0, 2.0) +
                     pow(s[2] - 3.0, 2.0) - 9.0;
    dconstraints[1] = 2.0 * (s[0] - 3.0);
    dconstraints[3] = 2.0 * (s[1] - 4.0);
    dconstraints[5] = 2.0 * (s[2] - 3.0);

    // cout << "f:\t" << objective << "\t";
    // cout << "g1:\t" << constraints[0] << "\t";
    // cout << "g2:\t" << constraints[1] << "\t";

    // for(auto si : s)
    // {
    //   cout << "\t" << si;
    // }
    //----------Update s----------
    mma.Update(s.data(), dobjective.data(), constraints.data(),
               dconstraints.data(), smin.data(), smax.data());
  }
  EXPECT_LT(objective, 8.78);
  EXPECT_LT(abs(constraints[0]), 1.0e-6);
  EXPECT_LT(abs(constraints[1]), 1.0e-6);
}

TEST(GCMMASolverTest, toy)
{
  // number of design variable & constraint
  int n = 3, m = 2;
  GCMMASolver gcmma(n, m);
  gcmma.SetMoveLimit(0.5);
  gcmma.SetOscillation(0.1);
  gcmma.SetAsymptotes(0.5, 0.7, 1.2);
  vector<double> s = vector<double>({4.0, 3.0, 2.0});
  vector<double> smin = vector<double>(n, 0.0);
  vector<double> smax = vector<double>(n, 5.0);
  double objective;                   // Function value of compliance
  vector<double> dobjective(n);       // Sensitivities of compliance
  vector<double> constraints(m);      // Function values of weight
  vector<double> dconstraints(m * n); // Sensitivities of weight

  //----------Optimize loop----------
  for(int k = 0; k < 10; k++)
  {
    //  cout << endl << "k = " << k << "\t";

    //  Objective
    objective = pow(s[0], 2.0) + pow(s[1], 2.0) + pow(s[2], 2.0);
    dobjective[0] = 2.0 * s[0];
    dobjective[1] = 2.0 * s[1];
    dobjective[2] = 2.0 * s[2];

    //  Constraint1
    constraints[0] = pow(s[0] - 5.0, 2.0) + pow(s[1] - 2.0, 2.0) +
                     pow(s[2] - 1.0, 2.0) - 9.0;
    dconstraints[0] = 2.0 * (s[0] - 5.0);
    dconstraints[2] = 2.0 * (s[1] - 2.0);
    dconstraints[4] = 2.0 * (s[2] - 1.0);

    //  Constraint2
    constraints[1] = pow(s[0] - 3.0, 2.0) + pow(s[1] - 4.0, 2.0) +
                     pow(s[2] - 3.0, 2.0) - 9.0;
    dconstraints[1] = 2.0 * (s[0] - 3.0);
    dconstraints[3] = 2.0 * (s[1] - 4.0);
    dconstraints[5] = 2.0 * (s[2] - 3.0);

    //----------Update s----------
    vector<double> snew(n); // new design variable
    double objectivenew =
        pow(snew[0], 2.0) + pow(snew[1], 2.0) + pow(snew[2], 2.0);
    vector<double> constraintsnew(2);

    gcmma.OuterUpdate(snew.data(), s.data(), objective, dobjective.data(),
                      constraints.data(), dconstraints.data(), smin.data(),
                      smax.data());

    constraintsnew[0] = pow(snew[0] - 5.0, 2.0) + pow(snew[1] - 2.0, 2.0) +
                        pow(snew[2] - 1.0, 2.0) - 9.0;
    constraintsnew[1] = pow(snew[0] - 3.0, 2.0) + pow(snew[1] - 4.0, 2.0) +
                        pow(snew[2] - 3.0, 2.0) - 9.0;

    bool conserv = gcmma.ConCheck(objectivenew, constraintsnew.data());
    for(int inneriter = 0; !conserv && inneriter < 15; ++inneriter)
    {
      // Inner iteration update
      gcmma.InnerUpdate(snew.data(), objectivenew, constraintsnew.data(),
                        s.data(), objective, dobjective.data(),
                        constraints.data(), dconstraints.data(), smin.data(),
                        smax.data());

      // Check conservativity
      objectivenew = pow(snew[0], 2.0) + pow(snew[1], 2.0) + pow(snew[2], 2.0);
      constraintsnew[0] = pow(snew[0] - 5.0, 2.0) + pow(snew[1] - 2.0, 2.0) +
                          pow(snew[2] - 1.0, 2.0) - 9.0;
      constraintsnew[1] = pow(snew[0] - 3.0, 2.0) + pow(snew[1] - 4.0, 2.0) +
                          pow(snew[2] - 3.0, 2.0) - 9.0;
      conserv = gcmma.ConCheck(objectivenew, constraintsnew.data());
    }
    s = snew;
    // cout << "f:\t" << objective << "\t";
    // cout << "g1:\t" << constraints[0] << "\t";
    // cout << "g2:\t" << constraints[1] << "\t";

    // for(auto si : s)
    // {
    //   cout << "\t" << si;
    // }
  }
  EXPECT_LT(objective, 8.78);
  EXPECT_LT(abs(constraints[0]), 1.0e-6);
  EXPECT_LT(abs(constraints[1]), 1.0e-6);
}

/// @todo OCSolver's test code