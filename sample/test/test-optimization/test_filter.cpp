#include <cmath>
#include <gtest/gtest.h>
#include <iostream>
#include <optimization/filter.hpp>
#include <optimization/filter_helm.hpp>
#include <problem/base.hpp>

using namespace std;
using namespace icarat;

namespace
{
// make 4 element mesh with quad4 for structure analysis
void make4Box(FEM &fem, vector<Element> &element, vector<Node> &node)
{
  // fem
  fem.ndim = 2;
  fem.nelx = 4;
  fem.voigt = 3;
  fem.dofnp = 2;
  fem.numnp = 9;
  fem.neq = fem.numnp * fem.dofnp;
  fem.numeq = fem.numeq;

  // element
  int xdiv = 2, ydiv = 2;
  for(int i = 0; i < ydiv; i++)
  {
    for(int j = 0; j < xdiv; j++)
    {
      Element actele;
      actele.ID = i;
      actele.eType = "quad4";
      actele.ne = 4;
      actele.ipmax = 4;
      actele.volume = 1.0;
      actele.nodeID.resize(4);
      actele.nodeID[0] = xdiv * i + j + i;
      actele.nodeID[1] = xdiv * i + j + i + 1;
      actele.nodeID[2] = xdiv * i + j + i + 1 + (xdiv + 1);
      actele.nodeID[3] = xdiv * i + j + i + (xdiv + 1);
      element.push_back(actele);
    }
  }

  // node
  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 3; j++)
    {
      Node actnode;
      actnode.ID = 3 * i + j;
      actnode.dof[0] = 2 * i;
      actnode.dof[1] = 2 * i + 1;
      actnode.x[0] = 1.0 * j;
      actnode.x[1] = 1.0 * i;
      node.push_back(actnode);
    }
  }
}
} // namespace

TEST(FilterTest, filter)
{
  FEM fem;
  vector<Element> element;
  vector<Node> node;
  make4Box(fem, element, node);
  double radius = 2.0;
  Filter<Element, Node> filter(fem, element, node, radius);
  vector<double> s(4, 0.5), dfds(4), answer(4);

  // sensitivity filter check
  dfds[0] = 1.0e-4;
  dfds[1] = 1.0;
  dfds[2] = 1.0e-4;
  dfds[3] = 1.0;

  auto dfds_new = filter.sensitivityFilter(s, dfds);

  // answer sensitivity
  answer[0] = 0.34587;
  answer[1] = 0.65423;
  answer[2] = 0.34587;
  answer[3] = 0.65423;

  double norm = 0.0;
  for(int i = 0; i < 4; i++)
    norm += pow(dfds_new[i] - answer[i], 2.0);
  EXPECT_LT(norm, 1.0e-10);

  // density filter check
  s[0] = 1.0e-4;
  s[1] = 1.0;
  s[2] = 1.0e-4;
  s[3] = 1.0;

  auto s_new = filter.densityFilter(s, false);
  dfds_new = filter.densityFilter(dfds, true);

  norm = 0.0;
  for(int i = 0; i < 4; i++)
    norm += pow(s_new[i] - answer[i], 2.0);
  EXPECT_LT(norm, 1.0e-10);

  norm = 0.0;
  for(int i = 0; i < 4; i++)
    norm += pow(dfds_new[i] - answer[i], 2.0);
  EXPECT_LT(norm, 1.0e-10);

  // threshold function check
  // beta=big value is almost step function
  double beta = 16.0;
  auto s_thre = filter.threshold(s_new, beta, 0.5);
  EXPECT_LT(s_thre[0], 0.01);
  EXPECT_GT(s_thre[1], 0.99);

  // beta=1.0 is almost linear function
  auto der_s = filter.thresholdDerivative(s_new, 1.0, 0.5);
  norm = 0.0;
  for(int i = 0; i < 4; i++)
    norm += pow(der_s[i] - 1.0, 2.0);
  EXPECT_LT(norm, 0.015);
}

TEST(FilterHelmholtzTest, filter)
{
  FEM fem;
  vector<Element> element;
  vector<Node> node;
  make4Box(fem, element, node);
  double radius = 2.0;
  FilterHelmholtz<Element, Node> filter(fem, element, node, radius);
  vector<double> s(4, 0.5), answer(4);

  // answer sensitivity
  answer[0] = 0.312569;
  answer[1] = 0.687531;
  answer[2] = 0.312569;
  answer[3] = 0.687531;

  // density filter check
  s[0] = 1.0e-4;
  s[1] = 1.0;
  s[2] = 1.0e-4;
  s[3] = 1.0;

  auto s_new = filter.densityFilter(s, "Eigen_CG", 0, false);
  auto deriv = filter.densityFilter(s, "Eigen_CG", 0, true);

  double norm = 0.0;
  for(int i = 0; i < 4; i++)
    norm += pow(s_new[i] - answer[i], 2.0);
  EXPECT_LT(norm, 1.0e-10);

  norm = 0.0;
  for(int i = 0; i < 4; i++)
    norm += pow(deriv[i] - answer[i], 2.0);
  EXPECT_LT(norm, 1.0e-10);

  // threshold function check
  // beta=big value is almost step function
  double beta = 16.0;
  auto s_thre = filter.threshold(s_new, beta, 0.5);
  EXPECT_LT(s_thre[0], 0.01);
  EXPECT_GT(s_thre[1], 0.99);

  // beta=1.0 is almost linear function
  auto der_s = filter.thresholdDerivative(s_new, 1.0, 0.5);
  norm = 0.0;
  for(int i = 0; i < 4; i++)
    norm += pow(der_s[i] - 1.0, 2.0);

  EXPECT_LT(norm, 0.015);
}
