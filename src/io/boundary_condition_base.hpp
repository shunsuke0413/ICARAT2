///  @file  boundary_condition_base.hpp
///  @author  Daiki Watanabe
///  @date  February 9, 2022.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <Eigen/Dense>
#include <io/mesh_icarat.hpp>
#include <numeric>
#include <problem/base.hpp>

namespace icarat
{

/// set Dirichlet object & value to node
template <class N>
void setDirichValue(int dofnp, N &node, std::vector<int> &flag,
                    std::vector<double> &val)
{
  if(node.dirich == nullptr)
  {
    std::shared_ptr<Dirichlet> dirich(new Dirichlet);
    node.dirich = dirich;
  }
  for(int k = 0; k < dofnp; k++)
  {
    if(node.dirich->flag[k] != 1)
      node.dirich->flag[k] = flag[k];
    node.dirich->val[k] = val[k];
  }
}

template <class N> int countDOF(FEM &fem, N &node)
{
  // initialization
  for(auto &n : node)
  {
    for(int j = 0; j < 6; j++)
      n.dof[j] = -1;
  }

  int counter = 0;
  for(auto &n : node)
  {
    // dirichlet is nothing
    if(n.dirich == nullptr)
    {
      for(int j = 0; j < fem.dofnp; j++)
      {
        n.dof[j] = counter;
        counter++;
      }
    }
    else
    {
      for(int j = 0; j < fem.dofnp; j++)
      {
        // dirichlet is not active
        if(n.dirich->flag[j] == 0)
        {
          n.dof[j] = counter;
          counter++;
        }
        // dirichlet is active... through
        else
          continue;
      }
    }
  }

  // constrainted dof
  int constCounter = counter;
  for(auto &n : node)
  {
    if(n.dirich != nullptr)
    {
      for(int j = 0; j < fem.dofnp; j++)
      {
        // dirichlet is active
        if(n.dirich->flag[j] != 0)
        {
          n.dof[j] = constCounter;
          constCounter++;
        }
      }
    }
  }

  assert(constCounter == fem.neq);

  std::cout << "numeq : " << counter << std::endl;
  return counter;
}

/// T is string
inline std::tuple<std::string, int> reducedOrderEtypeNe(std::string eType,
                                                        int ne)
{
  std::string eTyper;
  int ner;

  if(eType == "tria6")
  {
    eTyper = "tria3";
    ner = 3;
  }
  else if(eType == "quad8")
  {
    eTyper = "quad4";
    ner = 4;
  }
  else if(eType == "tetra10")
  {
    eTyper = "tria4";
    ner = 4;
  }
  else if(eType == "hexa20")
  {
    eTyper = "hexa8";
    ner = 8;
  }
  else
  {
    eTyper = eType;
    ner = ne;
  }

  return std::forward_as_tuple(eTyper, ner);
}

template <class E, class N> int countDOFFluid(FEM &fem, E &element, N &node)
{
  // initialization
  for(auto &n : node)
  {
    for(int j = 0; j < 6; j++)
      n.dof[j] = -1;
  }

  // count velocity part
  int counter = 0;
  for(auto &n : node)
  {
    // dirichlet is nothing
    if(n.dirich == nullptr)
    {
      for(int j = 0; j < fem.ndim; j++)
      {
        n.dof[j] = counter;
        counter++;
      }
    }
    else
    {
      for(int j = 0; j < fem.ndim; j++)
      {
        // dirichlet is not active
        if(n.dirich->flag[j] == 0)
        {
          n.dof[j] = counter;
          counter++;
        }
        // dirichlet is active... through
        else
          continue;
      }
    }
  }

  // count pressure part
  for(auto &e : element)
  {
    int nep;
    std::string eTypep;
    std::tie(eTypep, nep) = reducedOrderEtypeNe(e.eType, e.ne);
    for(int i = 0; i < nep; i++)
    {
      if(node[e.nodeID[i]].dirich == nullptr) // dirichlet is not set.
      {
        // When dof is not set
        if(node[e.nodeID[i]].dof[fem.ndim] == -1)
        {
          node[e.nodeID[i]].dof[fem.ndim] = counter;
          counter++;
        }
        else
          continue;
      }
      else // dirichlet is set.
      {
        if(node[e.nodeID[i]].dirich->flag[fem.ndim] == 0)
        {
          // When dof is not set
          if(node[e.nodeID[i]].dof[fem.ndim] == -1)
          {
            node[e.nodeID[i]].dof[fem.ndim] = counter;
            counter++;
          }
          else
            continue;
        }
        else // when value is preset.
          continue;
      }
    }
  }

  // constrainted dof
  int constCounter = counter;
  for(auto &n : node)
  {
    if(n.dirich != nullptr)
    {
      for(int j = 0; j < fem.dofnp; j++)
      {
        // dirichlet is active
        if(n.dirich->flag[j] != 0)
        {
          n.dof[j] = constCounter;
          constCounter++;
        }
      }
    }
  }

  assert(constCounter == fem.neq);

  std::cout << "numeq : " << counter << std::endl;
  return counter;
}

/// base class for constraint input
struct InputConstraint
{
  std::vector<int> flag;
  std::vector<double> val;
};

/// base class for load input
struct InputLoad
{
  /// true: distributed, false: concentrated load
  bool isDist = true;
  /// side ID of the element
  int sideID;
  /// load value. "moment" and "pressure" use this for one element
  std::vector<double> val;
  /// load type. Default is "xyz". "moment" and "pressure" are available.
  std::string type = "xyz";

  /// Reference axis origin, vector and sum of areadist where load is applied.
  /// These are used to compute "moment" loads.
  std::array<double, 3> axis_origin;
  std::array<double, 3> axis_vector;
  double sumareadist;
};

/// Class for creating distributed loads.
class SurfaceLoad
{
public:
  SurfaceLoad(FEM &fem, Element &actele) : fem_(fem), actele_(actele) {}

  /// Determine which line (surface) of element is subject to external loads.
  int judgeLine(int counter, std::vector<int> hitpoints);

  int setNumNode();

  int surfaceIpmax();

  ///  compute the equivalent nodal force at an element to Force.fext_org
  void generate(Force &force, Eigen::VectorXi &idof, Eigen::MatrixXd &X,
                InputLoad &inp_f);

private:
  /// Set and resize the parameters for the element's surface
  std::tuple<std::string, int, int>
  resizeSurfaceParameter(Eigen::VectorXi &idof, Eigen::MatrixXd &X,
                         std::string eType, int sideID);

  /// @brief make the pressure load value of the Neumann condition
  /// defined in local coordinates the value in global coordinates
  Eigen::VectorXd transLocalToGlobalPressure(Eigen::MatrixXd &X,
                                             Eigen::VectorXd &val, int iline);

  Eigen::VectorXd transLocalToGlobalMoment(Eigen::MatrixXd &X,
                                           Eigen::VectorXd &val,
                                           InputLoad &inp_f);

  FEM &fem_;
  Element &actele_;

  /// tria3,6
  /// 2
  /// |  |
  /// 5     4          y
  /// |        |      ↑
  /// 0----3-----1     →x
  ///
  /// side = 1...0-1
  /// side = 2...1-2
  /// side = 3...2-0

  /// tetra4,10
  /// 3
  /// |  9
  /// 7  8 2            z  y
  /// | 6    5      ↑↗
  /// 0---4----1     →x
  ///
  /// side = 1...0-1-2
  /// side = 2...0-1-3
  /// side = 3...1-2-3
  /// side = 4...2-0-3

  /// quad4,8
  /// 3-----6------2
  /// |            |
  /// 7            5       y
  /// |            |      ↑
  /// 0------4-----1       →x
  ///
  /// side = 1...0-1
  /// side = 2...1-2
  /// side = 3...2-3
  /// side = 4...3-0

  /// hexa8,20
  ///    7-----14------6
  ///  15           13 |
  /// 4-----12-----5  18
  /// | 19         |  |
  /// 16  3    10  17 2     z   y
  /// | 11         | 9      ↑↗
  /// 0------8-----1          →x
  ///
  /// side = 1...0-1-2-3
  /// side = 2...3-0-4-7
  /// side = 3...0-1-5-4
  /// side = 4...1-2-6-5
  /// side = 5...2-3-7-6
  /// side = 6...4-5-6-7
};

} // namespace icarat
