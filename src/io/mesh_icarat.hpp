///  @file  mesh_icarat.hpp
///  @author  Daiki Watanabe
///  @date  July 25, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include "mesh_base.hpp"
#include <problem/base.hpp>
#include <toml/toml.hpp>

namespace icarat
{
///  set box mesh data to element & node
/// element...set ID(ID) and number of node(ne), nodeIDs(nodeID), dof in the
/// element(numdof) node...set ID(ID) and coordinates(x)
template <class E, class N> class MeshIcarat
{
public:
  MeshIcarat(FEM &fem, const toml::value &);

  ///  set FEM data from toml file
  void setParameter(std::string probType);

  ///  main method of generating box mesh
  void generate(std::vector<E> &element, std::vector<N> &node);

  /// use for periodic boundary condition
  std::vector<std::vector<std::pair<int, int>>>
  generateMasterSlavePairs(std::vector<N> &node);

  /// set design variable for homogenized TO
  std::vector<double> setDesignHomo(std::vector<E> &element,
                                    std::vector<N> &node, double design_s0,
                                    double deltaH, double radius);

  /// getters
  int numeleTypes() { return 1; } // only one element type is available
  const int &nelx(int i) const { return nelx_; }
  const int &ne(int i) const { return ne_; }
  const int &ipmax(int i) const { return ipmax_; }
  const std::array<int, 3> &div() const { return div_; }
  const std::array<double, 3> &length() const { return length_; }

private:
  // generate quad4/8 or hexa8/20  box mesh (Element->ID,nodeID  Node->ID,x)
  void meshBox(std::vector<E> &element, std::vector<N> &node);

  void modifyMesh(std::vector<E> &element, std::vector<N> &node);

  FEM &fem_;
  const toml::value &config_;

  int nelx_;                     ///< number of element
  int ne_;                       ///< number of node in an element
  int ipmax_;                    ///< number of gauss point
  std::array<int, 3> div_;       ///< Number of divisions
  std::array<double, 3> length_; ///< length of sides
};

////////////////////////////////////////////////////////////////////
//////////////////////// public methods ///////////////////////////
///////////////////////////////////////////////////////////////////
template <class E, class N>
MeshIcarat<E, N>::MeshIcarat(FEM &fem, const toml::value &config)
    : fem_(fem), config_(config)
{
}

template <class E, class N>
void MeshIcarat<E, N>::setParameter(std::string probType)
{
  const auto &mesh = toml::find(config_, "mesh");
  // dimension
  fem_.ndim = toml::find<int>(mesh, "dimension");

  // dofnp, voigt
  std::tie(fem_.dofnp, fem_.voigt) = probTypeToDofnpVoigt(fem_.ndim, probType);

  // ne,ipmax
  const auto eleType = toml::find<std::string>(mesh, "type");

  ne_ = eleTypeToNe(eleType);
  ipmax_ = eleTypeToIpmax(eleType);

  // box division & length
  div_ = toml::find<std::array<int, 3>>(mesh, "division");
  length_ = toml::find<std::array<double, 3>>(mesh, "length");
  assert((int)div_.size() == 3);
  assert((int)length_.size() == 3);

  // numnp,neq,nelx
  if(eleType == "quad4" || eleType == "mitc4")
  {
    fem_.numnp = (div_[0] + 1) * (div_[1] + 1);
    fem_.neq = fem_.numnp * fem_.dofnp;
    fem_.nelx = div_[0] * div_[1];
  }
  else if(eleType == "quad8")
  {
    fem_.numnp = (2 * div_[0] + 1) * (div_[1] + 1) + (div_[0] + 1) * div_[1];
    fem_.neq = fem_.numnp * fem_.dofnp;
    fem_.nelx = div_[0] * div_[1];
  }
  else if(eleType == "hexa8")
  {
    fem_.numnp = (div_[0] + 1) * (div_[1] + 1) * (div_[2] + 1);
    fem_.neq = fem_.numnp * fem_.dofnp;
    fem_.nelx = div_[0] * div_[1] * div_[2];
  }
  else if(eleType == "hexa20")
  {
    fem_.numnp =
        4 * div_[0] * div_[1] * div_[2] +
        3 * (div_[0] * div_[1] + div_[1] * div_[2] + div_[2] * div_[0]) +
        2 * (div_[0] + div_[1] + div_[2]) + 1;
    fem_.neq = fem_.numnp * fem_.dofnp;
    fem_.nelx = div_[0] * div_[1] * div_[2];
  }
  else
    std::cerr << "Failed to read ne in mesh_icarat.hpp" << std::endl;

  if(probType == "fluid")
  {
    if(eleType == "quad8")
      fem_.neq -= (div_[0] + 1) * div_[1] + (div_[1] + 1) * div_[0];
    else if(eleType == "hexa20")
      fem_.neq -= (div_[0] + 1) * (div_[1] + 1) * (div_[2] + 1);
    else
    {
      std::cerr << "fluid problem need 2nd-order element in mesh_icarat.hpp"
                << std::endl;
      exit(1);
    }
  }

  nelx_ = fem_.nelx;
}

template <class E, class N>
void MeshIcarat<E, N>::generate(std::vector<E> &element, std::vector<N> &node)
{
  meshBox(element, node);

  const auto &mesh_key = toml::find(config_, "mesh");
  if(mesh_key.contains("selection"))
    modifyMesh(element, node);

  std::cout << "dimension : " << fem_.ndim << std::endl;
  std::cout << "number of node : " << fem_.numnp << std::endl;
  std::cout << "number of element : " << fem_.nelx << std::endl;
  std::cout << "total DOFs : " << fem_.neq << std::endl;
  std::cout << "check point at mesh generate           --->  ok" << std::endl;
}

template <class E, class N>
std::vector<std::vector<std::pair<int, int>>>
MeshIcarat<E, N>::generateMasterSlavePairs(std::vector<N> &node)
{
  // give tolerance
  double tol = 1.0e-7;

  double l_xmin = 0.0 - tol;
  double u_xmin = 0.0 + tol;
  double l_xmax = length_[0] - tol;
  double u_xmax = length_[0] + tol;

  double l_ymin = 0.0 - tol;
  double u_ymin = 0.0 + tol;
  double l_ymax = length_[1] - tol;
  double u_ymax = length_[1] + tol;

  double l_zmin = 0.0 - tol;
  double u_zmin = 0.0 + tol;
  double l_zmax = length_[2] - tol;
  double u_zmax = length_[2] + tol;

  std::vector<std::vector<std::pair<int, int>>> pairs(fem_.ndim);

  // 2D
  if(fem_.ndim == 2)
  {
    int countX = 0, countY = 0;

    for(int ii = 0; ii < (int)node.size(); ii++)
    {
      // x
      if(l_xmin <= node[ii].x[0] && node[ii].x[0] <= u_xmin)
      {
        std::pair<int, int> tmp;
        tmp.first = ii;
        pairs[0].push_back(tmp);
      }
      if(l_xmax <= node[ii].x[0] && node[ii].x[0] <= u_xmax)
      {
        pairs[0][countX].second = ii;
        countX++;
      }
      // y
      if(l_ymin <= node[ii].x[1] && node[ii].x[1] <= u_ymin)
      {
        std::pair<int, int> tmp;
        tmp.first = ii;
        pairs[1].push_back(tmp);
      }
      if(l_ymax <= node[ii].x[1] && node[ii].x[1] <= u_ymax)
      {
        pairs[1][countY].second = ii;
        countY++;
      }
    }
  }
  // 3D
  else
  {
    int countX = 0, countY = 0, countZ = 0;

    for(int ii = 0; ii < (int)node.size(); ii++)
    {
      // x
      if(l_xmin <= node[ii].x[0] && node[ii].x[0] <= u_xmin)
      {
        std::pair<int, int> tmp;
        tmp.first = ii;
        pairs[0].push_back(tmp);
      }
      if(l_xmax <= node[ii].x[0] && node[ii].x[0] <= u_xmax)
      {
        pairs[0][countX].second = ii;
        countX++;
      }
      // y
      if(l_ymin <= node[ii].x[1] && node[ii].x[1] <= u_ymin)
      {
        std::pair<int, int> tmp;
        tmp.first = ii;
        pairs[1].push_back(tmp);
      }
      if(l_ymax <= node[ii].x[1] && node[ii].x[1] <= u_ymax)
      {
        pairs[1][countY].second = ii;
        countY++;
      }
      // z
      if(l_zmin <= node[ii].x[2] && node[ii].x[2] <= u_zmin)
      {
        std::pair<int, int> tmp;
        tmp.first = ii;
        pairs[2].push_back(tmp);
      }
      if(l_zmax <= node[ii].x[2] && node[ii].x[2] <= u_zmax)
      {
        pairs[2][countZ].second = ii;
        countZ++;
      }
    }
  }

  return pairs;
}

template <class E, class N>
std::vector<double>
MeshIcarat<E, N>::setDesignHomo(std::vector<E> &element, std::vector<N> &node,
                                double design_s0, double deltaH, double radius)
{
  std::vector<double> design_s(fem_.nelx);

  if(fem_.ndim == 2)
  {
    assert(div_[0] % 2 == 0 && div_[1] % 2 == 0);

    double volume = length_[0] * length_[1];
    double deltaL =
        deltaH * pow(radius, 2.0) * M_PI / (volume - pow(radius, 2.0) * M_PI);

    double high = (1.0 + deltaH) * design_s0;
    double low = (1.0 - deltaL) * design_s0;

    Eigen::Vector2d centerP;
    centerP << length_[0] / 2.0, length_[1] / 2.0;

    for(int nel = 0; nel < fem_.nelx; nel++)
    {
      Eigen::Vector2d center = Eigen::Vector2d::Zero(2);

      for(int i = 0; i < element[nel].ne; i++)
      {
        center[0] += node[element[nel].nodeID[i]].x[0] / element[nel].ne;
        center[1] += node[element[nel].nodeID[i]].x[1] / element[nel].ne;
      }

      if((center - centerP).norm() <= radius)
        design_s[nel] = low;
      else
        design_s[nel] = high;
    }
  }
  else
  {
    assert(div_[0] % 2 == 0 && div_[1] % 2 == 0 && div_[2] % 2 == 0);

    double volume = length_[0] * length_[1] * length_[2];
    double deltaL = deltaH * 4.0 / 3.0 * pow(radius, 3.0) * M_PI /
                    (volume - 4.0 / 3.0 * pow(radius, 3.0) * M_PI);

    double high = (1.0 + deltaH) * design_s0;
    double low = (1.0 - deltaL) * design_s0;

    Eigen::Vector3d centerP;
    centerP << length_[0] / 2.0, length_[1] / 2.0, length_[2] / 2.0;

    for(int nel = 0; nel < fem_.nelx; nel++)
    {
      Eigen::Vector3d center = Eigen::Vector3d::Zero(3);

      for(int i = 0; i < element[nel].ne; i++)
      {
        center[0] += node[element[nel].nodeID[i]].x[0] / element[nel].ne;
        center[1] += node[element[nel].nodeID[i]].x[1] / element[nel].ne;
        center[2] += node[element[nel].nodeID[i]].x[2] / element[nel].ne;
      }

      if((center - centerP).norm() <= radius)
        design_s[nel] = low;
      else
        design_s[nel] = high;
    }
  }

  return design_s;
}

////////////////////////////////////////////////////////////////////
//////////////////////// private methods ///////////////////////////
////////////////////////////////////////////////////////////////////

template <class E, class N>
void MeshIcarat<E, N>::meshBox(std::vector<E> &element, std::vector<N> &node)
{
  auto eleType = toml::find<std::string>(config_, "mesh", "type");

  if(eleType == "quad4")
  {
    // coordinate
    for(int i = 0; i < div_[1] + 1; i++)
    {
      for(int j = 0; j < div_[0] + 1; j++)
      {
        int tmp = j + i + div_[0] * i;
        node[tmp].ID = j + i + div_[0] * i;

        node[tmp].x[0] = ((double)length_[0] / div_[0]) * j;
        node[tmp].x[1] = ((double)length_[1] / div_[1]) * i;
      }
    }
    // connectivity
    for(int i = 0; i < div_[1]; i++)
    {
      for(int j = 0; j < div_[0]; j++)
      {
        int tmp = div_[0] * i + j;
        element[tmp].eType = "quad4";
        element[tmp].ID = div_[0] * i + j;

        element[tmp].nodeID.resize(element[tmp].ne);
        element[tmp].nodeID[0] = div_[0] * i + j + i;
        element[tmp].nodeID[1] = div_[0] * i + j + i + 1;
        element[tmp].nodeID[2] = div_[0] * i + j + i + 1 + (div_[0] + 1);
        element[tmp].nodeID[3] = div_[0] * i + j + i + (div_[0] + 1);
      }
    }
  }
  else if(eleType == "quad8")
  {
    // coordinate
    for(int i = 0; i < div_[1] + 1; i++)
    {
      for(int j = 0; j < 2 * div_[0] + 1; j++)
      {
        int tmp = j + (3 * div_[0] + 2) * i;
        node[tmp].ID = j + (3 * div_[0] + 2) * i;
        node[tmp].x[0] = 0.5 * ((double)length_[0] / div_[0]) * j;
        node[tmp].x[1] = ((double)length_[1] / div_[1]) * i;
      }
    }
    for(int i = 0; i < div_[1]; i++)
    {
      for(int j = 0; j < div_[0] + 1; j++)
      {
        int tmp = j + (3 * div_[0] + 2) * i + 2 * div_[0] + 1;
        node[tmp].ID = j + (3 * div_[0] + 2) * i + 2 * div_[0] + 1;
        node[tmp].x[0] = ((double)length_[0] / div_[0]) * j;
        node[tmp].x[1] = 0.5 * ((double)length_[1] / div_[1]) +
                         ((double)length_[1] / div_[1]) * i;
      }
    }

    // connectivity
    for(int i = 0; i < div_[1]; i++)
    {
      for(int j = 0; j < div_[0]; j++)
      {
        int tmp = div_[0] * i + j;
        element[tmp].eType = "quad8";
        element[tmp].ID = div_[0] * i + j;

        element[tmp].nodeID.resize(element[tmp].ne);

        element[tmp].nodeID[0] = 2 * j + (3 * div_[0] + 2) * i;
        element[tmp].nodeID[1] = 2 * j + (3 * div_[0] + 2) * i + 2;
        element[tmp].nodeID[2] =
            2 * j + (3 * div_[0] + 2) * i + 3 * div_[0] + 4;
        element[tmp].nodeID[3] =
            2 * j + (3 * div_[0] + 2) * i + 3 * div_[0] + 2;
        element[tmp].nodeID[4] = 2 * j + (3 * div_[0] + 2) * i + 1;
        element[tmp].nodeID[5] =
            2 * j + (3 * div_[0] + 2) * i + 2 * (div_[0] + 1) - j;
        element[tmp].nodeID[6] =
            2 * j + (3 * div_[0] + 2) * i + 3 * div_[0] + 3;
        element[tmp].nodeID[7] =
            2 * j + (3 * div_[0] + 2) * i + 2 * (div_[0] + 1) - j - 1;
      }
    }
  }
  else if(eleType == "mitc4")
  {
    // coordinate
    for(int i = 0; i < div_[1] + 1; i++)
    {
      for(int j = 0; j < div_[0] + 1; j++)
      {
        int tmp = j + i + div_[0] * i;
        node[tmp].ID = j + i + div_[0] * i;

        node[tmp].x[0] = ((double)length_[0] / div_[0]) * j;
        node[tmp].x[1] = ((double)length_[1] / div_[1]) * i;
        node[tmp].x[2] = 0.0;
      }
    }
    // connectivity
    for(int i = 0; i < div_[1]; i++)
    {
      for(int j = 0; j < div_[0]; j++)
      {
        int tmp = div_[0] * i + j;
        element[tmp].ID = div_[0] * i + j;
        element[tmp].eType = "mitc4";

        element[tmp].nodeID.resize(element[tmp].ne);
        element[tmp].nodeID[0] = div_[0] * i + j + i;
        element[tmp].nodeID[1] = div_[0] * i + j + i + 1;
        element[tmp].nodeID[2] = div_[0] * i + j + i + 1 + (div_[0] + 1);
        element[tmp].nodeID[3] = div_[0] * i + j + i + (div_[0] + 1);
      }
    }
  }
  else if(eleType == "hexa8")
  {
    // coordinate
    for(int i = 0; i < div_[2] + 1; i++)
    {
      for(int j = 0; j < div_[1] + 1; j++)
      {
        for(int k = 0; k < div_[0] + 1; k++)
        {
          int tmp = k + j * (div_[0] + 1) + i * (div_[0] + 1) * (div_[1] + 1);
          node[tmp].ID =
              k + j * (1 + div_[0]) + i * (1 + div_[0]) * (1 + div_[1]);
          node[tmp].x[0] = ((double)length_[0] / div_[0]) * k;
          node[tmp].x[1] = ((double)length_[1] / div_[1]) * j;
          node[tmp].x[2] = ((double)length_[2] / div_[2]) * i;
        }
      }
    }
    // connectivity
    for(int i = 0; i < div_[2]; i++)
    {
      for(int j = 0; j < div_[1]; j++)
      {
        for(int k = 0; k < div_[0]; k++)
        {
          int tmp = div_[0] * div_[1] * i + div_[0] * j + k;
          element[tmp].eType = "hexa8";
          element[tmp].ID = div_[0] * div_[1] * i + div_[0] * j + k;

          element[tmp].nodeID.resize(element[tmp].ne);
          element[tmp].nodeID[0] =
              (div_[0] + 1) * (div_[1] + 1) * i + (div_[0] + 1) * j + k;
          element[tmp].nodeID[1] =
              (div_[0] + 1) * (div_[1] + 1) * i + (div_[0] + 1) * j + k + 1;
          element[tmp].nodeID[2] = (div_[0] + 1) * (div_[1] + 1) * i +
                                   (div_[0] + 1) * j + k + 1 + (div_[0] + 1);
          element[tmp].nodeID[3] = (div_[0] + 1) * (div_[1] + 1) * i +
                                   (div_[0] + 1) * j + k + (div_[0] + 1);

          element[tmp].nodeID[4] = (div_[0] + 1) * (div_[1] + 1) * i +
                                   (div_[0] + 1) * j + k +
                                   (div_[0] + 1) * (div_[1] + 1);
          element[tmp].nodeID[5] = (div_[0] + 1) * (div_[1] + 1) * i +
                                   (div_[0] + 1) * j + k + 1 +
                                   (div_[0] + 1) * (div_[1] + 1);
          element[tmp].nodeID[6] = (div_[0] + 1) * (div_[1] + 1) * i +
                                   (div_[0] + 1) * j + k + 1 + (div_[0] + 1) +
                                   (div_[0] + 1) * (div_[1] + 1);
          element[tmp].nodeID[7] = (div_[0] + 1) * (div_[1] + 1) * i +
                                   (div_[0] + 1) * j + k + (div_[0] + 1) +
                                   (div_[0] + 1) * (div_[1] + 1);
        }
      }
    }
  }
  else if(eleType == "hexa20")
  {
    // coordinate
    int counter = 0;
    for(int i = 0; i < 2 * div_[2] + 1; i++) // z
    {
      if(i % 2 == 0)
      {
        for(int j = 0; j < 2 * div_[1] + 1; j++) // y
        {
          if(j % 2 == 0)
          {
            for(int k = 0; k < 2 * div_[0] + 1; k++) // x
            {
              node[counter].ID = counter;
              node[counter].x[0] = k * (length_[0] / div_[0]) / 2.0;
              node[counter].x[1] = j * (length_[1] / div_[1]) / 2.0;
              node[counter].x[2] = i * (length_[2] / div_[2]) / 2.0;
              counter++;
            }
          }
          else
          {
            for(int k = 0; k < div_[0] + 1; k++) // x
            {
              node[counter].ID = counter;
              node[counter].x[0] = k * (length_[0] / div_[0]);
              node[counter].x[1] = j * (length_[1] / div_[1]) / 2.0;
              node[counter].x[2] = i * (length_[2] / div_[2]) / 2.0;
              counter++;
            }
          }
        }
      }
      else
      {
        for(int j = 0; j < div_[1] + 1; j++) // y
        {
          for(int k = 0; k < div_[0] + 1; k++) // x
          {
            node[counter].ID = counter;
            node[counter].x[0] = k * (length_[0] / div_[0]);
            node[counter].x[1] = j * (length_[1] / div_[1]);
            node[counter].x[2] = i * (length_[2] / div_[2]) / 2.0;
            counter++;
          }
        }
      }
    }

    // connectivity
    counter = 0;
    for(int i = 0; i < div_[2]; i++) // z
    {
      for(int j = 0; j < div_[1]; j++) // y
      {
        for(int k = 0; k < div_[0]; k++) // x
        {
          int tmp = counter;
          element[tmp].ID = tmp;
          element[tmp].eType = "hexa20";

          element[tmp].nodeID.resize(element[tmp].ne);

          element[tmp].nodeID[0] =
              i * (4 * div_[0] * div_[1] + 3 * div_[0] + 3 * div_[1] + 2) +
              j * (3 * div_[0] + 2) + k * 2;                   // x- y- z-
          element[tmp].nodeID[1] = element[tmp].nodeID[0] + 2; // x+ y- z-
          element[tmp].nodeID[2] =
              element[tmp].nodeID[1] + 3 * div_[0] + 2; // x+ y+ z-
          element[tmp].nodeID[3] =
              element[tmp].nodeID[0] + 3 * div_[0] + 2; // x- y+ z-
          element[tmp].nodeID[4] = element[tmp].nodeID[0] +
                                   4 * div_[0] * div_[1] + 3 * div_[0] +
                                   3 * div_[1] + 2; // x- y+ z+
          element[tmp].nodeID[5] = element[tmp].nodeID[1] +
                                   4 * div_[0] * div_[1] + 3 * div_[0] +
                                   3 * div_[1] + 2; // x- y+ z+
          element[tmp].nodeID[6] = element[tmp].nodeID[2] +
                                   4 * div_[0] * div_[1] + 3 * div_[0] +
                                   3 * div_[1] + 2; // x- y+ z+
          element[tmp].nodeID[7] = element[tmp].nodeID[3] +
                                   4 * div_[0] * div_[1] + 3 * div_[0] +
                                   3 * div_[1] + 2; // x- y+ z+
          element[tmp].nodeID[8] =
              element[tmp].nodeID[0] + 1; // between 0 and 1
          element[tmp].nodeID[9] = element[tmp].nodeID[8] + (div_[0] - k) * 2 +
                                   k + 1; // between 1 and 2
          element[tmp].nodeID[10] =
              element[tmp].nodeID[8] + 3 * div_[0] + 2; // between 2 and 3
          element[tmp].nodeID[11] = element[tmp].nodeID[0] + (div_[0] - k) * 2 +
                                    k + 1; // between 3 and 0
          element[tmp].nodeID[12] = element[tmp].nodeID[8] +
                                    4 * div_[0] * div_[1] + 3 * div_[0] +
                                    3 * div_[1] + 2; // between 4 and 5
          element[tmp].nodeID[13] = element[tmp].nodeID[9] +
                                    4 * div_[0] * div_[1] + 3 * div_[0] +
                                    3 * div_[1] + 2; // between 5 and 6
          element[tmp].nodeID[14] = element[tmp].nodeID[10] +
                                    4 * div_[0] * div_[1] + 3 * div_[0] +
                                    3 * div_[1] + 2; // between 6 and 7
          element[tmp].nodeID[15] = element[tmp].nodeID[11] +
                                    4 * div_[0] * div_[1] + 3 * div_[0] +
                                    3 * div_[1] + 2; // between 7 and 4
          element[tmp].nodeID[16] =
              element[tmp].nodeID[0] + (div_[1] - j) * (3 * div_[0] + 2) +
              j * (div_[0] + 1) + (div_[0] - k) * 2 + k + 1; // between 0 and 4
          element[tmp].nodeID[17] =
              element[tmp].nodeID[16] + 1; // between 1 and 5
          element[tmp].nodeID[18] =
              element[tmp].nodeID[17] + div_[0] + 1; // between 2 and 6
          element[tmp].nodeID[19] =
              element[tmp].nodeID[16] + div_[0] + 1; // between 3 and 7

          counter++;
        }
      }
    }
  }
}

template <class E, class N>
void MeshIcarat<E, N>::modifyMesh(std::vector<E> &element, std::vector<N> &node)
{
  using selection_type = std::array<std::array<double, 2>, 3>;
  const auto &select = toml::find<toml::array>(config_, "mesh", "selection");
  double tol = 1.0e-7;

  // 1. create isDesignNode & isDesignElement
  std::vector<bool> isDesignNode(fem_.numnp, false),
      isDesignElement(fem_.nelx, false);
  for(auto &i_select : select)
  {
    auto data = toml::get<selection_type>(i_select);
    double xmin = std::get<0>(std::get<0>(data)) - tol;
    double xmax = std::get<1>(std::get<0>(data)) + tol;
    double ymin = std::get<0>(std::get<1>(data)) - tol;
    double ymax = std::get<1>(std::get<1>(data)) + tol;
    double zmin = std::get<0>(std::get<2>(data)) - tol;
    double zmax = std::get<1>(std::get<2>(data)) + tol;
    if(fem_.ndim == 2)
    {
      for(int i = 0; i < fem_.numnp; i++)
        if(node[i].x[0] >= xmin && node[i].x[0] <= xmax && //
           node[i].x[1] >= ymin && node[i].x[1] <= ymax)
          isDesignNode[i] = true;
      for(int i = 0; i < fem_.nelx; i++)
      {
        double centerX = 0.0, centerY = 0.0;
        for(int j = 0; j < element[i].ne; j++)
        {
          centerX += node[element[i].nodeID[j]].x[0] / element[i].ne;
          centerY += node[element[i].nodeID[j]].x[1] / element[i].ne;
        }
        if(centerX >= xmin && centerX <= xmax && //
           centerY >= ymin && centerY <= ymax)
          isDesignElement[i] = true;
      }
    }
    else // 3D
    {
      for(int i = 0; i < fem_.numnp; i++)
        if(node[i].x[0] >= xmin && node[i].x[0] <= xmax && //
           node[i].x[1] >= ymin && node[i].x[1] <= ymax && //
           node[i].x[2] >= zmin && node[i].x[2] <= zmax)
          isDesignNode[i] = true;
      for(int i = 0; i < fem_.nelx; i++)
      {
        double centerX = 0.0, centerY = 0.0, centerZ = 0.0;
        for(int j = 0; j < element[i].ne; j++)
        {
          centerX += node[element[i].nodeID[j]].x[0] / element[i].ne;
          centerY += node[element[i].nodeID[j]].x[1] / element[i].ne;
          centerZ += node[element[i].nodeID[j]].x[2] / element[i].ne;
        }
        if(centerX >= xmin && centerX <= xmax && //
           centerY >= ymin && centerY <= ymax && //
           centerZ >= zmin && centerZ <= zmax)
          isDesignElement[i] = true;
      }
    }
  } // selection loop

  // 2. create numnp_new, nelx_new
  int numnp_new = std::count(isDesignNode.begin(), isDesignNode.end(), true);
  int nelx_new =
      std::count(isDesignElement.begin(), isDesignElement.end(), true);

  // 3. create node_new
  std::vector<int> nodeID_new(fem_.numnp);
  int counter = 0;
  std::vector<N> node_new;
  for(int i = 0; i < fem_.numnp; i++)
  {
    if(isDesignNode[i])
    {
      node[i].ID = counter;
      nodeID_new[i] = counter;
      counter++;
      node_new.push_back(node[i]);
    }
    else
    {
      nodeID_new[i] = -1;
    }
  }
  assert(node_new.size() == numnp_new);

  // 4. create element_new
  counter = 0;
  std::vector<E> element_new;
  for(int i = 0; i < fem_.nelx; i++)
  {
    if(isDesignElement[i])
    {
      element[i].ID = counter;
      counter++;
      for(int j = 0; j < element[i].ne; j++)
      {
        int newID = nodeID_new[element[i].nodeID[j]];
        if(newID == -1)
        {
          std::cerr << "error in mesh_icarat" << std::endl;
          exit(1);
        }
        element[i].nodeID[j] = newID;
      }
      element_new.push_back(element[i]);
    }
  }
  assert(element_new.size() == nelx_new);

  // 5. update each variables
  fem_.numnp = numnp_new;
  fem_.nelx = nelx_new;
  fem_.neq = fem_.numnp * fem_.dofnp;

  element = element_new;
  node = node_new;

  return;
}
} // namespace icarat