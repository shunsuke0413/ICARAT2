///  @file  mesh_base.hpp
///  @author  Daiki Watanabe
///  @date  February 8, 2022.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <iostream>
#include <string>

namespace icarat
{
/// T is string
template <typename T>
std::tuple<int, int> probTypeToDofnpVoigt(int ndim, T probType)
{
  int dofnp, voigt;
  if(probType == "structure")
  {
    dofnp = ndim;
    if(ndim == 2)
      voigt = 3;
    else
      voigt = 6;
  }
  else if(probType == "heat")
  {
    dofnp = 1;
    if(ndim == 2)
      voigt = 2;
    else
      voigt = 3;
  }
  else if(probType == "fluid")
  {
    dofnp = ndim + 1; // velocity + pressure
    if(ndim == 2)
      voigt = 3;
    else
      voigt = 6;
  }
  else if(probType == "shell")
  {
    dofnp = 6;
    voigt = 6;
  }
  else
  {
    std::cerr << "probType is incorrect in mesh_base.hpp." << std::endl;
    exit(1);
  }

  return std::forward_as_tuple(dofnp, voigt);
}

/// T is string
template <typename T> int eleTypeToNe(T eleType)
{
  int ne;

  if(eleType == "tria3")
    ne = 3;
  else if(eleType == "tria6")
    ne = 6;
  else if(eleType == "quad4")
    ne = 4;
  else if(eleType == "quad8")
    ne = 8;
  else if(eleType == "tetra4")
    ne = 4;
  else if(eleType == "tetra10")
    ne = 10;
  else if(eleType == "prism6")
    ne = 6;
  else if(eleType == "prism15")
    ne = 15;
  else if(eleType == "hexa8")
    ne = 8;
  else if(eleType == "hexa20")
    ne = 20;
  else
  {
    std::cerr << "Please choose element type in mesh_mdpa.hpp" << std::endl;
    exit(1);
  }

  return ne;
}

/// T is string
template <typename T> int eleTypeToIpmax(T eleType)
{
  int ipmax;

  if(eleType == "tria3")
    ipmax = 1;
  else if(eleType == "tria6")
    ipmax = 3;
  else if(eleType == "quad4")
    ipmax = 4;
  else if(eleType == "quad8")
    ipmax = 9;
  else if(eleType == "tetra4")
    ipmax = 1;
  else if(eleType == "tetra10")
    ipmax = 4;
  else if(eleType == "prism6")
    ipmax = 2;
  else if(eleType == "prism15")
    ipmax = 9;
  else if(eleType == "hexa8")
    ipmax = 8;
  else if(eleType == "hexa20")
    ipmax = 27;
  else
  {
    std::cerr << "Please choose element type in mesh_mdpa.hpp" << std::endl;
    exit(1);
  }

  return ipmax;
}

} // namespace icarat