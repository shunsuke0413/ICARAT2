///  @file  export_vtk.hpp
///  @author  Daiki Watanabe
///  @date  April 15, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <fstream>
#include <iostream>
#include <vector>

namespace icarat
{
template <class T> int cellTypes(T ndim, T ne)
{
  int cType;
  if(ndim == 2)
  {
    if(ne == 3)
      cType = 5;
    else if(ne == 6)
      cType = 22;
    else if(ne == 4)
      cType = 9;
    else if(ne == 8)
      cType = 23;
  }
  else if(ndim == 3)
  {
    if(ne == 4)
      cType = 10;
    else if(ne == 8)
      cType = 12;
    else if(ne == 6)
      cType = 13;
    else if(ne == 20)
      cType = 25;
    else if(ne == 15)
      cType = 26;
  }

  return cType;
}

//********************Make Headder********************
inline void makeHeaderVTK(std::ofstream &_fout)
{
  _fout << "# vtk DataFile Version 4.1\n";
  _fout << "vtk output\n";
  _fout << "ASCII\n";
  _fout << "DATASET UNSTRUCTURED_GRID\n";
}

//********************Add Fields ********************
inline void addFieldScalarVTK(std::ofstream &_fout,
                              std::vector<std::pair<std::string, double>> data)
{
  _fout << "\nFIELD FieldData " << data.size() << std::endl;
  for(auto &i_data : data)
  {
    std::replace(i_data.first.begin(), i_data.first.end(), ' ', '_');
    _fout << i_data.first << " 1 1 double" << std::endl;
    _fout << i_data.second << std::endl;
  }
}

//********************Add Points********************
template <class N>
void setPointVTK(std::vector<N> &_nodes, std::ofstream &_fout)
{
  _fout << "\nPOINTS\t" << _nodes.size() << "\tfloat\n";
  for(auto &node : _nodes)
  {
    for(int i = 0; i < 3; i++)
    {
      _fout << node.x[i] << "\t";
    }
    _fout << std::endl;
  }
}

//********************Add Elements*******************
template <class E>
void setElementVTK(std::vector<E> &_elements, std::ofstream &_fout)
{
  int datanum = 0;
  for(auto &element : _elements)
  {
    datanum += element.nodeID.size() + 1;
  }

  _fout << "\nCELLS " << _elements.size() << "\t" << datanum << "\n";
  for(auto &element : _elements)
  {
    _fout << element.nodeID.size() << "\t";
    for(int i = 0; i < (int)(element.nodeID.size()); i++)
    {
      _fout << element.nodeID[i] << "\t";
    }
    _fout << std::endl;
  }
}

//********************Add Element Types********************
template <class E>
void setElementTypeVTK(int ndim, std::vector<E> &element, std::ofstream &_fout)
{
  _fout << "\nCELL_TYPES\t" << element.size() << "\n";
  for(auto &e : element)
    _fout << cellTypes(ndim, e.ne) << "\n";
}

//********************Add Point Scalars********************
template <class N, class T>
void addPointScalarVTK(std::string _symbol, std::vector<N> &node,
                       std::ofstream &_fout, bool _isheader, T _value)
{
  if(_isheader)
  {
    _fout << "\nPOINT_DATA\t" << node.size() << "\n";
  }
  _fout << "SCALARS " << _symbol << " float\n";
  _fout << "LOOKUP_TABLE default\n";

  for(auto &n : node)
  {
    _fout << _value(n) << std::endl;
  }
}

//********************Add Point Vectors********************
template <class N, class T>
void addPointVectorVTK(std::string _symbol, std::vector<N> &node,
                       std::ofstream &_fout, bool _isheader, T _values)
{
  if(_isheader)
  {
    _fout << "\nPOINT_DATA\t" << node.size() << "\n";
  }
  _fout << "VECTORS " << _symbol << " float\n";
  for(auto &n : node)
  {
    for(int i = 0; i < 3; i++)
    {
      _fout << _values(n)[i] << "\t";
    }
    _fout << std::endl;
  }
}

//********************Add Element Scalers********************
template <class E, class T>
void addElementScalarVTK(std::string _symbol, std::vector<E> &element,
                         std::ofstream &_fout, bool _isheader, T _value)
{
  if(_isheader)
  {
    _fout << "\nCELL_DATA\t" << element.size() << "\n";
  }
  _fout << "SCALARS " << _symbol << " float\n";
  _fout << "LOOKUP_TABLE default\n";

  for(auto e : element)
  {
    _fout << _value(e) << std::endl;
  }
}

//********************Add Cell Vectors********************
template <class E, class T>
void addElementVectorVTK(std::string _symbol, std::vector<E> &element,
                         std::ofstream &_fout, bool _isheader, T _values)
{
  if(_isheader)
  {
    _fout << "\nCELL_DATA\t" << element.size() << "\n";
  }
  _fout << "VECTORS " << _symbol << " float\n";
  for(auto &e : element)
  {
    for(int i = 0; i < 3; i++)
    {
      _fout << _values(e)[i] << "\t";
    }
    _fout << std::endl;
  }
}

} // namespace icarat
