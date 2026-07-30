///  @file  mesh_vtk.hpp
///  @author  Daiki Watanabe
///  @date  October 15, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <iostream>
#include <problem/base.hpp>
#include <string>

namespace icarat
{
///  This class can set element's connectivity & node's coordinate
/// and generate element information (eType, ID, ne and nodeID)
/// from VTK file
/// This class also supports unstructured grids.
/// This class also supports restart feature of computation.(only ascii file)
template <class E, class N> class MeshVTK
{
public:
  /// constructor
  MeshVTK(FEM &fem, std::string &pathVTK, const toml::value &config);

  ///  set FEM data from dat file
  void setParameter(std::string probType);

  ///  main method of generating mesh
  void generate(std::vector<E> &element, std::vector<N> &node);

  /// import coordinate data
  ///@return size is [numnp][3]
  std::vector<std::vector<double>> importNodes();

  /// import point scalar data
  ///@return size is numnp
  std::vector<double> importPointScalars(std::string fieldName);

  /// import point vector data
  ///@return size is 3*numnp
  std::vector<double> importPointVectors(std::string fieldName);

  /// import element scalar data
  ///@return size is nelx
  std::vector<double> importElementScalars(std::string fieldName);

  /// getters
  /// number of element type
  int numeleTypes() { return (int)eleTypes_.size(); }
  /// ne of each element type
  int ne(int i) { return eleTypeToNe(eleTypes_[i].first); }
  /// ipmax of each element type
  int ipmax(int i) { return eleTypeToIpmax(eleTypes_[i].first); }
  /// nelx of each element type
  int nelx(int i) { return nelxs_[i]; }

private:
  FEM &fem_;
  std::string &pathVTK_;
  const toml::value &config_;
  std::vector<std::pair<std::string, std::string>>
      eleTypes_; ///<  element types (are selected from searchE_)
  std::vector<std::string> buffCType_; ///< buffer of cell type
  std::vector<int> nelxs_; ///< number of each element type's element
  std::vector<std::pair<std::string, std::string>>
      searchE_;              ///< element type and "CELL_TYPES"
  std::ifstream fin;         ///< VTK file
  bool withoffsets_ = false; ///< Determine if the file has "OFFSETS" set.
};

//////////////////////////////////////////////////
/////////////////////public///////////////////////
//////////////////////////////////////////////////

template <class E, class N>
MeshVTK<E, N>::MeshVTK(FEM &fem, std::string &pathVTK,
                       const toml::value &config)
    : fem_(fem), pathVTK_(pathVTK), config_(config), fin(pathVTK)
{
  typedef std::pair<std::string, std::string> pair;
  searchE_.push_back(pair("tria3", "5"));
  searchE_.push_back(pair("tria6", "22"));
  searchE_.push_back(pair("quad4", "9"));
  searchE_.push_back(pair("quad8", "23"));
  searchE_.push_back(pair("tetra4", "10"));
  searchE_.push_back(pair("tetra10", "24"));
  searchE_.push_back(pair("hexa8", "12"));
  searchE_.push_back(pair("hexa20", "25"));
}

template <class E, class N>
void MeshVTK<E, N>::setParameter(std::string probType)
{
  // check binary or ascii
  std::string tmp;
  bool judge_ascii = false;
  while(fin >> tmp)
  {
    if(tmp == "ASCII")
    {
      judge_ascii = true;
      break;
    }
  }
  if(judge_ascii == false)
  {
    std::cerr << "Please set ascii vtk file in mesh_vtk.hpp" << std::endl;
    exit(1);
  }

  // dimension
  fem_.ndim = toml::find<int>(config_, "mesh", "dimension");
  std::cout << "dimension : " << fem_.ndim << std::endl;

  // dofnp, voigt
  std::tie(fem_.dofnp, fem_.voigt) = probTypeToDofnpVoigt(fem_.ndim, probType);

  // element type
  auto eleType = toml::find<std::vector<std::string>>(config_, "mesh", "type");
  for(int i = 0; i < (int)eleType.size(); i++)
    for(int j = 0; j < (int)searchE_.size(); j++)
      if(eleType[i] == searchE_[j].first)
        eleTypes_.push_back(searchE_[j]);

  std::cout << "number of element type : " << (int)eleTypes_.size()
            << std::endl;

  // reset fin
  this->fin.clear();
  this->fin.seekg(0, std::ios_base::beg);

  // get nelx (total), numnp
  while(this->fin >> tmp)
  {
    if(tmp == "POINTS")
    {
      this->fin >> fem_.numnp;
    }
    if(tmp == "OFFSETS")
    {
      withoffsets_ = true;
    }
  }

  std::cout << "number of node : " << fem_.numnp << std::endl;

  fem_.neq = fem_.numnp * fem_.dofnp;
  std::cout << "total DOFs : " << fem_.neq << std::endl;

  // reset fin
  this->fin.clear();
  this->fin.seekg(0, std::ios_base::beg);

  // count cell types
  int nelxCT = 0;
  std::vector<int> eleCounter((int)eleTypes_.size());
  while(this->fin >> tmp)
  {
    if(tmp == "CELL_TYPES")
    {
      fin >> nelxCT;

      // element loop
      for(int nel = 0; nel < nelxCT; nel++)
      {
        fin >> tmp;
        // find cell type
        for(int type = 0; type < (int)eleTypes_.size(); type++)
          if(eleTypes_[type].second == tmp)
            eleCounter[type]++;

        buffCType_.push_back(tmp);
      }
      break;
    }
  }

  // set nelx
  nelxs_ = eleCounter;
  fem_.nelx = 0;
  for(int i = 0; i < eleCounter.size(); i++)
    fem_.nelx += eleCounter[i];

  std::cout << "number of element : " << fem_.nelx << std::endl;
}

template <class E, class N>
void MeshVTK<E, N>::generate(std::vector<E> &element, std::vector<N> &node)
{
  std::string tmp;

  // reset fin
  this->fin.clear();
  this->fin.seekg(0, std::ios_base::beg);

  while(this->fin >> tmp)
  {
    // node ID, coordinate
    if(tmp == "POINTS")
    {
      fin >> tmp >> tmp;
      for(int i = 0; i < fem_.numnp; i++)
      {
        node[i].ID = i;
        this->fin >> node[i].x[0] >> node[i].x[1] >> node[i].x[2];
      }
    }
    // element ID, nodeID (with offsets version)
    if(tmp == "CONNECTIVITY" && withoffsets_)
    {
      // Skip the extra parts.
      fin >> tmp;
      for(int nel = 0; nel < fem_.nelx; nel++)
      {
        element[nel].ID = nel;

        for(int type = 0; type < (int)eleTypes_.size(); type++)
        {
          // find a match for the cell type
          if(buffCType_[nel] == eleTypes_[type].second)
          {
            element[nel].eType = eleTypes_[type].first;
            int ne = eleTypeToNe(eleTypes_[type].first);
            // get connectivity
            element[nel].nodeID.resize(ne);
            for(int i = 0; i < ne; i++)
              fin >> element[nel].nodeID[i];
          }
        }
      }
      ////////////
      break;
      ////////////
    }
    // element ID, nodeID (without offsets version)
    if(tmp == "CELLS" && !withoffsets_)
    {
      // Skip the extra parts.
      fin >> tmp >> tmp;

      for(int nel = 0; nel < fem_.nelx; nel++)
      {
        element[nel].ID = nel;

        for(int type = 0; type < (int)eleTypes_.size(); type++)
        {
          // find a match for the cell type
          if(buffCType_[nel] == eleTypes_[type].second)
          {
            // Skip the extra parts.
            fin >> tmp;

            element[nel].eType = eleTypes_[type].first;
            int ne = eleTypeToNe(eleTypes_[type].first);
            // get connectivity
            element[nel].nodeID.resize(ne);
            for(int i = 0; i < ne; i++)
              fin >> element[nel].nodeID[i];
          }
        }
      }
      ////////////
      break;
      ////////////
    }
  }

  std::cout << "check point at mesh generate           --->  ok" << std::endl;
}

template <class E, class N>
std::vector<std::vector<double>> MeshVTK<E, N>::importNodes()
{
  std::vector<std::vector<double>> nodes(fem_.numnp);
  std::string tmp;
  bool flag = false;
  std::string fieldName = "POINTS";

  // reset fin
  this->fin.clear();
  this->fin.seekg(0, std::ios_base::beg);

  while(this->fin >> tmp)
  {
    // count node
    if(tmp == fieldName)
    {
      int nodesnum;
      fin >> nodesnum;
      fin >> tmp;
      for(int i = 0; i < nodesnum; i++)
      {
        double x, y, z;
        fin >> x >> y >> z;
        nodes.push_back({x, y, z});
      }
      flag = true;
      break;
    }
  }
  if(!flag)
  {
    std::cerr << fieldName + " couldn't load in mesh_vtk.hpp" << std::endl;
    exit(1);
  }

  return nodes;
}

template <class E, class N>
std::vector<double> MeshVTK<E, N>::importPointScalars(std::string fieldName)
{
  std::vector<double> poisca(fem_.numnp);
  std::string tmp;
  bool flag = false;

  // reset fin
  this->fin.clear();
  this->fin.seekg(0, std::ios_base::beg);

  while(this->fin >> tmp)
  {
    // count node
    if(tmp == fieldName)
    {
      // skip extra parts
      fin >> tmp >> tmp >> tmp;

      for(int np = 0; np < fem_.numnp; np++)
      {
        fin >> poisca[np];
      }
      flag = true;
      break;
    }
  }
  if(flag == false)
  {
    std::cerr << fieldName + " couldn't load in mesh_vtk.hpp" << std::endl;
    exit(1);
  }

  return poisca;
}

template <class E, class N>
std::vector<double> MeshVTK<E, N>::importPointVectors(std::string fieldName)
{
  std::vector<double> poivec(3 * fem_.numnp);
  std::string tmp;
  bool flag = false;

  // reset fin
  this->fin.clear();
  this->fin.seekg(0, std::ios_base::beg);

  while(this->fin >> tmp)
  {
    // count node
    if(tmp == fieldName)
    {
      // skip extra parts
      fin >> tmp;

      for(int np = 0; np < 3 * fem_.numnp; np++)
      {
        fin >> poivec[np];
      }
      flag = true;
      break;
    }
  }
  if(flag == false)
  {
    std::cerr << fieldName + " couldn't load in mesh_vtk.hpp" << std::endl;
    exit(1);
  }

  return poivec;
}

template <class E, class N>
std::vector<double> MeshVTK<E, N>::importElementScalars(std::string fieldName)
{
  std::vector<double> elesca(fem_.nelx);
  std::string tmp;
  bool flag = false;

  // reset fin
  this->fin.clear();
  this->fin.seekg(0, std::ios_base::beg);

  while(this->fin >> tmp)
  {
    // count node
    if(tmp == fieldName)
    {
      // skip extra parts
      fin >> tmp >> tmp >> tmp;

      for(int nel = 0; nel < fem_.nelx; nel++)
      {
        fin >> elesca[nel];
      }
      flag = true;
      break;
    }
  }
  if(flag == false)
  {
    std::cerr << fieldName + " couldn't load in mesh_vtk.hpp" << std::endl;
    exit(1);
  }

  return elesca;
}

} // namespace icarat
