///  @file  mesh_mdpa.hpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <io/boundary_condition_gid.hpp>
#include <io/file_reader.hpp>
#include <problem/base.hpp>
namespace icarat
{

///  This class sets element's connectivity & node's coordinate
/// and generate element information (eType, ID and nodeID)
/// from MDPA file
/// This class also supports unstructured grids.
template <class E, class N> class MeshMDPA
{
public:
  /// constructor
  MeshMDPA(FEM &fem, std::vector<std::string> &mdpa, const toml::value &config);

  ///  set FEM data from dat file
  void setParameter(std::string probType);

  ///  main method of generating mesh
  void generate(std::vector<E> &element, std::vector<N> &node);

  /// getters
  /// number of element type
  int numeleTypes() { return (int)components_.size(); }
  /// ne of each element type
  int ne(int i) { return eleTypeToNe(components_[i].first); }
  /// ipmax of each element type
  int ipmax(int i) { return eleTypeToIpmax(components_[i].first); }
  /// nelx of each element type
  int nelx(int i) { return nelxs_[i]; }

private:
  void inputMesh(int counter, std::string eType, int ne, std::string &line,
                 std::vector<E> &element, std::vector<N> &node);

  FEM &fem_;
  std::vector<std::string> &mdpa_;
  const toml::value &config_;
  std::vector<std::pair<std::string, std::string>>
      eleTypes_; ///<  element types (are selected from searchE_)
  std::vector<std::pair<std::string, std::string>> components_; ///<  components
  std::vector<int> nelxs_;           ///< each element type's nelx
  std::vector<std::string> searchN_; ///< keyword for node
  std::vector<std::pair<std::string, std::string>>
      searchE_; ///< keyword for element
};

//////////////////////////////////////////////////
/////////////////////public///////////////////////
//////////////////////////////////////////////////

template <class E, class N>
MeshMDPA<E, N>::MeshMDPA(FEM &fem, std::vector<std::string> &linesMDPA,
                         const toml::value &config)
    : fem_(fem), mdpa_(linesMDPA), config_(config)
{
  searchN_.push_back("Begin Nodes");
  searchN_.push_back("End Nodes");

  typedef std::pair<std::string, std::string> pair;
  searchE_.push_back(pair("tria3", "2D3"));
  searchE_.push_back(pair("tria6", "2D6"));
  searchE_.push_back(pair("quad4", "2D4"));
  searchE_.push_back(pair("quad8", "2D8"));
  searchE_.push_back(pair("tetra4", "3D4"));
  searchE_.push_back(pair("tetra10", "3D10"));
  searchE_.push_back(pair("prism6", "3D6"));
  searchE_.push_back(pair("prism15", "3D15"));
  searchE_.push_back(pair("hexa8", "3D8"));
  searchE_.push_back(pair("hexa20", "3D20"));

  searchE_.push_back(pair("end", "End Elements")); // 8
}

template <class E, class N>
void MeshMDPA<E, N>::setParameter(std::string probType)
{
  // dimension
  fem_.ndim = toml::find<int>(config_, "mesh", "dimension");
  std::cout << "dimension : " << fem_.ndim << std::endl;

  // dofnp, voigt
  std::tie(fem_.dofnp, fem_.voigt) = probTypeToDofnpVoigt(fem_.ndim, probType);

  // element type
  auto eleType = toml::find<toml::array>(config_, "mesh", "type");
  std::vector<std::pair<std::string, std::string>> eleTypes;
  for(int i = 0; i < (int)eleType.size(); i++)
  {
    auto tmp = toml::get<std::string>(eleType[i]);
    for(int j = 0; j < (int)searchE_.size(); j++)
      if(tmp == searchE_[j].first)
        eleTypes.push_back(searchE_[j]);
  }
  eleTypes_ = eleTypes;

  // nelx,numnp,neq
  int sumNelx = 0, nodeCounter = 0;
  std::vector<int> nelxs;
  std::vector<std::pair<std::string, std::string>> components;
  for(int line = 0; line < (int)mdpa_.size(); line++)
  {
    // count node
    if(mdpa_[line].find(searchN_[0]) != std::string::npos)
    {
      // get node information
      for(int i = 0; i <= 1e8; i++)
      {
        if(mdpa_[line + i + 1] == searchN_[1])
          break;
        else
          nodeCounter++;
      }
      if(nodeCounter == 1e8)
      {
        std::cerr << "Failed to load node information in mesh_mdpa.hpp"
                  << std::endl;
        exit(1);
      }
    }

    // count element
    for(int type = 0; type < (int)eleTypes_.size(); type++)
    {
      if(mdpa_[line].find(eleTypes_[type].second) != std::string::npos)
      {
        // get element information
        int nelxTmp = 0;
        for(int i = 0; i <= 1e8; i++)
        {
          if(mdpa_[line + i + 1] == searchE_.back().second)
            break;
          else
          {
            sumNelx++;
            nelxTmp++;
          }
        }
        if(sumNelx == 1e8)
        {
          std::cerr << "Failed to load element information in mesh_mdpa.hpp"
                    << std::endl;
          exit(1);
        }

        nelxs.push_back(nelxTmp);
        components.push_back(eleTypes_[type]);
      }
    }
  }

  nelxs_ = nelxs;
  components_ = components;

  std::cout << "number of element type: " << (int)eleTypes_.size() << std::endl;

  std::cout << "number of components: " << (int)components_.size() << std::endl;

  fem_.numnp = nodeCounter;
  std::cout << "number of node : " << fem_.numnp << std::endl;

  fem_.nelx = sumNelx;
  std::cout << "number of element : " << fem_.nelx << std::endl;

  fem_.neq = fem_.numnp * fem_.dofnp;
  std::cout << "total DOFs : " << fem_.neq << std::endl;
}

template <class E, class N>
void MeshMDPA<E, N>::generate(std::vector<E> &element, std::vector<N> &node)
{
  int elementID = 0, counterEle = 0;

  for(int line = 0; line < (int)mdpa_.size(); line++)
  {
    /// node
    if(mdpa_[line] == searchN_[0])
    {
      // get node information
      for(int i = 0; i < fem_.numnp; i++)
      {
        int nodeID;
        std::istringstream iss(mdpa_[line + i + 1]);
        iss >> nodeID >> node[i].x[0] >> node[i].x[1] >> node[i].x[2];
        node[i].ID = nodeID - 1;
      }
    }

    // count element
    for(int type = 0; type < (int)eleTypes_.size(); type++)
    {
      if(mdpa_[line].find(eleTypes_[type].second) != std::string::npos)
      {
        int ne = eleTypeToNe(eleTypes_[type].first);

        // get element information
        for(int i = 0; i < nelxs_[counterEle]; i++)
        {
          inputMesh(elementID, eleTypes_[type].first, ne, mdpa_[line + 1 + i],
                    element, node);
          elementID++;
        }
        counterEle++;
      }
    }
  }

  std::cout << "check point at mesh generate           --->  ok" << std::endl;
}

//////////////////////////////////////////////////
/////////////////////private///////////////////////
//////////////////////////////////////////////////

template <class E, class N>
void MeshMDPA<E, N>::inputMesh(int counter, std::string eType, int ne,
                               std::string &line, std::vector<E> &element,
                               std::vector<N> &node)
{
  std::istringstream iss(line);
  std::string str;

  element[counter].ID = counter;
  element[counter].eType = eType;

  // skip unnecessary part
  iss >> str >> str;

  // get connectivity
  int tmp;
  element[counter].nodeID.resize(ne);
  for(int i = 0; i < ne; i++)
  {
    iss >> tmp;
    element[counter].nodeID[i] = tmp - 1;
  }
}

} // namespace icarat
