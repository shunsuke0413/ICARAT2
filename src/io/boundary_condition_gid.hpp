///  @file  boundary_condition_gid.hpp
///  @author  Daiki Watanabe
///  @date  February 26, 2022.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include "boundary_condition_base.hpp"

namespace icarat
{

struct GIDConstraint
{
  std::string name; ///< group name in mdpa file
  std::vector<int> flag;
  std::vector<double> val;
};

struct GIDLoad : InputLoad
{
  std::string name; ///< group name in mdpa file
  // std::vector<double> val;
};

/// Class to process GID data files.
class GIDParser
{
public:
  GIDParser(std::string gidPath);

  /// read constraint settings from .gid folder
  ///@param [out] modelName model_part_name in ProjectParameters.json
  std::vector<GIDConstraint> readDirich(std::vector<std::string> &mdpa);

  std::vector<GIDLoad> readNeumann(std::vector<std::string> &mdpa);

  /// getter
  const std::string &mdpaName() const { return mdpaName_; }
  const std::string &modelName() const { return modelName_; }

private:
  std::string gidPath_;   ///< path to gid folder
  std::string mdpaName_;  ///< mdpa file name
  std::string modelName_; ///< model name
};

/// class of setting Dirichlet and Neumann boundary condition from .gid files
template <class E, class N> class BCGID
{
public:
  BCGID(FEM &fem, std::vector<E> &element, std::vector<N> &node,
        std::string gidPath)
      : fem(fem), element(element), node(node), gidP(gidPath){};

  ///  input Dirichlet condition
  void input_dirich(const toml::value &config, std::vector<std::string> &mdpa);

  /// set dof to each node
  /// @return fem.numeq(= total ODFs - dirichlet DOFs)
  int makeDOF();

  /// read & set Neumann condition's informations
  void input_neumann(const toml::value &config, std::vector<std::string> &mdpa,
                     Force &force);

  int counter_ = 0;

private:
  ///  main method of generating the equivalent nodal load
  void mainLoadLoop(int loadID, std::vector<int> &loadNode, Force &force);

  /// generating the equivalent nodal load with focusing on one element
  void loadOneElement(int nel, int loadID, std::vector<int> &loadNode,
                      Force &force);

  FEM &fem;
  std::vector<E> &element;
  std::vector<N> &node;
  std::vector<GIDLoad> inp_l;       ///< load structure
  std::vector<GIDConstraint> inp_c; ///< load structure
  GIDParser gidP;                   ///< GiD parser
};

/////////////////////////////////////
//////////////// public /////////////
/////////////////////////////////////

template <class E, class N>
void BCGID<E, N>::input_dirich(const toml::value &config,
                               std::vector<std::string> &mdpa)
{
  inp_c = gidP.readDirich(mdpa);
  for(int i = 0; i < (int)inp_c.size(); i++)
  {
    for(int line = 0; line < (int)mdpa.size(); line++)
    {
      // find constraint node's part
      if(mdpa[line].find(inp_c[i].name) != std::string::npos)
      {
        // init
        int counter = 0;
        bool stop = false;
        std::string tmp;
        while(stop == false)
        {
          std::istringstream iss(mdpa[line + 2 + counter]);
          iss >> tmp;
          if(tmp == "End")
            stop = true;
          else
          {
            setDirichValue(fem.dofnp, node[stoi(tmp) - 1], inp_c[i].flag,
                           inp_c[i].val);
            counter++;
          }
        }
        break;
      }
    } // mdpa file loop end
  }   // constraints loop end
#ifdef ICARAT_MPI
  if(fem.procid == 0)
#endif
    std::cout << "check point at dirich at macro         --->  ok" << std::endl;
}

template <class E, class N> int BCGID<E, N>::makeDOF()
{
  int numeq = countDOF(fem, node);

  return numeq;
}

template <class E, class N>
void BCGID<E, N>::input_neumann(const toml::value &config,
                                std::vector<std::string> &mdpa, Force &force)
{
  inp_l = gidP.readNeumann(mdpa);

  for(int constID = 0; constID < (int)inp_l.size(); constID++)
  {
    for(int line = 0; line < (int)mdpa.size(); line++)
    {
      // 1. find load's point
      if(mdpa[line].find(inp_l[constID].name) != std::string::npos)
      {
        // init
        std::vector<int> loadNode; // node where load is applied
        int counter = 0;
        bool stop = false;
        std::string tmp;
        while(stop == false)
        {
          std::istringstream iss(mdpa[line + 2 + counter]);
          iss >> tmp;
          if(tmp == "End")
            stop = true;
          else
          {
            loadNode.push_back(stoi(tmp) - 1);
            counter++;
          }
        }

        // 2. generate nodal loads for each element.
        mainLoadLoop(constID, loadNode, force);

        break;
      }
    } // mdpa file loop end
  }   // constraints loop end

#ifdef ICARAT_MPI
  if(fem.procid == 0)
#endif
    std::cout << "check point at neumann at macro        -->  ok" << std::endl;
}

/////////////////////////////////////
//////////////// private/////////////
/////////////////////////////////////

template <class E, class N>
void BCGID<E, N>::mainLoadLoop(int loadID, std::vector<int> &loadNode,
                               Force &force)
{
  for(int nel = 0; nel < fem.nelx; nel++)
    loadOneElement(nel, loadID, loadNode, force);

  // For the vtk's output, very small fe will be removed.
  bool display = true;
  for(int i = 0; i < (int)force.fext_org.size(); i++)
    if(abs(force.fext_org[i]) < 1.0e-10 && display == true)
    {
      force.fext_org[i] = 0.0;
      std::cerr << "Too small load's values will be skipped." << std::endl;
      display = false;
    }
}

template <class E, class N>
void BCGID<E, N>::loadOneElement(int nel, int loadID,
                                 std::vector<int> &loadNode, Force &force)
{
  Eigen::MatrixXd X(element[nel].ne, fem.ndim);
  Eigen::VectorXi idof(element[nel].numdof);
  for(int i = 0; i < element[nel].ne; i++)
  {
    for(int j = 0; j < fem.dofnp; j++)
      idof(fem.dofnp * i + j) = node[element[nel].nodeID[i]].dof[j];
    for(int j = 0; j < fem.ndim; j++)
      X(i, j) = node[element[nel].nodeID[i]].x[j];
  }

  SurfaceLoad surface(fem, element[nel]);
  int numnode = surface.setNumNode();
  int nodecounter = 0;
  std::vector<int> hitpoints;
  for(int j = 0; j < element[nel].ne; j++)
  {
    int tmp = element[nel].nodeID[j];
    // check node on this line
    for(int k = 0; k < (int)loadNode.size(); k++)
    {
      if(tmp == loadNode[k])
      {
        nodecounter++;
        hitpoints.push_back(j);
        if(nodecounter == numnode)
        {
          inp_l[loadID].sideID = surface.judgeLine(j, hitpoints);
          surface.generate(force, idof, X, inp_l[loadID]);
          counter_++;
          return;
        }
      }
    }
  }
}

} // namespace icarat
