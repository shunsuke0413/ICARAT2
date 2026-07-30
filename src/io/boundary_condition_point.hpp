///  @file  boundary_condition_points.hpp
///  @author  Daiki Watanabe
///  @date  February 26, 2022.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include "boundary_condition_base.hpp"
#include "file_reader.hpp"

namespace icarat
{

/// class of setting Dirichlet and Neumann boundary condition from .gid files
template <class E, class N> class BCPoint
{
public:
  BCPoint(FEM &fem, std::vector<E> &element, std::vector<N> &node,
          std::string project)
      : fem(fem), element(element), node(node), project(project)
  {
  }

  ///  input Dirichlet condition
  void input_dirich(const toml::value &config);

  /// set dof to each node
  /// @return fem.numeq(= total ODFs - dirichlet DOFs)
  int makeDOF();

  /// read & set Neumann condition's informations
  void input_neumann(const toml::value &config, Force &force);

private:
  ///  main method of generating the equivalent nodal load
  void mainLoadLoop(int loadID, std::vector<int> &loadNode, Force &force);

  /// generating the equivalent nodal load with focusing on one element
  void loadOneElement(int nel, int loadID, std::vector<int> &loadNode,
                      Force &force);

  double getLoadAreaDistance(int nel, int loadID, std::vector<int> &loadNode);

  FEM &fem;
  std::vector<E> &element;
  std::vector<N> &node;
  std::string project;
  std::vector<InputLoad> inp_l; ///< load structure
};

/////////////////////////////////////
//////////////// public /////////////
/////////////////////////////////////

template <class E, class N>
void BCPoint<E, N>::input_dirich(const toml::value &config)
{
  using dirich_type = std::tuple<std::array<int, 6>, std::array<double, 6>,
                                 std::array<std::string, 1>>;

  const auto &dirich =
      toml::find<toml::array>(config, "boundary_condition", "dirichlet");

  for(int dirichid = 0; dirichid < (int)dirich.size(); dirichid++)
  {
    dirich_type data = toml::get<dirich_type>(dirich[dirichid]);
    auto flag = std::get<0>(data);
    auto val = std::get<1>(data);
    auto csv_array = std::get<2>(data);

    auto flag_vec = std::vector<int>(flag.begin(), flag.end());
    auto val_vec = std::vector<double>(val.begin(), val.end());
    auto csv = read_file(project + "/" + csv_array[0]);

    std::string tmp;
    for(int i = 0; i < (int)csv.size(); i++)
    {
      std::istringstream iss(csv[i]);
      iss >> tmp;

      setDirichValue(fem.dofnp, node[stoi(tmp)], flag_vec, val_vec);
    }
  }

#ifdef ICARAT_MPI
  if(fem.procid == 0)
#endif
    std::cout << "check point at dirich at macro         --->  ok" << std::endl;

  return;
}

template <class E, class N> int BCPoint<E, N>::makeDOF()
{
  int numeq = countDOF(fem, node);

  return numeq;
}

template <class E, class N>
void BCPoint<E, N>::input_neumann(const toml::value &config, Force &force)
{
  using neum_type = std::tuple<std::array<double, 3>,
                               std::array<std::string, 1>, std::array<int, 1>>;
  using neum_mome_type =
      std::tuple<std::array<double, 1>, std::array<double, 3>,
                 std::array<double, 3>, std::array<std::string, 1>,
                 std::array<int, 1>>;
  using neum_pres_type = std::tuple<double, std::array<std::string, 1>>;
  const auto &neumann =
      toml::find<toml::array>(config, "boundary_condition", "neumann");
  const auto &neumann_mome =
      toml::find<toml::array>(config, "boundary_condition", "neumann_moment");
  const auto &neumann_pres =
      toml::find<toml::array>(config, "boundary_condition", "neumann_pressure");

  int counter = 0;

  // neumann loop
  for(int neumid = 0; neumid < (int)neumann.size(); neumid++)
  {
    auto data = toml::get<neum_type>(neumann[neumid]);
    // init
    std::vector<int> loadNode; // node where load is applied
    std::string tmp;
    auto csv_array = std::get<1>(data);
    auto csv = read_file(project + "/" + csv_array[0]);

    for(int i = 0; i < (int)csv.size(); i++)
    {
      std::istringstream iss(csv[i]);
      iss >> tmp;
      loadNode.push_back(stoi(tmp));
    }

    auto val = std::get<0>(data);
    std::vector<double> val_vec = std::vector<double>(val.begin(), val.end());
    InputLoad inp;
    inp.val = val_vec;
    inp_l.push_back(inp);

    // 2. generate nodal loads for each element.
    mainLoadLoop(counter, loadNode, force);

    counter++;
  }

  // moment loop
  for(int neumid = 0; neumid < (int)neumann_mome.size(); neumid++)
  {
    auto data = toml::get<neum_mome_type>(neumann_mome[neumid]);
    // init
    std::vector<int> loadNode; // node where load is applied
    std::string tmp;
    auto csv_array = std::get<3>(data);
    auto csv = read_file(project + "/" + csv_array[0]);

    for(int i = 0; i < (int)csv.size(); i++)
    {
      std::istringstream iss(csv[i]);
      iss >> tmp;
      loadNode.push_back(stoi(tmp));
    }

    auto val = std::get<0>(data);
    auto axis_origin = std::get<1>(data);
    auto axis_vector = std::get<2>(data);
    InputLoad inp;
    inp.val.resize(fem.ndim);
    inp.val[0] = val[0];
    inp.axis_origin = axis_origin;
    inp.axis_vector = axis_vector;
    inp.type = "moment";
    inp_l.push_back(inp);

    // 2. generate nodal loads for each element.
    mainLoadLoop(counter, loadNode, force);

    counter++;
  }

  // // pressure loop
  // for(int neumid = 0; neumid < (int)neumann_pres.size(); neumid++)
  // {
  //   auto data = toml::get<neum_pres_type>(neumann_pres[neumid]);
  //   // init
  //   std::vector<int> loadNode; // node where load is applied
  //   std::string tmp;
  //   auto csv = read_file(std::get<1>(data));

  //   for(int i = 0; i < (int)csv.size(); i++)
  //   {
  //     std::istringstream iss(csv[i]);
  //     iss >> tmp;
  //     loadNode.push_back(stoi(tmp));
  //   }

  //   auto val = std::get<0>(data);
  //   InputLoad inp;
  //   inp.val.resize(fem.ndim);
  //   inp.val[0] = val;
  //   inp.type = "pressure";
  //   inp_l.push_back(inp);

  //   // 2. generate nodal loads for each element.
  //   mainLoadLoop(counter, loadNode, force);

  //   counter++;
  // }

#ifdef ICARAT_MPI
  if(fem.procid == 0)
#endif
    std::cout << "check point at neumann at macro        -->  ok" << std::endl;

  return;
}

// /////////////////////////////////////
// //////////////// private/////////////
// /////////////////////////////////////

template <class E, class N>
void BCPoint<E, N>::mainLoadLoop(int loadID, std::vector<int> &loadNode,
                                 Force &force)
{
  // compute sum of load area * distance
  // if(inp_l[loadID].type == "moment")
  // {
  //   inp_l[loadID].sumareadist = 0.0;
  //   for(int nel = 0; nel < fem.nelx; nel++)
  //     inp_l[loadID].sumareadist += getLoadAreaDistance(nel, loadID,
  //     loadNode);

  //   inp_l[loadID].val[0] /= inp_l[loadID].sumareadist;
  // }

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

  return;
}

template <class E, class N>
void BCPoint<E, N>::loadOneElement(int nel, int loadID,
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
          return;
        }
      }
    }
  }

  return;
}

template <class E, class N>
double BCPoint<E, N>::getLoadAreaDistance(int nel, int loadID,
                                          std::vector<int> &loadNode)
{
}
} // namespace icarat
