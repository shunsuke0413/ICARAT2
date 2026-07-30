///  @file  boundary_condition_icarat.hpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include "boundary_condition_base.hpp"
#include <Eigen/Dense>
#include <io/mesh_icarat.hpp>
#include <problem/base.hpp>

namespace icarat
{

///  Structure used to input load condition in ICARAT.
struct IcaratLoad : InputLoad
{
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;
};

/// class of setting Dirichlet and Neumann boundary condition
///@attention This class supports only quadrilateral and hexahedral elements.
template <class E, class N> class BCIcarat
{
public:
  BCIcarat(FEM &fem, std::vector<E> &element, std::vector<N> &node)
      : fem(fem), element(element), node(node){};

  ///  input Dirichlet condition
  void input_dirich(const toml::value &config);

  /// set dof to each node
  /// @return fem.numeq(= total ODFs - dirichlet DOFs)
  int makeDOF();

  /// set dof to each node for fluid problem
  /// @return fem.numeq(= total ODFs - dirichlet DOFs)
  int makeDOFFluid();

  /// read & set Neumann condition's informations
  void input_neumann(const toml::value &config, Force &force);

  /// set dof to each node for periodic condition
  /// @return fem.numeq(= total DOFs - dirichlet DOFs)
  int makeDOF_periodic(std::vector<std::vector<std::pair<int, int>>> &pair);

private:
  ///  main method of making the equivalent nodal force
  void generateLoadMain(int numNeum, Force &force);

#ifdef ICARAT_MPI
  void generateLoadMainMPI(int numNeum, Force &force);
#endif

  FEM &fem;
  std::vector<E> &element;
  std::vector<N> &node;
  std::vector<IcaratLoad> inp_l; ///< load structure
};

/////////////////////////////////////
//////////////// public /////////////
/////////////////////////////////////

template <class E, class N>
void BCIcarat<E, N>::input_dirich(const toml::value &config)
{
  using dirich_type = std::tuple<std::array<int, 6>, std::array<double, 6>,
                                 std::array<std::array<double, 2>, 3>>;
  double tol = 1.0e-7;

  const auto &dirich =
      toml::find<toml::array>(config, "boundary_condition", "dirichlet");
  int numDirich = (int)dirich.size();

  for(auto &i_dirich : dirich)
  {
    dirich_type data = toml::get<dirich_type>(i_dirich);
    auto flag = std::get<0>(data);
    auto val = std::get<1>(data);
    auto range = std::get<2>(data);

    double xmin = std::get<0>(std::get<0>(range)) - tol;
    double xmax = std::get<1>(std::get<0>(range)) + tol;
    double ymin = std::get<0>(std::get<1>(range)) - tol;
    double ymax = std::get<1>(std::get<1>(range)) + tol;
    double zmin = std::get<0>(std::get<2>(range)) - tol;
    double zmax = std::get<1>(std::get<2>(range)) + tol;

    std::vector<double> val_vec = std::vector<double>(val.begin(), val.end());
    std::vector<int> flag_vec = std::vector<int>(flag.begin(), flag.end());

    // for 2D
    if(fem.ndim == 2)
    {
      for(int j = 0; j < fem.numnp; j++)
      {
        // check node on this line
        if(node[j].x[0] >= xmin && node[j].x[0] <= xmax &&
           node[j].x[1] >= ymin && node[j].x[1] <= ymax)
        {
          setDirichValue(fem.dofnp, node[j], flag_vec, val_vec);
        }
      }
    }
    // for 3D
    else if(fem.ndim == 3)
    {
      for(int j = 0; j < fem.numnp; j++)
      {
        // check node on this line
        if(node[j].x[0] >= xmin && node[j].x[0] <= xmax &&
           node[j].x[1] >= ymin && node[j].x[1] <= ymax &&
           node[j].x[2] >= zmin && node[j].x[2] <= zmax)
        {
          setDirichValue(fem.dofnp, node[j], flag_vec, val_vec);
        }
      }
    }
    else
    {
      std::cerr << "Error in boundary_condition_icarat.hpp" << std::endl;
      exit(1);
    }
  }

#ifdef ICARAT_MPI
  if(fem.procid == 0)
#endif
    std::cout << "check point at dirich at macro         --->  ok" << std::endl;
}

template <class E, class N> int BCIcarat<E, N>::makeDOF()
{
#ifdef ICARAT_MPI
  int numeq = countDOFMPI(fem, node);
#else
  int numeq = countDOF(fem, node);
#endif
  return numeq;
}

template <class E, class N> int BCIcarat<E, N>::makeDOFFluid()
{
#ifdef ICARAT_MPI
  int numeq = countDOFMPI(fem, node);
#else
  int numeq = countDOFFluid(fem, element, node);
#endif

  return numeq;
}

template <class E, class N>
int BCIcarat<E, N>::makeDOF_periodic(
    std::vector<std::vector<std::pair<int, int>>> &pair)
{
  // set master node's ID to slave ID
  for(const auto &pai : pair)
    for(const auto &pa : pai)
    {
      int sec = pa.second;
      if(node[sec].dirich == nullptr)
      {
        std::shared_ptr<Dirichlet> dirich(new Dirichlet);
        node[sec].dirich = dirich;
      }
    }

  // count dof
  int counter = 0;
  for(auto &n : node)
  {
    if(n.dirich != nullptr)
      continue;

    for(int j = 0; j < fem.dofnp; j++)
    {
      n.dof[j] = counter;
      counter++;
    }
  }

#ifdef ICARAT_MPI
  if(fem.procid == 0)
#endif
    std::cout << "numeq : " << counter << std::endl;

  // set master node's dof to slave node's dof
  for(const auto &pai : pair)
    for(const auto &pa : pai)
    {
      int sec = pa.second;
      for(int j = 0; j < fem.dofnp; j++)
      {
        node[sec].dof[j] = node[pa.first].dof[j];
      }
      // Dirichlet objects will be deleted.
      node[sec].dirich.reset();
    }

  return counter;
}

template <class E, class N>
void BCIcarat<E, N>::input_neumann(const toml::value &config, Force &force)
{
  using neum_type =
      std::tuple<std::array<double, 3>, std::array<std::array<double, 2>, 3>>;
  double tol = 1.0e-7;

  const auto neum =
      toml::find<toml::array>(config, "boundary_condition", "neumann");
  int numNeum = (int)neum.size();

  for(auto &i_neum : neum)
  {
    auto data = toml::get<neum_type>(i_neum);

    auto val = std::get<0>(data);
    auto range = std::get<1>(data);
    IcaratLoad inp_f_tmp;
    inp_f_tmp.isDist = true;
    if(fem.ndim == 2 &&
       std::get<0>(std::get<0>(range)) == std::get<1>(std::get<0>(range)) &&
       std::get<0>(std::get<1>(range)) == std::get<1>(std::get<1>(range)))
      inp_f_tmp.isDist = false;
    if(fem.ndim == 3 &&
       std::get<0>(std::get<0>(range)) == std::get<1>(std::get<0>(range)) &&
       std::get<0>(std::get<1>(range)) == std::get<1>(std::get<1>(range)) &&
       std::get<0>(std::get<2>(range)) == std::get<1>(std::get<2>(range)))
      inp_f_tmp.isDist = false;

    inp_f_tmp.xmin = std::get<0>(std::get<0>(range)) - tol;
    inp_f_tmp.xmax = std::get<1>(std::get<0>(range)) + tol;
    inp_f_tmp.ymin = std::get<0>(std::get<1>(range)) - tol;
    inp_f_tmp.ymax = std::get<1>(std::get<1>(range)) + tol;
    inp_f_tmp.zmin = std::get<0>(std::get<2>(range)) - tol;
    inp_f_tmp.zmax = std::get<1>(std::get<2>(range)) + tol;
    inp_f_tmp.val = std::vector<double>(
        val.begin(), val.end()); // TODO: check this can be compiled
    inp_l.push_back(inp_f_tmp);
  }

#ifdef ICARAT_MPI
  if(fem.procid == 0)
#endif
    std::cout << "check point at neumann at macro        --->  ok" << std::endl;

  // make external force for "force.fext_org"
  generateLoadMain(numNeum, force);

#ifdef ICARAT_MPI
  if(fem.procid == 0)
#endif
    std::cout << "check point at external force at macro --->  ok" << std::endl;
}

/////////////////////////////////////
//////////////// private/////////////
/////////////////////////////////////

template <class E, class N>
void BCIcarat<E, N>::generateLoadMain(int numNeum, Force &force)
{
  std::vector<int> force_flag(fem.numnp);

  for(int neum = 0; neum < numNeum; neum++)
  {
    if(fem.ndim == 2)
    {
      if(inp_l[neum].isDist) // 分布荷重
      {
        for(int nel = 0; nel < fem.nelx; nel++)
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
          for(int k = 0; k < element[nel].ne; k++)
          {
            int tmp = element[nel].nodeID[k];
            // check node on this line
            if(inp_l[neum].xmin <= node[tmp].x[0] &&
               node[tmp].x[0] <= inp_l[neum].xmax &&
               inp_l[neum].ymin <= node[tmp].x[1] &&
               node[tmp].x[1] <= inp_l[neum].ymax)
            {
              nodecounter++;
              hitpoints.push_back(k);
              if(nodecounter == numnode)
              {
                inp_l[neum].sideID = surface.judgeLine(k, hitpoints);
                surface.generate(force, idof, X, inp_l[neum]);
              }
            }
          }
        }
      }
      else // 集中荷重
      {
        for(int np = 0; np < fem.numnp; np++)
        {
          if(inp_l[neum].xmin <= node[np].x[0] &&
             node[np].x[0] <= inp_l[neum].xmax &&
             inp_l[neum].ymin <= node[np].x[1] &&
             node[np].x[1] <= inp_l[neum].ymax)
          {
            for(int dn = 0; dn < fem.dofnp; dn++)
            {
              int dof = node[np].dof[dn];
              if(dof < fem.numeq)
                force.fext_org[dof] += inp_l[neum].val[dn];
            }
          }
        }
      }
    }
    else if(fem.ndim == 3)
    {
      if(inp_l[neum].isDist) // 分布荷重
      {
        for(int nel = 0; nel < fem.nelx; nel++)
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
          for(int k = 0; k < element[nel].ne; k++)
          {
            int tmp = element[nel].nodeID[k];
            // check node on this line
            if(inp_l[neum].xmin <= node[tmp].x[0] &&
               node[tmp].x[0] <= inp_l[neum].xmax &&
               inp_l[neum].ymin <= node[tmp].x[1] &&
               node[tmp].x[1] <= inp_l[neum].ymax &&
               inp_l[neum].zmin <= node[tmp].x[2] &&
               node[tmp].x[2] <= inp_l[neum].zmax)
            {
              nodecounter++;
              hitpoints.push_back(k);
              if(nodecounter == numnode)
              {
                inp_l[neum].sideID = surface.judgeLine(k, hitpoints);
                surface.generate(force, idof, X, inp_l[neum]);
              }
            }
          }
        }
      }
      else // 集中荷重
      {
        for(int np = 0; np < fem.numnp; np++)
        {
          if(inp_l[neum].xmin <= node[np].x[0] &&
             node[np].x[0] <= inp_l[neum].xmax &&
             inp_l[neum].ymin <= node[np].x[1] &&
             node[np].x[1] <= inp_l[neum].ymax &&
             inp_l[neum].zmin <= node[np].x[2] &&
             node[np].x[2] <= inp_l[neum].zmax)
          {
            for(int dn = 0; dn < fem.dofnp; dn++)
            {
              int dof = node[np].dof[dn];
              if(dof < fem.numeq)
                force.fext_org[dof] += inp_l[neum].val[dn];
            }
          }
        }
      }
    }
  }

  // // For the vtk's output, very small fe will be removed.
  // bool display = true;
  // for(int i = 0; i < (int)force.fext_org.size(); i++)
  //   if(abs(force.fext_org[i]) < 1.0e-10 && display == true)
  //   {
  //     force.fext_org[i] = 0.0;
  //     std::cerr << "Too small load's values will be skipped." << std::endl;
  //     display = false;
  //   }
}

} // namespace icarat