///  @file  eigen_manipulation.cpp
///  @author  Daiki Watanabe
///  @date  November 20, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once

#include <Eigen/Core>
#include <problem/base.hpp>
#include <unsupported/Eigen/CXX11/Tensor>

namespace icarat
{
/// T...Eigen::Triplet<double>
/// assembling Ke & fe
template <typename T>
void assembling(int numdof, int numeq, Eigen::VectorXi &idof,
                Eigen::MatrixXd &Ke, Eigen::VectorXd &fe, std::vector<T> &tripK,
                std::vector<T> &tripF)
{
  for(int i = 0; i < numdof; i++)
  {
    // boundary condition check
    if(idof.coeff(i) >= numeq)
      continue;

    for(int j = 0; j < numdof; j++)
    {
      // boundary condition check
      if(idof.coeff(j) >= numeq)
        continue;

      tripK.push_back(T(idof.coeff(i), idof.coeff(j), Ke.coeff(i, j)));
    }
    tripF.push_back(T(idof.coeff(i), idof.coeff(i), fe.coeff(i)));
  }
}

template <typename T>
void assemblingFull(int numdof, int total_neq, const Eigen::VectorXi &idof,
                    const Eigen::MatrixXd &Ke, std::vector<T> &tripKfull)
{
  for(int i = 0; i < numdof; i++)
  {
    int global_row = idof.coeff(i);

    if(global_row < 0)
      continue;

    for(int j = 0; j < numdof; j++)
    {
      int global_col = idof.coeff(j);

      if(global_col < 0)
        continue;

      tripKfull.push_back(T(global_row, global_col, Ke.coeff(i, j)));
    }
  }
}

/// assembling only Ke
template <typename T>
void assembling(int numdof, int numeq, Eigen::VectorXi &idof,
                Eigen::MatrixXd &Ke, std::vector<T> &tripK)
{
  for(int i = 0; i < numdof; i++)
  {
    tripK.reserve(tripK.size() + numdof * numdof);

    // boundary condition check
    if(idof.coeff(i) >= numeq)
      continue;
    for(int j = 0; j < numdof; j++)
    {
      // boundary condition check
      if(idof.coeff(j) >= numeq)
        continue;
      tripK.push_back(T(idof.coeff(i), idof.coeff(j), Ke.coeff(i, j)));
    }
  }
}

/// assembling only fe
template <typename T>
void assembling(int numdof, int numeq, Eigen::VectorXi &idof,
                Eigen::VectorXd &fe, std::vector<T> &tripF)
{
  for(int i = 0; i < numdof; i++)
  {
    // boundary condition check
    if(idof.coeff(i) >= numeq)
      continue;

    tripF.push_back(T(idof.coeff(i), idof.coeff(i), fe.coeff(i)));
  }
}

/// assembling reaction force
template <typename T>
void assemblingReact(int numdof, Eigen::VectorXi &idof, Eigen::VectorXd &fe,
                     std::vector<T> &tripR)
{
  for(int i = 0; i < numdof; i++)
    tripR.push_back(T(idof.coeff(i), idof.coeff(i), fe.coeff(i)));
}

template <class T, class N>
void assemblingHomo(int ne, int numdof, Eigen::VectorXi &idof,
                    Eigen::MatrixXd &Ke, Eigen::VectorXd &fe,
                    std::vector<T> &tripletsK, std::vector<T> &tripletsF,
                    FEM &fem, std::vector<int> nodeID, std::vector<N> &nodes)
{
  // make cdof
  std::vector<std::vector<int>> cdof; ///< [numdof][numCnt]
  for(int i = 0; i < ne; i++)
  {
    // if not control point...cdof=idof
    if(nodes[nodeID[i]].mpc == nullptr)
    {
      for(int j = 0; j < fem.dofnp; j++)
      {
        std::vector<int> tmp(1, idof[fem.dofnp * i + j]);
        cdof.push_back(tmp);
      }
    }
    // if control point...cdof=idof + control point's dof
    else
    {
      for(int j = 0; j < fem.dofnp; j++)
      {
        std::vector<int> tmp = nodes[nodeID[i]].mpc->cdof[j];
        tmp.push_back(idof[fem.dofnp * i + j]);
        cdof.push_back(tmp);
      }
    }
  }

  // assembling
  for(int i = 0; i < numdof; i++)
  {
    for(int k = 0; k < (int)cdof[i].size(); k++)
    {
      for(int j = 0; j < numdof; j++)
      {
        for(int l = 0; l < (int)cdof[j].size(); l++)
        {
          T tripletK(cdof[i][k], cdof[j][l], Ke.coeff(i, j));
          tripletsK.push_back(tripletK);
        }
      }
      T tripletF(cdof[i][k], cdof[i][k], fe.coeff(i));
      tripletsF.push_back(tripletF);
    }
  }
}

template <class T, class N>
void assemblingReactHomo(int ne, int numdof, Eigen::VectorXi &idof,
                         Eigen::VectorXd &fe, std::vector<T> &tripletsR,
                         FEM &fem, std::vector<int> nodeID,
                         std::vector<N> &nodes)
{
  // make cdof
  std::vector<std::vector<int>> cdof; ///< [numdof][numCnt]
  for(int i = 0; i < ne; i++)
  {
    // if not control point...cdof=idof
    if(nodes[nodeID[i]].mpc == nullptr)
    {
      for(int j = 0; j < fem.dofnp; j++)
      {
        std::vector<int> tmp(1, idof[fem.dofnp * i + j]);
        cdof.push_back(tmp);
      }
    }
    // if control point...cdof=idof + control point's dof
    else
    {
      for(int j = 0; j < fem.dofnp; j++)
      {
        std::vector<int> tmp = nodes[nodeID[i]].mpc->cdof[j];
        tmp.push_back(idof[fem.dofnp * i + j]);
        cdof.push_back(tmp);
      }
    }
  }

  // assembling
  for(int i = 0; i < numdof; i++)
  {
    for(int k = 0; k < (int)cdof[i].size(); k++)
    {
      T tripletR(cdof[i][k], cdof[i][k], fe.coeff(i));
      tripletsR.push_back(tripletR);
    }
  }
}

#ifdef ICARAT_MPI
/// assembling Ke & fe used in MPI
/// T...Eigen::Triplet<double>
template <typename T>
void assembling(int numdof, Eigen::VectorXi &idof, Eigen::VectorXi &jdof,
                std::vector<bool> &isMasters, std::vector<int> &mprocids,
                Eigen::MatrixXd &Ke, Eigen::VectorXd &fe, std::vector<T> &Kp,
                std::vector<std::pair<int, T>> &Ks, std::vector<T> &Fp,
                std::vector<std::pair<int, T>> &Fs)
{
  for(int i = 0; i < numdof; i++) // row (local) index
  {
    // boundary condition check
    if(idof.coeff(i) < 0)
      continue;

    for(int j = 0; j < numdof; j++) // column (global) index
    {
      // boundary condition check
      if(jdof.coeff(j) < 0)
        continue;

      // slave
      if(!isMasters[i])
      {
        std::pair<int, T> pair;
        pair.first = mprocids[i];
        pair.second = T(idof.coeff(i), jdof.coeff(j), Ke.coeff(i, j));
        Ks.push_back(pair);
      }
      // master
      else
      {
        Kp.push_back(T(idof.coeff(i), jdof.coeff(j), Ke.coeff(i, j)));
      }
    }

    // slave
    if(!isMasters[i])
    {
      std::pair<int, T> pair;
      pair.first = mprocids[i];
      pair.second = T(idof.coeff(i), idof.coeff(i), fe.coeff(i));
      Fs.push_back(pair);
    }
    // master
    else
    {
      Fp.push_back(T(idof.coeff(i), idof.coeff(i), fe.coeff(i)));
    }
  }
}

/// assembling Ke used in MPI
/// T...Eigen::Triplet<double>
template <typename T>
void assembling(int numdof, Eigen::VectorXi &idof, Eigen::VectorXi &jdof,
                std::vector<bool> &isMasters, std::vector<int> &mprocids,
                Eigen::MatrixXd &Ke, std::vector<T> &Kp,
                std::vector<std::pair<int, T>> &Ks)
{
  for(int i = 0; i < numdof; i++) // row (local) index
  {
    // boundary condition check
    if(idof.coeff(i) < 0)
      continue;

    for(int j = 0; j < numdof; j++) // column (global) index
    {
      // boundary condition check
      if(jdof.coeff(j) < 0)
        continue;

      // slave
      if(!isMasters[i])
      {
        std::pair<int, T> pair;
        pair.first = mprocids[i];
        pair.second = T(idof.coeff(i), jdof.coeff(j), Ke.coeff(i, j));
        Ks.push_back(pair);
      }
      // master
      else
      {
        Kp.push_back(T(idof.coeff(i), jdof.coeff(j), Ke.coeff(i, j)));
      }
    }
  }
}
#endif

/// This function converts numeq values to node values within the element
/// considering dirichlet condition
template <class E, class N>
std::vector<Eigen::VectorXd>
numeqToElement(const FEM &fem, const std::vector<E> &element,
               const std::vector<N> &node, const Eigen::VectorXd &numeqValues)
{
  std::vector<Eigen::VectorXd> elementValues(fem.nelx);

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    int counter = 0;
    Eigen::VectorXd tmp(element[nel].numdof);
    for(int i = 0; i < element[nel].ne; i++)
    {
      for(int j = 0; j < fem.dofnp; j++)
      {
        int dof = node[element[nel].nodeID[i]].dof[j];
        // boundary conditioin check
        if(dof < fem.numeq)
          tmp.coeffRef(counter) = numeqValues.coeff(dof);
        else
          tmp.coeffRef(counter) = node[element[nel].nodeID[i]].dirich->val[j];

        counter++;
      }
    }
    elementValues[nel] = tmp;
  }
  return elementValues;
}

/// convert element values to global(numeq) values with Dirichlet conditiion
template <class E, class N>
Eigen::VectorXd
elementToNumeq(const FEM &fem, const std::vector<E> &element,
               const std::vector<N> &node,
               const std::vector<Eigen::VectorXd> &elementValues)
{
  Eigen::VectorXd numeqValues(fem.numeq);
  numeqValues.setZero();

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    int counter = 0;
    Eigen::VectorXi idof(element[nel].numdof);
    for(int i = 0; i < element[nel].ne; i++)
    {
      for(int j = 0; j < fem.dofnp; j++)
      {
        idof.coeffRef(counter) = node[element[nel].nodeID[i]].dof[j];
        counter++;
      }
    }

    for(int i = 0; i < counter; i++)
    {
      // boundary conditioin check
      if(idof.coeff(i) < fem.numeq)
        numeqValues[idof.coeff(i)] += elementValues[nel].coeff(i);
      else
        continue;
    }
  }
  return numeqValues;
}

} // namespace icarat