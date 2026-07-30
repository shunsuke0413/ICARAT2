///  @file  filter.hpp
///  @author  Daiki Watanabe
///  @date  November 15, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include "opt_base.hpp"
#include <problem/base.hpp>
#include <vector>

namespace icarat
{
/// Density or Sensitivity filter class in TO
/// This class can be only applied for uniform cube or square mesh
template <class E, class N> class Filter
{
public:
  /// @param [in] radius_ :filter radius
  Filter(FEM &fem_, std::vector<E> &element_, std::vector<N> &node_,
         double radius_);

  ///  apply sensitivity Filter
  std::vector<double> sensitivityFilter(const std::vector<double> &design_s,
                                        const std::vector<double> &dfds);

  /// apply density filter
  /// @returns filtered design_s and sensitivity derivative term of design_s
  /// @param  [in] isSens true...sensitivity derivative term,
  /// false...density filter procedure
  std::vector<double> densityFilter(const std::vector<double> &design_s,
                                    bool isSens);

  std::vector<double> densityFilter2(const std::vector<double> &design_s,
                                    bool isSens, int optstep, double mesh_length);

  /// update design variable using threshold function
  std::vector<double> threshold(const std::vector<double> &design_s,
                                double beta, double T);

  /// differential derivative of threshold function
  std::vector<double> thresholdDerivative(const std::vector<double> &design_s,
                                          double beta, double T);

private:
  /// make numele,ID_,dis_,dfilters_ds
  void preNeighbour();

  FEM &fem;
  std::vector<E> &element;
  std::vector<N> &node;
  double radius_;           ///< filter radius
  std::vector<int> numele_; //< number of the adjacent element for an element i
  std::vector<std::vector<int>>
      ID_; ///< ID[i][j] of the adjacent element j for an element i
  std::vector<std::vector<double>>
      dis_; ///< distance[i][j] from an element i to the adjacent element j
  std::vector<double> dfilters_ds; ///< it is used in "densityFilterSensitivity"
};

template <class E, class N>
Filter<E, N>::Filter(FEM &fem_, std::vector<E> &element_, std::vector<N> &node_,
                     double radius)
    : fem(fem_), element(element_), node(node_), radius_(radius),
      numele_(fem.nelx), ID_(fem.nelx), dis_(fem.nelx), dfilters_ds(fem.nelx)
{
  preNeighbour();
}

///////////////////////////////////////////////////////////
////////////////////////////public////////////////////////
//////////////////////////////////////////////////////////

template <class E, class N>
std::vector<double>
Filter<E, N>::sensitivityFilter(const std::vector<double> &design_s,
                                const std::vector<double> &dfds)
{
  std::vector<double> dfds_new(dfds.size()); ///< non-filtered sensitivity

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < fem.nelx; i++)
  {
    double aa = 0.0;
    double bb = 0.0;

    for(int j = 0; j < numele_[i]; j++)
    {
      double weight = radius_ - dis_[i][j];
      aa += weight;
      bb += (weight * design_s[ID_[i][j]] * dfds[ID_[i][j]]) /
            (design_s[ID_[i][j]] * element[ID_[i][j]].volume);
    }

    double cc = (design_s[i] * element[i].volume) / (design_s[i] * aa);
    dfds_new[i] = cc * bb;
  }

  return dfds_new;
}

template <class E, class N>
std::vector<double>
Filter<E, N>::densityFilter(const std::vector<double> &design_s, bool isSens)
{
  std::vector<double> filtered_s(fem.nelx);

  // density filter
  if(isSens == false)
  {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < fem.nelx; i++)
    {
      double aa = 0.0;
      double bb = 0.0;

      for(int j = 0; j < numele_[i]; j++)
      {
        double weight = radius_ - dis_[i][j];
        aa += weight;
        bb += weight * design_s[ID_[i][j]];
      }

      if(aa < 1.0e-10)
      {
        std::cerr << "coefficient aa is too small in filter.hpp." << std::endl;
        exit(1);
      }
      filtered_s[i] = bb / aa;
      dfilters_ds[i] = aa;
    }
  }
  // density filter's sensitivity
  else
  {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < (int)design_s.size(); i++)
    {
      for(int j = 0; j < numele_[i]; j++)
      {
        double weight = radius_ - dis_[i][j];
        filtered_s[i] += weight * design_s[ID_[i][j]] / dfilters_ds[ID_[i][j]];
      }
    }
  }
  return filtered_s;
}

template <class E, class N>
std::vector<double>
Filter<E, N>::densityFilter2(const std::vector<double> &design_s, bool isSens, int optstep, double mesh_length)
{
  std::vector<double> filtered_s(fem.nelx);
  double decay_amount = std::max((optstep - 46) * 0.05, 0.0);
  double radius_now = std::max(radius_ - decay_amount, mesh_length);

  // density filter
  if(isSens == false)
  {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < fem.nelx; i++)
    {
      double aa = 0.0;
      double bb = 0.0;

      for(int j = 0; j < numele_[i]; j++)
      {
        double weight = radius_now - dis_[i][j];
        aa += weight;
        bb += weight * design_s[ID_[i][j]];
      }

      if(aa < 1.0e-10)
      {
        std::cerr << "coefficient aa is too small in filter.hpp." << std::endl;
        exit(1);
      }
      filtered_s[i] = bb / aa;
      dfilters_ds[i] = aa;
    }
  }
  // density filter's sensitivity
  else
  {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < (int)design_s.size(); i++)
    {
      for(int j = 0; j < numele_[i]; j++)
      {
        double weight = radius_now - dis_[i][j];
        filtered_s[i] += weight * design_s[ID_[i][j]] / dfilters_ds[ID_[i][j]];
      }
    }
  }
  return filtered_s;
}

template <class E, class N>
std::vector<double>
Filter<E, N>::threshold(const std::vector<double> &d_filtered, double beta,
                        double T)
{
  std::vector<double> hat((int)d_filtered.size());

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < (int)hat.size(); i++)
    hat[i] = 0.5 * (tanh(T * beta) + tanh(beta * (d_filtered[i] - T))) /
             tanh(T * beta);

  return hat;
}

template <class E, class N>
std::vector<double>
Filter<E, N>::thresholdDerivative(const std::vector<double> &d_fil, double beta,
                                  double T)
{
  std::vector<double> dshat((int)d_fil.size());

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < (int)dshat.size(); i++)
    dshat[i] = (0.5 * (beta * (1.0 - pow(tanh(beta * (d_fil[i] - T)), 2.0))) /
                tanh(T * beta));

  return dshat;
}

///////////////////////////////////////////////////////////
////////////////////////////private////////////////////////
///////////////////////////////////////////////////////////

template <class E, class N> void Filter<E, N>::preNeighbour()
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < fem.nelx; i++)
  {
    int nei = element[i].ne;
    // coordinate of element's center
    Eigen::VectorXd centerI(fem.ndim);
    centerI.setZero();
    for(int dim = 0; dim < fem.ndim; dim++)
      for(int e = 0; e < nei; e++)
        centerI[dim] += node[element[i].nodeID[e]].x[dim] / nei;

    std::vector<int> IDi;
    std::vector<double> disi;
    int counter = 0;
    for(int j = 0; j < fem.nelx; j++)
    {
      int nej = element[j].ne;
      Eigen::VectorXd centerJ(fem.ndim);
      centerJ.setZero();

      for(int dim = 0; dim < fem.ndim; dim++)
        for(int e = 0; e < nej; e++)
          centerJ[dim] += node[element[j].nodeID[e]].x[dim] / nej;

      double disTmp = (centerI - centerJ).norm();

      if(disTmp <= radius_)
      {
        IDi.push_back(j);
        disi.push_back(disTmp);
        counter++;
      }
    } // j loop end

    // stock filter data of element i
    numele_[i] = counter;
    ID_[i] = IDi;
    dis_[i] = disi;
  }
}

} // namespace icarat
