///  @file  bmatrix_shell.hpp
///  @author  Daiki Watanabe
///  @date  June 5, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php
#pragma once
#include <problem/base.hpp>
#include <vector>

namespace icarat
{
class BmatrixShell
{
public:
  /// main function of making B matrix for shell element
  void make(Eigen::MatrixXd &Be, Eigen::MatrixXd &Ne, Eigen::MatrixXd &P,
            double &jac, int ndim, int ne, int ipmax, Eigen::MatrixXd &X,
            Eigen::Vector3d &length, double thickness, int ip, int ipz);

  void square4gauss4(Eigen::MatrixXd &Be, Eigen::MatrixXd &Ne,
                     Eigen::MatrixXd &P, double &jac, Eigen::MatrixXd &X,
                     Eigen::Vector3d &length, double thickness, int ip,
                     int ipz);
};
} // namespace icarat