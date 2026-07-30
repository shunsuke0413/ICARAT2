///  @file  bmatrix_total.hpp
///  @author  Daiki Watanabe
///  @date  June 3, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php
#pragma once
#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <vector>

namespace icarat
{
class BmatrixTotal
{
public:
  /// main function of making B matrix and deformation gradient tensor F
  /// BL...linealization of green lagrange strain (DE[u])
  /// BNL...general B matrix(\nabla u)
  void make(Eigen::MatrixXd &BL, Eigen::MatrixXd &BNL, Eigen::MatrixXd &Ne,
            double &jac, Eigen::MatrixXd &F, std::string eType, int ipmax,
            Eigen::MatrixXd &X, Eigen::MatrixXd &U, int ip);

  void makeDifferential(Eigen::MatrixXd &dBLdu, Eigen::MatrixXd &dFdu,
                        double &dJdu, std::string eType, int ipmax,
                        Eigen::MatrixXd &X, Eigen::MatrixXd &U, int ip,
                        int dofid);
};

template <class S, class G>
void bmatrix2dTotal(Eigen::MatrixXd &BL, Eigen::MatrixXd &BNL,
                    Eigen::MatrixXd &Ne, double &jac, Eigen::MatrixXd &F,
                    Eigen::MatrixXd &X, Eigen::MatrixXd &U, int ip)
{
  S shape;
  G gauss;

  Eigen::VectorXd NeArr = shape.N(gauss.points(ip));
  Eigen::MatrixXd dXdr = shape.dNdr(gauss.points(ip)) * X; ///<(d,d)
  Eigen::MatrixXd dNdX =
      dXdr.inverse() * shape.dNdr(gauss.points(ip)); ///<(d,n)
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(shape.d, shape.d);

  ///@ F=dx/dX = du/dX + I
  Eigen::MatrixXd dudX = (dNdX * U).transpose();
  F = I + dudX;

  if(Ne.rows() == 1)
    for(int n = 0; n < shape.n; n++)
      Ne.coeffRef(n) = NeArr.coeff(n);
  else if(Ne.rows() == 2)
    for(int n = 0; n < shape.n; n++)
    {
      Ne.coeffRef(0, shape.d * n) = NeArr.coeff(n);
      Ne.coeffRef(1, shape.d * n + 1) = NeArr.coeff(n);
    }
  else
  {
    std::cerr << "Invalid row in bmatrix_updated.hpp" << std::endl;
    exit(1);
  }

  for(int n = 0; n < shape.n; n++)
  {
    // xx
    BL.coeffRef(0, 2 * n) = F.coeff(0, 0) * dNdX.coeff(0, n);
    BL.coeffRef(0, 2 * n + 1) = F.coeff(1, 0) * dNdX.coeff(0, n);

    // yy
    BL.coeffRef(1, 2 * n) = F.coeff(0, 1) * dNdX.coeff(1, n);
    BL.coeffRef(1, 2 * n + 1) = F.coeff(1, 1) * dNdX.coeff(1, n);

    // xy
    BL.coeffRef(2, 2 * n) =
        F.coeff(0, 1) * dNdX.coeff(0, n) + F.coeff(0, 0) * dNdX.coeff(1, n);
    BL.coeffRef(2, 2 * n + 1) =
        F.coeff(1, 1) * dNdX.coeff(0, n) + F.coeff(1, 0) * dNdX.coeff(1, n);

    BNL.coeffRef(0, 2 * n) = dNdX.coeff(0, n);
    BNL.coeffRef(1, 2 * n + 1) = dNdX.coeff(0, n);
    BNL.coeffRef(2, 2 * n) = dNdX.coeff(1, n);
    BNL.coeffRef(3, 2 * n + 1) = dNdX.coeff(1, n);
  }

  // jacobian * gauss weight
  jac = dXdr.determinant() * gauss.weights(ip);
}

template <class S, class G>
void dBdu2dTotal(Eigen::MatrixXd &dBLdu, Eigen::MatrixXd &dFdu, double &dJdu,
                 Eigen::MatrixXd &X, Eigen::MatrixXd &U, int ip, int dofid)
{
  S shape;
  G gauss;

  Eigen::MatrixXd dXdr = shape.dNdr(gauss.points(ip)) * X; ///<(d,d)
  Eigen::MatrixXd dNdX =
      dXdr.inverse() * shape.dNdr(gauss.points(ip)); ///<(d,n)

  // F=dx/dX = du/dX + I
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(shape.d, shape.d);
  Eigen::MatrixXd dudX = (dNdX * U).transpose();
  Eigen::MatrixXd F = I + dudX;

  // dBLdu
  for(int n = 0; n < shape.n; n++)
  {
    int node_number = dofid / 2;
    if(dofid % 2 == 0) // x
    {
      for(int n = 0; n < shape.n; n++)
      {
        // xx
        dBLdu(0, 2 * n) = dNdX.coeff(0, node_number) * dNdX.coeff(0, n);
        dBLdu(0, 2 * n + 1) = 0.0;
        // yy
        dBLdu(1, 2 * n) = dNdX.coeff(1, node_number) * dNdX.coeff(1, n);
        dBLdu(1, 2 * n + 1) = 0.0;
        // xy
        dBLdu(2, 2 * n) = dNdX.coeff(0, n) * dNdX.coeff(1, node_number) +
                          dNdX.coeff(0, node_number) * dNdX.coeff(1, n);
        dBLdu(2, 2 * n + 1) = 0.0;
      }
    }
    else // y
    {
      for(int n = 0; n < shape.n; n++)
      {
        // xx
        dBLdu(0, 2 * n) = 0.0;
        dBLdu(0, 2 * n + 1) = dNdX.coeff(0, n) * dNdX.coeff(0, node_number);
        // yy
        dBLdu(1, 2 * n) = 0.0;
        dBLdu(1, 2 * n + 1) = dNdX.coeff(1, n) * dNdX.coeff(1, node_number);
        // xy
        dBLdu(2, 2 * n) = 0.0;
        dBLdu(2, 2 * n + 1) = dNdX.coeff(0, n) * dNdX.coeff(1, node_number) +
                              dNdX.coeff(0, node_number) * dNdX.coeff(1, n);
      }
    }
  }

  int node_number = dofid / 2;
  int xy = dofid % 2;

  // dFdu
  if(xy == 0)
    dFdu << dNdX.coeff(0, node_number), dNdX.coeff(1, node_number), 0.0, 0.0;
  else if(xy == 1)
    dFdu << 0.0, 0.0, dNdX.coeff(0, node_number), dNdX.coeff(1, node_number);
  else
  {
    std::cerr << "error in bmatrix" << std::endl;
    exit(1);
  }

  // dJdu
  if(xy == 0)
    dJdu = dNdX.coeff(0, node_number) * F(1, 1) -
           dNdX.coeff(1, node_number) * F(1, 0);
  else if(xy == 1)
    dJdu = dNdX.coeff(1, node_number) * F(0, 0) -
           dNdX.coeff(0, node_number) * F(0, 1);
  else
  {
    std::cerr << "error in bmatrix" << std::endl;
    exit(1);
  }
}

template <class S, class G>
void bmatrix3dTotal(Eigen::MatrixXd &BL, Eigen::MatrixXd &BNL,
                    Eigen::MatrixXd &Ne, double &jac, Eigen::MatrixXd &F,
                    Eigen::MatrixXd &X, Eigen::MatrixXd &U, int ip)
{
  S shape;
  G gauss;

  Eigen::VectorXd NeArr = shape.N(gauss.points(ip));
  Eigen::MatrixXd dXdr = shape.dNdr(gauss.points(ip)) * X; ///<(d,d)
  Eigen::MatrixXd dNdX =
      dXdr.inverse() * shape.dNdr(gauss.points(ip)); ///<(d,n)
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(shape.d, shape.d);

  ///@ F=dx/dX = du/dX + I
  F = I + (dNdX * U).transpose();

  if(Ne.rows() == 1)
    for(int n = 0; n < shape.n; n++)
      Ne.coeffRef(n) = NeArr.coeff(n);
  else if(Ne.rows() == 3)
    for(int n = 0; n < shape.n; n++)
    {
      Ne.coeffRef(0, shape.d * n) = NeArr.coeff(n);
      Ne.coeffRef(1, shape.d * n + 1) = NeArr.coeff(n);
      Ne.coeffRef(2, shape.d * n + 2) = NeArr.coeff(n);
    }
  else
  {
    std::cerr << "Invalid row in bmatrix_updated.hpp" << std::endl;
    exit(1);
  }

  for(int n = 0; n < shape.n; n++)
  {
    // xx
    BL.coeffRef(0, 3 * n) = F.coeff(0, 0) * dNdX.coeff(0, n);
    BL.coeffRef(0, 3 * n + 1) = F.coeff(1, 0) * dNdX.coeff(0, n);
    BL.coeffRef(0, 3 * n + 2) = F.coeff(2, 0) * dNdX.coeff(0, n);

    // yy
    BL.coeffRef(1, 3 * n) = F.coeff(0, 1) * dNdX.coeff(1, n);
    BL.coeffRef(1, 3 * n + 1) = F.coeff(1, 1) * dNdX.coeff(1, n);
    BL.coeffRef(1, 3 * n + 2) = F.coeff(2, 1) * dNdX.coeff(1, n);

    // zz
    BL.coeffRef(2, 3 * n) = F.coeff(0, 2) * dNdX.coeff(2, n);
    BL.coeffRef(2, 3 * n + 1) = F.coeff(1, 2) * dNdX.coeff(2, n);
    BL.coeffRef(2, 3 * n + 2) = F.coeff(2, 2) * dNdX.coeff(2, n);

    // xy
    BL.coeffRef(3, 3 * n) =
        F.coeff(0, 1) * dNdX.coeff(0, n) + F.coeff(0, 0) * dNdX.coeff(1, n);
    BL.coeffRef(3, 3 * n + 1) =
        F.coeff(1, 1) * dNdX.coeff(0, n) + F.coeff(1, 0) * dNdX.coeff(1, n);
    BL.coeffRef(3, 3 * n + 2) =
        F.coeff(2, 1) * dNdX.coeff(0, n) + F.coeff(2, 0) * dNdX.coeff(1, n);

    // yz
    BL.coeffRef(4, 3 * n) =
        F.coeff(0, 2) * dNdX.coeff(1, n) + F.coeff(0, 1) * dNdX.coeff(2, n);
    BL.coeffRef(4, 3 * n + 1) =
        F.coeff(1, 2) * dNdX.coeff(1, n) + F.coeff(1, 1) * dNdX.coeff(2, n);
    BL.coeffRef(4, 3 * n + 2) =
        F.coeff(2, 2) * dNdX.coeff(1, n) + F.coeff(2, 1) * dNdX.coeff(2, n);

    // zx
    BL.coeffRef(5, 3 * n) =
        F.coeff(0, 0) * dNdX.coeff(2, n) + F.coeff(0, 2) * dNdX.coeff(0, n);
    BL.coeffRef(5, 3 * n + 1) =
        F.coeff(1, 0) * dNdX.coeff(2, n) + F.coeff(1, 2) * dNdX.coeff(0, n);
    BL.coeffRef(5, 3 * n + 2) =
        F.coeff(2, 0) * dNdX.coeff(2, n) + F.coeff(2, 2) * dNdX.coeff(0, n);

    BNL.coeffRef(0, 3 * n) = dNdX.coeff(0, n);
    BNL.coeffRef(1, 3 * n + 1) = dNdX.coeff(0, n);
    BNL.coeffRef(2, 3 * n + 2) = dNdX.coeff(0, n);

    BNL.coeffRef(3, 3 * n) = dNdX.coeff(1, n);
    BNL.coeffRef(4, 3 * n + 1) = dNdX.coeff(1, n);
    BNL.coeffRef(5, 3 * n + 2) = dNdX.coeff(1, n);

    BNL.coeffRef(6, 3 * n) = dNdX.coeff(2, n);
    BNL.coeffRef(7, 3 * n + 1) = dNdX.coeff(2, n);
    BNL.coeffRef(8, 3 * n + 2) = dNdX.coeff(2, n);
  }

  // jacobian * gauss weight
  jac = dXdr.determinant() * gauss.weights(ip);
}

template <class S, class G>
void dBdu3dTotal(Eigen::MatrixXd &dBLdu, Eigen::MatrixXd &dFdu, double &dJdu,
                 Eigen::MatrixXd &X, Eigen::MatrixXd &U, int ip, int dofid)
{
  S shape;
  G gauss;

  Eigen::VectorXd NeArr = shape.N(gauss.points(ip));
  Eigen::MatrixXd dXdr = shape.dNdr(gauss.points(ip)) * X; ///<(d,d)
  Eigen::MatrixXd dNdX =
      dXdr.inverse() * shape.dNdr(gauss.points(ip)); ///<(d,n)
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(shape.d, shape.d);

  ///@ F=dx/dX = du/dX + I
  Eigen::MatrixXd dudX = (dNdX * U).transpose();
  Eigen::MatrixXd F = I + dudX;

  for(int n = 0; n < shape.n; n++)
  {
    int node_number = dofid / 3;
    if(dofid % 3 == 0)
    {
      for(int n = 0; n < shape.n; n++)
      {
        // xx
        dBLdu(0, 3 * n) = dNdX.coeff(0, node_number) * dNdX.coeff(0, n);
        dBLdu(0, 3 * n + 1) = 0.0;
        dBLdu(0, 3 * n + 2) = 0.0;
        // yy
        dBLdu(1, 3 * n) = dNdX.coeff(1, node_number) * dNdX.coeff(1, n);
        dBLdu(1, 3 * n + 1) = 0.0;
        dBLdu(1, 3 * n + 2) = 0.0;
        // xy
        dBLdu(2, 3 * n) = dNdX.coeff(0, n) * dNdX.coeff(1, node_number) +
                          dNdX.coeff(0, node_number) * dNdX.coeff(1, n);
        dBLdu(2, 3 * n + 1) = 0.0;
        dBLdu(2, 3 * n + 2) = 0.0;
      }
    }
    else
    {
      for(int n = 0; n < shape.n; n++)
      {
        // xx
        dBLdu(0, 2 * n) = 0.0;
        dBLdu(0, 2 * n + 1) = dNdX.coeff(0, n) * dNdX.coeff(0, node_number);
        // yy
        dBLdu(1, 2 * n) = 0.0;
        dBLdu(1, 2 * n + 1) = dNdX.coeff(1, n) * dNdX.coeff(1, node_number);
        // xy
        dBLdu(2, 2 * n) = 0.0;
        dBLdu(2, 2 * n + 1) = dNdX.coeff(0, n) * dNdX.coeff(1, node_number) +
                              dNdX.coeff(0, node_number) * dNdX.coeff(1, n);
      }
    }
  }
}

} // namespace icarat