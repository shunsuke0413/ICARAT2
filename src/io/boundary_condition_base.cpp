///  @file  boundary_condition_base.cpp
///  @author  Daiki Watanabe
///  @date  February 9, 2022.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "boundary_condition_base.hpp"
#include <basis/nmatrix_surface.hpp>

namespace icarat
{
using namespace std;
using namespace Eigen;

int SurfaceLoad::judgeLine(int counter, std::vector<int> hitpoints)
{
  int side;
  if(actele_.eType == "tria3" || actele_.eType == "tria6")
  {
    if(counter == 1)
      side = 1; // Ymin
    else
    {
      if(hitpoints[1] - hitpoints[0] == 1)
        side = 2; // Ymax
      else
        side = 3; // Xmin
    }
  }
  else if(actele_.eType == "quad4" || actele_.eType == "quad8")
  {
    if(counter == 1)
      side = 1; // Ymin
    else if(counter == 2)
      side = 2; // Xmax
    else
    {
      if(hitpoints[1] - hitpoints[0] == 1)
        side = 3; // Ymax
      else
        side = 4; // Xmin
    }
  }
  else if(actele_.eType == "tetra4" || actele_.eType == "tria10")
  {
    if(counter == 2)
      side = 1;
    else
    {
      if(hitpoints[2] + hitpoints[1] + hitpoints[0] == 4)
        side = 2;
      else if(hitpoints[2] + hitpoints[1] + hitpoints[0] == 6)
        side = 3;
      else
        side = 4;
    }
  }
  else if(actele_.eType == "hexa8" || actele_.eType == "hexa20")
  {
    if(counter == 3)
      side = 1; // Zmin
    else if(counter == 5)
      side = 3; // Ymin
    else if(counter == 6)
      side = 4; // Xmax
    else
    {
      if(hitpoints[2] - hitpoints[1] == 1)
        side = 2; // Xmin
      else if(hitpoints[1] - hitpoints[0] == 1)
        side = 5; // Ymax
      else
        side = 6; // Zmax
    }
  }
  else
  {
    cerr << "No load can be applied to this element type in "
            "boundary_condition_base.cpp "
         << endl;
    exit(1);
  }

  return side;
}

int SurfaceLoad::setNumNode()
{

  int numnode;
  if(actele_.eType == "tria3" || actele_.eType == "tria6")
    numnode = 2;
  else if(actele_.eType == "quad4" || actele_.eType == "quad8")
    numnode = 2;
  else if(actele_.eType == "tetra4" || actele_.eType == "tetra10")
    numnode = 3;
  else if(actele_.eType == "hexa8" || actele_.eType == "hexa20")
    numnode = 4;
  else
  {
    cerr << "No load can be applied to this element type in "
            "boundary_condition_base.cpp "
         << endl;
    exit(1);
  }

  return numnode;
}

int SurfaceLoad::surfaceIpmax()
{
  // set number of gauss point for appried surface
  int ipmax;
  // line2
  if(actele_.eType == "tria3" || actele_.eType == "quad4" ||
     actele_.eType == "mitc4")
    ipmax = 2;
  // line3
  else if(actele_.eType == "quad8" || actele_.eType == "tria6")
    ipmax = 3;
  // tria3
  else if(actele_.eType == "tetra4")
    ipmax = 1;
  // tria6
  else if(actele_.eType == "tetra10")
    ipmax = 3;
  // quad4
  else if(actele_.eType == "hexa8")
    ipmax = 4;
  // quad8
  else if(actele_.eType == "hexa20")
    ipmax = 9;
  else
  {
    std::cerr << "Error in boundary_condition_base.hpp" << std::endl;
    exit(1);
  }

  return ipmax;
}

void SurfaceLoad::generate(Force &force, VectorXi &idof, MatrixXd &X,
                           InputLoad &inp_f)
{
  // Convert val to transform global coordinates
  int min = std::min(3, fem_.dofnp);
  VectorXd val = VectorXd::Zero(min);
  for(int i = 0; i < min; i++)
    val[i] = inp_f.val[i];

  if(inp_f.type == "pressure")
    val = transLocalToGlobalPressure(X, val, inp_f.sideID);

  // get surface parameters
  string eTypeSur;
  int neSur, ipmaxSur;
  tie(eTypeSur, neSur, ipmaxSur) =
      resizeSurfaceParameter(idof, X, actele_.eType, inp_f.sideID);

  if(inp_f.type == "moment")
    val = transLocalToGlobalMoment(X, val, inp_f);

  // generate equivalent nodal load
  MatrixXd Ne = MatrixXd::Zero(min, idof.size());
  VectorXd fe = VectorXd::Zero(idof.size());
  SurfaceNmatrix nmatrix;
  for(int ip = 0; ip < ipmaxSur; ip++)
  {
    double jac = 0.0;

    // make shape function for surface
    nmatrix.make(Ne, jac, eTypeSur, ipmaxSur, X, ip);

    fe += Ne.transpose() * val * jac;
  }

  // assembling
  for(int i = 0; i < idof.size(); i++)
  {
    // check boundary condition
    if(idof.coeff(i) >= fem_.numeq)
      continue;

    force.fext_org(idof[i]) += fe[i];
  }
}

std::tuple<std::string, int, int>
SurfaceLoad::resizeSurfaceParameter(VectorXi &idof, MatrixXd &X,
                                    std::string eType, int sideID)
{
  int ne, ipmax;
  string eTypeSur; ///< element type of surface
  VectorXi neID;

  if(eType == "tria3")
  {
    eTypeSur = "line2";
    ne = 2;
    ipmax = 2;
    neID.resize(ne);
    if(sideID == 1)
      neID << 0, 1;
    else if(sideID == 2)
      neID << 1, 2;
    else if(sideID == 3)
      neID << 2, 0;
    else
    {
      cerr << "sideID is invalid in boundary_condition_base.cpp" << endl;
      exit(1);
    }
  }
  else if(eType == "tria6")
  {
    eTypeSur = "line3";
    ne = 3;
    ipmax = 3;
    neID.resize(ne);
    if(sideID == 1)
      neID << 0, 1, 3;
    else if(sideID == 2)
      neID << 1, 2, 4;
    else if(sideID == 3)
      neID << 2, 0, 5;
    else
    {
      cerr << "sideID is invalid in boundary_condition_base.cpp" << endl;
      exit(1);
    }
  }
  else if(eType == "quad4")
  {
    eTypeSur = "line2";
    ne = 2;
    ipmax = 2;
    neID.resize(ne);
    if(sideID == 1)
      neID << 0, 1;
    else if(sideID == 2)
      neID << 1, 2;
    else if(sideID == 3)
      neID << 2, 3;
    else if(sideID == 4)
      neID << 3, 0;
    else
    {
      cerr << "sideID is invalid in boundary_condition_base.cpp" << endl;
      exit(1);
    }
  }
  else if(eType == "quad8")
  {
    eTypeSur = "line3";
    ne = 3;
    ipmax = 3;
    neID.resize(ne);
    if(sideID == 1)
      neID << 0, 1, 4;
    else if(sideID == 2)
      neID << 1, 2, 5;
    else if(sideID == 3)
      neID << 2, 3, 6;
    else if(sideID == 4)
      neID << 3, 0, 7;
    else
    {
      cerr << "sideID is invalid in boundary_condition_base.cpp" << endl;
      exit(1);
    }
  }
  else if(eType == "tetra4")
  {
    eTypeSur = "tria3";
    ne = 3;
    ipmax = 1;
    neID.resize(ne);
    if(sideID == 1)
      neID << 0, 1, 2;
    else if(sideID == 2)
      neID << 0, 1, 3;
    else if(sideID == 3)
      neID << 1, 2, 3;
    else if(sideID == 4)
      neID << 2, 0, 3;
    else
    {
      cerr << "sideID is invalid in boundary_condition_base.cpp" << endl;
      exit(1);
    }
  }
  else if(eType == "tetra10")
  {
    eTypeSur = "tria6";
    ne = 6;
    ipmax = 3;
    neID.resize(ne);
    if(sideID == 1)
      neID << 0, 1, 2, 4, 5, 6;
    else if(sideID == 2)
      neID << 0, 1, 3, 4, 8, 7;
    else if(sideID == 3)
      neID << 1, 2, 3, 5, 9, 8;
    else if(sideID == 4)
      neID << 2, 0, 3, 6, 7, 9;
    else
    {
      cerr << "sideID is invalid in boundary_condition_base.cpp" << endl;
      exit(1);
    }
  }
  else if(eType == "hexa8")
  {
    eTypeSur = "quad4";
    ne = 4;
    ipmax = 4;
    neID.resize(ne);
    if(sideID == 1)
      neID << 0, 1, 2, 3;
    else if(sideID == 2)
      neID << 3, 0, 4, 7;
    else if(sideID == 3)
      neID << 0, 1, 5, 4;
    else if(sideID == 4)
      neID << 1, 2, 6, 5;
    else if(sideID == 5)
      neID << 2, 3, 7, 6;
    else if(sideID == 6)
      neID << 4, 5, 6, 7;
    else
    {
      cerr << "sideID is invalid in boundary_condition_base.cpp" << endl;
      exit(1);
    }
  }
  else if(eType == "hexa20")
  {
    eTypeSur = "quad8";
    ne = 8;
    ipmax = 9;
    neID.resize(ne);
    if(sideID == 1)
      neID << 0, 1, 2, 3, 8, 9, 10, 11;
    else if(sideID == 2)
      neID << 3, 0, 4, 7, 11, 16, 15, 19;
    else if(sideID == 3)
      neID << 0, 1, 5, 4, 8, 17, 12, 16;
    else if(sideID == 4)
      neID << 1, 2, 6, 5, 9, 18, 13, 17;
    else if(sideID == 5)
      neID << 2, 3, 7, 6, 10, 19, 14, 18;
    else if(sideID == 6)
      neID << 4, 5, 6, 7, 12, 13, 14, 15;
    else
    {
      cerr << "sideID is invalid in boundary_condition_base.cpp" << endl;
      exit(1);
    }
  }
  else
  {
    std::cerr << "Distributed loads cannot be applied to this element type in "
                 "boundary_condition_base.cpp."
              << endl;
    exit(1);
  }

  // convert to surface parameter
  MatrixXd Xtmp = X;
  VectorXi idoftmp = idof;
  X.resize(ne, fem_.ndim);
  idof.resize(fem_.dofnp * ne);

  for(int i = 0; i < ne; i++)
  {
    for(int j = 0; j < fem_.dofnp; j++)
      idof.coeffRef(fem_.dofnp * i + j) = idoftmp(fem_.dofnp * neID(i) + j);
    for(int j = 0; j < fem_.ndim; j++)
      X.coeffRef(i, j) = Xtmp(neID(i), j);
  }

  return forward_as_tuple(eTypeSur, ne, ipmax);
}

VectorXd SurfaceLoad::transLocalToGlobalPressure(MatrixXd &X, VectorXd &val,
                                                 int iline)
{
  assert(fem_.ndim == fem_.dofnp);

  VectorXd val_new = VectorXd::Zero(fem_.dofnp);
  // tria3 or tria6
  if(fem_.ndim == 2 && X.rows() == 3 || fem_.ndim == 2 && X.rows() == 6)
  {
    Matrix2d delv;
    Vector2d del1, del2;
    Matrix2d rotate;
    rotate << 0.0, -1.0, 1.0, 0.0;

    double tmp = val[0];
    val.setZero();
    val[1] = tmp;

    if(iline == 1)
    {
      del1(0) = X(1, 0) - X(0, 0);
      del1(1) = X(1, 1) - X(0, 1);
    }
    else if(iline == 2)
    {
      del1(0) = X(2, 0) - X(1, 0);
      del1(1) = X(2, 1) - X(1, 1);
    }
    else if(iline == 3)
    {
      del1(0) = X(0, 0) - X(2, 0);
      del1(1) = X(0, 1) - X(2, 1);
    }
    else
    {
      cerr << "error in boundary_condition_base.cpp" << endl;
      exit(1);
    }

    del2 = rotate * del1;
    del1.normalized();
    del2.normalized();

    delv << del1(0), del2(0), del1(1), del2(1);

    // rotation of val
    val_new = delv * val;
  }
  // tetra4 or tetra10
  else if(fem_.ndim == 3 && X.rows() == 4 || fem_.ndim == 3 && X.rows() == 10)
  {
    Vector3d del1, del2, del3;
    if(iline == 1)
    {
      del1(0) = X(1, 0) - X(0, 0);
      del1(1) = X(1, 1) - X(0, 1);
      del1(2) = X(1, 2) - X(0, 2);
      del2(0) = X(2, 0) - X(0, 0);
      del2(1) = X(2, 1) - X(0, 1);
      del2(2) = X(2, 2) - X(0, 2);
    }
    else if(iline == 2)
    {
      del1(0) = X(3, 0) - X(0, 0);
      del1(1) = X(3, 1) - X(0, 1);
      del1(2) = X(3, 2) - X(0, 2);
      del2(0) = X(1, 0) - X(0, 0);
      del2(1) = X(1, 1) - X(0, 1);
      del2(2) = X(1, 2) - X(0, 2);
    }
    else if(iline == 3)
    {
      del1(0) = X(3, 0) - X(1, 0);
      del1(1) = X(3, 1) - X(1, 1);
      del1(2) = X(3, 2) - X(1, 2);
      del2(0) = X(2, 0) - X(1, 0);
      del2(1) = X(2, 1) - X(1, 1);
      del2(2) = X(2, 2) - X(1, 2);
    }
    else if(iline == 4)
    {
      del1(0) = X(2, 0) - X(0, 0);
      del1(1) = X(2, 1) - X(0, 1);
      del1(2) = X(2, 2) - X(0, 2);
      del2(0) = X(3, 0) - X(0, 0);
      del2(1) = X(3, 1) - X(0, 1);
      del2(2) = X(3, 2) - X(0, 2);
    }
    else
    {
      cerr << "error in boundary_condition_base.cpp" << endl;
      exit(1);
    }
    del3 = del1.cross(del2);
    del3.normalize();

    val_new = val[0] * del3;
  }
  else
  {
    cerr << "Error in boundary_condition_base.cpp" << endl;
    exit(1);
  }

  // // 2D
  // if(fem_.ndim == 2)
  // {
  //   Matrix2d delv;
  //   Vector2d del1, del2;

  //   if(iline == 1)
  //   {
  //     del1(0) = X(1, 0) - X(0, 0);
  //     del1(1) = X(1, 1) - X(0, 1);
  //     del2(0) = X(3, 0) - X(0, 0);
  //     del2(1) = X(3, 1) - X(0, 1);
  //   }
  //   else if(iline == 2)
  //   {
  //     del1(0) = X(2, 0) - X(1, 0);
  //     del1(1) = X(2, 1) - X(1, 1);
  //     del2(0) = X(0, 0) - X(1, 0);
  //     del2(1) = X(0, 1) - X(1, 1);
  //   }
  //   else if(iline == 3)
  //   {
  //     del1(0) = X(3, 0) - X(2, 0);
  //     del1(1) = X(3, 1) - X(2, 1);
  //     del2(0) = X(1, 0) - X(2, 0);
  //     del2(1) = X(1, 1) - X(2, 1);
  //   }
  //   else if(iline == 4)
  //   {
  //     del1(0) = X(0, 0) - X(3, 0);
  //     del1(1) = X(0, 1) - X(3, 1);
  //     del2(0) = X(2, 0) - X(3, 0);
  //     del2(1) = X(2, 1) - X(3, 1);
  //   }

  //   del1.normalized();
  //   del2.normalized();

  //   delv << del1(0), del2(0), del1(1), del2(1);

  //   // rotation of val
  //   val_new = delv * valold;
  // }
  // // 3D solid
  // else if(fem_.ndim == 3 && fem_.dofnp != 6)
  // {
  //   Matrix3d delv;
  //   Vector3d del1, del2, del3;

  //   // Zmin
  //   if(iline == 1)
  //   {
  //     del1(0) = X(1, 0) - X(0, 0);
  //     del1(1) = X(1, 1) - X(0, 1);
  //     del1(2) = X(1, 2) - X(0, 2);
  //     del2(0) = X(3, 0) - X(0, 0);
  //     del2(1) = X(3, 1) - X(0, 1);
  //     del2(2) = X(3, 2) - X(0, 2);
  //     del3(0) = X(4, 0) - X(0, 0);
  //     del3(1) = X(4, 1) - X(0, 1);
  //     del3(2) = X(4, 2) - X(0, 2);
  //   }
  //   // Xmin
  //   else if(iline == 2)
  //   {
  //     del1(0) =
  //         node[element[nel].nodeID[0]].x[0] -
  //         node[element[nel].nodeID[4]].x[0];
  //     del1(1) =
  //         node[element[nel].nodeID[0]].x[1] -
  //         node[element[nel].nodeID[4]].x[1];
  //     del1(2) =
  //         node[element[nel].nodeID[0]].x[2] -
  //         node[element[nel].nodeID[4]].x[2];
  //     del2(0) =
  //         node[element[nel].nodeID[7]].x[0] -
  //         node[element[nel].nodeID[4]].x[0];
  //     del2(1) =
  //         node[element[nel].nodeID[7]].x[1] -
  //         node[element[nel].nodeID[4]].x[1];
  //     del2(2) =
  //         node[element[nel].nodeID[7]].x[2] -
  //         node[element[nel].nodeID[4]].x[2];
  //     del3(0) =
  //         node[element[nel].nodeID[5]].x[0] -
  //         node[element[nel].nodeID[4]].x[0];
  //     del3(1) =
  //         node[element[nel].nodeID[5]].x[1] -
  //         node[element[nel].nodeID[4]].x[1];
  //     del3(2) =
  //         node[element[nel].nodeID[5]].x[2] -
  //         node[element[nel].nodeID[4]].x[2];
  //   }
  //   // Ymin
  //   else if(iline == 3)
  //   {
  //     del1(0) =
  //         node[element[nel].nodeID[5]].x[0] -
  //         node[element[nel].nodeID[4]].x[0];
  //     del1(1) =
  //         node[element[nel].nodeID[5]].x[1] -
  //         node[element[nel].nodeID[4]].x[1];
  //     del1(2) =
  //         node[element[nel].nodeID[5]].x[2] -
  //         node[element[nel].nodeID[4]].x[2];
  //     del2(0) =
  //         node[element[nel].nodeID[0]].x[0] -
  //         node[element[nel].nodeID[4]].x[0];
  //     del2(1) =
  //         node[element[nel].nodeID[0]].x[1] -
  //         node[element[nel].nodeID[4]].x[1];
  //     del2(2) =
  //         node[element[nel].nodeID[0]].x[2] -
  //         node[element[nel].nodeID[4]].x[2];
  //     del3(0) =
  //         node[element[nel].nodeID[7]].x[0] -
  //         node[element[nel].nodeID[4]].x[0];
  //     del3(1) =
  //         node[element[nel].nodeID[7]].x[1] -
  //         node[element[nel].nodeID[4]].x[1];
  //     del3(2) =
  //         node[element[nel].nodeID[7]].x[2] -
  //         node[element[nel].nodeID[4]].x[2];
  //   }
  //   // Xmax
  //   else if(iline == 4)
  //   {
  //     del1(0) =
  //         node[element[nel].nodeID[5]].x[0] -
  //         node[element[nel].nodeID[1]].x[0];
  //     del1(1) =
  //         node[element[nel].nodeID[5]].x[1] -
  //         node[element[nel].nodeID[1]].x[1];
  //     del1(2) =
  //         node[element[nel].nodeID[5]].x[2] -
  //         node[element[nel].nodeID[1]].x[2];
  //     del2(0) =
  //         node[element[nel].nodeID[2]].x[0] -
  //         node[element[nel].nodeID[1]].x[0];
  //     del2(1) =
  //         node[element[nel].nodeID[2]].x[1] -
  //         node[element[nel].nodeID[1]].x[1];
  //     del2(2) =
  //         node[element[nel].nodeID[2]].x[2] -
  //         node[element[nel].nodeID[1]].x[2];
  //     del3(0) =
  //         node[element[nel].nodeID[0]].x[0] -
  //         node[element[nel].nodeID[1]].x[0];
  //     del3(1) =
  //         node[element[nel].nodeID[0]].x[1] -
  //         node[element[nel].nodeID[1]].x[1];
  //     del3(2) =
  //         node[element[nel].nodeID[0]].x[2] -
  //         node[element[nel].nodeID[1]].x[2];
  //   }
  //   // Ymax
  //   else if(iline == 5)
  //   {
  //     del1(0) =
  //         node[element[nel].nodeID[3]].x[0] -
  //         node[element[nel].nodeID[2]].x[0];
  //     del1(1) =
  //         node[element[nel].nodeID[3]].x[1] -
  //         node[element[nel].nodeID[2]].x[1];
  //     del1(2) =
  //         node[element[nel].nodeID[3]].x[2] -
  //         node[element[nel].nodeID[2]].x[2];
  //     del2(0) =
  //         node[element[nel].nodeID[1]].x[0] -
  //         node[element[nel].nodeID[2]].x[0];
  //     del2(1) =
  //         node[element[nel].nodeID[1]].x[1] -
  //         node[element[nel].nodeID[2]].x[1];
  //     del2(2) =
  //         node[element[nel].nodeID[1]].x[2] -
  //         node[element[nel].nodeID[2]].x[2];
  //     del3(0) =
  //         node[element[nel].nodeID[6]].x[0] -
  //         node[element[nel].nodeID[2]].x[0];
  //     del3(1) =
  //         node[element[nel].nodeID[6]].x[1] -
  //         node[element[nel].nodeID[2]].x[1];
  //     del3(2) =
  //         node[element[nel].nodeID[6]].x[2] -
  //         node[element[nel].nodeID[2]].x[2];
  //   }
  //   // Zmax
  //   else if(iline == 6)
  //   {
  //     del1(0) =
  //         node[element[nel].nodeID[4]].x[0] -
  //         node[element[nel].nodeID[5]].x[0];
  //     del1(1) =
  //         node[element[nel].nodeID[4]].x[1] -
  //         node[element[nel].nodeID[5]].x[1];
  //     del1(2) =
  //         node[element[nel].nodeID[4]].x[2] -
  //         node[element[nel].nodeID[5]].x[2];
  //     del2(0) =
  //         node[element[nel].nodeID[6]].x[0] -
  //         node[element[nel].nodeID[5]].x[0];
  //     del2(1) =
  //         node[element[nel].nodeID[6]].x[1] -
  //         node[element[nel].nodeID[5]].x[1];
  //     del2(2) =
  //         node[element[nel].nodeID[6]].x[2] -
  //         node[element[nel].nodeID[5]].x[2];
  //     del3(0) =
  //         node[element[nel].nodeID[1]].x[0] -
  //         node[element[nel].nodeID[5]].x[0];
  //     del3(1) =
  //         node[element[nel].nodeID[1]].x[1] -
  //         node[element[nel].nodeID[5]].x[1];
  //     del3(2) =
  //         node[element[nel].nodeID[1]].x[2] -
  //         node[element[nel].nodeID[5]].x[2];
  //   }

  //   del1.normalized();
  //   del2.normalized();
  //   del3.normalized();

  //   delv << del1(0), del2(0), del3(0), del1(1), del2(1), del3(1), del1(2),
  //       del2(2), del3(2);

  //   // rotation of val
  //   val_new = delv * val;
  // }
  // // 3D shell
  // else if(fem.dofPerNode == 6)
  // {
  //   Matrix3d delv;
  //   Vector3d del1, del2, del3;

  //   // Ymin
  //   if(iline == 1)
  //   {
  //     del1(0) =
  //         node[element[nel].nodeID[1]].x[0] -
  //         node[element[nel].nodeID[0]].x[0];
  //     del1(1) =
  //         node[element[nel].nodeID[1]].x[1] -
  //         node[element[nel].nodeID[0]].x[1];
  //     del1(2) =
  //         node[element[nel].nodeID[1]].x[2] -
  //         node[element[nel].nodeID[0]].x[2];
  //     del2(0) =
  //         node[element[nel].nodeID[3]].x[0] -
  //         node[element[nel].nodeID[0]].x[0];
  //     del2(1) =
  //         node[element[nel].nodeID[3]].x[1] -
  //         node[element[nel].nodeID[0]].x[1];
  //     del2(2) =
  //         node[element[nel].nodeID[3]].x[2] -
  //         node[element[nel].nodeID[0]].x[2];
  //   }
  //   // Xmax
  //   else if(iline == 2)
  //   {
  //     del1(0) =
  //         node[element[nel].nodeID[2]].x[0] -
  //         node[element[nel].nodeID[1]].x[0];
  //     del1(1) =
  //         node[element[nel].nodeID[2]].x[1] -
  //         node[element[nel].nodeID[1]].x[1];
  //     del1(2) =
  //         node[element[nel].nodeID[2]].x[2] -
  //         node[element[nel].nodeID[1]].x[2];
  //     del2(0) =
  //         node[element[nel].nodeID[0]].x[0] -
  //         node[element[nel].nodeID[1]].x[0];
  //     del2(1) =
  //         node[element[nel].nodeID[0]].x[1] -
  //         node[element[nel].nodeID[1]].x[1];
  //     del2(2) =
  //         node[element[nel].nodeID[0]].x[2] -
  //         node[element[nel].nodeID[1]].x[2];
  //   }
  //   // Xmin
  //   else if(iline == 3)
  //   {
  //     del1(0) =
  //         node[element[nel].nodeID[3]].x[0] -
  //         node[element[nel].nodeID[2]].x[0];
  //     del1(1) =
  //         node[element[nel].nodeID[3]].x[1] -
  //         node[element[nel].nodeID[2]].x[1];
  //     del1(2) =
  //         node[element[nel].nodeID[3]].x[2] -
  //         node[element[nel].nodeID[2]].x[2];
  //     del2(0) =
  //         node[element[nel].nodeID[1]].x[0] -
  //         node[element[nel].nodeID[2]].x[0];
  //     del2(1) =
  //         node[element[nel].nodeID[1]].x[1] -
  //         node[element[nel].nodeID[2]].x[1];
  //     del2(2) =
  //         node[element[nel].nodeID[1]].x[2] -
  //         node[element[nel].nodeID[2]].x[2];
  //   }
  //   // Xmin
  //   else if(iline == 4)
  //   {
  //     del1(0) =
  //         node[element[nel].nodeID[0]].x[0] -
  //         node[element[nel].nodeID[3]].x[0];
  //     del1(1) =
  //         node[element[nel].nodeID[0]].x[1] -
  //         node[element[nel].nodeID[3]].x[1];
  //     del1(2) =
  //         node[element[nel].nodeID[0]].x[2] -
  //         node[element[nel].nodeID[3]].x[2];
  //     del2(0) =
  //         node[element[nel].nodeID[2]].x[0] -
  //         node[element[nel].nodeID[3]].x[0];
  //     del2(1) =
  //         node[element[nel].nodeID[2]].x[1] -
  //         node[element[nel].nodeID[3]].x[1];
  //     del2(2) =
  //         node[element[nel].nodeID[2]].x[2] -
  //         node[element[nel].nodeID[3]].x[2];
  //   }

  //   del1.normalized();
  //   del2.normalized();
  //   del3 = (del1.cross(del2)).normalized();

  //   delv << del1(0), del2(0), del3(0), del1(1), del2(1), del3(1), del1(2),
  //       del2(2), del3(2);

  //   Vector3d val_3size;
  //   for(int i = 0; i < 3; i++)
  //     val_3size[i] = val[i];

  //   // rotation of val
  //   val_3size = delv * val_3size;

  //   for(int i = 0; i < 3; i++)
  //     val_new[i] = val_3size[i];
  // }
  // else
  // {
  //   std::cerr << "Error inboundary_condition_icarat.hpp" << std::endl;
  //   exit(1);
  // }

  //  cout << val_new << endl;

  return val_new;
}

VectorXd SurfaceLoad::transLocalToGlobalMoment(MatrixXd &X, VectorXd &val,
                                               InputLoad &inp_f)
{
  assert(fem_.dofnp == fem_.ndim);
  VectorXd val_new(fem_.dofnp);
  // 2D
  if(fem_.ndim == 2 && X.rows() == 2 || fem_.ndim == 2 && X.rows() == 3)
  {
    Vector2d origin, center;
    origin << inp_f.axis_origin[0], inp_f.axis_origin[1];
    center << (X(0, 0) + X(1, 0)) / 2.0, (X(0, 1) + X(1, 1)) / 2.0;
    Vector2d oricen = center - origin;
    Matrix2d rotate;
    rotate << 0.0, -1.0, 1.0, 0.0;
    val_new = val[0] * (rotate * oricen).normalized();
  }
  else if(fem_.ndim == 2 && X.rows() == 3 || fem_.ndim == 2 && X.rows() == 6)
  {
  }
  else
  {
    cerr << "error in boundary_condition_base.cpp" << endl;
    exit(1);
  }

  return val_new;
}
} // namespace icarat