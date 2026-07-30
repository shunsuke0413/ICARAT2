///  @file  base.hpp
///  @author  Daiki Watanabe
///  @date  May 31, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <memory>
#include <vector>
namespace icarat
{
///   problem setting related to FEM
struct FEM
{
  int neq;   ///< total dof
  int numeq; ///< number of equation(= neq - dirichletDOFs)
  int ndim;  ///< dimension
  int nelx;  ///< number of element
  int voigt; ///< the size of Voigt representation in a matrix
  int dofnp; ///< DOF per node(displacement=ndim,temperature=1,fluid=ndim+1)
  int numnp; ///< number of node
#ifdef ICARAT_MPI
  int procno;                 ///< number of MPI processes
  int procid;                 ///< ID of this process
  static const int block = 3; ///< Memory block size used in solver
  /// set in B.C. class
  int begdof; ///< beginning of dof
  int gnumeq; ///< global number of equation(= neq - dirichletDOFs)
  int dumdof; ///< Dummy dof resulting from blocking
#endif
};

/// basic element variable list class in FEM
struct Element
{
  /// set in element constructor
  int numdof; ///< dof in the element
  int ne;     ///< number of nodes in this element
  int ipmax;  ///< number of gauss point in this element

  /// set in mesh class
  std::string eType;       ///< element type
  int ID;                  ///< element ID
  std::vector<int> nodeID; ///< connectivity

  /// used in FEA
  double volume; ///< volume in this element
};

///  structure on Dirichlet conditions.
struct Dirichlet
{
  /// flag whether dirichlet is effective(=1) or not(=0).
  int flag[6];
  /// value of dirichlet(used only macro)
  double val[6];
};

///  structure on Neumann conditions.
struct Neumann
{
  int flag[3];   ///< dofflags of neumann node
  double val[3]; ///< value of load
};

/// multi point constraint structure
struct MPC
{
  /// DOF of control point(size: [dofnp][numCnt])
  /// numCnt is number of control point involved
  std::vector<std::vector<int>> cdof;

  /// value of control point (size is same of cdof)
  // std::vector<std::vector<double>> cval;
};

/// node variable list
struct Node
{
  Node(int totalTime) : valeach(6, totalTime), reaeach(3, totalTime)
  {
    valeach.setZero();
    reaeach.setZero();
  }
  // set in mesh class
  int ID;      ///< node ID
  int dof[6];  ///< (local) degree of freedom.
               ///<  When a dof is negative value, thedof is null.
  double x[3]; ///< coodinates of node

  // set in bc class
  std::shared_ptr<Dirichlet> dirich; ///< dirichlet BC structure
  std::shared_ptr<Neumann> neum;     ///< neumann  BC structure

  double val[6]; ///< solved value
  Eigen::MatrixXd valeach;

  double reaction[3]; ///< reaction force
  Eigen::MatrixXd reaeach;

  /// displacement increment in previous step (for increment-type material)
  double du[3];
  double du_old[3];

#ifdef ICARAT_MPI
  int gdof[6];   ///< global degree of freedom
  bool isMaster; ///< master/slave node (private node = true)
  int mprocid;   ///< master node's process ID
  int type;      ///< number of shared processes in this node.
  std::vector<std::pair<int, int>> sprocID; ///< share process & ID
#endif
};

struct Force
{
  Force(int numeq, int totalTime)
      : fext_org(numeq), fext(numeq), fint(numeq), fintarray(numeq, totalTime),
        fres(numeq), fexternal(numeq, totalTime), effec(numeq)
  {
    fext_org.setZero();
    fext.setZero();
    fint.setZero();
    fintarray.setZero();
    fres.setZero();
    fexternal.setZero();
    effec.setZero();
  }

  // set in B.C. class
  Eigen::VectorXd fext_org;  ///< orignal external force or heat
  Eigen::VectorXd fext;      ///< external force at the increment step
  Eigen::VectorXd fint;      ///< internal force
  Eigen::MatrixXd fintarray; ///< internal force
  Eigen::VectorXd fres;      ///< residual force (= fext - fint)
  Eigen::MatrixXd fexternal;
  Eigen::VectorXd effec;
};

struct Value
{
  Value(int numeq) : val(numeq) { val.setZero(); }
  Eigen::VectorXd val; ///< solved value
};

/// solution values of linear eigen problem
///(displacement & eigenvectors & eigenvalues)
struct ValueEigen : Value
{
  ValueEigen(int numeq, int numEigen) : Value(numeq), evalues(numEigen)
  {
    for(int i = 0; i < numEigen; i++)
    {
      Eigen::VectorXd evector = Eigen::VectorXd::Zero(numeq);
      evectors.push_back(evector);
    }
  }

  std::vector<double> evalues;           /// eigenvalues
  std::vector<Eigen::VectorXd> evectors; ///< eigenvectors
};

/// values of dynamic problem
struct ValueDyna : Value
{
  ValueDyna(int numeq)
      : Value(numeq), accel(numeq), velocity(numeq), du_incre(numeq)
  {
    accel.setZero();
    velocity.setZero();
    du_incre.setZero();
  }

  Eigen::VectorXd accel;    ///< acceleration
  Eigen::VectorXd velocity; ///< velocity
  Eigen::VectorXd du_incre; ///< time increment displacement
};

/// values of dynamic problem
struct ValueNR : Value
{
  ValueNR(int numeq) : Value(numeq), du(numeq) { du.setZero(); }

  Eigen::VectorXd du; ///< incremental displacement
};

/// solution of numerical material test of homogenization analysis
struct ValueHomo : Value
{
  ValueHomo(int numeq, int voigt) : Value(numeq)
  {
    for(int i = 0; i < voigt; i++)
    {
      Eigen::VectorXd value = Eigen::VectorXd::Zero(numeq);
      values.push_back(value);
    }
  }

  std::vector<Eigen::VectorXd> values; ///< response displacement
};

} // namespace icarat