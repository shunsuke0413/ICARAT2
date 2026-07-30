#include "make_element.hpp"
namespace icarat
{
namespace test
{
using namespace std;

void make1quad4(FEM &fem, vector<Element> &element, vector<Node> &node)
{
  // fem
  fem.ndim = 2;
  fem.nelx = 1;
  fem.voigt = 3;
  fem.dofnp = 2;
  fem.numnp = 4;
  fem.neq = fem.numnp * fem.dofnp;
  fem.numeq = fem.neq;

  // element
  Element actele;
  actele.ID = 0;
  actele.eType = "quad4";
  actele.ne = 4;
  actele.ipmax = 4;
  actele.numdof = fem.dofnp * actele.ne;
  actele.volume = 1.0;
  actele.nodeID.resize(4);
  actele.nodeID[0] = 0;
  actele.nodeID[1] = 1;
  actele.nodeID[2] = 3;
  actele.nodeID[3] = 2;
  element.push_back(actele);

  // node
  for(int i = 0; i < 2; i++)
  {
    for(int j = 0; j < 2; j++)
    {
      Node actnode;
      actnode.ID = 2 * i + j;
      actnode.dof[0] = 2 * actnode.ID;
      actnode.dof[1] = 2 * actnode.ID + 1;
      actnode.x[0] = 1.0 * j;
      actnode.x[1] = 1.0 * i;
      node.push_back(actnode);
    }
  }
}

void make1quad8(FEM &fem, vector<Element> &element, vector<Node> &node)
{
  // fem
  fem.ndim = 2;
  fem.nelx = 1;
  fem.voigt = 3;
  fem.dofnp = 2;
  fem.numnp = 8;
  fem.neq = fem.numnp * fem.dofnp;
  fem.numeq = fem.neq;

  // element
  Element actele;
  actele.ID = 0;
  actele.eType = "quad8";
  actele.ne = 8;
  actele.ipmax = 9;
  actele.numdof = fem.dofnp * actele.ne;
  actele.volume = 1.0;
  actele.nodeID.resize(8);
  actele.nodeID[0] = 0;
  actele.nodeID[1] = 2;
  actele.nodeID[2] = 7;
  actele.nodeID[3] = 5;
  actele.nodeID[4] = 1;
  actele.nodeID[5] = 4;
  actele.nodeID[6] = 6;
  actele.nodeID[7] = 3;
  element.push_back(actele);

  // node0
  Node actnode;
  actnode.ID = 0;
  actnode.dof[0] = 0;
  actnode.dof[1] = 1;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  node.push_back(actnode);
  // node1
  actnode.ID = 1;
  actnode.dof[0] = 2;
  actnode.dof[1] = 3;
  actnode.x[0] = 0.5;
  actnode.x[1] = 0.0;
  node.push_back(actnode);
  // node2
  actnode.ID = 2;
  actnode.dof[0] = 4;
  actnode.dof[1] = 5;
  actnode.x[0] = 1.0;
  actnode.x[1] = 0.0;
  node.push_back(actnode);
  // node3
  actnode.ID = 3;
  actnode.dof[0] = 6;
  actnode.dof[1] = 7;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.5;
  node.push_back(actnode);
  // node4
  actnode.ID = 4;
  actnode.dof[0] = 8;
  actnode.dof[1] = 9;
  actnode.x[0] = 1.0;
  actnode.x[1] = 0.5;
  node.push_back(actnode);
  // node5
  actnode.ID = 5;
  actnode.dof[0] = 10;
  actnode.dof[1] = 11;
  actnode.x[0] = 0.0;
  actnode.x[1] = 1.0;
  node.push_back(actnode);
  // node6
  actnode.ID = 6;
  actnode.dof[0] = 12;
  actnode.dof[1] = 13;
  actnode.x[0] = 0.5;
  actnode.x[1] = 1.0;
  node.push_back(actnode);
  // node7
  actnode.ID = 7;
  actnode.dof[0] = 14;
  actnode.dof[1] = 15;
  actnode.x[0] = 1.0;
  actnode.x[1] = 1.0;
  node.push_back(actnode);
}

void make1tria3(FEM &fem, vector<Element> &element, vector<Node> &node)
{
  // fem
  fem.ndim = 2;
  fem.nelx = 1;
  fem.voigt = 3;
  fem.dofnp = 2;
  fem.numnp = 3;
  fem.neq = fem.numnp * fem.dofnp;
  fem.numeq = fem.neq;

  // element
  Element actele;
  actele.ID = 0;
  actele.eType = "tria3";
  actele.ne = 3;
  actele.ipmax = 1;
  actele.numdof = fem.dofnp * actele.ne;
  actele.volume = 0.5;
  actele.nodeID.resize(3);
  actele.nodeID[0] = 0;
  actele.nodeID[1] = 1;
  actele.nodeID[2] = 2;

  element.push_back(actele);

  // node 0, 1, 2
  for(int i = 0; i < 2; i++)
  {
    for(int j = 0; j < 2; j++)
    {
      if(i == 1 && j == 1)
        continue;

      Node actnode;
      actnode.ID = 2 * i + j;
      actnode.dof[0] = 2 * actnode.ID;
      actnode.dof[1] = 2 * actnode.ID + 1;
      actnode.x[0] = 1.0 * j;
      actnode.x[1] = 1.0 * i;
      node.push_back(actnode);
    }
  }
}

void make1tria6(FEM &fem, vector<Element> &element, vector<Node> &node)
{
  // fem
  fem.ndim = 2;
  fem.nelx = 1;
  fem.voigt = 3;
  fem.dofnp = 2;
  fem.numnp = 6;
  fem.neq = fem.numnp * fem.dofnp;
  fem.numeq = fem.neq;

  // element
  Element actele;
  actele.ID = 0;
  actele.eType = "tria6";
  actele.ne = 6;
  actele.ipmax = 3;
  actele.numdof = fem.dofnp * actele.ne;
  actele.volume = 0.5;
  actele.nodeID.resize(6);
  actele.nodeID[0] = 0;
  actele.nodeID[1] = 2;
  actele.nodeID[2] = 5;
  actele.nodeID[3] = 1;
  actele.nodeID[4] = 4;
  actele.nodeID[5] = 3;

  element.push_back(actele);

  // node0
  Node actnode;
  actnode.ID = 0;
  actnode.dof[0] = 0;
  actnode.dof[1] = 1;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  node.push_back(actnode);

  // node1
  actnode.ID = 1;
  actnode.dof[0] = 2;
  actnode.dof[1] = 3;
  actnode.x[0] = 0.5;
  actnode.x[1] = 0.0;
  node.push_back(actnode);

  // node2
  actnode.ID = 2;
  actnode.dof[0] = 4;
  actnode.dof[1] = 5;
  actnode.x[0] = 1.0;
  actnode.x[1] = 0.0;
  node.push_back(actnode);

  // node3
  actnode.ID = 3;
  actnode.dof[0] = 6;
  actnode.dof[1] = 7;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.5;
  node.push_back(actnode);

  // node4
  actnode.ID = 4;
  actnode.dof[0] = 8;
  actnode.dof[1] = 9;
  actnode.x[0] = 0.5;
  actnode.x[1] = 0.5;
  node.push_back(actnode);

  // node5
  actnode.ID = 5;
  actnode.dof[0] = 10;
  actnode.dof[1] = 11;
  actnode.x[0] = 0.0;
  actnode.x[1] = 1.0;
  node.push_back(actnode);
}

void make1cubic8(FEM &fem, vector<Element> &element, vector<Node> &node)
{
  // fem
  fem.ndim = 3;
  fem.nelx = 1;
  fem.voigt = 6;
  fem.dofnp = 3;
  fem.numnp = 8;
  fem.neq = fem.numnp * fem.dofnp;
  fem.numeq = fem.neq;

  // element
  Element actele;
  actele.ID = 0;
  actele.eType = "hexa8";
  actele.ne = 8;
  actele.ipmax = 8;
  actele.numdof = fem.dofnp * actele.ne;
  actele.volume = 1.0;
  actele.nodeID.resize(8);
  actele.nodeID[0] = 0;
  actele.nodeID[1] = 1;
  actele.nodeID[2] = 3;
  actele.nodeID[3] = 2;
  actele.nodeID[4] = 4;
  actele.nodeID[5] = 5;
  actele.nodeID[6] = 7;
  actele.nodeID[7] = 6;

  element.push_back(actele);

  // node0
  Node actnode;
  actnode.ID = 0;
  actnode.dof[0] = 0;
  actnode.dof[1] = 1;
  actnode.dof[2] = 2;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node1
  actnode.ID = 1;
  actnode.dof[0] = 3;
  actnode.dof[1] = 4;
  actnode.dof[2] = 5;
  actnode.x[0] = 1.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node2
  actnode.ID = 2;
  actnode.dof[0] = 6;
  actnode.dof[1] = 7;
  actnode.dof[2] = 8;
  actnode.x[0] = 0.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node3
  actnode.ID = 3;
  actnode.dof[0] = 9;
  actnode.dof[1] = 10;
  actnode.dof[2] = 11;
  actnode.x[0] = 1.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node4
  actnode.ID = 4;
  actnode.dof[0] = 12;
  actnode.dof[1] = 13;
  actnode.dof[2] = 14;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node5
  actnode.ID = 5;
  actnode.dof[0] = 15;
  actnode.dof[1] = 16;
  actnode.dof[2] = 17;
  actnode.x[0] = 1.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node6
  actnode.ID = 6;
  actnode.dof[0] = 18;
  actnode.dof[1] = 19;
  actnode.dof[2] = 20;
  actnode.x[0] = 0.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node7
  actnode.ID = 7;
  actnode.dof[0] = 21;
  actnode.dof[1] = 22;
  actnode.dof[2] = 23;
  actnode.x[0] = 1.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);
}

void make1cubic20(FEM &fem, vector<Element> &element, vector<Node> &node)
{
  // fem
  fem.ndim = 3;
  fem.nelx = 1;
  fem.voigt = 6;
  fem.dofnp = 3;
  fem.numnp = 20;
  fem.neq = fem.numnp * fem.dofnp;
  fem.numeq = fem.neq;

  // element
  Element actele;
  actele.ID = 0;
  actele.eType = "hexa20";
  actele.ne = 20;
  actele.ipmax = 27;
  actele.numdof = fem.dofnp * actele.ne;
  actele.volume = 1.0;
  actele.nodeID.resize(20);
  actele.nodeID[0] = 0;
  actele.nodeID[1] = 1;
  actele.nodeID[2] = 2;
  actele.nodeID[3] = 4;
  actele.nodeID[4] = 4;
  actele.nodeID[5] = 5;
  actele.nodeID[6] = 6;
  actele.nodeID[7] = 7;
  actele.nodeID[8] = 8;
  actele.nodeID[9] = 9;
  actele.nodeID[10] = 10;
  actele.nodeID[11] = 11;
  actele.nodeID[12] = 12;
  actele.nodeID[13] = 13;
  actele.nodeID[14] = 14;
  actele.nodeID[15] = 15;
  actele.nodeID[16] = 16;
  actele.nodeID[17] = 17;
  actele.nodeID[18] = 18;
  actele.nodeID[19] = 19;

  element.push_back(actele);

  // node0
  Node actnode;
  actnode.ID = 0;
  actnode.dof[0] = 0;
  actnode.dof[1] = 1;
  actnode.dof[2] = 2;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node1
  actnode.ID = 1;
  actnode.dof[0] = 3;
  actnode.dof[1] = 4;
  actnode.dof[2] = 5;
  actnode.x[0] = 1.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node2
  actnode.ID = 2;
  actnode.dof[0] = 6;
  actnode.dof[1] = 7;
  actnode.dof[2] = 8;
  actnode.x[0] = 1.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node3
  actnode.ID = 3;
  actnode.dof[0] = 9;
  actnode.dof[1] = 10;
  actnode.dof[2] = 11;
  actnode.x[0] = 0.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node4
  actnode.ID = 4;
  actnode.dof[0] = 12;
  actnode.dof[1] = 13;
  actnode.dof[2] = 14;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node5
  actnode.ID = 5;
  actnode.dof[0] = 15;
  actnode.dof[1] = 16;
  actnode.dof[2] = 17;
  actnode.x[0] = 1.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node6
  actnode.ID = 6;
  actnode.dof[0] = 18;
  actnode.dof[1] = 19;
  actnode.dof[2] = 20;
  actnode.x[0] = 1.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node7
  actnode.ID = 7;
  actnode.dof[0] = 21;
  actnode.dof[1] = 22;
  actnode.dof[2] = 23;
  actnode.x[0] = 0.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node8
  actnode.ID = 8;
  actnode.dof[0] = 24;
  actnode.dof[1] = 25;
  actnode.dof[2] = 26;
  actnode.x[0] = 0.5;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node9
  actnode.ID = 9;
  actnode.dof[0] = 27;
  actnode.dof[1] = 28;
  actnode.dof[2] = 29;
  actnode.x[0] = 1.0;
  actnode.x[1] = 0.5;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node10
  actnode.ID = 10;
  actnode.dof[0] = 30;
  actnode.dof[1] = 31;
  actnode.dof[2] = 32;
  actnode.x[0] = 0.5;
  actnode.x[1] = 1.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node11
  actnode.ID = 11;
  actnode.dof[0] = 33;
  actnode.dof[1] = 34;
  actnode.dof[2] = 35;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.5;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node12
  actnode.ID = 12;
  actnode.dof[0] = 36;
  actnode.dof[1] = 37;
  actnode.dof[2] = 38;
  actnode.x[0] = 0.5;
  actnode.x[1] = 0.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node13
  actnode.ID = 13;
  actnode.dof[0] = 39;
  actnode.dof[1] = 40;
  actnode.dof[2] = 41;
  actnode.x[0] = 1.0;
  actnode.x[1] = 0.5;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node14
  actnode.ID = 14;
  actnode.dof[0] = 42;
  actnode.dof[1] = 43;
  actnode.dof[2] = 44;
  actnode.x[0] = 0.5;
  actnode.x[1] = 1.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node15
  actnode.ID = 15;
  actnode.dof[0] = 45;
  actnode.dof[1] = 46;
  actnode.dof[2] = 47;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.5;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node16
  actnode.ID = 16;
  actnode.dof[0] = 48;
  actnode.dof[1] = 49;
  actnode.dof[2] = 50;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.5;
  node.push_back(actnode);

  // node17
  actnode.ID = 17;
  actnode.dof[0] = 51;
  actnode.dof[1] = 52;
  actnode.dof[2] = 53;
  actnode.x[0] = 0.5;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.5;
  node.push_back(actnode);

  // node18
  actnode.ID = 18;
  actnode.dof[0] = 54;
  actnode.dof[1] = 55;
  actnode.dof[2] = 56;
  actnode.x[0] = 1.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 0.5;
  node.push_back(actnode);

  // node19
  actnode.ID = 19;
  actnode.dof[0] = 57;
  actnode.dof[1] = 58;
  actnode.dof[2] = 59;
  actnode.x[0] = 0.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 0.5;
  node.push_back(actnode);
}

void make1tetra4(FEM &fem, vector<Element> &element, vector<Node> &node)
{
  // fem
  fem.ndim = 3;
  fem.nelx = 1;
  fem.voigt = 6;
  fem.dofnp = 3;
  fem.numnp = 4;
  fem.neq = fem.numnp * fem.dofnp;
  fem.numeq = fem.neq;

  // element
  Element actele;
  actele.ID = 0;
  actele.eType = "tetra4";
  actele.ne = 4;
  actele.ipmax = 1;
  actele.numdof = fem.dofnp * actele.ne;
  actele.volume = 1.0;
  actele.nodeID.resize(4);
  actele.nodeID[0] = 0;
  actele.nodeID[1] = 2;
  actele.nodeID[2] = 1;
  actele.nodeID[3] = 3;

  element.push_back(actele);

  // node0
  Node actnode;
  actnode.ID = 0;
  actnode.dof[0] = 0;
  actnode.dof[1] = 1;
  actnode.dof[2] = 2;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node1
  actnode.ID = 1;
  actnode.dof[0] = 3;
  actnode.dof[1] = 4;
  actnode.dof[2] = 5;
  actnode.x[0] = 0.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node2
  actnode.ID = 2;
  actnode.dof[0] = 6;
  actnode.dof[1] = 7;
  actnode.dof[2] = 8;
  actnode.x[0] = 1.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node3
  actnode.ID = 3;
  actnode.dof[0] = 9;
  actnode.dof[1] = 10;
  actnode.dof[2] = 11;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);
}

void make1tetra10(FEM &fem, vector<Element> &element, vector<Node> &node)
{
  // fem
  fem.ndim = 3;
  fem.nelx = 1;
  fem.voigt = 6;
  fem.dofnp = 3;
  fem.numnp = 10;
  fem.neq = fem.numnp * fem.dofnp;
  fem.numeq = fem.neq;

  // element
  Element actele;
  actele.eType = "tetra10";
  actele.ID = 0;
  actele.ne = 10;
  actele.ipmax = 4;
  actele.numdof = fem.dofnp * actele.ne;
  actele.volume = 1.0;
  actele.nodeID.resize(10);
  actele.nodeID[0] = 0;
  actele.nodeID[1] = 1;
  actele.nodeID[2] = 2;
  actele.nodeID[3] = 3;
  actele.nodeID[4] = 4;
  actele.nodeID[5] = 5;
  actele.nodeID[6] = 6;
  actele.nodeID[7] = 7;
  actele.nodeID[8] = 8;
  actele.nodeID[9] = 9;

  element.push_back(actele);

  // node0
  Node actnode;
  actnode.ID = 0;
  actnode.dof[0] = 0;
  actnode.dof[1] = 1;
  actnode.dof[2] = 2;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node1
  actnode.ID = 1;
  actnode.dof[0] = 3;
  actnode.dof[1] = 4;
  actnode.dof[2] = 5;
  actnode.x[0] = 1.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node2
  actnode.ID = 2;
  actnode.dof[0] = 6;
  actnode.dof[1] = 7;
  actnode.dof[2] = 8;
  actnode.x[0] = 0.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node3
  actnode.ID = 3;
  actnode.dof[0] = 9;
  actnode.dof[1] = 10;
  actnode.dof[2] = 11;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node4
  actnode.ID = 4;
  actnode.dof[0] = 12;
  actnode.dof[1] = 13;
  actnode.dof[2] = 14;
  actnode.x[0] = 0.5;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node5
  actnode.ID = 5;
  actnode.dof[0] = 15;
  actnode.dof[1] = 16;
  actnode.dof[2] = 17;
  actnode.x[0] = 0.5;
  actnode.x[1] = 0.5;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node6
  actnode.ID = 6;
  actnode.dof[0] = 18;
  actnode.dof[1] = 19;
  actnode.dof[2] = 20;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.5;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node7
  actnode.ID = 7;
  actnode.dof[0] = 21;
  actnode.dof[1] = 22;
  actnode.dof[2] = 23;
  actnode.x[0] = 0.5;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.5;
  node.push_back(actnode);

  // node8
  actnode.ID = 8;
  actnode.dof[0] = 24;
  actnode.dof[1] = 25;
  actnode.dof[2] = 26;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.5;
  actnode.x[2] = 0.5;
  node.push_back(actnode);

  // node9
  actnode.ID = 9;
  actnode.dof[0] = 27;
  actnode.dof[1] = 28;
  actnode.dof[2] = 29;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.5;
  node.push_back(actnode);
}

void make1prism6(FEM &fem, vector<Element> &element, vector<Node> &node)
{
  // fem
  fem.ndim = 3;
  fem.nelx = 1;
  fem.voigt = 6;
  fem.dofnp = 3;
  fem.numnp = 6;
  fem.neq = fem.numnp * fem.dofnp;
  fem.numeq = fem.neq;

  // element
  Element actele;
  actele.ID = 0;
  actele.eType = "prism6";
  actele.ne = 6;
  actele.ipmax = 6;
  actele.numdof = fem.dofnp * actele.ne;
  actele.volume = 1.0;
  actele.nodeID.resize(6);
  actele.nodeID[0] = 0;
  actele.nodeID[1] = 1;
  actele.nodeID[2] = 2;
  actele.nodeID[3] = 3;
  actele.nodeID[4] = 4;
  actele.nodeID[5] = 5;

  element.push_back(actele);

  // node0
  Node actnode;
  actnode.ID = 0;
  actnode.dof[0] = 0;
  actnode.dof[1] = 1;
  actnode.dof[2] = 2;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node1
  actnode.ID = 1;
  actnode.dof[0] = 3;
  actnode.dof[1] = 4;
  actnode.dof[2] = 5;
  actnode.x[0] = 1.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node2
  actnode.ID = 2;
  actnode.dof[0] = 6;
  actnode.dof[1] = 7;
  actnode.dof[2] = 8;
  actnode.x[0] = 0.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node3
  actnode.ID = 3;
  actnode.dof[0] = 9;
  actnode.dof[1] = 10;
  actnode.dof[2] = 11;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node4
  actnode.ID = 4;
  actnode.dof[0] = 12;
  actnode.dof[1] = 13;
  actnode.dof[2] = 14;
  actnode.x[0] = 1.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node5
  actnode.ID = 5;
  actnode.dof[0] = 15;
  actnode.dof[1] = 16;
  actnode.dof[2] = 17;
  actnode.x[0] = 0.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);
}

void make1prism15(FEM &fem, vector<Element> &element, vector<Node> &node)
{
  // fem
  fem.ndim = 3;
  fem.nelx = 1;
  fem.voigt = 6;
  fem.dofnp = 3;
  fem.numnp = 15;
  fem.neq = fem.numnp * fem.dofnp;
  fem.numeq = fem.neq;

  // element
  Element actele;
  actele.ID = 0;
  actele.eType = "prism15";
  actele.ne = 15;
  actele.ipmax = 21;
  actele.numdof = fem.dofnp * actele.ne;
  actele.volume = 1.0;
  actele.nodeID.resize(15);
  actele.nodeID[0] = 0;
  actele.nodeID[1] = 1;
  actele.nodeID[2] = 2;
  actele.nodeID[3] = 3;
  actele.nodeID[4] = 4;
  actele.nodeID[5] = 5;
  actele.nodeID[6] = 6;
  actele.nodeID[7] = 7;
  actele.nodeID[8] = 8;
  actele.nodeID[9] = 9;
  actele.nodeID[10] = 10;
  actele.nodeID[11] = 11;
  actele.nodeID[12] = 12;
  actele.nodeID[13] = 13;
  actele.nodeID[14] = 14;

  element.push_back(actele);

  // node0
  Node actnode;
  actnode.ID = 0;
  actnode.dof[0] = 0;
  actnode.dof[1] = 1;
  actnode.dof[2] = 2;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node1
  actnode.ID = 1;
  actnode.dof[0] = 3;
  actnode.dof[1] = 4;
  actnode.dof[2] = 5;
  actnode.x[0] = 1.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node2
  actnode.ID = 2;
  actnode.dof[0] = 6;
  actnode.dof[1] = 7;
  actnode.dof[2] = 8;
  actnode.x[0] = 0.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node3
  actnode.ID = 3;
  actnode.dof[0] = 9;
  actnode.dof[1] = 10;
  actnode.dof[2] = 11;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node4
  actnode.ID = 4;
  actnode.dof[0] = 12;
  actnode.dof[1] = 13;
  actnode.dof[2] = 14;
  actnode.x[0] = 1.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node5
  actnode.ID = 5;
  actnode.dof[0] = 15;
  actnode.dof[1] = 16;
  actnode.dof[2] = 17;
  actnode.x[0] = 0.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node6
  actnode.ID = 6;
  actnode.dof[0] = 18;
  actnode.dof[1] = 19;
  actnode.dof[2] = 20;
  actnode.x[0] = 0.5;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node7
  actnode.ID = 7;
  actnode.dof[0] = 21;
  actnode.dof[1] = 22;
  actnode.dof[2] = 23;
  actnode.x[0] = 0.5;
  actnode.x[1] = 0.5;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node8
  actnode.ID = 8;
  actnode.dof[0] = 24;
  actnode.dof[1] = 25;
  actnode.dof[2] = 26;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.5;
  actnode.x[2] = 0.0;
  node.push_back(actnode);

  // node9
  actnode.ID = 9;
  actnode.dof[0] = 27;
  actnode.dof[1] = 28;
  actnode.dof[2] = 29;
  actnode.x[0] = 0.5;
  actnode.x[1] = 0.0;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node10
  actnode.ID = 10;
  actnode.dof[0] = 30;
  actnode.dof[1] = 31;
  actnode.dof[2] = 32;
  actnode.x[0] = 0.5;
  actnode.x[1] = 0.5;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node11
  actnode.ID = 11;
  actnode.dof[0] = 33;
  actnode.dof[1] = 34;
  actnode.dof[2] = 35;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.5;
  actnode.x[2] = 1.0;
  node.push_back(actnode);

  // node12
  actnode.ID = 12;
  actnode.dof[0] = 36;
  actnode.dof[1] = 37;
  actnode.dof[2] = 38;
  actnode.x[0] = 0.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.5;
  node.push_back(actnode);

  // node13
  actnode.ID = 13;
  actnode.dof[0] = 39;
  actnode.dof[1] = 40;
  actnode.dof[2] = 41;
  actnode.x[0] = 1.0;
  actnode.x[1] = 0.0;
  actnode.x[2] = 0.5;
  node.push_back(actnode);

  // node14
  actnode.ID = 14;
  actnode.dof[0] = 42;
  actnode.dof[1] = 43;
  actnode.dof[2] = 44;
  actnode.x[0] = 0.0;
  actnode.x[1] = 1.0;
  actnode.x[2] = 0.5;
  node.push_back(actnode);
}
} // namespace test
} // namespace icarat
