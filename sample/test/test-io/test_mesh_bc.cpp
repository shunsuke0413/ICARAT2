#include <cmath>
#include <gtest/gtest.h>
#include <io/boundary_condition_gid.hpp>
#include <io/boundary_condition_icarat.hpp>
#include <io/export_vtk.hpp>
#include <io/mesh_mdpa.hpp>
#include <io/mesh_vtk.hpp>
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace icarat;

TEST(MDPAIOTest, dirich_neumann_2d)
{
  // set gid folder, mdpa file toml file's path
  string gidPath = "test-io/2d.gid";
  GIDParser parser(gidPath);
  vector<string> mdpa = read_file(gidPath + "/" + parser.mdpaName() + ".mdpa");
  const auto config = toml::parse("test-io/2d.toml");

  FEM fem;
  MeshMDPA<Element, Node> mesh(fem, mdpa, config);
  mesh.setParameter("structure");

  vector<Node> node(fem.numnp);
  vector<Element> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      Element actele;
      actele.ne = mesh.ne(type);
      actele.ipmax = mesh.ipmax(type);
      actele.numdof = mesh.ne(type) * fem.dofnp;
      element.push_back(actele);
    }
  }

  mesh.generate(element, node);

  BCGID<Element, Node> bc(fem, element, node, gidPath);

  bc.input_dirich(config, mdpa);

  fem.numeq = bc.makeDOF();

  Force force(fem.numeq);

  bc.input_neumann(config, mdpa, force);

  VectorXd answer(fem.numeq); // =size is 12
  // equivalent load
  answer << 0, 0, 0, 0, 0.25, 0, 0, 0, 0.5, 0, 0.25, 0;

  EXPECT_LT((force.fext_org - answer).norm(), 1.0e-10);
}

TEST(MDPAIOTest, dirich_neumann_3d)
{
  // set gid folder, mdpa file config file's path
  string gidPath = "test-io/3d.gid";
  GIDParser parser(gidPath);
  vector<string> mdpa = read_file(gidPath + "/" + parser.mdpaName() + ".mdpa");
  const auto config = toml::parse("test-io/3d.toml");

  FEM fem;
  MeshMDPA<Element, Node> mesh(fem, mdpa, config);
  mesh.setParameter("structure");

  vector<Node> node(fem.numnp);
  vector<Element> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      Element actele;
      actele.ne = mesh.ne(type);
      actele.ipmax = mesh.ipmax(type);
      actele.numdof = mesh.ne(type) * fem.dofnp;
      element.push_back(actele);
    }
  }

  // 1.set mesh
  mesh.generate(element, node);

  BCGID<Element, Node> bc(fem, element, node, gidPath);

  bc.input_dirich(config, mdpa);

  fem.numeq = bc.makeDOF();

  Force force(fem.numeq);

  bc.input_neumann(config, mdpa, force);

  VectorXd answer(fem.numeq);
  // tetra element's equivalent nodal load
  answer << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.0 / 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 6.0, 0.0, 0.0,
      1.0 / 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 6.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.0 / 24.0, 0.0, 0.0, 1.0 / 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 6.0, 0.0,
      0.0, 1.0 / 6.0, 0.0, 0.0, 1.0 / 24.0, 0.0, 0.0;

  EXPECT_LT((force.fext_org - answer).norm(), 1.0e-10);
}

TEST(VTKIOTest, dirich_neumann_2d)
{
  const auto config = toml::parse("test-io/2d.toml");
  string vtk = "test-io/2d.vtk.test";
  FEM fem;
  MeshVTK<Element, Node> mesh(fem, vtk, config);
  mesh.setParameter("structure");

  vector<Node> node(fem.numnp);
  vector<Element> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      Element actele;
      actele.ne = mesh.ne(type);
      actele.ipmax = mesh.ipmax(type);
      actele.numdof = mesh.ne(type) * fem.dofnp;
      element.push_back(actele);
    }
  }

  // 1.set mesh
  mesh.generate(element, node);

  auto disp = mesh.importPointVectors("displacement");
  auto density = mesh.importElementScalars("density");
  vector<double> disp_answer{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4,
                             4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8};
  vector<double> density_answer{0, 1, 2, 3};

  double err1 = 0.0, err2 = 0.0;
  for(int i = 0; i < (int)disp.size(); i++)
    err1 += disp[i] - disp_answer[i];
  for(int i = 0; i < (int)density.size(); i++)
    err2 += density[i] - density_answer[i];

  EXPECT_LT(err1, 1.0e-10);
  EXPECT_LT(err2, 1.0e-10);

  BCIcarat<Element, Node> bc(fem, element, node);

  bc.input_dirich(config);

  fem.numeq = bc.makeDOF();

  Force force(fem.numeq);

  bc.input_neumann(config, force);
}

TEST(ICARATIOTest, dirich_neumann_2d)
{
  const auto config = toml::parse("test-io/quad8.toml");

  FEM fem;
  MeshIcarat<Element, Node> mesh(fem, config);
  mesh.setParameter("heat");

  vector<Node> node(fem.numnp);
  vector<Element> element;

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    Element actele;
    actele.ne = mesh.ne(0);
    actele.ipmax = mesh.ipmax(0);
    actele.numdof = mesh.ne(0) * fem.dofnp;
    element.push_back(actele);
  }

  // 1.set mesh
  mesh.generate(element, node);

  BCIcarat<Element, Node> bc(fem, element, node);

  // bc.input_dirich(config);

  fem.numeq = bc.makeDOF();

  Force force(fem.numeq);

  bc.input_neumann(config, force);

  // quadratic element's load
  VectorXd answer(fem.numeq);
  answer << 0, 0, 1.0 / 6.0, 0, 2.0 / 3.0, 0, 0, 1.0 / 6.0;
  EXPECT_LT((force.fext_org - answer).norm(), 1.0e-10);
}

TEST(ICARATIOTest, dirich_neumann_3d)
{
  const auto config = toml::parse("test-io/hexa20.toml");

  FEM fem;
  MeshIcarat<Element, Node> mesh(fem, config);
  mesh.setParameter("heat");

  vector<Node> node(fem.numnp);
  vector<Element> element;

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    Element actele;
    actele.ne = mesh.ne(0);
    actele.ipmax = mesh.ipmax(0);
    actele.numdof = mesh.ne(0) * fem.dofnp;
    element.push_back(actele);
  }

  // 1.set mesh
  mesh.generate(element, node);

  BCIcarat<Element, Node> bc(fem, element, node);

  // bc.input_dirich(config);

  fem.numeq = bc.makeDOF();

  Force force(fem.numeq);

  bc.input_neumann(config, force);

  // quadratic element's load
  VectorXd answer(fem.numeq);
  answer << 0, 0, -1.0 / 12.0, 0, 1.0 / 3.0, 0, 0, -1.0 / 12.0, 0, 1.0 / 3.0, 0,
      1.0 / 3.0, 0, 0, -1.0 / 12.0, 0, 1.0 / 3.0, 0, 0, -1.0 / 12.0;

  EXPECT_LT((answer - force.fext_org).norm(), 1.0e-10);
}
