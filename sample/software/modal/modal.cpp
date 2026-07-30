#include "modal.hpp"

using namespace std;
using namespace icarat;

void output(FEM &fem, vector<ModalLinearElastic> &element, vector<Node> &node,
            ValueEigen &value, int numeigen, string project);

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <project folder> \n";
    return 1;
  }

  // measure time
  auto t_start = chrono::system_clock::now();

  // set input file
  string project = argv[1];
  string tomlpath = project + "/TEXTfiles/tmp.toml";
  string vtkpath = project + "/VTKfiles/tmp.vtk";
  const auto config = toml::parse(tomlpath);

  auto density0 = toml::find<vector<double>>(config, "material", "density");

  // shift invert alnordi parameter
  double sigma = toml::find<double>(config, "this_problem", "shift");
  int numeigen = toml::find<int>(config, "this_problem", "numeigen");

  // set problem
  FEM fem;
  MeshVTK<ModalLinearElastic, Node> mesh(fem, vtkpath, config);
  mesh.setParameter("structure");
  auto design_s = mesh.importElementScalars("density");

  // set element(unstructured mesh version)
  vector<ModalLinearElastic> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      double density = 0.0;
      if(design_s[nel] >= 0.5)
        density = density0[0];
      else
        density = 1.0e-4 * density0[0];
      ModalLinearElastic actele(fem.voigt, mesh.ne(type), mesh.ipmax(type),
                                numdof, density);
      actele.setParameter(config);
      element.push_back(actele);
    }
  }

  // set node
  vector<Node> node(fem.numnp);

  // input mesh to element & node
  mesh.generate(element, node);

  // boundary condition class
  BCPoint<ModalLinearElastic, Node> BC(fem, element, node, project);

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  ValueEigen dis(fem.numeq, numeigen);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  auto IDs = getDesignNonDesignIDs(fem, element, node, config, project);

  int paraconv = 3 * numeigen;
  GeneralEigen<ModalLinearElastic, Node> compt(fem, element, node, dis,
                                               numeigen, paraconv, sigma);

  for(int i = 0; i < fem.nelx; i++)
    element[i].setDesignS(design_s[i]);

  compt.solve();

  output(fem, element, node, dis, numeigen, project);

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void output(FEM &fem, vector<ModalLinearElastic> &element, vector<Node> &node,
            ValueEigen &value, int numeigen, string project)
{
  string asciiname = project + "/VTKfiles/output_modal.vtk";

  ofstream fout(asciiname);
  makeHeaderVTK(fout);

  // make grid
  setPointVTK(node, fout);
  setElementVTK(element, fout);
  setElementTypeVTK(fem.ndim, element, fout);

  // export vector data to node
  addPointScalarVTK("ID", node, fout, true, [&](Node &n) { return n.ID; });
  // make eigen mode's output
  for(int ei = 0; ei < numeigen; ei++)
    addPointVectorVTK("eigenmode" + to_string(numeigen - ei), node, fout, false,
                      [&](Node &n)
                      {
                        vector<double> three(3);
                        for(int i = 0; i < fem.dofnp; i++)
                        {
                          // dirichlet is active.
                          if(n.dof[i] >= fem.numeq)
                            three[i] = n.dirich->val[i];
                          else
                            three[i] = value.evectors[ei](n.dof[i]);
                        }
                        return three;
                      });

  // export scalar data to element
  addElementScalarVTK("density", element, fout, true,
                      [&](ModalLinearElastic &e) { return e.density(); });

  addElementScalarVTK("design_s", element, fout, false,
                      [&](ModalLinearElastic &e) { return e.design_s(); });

  addElementScalarVTK("designable", element, fout, false,
                      [&](ModalLinearElastic &e) { return e.isDesignable; });

  addElementScalarVTK("young", element, fout, false,
                      [&](ModalLinearElastic &e) { return e.young(); });

  fout.close();
}