#include "modal.hpp"

using namespace std;
using namespace icarat;

void output(FEM &fem, vector<ModalLinearElastic> &element, vector<Node> &node,
            ValueEigen &value, int numeigen);

int main(int argc, char *argv[])
{
  if(argc != 3)
  {
    std::cout << "usage: " << argv[0] << " <.toml file> <.mdpa file>\n";
    return 1;
  }

  // measure time
  auto t_start = chrono::system_clock::now();

  // set input file
  const auto config = toml::parse(argv[1]);
  vector<string> mdpa = read_file(argv[2]);

  double density = toml::find<double>(config, "material", "density");

  // shift invert alnordi parameter
  double sigma = toml::find<double>(config, "this_problem", "shift");
  int numeigen = toml::find<int>(config, "this_problem", "numeigen");

  // set problem
  FEM fem;
  MeshMDPA<ModalLinearElastic, Node> mesh(fem, mdpa, config);
  mesh.setParameter("structure");

  // set element(unstructured mesh version)
  vector<ModalLinearElastic> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
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
  BCIcarat<ModalLinearElastic, Node> BC(fem, element, node);

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  ValueEigen dis(fem.numeq, numeigen);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  int paraconv = 3 * numeigen;
  GeneralEigen<ModalLinearElastic, Node> compt(fem, element, node, dis,
                                               numeigen, paraconv, sigma);
  compt.solve();

  output(fem, element, node, dis, numeigen);

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void output(FEM &fem, vector<ModalLinearElastic> &element, vector<Node> &node,
            ValueEigen &value, int numeigen)
{
  ofstream fout("res_mapped.vtk");
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
  addElementScalarVTK("ID", element, fout, true,
                      [&](ModalLinearElastic &e) { return e.ID; });

  addElementScalarVTK("design_s", element, fout, false,
                      [&](ModalLinearElastic &e) { return e.design_s(); });

  fout.close();
}