#include "explicit_method.hpp"
using namespace icarat;
using namespace std;

void output(FEM &fem, vector<DynamicElastic> &element, vector<Node> &node,
            int time);

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <.toml file>\n";
    return 1;
  }

  // measure time
  auto t_start = chrono::system_clock::now();

  // set input file
  const auto config = toml::parse(argv[1]);
  // parameters
  double deltaT = toml::find<double>(config, "this_problem", "delta_t");
  int timestep = toml::find<int>(config, "this_problem", "timestep");
  // material density (kg/mm3)
  double density = toml::find<double>(config, "material", "density");

  // set problem
  FEM fem;
  MeshIcarat<DynamicElastic, Node> mesh(fem, config);
  mesh.setParameter("structure");

  // set element
  vector<DynamicElastic> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      DynamicElastic actele(fem.voigt, mesh.ne(type), mesh.ipmax(type), numdof,
                            density);
      actele.setParameter(config);
      element.push_back(actele);
    }
  }

  // set node
  vector<Node> node(fem.numnp);

  // input mesh to element & node
  mesh.generate(element, node);

  // boundary condition class
  BCIcarat<DynamicElastic, Node> BC(fem, element, node);

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  ValueDyna dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  // ImplicitDyna<ElastDyna, Node> compt(fem, element, node, force, dis);
  ExplicitDyna<DynamicElastic, Node> compt(fem, element, node, force, dis);

  compt.initialization();

  for(int time = 1; time < timestep + 1; time++)
  {
    if((time - 1) % 10 == 0)
      force.fext = time * force.fext_org / timestep;

    cout << "------------ timestep no. " << time << " ------------" << endl;

    compt.solve(time, timestep, deltaT);

    for(auto &e : element)
      e.makeMisesStress(fem, node);

    output(fem, element, node, time);
  }

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void output(FEM &fem, vector<DynamicElastic> &element, vector<Node> &node,
            int time)
{
  ofstream fout("res" + to_string(time) + ".vtk");
  makeHeaderVTK(fout);

  // make grid
  setPointVTK(node, fout);
  setElementVTK(element, fout);
  setElementTypeVTK(fem.ndim, element, fout);

  // export vector data to node
  addPointVectorVTK("displacement", node, fout, true,
                    [&](Node &n) { return n.val; });

  // export scalar data to element
  addElementScalarVTK("design_s", element, fout, true,
                      [&](DynamicElastic &e) { return e.design_s(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](DynamicElastic &e) { return e.mises(); });

  fout.close();
}