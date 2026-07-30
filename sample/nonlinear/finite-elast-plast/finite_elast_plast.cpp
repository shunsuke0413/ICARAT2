#include "finite_elast_plast.hpp"

using namespace std;
using namespace icarat;

void output(FEM &fem, vector<FiniteMultipliEPUpdated> &element,
            vector<Node> &node, int loadStep);

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
  int numStep = toml::find<int>(config, "this_problem", "numstep");
  FEM fem;
  MeshIcarat<FiniteMultipliEPUpdated, Node> mesh(fem, config);
  mesh.setParameter("structure");

  // set element(unstructured mesh version)
  vector<FiniteMultipliEPUpdated> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      FiniteMultipliEPUpdated actele(fem.voigt, mesh.ne(type), mesh.ipmax(type),
                                     numdof);
      actele.setParameter(config);
      element.push_back(actele);
    }
  }

  // set node
  vector<Node> node(fem.numnp);
  // input mesh to element & node
  mesh.generate(element, node);

  // boundary condition class
  BCIcarat<FiniteMultipliEPUpdated, Node> BC(fem, element, node);

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  ValueNR dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);
  force.fext_org[6] = 1.0;
  force.fext_org[12] = 1.0;
  force.fext_org[18] = 1.0;
  force.fext_org[24] = 1.0;
  vector<double> factor = {0.0, 4000.0, 8000.0, 9000.0, 1.0e4, 1.1e4};
  NewtonRaphson<FiniteMultipliEPUpdated, Node> compt(fem, element, node, force,
                                                     dis);
  compt.initialization();

  /// Initialization of element internal variables
  for(int nel = 0; nel < fem.nelx; nel++)
    element[nel].initialization();

  for(int loadStep = 1; loadStep < numStep + 1; loadStep++)
  {
    force.fext = factor[loadStep] * force.fext_org;

    cout << "---------- load-step No. " << loadStep << " ----------" << endl;

    compt.solve(loadStep, numStep, 1.0e-8, 10, "Eigen_LDLT", 1);

    // make & stock stress
    for(auto &e : element)
      e.updateVariables(fem, node);

    output(fem, element, node, loadStep);
  }

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void output(FEM &fem, vector<FiniteMultipliEPUpdated> &element,
            vector<Node> &node, int loadStep)
{
  ofstream fout("res" + to_string(loadStep) + ".vtk");
  makeHeaderVTK(fout);

  // make grid
  setPointVTK(node, fout);
  setElementVTK(element, fout);
  setElementTypeVTK(fem.ndim, element, fout);

  // export vector data to node
  addPointVectorVTK("displacement", node, fout, true,
                    [&](Node &n) { return n.val; });
  addPointVectorVTK("reaction_force", node, fout, false,
                    [&](Node &n) { return n.reaction; });

  // export scalar data to element
  addElementScalarVTK("ID", element, fout, true,
                      [&](FiniteMultipliEPUpdated &e) { return e.ID; });

  addElementScalarVTK("density", element, fout, false,
                      [&](FiniteMultipliEPUpdated &e) { return e.design_s(); });

  addElementScalarVTK("ep", element, fout, false,
                      [&](FiniteMultipliEPUpdated &e)
                      {
                        double ep = 0.0;
                        for(int ip = 0; ip < e.ipmax; ip++)
                          ep += e.ep()[ip] / e.ipmax;
                        return ep;
                      });

  fout.close();
}