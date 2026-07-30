#include "finite_arc_length.hpp"
using namespace icarat;
using namespace std;

void output(FEM &fem, vector<FiniteMooneyTotal> &element, vector<Node> &node,
            int loadStep);

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
  int numstep = toml::find<int>(config, "this_problem", "numstep");

  double tolNR = 1.0e-6;
  int iterNR = 100;

  // set problem
  FEM fem;
  MeshIcarat<FiniteMooneyTotal, Node> mesh(fem, config);
  mesh.setParameter("structure");

  // set element
  vector<FiniteMooneyTotal> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      FiniteMooneyTotal actele(fem.voigt, mesh.ne(type), mesh.ipmax(type),
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
  BCIcarat<FiniteMooneyTotal, Node> BC(fem, element, node);

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  ValueNR dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  ArcLength<FiniteMooneyTotal, Node> compt(fem, element, node, force, dis);

  compt.initialization();
  int step = 0;
  while(1) // arc length method loop
  {
    step++;

    cout << "-------- increment-step No. " << step << " --------" << endl;

    compt.solve(tolNR, iterNR, "Eigen_LDLT", step, 0);

    for(auto &e : element)
      e.makeMisesStress(fem, node);

    output(fem, element, node, step);

    cout << "lambda: " << compt.lambda_ << endl;
    if(compt.lambda_ >= 1.0)
    {
      break;
    }
  }

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void output(FEM &fem, vector<FiniteMooneyTotal> &element, vector<Node> &node,
            int loadStep)
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

  addPointVectorVTK("reaction", node, fout, false,
                    [&](Node &n) { return n.reaction; });

  // export scalar data to element
  addElementScalarVTK("density", element, fout, true,
                      [&](FiniteMooneyTotal &e) { return e.design_s(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](FiniteMooneyTotal &e) { return e.mises(); });

  fout.close();
}