#include "contact.hpp"
using namespace icarat;
using namespace std;
using namespace Eigen;

void output(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
            Optimize &optim);

int main(int argc, char *argv[])
{
  if(argc != 3)
  {
    std::cout << "usage: " << argv[0] << " <master .toml> <slave .toml>\n";
    return 1;
  }

  // measure time
  auto t_start = chrono::system_clock::now();
  int convergence = 1; // convergence judge
  int numConst = 1;    // number of constraint
  double V0 = 0.0;     // first volume

  // set input file
  const auto config_m = toml::parse(argv[1]);
  const auto config_s = toml::parse(argv[2]);

  // parameters
  double deltaT = toml::find<double>(config_m, "this_problem", "delta_t");
  int timestep = toml::find<int>(config_m, "this_problem", "timestep");
  double tolNR = toml::find<double>(config_m, "this_problem", "tolNR");
  int maxIter = toml::find<int>(config_m, "this_problem", "maxIter");
  string solvName = toml::find<string>(config_m, "this_problem", "solvName");

  // master setting
  FEM fem_m;
  MeshIcarat<LinearElastic, Node> mesh_m(fem_m, config_m);
  mesh_m.setParameter("structure");
  vector<LinearElastic> element_m;
  for(int type = 0; type < mesh_m.numeleTypes(); type++)
  {
    int numdof = mesh_m.ne(type) * fem_m.dofnp;
    for(int nel = 0; nel < mesh_m.nelx(type); nel++)
    {
      LinearElastic actele(fem_m.voigt, mesh_m.ne(type), mesh_m.ipmax(type),
                           numdof);
      actele.setParameter(config_m);
      element_m.push_back(actele);
    }
  }
  vector<Node> node_m(fem_m.numnp);
  mesh_m.generate(element_m, node_m);
  BCIcarat<LinearElastic, Node> BCm(fem_m, element_m, node_m);
  BCm.input_dirich(config_m);
  fem_m.numeq = BCm.makeDOF();
  Force force_m(fem_m.numeq);
  Value dis_m(fem_m.numeq);

  // slave setting
  FEM fem_s;
  MeshIcarat<LinearElastic, Node> mesh_s(fem_s, config_s);
  mesh_s.setParameter("structure");
  vector<LinearElastic> element_s;
  for(int type = 0; type < mesh_s.numeleTypes(); type++)
  {
    int numdof = mesh_s.ne(type) * fem_s.dofnp;
    for(int nel = 0; nel < mesh_s.nelx(type); nel++)
    {
      LinearElastic actele(fem_s.voigt, mesh_s.ne(type), mesh_s.ipmax(type),
                           numdof);
      actele.setParameter(config_s);
      element_s.push_back(actele);
    }
  }
  vector<Node> node_s(fem_s.numnp);
  mesh_s.generate(element_s, node_s);
  BCIcarat<LinearElastic, Node> BCs(fem_s, element_s, node_s);
  BCs.input_dirich(config_s);
  fem_s.numeq = BCs.makeDOF();
  Force force_s(fem_s.numeq);
  Value dis_s(fem_s.numeq);

  ContactAnalysis<LinearElastic, Node> compt(fem_m, element_m, node_m, force_m,
                                             dis_m, fem_s, element_s, node_s,
                                             force_s, dis_s);

  for(int time = 1; time <= timestep; time++)
  {

    cout << "------------ timestep no. " << time << " ------------" << endl;

    compt.solve(time, timestep, tolNR, maxIter, solvName, 1);

    // for(auto &e : element_m)
    //   e.makeMisesStress(fem_m, node_m);
    // for(auto &e : element_s)
    //   e.makeMisesStress(fem_s, node_s);

    output(fem_m, element_m, node_m, time);
    output(fem_s, element_s, node_s, time);
    exit(1);
  }

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec " << endl;

  return 0;
}

void output(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
            int time)
{
  string outFile = "res" + to_string(time) + ".vtk";
  ofstream fout(outFile);
  makeHeaderVTK(fout);

  // make grid
  setPointVTK(node, fout);
  setElementVTK(element, fout);
  setElementTypeVTK(fem.ndim, element, fout);

  // export vector data to node
  addPointVectorVTK("displacement", node, fout, true,
                    [&](Node &n) { return n.val; });

  // export scalar data to element
  addElementScalarVTK("density", element, fout, true,
                      [&](LinearElastic &e) { return e.design_s(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](LinearElastic &e) { return e.mises(); });

  fout.close();
}
