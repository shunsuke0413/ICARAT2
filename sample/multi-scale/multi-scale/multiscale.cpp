#include "multiscale.hpp"

using namespace std;
using namespace icarat;
using namespace icarat::multiscale;
using namespace Eigen;

void outputmicro(string fname, FEM &fem, vector<mHomoElastic> &element,
                 vector<HomoNode> &node, ValueHomo &value);
void outputMacro(FEM &fem, vector<MHomoElastic> &element, vector<Node> &node);

int main(int argc, char *argv[])
{
  if(argc != 3)
  {
    std::cout << "usage: " << argv[0]
              << " <micro's .toml file> <macro's .toml file>\n";
    return 1;
  }

  auto t_start = chrono::system_clock::now(); // measure time
  int ID = 0; // macro ID used for localized analysis

  // set input file
  const auto mconfig = toml::parse(argv[1]);
  const auto Mconfig = toml::parse(argv[2]);

  //////////////////////////////////////////////////////////////
  ////////////////////////// set micro /////////////////////////
  //////////////////////////////////////////////////////////////
  Micro m;
  // set problem & mesh
  MeshIcarat<mHomoElastic, HomoNode> mmesh(m.fem, mconfig);
  mmesh.setParameter("structure");
  // set element
  for(int type = 0; type < mmesh.numeleTypes(); type++)
  {
    int numdof = mmesh.ne(type) * m.fem.dofnp;
    for(int nel = 0; nel < mmesh.nelx(type); nel++)
    {
      mHomoElastic actele(m.fem.voigt, mmesh.ne(type), mmesh.ipmax(type),
                          numdof);
      actele.setParameter(mconfig);
      m.element.push_back(actele);
    }
  }
  // set node
  m.node = vector<HomoNode>(m.fem.numnp);
  // input mesh to element & node
  mmesh.generate(m.element, m.node);
  auto pairs = mmesh.generateMasterSlavePairs(m.node);
  // boundary condition class
  BCIcarat<mHomoElastic, HomoNode> mbc(m.fem, m.element, m.node);
  // set periodic dirichlet condition
  m.fem.numeq = mbc.makeDOF_periodic(pairs);
  // set displacement
  ValueHomo mdis(m.fem.numeq, m.fem.voigt);
  array<double, 3> length = mmesh.length();
  Homogenization<mHomoElastic, HomoNode> mCompt(m.fem, m.element, m.node, mdis,
                                                pairs, length);

  //////////////////////////////////////////////////////////////
  ////////////////////////// set macro /////////////////////////
  //////////////////////////////////////////////////////////////
  Macro M;
  // set problem & mesh
  MeshIcarat<MHomoElastic, Node> Mmesh(M.fem, Mconfig);
  Mmesh.setParameter("structure");
  // set element
  for(int type = 0; type < Mmesh.numeleTypes(); type++)
  {
    int numdof = Mmesh.ne(type) * M.fem.dofnp;
    for(int nel = 0; nel < Mmesh.nelx(type); nel++)
    {
      MHomoElastic actele(M.fem.voigt, Mmesh.ne(type), Mmesh.ipmax(type),
                          numdof);
      actele.setParameter(Mconfig);
      M.element.push_back(actele);
    }
  }
  // set node
  M.node = vector<Node>(M.fem.numnp);
  // input mesh to element & node
  Mmesh.generate(M.element, M.node);
  // boundary condition class
  BCIcarat<MHomoElastic, Node> MBC(M.fem, M.element, M.node);
  // set dirichlet condition
  MBC.input_dirich(Mconfig);
  // get number of effective dofs
  M.fem.numeq = MBC.makeDOF();
  // set Macro force & displacement
  Force Mforce(M.fem.numeq);
  Value Mdis(M.fem.numeq);
  // set equivalent node load to Force class
  MBC.input_neumann(Mconfig, Mforce);
  LinearAnalysis<MHomoElastic, Node> MCompt(M.fem, M.element, M.node, Mforce,
                                            Mdis);

  //////////////////////start analysis//////////////////////

  // make square ring
  double half = 0.5;
  for(int i = 30; i < 70; i++)
    for(int j = 30; j < 70; j++)
      m.element[100 * i + j].setDesignS(half);

  // numerical material testing
  MatrixXd CH = mCompt.solveNMT(true);

  // get micro values
  for(auto &e : m.element)
    e.makeMicroValues(m.fem, m.node);

  string fname = "res_nmt.vtk";

  outputmicro(fname, m.fem, m.element, m.node, mdis);

  //  set CH to macro element
  for(auto &Me : M.element)
    Me.setDe(CH);

  // macro analysis
  MCompt.solve("Eigen_LDLT", 1);

  // get stress
  for(auto &e : M.element)
    e.makeMisesStress(M.fem, M.node);

  outputMacro(M.fem, M.element, M.node);

  // get macro strain & stress
  VectorXd Mstrain = M.element[ID].strain();
  VectorXd Mstress(M.fem.voigt);
  Mstress.setZero();

  // localization analysis
  mCompt.solveLocal(Mstrain, Mstress, true);

  // get micro values
  for(auto &e : m.element)
    e.makeMicroValues(m.fem, m.node);

  fname = "res_local.vtk";

  outputmicro(fname, m.fem, m.element, m.node, mdis);

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void outputmicro(string fname, FEM &fem, vector<mHomoElastic> &element,
                 vector<HomoNode> &node, ValueHomo &value)
{
  ofstream fout(fname);
  makeHeaderVTK(fout);

  // make grid
  setPointVTK(node, fout);
  setElementVTK(element, fout);
  setElementTypeVTK(fem.ndim, element, fout);

  // export vector data to node
  addPointVectorVTK("displacement", node, fout, true,
                    [&](Node &n) { return n.val; });

  // micro strain in each direction
  for(int d = 0; d < fem.voigt; d++)
    addPointVectorVTK("direction" + to_string(d), node, fout, false,
                      [&](HomoNode &n) { return n.nmtval[d]; });

  // export scalar data to element
  addElementScalarVTK("density", element, fout, true,
                      [&](mHomoElastic &e) { return e.design_s(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](mHomoElastic &e) { return e.mises(); });

  fout.close();
}

void outputMacro(FEM &fem, vector<MHomoElastic> &element, vector<Node> &node)
{
  ofstream fout("res_macro.vtk");
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
                      [&](MHomoElastic &e) { return e.design_s(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](MHomoElastic &e) { return e.mises(); });

  fout.close();
}
