#include "opt_multiscale.hpp"
using namespace std;
using namespace icarat;
using namespace icarat::opt_multiscale;
using namespace Eigen;

void outputmicro(FEM &fem, vector<mHomoElastic> &element,
                 vector<HomoNode> &node, ValueHomo &value, int optstep);
void outputMacro(FEM &fem, vector<MHomoElastic> &element, vector<Node> &node,
                 int optstep);

int main(int argc, char *argv[])
{
  if(argc != 3)
  {
    std::cout << "usage: " << argv[0]
              << " <micro's .toml file> <macro's .toml file>\n";
    return 1;
  }

  auto t_start = chrono::system_clock::now(); // measure time
  int numConst = 1;                           // number of constraint
  double xmin = 1.0e-4, xmax = 1.0; // upper & lower bound of design variable
  double V0 = 0.0;                  // initial volume
  int convergence = 1;              // convergence checker

  // set input file
  const auto mconfig = toml::parse(argv[1]);
  const auto Mconfig = toml::parse(argv[2]);

  //////////////////////////////////////////////////////////////
  ////////////////////////// micro /////////////////////////////
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
  ////////////////////////// Macro /////////////////////////////
  //////////////////////////////////////////////////////////////
  Macro M;
  // set problem & mesh
  MeshIcarat<MHomoElastic, Node> MMesh(M.fem, Mconfig);
  MMesh.setParameter("structure");
  // set element
  for(int type = 0; type < MMesh.numeleTypes(); type++)
  {
    int numdof = MMesh.ne(type) * M.fem.dofnp;
    for(int nel = 0; nel < MMesh.nelx(type); nel++)
    {
      MHomoElastic actele(M.fem.voigt, MMesh.ne(type), MMesh.ipmax(type),
                          numdof);
      actele.setParameter(Mconfig);
      M.element.push_back(actele);
    }
  }
  // set node
  M.node = vector<Node>(M.fem.numnp);
  // input mesh to element & node
  MMesh.generate(M.element, M.node);
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

  // optimization setting
  Optimize optim(m.fem.nelx, numConst, xmin, xmax);
  optim.optstep = 0;
  for(int i = 0; i < m.fem.nelx; i++)
    optim.design_s[i] = m.element[i].design_s();
  optim.design_s =
      mmesh.setDesignHomo(m.element, m.node, m.element[0].design_s(), 0.2, 0.3);

  MMASolver mma(m.fem.nelx, numConst);
  double filterR = 3.0 * abs(m.node[0].x[0] - m.node[1].x[0]);
  Filter<mHomoElastic, HomoNode> filter(m.fem, m.element, m.node, filterR);

  //////////////////////////////////////////////////////////
  /////////////////////// start TO /////////////////////////
  //////////////////////////////////////////////////////////
  while(convergence == 1)
  {
    optim.optstep++;

    cout << "============== optimization step ==============" << endl;
    cout << "step no. " << optim.optstep << endl;

    // filtering
    auto design_filtered = filter.densityFilter(optim.design_s, false);
    double beta = min(pow(2, (int)(optim.optstep / 40)), 16.0);
    cout << "beta" << beta << endl;
    double T = 0.5;
    auto design_threshold = filter.threshold(design_filtered, beta, T);
    auto threDeriv = filter.thresholdDerivative(design_filtered, beta, T);

    for(int i = 0; i < m.fem.nelx; i++)
      m.element[i].setDesignS(design_threshold[i]);

    // numerical material testing
    MatrixXd CH = mCompt.solveNMT(true);
    for(int i = 0; i < m.fem.nelx; i++)
      m.element[i].makeMicroValues(m.fem, m.node);

    outputmicro(m.fem, m.element, m.node, mdis, optim.optstep);

    // set CH to macro element
    for(int i = 0; i < M.fem.nelx; i++)
      M.element[i].setDe(CH);

    // macro analysis
    MCompt.solve("Eigen_LDLT", 1);

    // get stress
    for(int i = 0; i < M.fem.nelx; i++)
      M.element[i].makeMisesStress(M.fem, M.node);

    // const func
    if(optim.optstep == 1)
    {
      V0 = 0.0;
      for(int i = 0; i < m.fem.nelx; i++)
        V0 += optim.design_s[i] * m.element[i].volume;
      if(V0 == 0.0)
      {
        cerr << "V0 = 0. in optimization" << endl;
        exit(1);
      }

      optim.const_h[0] = 0.0;
      cout << "first volume: " << V0 << endl;
    }
    else
    {
      optim.const_h[0] = 0.0;
      for(int i = 0; i < m.fem.nelx; i++)
        optim.const_h[0] += optim.design_s[i] * m.element[i].volume;

      optim.const_h[0] = optim.const_h[0] / V0 - 1.0;

      cout << "const value: " << optim.const_h[0] << endl;
    }

    // get objective function
    optim.object_f = Mdis.val.dot(Mforce.fext_org);
    cout << "function value: " << optim.object_f << endl;

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);
    outputMacro(M.fem, M.element, M.node, optim.optstep);

    convergence = optim.isConvergence(1.0e-6, 200);

    // analytical
    sensitivity(M, Mforce, Mdis, m, optim, V0);

    // filtering sensitivity
    for(int i = 0; i < m.fem.nelx; i++)
    {
      optim.dfds[i] = threDeriv[i] * optim.dfds[i];
      optim.dgds[i] = threDeriv[i] * optim.dgds[i];
    }
    optim.dfds = filter.densityFilter(optim.dfds, true);
    optim.dgds = filter.densityFilter(optim.dgds, true);

    // output sensivity graph
    exportGraph(optim.dfds, "objective.csv");

    // fdm
    if(toml::find<bool>(Mconfig, "this_problem", "fdm"))
      fdm(M, m, Mforce, Mdis, mdis, optim, filter, pairs, length);

    // Get updated design variables with OC
    mma.Update(optim.design_s.data(), optim.dfds.data(), optim.const_h.data(),
               optim.dgds.data(), optim.xmin.data(), optim.xmax.data());

    // not convergence
    if(optim.optstep == 1000)
    {
      cerr << "not convergence in this optimization" << endl;
      convergence = 0;
    }
  }

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void outputmicro(FEM &fem, vector<mHomoElastic> &element,
                 vector<HomoNode> &node, ValueHomo &value, int optstep)
{
  ofstream fout("res_micro" + to_string(optstep) + ".vtk");
  makeHeaderVTK(fout);

  // make grid
  setPointVTK(node, fout);
  setElementVTK(element, fout);
  setElementTypeVTK(fem.ndim, element, fout);

  // export vector data to node
  addPointScalarVTK("ID", node, fout, true, [&](HomoNode &n) { return n.ID; });

  // micro strain in each direction
  for(int d = 0; d < fem.voigt; d++)
    addPointVectorVTK("direction" + to_string(d), node, fout, false,
                      [&](HomoNode &n) { return n.nmtval[d]; });

  // export scalar data to element
  addElementScalarVTK("density", element, fout, true,
                      [&](mHomoElastic &e) { return e.design_s(); });

  addElementScalarVTK("s_energy", element, fout, false,
                      [&](mHomoElastic &e) { return e.sEnergy().sum(); });

  fout.close();
}

void outputMacro(FEM &fem, vector<MHomoElastic> &element, vector<Node> &node,
                 int optstep)
{
  ofstream fout("res_macro" + to_string(optstep) + ".vtk");
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

  fout.close();
}
