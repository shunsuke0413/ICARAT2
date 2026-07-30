#include "opt_mma_Helmholtz.hpp"
using namespace icarat;
using namespace icarat::helmholtz;
using namespace std;
using namespace Eigen;

void sensitivity(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
                 Optimize &optim, double V0);
void output(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
            Optimize &optim);

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <.toml file>\n";
    return 1;
  }

  // measure time
  auto t_start = chrono::system_clock::now();
  int convergence = 1; // convergence judge
  int numConst = 1;    // number of constraint
  double V0 = 0.0;     // first volume

  // set input file
  const auto config = toml::parse(argv[1]);

  // set problem
  FEM fem;
  MeshIcarat<LinearElastic, Node> mesh(fem, config);
  mesh.setParameter("structure");

  // set element
  vector<LinearElastic> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      LinearElastic actele(fem.voigt, mesh.ne(type), mesh.ipmax(type), numdof);
      actele.setParameter(config);
      element.push_back(actele);
    }
  }
  // set node
  vector<Node> node(fem.numnp);

  // input mesh to element & node
  mesh.generate(element, node);

  // boundary condition class
  BCIcarat<LinearElastic, Node> BC(fem, element, node);
  BC.input_dirich(config);

  // set dirichlet condition
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  Value dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  // set optimization's classes
  Optimize optim(fem.nelx, numConst, 1.0e-4, 1.0);
  MMASolver mma(fem.nelx, numConst);
  // filter radius
  double filterR = 2.0 * abs(node[0].x[0] - node[1].x[0]);
  FilterHelmholtz<LinearElastic, Node> filter(fem, element, node, filterR);

  // initial design field
  optim.optstep = 0;
  for(int i = 0; i < fem.nelx; i++)
    optim.design_s[i] = element[i].design_s();

  LinearAnalysis<LinearElastic, Node> compt(fem, element, node, force, dis);

  // optimization start
  while(convergence == 1)
  {
    optim.optstep++;

    cout << "============== optimization step ==============" << endl;
    cout << "step no. " << optim.optstep << endl;

    // filtering procedure
    auto fil_s = filter.densityFilter(optim.design_s, "Eigen_CG", 0, false);
    double beta = min(pow(2, (int)(optim.optstep / 40)), 16.0);
    cout << "beta" << beta << endl;
    auto thre_s = filter.threshold(fil_s, beta, 0.5);

    // set threshold design variable
    for(int i = 0; i < fem.nelx; i++)
      element[i].setDesignS(thre_s[i]);

    // structural analysis
    compt.solve("Eigen_LDLT", 1);

    // get  stress for output
    for(auto &e : element)
      e.makeMisesStress(fem, node);

    if(optim.optstep == 1)
    {
      V0 = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        V0 += element[i].design_s() * element[i].volume;
      if(V0 == 0.0)
        cerr << "V0 = 0. in optimization" << endl;
    }
    else
    {
      // get const_func
      optim.const_h[0] = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        optim.const_h[0] += element[i].design_s() * element[i].volume / V0;
      optim.const_h[0] -= 1;
      cout << "Constraint value:" << optim.const_h[0] << endl;
    }

    // get objective function
    optim.object_f = dis.val.dot(force.fext);
    cout << "function value = " << optim.object_f << endl;

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);

    output(fem, element, node, optim);

    convergence = optim.isConvergence(1.0e-7, 50);

    // analytical
    sensitivity(fem, element, node, optim, V0);
    auto thrediff = filter.thresholdDerivative(fil_s, beta, 0.5);
    for(int i = 0; i < fem.nelx; i++)
    {
      optim.dfds[i] *= thrediff[i];
      optim.dgds[i] *= thrediff[i];
    }
    optim.dfds = filter.densityFilter(optim.dfds, "Eigen_CG", 0, true);
    optim.dgds = filter.densityFilter(optim.dgds, "Eigen_CG", 0, true);

    // output sensivity graph
    if(optim.optstep == 1)
      exportGraph(optim.dfds, "objective.csv");

    // finite difference method
    if(toml::find<bool>(config, "this_problem", "fdm"))
      fdm(fem, element, node, force, dis, optim, filter);

    // Get updated design variables with MMA
    mma.Update(optim.design_s.data(), optim.dfds.data(), optim.const_h.data(),
               optim.dgds.data(), optim.xmin.data(), optim.xmax.data());

    // not convergence
    if(optim.optstep == 1000)
    {
      cerr << "not convergence in this optimization" << endl;
      exit(1);
    }
  }

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec " << endl;

  return 0;
}

void sensitivity(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
                 Optimize &optim, double V0)
{
  for(int nel = 0; nel < fem.nelx; nel++)
  {
    double dens = element[nel].design_s();
    double young1 = element[nel].young1();
    double young2 = element[nel].young2();
    double p = element[nel].pp();
    double young, dyoung;
    if(young2 >= young1)
    {
      young = (1.0 - pow(dens, p)) * young1 + pow(dens, p) * young2;
      dyoung = p * pow(dens, p - 1.0) * (young2 - young1);
    }
    else if(young2 < young1)
    {
      young = pow(1.0 - dens, p) * young1 + (1.0 - pow(1.0 - dens, p)) * young2;
      dyoung = (young2 - young1) * p * pow(1.0 - dens, p - 1.0);
    }

    double E = dyoung / young;

    Eigen::VectorXd disp = Eigen::VectorXd::Zero(element[nel].numdof);
    for(int i = 0; i < element[nel].ne; i++)
      for(int j = 0; j < fem.ndim; j++)
        disp.coeffRef(fem.ndim * i + j) = node[element[nel].nodeID[i]].val[j];

    optim.dfds[nel] = -E * disp.transpose() * element[nel].Ke() * disp;

    optim.dgds[nel] = element[nel].volume / V0;
  }
}

void output(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
            Optimize &optim)
{
  string outFile = "res" + to_string(optim.optstep) + ".vtk";
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

  // convert ascii to binary file using "meshio"
  string command = "meshio-binary " + outFile;
  auto tmp = system(command.c_str());
}
