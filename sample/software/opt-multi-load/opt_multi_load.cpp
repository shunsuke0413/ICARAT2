#include "opt_multi_load.hpp"

using namespace std;
using namespace icarat;
using namespace Eigen;

void sensitivity(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
                 Optimize &optim, double V0);
void output(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
            Optimize &optim, vector<vector<array<double, 3>>> &disps,
            vector<vector<double>> &misess, string project);

int main(int argc, char *argv[])
{

  if(argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <project folder's path>\n";
    std::cout << "This program assumes the use of TOP-MVP software.\n";
    return 1;
  }

  // measure time
  auto t_start = chrono::system_clock::now();
  int convergence = 1; // convergence judge
  int numConst = 1;    // number of constraint
  double V0 = 0.0;     // first volume
  double T = 0.5;

  // set input file
  string project = argv[1];
  string tomlpath = project + "/TEXTfiles/tmp.toml";
  string vtkpath = project + "/VTKfiles/tmp.vtk";
  const auto config = toml::parse(tomlpath);
  // set problem
  FEM fem;
  // set problem parameters
  int numlc = toml::find<int>(config, "boundary_condition", "num_loadcase");

  MeshVTK<LinearElastic, Node> mesh(fem, vtkpath, config);
  // MeshMDPA<LinearElastic, Node> mesh(fem, inputMDPA, config);
  mesh.setParameter("structure");

  // set element(unstructured mesh version)
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

  //  boundary condition class
  //  BCIcarat<LinearElastic, Node> BC(fem, element, node);
  BCPoint<LinearElastic, Node> BC(fem, element, node, project);

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force forcetmp(fem.numeq);
  vector<Force> force(numlc, forcetmp);
  Value dis(fem.numeq);

  // set equivalent node load to Force class
  for(int i = 0; i < numlc; i++)
  {
    std::cout << "generating case " + to_string(i + 1) + "'s load..." << endl;
    BC.input_neumann(config, force[i]);
  }

  // get pair of design area's ID & non-design area's ID
  auto IDs = getDesignNonDesignIDs(fem, element, node, config, project);

  // set optimization's classes
  Optimize optim(fem.nelx, numConst, 1.0e-4, 1.0);
  MMASolver mma(fem.nelx, 1);
  double filterR = toml::find<double>(config, "this_problem", "filterR");
  FilterHelmholtz<LinearElastic, Node> filter(fem, element, node, filterR);

  // initialization
  optim.optstep = 0;
  for(int i = 0; i < fem.nelx; i++)
    optim.design_s[i] = element[i].design_s();

  // optimization start
  while(convergence == 1)
  {
    optim.optstep++;

    std::cout << "============== optimization step ==============" << endl;
    std::cout << "step no. " << optim.optstep << endl;

    auto design_filtered =
        filter.densityFilter(optim.design_s, "Eigen_CG", 1, false);

    /// set threshold function
    double beta = min(pow(2, (int)(optim.optstep / 40)), 16.0);
    std::cout << "beta: " << beta << endl;
    auto design_threshold = filter.threshold(design_filtered, beta, T);
    auto threDeriv = filter.thresholdDerivative(design_filtered, beta, T);

    // set threshold design variable
    for(int i = 0; i < fem.nelx; i++)
      element[i].setDesignS(design_threshold[i]);

    // initialization
    optim.object_f = 0.0;
    fill(optim.dfds.begin(), optim.dfds.end(), 0.0);
    vector<vector<double>> misess(numlc);
    vector<vector<array<double, 3>>> disps(numlc);
    for(int lc = 0; lc < numlc; lc++)
    {
      std::cout << "computing load case " + to_string(lc + 1) + "..." << endl;
      LinearAnalysis<LinearElastic, Node> compt(fem, element, node, force[lc],
                                                dis);

      compt.solve("Eigen_CG", 1);

      // restore displacement & stress
      misess[lc].resize(fem.nelx);
      for(int i = 0; i < fem.nelx; i++)
      {
        element[i].makeMisesStress(fem, node);
        misess[lc][i] = element[i].mises();
      }

      disps[lc].resize(fem.numnp);
      for(int np = 0; np < fem.numnp; np++)
        for(int j = 0; j < 3; j++)
          disps[lc][np][j] = node[np].val[j];

      // get objective function
      optim.object_f += dis.val.dot(force[lc].fext);
      std::cout << "function value = " << optim.object_f << endl;

      // get const_func
      if(optim.optstep == 1)
      {
        V0 = 0.0;
        for(int i = 0; i < (int)IDs.first.size(); i++)
          V0 += element[IDs.first[i]].design_s() * element[IDs.first[i]].volume;
        if(V0 == 0.0)
          cerr << "V0 = 0. in optimization" << endl;

        optim.const_h[0] = 0.0;
      }
      else
      {
        optim.const_h[0] = 0.0;
        for(int i = 0; i < (int)IDs.first.size(); i++)
          optim.const_h[0] += element[IDs.first[i]].design_s() *
                              element[IDs.first[i]].volume / V0;

        optim.const_h[0] += -1.0;

        cout << "const value: " << optim.const_h[0] << endl;
      }

      // analytical
      sensitivity(fem, element, node, optim, V0);
    }

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);
    if((optim.optstep - 1) % 10 == 0)
      output(fem, element, node, optim, disps, misess, project);

    convergence = optim.isConvergence(1.0e-6, 50);

    for(int i = 0; i < fem.nelx; i++)
    {
      optim.dfds[i] *= threDeriv[i];
      optim.dgds[i] *= threDeriv[i];
    }

    // filtering sensitivity
    optim.dfds = filter.densityFilter(optim.dfds, "Eigen_CG", 0, true);
    optim.dgds = filter.densityFilter(optim.dgds, "Eigen_CG", 0, true);

    // In the non-design domain, the sensitivity is 0
    for(int i = 0; i < (int)IDs.second.size(); i++)
    {
      optim.dfds[IDs.second[i]] = 0.0;
      optim.dgds[IDs.second[i]] = 0.0;
    }

    // Get updated design variables with MMA
    mma.Update(optim.design_s.data(), optim.dfds.data(), optim.const_h.data(),
               optim.dgds.data(), optim.xmin.data(), optim.xmax.data());

    // not convergence
    if(optim.optstep == 200)
    {
      cerr << "optimization step reached 200 step." << endl;
      convergence = 0;
    }
  }

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  std::cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void sensitivity(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
                 Optimize &optim, double V0)
{
  assert(element[0].young2() > element[0].young1());

  for(int nel = 0; nel < fem.nelx; nel++)
  {
    double dens = element[nel].design_s();
    double young1 = element[nel].young1();
    double young2 = element[nel].young2();
    double p = element[nel].pp();

    double young = (1.0 - pow(dens, p)) * young1 + pow(dens, p) * young2;
    double dyoung = p * pow(dens, p - 1.0) * (young2 - young1);

    double E = dyoung / young;

    Eigen::VectorXd disp = Eigen::VectorXd::Zero(element[nel].numdof);
    for(int i = 0; i < element[nel].ne; i++)
      for(int j = 0; j < fem.ndim; j++)
        disp.coeffRef(fem.ndim * i + j) = node[element[nel].nodeID[i]].val[j];

    optim.dfds[nel] += -E * disp.transpose() * element[nel].Ke() * disp;

    optim.dgds[nel] = element[nel].volume / V0;
  }

  // output sensivity graph
  if(optim.optstep == 1)
    exportGraph(optim.dfds, "objective.csv");
}

void output(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
            Optimize &optim, vector<vector<array<double, 3>>> &disps,
            vector<vector<double>> &misess, string project)
{
  string asciiname =
      project + "/VTKfiles/output" + to_string(optim.optstep) + ".vtk";
  ofstream fout(asciiname);
  makeHeaderVTK(fout);

  // make grid
  setPointVTK(node, fout);
  setElementVTK(element, fout);
  setElementTypeVTK(fem.ndim, element, fout);

  // export vector data to node
  for(int i = 0; i < disps.size(); i++)
  {
    bool isFirst;
    if(i == 0)
      isFirst = true;
    else
      isFirst = false;
    addPointVectorVTK("displacement_case" + to_string(i), node, fout, isFirst,
                      [&](Node &n) { return disps[i][n.ID]; });
  }
  // export scalar data to element
  addElementScalarVTK("density", element, fout, true,
                      [&](LinearElastic &e) { return e.design_s(); });

  for(int i = 0; i < misess.size(); i++)
    addElementScalarVTK("mises_case" + to_string(i), element, fout, false,
                        [&](LinearElastic &e) { return misess[i][e.ID]; });

  fout.close();
}