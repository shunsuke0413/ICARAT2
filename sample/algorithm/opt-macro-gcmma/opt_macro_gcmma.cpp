#include "opt_macro_gcmma.hpp"

using namespace std;
using namespace icarat;

void sensitivity(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
                 Optimize &optim, vector<double> &threDiff);
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

  // set problem & mesh
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

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  Value dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  // optimization class
  Optimize optim(fem.nelx, numConst, 1.0e-4, 1.0);
  GCMMASolver gcmma(fem.nelx, numConst);
  // filter radius (is 2-cell length)
  double filterR = 2.0 * mesh.length()[0] / mesh.div()[0];
  Filter<LinearElastic, Node> filter(fem, element, node, filterR);

  // fisrt design field
  for(int nel = 0; nel < fem.nelx; nel++)
    optim.design_s[nel] = element[nel].design_s();

  LinearAnalysis<LinearElastic, Node> compt(fem, element, node, force, dis);

  optim.optstep = 0;

  // optimization start
  while(convergence == 1)
  {
    optim.optstep++;

    std::cout << "============== optimization step ==============" << endl;
    std::cout << "step no. " << optim.optstep << endl;

    // filtering procedure
    int beta_tmp = max((optim.optstep - 100) / 40.0, 1.0);
    double beta = min((double)beta_tmp, 8.0);
    double T = 0.5;
    std::cout << "beta: " << beta << endl;

    auto design_filtered = filter.densityFilter(optim.design_s, false);
    auto design_threshold = filter.threshold(design_filtered, beta, T);
    auto threDifferent = filter.thresholdDerivative(design_filtered, beta, T);

    for(int nel = 0; nel < fem.nelx; nel++)
      element[nel].setDesignS(design_threshold[nel]);

    compt.solve("Eigen_CG", 1);

    /// get  stress
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

    // get objective function
    optim.object_f = dis.val.dot(force.fext);

    // get const_func
    optim.const_h[0] = 0.0;
    for(int i = 0; i < fem.nelx; i++)
      optim.const_h[0] += element[i].design_s() * element[i].volume;
    optim.const_h[0] -= V0;

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);

    std::cout << "----------------- gcmma start -----------------" << endl;
    std::cout << "function value : " << optim.object_f << endl;
    std::cout << "constraint value : " << optim.const_h[0] << endl;

    // export vtk
    output(fem, element, node, optim);

    // convergence check
    convergence = optim.isConvergence(1.0e-5, 50);

    // analytical
    sensitivity(fem, element, node, optim, threDifferent);

    optim.dfds = filter.densityFilter(optim.dfds, true);
    optim.dgds = filter.densityFilter(optim.dgds, true);

    /////////////////////////////////////////////////////////
    ////////////////////// GCMMA start //////////////////////
    /////////////////////////////////////////////////////////
    // temporary design variable for GCMMA scheme
    vector<double> design_s_new(fem.nelx);
    vector<double> const_h_new(numConst);
    double object_f_new;

    gcmma.OuterUpdate(design_s_new.data(), optim.design_s.data(),
                      optim.object_f, optim.dfds.data(), optim.const_h.data(),
                      optim.dgds.data(), optim.xmin.data(), optim.xmax.data());

    design_filtered = filter.densityFilter(design_s_new, false);
    design_threshold = filter.threshold(design_filtered, beta, T);

    for(int nel = 0; nel < fem.nelx; nel++)
      element[nel].setDesignS(design_threshold[nel]);

    compt.solve("Eigen_CG", 0);

    // get objective function
    object_f_new = dis.val.dot(force.fext);
    std::cout << "function value : " << optim.object_f << endl;

    // get const_func
    const_h_new[0] = 0.0;
    for(int i = 0; i < fem.nelx; i++)
      const_h_new[0] += element[i].design_s() * element[i].volume;
    const_h_new[0] -= V0;
    std::cout << "constraint value : " << optim.const_h[0] << endl;

    // first gcmma check
    bool conserv = gcmma.ConCheck(object_f_new, const_h_new.data());
    if(conserv)
      std::cout << "GCMMA inner loop: 0" << endl;

    // GCMMA inner loop
    for(int inneriter = 0; !conserv && inneriter < 15; ++inneriter)
    {
      // Inner iteration update
      gcmma.InnerUpdate(design_s_new.data(), object_f_new, const_h_new.data(),
                        optim.design_s.data(), optim.object_f,
                        optim.dfds.data(), optim.const_h.data(),
                        optim.dgds.data(), optim.xmin.data(),
                        optim.xmax.data());

      design_filtered = filter.densityFilter(design_s_new, false);
      design_threshold = filter.threshold(design_filtered, beta, T);

      for(int nel = 0; nel < fem.nelx; nel++)
        element[nel].setDesignS(design_threshold[nel]);

      compt.solve("Eigen_CG", 0);

      // get objective function
      object_f_new = dis.val.dot(force.fext);

      // get const_func
      const_h_new[0] = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        const_h_new[0] += element[i].design_s() * element[i].volume;
      const_h_new[0] -= V0;

      // gcmma check (if conserv=1, it'll break out of the loop.)
      conserv = gcmma.ConCheck(object_f_new, const_h_new.data());
      std::cout << "GCMMA inner loop: " << inneriter << endl;
    }

    // finally update design variables
    optim.design_s = design_s_new;

    // not convergence
    if(optim.optstep == 1000)
    {
      cerr << "not convergence in this optimization" << endl;
      convergence = 1;
    }
  }

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  std::cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void sensitivity(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
                 Optimize &optim, vector<double> &threDiff)
{
  assert(element[0].young1() < element[0].young2());

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

    optim.dfds[nel] =
        -threDiff[nel] * E * disp.transpose() * element[nel].Ke() * disp;

    optim.dgds[nel] = threDiff[nel] * element[nel].volume;
  }

  // output sensivity graph
  if(optim.optstep == 1)
    exportGraph(optim.dfds, "objective.csv");
}

void output(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
            Optimize &optim)
{
  ofstream fout("res" + to_string(optim.optstep) + ".vtk");
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
                      [&](LinearElastic &e) { return e.design_s(); });

  addElementScalarVTK("young", element, fout, false,
                      [&](LinearElastic &e) { return e.young(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](LinearElastic &e) { return e.mises(); });

  fout.close();
}