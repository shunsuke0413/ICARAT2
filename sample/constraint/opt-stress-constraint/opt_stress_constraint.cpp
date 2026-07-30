#include "opt_stress_constraint.hpp"

using namespace std;
using namespace Eigen;
using namespace icarat;
using namespace icarat::stress;

void output(FEMStressConst &fem, std::vector<LinearElastic> &element,
            std::vector<Node> &node, Optimize &optim);

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <.toml file>\n";
    return 1;
  }

  auto t_start = chrono::system_clock::now(); // measure time
  int convergence = 1;                        // convergence judge
  int numConst = 1;                           // number of constraint
  double V0 = 0.0;                            // first volume

  // set input file
  const auto config = toml::parse(argv[1]);

  // set problem
  FEMStressConst fem;
  fem.stressLim = toml::find<double>(config, "this_problem",
                                     "stress_limit"); // stress limit
  fem.pNorm =
      toml::find<double>(config, "this_problem", "p_norm"); // p-norm parameter
  fem.qp = toml::find<double>(config, "this_problem",
                              "qp"); // parameter for qp-relaxation

  MeshIcarat<LinearElastic, Node> mesh(fem, config);
  mesh.setParameter("structure");

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

  Optimize optim(fem.nelx, numConst, 1.0e-4, 1.0);
  GCMMASolver gcmma(fem.nelx, numConst);
  gcmma.SetAsymptotes(0.1, 0.7, 1.2);

  double filterR = 2.0 * mesh.length()[0] / mesh.div()[0];
  Filter<LinearElastic, Node> filter(fem, element, node, filterR);

  // init
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

    // get filtered density
    auto design_filtered = filter.densityFilter(optim.design_s, false);
    double beta = min(pow(2.0, (int)(optim.optstep / 40)), 8.0);
    double T = 0.5;
    auto design_threshold = filter.threshold(design_filtered, beta, T);
    auto threDifferent = filter.thresholdDerivative(design_filtered, beta, T);

    // set threshold design variable
    for(int i = 0; i < fem.nelx; i++)
      element[i].setDesignS(design_threshold[i]);

    compt.solve("Eigen_CG", 0);

    // get stress
    for(auto &e : element)
      e.makeMisesStress(fem, node);

    // get volume const
    double V = 0.0;
    if(optim.optstep == 1)
    {
      V0 = 0.0;
      for(int i = 0; i < fem.nelx; i++)
      {
        V0 += element[i].volume;
        V += element[i].design_s() * element[i].volume;
      }
    }
    else
    {
      for(int i = 0; i < fem.nelx; i++)
        V += element[i].design_s() * element[i].volume;
    }
    optim.object_f = V / V0;
    cout << "function value: " << optim.object_f << endl;

    // get stress const
    optim.const_h[0] = calConstraintValue(fem, element);
    cout << "stress constraint value:" << optim.const_h[0] << endl;

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);
    output(fem, element, node, optim);

    convergence = optim.isConvergence(1.0e-5, 250);

    //  analytical sensitivity
    optim.dgds = sensStressConst(fem, element, node, compt.K());

    // add threshold derivative
    for(int i = 0; i < fem.nelx; i++)
    {
      optim.dfds[i] = element[i].volume * threDifferent[i] / V0;
      optim.dgds[i] *= threDifferent[i];
    }

    // add filter's derivative
    optim.dfds = filter.densityFilter(optim.dfds, true);
    optim.dgds = filter.densityFilter(optim.dgds, true);

    // output stress constraint's sensivity graph
    if(optim.optstep == 1)
      exportGraph(optim.dgds, "constraint_sens.csv");

    // FDM
    if(toml::find<bool>(config, "this_problem", "fdm") == true)
      fdm(fem, element, node, force, dis, optim, filter);

    /////////////////////////////////////////////////////////
    ////////////////////// GCMMA start //////////////////////
    /////////////////////////////////////////////////////////
    // temporary design variable for GCMMA scheme
    vector<double> design_s_new(fem.nelx);
    vector<double> const_h_new(numConst);
    double object_f_new;

    // get new desgn variable
    gcmma.OuterUpdate(design_s_new.data(), optim.design_s.data(),
                      optim.object_f, optim.dfds.data(), optim.const_h.data(),
                      optim.dgds.data(), optim.xmin.data(), optim.xmax.data());

    // recompute objective function
    design_filtered = filter.densityFilter(design_s_new, false);
    design_threshold = filter.threshold(design_filtered, beta, T);
    for(int nel = 0; nel < fem.nelx; nel++)
      element[nel].setDesignS(design_threshold[nel]);

    compt.solve("Eigen_CG", 0);
    for(int nel = 0; nel < fem.nelx; nel++)
      element[nel].makeMisesStress(fem, node);

    V = 0.0;
    for(int i = 0; i < fem.nelx; i++)
      V += element[i].design_s() * element[i].volume;
    object_f_new = V / V0;
    std::cout << "function value : " << object_f_new << endl;

    // get const_func
    const_h_new[0] = calConstraintValue(fem, element);
    std::cout << "constraint value : " << const_h_new[0] << endl;

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
      for(int nel = 0; nel < fem.nelx; nel++)
        element[nel].makeMisesStress(fem, node);

      // get objective function
      V = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        V += element[i].design_s() * element[i].volume;
      object_f_new = V / V0;
      cout << "function value: " << object_f_new << endl;

      // get const_func
      const_h_new[0] = calConstraintValue(fem, element);
      std::cout << "constraint value : " << const_h_new[0] << endl;

      // gcmma check (if conserv=1, it'll break out of the loop.)
      conserv = gcmma.ConCheck(object_f_new, const_h_new.data());
      std::cout << "GCMMA inner loop: " << inneriter << endl;
    }

    // finally update design variables
    optim.design_s = design_s_new;

    // not convergence
    if(optim.optstep == 500)
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

void output(FEMStressConst &fem, vector<LinearElastic> &element,
            vector<Node> &node, Optimize &optim)
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

  // export scalar data to element
  addElementScalarVTK("density", element, fout, true,
                      [&](LinearElastic &e) { return e.design_s(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](LinearElastic &e) { return e.mises(); });

  fout.close();
}
