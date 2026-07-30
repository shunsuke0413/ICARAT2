#include "opt_macro_mma_nondesign.hpp"

using namespace std;
using namespace icarat;

void sensitivity(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
                 Optimize &optim, double V0);
void output(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
            Optimize &optim);
pair<vector<int>, vector<int>>
getDesignNonDesignIDs(FEM &fem, vector<LinearElastic> &element,
                      vector<Node> &node, const toml::value &config);

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <.toml file> <.mdpa file>\n";
    return 1;
  }

  // measure time
  auto t_start = chrono::system_clock::now();
  int convergence = 1; // convergence judge
  int numConst = 1;    // number of constraint
  double V0 = 0.0;     // first volume
  double T = 0.5;
  double xmin = 1.0e-4;
  double xmax = 1.0;

  // set input file
  const auto config = toml::parse(argv[1]);
  double filterR = toml::find<double>(config, "this_problem", "filterR");

  // set problem
  FEM fem;
  MeshIcarat<LinearElastic, Node> mesh(fem, config);
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

  // get pair of design area's ID & non-design area's ID
  auto IDs = getDesignNonDesignIDs(fem, element, node, config);

  // set optimization's classes
  Optimize optim(fem.nelx, numConst, xmin, xmax);
  MMASolver mma(fem.nelx, numConst);
  Filter<LinearElastic, Node> filter(fem, element, node, filterR);

  // initialization
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

    auto design_filtered = filter.densityFilter(optim.design_s, false);

    /// set threshold function
    double beta = min(pow(2, (int)(optim.optstep / 40)), 16.0);
    cout << "beta: " << beta << endl;
    auto design_threshold = filter.threshold(design_filtered, beta, T);
    auto threDeriv = filter.thresholdDerivative(design_filtered, beta, T);

    // set threshold design variable
    for(int i = 0; i < fem.nelx; i++)
      element[i].setDesignS(design_threshold[i]);

    compt.solve("Eigen_CG", 1);

    /// get  stress
    for(auto &e : element)
      e.makeMisesStress(fem, node);

    // get objective function
    optim.object_f = dis.val.dot(force.fext);
    cout << "function value = " << optim.object_f << endl;

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

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);
    output(fem, element, node, optim);

    convergence = optim.isConvergence(1.0e-5, 50);

    // analytical
    sensitivity(fem, element, node, optim, V0);

    for(int i = 0; i < fem.nelx; i++)
    {
      optim.dfds[i] *= threDeriv[i];
      optim.dgds[i] *= threDeriv[i];
    }

    // filtering sensitivity
    optim.dfds = filter.densityFilter(optim.dfds, true);
    optim.dgds = filter.densityFilter(optim.dgds, true);

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
    if(optim.optstep == 1000)
    {
      cerr << "not convergence in this optimization" << endl;
      convergence = 1;
    }
  }

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec \n";

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

    optim.dfds[nel] = -E * disp.transpose() * element[nel].Ke() * disp;

    optim.dgds[nel] = element[nel].volume / V0;
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

  // export scalar data to element
  addElementScalarVTK("density", element, fout, true,
                      [&](LinearElastic &e) { return e.design_s(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](LinearElastic &e) { return e.mises(); });

  fout.close();
}

using Type =
    std::tuple<std::array<double, 1>, std::array<std::array<double, 2>, 3>>;
pair<vector<int>, vector<int>>
getDesignNonDesignIDs(FEM &fem, vector<LinearElastic> &element,
                      vector<Node> &node, const toml::value &config)
{
  vector<int> designIDs, nondesignIDs;

  const auto non_designs =
      toml::find<toml::array>(config, "this_problem", "non_design");

  for(auto &nonID : non_designs)
  {
    Type data = toml::get<Type>(nonID);
    double value = std::get<0>(data)[0];
    auto range = std::get<1>(data);
    for(int i = 0; i < fem.nelx; i++)
    {
      // get coordinates of cell center
      double centerX = 0.0, centerY = 0.0, centerZ = 0.0;
      for(int j = 0; j < element[i].ne; j++)
      {
        centerX += node[element[i].nodeID[j]].x[0] / element[i].ne;
        centerY += node[element[i].nodeID[j]].x[1] / element[i].ne;
        centerZ += node[element[i].nodeID[j]].x[2] / element[i].ne;
      }

      if(fem.ndim == 3)
      {
        if(range[0][0] < centerX && centerX < range[0][1] && //
           range[1][0] < centerY && centerY < range[1][1] && //
           range[2][0] < centerZ && centerZ < range[2][1])
        {
          nondesignIDs.push_back(i);
          element[i].setDesignS(value);
        }
        else
        {
          designIDs.push_back(i);
        }
      }
      else
      {
        if(range[0][0] < centerX && centerX < range[0][1] && //
           range[1][0] < centerY && centerY < range[1][1])
        {
          nondesignIDs.push_back(i);
          element[i].setDesignS(value);
        }
        else
        {
          designIDs.push_back(i);
        }
      }
    }
  }

  // 非設計領域が無い場合
  if(non_designs.size() == 0)
  {
    std::iota(designIDs.begin(), designIDs.end(), 0);
    designIDs.resize(fem.nelx);
  }

  // Delete duplicate IDs.
  std::sort(designIDs.begin(), designIDs.end());
  designIDs.erase(std::unique(designIDs.begin(), designIDs.end()),
                  designIDs.end());

  std::sort(nondesignIDs.begin(), nondesignIDs.end());
  nondesignIDs.erase(std::unique(nondesignIDs.begin(), nondesignIDs.end()),
                     nondesignIDs.end());

  return make_pair(designIDs, nondesignIDs);
}
