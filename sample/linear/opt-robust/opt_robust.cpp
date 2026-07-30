///@sa Zhao J. et al., Robust topology optimization under loading uncertainty
/// based on linear elastic theory and orthogonal diagonalization of symmetric
/// matrices, CMAME, 2014.
///@sa Xu Y. et al., Robust topology optimization for multiple fiber-reinforced
/// plastic (FRP) composites under loading uncertainties, SMO, 2019.
#include "opt_robust.hpp"
#include "uncertainty.hpp"

using namespace std;
using namespace icarat;
using namespace icarat::uncertainty;

void output(FEM &fem, vector<LinearElastic> &element, vector<Node> &node,
            int optstep);
pair<vector<int>, vector<int>>
getDesignNonDesignIDs(FEM &fem, vector<LinearElastic> &element,
                      vector<Node> &node, const toml::value &config);

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
  string tomlname = argv[1];
  tomlname.erase(tomlname.size() - 5, 5);
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
  int numloadcase =
      toml::find<int>(config, "boundary_condition", "num_loadcase");
  int numxi = fem.ndim * numloadcase;
  Force force(fem.numeq);
  Value dis(fem.numeq);

  vector<Force> force_unc(numxi, force);
  vector<Value> dis_unc(numxi, dis);

  // set equivalent node load to Force class
  // とりあえず，全force_uncに同じforceの値を入れる．
  for(int i = 0; i < numxi; i++)
    BC.input_neumann(config, force_unc[i]);

  // ここでforce_uncのお化粧をする。（x方向とy方向それぞれのforceに分解（偶数番目がx方向、奇数番目がy方向となるように場合わけを行う））
  // ここでは，y座標に値が入っていることを前提とする．
  for(int loadcase = 0; loadcase < numloadcase; loadcase++)
  {
    for(int np = 0; np < fem.numnp; np++)
    {
      // 各節点で，荷重が値を持っているかを判定
      bool haveValue = false;
      for(int dim = 0; dim < fem.ndim; dim++)
      {
        int dof = node[np].dof[dim];
        if(dof < fem.numeq &&
           force_unc[fem.ndim * loadcase].fext_org[dof] != 0.0)
          haveValue = true;
      }

      if(haveValue) // その節点の荷重が値を持っていた場合
      {
        // 大きさを計算
        double rr = 0.0;
        for(int dim = 0; dim < fem.ndim; dim++)
        {
          int dof = node[np].dof[dim];
          if(dof < fem.numeq)
            rr += pow(force_unc[fem.ndim * loadcase].fext_org[dof], 2.0);
        }
        rr = sqrt(rr);

        for(int dim1 = 0; dim1 < fem.ndim; dim1++)
        {
          for(int dim2 = 0; dim2 < fem.ndim; dim2++)
          {
            int dof = node[np].dof[dim2];
            if(dof < fem.numeq)
            {
              if(dim1 == dim2)
                force_unc[fem.ndim * loadcase + dim1].fext_org[dof] = rr;
              else
                force_unc[fem.ndim * loadcase + dim1].fext_org[dof] = 0.0;
            }
          }
        }
      }
    }
  }

  // set uncertainty information
  int numsam = toml::find<int>(config, "uncertainty", "num_sample");
  Uncertainty unc(fem, numsam, numloadcase);
  unc.setParameter(config);

  // get pair of design area's ID & non-design area's ID
  auto IDs = getDesignNonDesignIDs(fem, element, node, config);

  // optimization class
  Optimize optim(fem.nelx, numConst, 1.0e-4, 1.0);
  MMASolver mma(fem.nelx, numConst);

  // filter radius (is 2-cell length)
  double filterR = toml::find<double>(config, "this_problem", "filterR");
  FilterHelmholtz<LinearElastic, Node> filter(fem, element, node, filterR);

  // fisrt design field
  for(int nel = 0; nel < fem.nelx; nel++)
    optim.design_s[nel] = element[nel].design_s();

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

    auto design_filtered =
        filter.densityFilter(optim.design_s, "Eigen_CG", 0, false);
    auto design_threshold = filter.threshold(design_filtered, beta, T);
    auto threDifferent = filter.thresholdDerivative(design_filtered, beta, T);

    for(int nel = 0; nel < fem.nelx; nel++)
      element[nel].setDesignS(design_threshold[nel]);

    for(int i = 0; i < numxi; i++)
    {
      LinearAnalysis<LinearElastic, Node> compt(fem, element, node,
                                                force_unc[i], dis_unc[i]);
      compt.solve("Eigen_CG", 1);

      /// get  stress
      for(auto &e : element)
        e.makeMisesStress(fem, node);
    }

    // export vtk
    output(fem, element, node, optim.optstep);
    // Calculating the expected value and variance of a random variable
    unc.make_normal_dis();

    // Calculate expected value of mean compliance and variance
    // get objective function
    unc.print_rcomp(force_unc, dis_unc);
    optim.object_f = unc.func_new;

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
      optim.const_h[0] -= 1.0;
    }

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);

    std::cout << "function value : " << optim.object_f << endl;
    std::cout << "constraint value : " << optim.const_h[0] << endl;

    // convergence check
    convergence = optim.isConvergence(1.0e-7, 50);

    // uncertain sensitivity analysis
    unc.unc_sensitivity_macro(force_unc, dis_unc, element, node, V0);

    for(int i = 0; i < fem.nelx; i++)
    {
      optim.dfds[i] = threDifferent[i] * unc.dfds[i];
      optim.dgds[i] = threDifferent[i] * unc.dgds[i];
    }

    // analytical
    optim.dfds = filter.densityFilter(optim.dfds, "Eigen_CG", 0, true);
    optim.dgds = filter.densityFilter(optim.dgds, "Eigen_CG", 0, true);

    // output sensivity graph
    if(optim.optstep == 1)
      exportGraph(optim.dfds, "objective.csv");

    // finite difference method
    if(toml::find<bool>(config, "this_problem", "fdm"))
      unc.fdm(element, node, force_unc, dis_unc, optim, filter);

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
      exit(1);
    }
  }

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  std::cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void output(FEM &fem, std::vector<LinearElastic> &element, vector<Node> &node,
            int optstep)
{
  ofstream fout("res" + to_string(optstep) + ".vtk");
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

using Type = std::tuple<std::array<double, 1>, std::array<double, 3>,
                        std::array<double, 1>, std::array<string, 2>,
                        std::array<double, 1>>;
pair<vector<int>, vector<int>>
getDesignNonDesignIDs(FEM &fem, vector<LinearElastic> &element,
                      vector<Node> &node, const toml::value &config)
{
  vector<int> designIDs, nondesignIDs;

  const auto &non_designs =
      toml::find<toml::array>(config, "this_problem", "non_design");

  if(non_designs.size() == 0)
  {
    designIDs.resize(fem.nelx);
    iota(designIDs.begin(), designIDs.end(), 0);
  }
  else
  {
    for(auto &nonID : non_designs)
    {
      Type data = toml::get<Type>(nonID);
      double value = std::get<0>(data)[0];
      auto origin = std::get<1>(data);
      double radius = std::get<2>(data)[0];
      string push = std::get<3>(data)[0];
      double thick = std::get<4>(data)[0];

      for(int i = 0; i < fem.nelx; i++)
      {
        // get coordinates of cell center
        vector<double> center(3);
        for(int j = 0; j < element[i].ne; j++)
        {
          center[0] += node[element[i].nodeID[j]].x[0] / element[i].ne;
          center[1] += node[element[i].nodeID[j]].x[1] / element[i].ne;
          center[2] += node[element[i].nodeID[j]].x[2] / element[i].ne;
        }

        if(fem.ndim == 3)
        {
          if(sqrt(pow(origin[0] - center[0], 2.0) +
                  pow(origin[1] - center[1], 2.0) +
                  pow(origin[2] - center[2], 2.0)) <= radius)
          {
            nondesignIDs.push_back(i);
            element[i].setDesignS(value);
          }
          else
          {
            designIDs.push_back(i);
          }
        }
        else // 2D
        {
          if(sqrt(pow(origin[0] - center[0], 2.0) +
                  pow(origin[1] - center[1], 2.0)) <= radius)
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
  }

  return make_pair(designIDs, nondesignIDs);
}