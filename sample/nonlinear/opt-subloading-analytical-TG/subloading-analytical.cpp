#include "subloading-analytical.hpp"

using namespace std;
using namespace Eigen;
using namespace icarat;
using namespace icarat::opt_elast_plast;

void output(FEM &fem, vector<OptLinearElastPlast> &element, vector<Node> &node,
            Optimize &optim);

int main(int argc, char *argv[])
{
  if(argc == 1)
  {
    std::cout << "usage1: " << argv[0] << " <.toml file>\n";
    // std::cout << "usage2: " << argv[0] << " <.toml file> <.vtk file>\n";//vtk
    std::cout << "usage2: " << argv[0]
              << " <.toml file> <.mdpa file>\n"; // Icarat
    return 1;
  }

  // measure time
  auto t_start = chrono::system_clock::now();
  int convergence = 1; // convergence judge
  int numConst = 1;    // number of constraint
  double V0 = 0.0;     // first volume
  double T = 0.5;
  // set input file
  const auto config = toml::parse(argv[1]);
  int numStep = toml::find<int>(config, "this_problem", "numstep");
  double filterR = toml::find<double>(config, "this_problem", "filterR");
  std::array<int, 4> covers = {
      toml::find<int>(config, "boundary_condition", "cover-L"),
      toml::find<int>(config, "boundary_condition", "cover-R"),
      toml::find<int>(config, "boundary_condition", "cover-B"),
      toml::find<int>(config, "boundary_condition", "cover-T")};
  // string pathVTK=argv[2];//vtk
  FEM fem;
  MeshIcarat<OptLinearElastPlast, Node> mesh(fem, config); // Icarat
  // MeshVTK<OptLinearElastPlast, Node> mesh(fem, pathVTK,config);//vtk
  mesh.setParameter("structure");

  // set element(unstructured mesh version)
  vector<OptLinearElastPlast> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      OptLinearElastPlast actele(fem.voigt, mesh.ne(type), mesh.ipmax(type),
                                 numdof);
      actele.setParameter(config);
      element.push_back(actele);
    }
  }

  int totalTime = numStep;
  // auto designs=mesh.importElementScalars("density");//vtk
  //  set node
  vector<Node> node(fem.numnp, icarat::Node(totalTime));
  // input mesh to element & node
  mesh.generate(element, node);

  // boundary condition class
  BCIcarat<OptLinearElastPlast, Node> BC(fem, element, node);

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq, totalTime);
  ValueNR dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  // // set control node for load-displacement
  // int direction;
  // int control_ID;
  // std::string dire =
  //     toml::find<std::string>(config, "load_displacement", "direction");
  // std::array<double, 3> coordinate = toml::find<std::array<double, 3>>(
  //     config, "load_displacement", "observation_point");

  // if(dire == "X")
  //   direction = 0;
  // else if(dire == "Y")
  //   direction = 1;
  // else if(dire == "Z")
  //   direction = 2;

  // for(int i = 0; i < fem.numnp; i++)
  // {
  //   if(node[i].x[0] == coordinate[0] && node[i].x[1] == coordinate[1] &&
  //      node[i].x[2] == coordinate[2])
  //   {
  //     control_ID = i;
  //     break;
  //   }
  // }

  int element_ID = toml::find<int>(config, "mises_eqstn", "element_ID");
  // int element_ID_sens =
  //     toml::find<int>(config, "sensitivity_history", "element_ID_sens");

  Optimize optim(fem.nelx, numConst, 1.0e-4, 1.0);
  OC oc(fem.nelx);
  MMASolver mma(fem.nelx, numConst);
  Filter<OptLinearElastPlast, Node> filter(fem, element, node, filterR);
  double mesh_length = mesh.length()[0] / mesh.div()[0];

  // initialization
  optim.optstep = 0;
  for(int i = 0; i < fem.nelx; i++)
    optim.design_s[i] = element[i].design_s();

  for(int el = 0; el < fem.nelx; el++)
  {
    int ix = el % mesh.div()[0];
    int iy = el / mesh.div()[0];

    if(ix < covers[0] || ix >= (mesh.div()[0] - covers[1]) || iy < covers[2] ||
       iy >= (mesh.div()[1] - covers[3]))
    {
      optim.xmin[el] = 1.0 - 1.0e-4;
      optim.xmax[el] = 1.0;
      optim.design_s[el] = 1.0;
    }
  }
  // ==========================================
  // 1. メッシュの分割数を取得
  // ==========================================
  int nx = mesh.div()[0];
  int ny = mesh.div()[1];

  // ==========================================
  // 2. TOMLから非設計領域のボックス（インデックス範囲）を読み込む
  // ==========================================
  struct BoxIndex3D
  {
    int xmin, xmax;
    int ymin, ymax;
    int zmin, zmax;
  };

  std::vector<BoxIndex3D> non_design_list;
  auto boxes = toml::find<std::vector<std::vector<int>>>(config, "non_design",
                                                         "non_design_boxes");

  for(const auto &b : boxes)
  {
    BoxIndex3D box;
    box.xmin = b[0];
    box.xmax = b[1];
    box.ymin = b[2];
    box.ymax = b[3];
    box.zmin = b[4];
    box.zmax = b[5];
    non_design_list.push_back(box);
  }

  // ==========================================
  // 3. 全要素分のフラグ配列を用意し、ボックス内か判定して保存する
  // ==========================================
  std::vector<bool> is_non_design_flag(fem.nelx, false);

  for(int el = 0; el < fem.nelx; el++)
  {
    int iz = el / (nx * ny);
    int iy = (el / nx) % ny;
    int ix = el % nx;

    bool is_in_non_design = false;
    for(const auto &box : non_design_list)
    {
      if(ix >= box.xmin && ix <= box.xmax && iy >= box.ymin && iy <= box.ymax &&
         iz >= box.zmin && iz <= box.zmax)
      {
        is_in_non_design = true;
        break;
      }
    }
    is_non_design_flag[el] = is_in_non_design;
  }

  // ==========================================
  // 4. 判定結果をもとに、最適化変数（optim）の境界条件を固定する
  // ==========================================
  for(int el = 0; el < fem.nelx; el++)
  {
    if(is_non_design_flag[el])
    {
      optim.xmin[el] = 1.0e-4;
      optim.xmax[el] = 1.0e-4 + 1.0e-4;
      optim.design_s[el] = 1.0e-4;
    }
  }

  NewtonRaphson<OptLinearElastPlast, Node> compt(fem, element, node, force,
                                                 dis);

  // compute time of preprocess
  auto t_now = chrono::system_clock::now();
  auto t_dur = t_now - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  // addPlot(optim.optstep, t_sec, "computetime.csv", 1);

  // optimization start
  while(convergence == 1)
  {
    optim.optstep++;

    cout << "============== optimization step ==============" << endl;
    cout << "step no. " << optim.optstep << endl;

    auto design_filtered = filter.densityFilter2(optim.design_s, false,
                                                 optim.optstep, mesh_length);

    for(int el = 0; el < fem.nelx; el++)
    {
      int ix = el % mesh.div()[0];
      int iy = el / mesh.div()[0];

      if(ix < covers[0] || ix >= (mesh.div()[0] - covers[1]) ||
         iy < covers[2] || iy >= (mesh.div()[1] - covers[3]))
      {
        design_filtered[el] = 1.0;
      }

      if(is_non_design_flag[el])
      {
        design_filtered[el] = 1.0e-4;
      }
    }

    // /// set threshold function
    // double beta = 1.0;
    // cout << "beta: " << beta << endl;
    // //auto design_threshold = filter.threshold(design_filtered, beta, T);
    // //for optimization auto design_threshold = design_filtered; //for
    // sensitivity analysis auto threDeriv =
    // filter.thresholdDerivative(design_filtered, beta, T);

    // set threshold design variable
    for(int i = 0; i < fem.nelx; i++)
    {
      element[i].setDesignS(design_filtered[i]); // Icarat
      // element[i].setDesignS(optim.design_s[i]);//Icarat
      // element[i].setDesignS(designs[i]);//vtk
    }

    // Initialization of element internal variables
    compt.initialization();
    for(int nel = 0; nel < fem.nelx; nel++)
    {
      element[nel].initialization();
      element[nel].init_opt();
    }

    // output_str(fem, element, node, optim, 0);//

    // addPlot(node[control_ID].val[direction],
    // node[control_ID].reaction[direction],
    //         "load_displacement.csv", 1);

    addPlot(element[element_ID].eqstn(), element[element_ID].mises(),
            "mises_eqstn.csv", 1);
    addPlot(element[element_ID].eqstn_gauss(),
            element[element_ID].mises_gauss(), "mises_eqstn_gauss.csv", 1);

    for(int loadStep = 1; loadStep < numStep + 1; loadStep++)
    {
      for(auto &e : element)
        e.iteration = 0;
      force.fext = loadStep * force.fext_org / numStep;

      cout << "---------- load-step No. " << loadStep << " ----------" << endl;

      compt.solve_M2(loadStep, numStep, 1.0e-6, 10, "Eigen_LDLT", 0);

      vector<Triplet<double>> tripletsR;
      for(int nel = 0; nel < fem.nelx; nel++)
      {
        element[nel].stress_plastic(fem, node, tripletsR);
        element[nel].sensitivity(fem, node);
      }
      // output_str(fem, element, node, optim, loadStep);//

      // addPlot(node[control_ID].val[direction],
      //         node[control_ID].reaction[direction], "load_displacement.csv",
      //         0);
      addPlot(element[element_ID].eqstn(), element[element_ID].mises(),
              "mises_eqstn.csv", 0);
      addPlot(element[element_ID].eqstn_gauss(),
              element[element_ID].mises_gauss(), "mises_eqstn_gauss.csv", 0);
    }

    // get objective function
    optim.object_f = 0.0;
    for(int i = 0; i < fem.nelx; i++)
      optim.object_f += element[i].func_new();

    cout << "object function: " << optim.object_f << endl;

    // get const_func
    if(optim.optstep == 1)
    {
      V0 = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        V0 += element[i].design_s() * element[i].volume;
      if(V0 == 0.0)
        cerr << "V0 = 0. in optimization" << endl;

      // optim.const_h[0] = V0; //OC
      optim.const_h[0] = 0.0; // MMA
    }
    else
    {
      optim.const_h[0] = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        optim.const_h[0] += element[i].design_s() * element[i].volume;
      optim.const_h[0] -= V0; // MMA

      cout << "const value: " << optim.const_h[0] << ",  "
           << "first volume: " << V0 << endl;
    }

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);
    output(fem, element, node, optim);

    convergence = optim.isConvergence(1.0e-7, 50);

    // for(int i = 0; i < fem.nelx; i++)
    // {
    //   optim.dfds[i] = element[i].objective() * threDeriv[i];
    //   optim.dgds[i] = element[i].volume * threDeriv[i];
    // }
    for(int i = 0; i < fem.nelx; i++)
    {
      optim.dfds[i] = element[i].objective();
      optim.dgds[i] = element[i].volume;
    }
    // filtering sensitivity
    optim.dfds =
        filter.densityFilter2(optim.dfds, true, optim.optstep, mesh_length);
    optim.dgds =
        filter.densityFilter2(optim.dgds, true, optim.optstep, mesh_length);

    for(int el = 0; el < fem.nelx; el++)
    {
      int ix = el % mesh.div()[0];
      int iy = el / mesh.div()[0];

      if(ix < covers[0] || ix >= (mesh.div()[0] - covers[1]) ||
         iy < covers[2] || iy >= (mesh.div()[1] - covers[3]))
      {
        optim.dfds[el] = 0.0; // 目的関数の感度を強制リセット
        optim.dgds[el] = 0.0;
        // 制約関数の感度（dgdxなど）配列があればそれも 0.0 にする
      }

      if(is_non_design_flag[el])
      {
        optim.dfds[el] = 0.0;
        optim.dgds[el] = 0.0;
      }
    }

    exportGraph(optim.dfds, "objective.csv");

    // exit(1);

    // if(toml::find<bool>(config, "this_problem", "taylor_test"))
    //   taylor_test(fem, element, node, force, dis, optim, numStep, filterR);
    if(toml::find<bool>(config, "this_problem", "fdm"))
      fdm(fem, element, node, force, dis, optim, numStep, filterR);

    // Get updated design variables with MMA
    mma.Update(optim.design_s.data(), optim.dfds.data(), optim.const_h.data(),
               optim.dgds.data(), optim.xmin.data(), optim.xmax.data());

    for(int el = 0; el < fem.nelx; el++)
    {
      if(is_non_design_flag[el])
      {
        optim.design_s[el] = 1.0e-4;
      }
    }

    // not convergence
    if(optim.optstep == 1000)
    {
      cerr << "not convergence in this optimization" << endl;
      convergence = 1;
    }

    // compute time at each step
    // t_now = chrono::system_clock::now();
    // t_dur = t_now - t_start;
    // t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
    // addPlot(optim.optstep, t_sec, "computetime.csv", 0);
  }

  auto t_end = chrono::system_clock::now();
  t_dur = t_end - t_start;
  t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void output(FEM &fem, vector<OptLinearElastPlast> &element, vector<Node> &node,
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
                      [&](OptLinearElastPlast &e) { return e.design_s(); });

  addElementScalarVTK("epstn", element, fout, false,
                      [&](OptLinearElastPlast &e) { return e.eqstn(); });

  addElementScalarVTK("epsilon", element, fout, false,
                      [&](OptLinearElastPlast &e) { return e.epsilon(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](OptLinearElastPlast &e) { return e.mises(); });

  fout.close();
}
