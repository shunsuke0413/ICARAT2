#include "struct-subloading-analytical.hpp"

using namespace std;
using namespace Eigen;
using namespace icarat;
using namespace icarat::opt_elast_plast;

void output(FEM &fem, vector<OptLinearElastPlast> &element, vector<Node> &node,
            Optimize &optim);

void output_str(FEM &fem, vector<OptLinearElastPlast> &element,
                vector<Node> &node, Optimize &optim, int loadstep);

int main(int argc, char *argv[])
{
  if(argc < 3)
  {
    std::cout << "usage2: " << argv[0] << " <.toml file> <.vtk file>\n"; // vtk
    return 1;
  }

  // measure time
  auto t_start = chrono::system_clock::now();
  int numConst = 1; // number of constraint
  double T = 0.5;
  // set input file
  const auto config = toml::parse(argv[1]);
  int numStep = toml::find<int>(config, "this_problem", "numstep");
  string pathVTK = argv[2]; // vtk
  FEM fem;

  // set control node for load-displacement
  int direction;
  int control_ID;
  std::string dire =
      toml::find<std::string>(config, "load_displacement", "direction");
  std::array<double, 3> coordinate = toml::find<std::array<double, 3>>(
      config, "load_displacement", "observation_point");

  if(dire == "X")
    direction = 0;
  else if(dire == "Y")
    direction = 1;
  else if(dire == "Z")
    direction = 2;

  Optimize optim(fem.nelx, numConst, 1.0e-4, 1.0);

  int element_ID = toml::find<int>(config, "mises_eqstn", "element_ID");

  optim.optstep = 1;

  for(int filenum = 51; filenum <= 51;
      filenum++) // ループ数変えたらここも変える!!!!!!
  {
    ostringstream fileNameStream;
    fileNameStream << "res" << filenum << ".vtk";
    string fileName = fileNameStream.str();

    MeshVTK<OptLinearElastPlast, Node> mesh(fem, fileName, config); // vtk
    mesh.setParameter("structure");

    // set element(unstructured mesh version)
    vector<OptLinearElastPlast> element;
    element.clear();
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
    auto designs = mesh.importElementScalars("density"); // vtk
    // set threshold design variable
    for(int i = 0; i < fem.nelx; i++)
    {
      element[i].setDesignS(designs[i]); // vtk
    }

    int totalTime = numStep;
    // set node
    vector<Node> node(fem.numnp, totalTime);
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

    for(int i = 0; i < fem.numnp; i++)
    {
      if(node[i].x[0] == coordinate[0] && node[i].x[1] == coordinate[1] &&
         node[i].x[2] == coordinate[2])
      {
        control_ID = i;
        break;
      }
    }

    // optim.design_s = designs;

    NewtonRaphson<OptLinearElastPlast, Node> compt(fem, element, node, force,
                                                   dis);

    // Initialization of element internal variables
    compt.initialization();
    for(int nel = 0; nel < fem.nelx; nel++)
    {
      element[nel].initialization();
      element[nel].init_opt();
    }

    // --- TOMLの読み込み部分 ---
    // ターゲットのリストを配列として取得
    auto boundary_targets =
        toml::find<std::vector<std::map<std::string, double>>>(
            config, "reaction_boundary", "targets");

    // --- 荷重-変位プロットの計算部分（mainのループ内） ---
    double total_R = 0.0;

    // TOMLで指定されたすべての固定端・支点をループで回して足し合わせる
    for(const auto &target : boundary_targets)
    {
      // デフォルトは -1.0（チェックしない）にしておく
      double tx = -1.0;
      double ty = -1.0;
      double tz = -1.0;

      // mapの中に "X", "Y", "Z" キーがあれば、その値を代入
      if(target.count("X"))
        tx = target.at("X");
      if(target.count("Y"))
        ty = target.at("Y");
      if(target.count("Z"))
        tz = target.at("Z");

      // 1箇所（または1線・1面）ずつの反力を計算して、全体の合計に累積する
      total_R += compt.calculateTotalReaction(direction, tx, ty, tz);
    }

    // output_str(fem, element, node, optim, 0);
    addPlot(node[control_ID].val[direction], total_R, "load_displacement.csv",
            1);
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

      compt.solve_M2(loadStep, numStep, 1.0e-6, 10, "LDLT", 0);

      vector<Triplet<double>> tripletsR;
      for(int nel = 0; nel < fem.nelx; nel++)
      {
        element[nel].stress_plastic(fem, node, tripletsR);
        element[nel].sensitivity(fem, node);
      }

      total_R = 0.0;
      // TOMLで指定されたすべての固定端・支点をループで回して足し合わせる
      for(const auto &target : boundary_targets)
      {
        // デフォルトは -1.0（チェックしない）にしておく
        double tx = -1.0;
        double ty = -1.0;
        double tz = -1.0;

        // mapの中に "X", "Y", "Z" キーがあれば、その値を代入
        if(target.count("X"))
          tx = target.at("X");
        if(target.count("Y"))
          ty = target.at("Y");
        if(target.count("Z"))
          tz = target.at("Z");

        // 1箇所（または1線・1面）ずつの反力を計算して、全体の合計に累積する
        total_R += compt.calculateTotalReaction(direction, tx, ty, tz);
      }

      // output_str(fem, element, node, optim, loadStep);
      addPlot(node[control_ID].val[direction], total_R, "load_displacement.csv",
              0);
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

    // output result
    addPlot(filenum, optim.object_f, "function.csv", filenum);
    output(fem, element, node, optim);
    optim.optstep++;
  }

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

void output(FEM &fem, vector<OptLinearElastPlast> &element, vector<Node> &node,
            Optimize &optim)
{
  ofstream fout("str" + to_string(optim.optstep) + ".vtk");
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

void output_str(FEM &fem, vector<OptLinearElastPlast> &element,
                vector<Node> &node, Optimize &optim, int loadstep)
{
  ofstream fout("Str" + to_string(loadstep) + ".vtk");
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
