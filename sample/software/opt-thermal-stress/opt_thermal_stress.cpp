///@sa 市川智ら, 熱応力を考慮したトポロジー最適化の基礎的検討, 土木学会論文集,
/// Vol.71, No.2, 2015.

#include "opt_thermal_stress.hpp"
#include <misc/assembling.hpp>
using namespace icarat;
using namespace std;
using namespace Eigen;

typedef std::vector<std::vector<std::vector<double>>> ijk;
typedef std::vector<std::vector<double>> ij;
pair<ij, ijk> sensitivity(FEM &fem,
                          vector<LinearThermalStressElastMultiMat> &element,
                          vector<Node> &node, Force &force, ValueNR &value,
                          const SparseMatrix<double> &K, VectorXd &dualforce,
                          int numMat, vector<double> &VV);

void output(FEM &fem, vector<LinearThermalStressElastMultiMat> &element,
            vector<Node> &node, Optimize &optim, string project, int numMat);

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <project folder path>\n";
    return 1;
  }

  // measure time
  auto t_start = chrono::system_clock::now();
  int convergence = 1; // convergence judge
  double xmin = 1.0e-4;
  double xmax = 1.0;
  double eta = 0.5;
  double tolNR = 1.0e-6;
  int iterNR = 100;
  double V0;

  // set input file
  string project = argv[1];
  string tomlpath = project + "/TEXTfiles/tmp.toml";
  string vtkpath = project + "/VTKfiles/tmp.vtk";
  const auto config = toml::parse(tomlpath);

  // set problem
  FEM fem;

  double deltaT = toml::find<double>(config, "this_problem", "deltaT");
  int numMat = toml::find<int>(config, "material", "numMat");
  int numConst = numMat;       // number of constraint
  vector<double> VV(numConst); // initial volumes

  MeshVTK<LinearThermalStressElastMultiMat, Node> mesh(fem, vtkpath, config);
  mesh.setParameter("structure");

  // set element
  vector<LinearThermalStressElastMultiMat> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      LinearThermalStressElastMultiMat actele(fem.voigt, mesh.ne(type),
                                              mesh.ipmax(type), numdof, numMat);
      actele.setParameter(config);
      element.push_back(actele);
    }
  }

  // set node
  vector<Node> node(fem.numnp);

  // input mesh to element & node
  mesh.generate(element, node);

  // boundary condition class
  BCPoint<LinearThermalStressElastMultiMat, Node> BC(fem, element, node,
                                                     project);

  // set dirichlet condition
  BC.input_dirich(config);
  // get number of effective dofs
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  ValueNR dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  // get pair of design area's ID & non-design area's ID
  auto IDs = getDesignNonDesignIDs(fem, element, node, config, project);

  // set optimization's classes
  Optimize optim(numMat * fem.nelx, numConst, xmin, xmax);
  MMASolver mma(numMat * fem.nelx, numConst);
  double filterR = toml::find<double>(config, "this_problem", "filterR");
  FilterHelmholtz<LinearThermalStressElastMultiMat, Node> filter(fem, element,
                                                                 node, filterR);

  // optimization settings
  optim.optstep = 0;
  for(int matID = 0; matID < numMat; matID++)
    for(int i = 0; i < fem.nelx; i++)
      optim.design_s[i + matID * fem.nelx] = element[i].design_s(matID);

  auto nodeIDstr = read_file(project + "/surfaces/surface_7.csv");
  vector<int> nodeID;
  for(auto &n : nodeIDstr)
    nodeID.push_back(stoi(n));

  int direction =
      toml::find<int>(config, "this_problem", "direction"); // x=0,y=1,z=2

  VectorXd phi = VectorXd::Zero(fem.numeq);
  for(auto &n : nodeID)
    phi(node[n].dof[direction]) = 1.0;

  NewtonRaphson<LinearThermalStressElastMultiMat, Node> compt(fem, element,
                                                              node, force, dis);

  // cout << IDs.first.size() << "   " << IDs.second.size() << endl;

  // optimization start
  while(convergence == 1)
  {
    optim.optstep++;

    std::cout << "============== optimization step ==============" << endl;
    std::cout << "step no. " << optim.optstep << endl;

    // set threshold function
    double beta = min(pow(max(1, (int)(optim.optstep - 50) / 50), 2.0), 8.0);
    std::cout << "beta" << beta << std::endl;

    // filtering procedure
    vector<vector<double>> threDiffs;
    for(int matID = 0; matID < numMat; matID++)
    {
      vector<double> design_s_tmp(fem.nelx);
      for(int i = 0; i < fem.nelx; i++)
        design_s_tmp[i] = optim.design_s[i + fem.nelx * matID];
      auto design_filtered =
          filter.densityFilter(design_s_tmp, "Eigen_CG", 0, false);
      auto threDiff = filter.thresholdDerivative(design_filtered, beta, eta);
      threDiffs.push_back(threDiff);
      auto design_threshold = filter.threshold(design_filtered, beta, eta);
      for(int i = 0; i < fem.nelx; i++)
        element[i].setDesignS(matID,
                              std::clamp(design_threshold[i], xmin, xmax));
    }

    // generate thermal load
    VectorXd fth = VectorXd::Zero(fem.numeq);
    for(int i = 0; i < fem.nelx; i++)
    {
      VectorXd fthBuff = element[i].generateThermalLoad(fem, node, deltaT);
      fth += fthBuff;
    }

    // FE with thermal stress
    compt.initialization();
    force.fext = force.fext_org + fth;
    compt.solve(1, 1, tolNR, iterNR, "Eigen_LDLT", 0);

    for(auto &e : element)
      e.makeMisesStress(fem, node);

    // get const_func
    if(optim.optstep == 1)
    {
      V0 = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        V0 += element[i].volume;

      assert(V0 != 0.0);

      cout << "Volume: " << V0 << endl;
    }

    for(int matID = 0; matID < numMat; matID++)
    {
      VV[matID] = 0.0;
      for(int i = 0; i < IDs.first.size(); i++)
      {
        VV[matID] += element[IDs.first[i]].design_s(matID) *
                     element[IDs.first[i]].volume;
      }
      std::cout << "Percentage of material " << to_string(matID)
                << " present (%):" << 100.0 * VV[matID] / V0 << endl;
    }

    // get objective function
    optim.object_f = 0.0;
    for(auto &n : nodeID)
      optim.object_f += pow(node[n].val[direction], 2.0);
    std::cout << "function value " << optim.object_f << endl;

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);
    output(fem, element, node, optim, project, numMat);

    convergence = optim.isConvergence(1.0e-5, 350);

    // analytical sensitivity
    VectorXd dualforce = -2.0 * dis.val.array() * phi.array();
    auto sens = sensitivity(fem, element, node, force, dis, compt.K(),
                            dualforce, numMat, VV);

    for(int matID = 0; matID < numMat; matID++)
    {
      // threshold derivative
      for(int i = 0; i < fem.nelx; i++)
      {
        sens.first[matID][i] *= threDiffs[matID][i];

        for(int constID = 0; constID < numMat; constID++)
          sens.second[matID][constID][i] *= threDiffs[matID][i];
      }

      // filter derivative
      sens.first[matID] =
          filter.densityFilter(sens.first[matID], "Eigen_CG", 0, true);
      for(int constID = 0; constID < numMat; constID++)
        sens.second[matID][constID] = filter.densityFilter(
            sens.second[matID][constID], "Eigen_CG", 0, true);

      // output sensivity graph
      if(optim.optstep == 1)
        exportGraph(sens.first[matID], "objective" + to_string(matID) + ".csv");

      // set for MMA
      for(int i = 0; i < fem.nelx; i++)
      {
        optim.dfds[matID * fem.nelx + i] = sens.first[matID][i];
        for(int constID = 0; constID < numMat; constID++)
          optim.dgds[numMat * fem.nelx * matID + numMat * i + constID] =
              sens.second[matID][constID][i];
      }
    }

    // In the non-design domain, the sensitivity is 0
    for(int matID = 0; matID < numMat; matID++)
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
      cerr << "not convergence in this optimization" << endl;
      convergence = 0;
    }

  } // optimization loop

  auto t_end = chrono::system_clock::now();
  auto t_dur = t_end - t_start;
  auto t_sec = chrono::duration_cast<chrono::seconds>(t_dur).count();
  std::cout << "Computation time: " << t_sec << " sec \n";

  return 0;
}

pair<ij, ijk> sensitivity(FEM &fem,
                          vector<LinearThermalStressElastMultiMat> &element,
                          vector<Node> &node, Force &force, ValueNR &value,
                          const SparseMatrix<double> &K, VectorXd &dualforce,
                          int numMat, vector<double> &VV)
{
  ij dfds(numMat);  //[matID][elementID]
  ijk dgds(numMat); //[matID][constID][elementID]

  std::cout << "---- computing adjoint vector ----" << endl;
  VectorXd lambda = VectorXd::Zero(fem.numeq);
  Solvers solver(K, dualforce, lambda);
  solver.solve("Eigen_LDLT", 0);
  std::cout << "------- finish to compute -------" << endl;

  vector<VectorXd> lambdaEle =
      numeqToElement<LinearThermalStressElastMultiMat, Node>(fem, element, node,
                                                             lambda);

  for(int matID = 0; matID < numMat; matID++)
  {
    // allocation
    dfds[matID].resize(fem.nelx);
    dgds[matID].resize(numMat);
    for(int constID = 0; constID < numMat; constID++)
      dgds[matID][constID].resize(fem.nelx);

    for(int nel = 0; nel < fem.nelx; nel++)
    {
      // Young derivative
      double dYoung = 0.0;
      for(int matID1 = 0; matID1 < numMat; matID1++)
        dYoung += element[nel].dwds(matID1, matID, element[nel].pp()) *
                  element[nel].eachYoung(matID1);
      dYoung /= element[nel].young();

      // thermal stress derivative
      double dBeta = 0.0;
      for(int matID1 = 0; matID1 < numMat; matID1++)
        dBeta += element[nel].dwds(matID1, matID, 0.0) *
                 element[nel].alpha(matID1) * element[nel].eachYoung(matID1);
      dBeta /= element[nel].beta();

      // Allocate final displacement and internal force to each element
      VectorXd valEle = VectorXd::Zero(element[nel].numdof);
      for(int i = 0; i < element[nel].ne; i++)
        for(int j = 0; j < fem.ndim; j++)
          valEle.coeffRef(fem.ndim * i + j) =
              node[element[nel].nodeID[i]].val[j];

      VectorXd dfthds = dBeta * element[nel].fthi();
      dfds[matID][nel] = lambdaEle[nel].transpose() *
                         (dYoung * element[nel].Ke() * valEle - dfthds);

      for(int constID = 0; constID < numMat; constID++)
      {
        if(constID == matID)
          dgds[matID][constID][nel] = element[nel].volume / VV[matID];
        else
          dgds[matID][constID][nel] = 0.0;
      }
    }
  }

  pair<ij, ijk> pair;
  pair.first = dfds;
  pair.second = dgds;

  return pair;
}

void output(FEM &fem, vector<LinearThermalStressElastMultiMat> &element,
            vector<Node> &node, Optimize &optim, string project, int numMat)
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
  addPointVectorVTK("displacement", node, fout, true,
                    [&](Node &n) { return n.val; });

  // export scalar data to element
  addElementScalarVTK("young", element, fout, true,
                      [&](LinearThermalStressElastMultiMat &e)
                      { return e.young(); });

  addElementScalarVTK("binary", element, fout, false,
                      [&](LinearThermalStressElastMultiMat &e)
                      {
                        double tmp = 0.0;
                        if(e.design_s(1) > e.design_s(0))
                          tmp = 1.0;
                        if(!e.isDesignable)
                          tmp = 2.0;
                        return tmp;
                      });

  for(int matID = 0; matID < numMat; matID++)
    addElementScalarVTK("density" + to_string(matID), element, fout, false,
                        [&](LinearThermalStressElastMultiMat &e)
                        { return e.design_s(matID); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](LinearThermalStressElastMultiMat &e)
                      { return e.mises(); });

  addElementScalarVTK("designable", element, fout, false,
                      [&](LinearThermalStressElastMultiMat &e)
                      { return e.isDesignable; });

  fout.close();
}
