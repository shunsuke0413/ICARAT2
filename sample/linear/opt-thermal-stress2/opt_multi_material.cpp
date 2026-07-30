#include "opt_multi_material.hpp"
using namespace icarat;
using namespace icarat::opt_multi_material;
using namespace std;
using namespace Eigen;

std::pair<ij, ijk>
sensitivity(FEM &fem, vector<LinearThermalStressElastMultiMat> &element,
            vector<Node> &node, Force &force, ValueNR &value,
            const SparseMatrix<double> &K, VectorXd &dualforce, int numMat,
            vector<double> &VV);

void output(FEM &fem, vector<LinearThermalStressElastMultiMat> &element,
            vector<Node> &node, Optimize &optim, int numMat);

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "usage: " << argv[0] << " <.toml file>\n";
    return 1;
  }

  // measure time
  auto t_start = chrono::system_clock::now();

  // set input file
  const auto config = toml::parse(argv[1]);

  int convergence = 1; // convergence judge
  double deltaT = toml::find<double>(config, "this_problem", "deltaT");
  int numMat = toml::find<int>(config, "material", "num_mat");
  int numConst = numMat;                     // number of constraint
  vector<double> V0(numConst), VV(numConst); //  volumes
  double fullV;                              // volume of structural model
  double tolNR = 1.0e-6;
  int iterNR = 100;

  // set problem
  FEM fem;
  MeshIcarat<LinearThermalStressElastMultiMat, Node> mesh(fem, config);
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
  BCIcarat<LinearThermalStressElastMultiMat, Node> BC(fem, element, node);
  BC.input_dirich(config);

  // set dirichlet condition
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  ValueNR dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  // set optimization's classes
  Optimize optim(numMat * fem.nelx, numConst, 1.0e-4, 1.0);
  MMASolver mma(numMat * fem.nelx, numConst);
  // filter radius
  double filterR = toml::find<double>(config, "this_problem", "filterR");
  FilterHelmholtz<LinearThermalStressElastMultiMat, Node> filter(fem, element,
                                                                 node, filterR);

  // initial design field
  optim.optstep = 0;
  for(int matID = 0; matID < numMat; matID++)
    for(int i = 0; i < fem.nelx; i++)
      optim.design_s[i + matID * fem.nelx] = element[i].design_s(matID);

  vector<int> nodeID = {59, 60, 61};
  int direction = 1;
  VectorXd phi = VectorXd::Zero(fem.numeq);
  for(auto &n : nodeID)
    phi(node[n].dof[direction]) = 1.0;

  NewtonRaphson<LinearThermalStressElastMultiMat, Node> compt(fem, element,
                                                              node, force, dis);

  output(fem, element, node, optim, numMat);

  // optimization start
  while(convergence == 1)
  {
    optim.optstep++;

    cout << "============== optimization step ==============" << endl;
    cout << "step no. " << optim.optstep << endl;

    // filtering procedure
    double beta = min(pow(2, (int)(optim.optstep / 40)), 16.0);
    cout << "beta" << beta << endl;
    vector<vector<double>> threDiffs;
    for(int matID = 0; matID < numMat; matID++)
    {
      vector<double> design_s_tmp(fem.nelx);
      for(int i = 0; i < fem.nelx; i++)
        design_s_tmp[i] = optim.design_s[i + fem.nelx * matID];
      auto design_filtered =
          filter.densityFilter(design_s_tmp, "Eigen_CG", 0, false);
      auto threDiff = filter.thresholdDerivative(design_filtered, beta, 0.5);
      threDiffs.push_back(threDiff);
      // set threshold function
      auto design_threshold = filter.threshold(design_filtered, beta, 0.5);
      // set threshold design variable
      for(int i = 0; i < fem.nelx; i++)
        element[i].setDesignS(matID, design_threshold[i]);
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

    /// get  stress for output
    for(auto &e : element)
      e.makeMisesStress(fem, node);

    // get objective function
    optim.object_f = 0.0;
    for(auto &n : nodeID)
      optim.object_f += pow(node[n].val[direction], 2.0);
    // optim.object_f = dis.val.dot(force.fext);
    std::cout << "function value = " << optim.object_f << endl;

    // get const_func
    if(optim.optstep == 1)
    {
      fullV = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        fullV += element[i].volume;

      for(int matID = 0; matID < numMat; matID++)
      {
        V0[matID] = 0.0;
        for(int i = 0; i < fem.nelx; i++)
          V0[matID] += element[i].design_s(matID) * element[i].volume;
      }
    }

    std::cout << "Full Volume " << fullV << endl;

    for(int matID = 0; matID < numMat; matID++)
    {
      VV[matID] = 0.0;
      for(int i = 0; i < fem.nelx; i++)
        VV[matID] += element[i].design_s(matID) * element[i].volume;
      std::cout << "Percentage of material " << to_string(matID)
                << " present (%):" << 100.0 * VV[matID] / fullV << endl;
    }

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);

    output(fem, element, node, optim, numMat);

    convergence = optim.isConvergence(1.0e-5, 50);

    // analytical
    VectorXd dualforce = -2.0 * dis.val.array() * phi.array();
    //    auto sens = sensitivity(numMat, fem, element, node, VV, filter);
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

    // finite difference method
    if(toml::find<bool>(config, "this_problem", "fdm"))
      fdm(fem, element, node, force, dis, optim, filter, numMat, nodeID,
          direction, deltaT);

    // Get updated design variables with MMA
    mma.Update(optim.design_s.data(), optim.dfds.data(), optim.const_h.data(),
               optim.dgds.data(), optim.xmin.data(), optim.xmax.data());

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
  cout << "Computation time: " << t_sec << " sec " << endl;

  return 0;
}

std::pair<ij, ijk>
sensitivity(FEM &fem, vector<LinearThermalStressElastMultiMat> &element,
            vector<Node> &node, Force &force, ValueNR &value,
            const SparseMatrix<double> &K, VectorXd &dualforce, int numMat,
            vector<double> &VV)
{
  ij dfds(numMat);  //[matID][elementID]
  ijk dgds(numMat); //[matID][constID][elementID]

  std::cout << "---- computing adjoint vector ----" << endl;
  VectorXd lambda = VectorXd::Zero(fem.numeq);
  Solvers solver(K, dualforce, lambda);
  solver.solve("Eigen_CG", 0);
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
                  element[nel].eachYoung(matID1) / element[nel].young();

      // thermal stress derivative
      double dBeta = 0.0;
      for(int matID1 = 0; matID1 < numMat; matID1++)
        dBeta += element[nel].dwds(matID1, matID, 0.0) *
                 element[nel].alpha(matID1) * element[nel].eachYoung(matID1);
      dBeta /= element[nel].beta();

      Eigen::VectorXd valEle = Eigen::VectorXd::Zero(element[nel].numdof);
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
            vector<Node> &node, Optimize &optim, int numMat)
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
  addElementScalarVTK("young", element, fout, true,
                      [&](LinearThermalStressElastMultiMat &e)
                      { return e.young(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](LinearThermalStressElastMultiMat &e)
                      { return e.mises(); });

  for(int i = 0; i < numMat; i++)
  {
    addElementScalarVTK("volume_frac" + to_string(i + 1), element, fout, false,
                        [&](LinearThermalStressElastMultiMat &e)
                        { return e.makeCoeff(1.0)[i + 1]; });
  }

  fout.close();
}
