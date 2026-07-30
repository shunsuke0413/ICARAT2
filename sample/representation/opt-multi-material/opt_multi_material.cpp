#include "opt_multi_material.hpp"
using namespace icarat;
using namespace icarat::opt_multi_material;
using namespace std;
using namespace Eigen;

pair<ij, ijk> sensitivity(int numMat, FEM &fem,
                          vector<LinearElastMultiMat> &element,
                          vector<Node> &node, vector<double> V0,
                          Filter<LinearElastMultiMat, Node> &filter);
void output(FEM &fem, vector<LinearElastMultiMat> &element, vector<Node> &node,
            Optimize &optim, int numMat);

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
  // number of materials
  int numMat = toml::find<int>(config, "material", "num_mat");
  int numConst = numMat;       // number of constraint
  vector<double> V0(numConst); // initial volumes
  string optScheme = toml::find<string>(config, "material", "func_type");

  // set problem
  FEM fem;
  MeshIcarat<LinearElastMultiMat, Node> mesh(fem, config);
  mesh.setParameter("structure");

  // set element
  vector<LinearElastMultiMat> element;
  for(int type = 0; type < mesh.numeleTypes(); type++)
  {
    int numdof = mesh.ne(type) * fem.dofnp;
    for(int nel = 0; nel < mesh.nelx(type); nel++)
    {
      LinearElastMultiMat actele(fem.voigt, mesh.ne(type), mesh.ipmax(type),
                                 numdof, numMat, optScheme);
      actele.setParameter(config);
      element.push_back(actele);
    }
  }

  // set node
  vector<Node> node(fem.numnp);

  // input mesh to element & node
  mesh.generate(element, node);

  // boundary condition class
  BCIcarat<LinearElastMultiMat, Node> BC(fem, element, node);
  BC.input_dirich(config);

  // set dirichlet condition
  fem.numeq = BC.makeDOF();

  // set force & displacement
  Force force(fem.numeq);
  Value dis(fem.numeq);

  // set equivalent node load to Force class
  BC.input_neumann(config, force);

  // set optimization's classes
  Optimize optim(numMat * fem.nelx, numConst, 1.0e-4, 1.0);
  MMASolver mma(numMat * fem.nelx, numConst);
  // filter radius
  double filterR = toml::find<double>(config, "this_problem", "filterR");
  Filter<LinearElastMultiMat, Node> filter(fem, element, node, filterR);

  // initial design field
  optim.optstep = 0;
  for(int matID = 0; matID < numMat; matID++)
    for(int i = 0; i < fem.nelx; i++)
      optim.design_s[i + matID * fem.nelx] = element[i].design_ss(matID);

  LinearAnalysis<LinearElastMultiMat, Node> compt(fem, element, node, force,
                                                  dis);

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
      auto design_filtered = filter.densityFilter(design_s_tmp, false);
      auto threDiff = filter.thresholdDerivative(design_filtered, beta, 0.5);
      threDiffs.push_back(threDiff);
      // set threshold function
      auto design_threshold = filter.threshold(design_filtered, beta, 0.5);
      // set threshold design variable
      for(int i = 0; i < fem.nelx; i++)
        element[i].setDesignSS(matID, design_threshold[i]);
    }

    // structural analysis
    compt.solve("Eigen_CG", 1);

    /// get  stress for output
    for(auto &e : element)
      e.makeMisesStress(fem, node);

    // get objective function
    optim.object_f = dis.val.dot(force.fext);
    cout << "function value = " << optim.object_f << endl;

    // get const_func
    if(optim.optstep == 1)
    {
      for(int matID = 0; matID < numMat; matID++)
      {
        V0[matID] = 0.0;
        // DMO
        if(optScheme == "DMO")
          for(int i = 0; i < fem.nelx; i++)
            V0[matID] += element[i].design_ss(matID) * element[i].volume;
        // extended SIMP
        else
          for(int i = 0; i < fem.nelx; i++)
          {
            V0[matID] +=
                element[i].makeCoeff(1.0)[matID + 1] * element[i].volume;
          }

        assert(V0[matID] != 0.0);

        optim.const_h[matID] = 0.0;
        cout << "Constraint value" + to_string(matID) << ": "
             << optim.const_h[matID] << endl;
      }
    }
    else
    {
      for(int matID = 0; matID < numMat; matID++)
      {
        optim.const_h[matID] = 0.0;
        // DMO
        if(optScheme == "DMO")
          for(int i = 0; i < fem.nelx; i++)
            optim.const_h[matID] +=
                element[i].design_ss(matID) * element[i].volume / V0[matID];
        // extended SIMP
        else
          for(int i = 0; i < fem.nelx; i++)
            optim.const_h[matID] += element[i].makeCoeff(1.0)[matID + 1] *
                                    element[i].volume / V0[matID];

        optim.const_h[matID] -= 1.0;
        cout << "Constraint value" + to_string(matID) << ": "
             << optim.const_h[matID] << endl;
      }
    }

    // output result
    addPlot(optim.optstep, optim.object_f, "function.csv", optim.optstep);

    output(fem, element, node, optim, numMat);

    convergence = optim.isConvergence(1.0e-5, 50);

    // analytical
    auto tmp = sensitivity(numMat, fem, element, node, V0, filter);

    for(int matID = 0; matID < numMat; matID++)
    {
      // threshold derivative
      for(int i = 0; i < fem.nelx; i++)
      {
        tmp.first[matID][i] *= threDiffs[matID][i];

        for(int constID = 0; constID < numMat; constID++)
          tmp.second[matID][constID][i] *= threDiffs[matID][i];
      }

      // filter derivative
      tmp.first[matID] = filter.densityFilter(tmp.first[matID], true);
      for(int constID = 0; constID < numMat; constID++)
        tmp.second[matID][constID] =
            filter.densityFilter(tmp.second[matID][constID], true);

      // output sensivity graph
      if(optim.optstep == 1)
        exportGraph(tmp.first[matID], "objective" + to_string(matID) + ".csv");

      // set for MMA
      for(int i = 0; i < fem.nelx; i++)
      {
        optim.dfds[matID * fem.nelx + i] = tmp.first[matID][i];
        for(int constID = 0; constID < numMat; constID++)
          optim.dgds[numMat * fem.nelx * matID + numMat * i + constID] =
              tmp.second[matID][constID][i];
      }
    }

    // finite difference method
    if(toml::find<bool>(config, "this_problem", "fdm"))
      fdm(fem, element, node, force, dis, optim, filter, numMat);

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

pair<ij, ijk> sensitivity(int numMat, FEM &fem,
                          vector<LinearElastMultiMat> &element,
                          vector<Node> &node, vector<double> V0,
                          Filter<LinearElastMultiMat, Node> &filter)
{
  ij dfds(numMat);  //[matID][elementID]
  ijk dgds(numMat); //[matID][constID][elementID]
  for(int matID = 0; matID < numMat; matID++)
  {
    vector<double> dfds_tmp;
    ij dgds_tmp(numMat);
    for(int nel = 0; nel < fem.nelx; nel++)
    {
      // get young derivative
      double E = element[nel].dYoung(matID) / element[nel].young();

      Eigen::VectorXd disp = Eigen::VectorXd::Zero(element[nel].numdof);
      for(int i = 0; i < element[nel].ne; i++)
        for(int j = 0; j < fem.ndim; j++)
          disp.coeffRef(fem.ndim * i + j) = node[element[nel].nodeID[i]].val[j];

      dfds_tmp.push_back(-E * disp.transpose() * element[nel].Ke() * disp);

      for(int constID = 0; constID < numMat; constID++)
      {
        // DMO
        if(element[nel].type() == "DMO")
        {
          if(constID == matID)
            dgds_tmp[constID].push_back(element[nel].volume / V0[matID]);
          else
            dgds_tmp[constID].push_back(0.0);
        }
        // extended SIMP
        else
        {
          dgds_tmp[constID].push_back(
              element[nel].dwds(constID + 1, matID, 1.0) * element[nel].volume /
              V0[matID]);
        }
      }
    }

    dfds[matID] = dfds_tmp;
    dgds[matID] = dgds_tmp;
  }

  pair<ij, ijk> pair;
  pair.first = dfds;
  pair.second = dgds;

  return pair;
}

void output(FEM &fem, vector<LinearElastMultiMat> &element, vector<Node> &node,
            Optimize &optim, int numMat)
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
                      [&](LinearElastMultiMat &e) { return e.young(); });

  addElementScalarVTK("mises", element, fout, false,
                      [&](LinearElastMultiMat &e) { return e.mises(); });

  addElementScalarVTK("binary", element, fout, false,
                      [&](LinearElastMultiMat &e)
                      {
                        double max = 0.0;
                        int maxID;
                        for(int i = 0; i < e.numMat(); i++)
                        {
                          if(e.design_ss(i) > max)
                          {
                            max = e.design_ss(i);
                            maxID = i;
                          }
                        }
                        if(max < 0.5)
                          return -1;
                        else
                          return maxID;
                      });

  for(int i = 0; i < numMat; i++)
  {
    addElementScalarVTK("volume_frac" + to_string(i + 1), element, fout, false,
                        [&](LinearElastMultiMat &e)
                        { return e.makeCoeff(1.0)[i + 1]; });
  }

  fout.close();
}
