///"g++ -std=c++14 vtk_to_in.cpp -o vtk_to_in.out"

#include "vtk_to_in.hpp"
#include <sys/stat.h>

using namespace std;

void makeSubMesh(int procNumber, string &metisFile, string &inName,
                 VTKImport &vtk);

int main(int argc, char **argv)
{
  //////////////input zone///////////////
  // 0. store vtk data
  cout << "import vtk file..." << endl;
  fflush(stdout);
  VTKImport vtk("output.vtk", "HEXA8");
  cout << "done" << endl;
  fflush(stdout);

  // process number
  int procNumber = 10;
  // command to create a folder to store in files
  string dir_in = "../sample/mpi/opt-linear/import-big";
  // .in file name
  string inName = dir_in + "/mesh";
  // metis input file name
  string metisFile = "metis.mesh";
  // run metis command
  string metis_exe_str = "./../vendor/metis/bin/mpmetis " + metisFile + " " +
                         to_string(procNumber);
  const char *metis_exe = metis_exe_str.c_str();
  ////////////////////////////////////////

  // 1. make mesh file for metis
  cout << "export inputfile for metis..." << endl;
  fflush(stdout);
  ofstream fout(metisFile);

  fout << vtk.nelx << "\n";
  for(int i = 0; i < vtk.nelx; i++)
  {
    for(int j = 0; j < (int)vtk.element[i].nodeID.size(); j++)
      fout << vtk.element[i].nodeID[j] + 1 << " ";

    if(i != vtk.nelx - 1)
      fout << "\n";
  }

  fout.close();
  cout << "done" << endl;
  fflush(stdout);

  // 2. run metis
  int res = system(metis_exe);

  // 3. make in file
  res = mkdir(dir_in.c_str(), 0777);
  makeSubMesh(procNumber, metisFile, inName, vtk);
  fflush(stdout);

  return 0;
}

void makeSubMesh(int procNumber, string &metisFile, string &inName,
                 VTKImport &vtk)
{
  // 1. read metis output file
  ifstream fin(metisFile + ".epart." + to_string(procNumber));
  vector<int> localNelx(procNumber);
  vector<vector<int>> localElementID(procNumber, vector<int>(vtk.nelx));
  int procid;
  for(int nel = 0; nel < vtk.nelx; nel++)
  {
    fin >> procid;
    localElementID[procid][localNelx[procid]++] = nel;
  }

  // 2. check shared nodes
  vector<vector<int>> sharedNodeFlag(procNumber, vector<int>(vtk.numnp));
  vector<int> localNumnp(procNumber);
  // each process loop
  for(int processID = 0; processID < procNumber; processID++)
  {
    localNumnp[processID] = 0;
    for(int i = 0; i < vtk.numnp; i++)
    {
      vtk.node[i].localNodeID = -1;
      sharedNodeFlag[processID][i] = 0;
    }

    // element loop
    for(int nel = 0; nel < localNelx[processID]; nel++)
    {
      int ne = vtk.element[localElementID[processID][nel]].nodeID.size();

      for(int k = 0; k < ne; k++)
      {
        int nIDtmp = vtk.element[localElementID[processID][nel]].nodeID[k];

        // when this node is shared
        if(vtk.node[nIDtmp].localNodeID == -1)
          vtk.node[nIDtmp].localNodeID = localNumnp[processID]++;

        // when this node is shared & meet this node for the first time
        if(vtk.node[nIDtmp].globalID > -1 &&
           sharedNodeFlag[processID][nIDtmp] == 0)
        {
          vtk.node[nIDtmp].type++;
          sharedNodeFlag[processID][nIDtmp] = 1;
          std::pair<int, int> pair;
          pair.first = processID;
          pair.second = vtk.node[nIDtmp].localNodeID;
          vtk.node[nIDtmp].procNode.push_back(pair);
        }
      } // ne loop end
    }   // element loop end
  }     // process loop end

  // 3. make in file
  cout << "make in file..." << endl;
  fflush(stdout);
  for(int procid = 0; procid < procNumber; procid++)
  {
    cout << "\r" << procid << "/" << procNumber - 1;
    fflush(stdout);

    string inFile;
    if(procid < 10)
      inFile = inName + "0" + to_string(procid) + ".in";
    else
      inFile = inName + to_string(procid) + ".in";

    ofstream fout(inFile);

    // number of node & element (local)
    fout << localNumnp[procid] << " " << localNelx[procid] << endl;

    // coordinate
    localNumnp[procid] = 0;
    for(int j = 0; j < vtk.numnp; j++)
      vtk.node[j].localNodeID = -1;

    for(int nel = 0; nel < localNelx[procid]; nel++)
    {
      int ne = vtk.element[localElementID[procid][nel]].nodeID.size();

      for(int k = 0; k < ne; k++)
      {
        int nIDtmp = vtk.element[localElementID[procid][nel]].nodeID[k];

        if(vtk.node[nIDtmp].localNodeID == -1)
        {
          vtk.node[nIDtmp].localNodeID = localNumnp[procid]++;
          fout << vtk.node[nIDtmp].localNodeID << " " << vtk.node[nIDtmp].x[0]
               << " " << vtk.node[nIDtmp].x[1] << " " << vtk.node[nIDtmp].x[2]
               << " " << vtk.node[nIDtmp].type << " ";

          for(int p = 0; p < (int)vtk.node[nIDtmp].procNode.size(); p++)
            fout << vtk.node[nIDtmp].procNode[p].first << " "
                 << vtk.node[nIDtmp].procNode[p].second << " ";

          fout << std::endl;
        }
      } // ne loop
    }   // element loop

    // connectivity
    for(int nel = 0; nel < localNelx[procid]; nel++)
    {
      int ne = vtk.element[localElementID[procid][nel]].nodeID.size();

      for(int k = 0; k < ne; k++)
      {
        int nIDtmp = vtk.element[localElementID[procid][nel]].nodeID[k];
        fout << vtk.node[nIDtmp].localNodeID << " ";
      }
      fout << std::endl;
    }

    // close .in file
    fout.close();
  }

  cout << endl << "end" << endl;
  fflush(stdout);

  // close metis output file
  fin.close();
}