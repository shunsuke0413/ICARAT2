#pragma once
#include <fstream>
#include <iostream>
#include <vector>

/// node class
struct Node
{
  double x[3];  ///< coordinate x
  int globalID; /// global ID
  int localNodeID;

  /// @note Type have more than 1
  ///       type=1...no sharing
  ///       type=n>1...sharing with n nodes of other process
  int type = 0;
  std::vector<std::pair<int, int>> procNode; ///< process,Node
};

/// element class
struct Element
{
  std::vector<int> nodeID;
};

/// vtk import class
///@sa https://github.com/PANFACTORY/PANSLBM2
class VTKImport
{
public:
  VTKImport() = delete;
  VTKImport(std::string _fname, std::string _ctype)
  {
    this->fin.open(_fname);
    std::string tmp;

    int cTypeInt = getCellType(_ctype);

    // 1. read cell type
    std::vector<int> CTarray;
    while(this->fin >> tmp)
    {
      if(tmp == "CELL_TYPES")
      {
        int numCtypes, tNumber;
        this->fin >> numCtypes;

        int counter = 0;
        for(int i = 0; i < numCtypes; i++)
        {
          this->fin >> tNumber;
          if(tNumber == cTypeInt)
            CTarray.push_back(counter);

          counter++;
        }

        nelx = counter;

        break;
      }
    }

    // 2. read points and cells
    this->fin.clear();
    this->fin.seekg(0, std::ios_base::beg);
    while(this->fin >> tmp)
    {
      if(tmp == "POINTS")
      {
        this->fin >> numnp >> tmp;
        node.resize(numnp);

        for(int i = 0; i < numnp; i++)
        {
          node[i].globalID = i;
          this->fin >> node[i].x[0] >> node[i].x[1] >> node[i].x[2];
        }
      }
      else if(tmp == "CELLS")
      {
        int Tnelx;
        this->fin >> Tnelx >> tmp;

        element.resize(nelx);

        int ne;
        int counter = 0;
        for(int i = 0; i < Tnelx; i++)
        {
          // If it matches the celltype
          if(i == CTarray[counter])
          {
            this->fin >> ne;
            element[counter].nodeID.resize(ne);
            for(int j = 0; j < ne; j++)
              this->fin >> element[counter].nodeID[j];

            counter++;
          }
          // If it doesn't match the celltype
          else
          {
            this->fin >> ne;
            for(int j = 0; j < ne; j++)
              this->fin >> tmp;
          }
        }
        ////////////////////
        break;
        ////////////////////
      }
    }
  }
  VTKImport(const VTKImport &) = delete;
  ~VTKImport() {}

  std::ifstream fin;
  int nelx;
  int numnp;
  std::vector<Element> element;
  std::vector<Node> node;

private:
  int getCellType(std::string name)
  {
    if(name == "TRIA3")
      return 5;
    else if(name == "TETRA4")
      return 10;
    else if(name == "HEXA8")
      return 12;
    else if(name == "QUAD4")
      return 9;
    else
    {
      std::cerr << "Element type is invalid." << std::endl;
      exit(1);
    }
  }
};
