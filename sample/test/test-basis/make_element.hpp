#pragma once
#include <problem/base.hpp>
#include <vector>

namespace icarat
{
namespace test
{
// make 1 element mesh with quad4 for structure analysis
void make1quad4(FEM &fem, std::vector<Element> &element,
                std::vector<Node> &node);

// make 1 element mesh with quad4 for structure analysis
void make1quad8(FEM &fem, std::vector<Element> &element,
                std::vector<Node> &node);

// make 1 element mesh with tria3 for structure analysis
void make1tria3(FEM &fem, std::vector<Element> &element,
                std::vector<Node> &node);

// make 1 element mesh with tria6 for structure analysis
void make1tria6(FEM &fem, std::vector<Element> &element,
                std::vector<Node> &node);

// make 1 element mesh with cubic8 for structure analysis
void make1cubic8(FEM &fem, std::vector<Element> &element,
                 std::vector<Node> &node);

// make 1 element mesh with cubic20 for structure analysis
void make1cubic20(FEM &fem, std::vector<Element> &element,
                  std::vector<Node> &node);

// make 1 element mesh with cubic8 for structure analysis
void make1cubic8(FEM &fem, std::vector<Element> &element,
                 std::vector<Node> &node);

// make 1 element mesh with tetra4 for structure analysis
void make1tetra4(FEM &fem, std::vector<Element> &element,
                 std::vector<Node> &node);

// make 1 element mesh with tetra10 for structure analysis
void make1tetra10(FEM &fem, std::vector<Element> &element,
                  std::vector<Node> &node);

// make 1 element mesh with prism6 for structure analysis
void make1prism6(FEM &fem, std::vector<Element> &element,
                 std::vector<Node> &node);

// make 1 element mesh with prism15 for structure analysis
void make1prism15(FEM &fem, std::vector<Element> &element,
                  std::vector<Node> &node);
} // namespace test
} // namespace icarat