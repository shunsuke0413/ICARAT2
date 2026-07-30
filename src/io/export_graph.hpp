/// @file  export_graph.hpp
///  @author	Daiki Watanabe
///  @date		November 27, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <cassert>
#include <fstream>
#include <string>
namespace icarat
{
/// 2D graph ordered on the x-axis (T1 is vector data)
template <typename T> void exportGraph(T &yy, std::string filename)
{
  std::ofstream writing_file;
  writing_file.open(filename, std::ios::out);

  for(int i = 0; i < (int)yy.size(); i++)
    writing_file << i << " " << yy[i] << std::endl;

  writing_file.close();
}

/// 2D graph export (T1,T2 is vector data)
template <typename T1, typename T2>
void exportGraph(T1 &xx, T2 &yy, std::string filename)
{
  assert((int)xx.size() == (int)yy.size());

  std::ofstream writing_file;
  writing_file.open(filename, std::ios::out);

  for(int i = 0; i < (int)xx.size(); i++)
    writing_file << xx[i] << " " << yy[i] << std::endl;

  writing_file.close();
}

/// 3D graph export (T1,T2 is vector data)
template <typename T1, typename T2, typename T3>
void exportGraph(T1 &xx, T2 &yy, T2 &zz, std::string filename)
{
  std::ofstream writing_file;
  writing_file.open(filename, std::ios::out);

  for(int i = 0; i < (int)xx.size(); i++)
    writing_file << xx[i] << " " << yy[i] << " " << zz[i] << std::endl;

  writing_file.close();
}

/// 2D add graph plot (T1,T2 is scalar data)
template <typename T1, typename T2>
void addPlot(T1 xx, T2 yy, std::string filename, int isFirst)
{
  std::ofstream writing_file;
  if(isFirst == 1) // step1...recreate file
    writing_file.open(filename, std::ios::out);
  else // other...add data to file
    writing_file.open(filename, std::ios::app);

  writing_file << xx << " " << yy << std::endl;

  writing_file.close();
}

/// 3D add graph plot (T1,T2,T3 is scalar data)
/// @param [in] isFirst ... true=1, false=other value
template <typename T1, typename T2, typename T3>
void addPlot(T1 &xx, T2 &yy, T3 &zz, std::string filename, int isFirst)
{
  std::ofstream writing_file;
  if(isFirst == 1)
    writing_file.open(filename, std::ios::out);
  else
    writing_file.open(filename, std::ios::app);

  writing_file << xx << " " << yy << " " << zz << std::endl;

  writing_file.close();
}
} // namespace icarat
