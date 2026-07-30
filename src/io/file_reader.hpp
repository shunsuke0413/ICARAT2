///  @file  file_reader.hpp
///  @author	Daiki Watanabe
///  @author	Takumi Sugiura
///  @date		January 22, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace icarat
{
//  read text and convert to vector of strings
template <class T> std::vector<std::string> read_file(T filePass)
{
  std::string line;
  std::vector<std::string> lines;

  std::ifstream input_file(filePass);
  if(!input_file.is_open())
  {
    std::cerr << "Could not open the file - '" << filePass << "'" << std::endl;
    exit(1);
  }

  while(getline(input_file, line))
    lines.push_back(line);

  input_file.close();

  return lines;
}

/**
 * @fn
 * @brief this funtion reads csv file as std::vector<std::vector<T>>. All colums
 * are read as type T.
 * @param path file path
 * @param hasHeader If file has hader in the first line, give true to skip it.
 * @param separator The default separator is ',' but you can give any character
 * as the separator
 */
template <typename T>
std::vector<std::vector<T>> read_csv(const std::string path, bool hasHeader,
                                     char separator = ',')
{

  std::ifstream fp(path, std::ios::in);
  if(!fp)
  {
    std::cout << "no such file or directory: " << path << std::endl;
    std::abort();
  }

  std::string buff_i;
  if(hasHeader)
  {
    std::getline(fp, buff_i);
  }

  // read first data
  std::getline(fp, buff_i);
  std::stringstream ss_first(buff_i);
  std::string buff_first;
  std::vector<T> first_tmp;
  while(std::getline(ss_first, buff_first, separator))
  {
    first_tmp.push_back(std::stod(buff_first));
  }

  std::vector<std::vector<T>> tmp(first_tmp.size());
  for(size_t i = 0, size = first_tmp.size(); i < size; i++)
  {
    tmp[i].push_back(first_tmp[i]);
  }
  while(std::getline(fp, buff_i))
  {
    std::stringstream ss(buff_i);
    std::string buff_j;
    std::vector<T> tmp_i;
    int counter_j = 0;
    while(std::getline(ss, buff_j, separator))
    {
      tmp[counter_j].push_back(std::stod(buff_j));
      counter_j++;
    }
  }
  return tmp;
}

} // namespace icarat
