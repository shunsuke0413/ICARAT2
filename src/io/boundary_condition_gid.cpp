///  @file  boundary_condition_gid.cpp
///  @author  Daiki Watanabe
///  @date  February 26, 2022.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#include "boundary_condition_gid.hpp"
#include <fstream>
#include <iostream>
#include <misc/json.hpp>

using namespace std;
using json = nlohmann::json;

namespace icarat
{

GIDParser::GIDParser(string gidPath) : gidPath_(gidPath)
{
  // read json
  ifstream reading(gidPath_ + "/ProjectParameters.json", ios::in);
  json j;
  reading >> j;

  mdpaName_ = j["solver_settings"]["model_import_settings"]["input_filename"]
                  .get<string>();

  modelName_ = j["solver_settings"]["model_part_name"].get<string>();
}

vector<GIDConstraint> GIDParser::readDirich(std::vector<std::string> &mdpa)
{
  // read json
  ifstream reading(gidPath_ + "/ProjectParameters.json", ios::in);
  json j;
  reading >> j;

  vector<GIDConstraint> gconsts;
  int size = (int)j["processes"]["constraints_process_list"].size();
  for(int i = 0; i < size; i++)
  {
    json consti = j["processes"]["constraints_process_list"][i];
    GIDConstraint gconst;

    // name...structure name is deleted.
    string fullName = consti["Parameters"]["model_part_name"].get<string>();
    gconst.name = fullName.erase(0, modelName_.length() + 1);

    // flag
    vector<bool> flag = consti["Parameters"]["constrained"].get<vector<bool>>();
    for(int j = 0; j < (int)flag.size(); j++)
      gconst.flag.push_back(flag[j]);

    // value
    auto val = consti["Parameters"]["value"];
    for(int j = 0; j < (int)val.size(); j++)
      if(val[j].is_null())
        gconst.val.push_back(0.0);
      else
        gconst.val.push_back(val[j].get<double>());

    gconsts.push_back(gconst);
  }

  return gconsts;
}

vector<GIDLoad> GIDParser::readNeumann(std::vector<std::string> &mdpa)
{
  // read json
  ifstream reading(gidPath_ + "/ProjectParameters.json", ios::in);
  json j;
  reading >> j;

  string modelName = j["solver_settings"]["model_part_name"].get<string>();

  vector<GIDLoad> gloads;
  int size = (int)j["processes"]["loads_process_list"].size();
  for(int i = 0; i < size; i++)
  {
    json loadi = j["processes"]["loads_process_list"][i];
    GIDLoad gload;

    // name...structure name is deleted.
    string fullName = loadi["Parameters"]["model_part_name"].get<string>();
    gload.name = fullName.erase(0, modelName.length() + 1);

    // value
    double modulus = loadi["Parameters"]["modulus"].get<double>();
    gload.val = loadi["Parameters"]["direction"].get<vector<double>>();
    for(int j = 0; j < (int)gload.val.size(); j++)
      gload.val[j] *= modulus;

    gloads.push_back(gload);
  }

  return gloads;
}

} // namespace icarat

// struct GIDConstraint
// {
//   std::string name; ///< group name in mdpa file
//   std::vector<int> flag;
//   std::vector<double> val;
// };

// struct GIDLoad
// {
//   std::string name; ///< group name in mdpa file
//   std::vector<double> val;
// };

// int main()
// {
//   // read json
//   ifstream reading("ProjectParameters.json", ios::in);
//   json j;
//   reading >> j;

//   vector<GIDConstraint> gconsts;
//   int size = (int)j["processes"]["constraints_process_list"].size();
//   for(int i = 0; i < size; i++)
//   {
//     json constraint = j["processes"]["constraints_process_list"][i];
//     GIDConstraint gconst;

//     gconst.name = constraint["Parameters"]["model_part_name"].get<string>();

//     vector<bool> flag =
//         constraint["Parameters"]["constrained"].get<vector<bool>>();
//     for(int j = 0; j < (int)flag.size(); j++)
//       gconst.flag.push_back(flag[j]);

//     gconst.val = constraint["Parameters"]["value"].get<vector<double>>();

//     gconsts.push_back(gconst);
//   }

//   return 0;
// }