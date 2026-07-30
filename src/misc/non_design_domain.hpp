#pragma once

#include <io/file_reader.hpp>
#include <problem/base.hpp>
#include <toml/toml.hpp>

namespace icarat
{
template <class E, class N>
std::pair<std::vector<int>, std::vector<int>>
getDesignNonDesignIDs(FEM &fem, std::vector<E> &element, std::vector<N> &node,
                      const toml::value &config, std::string project)
{
  using NonDesignType =
      std::tuple<std::array<double, 1>, std::array<std::string, 1>>;

  const auto non_designs =
      toml::find<toml::array>(config, "this_problem", "non_design");

  // 1. nondesignIDsの作成
  std::vector<int> nondesignIDs;
  for(auto &non : non_designs)
  {
    auto data = toml::get<NonDesignType>(non);
    double value = std::get<0>(data)[0];
    std::string nonIDs_path = project + "/" + std::get<1>(data)[0];
    auto nonIDs_str = read_file(nonIDs_path);

    for(auto &p : nonIDs_str)
    {
      int tmp = std::stoi(p);
      nondesignIDs.push_back(tmp);

      element[tmp].isDesignable = false;
    }
  }

  // 昇順にし，被りを削除
  std::sort(nondesignIDs.begin(), nondesignIDs.end());
  nondesignIDs.erase(std::unique(nondesignIDs.begin(), nondesignIDs.end()),
                     nondesignIDs.end());

  // 2. 全数からnondesignIDsの番号を消す＝designIDs
  std::vector<int> allNumbers(fem.nelx);
  std::iota(allNumbers.begin(), allNumbers.end(), 0);

  std::vector<int> designIDs;
  std::set_difference(allNumbers.begin(), allNumbers.end(),
                      nondesignIDs.begin(), nondesignIDs.end(),
                      std::back_inserter(designIDs));

  return std::make_pair(designIDs, nondesignIDs);
}
} // namespace icarat