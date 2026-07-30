#pragma once
#include <io/export_graph.hpp>
#include <io/file_reader.hpp>
#include <optimization/filter_helm.hpp>
#include <optimization/opt_base.hpp>
#include <problem/base.hpp>
#include <problem/elasticity.hpp>
#include <vector>

namespace icarat
{
namespace uncertainty
{
/// class of robust TO
class Uncertainty
{
public:
  Uncertainty(FEM &fem, int numsam, int numloadcase)
      : dfds(fem.nelx), dgds(fem.nelx), fem_(fem), numsam_(numsam),
        numloadcase_(numloadcase), mean_val_(numloadcase),
        mean_dir_(numloadcase), dev_val_(numloadcase), dev_dir_(numloadcase)
  {
  }

  void setParameter(const toml::value &config);
  void make_normal_dis();
  void print_rcomp(std::vector<Force> &force_unc, std::vector<Value> &dis_unc);
  void unc_sensitivity_macro(std::vector<Force> &force_unc,
                             std::vector<Value> &dis_unc,
                             std::vector<LinearElastic> &element,
                             std::vector<Node> &node, double V0);
  void fdm(std::vector<LinearElastic> &element, std::vector<Node> &node,
           std::vector<Force> &force_unc, std::vector<Value> &dis_unc,
           Optimize &optim, FilterHelmholtz<LinearElastic, Node> &filter);

  std::vector<double> dfds; // 目的関数の感度
  std::vector<double> dgds; // 制約関数の感度
  double func_new;

private:
  FEM &fem_;

  int numsam_; // サンプルサイズ

  // 荷重ケース数(実際の独立な確率変数xiとの関係は， xi=ndim*numloadcase_)
  int numloadcase_;

  // 最大値の要素番号（3次元問題でのみ使用）
  int max_index_;

  // 荷重の大きさの平均値と偏差
  std::vector<double> mean_val_, dev_val_;
  // 荷重の方向の平均値と偏差(ラジアン)
  std::vector<double> mean_dir_, dev_dir_;

  double alpha_; // 重み係数

  double Cmean_; // 平均コンプライアンスの平均値
  double Cvar_;  // 平均コンプライアンスの分散値

  // void print_rcompの変数
  Eigen::MatrixXd E_; // 平均コンプライアンスが従う平均値のデータ
  Eigen::MatrixXd V_; // 平均コンプライアンスが分散値のデータ

  Eigen::MatrixXd beta_ij_;
};
} // namespace uncertainty
} // namespace icarat
