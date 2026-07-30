#include "uncertainty.hpp"
#include <analysis/linear.hpp>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <io/export_graph.hpp>
#include <iostream>
#include <optimization/filter.hpp>
#include <random>

namespace icarat
{
namespace uncertainty
{
using namespace std;

void Uncertainty::setParameter(const toml::value &config)
{
  const auto &config_unc = toml::find(config, "uncertainty");
  numsam_ = toml::find<int>(config_unc, "num_sample");

  // インプットする荷重の大きさと荷重の方向の平均値と偏差
  double inp_vval = toml::find<double>(config_unc, "inp_vval");
  double inp_vdir = toml::find<double>(config_unc, "inp_vdir");

  alpha_ = toml::find<double>(config_unc, "alpha");

  const auto neum =
      toml::find<toml::array>(config, "boundary_condition", "neumann");
  using neum_type =
      std::tuple<std::array<double, 3>, std::array<std::array<double, 2>, 3>>;

  int counter = 0;
  for(auto &i_neum : neum)
  {
    auto data = toml::get<neum_type>(i_neum);
    auto val = std::get<0>(data);
    double inp_mdir, inp_mval;

    // 極座標表示
    if(fem_.ndim == 2)
    {
      inp_mval = sqrt(pow(val[0], 2.0) + pow(val[1], 2.0));
      inp_mdir = acos(val[0] / inp_mval);
      assert(!isnan(inp_mdir));
    }
    else
    {
      inp_mval = sqrt(pow(val[0], 2.0) + pow(val[1], 2.0) + pow(val[2], 2.0));

      // 絶対値の最大値のイテレータを取得
      auto max_val = std::max_element(val.begin(), val.end(),
                                      [](double a, double b)
                                      { return std::abs(a) < std::abs(b); });

      // 最大値の要素番号を取得
      max_index_ = std::distance(val.begin(), max_val);

      inp_mdir = acos(*max_val / inp_mval);
      assert(!isnan(inp_mdir));
    }
    mean_val_[counter] = inp_mval;
    mean_dir_[counter] = inp_mdir;
    dev_val_[counter] = inp_vval;
    dev_dir_[counter] = M_PI * inp_vdir / 180;
    cout << "mean_val " << mean_val_[counter] << endl;
    cout << "mean_dir " << mean_dir_[counter] << endl;
    cout << "dev_val " << dev_val_[counter] << endl;
    cout << "dev_dir " << dev_dir_[counter] << endl;
  }

  // 2階のテンソル
  E_.resize(numloadcase_ * fem_.ndim, numloadcase_ * fem_.ndim);
  // 4階のテンソル
  V_.resize(pow(numloadcase_ * fem_.ndim, 2), pow(numloadcase_ * fem_.ndim, 2));
}

void Uncertainty::make_normal_dis()
{
  // xi: 荷重の確率変数ベクトル
  //[荷重ケース数][サンプルデータ数(=numsam_)][次元数(=ndim)]
  vector<vector<vector<double>>> xi(numloadcase_);

  for(int ranvar = 0; ranvar < numloadcase_; ranvar++)
  {
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());

    // 平均0.0、標準偏差1.0のガウス分布
    std::normal_distribution<> normal(0.0, 1.0);
    // [0.0, 1.0]の一様分布
    std::uniform_real_distribution<> uniform(0.0, 1.0);

    ///////////////////////////////////////////////////
    ///////////////////荷重の大きさ////////////////////
    ///////////////////////////////////////////////////
    vector<double> val(numsam_);
    double mu_val = 0.0;
    for(int i = 0; i < numsam_; i++)
    {
      double x = max(min(normal(engine), 10000.0), -10000.0);
      val[i] = mean_val_[ranvar] + dev_val_[ranvar] * x;
      mu_val += val[i];
    }
    // 平均値
    mu_val = mu_val / numsam_;
    // 偏差
    double sigma_val = 0.0;
    for(int i = 0; i < numsam_; i++)
      sigma_val += pow((val[i] - mu_val), 2);
    sigma_val = sigma_val / numsam_;
    sigma_val = sqrt(sigma_val);
    assert(!isnan(sigma_val));

    cout << "mean value of load magnitude: " << mu_val << endl
         << "deviation of load magnitude: " << sigma_val << endl;

    ///////////////////////////////////////////////////
    ///////////////////荷重の方向////////////////////
    ///////////////////////////////////////////////////
    double mu_dir = 0.0;
    std::vector<double> theta(numsam_), phi(numsam_);
    for(int i = 0; i < numsam_; i++)
    {
      double x1 = max(min(normal(engine), 10000.0), -10000.0);
      double x2 = uniform(engine);
      theta[i] = mean_dir_[ranvar] + dev_dir_[ranvar] * x1;
      phi[i] = 2.0 * M_PI * x2;
      mu_dir += theta[i];
    }
    // 平均値
    mu_dir = mu_dir / numsam_;
    // 偏差
    double sigma_dir = 0.0;
    for(int i = 0; i < numsam_; i++)
      sigma_dir += pow((theta[i] - mu_dir), 2);
    sigma_dir = sqrt(sigma_dir / numsam_);

    assert(!isnan(sigma_dir));

    cout << "mean value of load direction(°): " << mu_dir * 180.0 / M_PI << endl
         << "deviation of load direction(°): " << sigma_dir * 180.0 / M_PI
         << endl;

    vector<vector<double>> xi_tmp(numsam_, vector<double>(fem_.ndim));
    if(fem_.ndim == 2)
    {
      for(int i = 0; i < numsam_; i++)
      {
        xi_tmp[i][0] = val[i] * cos(theta[i]);
        xi_tmp[i][1] = val[i] * sin(theta[i]);
      }
    }
    else // 3D
    {
      for(int i = 0; i < numsam_; i++)
      {
        if(max_index_ == 0)
        {
          xi_tmp[i][1] = val[i] * sin(theta[i]) * sin(phi[i]);
          xi_tmp[i][2] = val[i] * sin(theta[i]) * cos(phi[i]);
          xi_tmp[i][0] = val[i] * cos(theta[i]);
        }
        else if(max_index_ == 2)
        {
          xi_tmp[i][0] = val[i] * sin(theta[i]) * sin(phi[i]);
          xi_tmp[i][1] = val[i] * sin(theta[i]) * cos(phi[i]);
          xi_tmp[i][2] = val[i] * cos(theta[i]);
        }
        else
        {
          xi_tmp[i][2] = val[i] * sin(theta[i]) * sin(phi[i]);
          xi_tmp[i][0] = val[i] * sin(theta[i]) * cos(phi[i]);
          xi_tmp[i][1] = val[i] * cos(theta[i]);
        }
      }
    }

    xi[ranvar] = xi_tmp;
  }

  if(fem_.ndim == 2)
  {
    // mean of xi(2乗)
    double Exx = 0.0, Eyy = 0.0, Exxxx = 0.0, Eyyyy = 0.0;
    for(int j = 0; j < numloadcase_; j++)
      for(int i = 0; i < numsam_; i++)
      {
        Exx += pow(xi[j][i][0], 2.0) / (numsam_ * numloadcase_);
        Eyy += pow(xi[j][i][1], 2.0) / (numsam_ * numloadcase_);
        Exxxx += pow(xi[j][i][0], 4.0) / (numsam_ * numloadcase_);
        Eyyyy += pow(xi[j][i][1], 4.0) / (numsam_ * numloadcase_);
      }

    double Vxxxx = Exxxx - Exx * Exx;
    double Vyyyy = Eyyyy - Eyy * Eyy;

    double Exy = 0.0, Exxxy = 0.0, Exxyy = 0.0, Exyyy = 0.0;
    for(int j = 0; j < numloadcase_; j++)
      for(int k = 0; k < numloadcase_; k++)
        for(int i = 0; i < numsam_; i++)
        {
          Exy += xi[j][i][0] * xi[k][i][1] /
                 (numsam_ * numloadcase_ * numloadcase_);
          Exxxy += pow(xi[j][i][0], 3.0) * xi[k][i][1] /
                   (numsam_ * numloadcase_ * numloadcase_);
          Exxyy += pow(xi[j][i][0], 2.0) * pow(xi[k][i][1], 2.0) /
                   (numsam_ * numloadcase_ * numloadcase_);
          Exyyy += xi[j][i][0] * pow(xi[k][i][1], 3.0) /
                   (numsam_ * numloadcase_ * numloadcase_);
        }

    double Vxxxy = Exxxy - Exx * Exy;
    double Vxxyy = Exxyy - Exx * Eyy;
    double Vxyyy = Exyyy - Eyy * Exy;

    // element of array E
    E_ << Exx, Exy, //
        Exy, Eyy;

    // element of array V
    V_ << Vxxxx, Vxxxy, Vxxxy, Vxxyy, // xxxx,xxxy,xxyx,xxyy
        Vxxxy, Vxxyy, Vxxyy, Vxyyy,   // xyxx,xyxy,xyyx,xyyy
        Vxxxy, Vxxyy, Vxxyy, Vxyyy,   // yxxx,yxxy,yxyx,yxyy
        Vxxyy, Vxyyy, Vxyyy, Vyyyy;   // yyxx,yyxy,yyyx,yyyy
  }
  else // 3D
  {
    // mean of xi(2乗)
    double Exx = 0.0, Eyy = 0.0, Ezz = 0.0;
    double Exxxx = 0.0, Eyyyy = 0.0, Ezzzz = 0.0;
    for(int j = 0; j < numloadcase_; j++)
    {
      for(int i = 0; i < numsam_; i++)
      {
        Exx += pow(xi[j][i][0], 2.0) / (numsam_ * numloadcase_);
        Eyy += pow(xi[j][i][1], 2.0) / (numsam_ * numloadcase_);
        Ezz += pow(xi[j][i][2], 2.0) / (numsam_ * numloadcase_);
        Exxxx += pow(xi[j][i][0], 4.0) / (numsam_ * numloadcase_);
        Eyyyy += pow(xi[j][i][1], 4.0) / (numsam_ * numloadcase_);
        Ezzzz += pow(xi[j][i][2], 4.0) / (numsam_ * numloadcase_);
      }
    }

    double Vxxxx = Exxxx - Exx * Exx;
    double Vyyyy = Eyyyy - Eyy * Eyy;
    double Vzzzz = Ezzzz - Ezz * Ezz;

    double Exy = 0.0, Exz = 0.0, Eyz = 0.0;
    double Exxxy = 0.0, Exxyy = 0.0, Exyyy = 0.0;
    double Exxxz = 0.0, Exxzz = 0.0, Exzzz = 0.0;
    double Eyyyz = 0.0, Eyyzz = 0.0, Eyzzz = 0.0;
    for(int j = 0; j < numloadcase_; j++)
      for(int k = 0; k < numloadcase_; k++)
        for(int i = 0; i < numsam_; i++)
        {
          Exy += xi[j][i][0] * xi[k][i][1] /
                 (numsam_ * numloadcase_ * numloadcase_);
          Eyz += xi[j][i][1] * xi[k][i][2] /
                 (numsam_ * numloadcase_ * numloadcase_);
          Exz += xi[j][i][0] * xi[k][i][2] /
                 (numsam_ * numloadcase_ * numloadcase_);

          Exxxy += pow(xi[j][i][0], 3.0) * xi[k][i][1] /
                   (numsam_ * numloadcase_ * numloadcase_);
          Exxyy += pow(xi[j][i][0], 2.0) * pow(xi[k][i][1], 2.0) /
                   (numsam_ * numloadcase_ * numloadcase_);
          Exyyy += xi[j][i][0] * pow(xi[k][i][1], 3.0) /
                   (numsam_ * numloadcase_ * numloadcase_);

          Eyyyz += pow(xi[j][i][1], 3.0) * xi[k][i][2] /
                   (numsam_ * numloadcase_ * numloadcase_);
          Eyyzz += pow(xi[j][i][1], 2.0) * pow(xi[k][i][2], 2.0) /
                   (numsam_ * numloadcase_ * numloadcase_);
          Eyzzz += xi[j][i][1] * pow(xi[k][i][2], 3.0) /
                   (numsam_ * numloadcase_ * numloadcase_);

          Exxxz += pow(xi[j][i][0], 3.0) * xi[k][i][2] /
                   (numsam_ * numloadcase_ * numloadcase_);
          Exxzz += pow(xi[j][i][0], 2.0) * pow(xi[k][i][2], 2.0) /
                   (numsam_ * numloadcase_ * numloadcase_);
          Exzzz += xi[j][i][0] * pow(xi[k][i][2], 3.0) /
                   (numsam_ * numloadcase_ * numloadcase_);
        }

    double Vxxxy = Exxxy - Exx * Exy;
    double Vxxyy = Exxyy - Exx * Eyy;
    double Vxyyy = Exyyy - Eyy * Exy;

    double Vyyyz = Eyyyz - Eyy * Eyz;
    double Vyyzz = Eyyzz - Eyy * Ezz;
    double Vyzzz = Eyzzz - Eyz * Ezz;

    double Vxxxz = Exxxz - Exx * Exz;
    double Vxxzz = Exxzz - Exx * Ezz;
    double Vxzzz = Exzzz - Exz * Ezz;

    double Exxyz = 0.0, Exyyz = 0.0, Exyzz = 0.0;
    for(int j = 0; j < numloadcase_; j++)
      for(int k = 0; k < numloadcase_; k++)
        for(int l = 0; l < numloadcase_; l++)
          for(int i = 0; i < numsam_; i++)
          {
            Exxyz += pow(xi[j][i][0], 2.0) * xi[k][i][1] * xi[l][i][2] /
                     (numsam_ * numloadcase_ * numloadcase_ * numloadcase_);
            Exyyz += xi[j][i][0] * pow(xi[k][i][1], 2.0) * xi[l][i][2] /
                     (numsam_ * numloadcase_ * numloadcase_ * numloadcase_);
            Exyzz += xi[j][i][0] * xi[k][i][1] * pow(xi[l][i][2], 2.0) /
                     (numsam_ * numloadcase_ * numloadcase_ * numloadcase_);
          }

    double Vxxyz = Exxyz - Exx * Eyz;
    double Vxyyz = Exyyz - Exy * Eyz;
    double Vxyzz = Exyzz - Exy * Ezz;

    // element of array E
    E_ << Exx, Exy, Exz, // xx,xy,xz
        Exy, Eyy, Eyz,   // yx,yy,yz
        Exz, Eyz, Ezz;   // zx,zy,zz

    // element of array V
    V_ << Vxxxx, Vxxxy, Vxxxz, Vxxxy, Vxxyy, Vxxyz, Vxxxz, Vxxyz,
        Vxxzz, // xxxx,xxxy,xxxz,xxyx,xxyy,xxyz,xxzx,xxzy,xxzz
        Vxxxy, Vxxyy, Vxxyz, Vxxyy, Vxyyy, Vxyyz, Vxxyz, Vxyyz,
        Vxyzz, // xyxx,xyxy,xyxz,xyyx,xyyy,xyyz,xyzx,xyzy,xyzz
        Vxxxz, Vxxyz, Vxxzz, Vxxyz, Vxyyz, Vxyzz, Vxxzz, Vxyzz,
        Vxzzz, // xzxx,xzxy,xzxz,xzyx,xzyy,xzyz,xzzx,xzzy,xzzz
        Vxxxy, Vxxyy, Vxxyz, Vxxyy, Vxyyy, Vxyyz, Vxxyz, Vxyyz,
        Vxyzz, // yxxx,yxxy,yxxz,yxyx,yxyy,yxyz,yxzx,yxzy,yxzz
        Vxxyy, Vxyyy, Vxyyz, Vxyyy, Vyyyy, Vyyyz, Vxyyz, Vyyyz,
        Vyyzz, // yyxx,yyxy,yyxz,yyyx,yyyy,yyyz,yyzx,yyzy,yyzz
        Vxxyz, Vxyyz, Vxyzz, Vxyyz, Vyyyz, Vyyzz, Vxyzz, Vyyzz,
        Vyzzz, // yzxx,yzxy,yzxz,yzyx,yzyy,yzyz,yzzx,yzzy,yzzz
        Vxxxz, Vxxyz, Vxxzz, Vxxyz, Vxyyz, Vxyzz, Vxxzz, Vxyzz,
        Vxzzz, // zxxx,zxxy,zxxz,zxyx,zxyy,zxyz,zxzx,zxzy,zxzz
        Vxxyz, Vxyyz, Vxyzz, Vxyyz, Vyyyz, Vyyzz, Vxyzz, Vyyzz,
        Vyzzz, // yzxx,zyxy,zyxz,zyyx,zyyy,zyyz,zyzx,zyzy,zyzz
        Vxxzz, Vxyzz, Vxzzz, Vxyzz, Vyyzz, Vyzzz, Vxzzz, Vyzzz,
        Vzzzz; // zzxx,zzxy,zzxz,zzyx,zzyy,zzyz,zzzx,zzzy,zzzz
  }
}

void Uncertainty::print_rcomp(vector<Force> &force_unc, vector<Value> &dis_unc)
{
  int numxi = fem_.ndim * numloadcase_;
  int numxi2 = numxi * numxi;
  int numxi3 = numxi2 * numxi;

  // mean part
  Cmean_ = 0.0;
  for(int i = 0; i < numxi; i++)
    for(int j = 0; j < numxi; j++)
      Cmean_ += dis_unc[j].val.dot(force_unc[i].fext) * E_(i, j);

  // variance part
  Cvar_ = 0.0;
  for(int i = 0; i < numxi; i++)
    for(int j = 0; j < numxi; j++)
      for(int k = 0; k < numxi; k++)
        for(int l = 0; l < numxi; l++)
          Cvar_ += dis_unc[j].val.dot(force_unc[i].fext) *
                   dis_unc[l].val.dot(force_unc[k].fext) *
                   V_(numxi * i + k, numxi * j + l);

  // 感度でCvar_の逆数をとるため，最小値は0より大きい必要がある．
  Cvar_ = max(1.0e-7, Cvar_);
  func_new = Cmean_ + alpha_ * sqrt(Cvar_);
  std::cout << "Cmean_ sum:" << Cmean_ << endl;
  std::cout << "Cvar_ sum:" << Cvar_ << endl;
  std::cout << "uncertain compliance:" << func_new << endl;
}

void Uncertainty::unc_sensitivity_macro(std::vector<Force> &force_unc,
                                        std::vector<Value> &dis_unc,
                                        std::vector<LinearElastic> &element,
                                        std::vector<Node> &node, double V0)
{
  int numxi = fem_.ndim * numloadcase_;
  int numxi2 = numxi * numxi;
  int numxi3 = numxi2 * numxi;

  beta_ij_.resize(numxi, numxi);
  beta_ij_.setZero();
  // first term
  beta_ij_ += E_;

  // second term
  for(int i = 0; i < numxi; i++)
    for(int j = 0; j < numxi; j++)
      for(int k = 0; k < numxi; k++)
        for(int l = 0; l < numxi; l++)
          beta_ij_(i, j) += alpha_ / sqrt(Cvar_) *
                            dis_unc[k].val.dot(force_unc[l].fext) *
                            V_(numxi * i + k, numxi * j + l);

  std::cout << "beta_ij: " << beta_ij_ << endl;

  for(int nel = 0; nel < fem_.nelx; nel++)
  {
    dfds[nel] = 0.0;
    for(int i = 0; i < numxi; i++)
    {
      for(int j = 0; j < numxi; j++)
      {
        double dens = element[nel].design_s();
        double young1 = element[nel].young1();
        double young2 = element[nel].young2();
        double p = element[nel].pp();

        // 要素節点変位ベクトルの作成
        Eigen::VectorXd dispi(element[nel].numdof);
        Eigen::VectorXd dispj(element[nel].numdof);
        for(int aa = 0; aa < element[nel].ne; aa++)
        {
          for(int bb = 0; bb < fem_.ndim; bb++)
          {
            int dof = node[element[nel].nodeID[aa]].dof[bb];
            if(dof < fem_.numeq)
            {
              dispi.coeffRef(fem_.ndim * aa + bb) =
                  dis_unc[i].val[node[element[nel].nodeID[aa]].dof[bb]];
              dispj.coeffRef(fem_.ndim * aa + bb) =
                  dis_unc[j].val[node[element[nel].nodeID[aa]].dof[bb]];
            }
            else // fixされている部分は0とする．
            {
              dispi.coeffRef(fem_.ndim * aa + bb) = 0.0;
              dispj.coeffRef(fem_.ndim * aa + bb) = 0.0;
            }
          }
        }

        double young = (1.0 - pow(dens, p)) * young1 + pow(dens, p) * young2;
        double dyoung = p * pow(dens, p - 1.0) * (young2 - young1);
        dfds[nel] -= beta_ij_(i, j) * dyoung / young * dispj.transpose() *
                     element[nel].Ke() * dispi;
      }
    }

    dgds[nel] = element[nel].volume / V0;
  }
}

void Uncertainty::fdm(std::vector<LinearElastic> &element,
                      std::vector<Node> &node, std::vector<Force> &force_unc,
                      std::vector<Value> &dis_unc, Optimize &optim,
                      FilterHelmholtz<LinearElastic, Node> &filter)
{

  cout << "          Sensitivity analysis by FDM          " << endl;

  // init
  double delta = 1.0e-07;
  int numxi = fem_.ndim * numloadcase_;
  vector<double> object(fem_.nelx);
  double func_old = func_new;
  auto design_s_old = optim.design_s;

  /// cal df / ds
  for(int nel = 0; nel < fem_.nelx; nel++)
  {
    optim.design_s[nel] += delta;

    auto design_filtered =
        filter.densityFilter(optim.design_s, "Eigen_CG", 0, false);
    double beta = 1.0, T = 0.5;
    auto design_threshold = filter.threshold(design_filtered, beta, T);

    // set threshold design variable
    for(int i = 0; i < fem_.nelx; i++)
      element[i].setDesignS(design_threshold[i]);

    /// cal f(s + delta)
    for(int i = 0; i < numxi; i++)
    {
      LinearAnalysis<LinearElastic, Node> compt(fem_, element, node,
                                                force_unc[i], dis_unc[i]);
      compt.solve("Eigen_CG", 0);
    }

    /// get  stress
    for(auto &e : element)
      e.makeMisesStress(fem_, node);

    print_rcomp(force_unc, dis_unc);
    double func_new1 = func_new;

    /// dfds
    object[nel] = (func_new1 - func_old) / delta;

    /// init
    optim.design_s = design_s_old;

    cout << "Element ID : " << nel << " dfds macro = " << object[nel] << endl;
  }

  // output sensivity graph
  exportGraph(object, "objective_fdm.csv");

  cout << "========================================================" << endl;
  cout << "Finish computation succecfully at fdm.cpp               " << endl;
  cout << "========================================================" << endl;
  exit(0);
}
} // namespace uncertainty
} // namespace icarat
