#include "opt_elast_plast.hpp"
namespace icarat
{
namespace opt_elast_plast
{
using namespace std;
using namespace Eigen;

namespace
{
double dSIMP(double density, double young1, double young2, double pp)
{
  assert(young2 >= young1);
  return (young2 - young1) * pp * pow(density, pp - 1.0);
}
} // namespace

void OptLinearElastPlast::sensitivity(FEM &fem, vector<Node> &node)
{
  if(mattype_ == "plane_strain")
    sensitivity_plane_strain(fem, node);
  else if(mattype_ == "3D")
    sensitivity_3d(fem, node);
  else
  {
    cerr << "mattype is incorrect in elasto_plast_iso.cpp" << endl;
    exit(1);
  }
}

////////////////////////////////////////////////////////////
////////////////////////////private/////////////////////////
////////////////////////////////////////////////////////////

void OptLinearElastPlast::sensitivity_plane_strain(FEM &fem, vector<Node> &node)
{
  // get jacobian
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
  MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  for(int i = 0; i < ne; i++)
    for(int j = 0; j < fem.ndim; j++)
      X.coeffRef(i, j) = node[nodeID[i]].x[j];

  Bmatrix bmatrix;
  double jac = 0.0;

  for(int ip = 0; ip < ipmax; ip++)
  {
    // make strain increment
    VectorXd strain_n(4);
    strain_n << strain_gp(ip)[0], strain_gp(ip)[1], 0.5 * strain_gp(ip)[2], 0.0;

    VectorXd strain_n0(4);
    strain_n0[0] = strain_old_[ip][0];
    strain_n0[1] = strain_old_[ip][1];
    strain_n0[2] = 0.5 * strain_old_[ip][2];
    strain_n0[3] = 0.0;

    VectorXd strain_incre(4);
    strain_incre[0] = strain_n[0] - strain_n0[0];
    strain_incre[1] = strain_n[1] - strain_n0[1];
    strain_incre[2] = strain_n[2] - strain_n0[2];
    strain_incre[3] = 0.0;

    // get jacobian
    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);

    // objective function
    double func =
        0.5 * jac *
        ((stress_gp(ip)[0] + stress_old_[ip][0]) * strain_incre[0] +
         (stress_gp(ip)[1] + stress_old_[ip][1]) * strain_incre[1] +
         (sig3(ip) + stress_old_[ip][3]) * strain_incre[3] +
         2.0 * (stress_gp(ip)[2] + stress_old_[ip][2]) * strain_incre[2]);

    func_new_ = func_new_ - func;

    double dyoung = dSIMP(design_s(), young1(), young2(), pp());

    double dEh = dSIMP(design_s(), Ehs()[0], Ehs()[1], pp());

    double dsigma_y = dSIMP(design_s(), sigma_ys()[0], sigma_ys()[1], pp());

    double hards = young() * Eh() / (young() - Eh());
    double dhards =
        (dyoung * Eh() + young() * dEh) / (young() - Eh()) -
        (young() * Eh() * (dyoung - dEh)) / pow(young() - Eh(), 2.0);

    double Ehi = betaH() * hards;
    double Ehk = (1.0 - betaH()) * hards;
    double G = (0.5 * young()) / (1.0 + poisson());
    double K = young() / (3.0 * (1.0 - 2.0 * poisson()));
    double dK = dyoung / (3.0 * (1.0 - 2.0 * poisson()));
    double strain_p_bar = sqrt(2.0 / 3.0) * gamma(ip);

    double xi_tri_bar = sqrt(3.0 / 2.0) * xi_tri_norm(ip);
    // No.1 G_object
    double G_obj = (0.5 * dyoung) / (1.0 + poisson());

    double tr_strain = strain_incre[0] + strain_incre[1] + strain_incre[3];

    VectorXd strain_incre1(4);
    strain_incre1[0] = strain_incre[0] - (1.0 / 3.0) * tr_strain;
    strain_incre1[1] = strain_incre[1] - (1.0 / 3.0) * tr_strain;
    strain_incre1[2] = strain_incre[2];
    strain_incre1[3] = strain_incre[3] - (1.0 / 3.0) * tr_strain;

    // No.7 cal deviation trial stress obj
    VectorXd stress_tri1_obj(4);
    stress_tri1_obj[0] = stress1_obj_[ip][0] + 2.0 * G_obj * strain_incre1[0];
    stress_tri1_obj[1] = stress1_obj_[ip][1] + 2.0 * G_obj * strain_incre1[1];
    stress_tri1_obj[2] = stress1_obj_[ip][2] + 2.0 * G_obj * strain_incre1[2];
    stress_tri1_obj[3] = -(stress_tri1_obj[0] + stress_tri1_obj[1]);

    // No.6 trial stress obj
    VectorXd stress_tri_obj(4);
    stress_tri_obj[0] =
        stress_obj_[ip][0] + dK * tr_strain + 2.0 * G_obj * strain_incre1[0];
    stress_tri_obj[1] =
        stress_obj_[ip][1] + dK * tr_strain + 2.0 * G_obj * strain_incre1[1];
    stress_tri_obj[2] = stress_obj_[ip][2] + 2.0 * G_obj * strain_incre1[2];
    stress_tri_obj[3] = poisson() * (stress_tri_obj[0] + stress_tri_obj[1]);

    // get stress_fin_obj
    VectorXd stress_fin_obj(4), stress_fin1_obj(4);
    double AA;
    // Elastic
    if(yield(ip) <= 0.0)
    {
      // No.8 cal deviation final stress obj
      stress_fin1_obj[0] = stress_tri1_obj[0];
      stress_fin1_obj[1] = stress_tri1_obj[1];
      stress_fin1_obj[2] = stress_tri1_obj[2];
      stress_fin1_obj[3] = stress_tri1_obj[3];

      // final stress objective
      stress_fin_obj[0] = stress_tri_obj[0];
      stress_fin_obj[1] = stress_tri_obj[1];
      stress_fin_obj[2] = stress_tri_obj[2];
      stress_fin_obj[3] = stress_tri_obj[3];
    }
    // Plastic */
    else
    {
      // No.2 stress_try_ver_obj
      AA = xi_tri(ip)[0] * stress_tri1_obj[0] +
           xi_tri(ip)[1] * stress_tri1_obj[1] +
           xi_tri(ip)[3] * stress_tri1_obj[3] +
           2.0 * xi_tri(ip)[2] * stress_tri1_obj[2];

      double stress_tri_bar_obj = (1.5 * AA) / xi_tri_bar;

      // No.3 incremental plastic strain obj
      double dk_dpstrain = Ehi;

      // for hards objective
      double Ehi_obj = dhards * betaH();

      double dk_dsi = dsigma_y + Ehi_obj * this->strain_p_bar(ip);

      AA = stress_tri_bar_obj - dk_dpstrain * strain_p_bar_obj_[ip] - dk_dsi -
           3.0 * G_obj * strain_p_bar;
      double BB = dk_dpstrain + 3.0 * G;

      double strain_p_bar_obj = AA / BB;

      // No.4 stress final bar obj
      double stress_fin_bar_obj = stress_tri_bar_obj -
                                  3.0 * G_obj * strain_p_bar -
                                  3.0 * G * strain_p_bar_obj;

      // No.8 cal deviation final stress obj
      for(int i = 0; i < 3; i++)
      {
        stress_fin1_obj[i] =
            (stress_fin_bar_obj * xi_tri(ip)[i] +
             stress_fin_bar(ip) * stress_tri1_obj[i]) /
                xi_tri_bar -
            (stress_fin_bar(ip) * stress_tri_bar_obj * xi_tri(ip)[i]) /
                pow(xi_tri_bar, 2.0);
      }
      stress_fin1_obj[3] = -(stress_fin1_obj[0] + stress_fin1_obj[1]);

      // No.9 cal Hydrostatic p
      double tr_stress_obj =
          stress_tri_obj[0] + stress_tri_obj[1] + stress_tri_obj[3];

      VectorXd P_obj(4);
      P_obj[0] = (1.0 / 3.0) * tr_stress_obj;
      P_obj[1] = (1.0 / 3.0) * tr_stress_obj;
      P_obj[2] = 0.0;
      P_obj[3] = (1.0 / 3.0) * tr_stress_obj;

      // final stress objective
      stress_fin_obj[0] = stress_fin1_obj[0] + P_obj[0];
      stress_fin_obj[1] = stress_fin1_obj[1] + P_obj[1];
      stress_fin_obj[2] = stress_fin1_obj[2] + P_obj[2];
      stress_fin_obj[3] = poisson() * (stress_fin_obj[0] + stress_fin_obj[1]);

      // update plastic strain objective
      strain_p_bar_obj_[ip] += strain_p_bar_obj;
    }

    AA =
        0.5 * (stress_fin_obj[0] + stress_obj_[ip][0]) * strain_incre[0] * jac +
        0.5 * (stress_fin_obj[1] + stress_obj_[ip][1]) * strain_incre[1] * jac +
        0.5 * (stress_fin_obj[3] + stress_obj_[ip][3]) * strain_incre[3] * jac +
        2.0 * 0.5 * (stress_fin_obj[2] + stress_obj_[ip][2]) * strain_incre[2] *
            jac;

    objective_ -= AA;

    // update final deviation stress objective
    stress1_obj_[ip][0] = stress_fin1_obj[0];
    stress1_obj_[ip][1] = stress_fin1_obj[1];
    stress1_obj_[ip][2] = stress_fin1_obj[2];
    stress1_obj_[ip][3] = stress_fin1_obj[3];

    // update final stress objective
    stress_obj_[ip][0] = stress_fin_obj[0];
    stress_obj_[ip][1] = stress_fin_obj[1];
    stress_obj_[ip][2] = stress_fin_obj[2];
    stress_obj_[ip][3] = stress_fin_obj[3];

    // stock old strain & stress
    for(int i = 0; i < 3; i++)
    {
      strain_old_[ip](i) = strain_gp(ip)(i);
      stress_old_[ip](i) = stress_gp(ip)(i);
    }
    stress_old_[ip](3) = sig3(ip);
  }
}

void OptLinearElastPlast::sensitivity_3d(FEM &fem, vector<Node> &node)
{
  // get jacobian
  MatrixXd Ne = MatrixXd::Zero(fem.ndim, fem.ndim * ne);
  MatrixXd Be = MatrixXd::Zero(fem.voigt, numdof);
  MatrixXd X = MatrixXd::Zero(ne, fem.ndim);
  for(int i = 0; i < ne; i++)
    for(int j = 0; j < fem.ndim; j++)
      X.coeffRef(i, j) = node[nodeID[i]].x[j];

  Bmatrix bmatrix;
  double jac = 0.0;

  for(int ip = 0; ip < ipmax; ip++)
  {
    // SIMP derivative
    double dyoung = dSIMP(design_s(), young1(), young2(), pp());
    double dEh = dSIMP(design_s(), Ehs()[0], Ehs()[1], pp());
    double dsigma_y = dSIMP(design_s(), sigma_ys()[0], sigma_ys()[1], pp());

    double hards = young() * Eh() / (young() - Eh());
    double dhards =
        (dyoung * Eh() + young() * dEh) / (young() - Eh()) -
        (young() * Eh() * (dyoung - dEh)) / pow(young() - Eh(), 2.0);

    VectorXd onefive(6);
    onefive << 1.0, 1.0, 1.0, 0.5, 0.5, 0.5;

    // strain increment
    VectorXd strain_n(6);
    strain_n[0] = strain_gp(ip)[0];
    strain_n[1] = strain_gp(ip)[1];
    strain_n[2] = strain_gp(ip)[2];
    strain_n[3] = 0.5 * strain_gp(ip)[3];
    strain_n[4] = 0.5 * strain_gp(ip)[4];
    strain_n[5] = 0.5 * strain_gp(ip)[5];

    VectorXd strain_n0(6);
    strain_n0[0] = strain_old_[ip][0];
    strain_n0[1] = strain_old_[ip][1];
    strain_n0[2] = strain_old_[ip][2];
    strain_n0[3] = 0.5 * strain_old_[ip][3];
    strain_n0[4] = 0.5 * strain_old_[ip][4];
    strain_n0[5] = 0.5 * strain_old_[ip][5];

    double xi_tri_bar = sqrt(3.0 / 2.0) * xi_tri_norm(ip);

    double G_obj = (0.5 * dyoung) / (1.0 + poisson());

    VectorXd strain_incre = strain_n - strain_n0;

    bmatrix.make(Be, Ne, jac, eType, ipmax, X, ip);

    // objective function
    double func = 0.0;
    for(int i = 0; i < 3; i++)
      func += 0.5 * (stress_gp(ip)[i] + stress_old_[ip][i]) *
              (strain_incre[i]) * jac;
    for(int i = 3; i < 6; i++)
      func += 2.0 * 0.5 * (stress_gp(ip)[i] + stress_old_[ip][i]) *
              (strain_incre[i]) * jac;

    func_new_ = func_new_ - func;

    double Ehi = betaH() * hards;
    double Ehk = (1.0 - betaH()) * hards;
    double G = (0.5 * young()) / (1.0 + poisson());
    double K = young() / (3.0 * (1.0 - 2.0 * poisson()));
    double dK = dyoung / (3.0 * (1.0 - 2.0 * poisson()));
    double strain_p_bar = sqrt(2.0 / 3.0) * gamma(ip);

    double tr_strain = strain_incre[0] + strain_incre[1] + strain_incre[2];

    // No.6 trial stress obj
    MatrixXd Detmp(6, 6);
    // elastic modulus tensor
    De(Detmp, young(), poisson());
    VectorXd st_inc = dyoung / young() * Detmp * strain_incre;

    for(int i = 3; i < 6; i++)
      st_inc[i] = 2.0 * st_inc[i];

    VectorXd stress_tri_obj = stress_obj_[ip] + st_inc;

    VectorXd onezero(6);
    onezero << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;
    VectorXd strain_incre1 = strain_incre - (1.0 / 3.0) * tr_strain * onezero;

    // No.7 cal deviation trial stress obj
    VectorXd stress_tri1_obj = stress1_obj_[ip] + 2.0 * G_obj * strain_incre1;

    // get stress_fin_obj
    VectorXd stress_fin_obj(6), stress_fin1_obj(6);
    double AA;
    // Elastic
    if(yield(ip) <= 0.0)
    {
      // No.8 cal deviation final stress obj
      stress_fin1_obj = stress_tri1_obj;

      // final stress objective
      stress_fin_obj = stress_tri_obj;
    }
    // Plastic */
    else
    {
      // No.2 stress_try_ver_obj
      AA = 0.0;
      for(int i = 0; i < 3; i++)
        AA += xi_tri(ip)[i] * stress_tri1_obj[i];
      for(int i = 3; i < 6; i++)
        AA += 2.0 * xi_tri(ip)[i] * stress_tri1_obj[i];
      double stress_tri_bar_obj = (1.5 * AA) / xi_tri_bar;

      // No.3 incremental plastic strain obj
      double dk_dpstrain = Ehi;

      // for hards objective
      double Ehi_obj = dhards * betaH();

      double dk_dsi = dsigma_y + Ehi_obj * this->strain_p_bar(ip);

      AA = stress_tri_bar_obj - dk_dpstrain * strain_p_bar_obj_[ip] - dk_dsi -
           3.0 * G_obj * strain_p_bar;
      double BB = dk_dpstrain + 3.0 * G;

      double strain_p_bar_obj = AA / BB;

      // No.4 stress final bar obj
      double stress_fin_bar_obj = stress_tri_bar_obj -
                                  3.0 * G_obj * strain_p_bar -
                                  3.0 * G * strain_p_bar_obj;

      // No.8 cal deviation final stress obj
      stress_fin1_obj = (stress_fin_bar_obj * xi_tri(ip) +
                         stress_fin_bar(ip) * stress_tri1_obj) /
                            xi_tri_bar -
                        (stress_fin_bar(ip) * stress_tri_bar_obj * xi_tri(ip)) /
                            pow(xi_tri_bar, 2.0);

      // No.9 cal Hydrostatic p
      double tr_stress_obj =
          stress_tri_obj[0] + stress_tri_obj[1] + stress_tri_obj[2];

      VectorXd P_obj = (1.0 / 3.0) * tr_stress_obj * onezero;

      // final stress objective
      stress_fin_obj = stress_fin1_obj + P_obj;

      // update plastic strain objective
      strain_p_bar_obj_[ip] += strain_p_bar_obj;
    }

    VectorXd onetwo(6);
    onetwo << 1.0, 1.0, 1.0, 2.0, 2.0, 2.0;
    VectorXd AAtmp = onetwo.array() *
                     (stress_fin_obj + stress_obj_[ip]).array() *
                     strain_incre.array();

    AA = (0.5 * AAtmp * jac).sum();

    objective_ -= AA;

    // update final deviation stress objective
    stress1_obj_[ip] = stress_fin1_obj;

    // update final stress objective
    stress_obj_[ip] = stress_fin_obj;

    // stock old strain & stress
    strain_old_[ip] = strain_gp(ip);
    stress_old_[ip] = stress_gp(ip);
  }
}
} // namespace opt_elast_plast
} // namespace icarat