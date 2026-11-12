#include "opt_finite_strain.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace icarat
{
namespace opt_finite_strain
{
using namespace Eigen;
using namespace std;

void fdm(FEM &fem, vector<FiniteNeoHookedamage> &element, vector<Node> &node,
         Force &force, ValueNR &dis, Optimize &optim,
         Filter<FiniteNeoHookedamage, Node> &filter, int numstep,
         const SparseMatrix<double> &weights, double beta, double V0)
{
  // All perturbations start from zero displacement and zero damage. Restore
  // primal data even when a perturbed solve or derivative check fails.
  const auto originalElement = element;
  const auto originalNode = node;
  const auto originalForce = force;
  const auto originalDis = dis;
  const auto originalOptim = optim;
  auto restore = [&]() {
    element = originalElement;
    node = originalNode;
    force = originalForce;
    dis = originalDis;
    optim = originalOptim;
    filter.densityFilter(optim.design_s, false);
  };
  struct Result
  {
    double objective, volume;
    vector<unsigned char> branches;
  };
  auto reanalyze = [&](const vector<double> &design) {
    auto filtered = filter.densityFilter(design, false);
    auto rho = filter.threshold(filtered, beta, 0.5);
    Result result{0.0, -1.0, {}};
    for(int e = 0; e < fem.nelx; ++e)
    {
      result.branches.push_back(rho[e] < 1.e-4 ? 0 : rho[e] > 1.0 ? 2 : 1);
      rho[e] = std::clamp(rho[e], 1.e-4, 1.0);
      element[e].setDesignS(rho[e]);
    }
    auto gamma = filter.threshold(rho, 50.0, 0.1);
    for(int e = 0; e < fem.nelx; ++e)
    {
      result.branches.push_back(gamma[e] < 1.e-4 ? 0 : gamma[e] > 1.0 ? 2 : 1);
      gamma[e] = std::clamp(gamma[e], 1.e-4, 1.0);
      element[e].setGamma(gamma[e]);
    }
    NewtonRaphson<FiniteNeoHookedamage, Node> computation(
        fem, element, node, force, dis, numstep, {});
    computation.initialization();
    initializeDamagePath(fem, element, node);
    DamageHistory history{captureDamageFrame(fem, element, node)};
    for(int n = 1; n <= numstep; ++n)
    {
      const auto tau = computation.tau();
      const auto input = filter.tauFilter(tau);
      const VectorXd expected = weights * Map<const VectorXd>(tau.data(), tau.size());
      if((expected - Map<const VectorXd>(input.data(), input.size())).norm() >
         1.e-10 * max(1.0, expected.norm()))
        throw runtime_error("Damage filter reconstruction differs from tauFilter.");
      force.fext = static_cast<double>(n) * force.fext_org / numstep;
      const int iterations = computation.solve(n, numstep, 1.e-10, 100,
                                                "Eigen_LDLT", 0,
                                                input);
      if(iterations >= 100)
        throw runtime_error("Finite-difference forward analysis did not converge.");
      history.push_back(captureDamageFrame(fem, element, node, input));
      for(int e = 0; e < fem.nelx; ++e)
      {
        const double trial = element[e].calculateTensileDamage(input[e]);
        result.branches.push_back(trial > history[n - 1].damage(e) ? 1 : 0);
        result.branches.push_back(trial <= 0.0 ? 0 : trial >= 0.999 ? 2 : 1);
      }
    }
    result.objective = -damageWork(fem, element, node, history);
    for(const auto &e : element)
      result.volume += e.design_s() * e.volume / V0;
    return result;
  };
  try
  {
    ofstream comparison("sensitivity_comparison.csv");
    if(!comparison)
      throw runtime_error("Cannot open sensitivity_comparison.csv.");
    comparison << "element,h,adjoint,finite_difference,relative_error,volume_adjoint,volume_finite_difference,difference_scheme,branch_status\n";
    comparison << setprecision(17);
    const auto baseline = reanalyze(originalOptim.design_s);
    if(abs(baseline.objective - originalOptim.object_f) >
       1.e-7 * max(1.0, abs(baseline.objective)))
      throw runtime_error("Unperturbed full-path reanalysis differs from the objective.");
    cout << "Full-path finite differences (negative trapezoidal work)\n";
    for(int e = 0; e < fem.nelx; ++e)
      for(double h : {1.e-5, 1.e-6})
      {
        auto plusDesign = originalOptim.design_s;
        auto minusDesign = plusDesign;
        const bool canPlus = plusDesign[e] + h <= originalOptim.xmax[e];
        const bool canMinus = minusDesign[e] - h >= originalOptim.xmin[e];
        if(!canPlus && !canMinus)
          throw runtime_error("Finite-difference step exceeds both design bounds.");
        Result plus = baseline, minus = baseline;
        if(canPlus)
        {
          plusDesign[e] += h;
          plus = reanalyze(plusDesign);
        }
        if(canMinus)
        {
          minusDesign[e] -= h;
          minus = reanalyze(minusDesign);
        }
        const double span = (canPlus && canMinus ? 2.0 : 1.0) * h;
        const double fd = (plus.objective - minus.objective) / span;
        const double fdVolume = (plus.volume - minus.volume) / span;
        const double analytic = originalOptim.dfds[e];
        const double error = abs(fd - analytic) / max(1.e-10, max(abs(fd), abs(analytic)));
        const bool smooth = plus.branches == minus.branches;
        comparison << e << ',' << h << ',' << analytic << ',' << fd << ','
                   << error << ',' << originalOptim.dgds[e] << ',' << fdVolume << ','
                   << (canPlus && canMinus ? "central" : "one-sided") << ','
                   << (smooth ? "smooth-branches" : "BRANCH-CHANGE") << endl;
        if(!comparison)
          throw runtime_error("Failed writing sensitivity_comparison.csv.");
        cout << setprecision(12) << "element=" << e << " h=" << h
             << " adjoint=" << analytic << " FD=" << fd
             << " relative_error=" << error
             << " volume_FD=" << fdVolume
             << " volume_adjoint=" << originalOptim.dgds[e]
             << (canPlus && canMinus ? " central" : " one-sided")
             << (smooth ? " smooth-branches" : " BRANCH-CHANGE") << '\n';
        // Report rather than assert at nonsmooth switches. At design bounds
        // the selected derivative is the interior-side derivative.
        if(smooth && h == 1.e-6 &&
           abs(fd - analytic) > 2.e-3 * max(1.e-6, max(abs(fd), abs(analytic))))
          throw runtime_error("Full-path objective finite difference did not match the adjoint.");
        if(smooth && h == 1.e-6 &&
           abs(fdVolume - originalOptim.dgds[e]) > 2.e-4 * max(1.e-6, abs(fdVolume)))
          throw runtime_error("Filtered volume finite difference did not match.");
      }
    restore();
  }
  catch(...)
  {
    restore();
    throw;
  }
}
} // namespace opt_finite_strain
} // namespace icarat
