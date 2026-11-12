#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace icarat
{
struct DamageLawValue
{
  double trial, derivative;
};
// One law for the primal model and the history adjoint. The derivative is
// zero on the clipped branches, including the clipping boundaries.
inline DamageLawValue evaluateDamageLaw(double tau)
{
  if(!std::isfinite(tau) || tau < 0.0)
    throw std::runtime_error("Invalid filtered damage energy norm.");
  if(tau < 0.35)
    return {0.0, 0.0};
  const double a = -1.22005, b = 52.3473, c = 65.1249, d = 16.550;
  const double denominator = c * tau + d;
  const double A = 1.0 - (a * tau + b) / denominator;
  const double g = 1.0 / (1.0 + std::exp(-5.0 * (tau - 0.35)));
  const double raw = A * g;
  return {std::clamp(raw, 0.0, 0.999),
          raw > 0.0 && raw < 0.999
              ? (b * c - a * d) / (denominator * denominator) * g +
                    5.0 * A * g * (1.0 - g)
              : 0.0};
}
} // namespace icarat
