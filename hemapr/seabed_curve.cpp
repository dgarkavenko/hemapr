#include "seabed_curve.h"

#include <algorithm>
#include <cmath>

#include "helpers.h"

namespace {
double smoothstep(double t) {
  t = clamp01(t);
  return t * t * (3.0 - 2.0 * t);
}

double smootherstep(double t) {
  t = clamp01(t);
  return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
}

double sigmoid01(double t, double k) {
  t = clamp01(t);
  auto s = [&](double x) { return 1.0 / (1.0 + std::exp(-k * (x - 0.5))); };
  double a = s(0.0), b = s(1.0);
  double y = (s(t) - a) / (b - a);
  return clamp01(y);
}
}  // namespace

double seabedEase(double t, const Config& cfg) {
  t = clamp01(t);

  if (cfg.seabedCurve == "linear") return t;
  if (cfg.seabedCurve == "smoothstep") return smoothstep(t);
  if (cfg.seabedCurve == "smootherstep") return smootherstep(t);

  if (cfg.seabedCurve == "smootherstep_pow") {
    double g = std::max(0.05, std::min(10.0, cfg.seabedGamma));
    double w = std::pow(t, g);
    return smootherstep(w);
  }

  if (cfg.seabedCurve == "sigmoid") {
    double g = std::max(0.05, std::min(10.0, cfg.seabedGamma));
    double k = std::max(0.5, std::min(20.0, 8.0 / g));
    return sigmoid01(t, k);
  }

  return smootherstep(t);
}
