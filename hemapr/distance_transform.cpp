#include "distance_transform.h"

#include <algorithm>

namespace {
void distanceTransform1D(const double* f, double* d, int n) {
  const double INF = 1e30;

  std::vector<int> v(n);
  std::vector<double> z(n + 1);

  int k = 0;
  v[0] = 0;
  z[0] = -INF;
  z[1] = INF;

  for (int q = 1; q < n; ++q) {
    double s = ((f[q] + (double)q * q) - (f[v[k]] + (double)v[k] * v[k])) / (2.0 * (q - v[k]));
    while (s <= z[k]) {
      --k;
      s = ((f[q] + (double)q * q) - (f[v[k]] + (double)v[k] * v[k])) / (2.0 * (q - v[k]));
    }
    ++k;
    v[k] = q;
    z[k] = s;
    z[k + 1] = INF;
  }

  k = 0;
  for (int q = 0; q < n; ++q) {
    while (z[k + 1] < q) ++k;
    double dx = (double)q - (double)v[k];
    d[q] = dx * dx + f[v[k]];
  }
}
}  // namespace

std::vector<double> distanceToLandSquared(
  const std::vector<uint8_t>& landMask, int w, int h,
  Progress* progCols, Progress* progRows)
{
  const double INF = 1e30;

  std::vector<double> tmp((size_t)w * (size_t)h, INF);
  std::vector<double> dist2((size_t)w * (size_t)h, INF);

  std::vector<double> f(std::max(w, h));
  std::vector<double> d(std::max(w, h));

  if (progCols) progCols->update(0);
  for (int x = 0; x < w; ++x) {
    if (progCols) progCols->update(x);
    for (int y = 0; y < h; ++y) {
      f[y] = landMask[(size_t)y * (size_t)w + (size_t)x] ? 0.0 : INF;
    }
    distanceTransform1D(f.data(), d.data(), h);
    for (int y = 0; y < h; ++y) {
      tmp[(size_t)y * (size_t)w + (size_t)x] = d[y];
    }
  }
  if (progCols) progCols->finish();

  if (progRows) progRows->update(0);
  for (int y = 0; y < h; ++y) {
    if (progRows) progRows->update(y);
    for (int x = 0; x < w; ++x) {
      f[x] = tmp[(size_t)y * (size_t)w + (size_t)x];
    }
    distanceTransform1D(f.data(), d.data(), w);
    for (int x = 0; x < w; ++x) {
      dist2[(size_t)y * (size_t)w + (size_t)x] = d[x];
    }
  }
  if (progRows) progRows->finish();

  return dist2;
}
