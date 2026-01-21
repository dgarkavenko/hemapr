#include "resample.h"

#include <cmath>

#include "helpers.h"

namespace {
float cubicCatmullRom(float p0, float p1, float p2, float p3, float t) {
  float a0 = -0.5f * p0 + 1.5f * p1 - 1.5f * p2 + 0.5f * p3;
  float a1 =  p0 - 2.5f * p1 + 2.0f * p2 - 0.5f * p3;
  float a2 = -0.5f * p0 + 0.5f * p2;
  float a3 =  p1;
  return ((a0 * t + a1) * t + a2) * t + a3;
}
}  // namespace

float sampleBilinear(const std::vector<float>& src, int sw, int sh, float x, float y) {
  int x0 = (int)std::floor(x);
  int y0 = (int)std::floor(y);
  int x1 = x0 + 1, y1 = y0 + 1;
  float tx = x - x0, ty = y - y0;

  x0 = clampi(x0, 0, sw - 1); x1 = clampi(x1, 0, sw - 1);
  y0 = clampi(y0, 0, sh - 1); y1 = clampi(y1, 0, sh - 1);

  float p00 = src[y0 * sw + x0];
  float p10 = src[y0 * sw + x1];
  float p01 = src[y1 * sw + x0];
  float p11 = src[y1 * sw + x1];

  float a = p00 + (p10 - p00) * tx;
  float b = p01 + (p11 - p01) * tx;
  return a + (b - a) * ty;
}

float sampleBicubic(const std::vector<float>& src, int sw, int sh, float x, float y) {
  int ix = (int)std::floor(x);
  int iy = (int)std::floor(y);
  float tx = x - ix;
  float ty = y - iy;

  float row[4];
  for (int m = -1; m <= 2; ++m) {
    float col[4];
    int yy = clampi(iy + m, 0, sh - 1);
    for (int n = -1; n <= 2; ++n) {
      int xx = clampi(ix + n, 0, sw - 1);
      col[n + 1] = src[yy * sw + xx];
    }
    row[m + 1] = cubicCatmullRom(col[0], col[1], col[2], col[3], tx);
  }
  return cubicCatmullRom(row[0], row[1], row[2], row[3], ty);
}
