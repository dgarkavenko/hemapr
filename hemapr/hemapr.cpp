// main.cpp
// Windows console app: GeoTIFF -> 16-bit grayscale PNG for Unreal Landscape
// Dependencies: GDAL, fmt (vcpkg), lodepng (vendored or adjust include accordingly)

#define FMT_HEADER_ONLY
#define FMT_UNICODE 0
#include <fmt/core.h>
#include <fmt/core.h>
#include <fmt/color.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>   // fflush
#include <filesystem>
#include <fstream>
#include <iostream> // only for std::getline
#include <limits>
#include <string>
#include <vector>

#include "gdal_priv.h"
#include "cpl_conv.h"

// If you installed lodepng via vcpkg, change this to: #include <lodepng.h>
#include "lodepng.h"

namespace fs = std::filesystem;

// ------------------------ small helpers ------------------------

static inline int clampi(int v, int lo, int hi) { return std::max(lo, std::min(hi, v)); }
static inline double clamp01(double v) { return std::max(0.0, std::min(1.0, v)); }
static inline double lerp(double a, double b, double t) { return a + (b - a) * t; }

static inline uint16_t to16From8(int v8) {
  v8 = clampi(v8, 0, 255);
  return static_cast<uint16_t>(v8 * 257); // exact 8->16 expansion
}

static std::string trim(const std::string& s) {
  size_t a = 0, b = s.size();
  while (a < b && std::isspace((unsigned char)s[a])) ++a;
  while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
  return s.substr(a, b - a);
}

static std::string lower(std::string s) {
  for (auto& c : s) c = (char)std::tolower((unsigned char)c);
  return s;
}

static bool hasTifExt(const fs::path& p) {
  auto ext = lower(p.extension().string());
  return ext == ".tif" || ext == ".tiff";
}

static int readIntWithDefault(const std::string& prompt, int def) {
  fmt::print("{} [{}]: ", prompt, def);
  std::string s;
  std::getline(std::cin, s);
  s = trim(s);
  if (s.empty()) return def;
  try { return std::stoi(s); } catch (...) { return def; }
}

// ------------------------ progress ------------------------

struct Progress {
  std::string label;
  int total = 1;
  int lastPct = -1;

  Progress(std::string l, int t) : label(std::move(l)), total(std::max(1, t)) {}

  void update(int done) {
    done = std::max(0, std::min(total, done));
    int pct = (int)std::llround(100.0 * (double)done / (double)total);
    if (pct != lastPct) {
      lastPct = pct;
      fmt::print("\r{:<26} {:3d}%", label, pct);
      std::fflush(stdout);
    }
  }

  void finish() {
    update(total);
    fmt::print("\n");
  }
};

// ------------------------ config ------------------------

enum class InterpMode { Bilinear = 1, Bicubic = 2 };

struct Config {
  InterpMode interp = InterpMode::Bicubic;
  int maxProximityPixels = 4096;

  double seaLevel = 0.0;          // <= seaLevel is sea
  double landMinHeight = 0.01;    // maps to landMin
  double landMidHeight = 2048.0;  // maps to landMid

  uint16_t seaShore = to16From8(128); // distance 0 -> #808080
  uint16_t seaFar   = to16From8(0);   // distance max -> #000000

  uint16_t landZero = to16From8(0);   // height seaLevel -> black-ish
  uint16_t landMin  = to16From8(128); // height landMinHeight -> #808080
  uint16_t landMid  = to16From8(191); // height landMidHeight -> #bfbfbf

  // Seabed curve shaping: "steep but smooth in-out"
  std::string seabedCurve = "smootherstep_pow"; // linear|smoothstep|smootherstep|smootherstep_pow|sigmoid
  double seabedGamma = 0.55;                    // for smootherstep_pow (0.4..0.9 usually nice)
  int shoreOffsetPixels = 1;                    // subtract to make first sea pixel match shore
};

static bool parseHexColorToGray16(const std::string& s, uint16_t& out) {
  std::string t = trim(s);
  if (!t.empty() && t[0] == '#') t = t.substr(1);
  if (t.size() != 6) return false;

  auto hex2 = [&](int i) -> int {
    auto h = [](char c)->int {
      if (c >= '0' && c <= '9') return c - '0';
      if (c >= 'a' && c <= 'f') return 10 + (c - 'a');
      if (c >= 'A' && c <= 'F') return 10 + (c - 'A');
      return -1;
    };
    int a = h(t[i]);
    int b = h(t[i + 1]);
    if (a < 0 || b < 0) return -1;
    return a * 16 + b;
  };

  int r = hex2(0), g = hex2(2), b = hex2(4);
  if (r < 0 || g < 0 || b < 0) return false;

  int gray8 = (r == g && g == b)
    ? r
    : (int)std::llround(0.2126 * r + 0.7152 * g + 0.0722 * b);

  out = to16From8(gray8);
  return true;
}

static void writeDefaultConfig(const fs::path& path) {
  std::ofstream f(path);
  f <<
    "# geotiff2png.ini\n\n"
    "interpolation = bicubic          # bilinear | bicubic\n"
    "max_proximity_pixels = 4096\n"
    "shore_offset_pixels = 1\n\n"
    "sea_level = 0.0                  # <= sea_level is sea\n"
    "land_min_height = 0.01\n"
    "land_mid_height = 2048.0\n\n"
    "seabed_curve = smootherstep_pow  # linear | smoothstep | smootherstep | smootherstep_pow | sigmoid\n"
    "seabed_gamma = 0.55\n\n"
    "sea_shore_color = #808080\n"
    "sea_far_color   = #000000\n\n"
    "land_zero_color = #000000\n"
    "land_min_color  = #808080\n"
    "land_mid_color  = #bfbfbf\n";
}

static bool loadConfig(const fs::path& path, Config& cfg) {
  if (!fs::exists(path)) return false;
  std::ifstream f(path);
  if (!f) return false;

  std::string line;
  while (std::getline(f, line)) {
    std::string s = trim(line);
    if (s.empty()) continue;
    if (s[0] == '#' || s[0] == ';') continue;

    auto eq = s.find('=');
    if (eq == std::string::npos) continue;

    std::string key = lower(trim(s.substr(0, eq)));
    std::string val = trim(s.substr(eq + 1));

    if (key == "interpolation") {
      std::string v = lower(val);
      if (v == "bilinear") cfg.interp = InterpMode::Bilinear;
      else if (v == "bicubic") cfg.interp = InterpMode::Bicubic;
    } else if (key == "max_proximity_pixels") {
      try { cfg.maxProximityPixels = std::stoi(val); } catch (...) {}
      cfg.maxProximityPixels = std::max(1, cfg.maxProximityPixels);
    } else if (key == "shore_offset_pixels") {
      try { cfg.shoreOffsetPixels = std::stoi(val); } catch (...) {}
      cfg.shoreOffsetPixels = std::max(0, cfg.shoreOffsetPixels);
    } else if (key == "sea_level") {
      try { cfg.seaLevel = std::stod(val); } catch (...) {}
    } else if (key == "land_min_height") {
      try { cfg.landMinHeight = std::stod(val); } catch (...) {}
    } else if (key == "land_mid_height") {
      try { cfg.landMidHeight = std::stod(val); } catch (...) {}
    } else if (key == "seabed_curve") {
      cfg.seabedCurve = lower(val);
    } else if (key == "seabed_gamma") {
      try { cfg.seabedGamma = std::stod(val); } catch (...) {}
      cfg.seabedGamma = std::max(0.05, std::min(10.0, cfg.seabedGamma));
    } else if (key == "sea_shore_color") {
      uint16_t c; if (parseHexColorToGray16(val, c)) cfg.seaShore = c;
    } else if (key == "sea_far_color") {
      uint16_t c; if (parseHexColorToGray16(val, c)) cfg.seaFar = c;
    } else if (key == "land_zero_color") {
      uint16_t c; if (parseHexColorToGray16(val, c)) cfg.landZero = c;
    } else if (key == "land_min_color") {
      uint16_t c; if (parseHexColorToGray16(val, c)) cfg.landMin = c;
    } else if (key == "land_mid_color") {
      uint16_t c; if (parseHexColorToGray16(val, c)) cfg.landMid = c;
    }
  }

  // Basic sanity
  if (cfg.landMidHeight < cfg.landMinHeight) std::swap(cfg.landMidHeight, cfg.landMinHeight);
  if (cfg.maxProximityPixels < 1) cfg.maxProximityPixels = 1;

  return true;
}

// ------------------------ seabed curve ------------------------

static double smoothstep(double t) {
  t = clamp01(t);
  return t * t * (3.0 - 2.0 * t);
}

static double smootherstep(double t) {
  t = clamp01(t);
  return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
}

static double sigmoid01(double t, double k) {
  t = clamp01(t);
  auto s = [&](double x) { return 1.0 / (1.0 + std::exp(-k * (x - 0.5))); };
  double a = s(0.0), b = s(1.0);
  double y = (s(t) - a) / (b - a);
  return clamp01(y);
}

// “steep but smooth in-out” lives here
static double seabedEase(double t, const Config& cfg) {
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
    // map gamma to a usable steepness range
    double g = std::max(0.05, std::min(10.0, cfg.seabedGamma));
    double k = std::max(0.5, std::min(20.0, 8.0 / g));
    return sigmoid01(t, k);
  }

  return smootherstep(t);
}

// ------------------------ resampling ------------------------

static float cubicCatmullRom(float p0, float p1, float p2, float p3, float t) {
  float a0 = -0.5f * p0 + 1.5f * p1 - 1.5f * p2 + 0.5f * p3;
  float a1 =  p0 - 2.5f * p1 + 2.0f * p2 - 0.5f * p3;
  float a2 = -0.5f * p0 + 0.5f * p2;
  float a3 =  p1;
  return ((a0 * t + a1) * t + a2) * t + a3;
}

static float sampleBilinear(const std::vector<float>& src, int sw, int sh, float x, float y) {
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

static float sampleBicubic(const std::vector<float>& src, int sw, int sh, float x, float y) {
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

// ------------------------ distance transform ------------------------
// Felzenszwalb & Huttenlocher, 1D squared Euclidean distance transform

static void distanceTransform1D(const double* f, double* d, int n) {
  const double INF = 1e30;

  std::vector<int> v(n);
  std::vector<double> z(n + 1);

  int k = 0;
  v[0] = 0;
  z[0] = -INF;
  z[1] =  INF;

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

static std::vector<double> distanceToLandSquared(
  const std::vector<uint8_t>& landMask, int w, int h,
  Progress* progCols, Progress* progRows)
{
  const double INF = 1e30;

  std::vector<double> tmp((size_t)w * (size_t)h, INF);
  std::vector<double> dist2((size_t)w * (size_t)h, INF);

  std::vector<double> f(std::max(w, h));
  std::vector<double> d(std::max(w, h));

  // Columns
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

  // Rows
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

static void writeU16BE(std::vector<unsigned char>& buf, size_t idx, uint16_t v) {
  buf[idx * 2 + 0] = (unsigned char)((v >> 8) & 0xFF);
  buf[idx * 2 + 1] = (unsigned char)(v & 0xFF);
}

// ------------------------ main ------------------------

int main() {
  GDALAllRegister();

  // Config
  fs::path cfgPath = fs::current_path() / "geotiff2png.ini";
  Config cfg;

  if (!fs::exists(cfgPath)) {
    writeDefaultConfig(cfgPath);
    fmt::print("Created default config: {}\n", cfgPath.string());
  }
  if (loadConfig(cfgPath, cfg)) {
    fmt::print("Loaded config: {}\n", cfgPath.string());
  } else {
    fmt::print("Config read failed; using defaults.\n");
  }

  fmt::print(fg(fmt::color::medium_sea_green) | fmt::emphasis::bold,
  "Config: interp={}, maxProx={}px, seabedCurve={}, gamma={}\n",
             (cfg.interp == InterpMode::Bicubic ? "bicubic" : "bilinear"),
             cfg.maxProximityPixels, cfg.seabedCurve, cfg.seabedGamma);

  fs::path folder = fs::current_path().concat("/data/");

  std::vector<fs::path> tifs;
  for (auto& e : fs::directory_iterator(folder)) {
    if (e.is_regular_file() && hasTifExt(e.path())) tifs.push_back(e.path());
  }
  std::sort(tifs.begin(), tifs.end());

  if (tifs.empty()) {
    fmt::print(stderr, "No .tif/.tiff files found in: {}\n", folder.string());
    return 1;
  }

  fmt::print("\nFound GeoTIFFs:\n");
  int showN = std::min<int>(9, (int)tifs.size());
  
  for (int i = 0; i < showN; ++i) {
    auto fn = tifs[i].filename().string();
    fmt::print("[{}] {}\n", i + 1, fmt::styled(fn, fmt::fg(fmt::color::light_golden_rod_yellow)));
  }

  int choice = readIntWithDefault("\nChoose (1-9)", 1);
  choice = clampi(choice, 1, showN);
  fs::path inPath = tifs[choice - 1];

  fmt::print("\nOpening: {}\n", inPath.string());

  GDALDataset* ds = (GDALDataset*)GDALOpen(inPath.string().c_str(), GA_ReadOnly);
  if (!ds) {
    fmt::print(stderr, "Failed to open GeoTIFF.\n");
    return 1;
  }

  GDALRasterBand* band = ds->GetRasterBand(1);
  if (!band) {
    fmt::print(stderr, "Missing raster band 1.\n");
    GDALClose(ds);
    return 1;
  }

  const int sw = ds->GetRasterXSize();
  const int sh = ds->GetRasterYSize();
  fmt::print("Original size: {} x {}\n", sw, sh);

  // Keep proportions: ask for bigger side
  const int biggerDefault = std::max(sw, sh);
  int biggerSide = readIntWithDefault("Target bigger side (keeps aspect)", biggerDefault);
  if (biggerSide < 1) biggerSide = biggerDefault;

  const double scale = (double)biggerSide / (double)std::max(sw, sh);
  const int ow = std::max(1, (int)std::llround(sw * scale));
  const int oh = std::max(1, (int)std::llround(sh * scale));
  fmt::print("Output size: {} x {}\n", ow, oh);

  // Read full-res float32
  std::vector<float> src((size_t)sw * (size_t)sh);
  CPLErr err = band->RasterIO(GF_Read, 0, 0, sw, sh,
                              src.data(), sw, sh, GDT_Float32,
                              0, 0, nullptr);
  if (err != CE_None) {
    fmt::print(stderr, "RasterIO failed.\n");
    GDALClose(ds);
    return 1;
  }

  int hasNoData = 0;
  double noData = band->GetNoDataValue(&hasNoData);

  for (auto& v : src) {
    if (!std::isfinite(v)) { v = 0.0f; continue; }
    if (hasNoData && std::fabs((double)v - noData) < 1e-12) v = 0.0f;
  }

  // Resample
  std::vector<float> elev((size_t)ow * (size_t)oh, 0.0f);
  Progress progResample("Resampling", oh);

  for (int y = 0; y < oh; ++y) {
    progResample.update(y);
    float sy = ((y + 0.5f) * (float)sh / (float)oh) - 0.5f;
    for (int x = 0; x < ow; ++x) {
      float sx = ((x + 0.5f) * (float)sw / (float)ow) - 0.5f;
      float v = (cfg.interp == InterpMode::Bicubic)
        ? sampleBicubic(src, sw, sh, sx, sy)
        : sampleBilinear(src, sw, sh, sx, sy);
      elev[(size_t)y * (size_t)ow + (size_t)x] = std::isfinite(v) ? v : 0.0f;
    }
  }
  progResample.finish();

  // Land mask + max elevation
  std::vector<uint8_t> land((size_t)ow * (size_t)oh, 0);
  float maxElev = 0.0f;

  for (size_t i = 0; i < elev.size(); ++i) {
    if ((double)elev[i] > cfg.seaLevel) {
      land[i] = 1;
      maxElev = std::max(maxElev, elev[i]);
    }
  }
  fmt::print("Max land elevation (after resample): {:.3f}\n", (double)maxElev);

  // Distance-to-land
  Progress progCols("Distance-to-land (cols)", ow);
  Progress progRows("Distance-to-land (rows)", oh);
  std::vector<double> dist2 = distanceToLandSquared(land, ow, oh, &progCols, &progRows);

  // Map to 16-bit grayscale PNG buffer (big-endian)
  std::vector<unsigned char> pngBytes((size_t)ow * (size_t)oh * 2);
  Progress progMap("Mapping to 16-bit", oh);

  const double denomMin = std::max(1e-12, (cfg.landMinHeight - cfg.seaLevel));
  const double denomMid = std::max(1e-12, (cfg.landMidHeight - cfg.landMinHeight));
  const double maxProx  = (double)std::max(1, cfg.maxProximityPixels);

  for (int y = 0; y < oh; ++y) {
    progMap.update(y);
    for (int x = 0; x < ow; ++x) {
      size_t i = (size_t)y * (size_t)ow + (size_t)x;
      double h = (double)elev[i];

      uint16_t outV = 0;

      if (h > cfg.seaLevel) {
        // Land mapping:
        // seaLevel..landMinHeight -> landZero..landMin
        // landMinHeight..landMidHeight -> landMin..landMid
        // > landMidHeight -> landMid..white based on maxElev
        if (h <= cfg.landMinHeight) {
          double t = clamp01((h - cfg.seaLevel) / denomMin);
          outV = (uint16_t)std::llround(lerp((double)cfg.landZero, (double)cfg.landMin, t));
        } else if (h <= cfg.landMidHeight) {
          double t = clamp01((h - cfg.landMinHeight) / denomMid);
          outV = (uint16_t)std::llround(lerp((double)cfg.landMin, (double)cfg.landMid, t));
        } else {
          if ((double)maxElev <= cfg.landMidHeight + 1e-6) {
            outV = cfg.landMid;
          } else {
            double t = clamp01((h - cfg.landMidHeight) / std::max(1e-12, ((double)maxElev - cfg.landMidHeight)));
            outV = (uint16_t)std::llround(lerp((double)cfg.landMid, 65535.0, t));
          }
        }
      } else {
        // Sea/seabed mapping by proximity:
        // distance 0 -> seaShore (typically #808080)
        // distance max -> seaFar (typically black)
        // plus "steep but smooth in-out" curve shaping
        double dpx = std::sqrt(std::max(0.0, dist2[i]));
        dpx = std::max(0.0, dpx - (double)cfg.shoreOffsetPixels);

        double t = clamp01(dpx / maxProx);
        t = seabedEase(t, cfg);

        double v = lerp((double)cfg.seaShore, (double)cfg.seaFar, t);
        outV = (uint16_t)std::llround(std::max(0.0, std::min(65535.0, v)));
      }

      writeU16BE(pngBytes, i, outV);
    }
  }
  progMap.finish();

  fs::path outPath = inPath;
  outPath.replace_extension(".png");

  fmt::print("Writing 16-bit PNG: {}\n", outPath.string());
  unsigned error = lodepng::encode(outPath.string(), pngBytes, (unsigned)ow, (unsigned)oh, LCT_GREY, 16);
  if (error) {
    fmt::print(stderr, "PNG encode error {}: {}\n", error, lodepng_error_text(error));
    GDALClose(ds);
    return 1;
  }

  GDALClose(ds);
  fmt::print("Done.\n");
  return 0;
}
