#define FMT_HEADER_ONLY
#define FMT_UNICODE 0

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <fmt/color.h>
#include <fmt/core.h>

#include "cpl_conv.h"
#include "gdal_priv.h"

#include "config.h"
#include "distance_transform.h"
#include "helpers.h"
#include "progress.h"
#include "resample.h"
#include "seabed_curve.h"

#include "lodepng.h"

namespace fs = std::filesystem;

int main() {
  GDALAllRegister();

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

  const int biggerDefault = std::max(sw, sh);
  int biggerSide = readIntWithDefault("Target bigger side (keeps aspect)", biggerDefault);
  if (biggerSide < 1) biggerSide = biggerDefault;

  const double scale = (double)biggerSide / (double)std::max(sw, sh);
  const int ow = std::max(1, (int)std::llround(sw * scale));
  const int oh = std::max(1, (int)std::llround(sh * scale));
  fmt::print("Output size: {} x {}\n", ow, oh);

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

  std::vector<uint8_t> land((size_t)ow * (size_t)oh, 0);
  float maxElev = 0.0f;

  for (size_t i = 0; i < elev.size(); ++i) {
    if ((double)elev[i] > cfg.seaLevel) {
      land[i] = 1;
      maxElev = std::max(maxElev, elev[i]);
    }
  }
  fmt::print("Max land elevation (after resample): {:.3f}\n", (double)maxElev);

  Progress progCols("Distance-to-land (cols)", ow);
  Progress progRows("Distance-to-land (rows)", oh);
  std::vector<double> dist2 = distanceToLandSquared(land, ow, oh, &progCols, &progRows);

  std::vector<unsigned char> pngBytes((size_t)ow * (size_t)oh * 2);
  Progress progMap("Mapping to 16-bit", oh);

  const double denomMin = std::max(1e-12, (cfg.landMinHeight - cfg.seaLevel));
  const double denomMid = std::max(1e-12, (cfg.landMidHeight - cfg.landMinHeight));
  const double maxProx = (double)std::max(1, cfg.maxProximityPixels);

  for (int y = 0; y < oh; ++y) {
    progMap.update(y);
    for (int x = 0; x < ow; ++x) {
      size_t i = (size_t)y * (size_t)ow + (size_t)x;
      double h = (double)elev[i];

      uint16_t outV = 0;

      if (h > cfg.seaLevel) {
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
            double t = clamp01((h - cfg.landMidHeight) /
                               std::max(1e-12, ((double)maxElev - cfg.landMidHeight)));
            outV = (uint16_t)std::llround(lerp((double)cfg.landMid, 65535.0, t));
          }
        }
      } else {
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
  unsigned error = lodepng::encode(outPath.string(), pngBytes,
                                   (unsigned)ow, (unsigned)oh, LCT_GREY, 16);
  if (error) {
    fmt::print(stderr, "PNG encode error {}: {}\n", error, lodepng_error_text(error));
    GDALClose(ds);
    return 1;
  }

  GDALClose(ds);
  fmt::print("Done.\n");
  return 0;
}
