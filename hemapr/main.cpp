#define FMT_HEADER_ONLY
#define FMT_UNICODE 0

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <optional>
#include <regex>
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

struct ExpectedBounds {
  double xmin;
  double xmax;
  double ymin;
  double ymax;
};

std::optional<ExpectedBounds> parseExpectedBounds(const fs::path& path) {
  std::string name = path.stem().string();
  std::replace(name.begin(), name.end(), ',', '.');
  std::regex re(R"(.*_(-?\d+(?:\.\d+)?)_(-?\d+(?:\.\d+)?)_(-?\d+(?:\.\d+)?)_(-?\d+(?:\.\d+)?)$)");
  std::smatch match;
  if (!std::regex_match(name, match, re)) return std::nullopt;
  try {
    ExpectedBounds bounds{
      std::stod(match[1].str()),
      std::stod(match[2].str()),
      std::stod(match[3].str()),
      std::stod(match[4].str())
    };
    return bounds;
  } catch (...) {
    return std::nullopt;
  }
}

bool getGeoTransformBounds(GDALDataset* ds, double& xmin, double& xmax, double& ymin, double& ymax) {
  double gt[6] = {0};
  if (ds->GetGeoTransform(gt) != CE_None) return false;
  const int sw = ds->GetRasterXSize();
  const int sh = ds->GetRasterYSize();
  double cornersX[4] = {
    gt[0],
    gt[0] + gt[1] * sw,
    gt[0] + gt[2] * sh,
    gt[0] + gt[1] * sw + gt[2] * sh
  };
  double cornersY[4] = {
    gt[3],
    gt[3] + gt[4] * sw,
    gt[3] + gt[5] * sh,
    gt[3] + gt[4] * sw + gt[5] * sh
  };
  xmin = *std::min_element(cornersX, cornersX + 4);
  xmax = *std::max_element(cornersX, cornersX + 4);
  ymin = *std::min_element(cornersY, cornersY + 4);
  ymax = *std::max_element(cornersY, cornersY + 4);
  return true;
}

void expandGeoTiffToBounds(const fs::path& inPath,
                           GDALDataset* ds,
                           double expXMin,
                           double expXMax,
                           double expYMin,
                           double expYMax) {
  double gt[6] = {0};
  if (ds->GetGeoTransform(gt) != CE_None) {
    fmt::print(stderr, "GeoTransform unavailable; cannot expand.\n");
    return;
  }

  const int sw = ds->GetRasterXSize();
  const int sh = ds->GetRasterYSize();
  const double px = gt[1];
  const double py = gt[5];

  if (px == 0.0 || py == 0.0) {
    fmt::print(stderr, "Invalid pixel size; cannot expand.\n");
    return;
  }

  double xmin;
  double xmax;
  double ymin;
  double ymax;
  if (!getGeoTransformBounds(ds, xmin, xmax, ymin, ymax)) {
    fmt::print(stderr, "Failed to read GeoTransform bounds.\n");
    return;
  }

  const double tol = std::max(std::abs(px), std::abs(py)) * 0.01;
  if (std::abs(xmin - expXMin) <= tol &&
      std::abs(xmax - expXMax) <= tol &&
      std::abs(ymin - expYMin) <= tol &&
      std::abs(ymax - expYMax) <= tol) {
    fmt::print("GeoTIFF already matches expected bounds.\n");
    return;
  }

  if (px <= 0.0 || py >= 0.0) {
    fmt::print(stderr, "Unsupported GeoTransform orientation; cannot expand safely.\n");
    return;
  }

  const double colF = (expXMin - gt[0]) / px;
  const double rowF = (expYMax - gt[3]) / py;
  const double colMaxF = (expXMax - gt[0]) / px;
  const double rowMaxF = (expYMin - gt[3]) / py;

  const int minCol = (int)std::floor(colF + 1e-9);
  const int minRow = (int)std::floor(rowF + 1e-9);
  int maxCol = (int)std::ceil(colMaxF - 1e-9);
  int maxRow = (int)std::ceil(rowMaxF - 1e-9);

  const int expectedWpx = (int)std::llround((expXMax - expXMin) / px);
  const int expectedHpx = (int)std::llround((expYMax - expYMin) / std::abs(py));
  if (expectedWpx <= 0 || expectedHpx <= 0) {
    fmt::print(stderr, "Computed invalid expected size; cannot expand.\n");
    return;
  }

  maxCol = minCol + expectedWpx;
  maxRow = minRow + expectedHpx;

  const int newW = expectedWpx;
  const int newH = expectedHpx;

  if (newW <= 0 || newH <= 0) {
    fmt::print(stderr, "Computed invalid target size; cannot expand.\n");
    return;
  }

  std::vector<float> src((size_t)sw * (size_t)sh);
  GDALRasterBand* band = ds->GetRasterBand(1);
  if (!band) {
    fmt::print(stderr, "Missing raster band 1; cannot expand.\n");
    return;
  }
  CPLErr err = band->RasterIO(GF_Read, 0, 0, sw, sh,
                              src.data(), sw, sh, GDT_Float32,
                              0, 0, nullptr);
  if (err != CE_None) {
    fmt::print(stderr, "RasterIO failed; cannot expand.\n");
    return;
  }

  int hasNoData = 0;
  double noData = band->GetNoDataValue(&hasNoData);
  const float fillValue = hasNoData ? (float)noData : 0.0f;
  std::vector<float> out((size_t)newW * (size_t)newH, fillValue);
  for (int y = 0; y < sh; ++y) {
    int dy = y - minRow;
    if (dy < 0 || dy >= newH) continue;
    for (int x = 0; x < sw; ++x) {
      int dx = x - minCol;
      if (dx < 0 || dx >= newW) continue;
      out[(size_t)dy * (size_t)newW + (size_t)dx] = src[(size_t)y * (size_t)sw + (size_t)x];
    }
  }

  fs::path outPath = inPath;
  outPath.replace_filename(inPath.stem().string() + "_expanded.tif");

  GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
  if (!driver) {
    fmt::print(stderr, "GTiff driver not available; cannot expand.\n");
    return;
  }

  GDALDataset* outDs = driver->Create(outPath.string().c_str(), newW, newH, 1, GDT_Float32, nullptr);
  if (!outDs) {
    fmt::print(stderr, "Failed to create expanded GeoTIFF.\n");
    return;
  }

  double newGt[6] = {gt[0] + (double)minCol * px, px, gt[2],
                     gt[3] + (double)minRow * py, gt[4], py};
  outDs->SetGeoTransform(newGt);
  const char* proj = ds->GetProjectionRef();
  if (proj && proj[0] != '\0') {
    outDs->SetProjection(proj);
  }

  GDALRasterBand* outBand = outDs->GetRasterBand(1);
  if (!outBand) {
    fmt::print(stderr, "Failed to access output band.\n");
    GDALClose(outDs);
    return;
  }
  if (hasNoData) {
    outBand->SetNoDataValue(noData);
  }

  err = outBand->RasterIO(GF_Write, 0, 0, newW, newH,
                          out.data(), newW, newH, GDT_Float32,
                          0, 0, nullptr);
  if (err != CE_None) {
    fmt::print(stderr, "Failed to write expanded GeoTIFF.\n");
    GDALClose(outDs);
    return;
  }

  GDALClose(outDs);

  fmt::print("Expanded GeoTIFF written: {}\n", outPath.string());
}

int main() {
  GDALAllRegister();

  fs::path cfgPath = fs::current_path() / "geotiff2png.ini";

  if (!fs::exists(cfgPath)) {
    writeDefaultConfig(cfgPath);
    fmt::print("Created default config: {}\n", cfgPath.string());
  }
  Config cfg;

  fs::path folder = fs::current_path().concat("/data/");
  bool processAnother = true;

  while (processAnother) {
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

    cfg = Config{};
    if (loadConfig(cfgPath, cfg)) {
      fmt::print("Loaded config: {}\n", cfgPath.string());
    } else {
      fmt::print("Config read failed; using defaults.\n");
    }

    std::string tilesLabel = (cfg.outputTilesX > 0 && cfg.outputTilesY > 0)
      ? fmt::format("{}x{}", cfg.outputTilesX, cfg.outputTilesY)
      : "0";
    std::string targetLabel = (cfg.targetWidth > 0 && cfg.targetHeight > 0)
      ? fmt::format("{}x{}", cfg.targetWidth, cfg.targetHeight)
      : "0";
    fmt::print(fg(fmt::color::medium_sea_green) | fmt::emphasis::bold,
               "Config: interp={}, maxProx={}px, tiles={}, target={}, seabedCurve={}, gamma={}\n",
               (cfg.interp == InterpMode::Bicubic ? "bicubic" : "bilinear"),
               cfg.maxProximityPixels, tilesLabel, targetLabel, cfg.seabedCurve, cfg.seabedGamma);

    GDALDataset* ds = (GDALDataset*)GDALOpen(inPath.string().c_str(), GA_ReadOnly);
    if (!ds) {
      fmt::print(stderr, "Failed to open GeoTIFF.\n");
      return 1;
    }

    auto expected = parseExpectedBounds(inPath);
    if (expected) {
      double xmin;
      double xmax;
      double ymin;
      double ymax;
      if (getGeoTransformBounds(ds, xmin, xmax, ymin, ymax)) {
        double gt[6] = {0};
        const bool hasGeoTransform = (ds->GetGeoTransform(gt) == CE_None);
        const double px = hasGeoTransform ? gt[1] : 0.0;
        const double py = hasGeoTransform ? gt[5] : 0.0;
        fmt::print("Expected bounds from filename: xmin={:.6f} xmax={:.6f} ymin={:.6f} ymax={:.6f}\n",
                   expected->xmin, expected->xmax, expected->ymin, expected->ymax);
        fmt::print("Actual bounds from GeoTIFF:     xmin={:.6f} xmax={:.6f} ymin={:.6f} ymax={:.6f}\n",
                   xmin, xmax, ymin, ymax);
        const double expectedWidth = expected->xmax - expected->xmin;
        const double expectedHeight = expected->ymax - expected->ymin;
        const double actualWidth = xmax - xmin;
        const double actualHeight = ymax - ymin;
        const double sizeTol = hasGeoTransform
          ? std::max(std::abs(px), std::abs(py)) * 0.01
          : 1e-9;
        const bool sizeMatches = std::abs(actualWidth - expectedWidth) <= sizeTol &&
                                 std::abs(actualHeight - expectedHeight) <= sizeTol;

        if (sizeMatches) {
          fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "Valid\n");
        } else {
          if (hasGeoTransform && px != 0.0 && py != 0.0) {
            const double extendLeft = std::max(0.0, xmin - expected->xmin);
            const double extendRight = std::max(0.0, expected->xmax - xmax);
            const double extendTop = std::max(0.0, expected->ymax - ymax);
            const double extendBottom = std::max(0.0, ymin - expected->ymin);

            const auto pxCount = [](double dist, double pixelSize) {
              return (int)std::llround(dist / std::abs(pixelSize));
            };

            const bool needsExtension = extendLeft > sizeTol || extendRight > sizeTol ||
                                        extendTop > sizeTol || extendBottom > sizeTol;
            if (needsExtension) {
              fmt::print("Expected bounds extend beyond actual by:\n");
              fmt::print("  left  {:.6f} ({} px)\n", extendLeft, pxCount(extendLeft, px));
              fmt::print("  right {:.6f} ({} px)\n", extendRight, pxCount(extendRight, px));
              fmt::print("  top   {:.6f} ({} px)\n", extendTop, pxCount(extendTop, py));
              fmt::print("  bottom {:.6f} ({} px)\n", extendBottom, pxCount(extendBottom, py));
            } else {
              fmt::print("Expected bounds are smaller than actual; expansion would crop.\n");
            }
          }

          bool modify = readYesNo("Modify GeoTIFF to expected bounds?", false);
          if (modify) {
            expandGeoTiffToBounds(inPath, ds,
                                  expected->xmin, expected->xmax,
                                  expected->ymin, expected->ymax);
          }
        }
      } else {
        fmt::print(stderr, "GeoTransform unavailable; cannot validate bounds.\n");
      }
    } else {
      fmt::print("Filename does not contain expected bounds (xmin_xmax_ymin_ymax).\n");
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

    const int targetW = cfg.targetWidth;
    const int targetH = cfg.targetHeight;
    const double scaleW = (targetW > 0) ? (double)targetW / (double)sw : 1.0;
    const double scaleH = (targetH > 0) ? (double)targetH / (double)sh : 1.0;
    const double scale = std::min(1.0, std::min(scaleW, scaleH));
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

    const bool hasTiles = cfg.outputTilesX > 0 && cfg.outputTilesY > 0 &&
                          (cfg.outputTilesX > 1 || cfg.outputTilesY > 1);
    if (!hasTiles) {
      fmt::print("Writing 16-bit PNG: {}\n", outPath.string());
      unsigned error = lodepng::encode(outPath.string(), pngBytes,
                                       (unsigned)ow, (unsigned)oh, LCT_GREY, 16);
      if (error) {
        fmt::print(stderr, "PNG encode error {}: {}\n", error, lodepng_error_text(error));
        GDALClose(ds);
        return 1;
      }
    } else {
      if (ow % cfg.outputTilesX != 0 || oh % cfg.outputTilesY != 0) {
        fmt::print(stderr, "Output size {}x{} not divisible by tiles {}x{}.\n",
                   ow, oh, cfg.outputTilesX, cfg.outputTilesY);
        GDALClose(ds);
        return 1;
      }
      const int tileW = ow / cfg.outputTilesX;
      const int tileH = oh / cfg.outputTilesY;
      fmt::print("Writing tiled PNGs ({}x{} tiles of {}x{} each).\n",
                 cfg.outputTilesX, cfg.outputTilesY, tileW, tileH);
      for (int ty = 0; ty < cfg.outputTilesY; ++ty) {
        for (int tx = 0; tx < cfg.outputTilesX; ++tx) {
          std::vector<unsigned char> tileBytes((size_t)tileW * (size_t)tileH * 2);
          for (int y = 0; y < tileH; ++y) {
            size_t srcRow = (size_t)(ty * tileH + y) * (size_t)ow + (size_t)(tx * tileW);
            size_t dstRow = (size_t)y * (size_t)tileW;
            const unsigned char* srcPtr = pngBytes.data() + srcRow * 2;
            unsigned char* dstPtr = tileBytes.data() + dstRow * 2;
            std::copy(srcPtr, srcPtr + (size_t)tileW * 2, dstPtr);
          }
          fs::path tilePath = outPath;
          tilePath.replace_filename(
            fmt::format("{}_{}x{}_{}_{}.png",
                        outPath.stem().string(),
                        cfg.outputTilesX, cfg.outputTilesY,
                        tx + 1, ty + 1));
          unsigned error = lodepng::encode(tilePath.string(), tileBytes,
                                           (unsigned)tileW, (unsigned)tileH, LCT_GREY, 16);
          if (error) {
            fmt::print(stderr, "PNG encode error {}: {}\n", error, lodepng_error_text(error));
            GDALClose(ds);
            return 1;
          }
        }
      }
    }

    GDALClose(ds);
    fmt::print("Done.\n");
    processAnother = readYesNo("Process another GeoTIFF?", true);
  }

  return 0;
}
