#pragma once

#include <cstdint>
#include <filesystem>
#include <string>

#include "helpers.h"

enum class InterpMode { Bilinear = 1, Bicubic = 2 };

struct Config {
  InterpMode interp = InterpMode::Bicubic;
  int maxProximityPixels = 4096;
  int outputTilesX = 0;
  int outputTilesY = 0;

  double seaLevel = 0.0;
  double landMinHeight = 0.01;
  double landMidHeight = 2048.0;

  uint16_t seaShore = to16From8(128);
  uint16_t seaFar = to16From8(0);

  uint16_t landZero = to16From8(0);
  uint16_t landMin = to16From8(128);
  uint16_t landMid = to16From8(191);

  std::string seabedCurve = "smootherstep_pow";
  double seabedGamma = 0.55;
  int shoreOffsetPixels = 1;
};

bool parseHexColorToGray16(const std::string& s, uint16_t& out);
void writeDefaultConfig(const fs::path& path);
bool loadConfig(const fs::path& path, Config& cfg);
