#include "config.h"

#include <cmath>
#include <fstream>
#include <utility>

bool parseOutputTiles(const std::string& s, int& outX, int& outY) {
  std::string v = lower(trim(s));
  if (v.empty()) return false;
  if (v == "0") {
    outX = 0;
    outY = 0;
    return true;
  }
  auto xPos = v.find('x');
  if (xPos == std::string::npos) return false;
  std::string lhs = v.substr(0, xPos);
  std::string rhs = v.substr(xPos + 1);
  if (lhs.empty() || rhs.empty()) return false;
  try {
    outX = std::stoi(lhs);
    outY = std::stoi(rhs);
  } catch (...) {
    return false;
  }
  if (outX < 1 || outY < 1) return false;
  return true;
}

bool parseHexColorToGray16(const std::string& s, uint16_t& out) {
  std::string t = trim(s);
  if (!t.empty() && t[0] == '#') t = t.substr(1);
  if (t.size() != 6) return false;

  auto hex2 = [&](int i) -> int {
    auto h = [](char c) -> int {
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

void writeDefaultConfig(const fs::path& path) {
  std::ofstream f(path);
  f <<
    "# geotiff2png.ini\n\n"
    "interpolation = bicubic          # bilinear | bicubic\n"
    "max_proximity_pixels = 4096\n"
    "output_tiles = 0                 # 0 | NxM (e.g. 4x2)\n"
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

bool loadConfig(const fs::path& path, Config& cfg) {
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
    } else if (key == "output_tiles") {
      int tilesX = cfg.outputTilesX;
      int tilesY = cfg.outputTilesY;
      if (parseOutputTiles(val, tilesX, tilesY)) {
        cfg.outputTilesX = tilesX;
        cfg.outputTilesY = tilesY;
      }
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

  if (cfg.landMidHeight < cfg.landMinHeight) std::swap(cfg.landMidHeight, cfg.landMinHeight);
  if (cfg.maxProximityPixels < 1) cfg.maxProximityPixels = 1;

  return true;
}
