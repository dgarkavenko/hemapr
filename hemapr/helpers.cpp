#include "helpers.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <iostream>

#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif
#ifndef FMT_UNICODE
#define FMT_UNICODE 0
#endif

#include <fmt/core.h>

int clampi(int v, int lo, int hi) { return std::max(lo, std::min(hi, v)); }
double clamp01(double v) { return std::max(0.0, std::min(1.0, v)); }
double lerp(double a, double b, double t) { return a + (b - a) * t; }

uint16_t to16From8(int v8) {
  v8 = clampi(v8, 0, 255);
  return static_cast<uint16_t>(v8 * 257);
}

std::string trim(const std::string& s) {
  size_t a = 0, b = s.size();
  while (a < b && std::isspace((unsigned char)s[a])) ++a;
  while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
  return s.substr(a, b - a);
}

std::string lower(std::string s) {
  for (auto& c : s) c = (char)std::tolower((unsigned char)c);
  return s;
}

bool hasTifExt(const fs::path& p) {
  auto ext = lower(p.extension().string());
  return ext == ".tif" || ext == ".tiff";
}

int readIntWithDefault(const std::string& prompt, int def) {
  fmt::print("{} [{}]: ", prompt, def);
  std::string s;
  std::getline(std::cin, s);
  s = trim(s);
  if (s.empty()) return def;
  try { return std::stoi(s); } catch (...) { return def; }
}

bool readYesNo(const std::string& prompt, bool def) {
  fmt::print("{} [{}]: ", prompt, def ? "Y/n" : "y/N");
  std::string s;
  std::getline(std::cin, s);
  s = lower(trim(s));
  if (s.empty()) return def;
  if (s == "y" || s == "yes") return true;
  if (s == "n" || s == "no") return false;
  return def;
}

void writeU16BE(std::vector<unsigned char>& buf, size_t idx, uint16_t v) {
  buf[idx * 2 + 0] = (unsigned char)((v >> 8) & 0xFF);
  buf[idx * 2 + 1] = (unsigned char)(v & 0xFF);
}
