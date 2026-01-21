#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif
#ifndef FMT_UNICODE
#define FMT_UNICODE 0
#endif

#include <fmt/core.h>

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
