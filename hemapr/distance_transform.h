#pragma once

#include <cstdint>
#include <vector>

#include "progress.h"

std::vector<double> distanceToLandSquared(
  const std::vector<uint8_t>& landMask, int w, int h,
  Progress* progCols, Progress* progRows);
