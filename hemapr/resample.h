#pragma once

#include <vector>

float sampleBilinear(const std::vector<float>& src, int sw, int sh, float x, float y);
float sampleBicubic(const std::vector<float>& src, int sw, int sh, float x, float y);
