#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace fs = std::filesystem;

int clampi(int v, int lo, int hi);
double clamp01(double v);
double lerp(double a, double b, double t);
uint16_t to16From8(int v8);
std::string trim(const std::string& s);
std::string lower(std::string s);
bool hasTifExt(const fs::path& p);
int readIntWithDefault(const std::string& prompt, int def);
bool readYesNo(const std::string& prompt, bool def);
void writeU16BE(std::vector<unsigned char>& buf, size_t idx, uint16_t v);
