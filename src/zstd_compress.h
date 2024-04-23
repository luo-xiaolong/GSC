#pragma once

#include <iostream>
#include <vector>
// #include <zstd.h>

// #include <thread>
#include "../include/zstd-1.5.2/lib/zstd.h"
using namespace std;
namespace zstd {
    bool zstd_compress(const std::vector<uint8_t>& srcContent, std::vector<uint8_t>& cBuff);
    bool zstd_compress(vector<uint32_t>& v_text, std::vector<uint8_t>& cBuff, size_t& cSizeActual);
    bool zstd_decompress(const std::vector<uint8_t>& cBuff, std::vector<uint8_t>& dBuff);
    bool zstd_decompress(const std::vector<uint8_t>& cBuff, std::vector<uint32_t>& v_text, size_t& dSizeActual);
}