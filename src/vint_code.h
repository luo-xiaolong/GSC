#pragma once
#include <vector>
#include <cstdint>
#include "defs.h"
namespace vint_code {
    uint32_t ReadVint(std::vector<uint8_t>& buffer, size_t& pos);
    size_t WriteVint(uint32_t value, std::vector<uint8_t>& buffer);
    
    std::vector<uint8_t> EncodeArray(const std::vector<uint32_t>& arr);
    std::vector<uint32_t> DecodeArray(std::vector<uint8_t>& buffer);
};


