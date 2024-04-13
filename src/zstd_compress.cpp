#include "zstd_compress.h"

namespace zstd {


 bool zstd_compress(const std::vector<uint8_t>& srcContent, std::vector<uint8_t>& cBuff) {
    size_t cSizeActual = 0;

    size_t cSize = ZSTD_compressBound(srcContent.size());

 
    cBuff.resize(cSize);


    cSizeActual = ZSTD_compress(cBuff.data(), cSize, srcContent.data(), srcContent.size(), 10);
    cBuff.resize(cSizeActual);

    return !ZSTD_isError(cSizeActual);
  }

  bool zstd_decompress(const std::vector<uint8_t>& cBuff, std::vector<uint8_t>& dBuff) {

    size_t dSizeActual = 0;

    size_t dSize = ZSTD_getDecompressedSize(cBuff.data(), cBuff.size());


    dBuff.resize(dSize);


    dSizeActual = ZSTD_decompress(dBuff.data(), dSize, cBuff.data(), cBuff.size());

    dBuff.resize(dSizeActual);

    return !ZSTD_isError(dSizeActual);
  }
}
