#include "zstd_compress.h"

namespace zstd {

//  bool zstd_compress(vector<uint32_t>& v_text, std::vector<uint8_t>& cBuff, size_t& cSizeActual) {
        
//         vector<uint8_t> srcContent;
//         srcContent.resize(v_text.size()*4,0);

//         for(size_t i =0;i<v_text.size();++i){
//             uint32_t x = v_text[i];
//             srcContent[i*4] = x & 0xffu;
//             x >>= 8;
//             srcContent[i*4 + 1] = x & 0xffu;
//             x >>= 8;
//             srcContent[i*4 + 2] = x & 0xffu;
//             x >>= 8;
//             srcContent[i*4 + 3] = x & 0xffu;
//         }
//     // 获取压缩后的最大可能大小
//     size_t cSize = ZSTD_compressBound(srcContent.size());

//     // 将压缩缓冲区设置为压缩后的最大可能大小
//     cBuff.resize(cSize);

//     // 压缩字符串
//     cSizeActual = ZSTD_compress(cBuff.data(), cSize, srcContent.data(), srcContent.size(), 12);
//     cBuff.resize(cSizeActual);
//     // 如果压缩失败，则返回false
//     return !ZSTD_isError(cSizeActual);
//   }
 bool zstd_compress(const std::vector<uint8_t>& srcContent, std::vector<uint8_t>& cBuff) {
    size_t cSizeActual = 0;
    // 获取压缩后的最大可能大小
    size_t cSize = ZSTD_compressBound(srcContent.size());

    // 将压缩缓冲区设置为压缩后的最大可能大小
    cBuff.resize(cSize);

    // 压缩字符串
    cSizeActual = ZSTD_compress(cBuff.data(), cSize, srcContent.data(), srcContent.size(), 10);
    cBuff.resize(cSizeActual);
    // 如果压缩失败，则返回false
    return !ZSTD_isError(cSizeActual);
  }

  bool zstd_decompress(const std::vector<uint8_t>& cBuff, std::vector<uint8_t>& dBuff) {

    size_t dSizeActual = 0;
    // 获取解压后的大小
    size_t dSize = ZSTD_getDecompressedSize(cBuff.data(), cBuff.size());

    // 将解压缓冲区设置为解压后的大小
    dBuff.resize(dSize);

    // 解压字符串
    dSizeActual = ZSTD_decompress(dBuff.data(), dSize, cBuff.data(), cBuff.size());

    dBuff.resize(dSizeActual);

    // 如果解压失败，则返回false
    return !ZSTD_isError(dSizeActual);
  }
//     bool zstd_decompress(const std::vector<uint8_t>& cBuff, std::vector<uint32_t>& v_text, size_t& dSizeActual) {
    
//     vector<uint8_t> dBuff;
//     // 获取解压后的大小
//     size_t dSize = ZSTD_getDecompressedSize(cBuff.data(), cBuff.size());

//     // 将解压缓冲区设置为解压后的大小
//     dBuff.resize(dSize);

//     // 解压字符串
//     dSizeActual = ZSTD_decompress(dBuff.data(), dSize, cBuff.data(), cBuff.size());

//     dBuff.resize(dSizeActual);

//     v_text.resize(dSizeActual>>2,0);
        
//     for(size_t i = 0;i<v_text.size();i++){
// 		uint64_t shift = 0;
//         for (int j = 0; j < 4; ++j)
// 		{
// 			uint32_t c = dBuff[i*4+j];
                
//             v_text[i] += c << shift;
// 			shift += 8;

// 		}
        
//     }
//     // 如果解压失败，则返回false
//     return !ZSTD_isError(dSizeActual);
//   }
}


// using namespace std;

// const size_t CHUNK_SIZE = 4096;

// // 压缩函数
// void zstd_compress_thread(const vector<uint8_t>& srcContent, vector<uint8_t>& cBuff, size_t offset, size_t chunkSize) {
//     size_t cSizeActual = 0;
//     // 获取压缩后的最大可能大小
//     size_t cSize = ZSTD_compressBound(chunkSize);

//     // 将压缩缓冲区设置为压缩后的最大可能大小
//     cBuff.resize(cSize, 0);

//     // 压缩字符串
//     cSizeActual = ZSTD_compress(cBuff.data(), cSize, srcContent.data() + offset, chunkSize, 3);
//     cBuff.resize(cSizeActual);
// }

// // 多线程压缩
// bool zstd_compress(const vector<uint8_t>& srcContent, vector<uint8_t>& cBuff) {
//     size_t dataSize = srcContent.size();
//     size_t numChunks = dataSize / CHUNK_SIZE + ((dataSize % CHUNK_SIZE) ? 1 : 0);
//     vector<thread> threads(numChunks);

//     cBuff.clear();

//     for (size_t i = 0; i < numChunks; ++i) {
//         size_t offset = i * CHUNK_SIZE;
//         size_t chunkSize = min(CHUNK_SIZE, dataSize - offset);

//         cBuff.resize(cBuff.size() + ZSTD_compressBound(chunkSize));

//         // 为每个chunk启动一个线程压缩
//         threads[i] = thread(zstd_compress_thread, std::cref(srcContent), std::ref(cBuff), offset, chunkSize);
//     }

//     // 等待所有线程完成
//     for (auto& thread : threads) {
//         thread.join();
//     }

//     // 如果压缩失败，则返回false
//     return !ZSTD_isError(ZSTD_getErrorCode(0));
// }

// int main() {
//     vector<uint8_t> srcContent = {'h', 'e', 'l', 'l', 'o', ',', ' ', 'w', 'o', 'r', 'l', 'd', '!'};
//     vector<uint8_t> compressedContent;

//     if (zstd_compress(srcContent, compressedContent)) {
//         cout << "压缩成功！" << endl;
//         cout << "压缩前大小：" << srcContent.size() << endl;
//         cout << "压缩后大小：" << compressedContent.size() << endl;
//     } else {
//         cout << "压缩失败！" << endl;
//     }

//     return 0;
// }
