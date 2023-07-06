#include "vint_code.h"
#include <iostream>

namespace vint_code {
// 从输入缓冲区中读取一个VINT
// buffer: 输入缓冲区
// pos: 当前读取位置，函数结束后会修改为下一个读取位置
// 返回值: 解码得到的整数

uint32_t ReadVint(std::vector<uint8_t>& buffer, size_t& pos)
{
    if(buffer[pos] == '\0')
    {
        pos++;
        return 0;
    }
    uint8_t firstByte = buffer[pos++];
    uint8_t mask = 0x80;
    uint32_t value = firstByte & 0x7f;
    int shift = 7;
    while (firstByte & mask)
    {
        firstByte = buffer[pos++];
        // value = (value << 7) | (firstByte & 0x7f);
        value = ((firstByte & 0x7f) << shift) | value ;
        shift += 7;
        // mask <<= 7;
    }

    return value;
}

// 将一个VINT写入输出缓冲区
// value: 待写入的整数
// buffer: 输出缓冲区
// 返回值: 写入的字节数
size_t WriteVint(uint32_t value, std::vector<uint8_t>& buffer)
{
    size_t size = 0;

    while (value > 0x7f)
    {
        buffer.push_back((value & 0x7f) | 0x80);
        value >>= 7;
        size++;
    }
    if(value)
        buffer.push_back(value & 0x7f);
    else
        buffer.push_back('\0');
    size++;

    return size;
}
std::vector<uint8_t> EncodeArray(const std::vector<uint32_t>& arr)
{
    std::vector<uint8_t> buffer;

    buffer.reserve(arr.size() * 4);

    for (const auto& value : arr)
    {
        WriteVint(value, buffer);
    }

    return buffer;
}
std::vector<uint32_t> DecodeArray(std::vector<uint8_t> &buffer)
{
    size_t size = buffer.size();
    std::vector<uint32_t> arr;
    arr.reserve(size / 4);
    size_t pos = 0;

    while (pos < size)
    {
        arr.push_back(ReadVint(buffer, pos));
    }

    return arr;
}
};