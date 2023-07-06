/*
 This file is a part of GTC software distributed under GNU GPL 3 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */
#ifndef _BIT_MEMORY_H
#define _BIT_MEMORY_H

#include "defs.h"
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

typedef enum {mode_none, mode_file_read, mode_file_write, mode_file_read_ra, mode_mem_read, mode_mem_write, mode_mem_read_ra} t_mode;

// ********************************************************************************************
class CBitMemory {
public:
    uint8_t *mem_buffer;
    bool mem_buffer_ownership;
    int64_t mem_buffer_pos;
    int32_t word_buffer_pos;

    // int32_t byte_buffer_pos;
private:
    int64_t mem_buffer_size;
    uint32_t word_buffer;
    int32_t word_buffer_size;

    // uint32_t byte_buffer;
    // int32_t byte_buffer_size;
    t_mode mode;
    
public:
    uint32_t n_bit_mask[32];
    
    inline bool FlushPartialWordBuffer();
    inline bool FlushInputWordBuffer();
    
    CBitMemory();
    CBitMemory(const CBitMemory &y);
    ~CBitMemory();
    
    bool Open(uint8_t *p, int64_t size, bool force_open = false);
    bool Create(int64_t size = 1);
    bool Complete();
    bool Close();
    bool Restart();

    bool TakeOwnership();
    
    inline uint64_t GetPos(void);
    bool SetPos(int64_t pos);
    
    inline int32_t GetWordPos(void);
    
    inline bool PutBit(const uint32_t word);

    inline bool PutBits(uint32_t word, int32_t n_bits);
    inline bool PutBytes(const unsigned char *data, int64_t n_bytes);
    inline bool PutByte(const unsigned char byte);

    inline bool PutWord(const uint32_t data);

    
    inline bool GetBit(uint32_t &word);

    inline bool GetBits(uint32_t &word, uint32_t n_bits);

    inline bool GetByte(uint32_t &byte);
    inline bool GetWord(uint32_t &data);
    inline bool GetWord(int32_t &data);

    inline bool discardBits(uint32_t word);  //AD //while putting
    inline bool GetBitsAndDiscard(uint32_t n_bits); //while getting

};


// ********************************************************************************************
bool CBitMemory::GetBitsAndDiscard(uint32_t n_bits)
{
    uint32_t count;
   
    while(n_bits)
    {
        if(word_buffer_pos == 0)
        {
            count = n_bits/8;
            if(count > 1)
            {
                if(mem_buffer_pos + count - 1 >= mem_buffer_size)
                    return false;
                //  byte = mem_buffer[mem_buffer_pos++];
                mem_buffer_pos += count - 1;
                n_bits -= (count - 1) << 3;  //*8
            }
            if(!GetByte(word_buffer))
                return false;
            word_buffer_pos = 8;
        }
        
        if((int32_t) n_bits > word_buffer_pos)
        {
          //  word <<= word_buffer_pos;
          //  word += word_buffer & n_bit_mask[word_buffer_pos];
            n_bits -= word_buffer_pos;
            word_buffer_pos = 0;
        }
        else
        {
          //  word <<= n_bits;
            word_buffer_pos -= n_bits;
          //  word += (word_buffer >> word_buffer_pos) & n_bit_mask[n_bits];
            return true;
        }
    }
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::discardBits(uint32_t n_bits)  //while putting
{
    uint32_t word;
    while(n_bits)
    {
        if((int32_t) n_bits <= word_buffer_pos)
        {
            word_buffer_pos -= n_bits;
            word_buffer = word_buffer >> n_bits ;
            if(!word_buffer_pos)
                word_buffer = 0;
            return true;
        }
        else //(n_bits > word_buf_pos)
        {
            n_bits -= word_buffer_pos;
            word_buffer_pos = 0;
            word_buffer = 0;
        }
        
        mem_buffer_pos = mem_buffer_pos - 4;
        //memcpy(&word_buffer, &mem_buffer[mem_buffer_pos], 4);
        GetWord(word);
        mem_buffer_pos = mem_buffer_pos - 4;
        word_buffer = word;
        word_buffer_pos = 32;

    }
    return false;
}


// ********************************************************************************************
uint64_t CBitMemory::GetPos(void)
{
    return mem_buffer_pos;
}

// ********************************************************************************************
int32_t CBitMemory::GetWordPos(void)
{
    return word_buffer_pos;
}


// ********************************************************************************************
bool CBitMemory::PutBit(const uint32_t word)
{
    if(word_buffer_pos < word_buffer_size)
    {
        word_buffer <<= 1;
        word_buffer += word;
        ++word_buffer_pos;
    }
    else
    {
        PutWord(word_buffer);
        word_buffer_pos = 1;
        word_buffer = word;
    }
    
    return true;
};
// bool CBitMemory::PutBit(const uint32_t word)
// {
//     if(byte_buffer_pos < byte_buffer_size)
//     {
//         byte_buffer <<= 1;
//         byte_buffer += word;
//         ++byte_buffer_pos;
//     }
//     else
//     {
//         // PutWord(word_buffer);
//         PutByte(byte_buffer)
//         byte_buffer_pos = 1;
//         byte_buffer = word;
//     }
    
//     return true;
// };


// ********************************************************************************************
bool CBitMemory::PutBits(uint32_t word, int32_t n_bits)
{
    int32_t rest_bits = word_buffer_size - word_buffer_pos;
    if(n_bits >= rest_bits)
    {
        n_bits -= rest_bits;
        word_buffer <<= rest_bits;
        word_buffer += word >> n_bits;
        word &= n_bit_mask[n_bits];
        word_buffer_pos = 0;
        PutWord(word_buffer);
        word_buffer = 0;
    }
    
    word_buffer     <<= n_bits;
    word_buffer     += word;
    word_buffer_pos += n_bits;
    
    return true;
}

// ********************************************************************************************

bool CBitMemory::FlushPartialWordBuffer()
{
    word_buffer <<= (32 - word_buffer_pos) & 7;
    
    if(word_buffer_pos > 24)
        PutByte(word_buffer >> 24);
    if(word_buffer_pos > 16)
        PutByte((word_buffer >> 16) & 0xFF);
    if(word_buffer_pos > 8)
        PutByte((word_buffer >> 8) & 0xFF);
    if(word_buffer_pos > 0)
        PutByte(word_buffer & 0xFF);
    
    word_buffer     = 0;
    word_buffer_pos = 0;
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::FlushInputWordBuffer()
{
    word_buffer_pos = 0;
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::PutByte(const unsigned char byte)
{
    if(mem_buffer_pos + 1 > mem_buffer_size)
    {
        mem_buffer_size = (uint64_t) ((mem_buffer_pos + 1) * 1.5);
        uint8_t *new_mem_buffer = new uint8_t[mem_buffer_size];
        if(mem_buffer)
        {
            copy_n(mem_buffer, mem_buffer_pos, new_mem_buffer);
            delete[] mem_buffer;
        }
        mem_buffer = new_mem_buffer;
    }
    
    mem_buffer[mem_buffer_pos++] = byte;
    
    return true;
}


// ********************************************************************************************
bool CBitMemory::PutBytes(const unsigned char *data, int64_t n_bytes)
{
   
    if(mem_buffer_pos + n_bytes > mem_buffer_size)
    {
        mem_buffer_size = (uint64_t) ((mem_buffer_pos + n_bytes) * 1.5);
        uint8_t *new_mem_buffer = new uint8_t[mem_buffer_size];
        if(mem_buffer)
        {
            copy_n(mem_buffer, mem_buffer_pos, new_mem_buffer);
            delete[] mem_buffer;
        }
        mem_buffer = new_mem_buffer;
    }
    copy_n(data, n_bytes, mem_buffer+mem_buffer_pos);
    mem_buffer_pos += n_bytes;
    
    return true;
}


// ********************************************************************************************
bool CBitMemory::PutWord(const uint32_t data)
{
    PutByte(data >> 24);
    PutByte((data >> 16) & 0xFF);
    PutByte((data >> 8) & 0xFF);
    PutByte(data & 0xFF);
    
    return true;
}
// ********************************************************************************************
bool CBitMemory::GetBit(uint32_t &word)
{
    if(word_buffer_pos == 0)
    {
        if(!GetByte(word_buffer))
            return false;
        word_buffer_pos = 7;
        word = word_buffer >> 7;
    }
    else
        word = (word_buffer >> (--word_buffer_pos)) & 1;
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::GetBits(uint32_t &word, uint32_t n_bits)
{
    word = 0;
    while(n_bits)
    {
        if(word_buffer_pos == 0)
        {
            if(!GetByte(word_buffer))
                return false;
            word_buffer_pos = 8;
        }
        
        if((int32_t) n_bits > word_buffer_pos)
        {
            word <<= word_buffer_pos;
            word += word_buffer & n_bit_mask[word_buffer_pos];
            n_bits -= word_buffer_pos;
            word_buffer_pos = 0;
        }
        else
        {
            word <<= n_bits;
            word_buffer_pos -= n_bits;
            word += (word_buffer >> word_buffer_pos) & n_bit_mask[n_bits];
            return true;
        }
    }
    
    return true;
}





// ********************************************************************************************
bool CBitMemory::GetByte(uint32_t &byte)
{
    if(mem_buffer_pos >= mem_buffer_size)
        return false;
    byte = mem_buffer[mem_buffer_pos++];
    return true;
};


// ********************************************************************************************
bool CBitMemory::GetWord(uint32_t &data)
{
    uint32_t c = 0;
    bool r;
    
    r = GetByte(c);
    data = c;
    r &= GetByte(c);
    data = (data << 8) + c;
    r &= GetByte(c);
    data = (data << 8) + c;
    r &= GetByte(c);
    data = (data << 8) + c;
    
    return r;
}

// ********************************************************************************************
bool CBitMemory::GetWord(int32_t &data)
{
    uint32_t c;
    bool r;
    
    r = GetByte(c);
    data = c;
    r &= GetByte(c);
    data = (data << 8) + c;
    r &= GetByte(c);
    data = (data << 8) + c;
    r &= GetByte(c);
    data = (data << 8) + c;
    
    return r;
}
class CBuffer {
public:
	enum class buffer_t {none, flag, integer, real, text};

private:
	uint32_t max_size;
#ifndef IGNORE_OFFSET_IN_FIRST_BLOCK
	uint32_t offset;
#endif
	buffer_t type;

	vector<uint32_t> v_size;
	vector<uint8_t> v_data;

	bool is_function;
	bool is_no_data;		// buffer was set to empty vector

	// function_data_item_t fun;

	uint32_t v_size_pos;
	uint32_t v_data_pos;
public:
	CBuffer();
	~CBuffer();

	void SetMaxSize(uint32_t _max_size, uint32_t _offset);

	// Output buffer methods
	void WriteFlag(uint8_t f)
	{
		v_size.emplace_back(f);

		type = buffer_t::flag;
	}

	// void WriteInt(char* p, uint32_t size)
	// {
	// 	v_size.emplace_back(size);

	// 	v_data.insert(v_data.end(), p, p + 4 * size);

	// 	type = buffer_t::integer;
	// }
    void WriteInt(char* p, uint32_t size)
	{
		v_size.emplace_back(size>>2);

		v_data.insert(v_data.end(), p, p + size);

		type = buffer_t::integer;
	}

	void WriteInt64(int64_t x)
	{
		int sign = 0;

		if (x < 0)
		{
			sign = 1;
			x = -x;
		}

		uint8_t bytes[8];
		int no_bytes = 0;
		size_t tmp = (size_t)x;

		for (; tmp; ++no_bytes)
		{
			bytes[no_bytes] = tmp & 0xff;
			tmp >>= 8;
		}

		v_size.emplace_back(sign + no_bytes * 2);

		for (int i = no_bytes - 1; i >= 0; --i)
			v_data.emplace_back(bytes[i]);
	}
    void WriteReal(char* p, uint32_t size)
	{
		v_size.emplace_back(size>>2);

		v_data.insert(v_data.end(), p, p +size);

		type = buffer_t::real;
	}
	// void WriteReal(char* p, uint32_t size)
	// {
	// 	v_size.emplace_back(size);

	// 	v_data.insert(v_data.end(), p, p + 4 * size);

	// 	type = buffer_t::real;
	// }

	void WriteText(char* p, uint32_t size)
	{
		v_size.push_back(size);

		if (size)
			v_data.insert(v_data.end(), p, p + size);

		type = buffer_t::text;
	}

	void GetBuffer(vector<uint32_t>& _v_size, vector<uint8_t>& _v_data);
	
	bool IsFull(void)
	{
#ifndef IGNORE_OFFSET_IN_FIRST_BLOCK
		return v_data.size() + 4 * v_size.size() >= max_size + offset;
#else
		return v_data.size() + 4 * v_size.size() >= max_size;
#endif
	}

	// Input buffer methods
	void ReadFlag(uint8_t &flag)
	{
		if (v_size.empty())
		{
			flag = 0;

			return;
		}

		flag = (uint8_t)v_size[v_size_pos++];
	}
    void ReadInt(char* &p, uint32_t& size)
	{
        
		if (v_size.empty())
		{
			p = nullptr;

			return;
		}
		size = v_size[v_size_pos++]*4;
        
		if (size)
		{
			p = new char[size];
            // cout<<"ReadInt: "<<v_data.size()<<endl;
			copy_n(v_data.begin() + v_data_pos, size, p);
			v_data_pos += size;

            
		}
		else
			p = nullptr;
	}
	// void ReadInt(char* &p, uint32_t& size)
	// {
        
	// 	if (v_size.empty())
	// 	{
	// 		p = nullptr;

	// 		return;
	// 	}
	// 	size = v_size[v_size_pos++];
        
	// 	if (size)
	// 	{
	// 		p = new char[size * 4];
    //         // cout<<"ReadInt: "<<v_data.size()<<endl;
	// 		copy_n(v_data.begin() + v_data_pos, 4 * size, p);
	// 		v_data_pos += 4 * size;

            
	// 	}
	// 	else
	// 		p = nullptr;
	// }

	void ReadInt64(int64_t &x)
	{
		int sign = v_size[v_size_pos++];

		int no_bytes = sign / 2;
		sign &= 1;

		x = 0;

		for (int i = 0; i < no_bytes; ++i)
			x += ((int64_t)v_data[v_data_pos++]) << (8 * (no_bytes - i - 1));

		if (sign)
			x = -x;
	}
    void ReadReal(char* &p, uint32_t& size)
	{
		if (v_size.empty())
		{
			p = nullptr;

			return;
		}

		size = v_size[v_size_pos++]*4;

		if (size)
		{
			p = new char[size];

			copy_n(v_data.begin() + v_data_pos, size, p);
			v_data_pos += size;
		}
		else
			p = nullptr;
	}
	// void ReadReal(char* &p, uint32_t& size)
	// {
	// 	if (v_size.empty())
	// 	{
	// 		p = nullptr;

	// 		return;
	// 	}

	// 	size = v_size[v_size_pos++];

	// 	if (size)
	// 	{
	// 		p = new char[size * 4];

	// 		copy_n(v_data.begin() + v_data_pos, 4 * size, p);
	// 		v_data_pos += 4 * size;
	// 	}
	// 	else
	// 		p = nullptr;
	// }

	void ReadText(char* &p, uint32_t& size)
	{
		if (v_size.empty())
		{
			p = nullptr;

			return;
		}

		size = v_size[v_size_pos++];

		if (size)
		{
			p = new char[size + 1];

			copy_n(v_data.begin() + v_data_pos, size, p);
			p[size] = 0;
			v_data_pos += size;
		}
		else
			p = nullptr;
	}

	void SetBuffer(vector<uint32_t>& _v_size, vector<uint8_t>& _v_data);
    bool Restart();
	// void SetFunction(function_data_item_t& _fun);
	void FuncInt(char*& p, uint32_t& size, char* src_p, uint32_t src_size);
	void FuncReal(char*& p, uint32_t& size, char* src_p, uint32_t src_size);

	bool IsEmpty()
	{
		return v_size_pos >= v_size.size() && !is_no_data;
	}
};



class CVariantsBuffer {
public:

private:
	uint32_t max_size;
#ifndef IGNORE_OFFSET_IN_FIRST_BLOCK
	uint32_t offset;
#endif

	vector<uint8_t> v_data;


	// function_data_item_t fun;

	uint32_t v_data_pos;
public:
	CVariantsBuffer();
	~CVariantsBuffer();

	void SetMaxSize(uint32_t _max_size, uint32_t _offset);

	// Output buffer methods
    void append_str(vector<uint8_t>& v_comp, const string& x) {
        size_t old_size = v_comp.size();
        size_t new_size = old_size + x.size() + 1; // 加 1 表示空字符
        v_comp.reserve(new_size);
        copy(x.begin(), x.end(), std::back_inserter(v_comp));
        v_comp.emplace_back(0); // 在末尾添加空字符
    }
    
    template <typename T>
	void append(vector<uint8_t> &v_comp, T x)
	{
    	append_str(v_comp, std::to_string(x));
	};
	

	// void WriteInt64(int64_t x)
	// {
	// 	int sign = 0;

	// 	if (x < 0)
	// 	{
	// 		sign = 1;
	// 		x = -x;
	// 	}

	// 	uint8_t bytes[8];
	// 	int no_bytes = 0;
	// 	size_t tmp = (size_t)x;

	// 	for (; tmp; ++no_bytes)
	// 	{
	// 		bytes[no_bytes] = tmp & 0xff;
	// 		tmp >>= 8;
	// 	}

	// 	v_size.emplace_back(sign + no_bytes * 2);

	// 	for (int i = no_bytes - 1; i >= 0; --i)
	// 		v_data.emplace_back(bytes[i]);
	// }

	void GetBuffer(vector<uint8_t>& _v_data);
	
	bool IsFull(void)
	{
#ifndef IGNORE_OFFSET_IN_FIRST_BLOCK
		return v_data.size() >= max_size + offset;
#else
		return v_data.size() >= max_size;
#endif
	}
    // *******************************************************************************************************************************
    void read_str(const vector<uint8_t>& v_comp, size_t& pos, string& x)
    {
	    x.clear();
        auto null_pos = find(v_comp.begin() + pos, v_comp.end(), 0);
        x.assign(v_comp.begin() + pos, null_pos);
        pos = null_pos - v_comp.begin() + 1;
    }
    template<typename T>
	void read(vector<uint8_t>& v_comp, size_t& pos, T& x) {
    	string t;
    	read_str(v_comp, pos, t);
    	x = stoll(t);
	};

	// void ReadInt64(int64_t &x)
	// {
	// 	int sign = v_size[v_size_pos++];

	// 	int no_bytes = sign / 2;
	// 	sign &= 1;

	// 	x = 0;

	// 	for (int i = 0; i < no_bytes; ++i)
	// 		x += ((int64_t)v_data[v_data_pos++]) << (8 * (no_bytes - i - 1));

	// 	if (sign)
	// 		x = -x;
	// }



	void SetBuffer(vector<uint8_t>& _v_data);
    bool Restart();



};
#endif
