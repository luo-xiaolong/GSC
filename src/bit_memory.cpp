/*
 This file is a part of GTC software distributed under GNU GPL 3 licence.
 
 Authors: Agnieszka Danek and Sebastian Deorowicz
 
 Version: 1
 Date   : 2017-April
 */
#include "defs.h"
#include "bit_memory.h"
#include <algorithm>
using namespace std;
// ********************************************************************************************
CBitMemory::CBitMemory()
{
    mem_buffer	     = NULL;
    mem_buffer_size  = 0;
    mem_buffer_pos   = 0;
    word_buffer_pos  = 0;
    word_buffer		 = 0;
    word_buffer_size = 32;
    mode			 = mode_none;
    mem_buffer_ownership = false;
    
    for(int32_t i = 0; i < 32; ++i)
        n_bit_mask[i] = (1u << i) - 1;
}
// ********************************************************************************************
CBitMemory::CBitMemory(const CBitMemory &y)
{
    CBitMemory &x = const_cast<CBitMemory &>(y);
    
    mem_buffer = x.mem_buffer;
    x.mem_buffer = NULL;
    mem_buffer_ownership = x.mem_buffer_ownership;
    x.mem_buffer_ownership = false;
    
    mem_buffer_pos = x.mem_buffer_pos;
    x.mem_buffer_pos = 0;
    
    mem_buffer_size = x.mem_buffer_size;
    x.mem_buffer_size = 0;
    
    word_buffer = x.word_buffer;
    word_buffer = 0;
    
    word_buffer_pos = x.word_buffer_pos;
    x.word_buffer_pos = 0;
    
    word_buffer_size = x.word_buffer_size;
    
    mode = x.mode;
    mode = mode_none;
    
    copy_n(x.n_bit_mask, 32, n_bit_mask);
}
// ********************************************************************************************
CBitMemory::~CBitMemory()
{
    if(mem_buffer && mem_buffer_ownership)
        delete[] mem_buffer;
}
// ********************************************************************************************
bool CBitMemory::TakeOwnership()
{
    if(mem_buffer && mem_buffer_ownership)
    {
        mem_buffer_ownership = false;
        return true;
    }
    
    return false;
}
// ********************************************************************************************
bool CBitMemory::Open(uint8_t *p, int64_t size, bool force_open)
{
    if(!force_open)
    {
        if(mode != mode_none)
            return false;
    }
    
    if(mem_buffer && mem_buffer_ownership)
        delete[] mem_buffer;
    
    mem_buffer_size = size;
    /*	if(!mem_buffer_size)
     mem_buffer_size = 1;
     mem_buffer = new uint8_t[mem_buffer_size];
     copy_n(p, size, mem_buffer);*/
    mem_buffer = p;
    mem_buffer_ownership = false;
    if(!mem_buffer_size)
    {
        mem_buffer_size = 1;
        mem_buffer = new uint8_t[mem_buffer_size];
        mem_buffer_ownership = true;
        copy_n(p, size, mem_buffer);
    }
    
    mode = mode_mem_read;
    
    word_buffer_size = 8;
    mem_buffer_pos   = 0;
    
    return mode == mode_mem_read;
}
// ********************************************************************************************
bool CBitMemory::Restart()
{
    mode = mode_mem_read;
    
    word_buffer_size = 32;
    mem_buffer_pos   = 0;
    word_buffer_pos  = 0;
    word_buffer		 = 0;
    
    return true;
}

// ********************************************************************************************
bool CBitMemory::Create(int64_t size)
{
    if(mode != mode_none)
        return false;
    
    if(mem_buffer && mem_buffer_ownership)
        delete[] mem_buffer;
    
    if(!size)
        size = 1;
    mem_buffer_size = size;
    mem_buffer = new uint8_t[mem_buffer_size];
    mem_buffer_ownership = true;
    
    mode = mode_mem_write;
    mem_buffer_pos = 0;
    
    word_buffer_size = 32;

    // byte_buffer_size = 8;
    
    return mode == mode_mem_write;
}
// ********************************************************************************************
bool CBitMemory::Close()
{
    if(mode != mode_mem_write && mode != mode_mem_read)
        return false;
    
    if(mode == mode_mem_write)
    {
        if(word_buffer_pos)
            FlushPartialWordBuffer();
    }
    
    if(mem_buffer && mem_buffer_ownership)
        delete[] mem_buffer;
    
    mem_buffer = NULL;
    mem_buffer_ownership = false;
    mode = mode_none;
    
    return true;
}
// ********************************************************************************************
bool CBitMemory::Complete()
{
    if(mode != mode_mem_write && mode != mode_mem_read)
        return false;
    
    if(mode == mode_mem_write)
    {
        if(word_buffer_pos)
            FlushPartialWordBuffer();
    }
    
    return true;
}
// ********************************************************************************************
bool CBitMemory::SetPos(int64_t pos)
{
    if(mode != mode_file_read && mode != mode_mem_read)
        return false;
    
    if((int64_t) pos > mem_buffer_size)
        return false;
    mem_buffer_pos = pos;
    
    word_buffer_pos = 0;
    word_buffer     = 0;
    
    return true;
}
// ********************************************************************************************
// ********************************************************************************************
CBuffer::CBuffer()
{
	type = buffer_t::none;
	max_size = 0;

	v_size_pos = 0;
	v_data_pos = 0;

	is_function = false;
	is_no_data = false;
}

// ********************************************************************************************
CBuffer::~CBuffer()
{
}

// ********************************************************************************************
void CBuffer::SetMaxSize(uint32_t _max_size, uint32_t _offset)
{
	max_size = _max_size;
#ifndef IGNORE_OFFSET_IN_FIRST_BLOCK
	offset = _offset;
#endif
}

// ********************************************************************************************
void CBuffer::GetBuffer(vector<uint32_t>& _v_size, vector<uint8_t>& _v_data)
{
    
	v_size.shrink_to_fit();
	v_data.shrink_to_fit();

	swap(_v_size, v_size);
	swap(_v_data, v_data);

	v_size.clear();
	v_data.clear();

/*	v_size.shrink_to_fit();
	v_data.shrink_to_fit();

	v_size.reserve(_v_size.size() / 4);
	v_data.reserve(_v_data.size() / 4);*/

	type = buffer_t::none;

#ifndef IGNORE_OFFSET_IN_FIRST_BLOCK
	offset = 0;
#endif
}
// ********************************************************************************************
void CBuffer::SetBuffer(vector<uint32_t>& _v_size, vector<uint8_t>& _v_data)
{
	v_size = move(_v_size);
	v_data = move(_v_data);

	_v_size.clear();
	_v_data.clear();

	v_size_pos = 0;
	v_data_pos = 0;

	is_no_data = v_size.empty();
    
}




// ********************************************************************************************
// ********************************************************************************************
CVariantsBuffer::CVariantsBuffer()
{
	max_size = 0;
	v_data_pos = 0;

}

// ********************************************************************************************
CVariantsBuffer::~CVariantsBuffer()
{
}

// ********************************************************************************************
void CVariantsBuffer::SetMaxSize(uint32_t _max_size, uint32_t _offset)
{
	max_size = _max_size;
#ifndef IGNORE_OFFSET_IN_FIRST_BLOCK
	offset = _offset;
#endif
}

// ********************************************************************************************
void CVariantsBuffer::GetBuffer(vector<uint8_t>& _v_data)
{
    

	v_data.shrink_to_fit();

	swap(_v_data, v_data);

	v_data.clear();

/*	v_size.shrink_to_fit();
	v_data.shrink_to_fit();

	v_size.reserve(_v_size.size() / 4);
	v_data.reserve(_v_data.size() / 4);*/


#ifndef IGNORE_OFFSET_IN_FIRST_BLOCK
	offset = 0;
#endif
}
// ********************************************************************************************
void CVariantsBuffer::SetBuffer(vector<uint8_t>& _v_data)
{

	v_data = move(_v_data);


	_v_data.clear();

	v_data_pos = 0;


}
