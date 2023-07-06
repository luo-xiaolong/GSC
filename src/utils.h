#pragma once
// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include "defs.h"
#include <mutex>
#include <condition_variable>
#include <nmmintrin.h>
#include <string>
#include <vector>
#include <iostream>
#include "htslib/kstring.h"
#include "htslib/bgzf.h"
#include "htslib/vcf.h"
using namespace std;

#ifdef OUR_STRTOL
#ifdef __APPLE__
long int strtol(const char* str, char** endptr, int base) ;
#else
long int strtol(const char* str, char** endptr, int base) noexcept;
#endif
#endif
struct chunk_stream {
    uint32_t cur_chunk_actual_pos;
    size_t offset;
    chunk_stream() : cur_chunk_actual_pos(0), offset(0)
	{};
	chunk_stream(int64_t _cur_chunk_actual_pos, size_t _offset) : cur_chunk_actual_pos(_cur_chunk_actual_pos), offset(_offset)
	{};
};
typedef struct sortblock {
	int64_t val;
	int p;
	string s_ref;
    sortblock(int64_t a, int b, string str_ref) : val(a), p(b), s_ref(str_ref){}

} sblock;
// *****************************************************************************************

void append_str(vector<uint8_t>& v_comp, const string& x);
void read_str(const vector<uint8_t>& v_comp, size_t& pos, string& x);
void my_merge(vector<sblock> &a, int start, int mid, int end);
void my_merge_sort_r(vector<sblock> &a, int start, int end);
void my_merge_sort(vector<sblock> &a, int len);  
// void insertion_sort(vector<sblock>& arr);
// void merge(vector<sblock>& arr, int l, int m, int r);
// void merge_sort(vector<sblock>& arr, int l, int r);
// void my_merge_sort(vector<sblock>& arr);
template <typename T>
inline void append(vector<uint8_t> &v_comp, T x)
{
    append_str(v_comp, std::to_string(x));
};

template<typename T>
inline void read(vector<uint8_t>& v_comp, size_t& pos, T& x) {
    string t;
    read_str(v_comp, pos, t);
    try{
        x = stoll(t);
    }catch(const std::invalid_argument& e){
        
        std::cerr << "Invalid argument: " << e.what()<<":"<<v_comp.size()<<":"<<pos<<":"<<(t.size())<< std::endl;
    }
};

class MyBarrier {
public:
    MyBarrier(unsigned int count) :
        m_count(count), m_generation(0),
        m_count_reset_value(count)
    {
    }

    void count_down_and_wait()
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        unsigned int gen = m_generation;
        if (--m_count == 0)
        {
            m_generation++;
            m_count = m_count_reset_value;
            m_cond.notify_all();
            return;
        }
        while (gen == m_generation)
            m_cond.wait(lock);
    }

private:
    std::mutex m_mutex;
    std::condition_variable m_cond;
    unsigned int m_count;
    unsigned int m_generation;
    unsigned int m_count_reset_value;
};

// *****************************************************************************************
#define hts_expand00(type_t, n, m, ptr) if ((n) > (m)) { \
int t = (m); (m) = (n); kroundup32(m); \
(ptr) = (type_t*)realloc((ptr), (m) * sizeof(type_t)); \
memset(((type_t*)ptr)+t,0,sizeof(type_t)*((m)-t)); \
}


int bcf_update_genotypes_fast(const bcf_hdr_t *hdr, bcf1_t *line,kstring_t &str);
