#pragma once

#include <iostream>
#include "gsc_params.h"
#include "defs.h"
#include "bit_memory.h"
// #include <array>
#include "nmmintrin.h"
#include <immintrin.h>
// #include <random>
#include <map>
#include <cstring>
#include "variant.h"
#include "utils.h"
#include "vint_code.h"
#include <mutex>
#include <condition_variable>
typedef vector<uint64_t> mc_vec_t;        
template <typename T>
inline uint64_t bit_cost(const T &x, const T &y)
{
    uint64_t r = 0;
    int max_i = x.size();
    // cout << "max_i:" << max_i << endl;
    switch (max_i % 4)
    {
    case 3:
        --max_i;
        r += _mm_popcnt_u64(x[max_i] ^ y[max_i]);
    case 2:
        --max_i;
        r += _mm_popcnt_u64(x[max_i] ^ y[max_i]);
    case 1:
        --max_i;
        // cout << "max_i:" << max_i << endl;
        r += _mm_popcnt_u64(x[max_i] ^ y[max_i]);
    }
    for (int i = max_i - 1; i >= 0;)
    {
        r += _mm_popcnt_u64(x[i] ^ y[i]);
        --i;
        r += _mm_popcnt_u64(x[i] ^ y[i]);
        --i;
        r += _mm_popcnt_u64(x[i] ^ y[i]);
        --i;
        r += _mm_popcnt_u64(x[i] ^ y[i]);
        --i;
    }
    return r;
}

#ifdef __AVX2__
template <typename T>
inline uint64_t bit_cost(const T& x, const T& y, uint64_t best_cost, uint64_t* new_xor)
{
    const int n = x.size();
    const int m = n / 4 * 4;
    uint64_t r = 0;

    __m256i* x_ptr = (__m256i*)x.data();
    __m256i* y_ptr = (__m256i*)y.data();
    __m256i* new_xor_ptr = (__m256i*)new_xor;

    for (int i = 0; i < m / 4; ++i)
    {
        const __m256i xx = _mm256_loadu_si256(x_ptr + i);
        const __m256i yy = _mm256_loadu_si256(y_ptr + i);
        const __m256i xor_val = _mm256_xor_si256(xx, yy);
        new_xor_ptr[i] = xor_val;
        __m256i v2 = _mm256_popcnt_epi64(xor_val);
        r += _mm256_reduce_add_epi64(v2);

    }

    for (int i = m; i < n; ++i)
    {
        new_xor[i] = x[i] ^ y[i];
        r += _mm_popcnt_u64(new_xor[i]);
    }

    return r;
}
#else
template <typename T>
inline uint64_t bit_cost(const T &x, const T &y, uint64_t best_cost,uint64_t* new_xor) 
{
    uint64_t r = 0;
    int max_i = x.size();
    // cout<<"x[max_i]:"<< x[0]<<endl;
    switch (max_i % 4)
    {
    case 3:
        --max_i;
        new_xor[max_i]=(x[max_i] ^ y[max_i]);           
        r += _mm_popcnt_u64(new_xor[max_i]);
        
    case 2:
        --max_i;
        new_xor[max_i]=(x[max_i] ^ y[max_i]);        
        r += _mm_popcnt_u64(new_xor[max_i]);
        
    case 1:
        --max_i;
        new_xor[max_i]=(x[max_i] ^ y[max_i]);       
        r += _mm_popcnt_u64(new_xor[max_i]);
        
        // cout << "R:" << r << endl;
    }
    for (int i = max_i - 1; i >= 0 && r < best_cost;)
    {
        // cout << "i:" << i << endl;
        new_xor[i]=x[i] ^ y[i];                         
        r += _mm_popcnt_u64( new_xor[i]);               
        --i;

        new_xor[i]=x[i] ^ y[i];                         
        r += _mm_popcnt_u64( new_xor[i]);
        --i;

        new_xor[i]=x[i] ^ y[i];                        
        r += _mm_popcnt_u64( new_xor[i]);
        --i;

        new_xor[i]=x[i] ^ y[i];                        
        r += _mm_popcnt_u64( new_xor[i]);
        --i;
    }
    return r;
}
#endif
class BlockProcess
{

    GSC_Params params ;
    uint64_t cur_no_vec;
    uint64_t cur_vec_id;
    uint32_t no_copy;
    uint32_t no_samples_index;
    uint8_t *data = nullptr;
    uint32_t perm_lut8[8];
    uint64_t perm_lut64[64];
    uint64_t start = 0;
    uint32_t cur_block_id = 0;
    // int64_t prev_pos;
    void permute_range_vec(uint64_t id_start, uint64_t id_stop, vector<uint32_t> &v_perm,vector<bool> &zeros, vector<bool> &copies, vector<uint32_t> &origin_of_copy, vector<uint8_t> &samples_indexes);
    inline void get_perm(vector<uint32_t> perm, int n,vector<variant_desc_t> &v_vcf_data_compress);
    mutex mtx_v_part1;
	condition_variable cv_v_part1;

public:
    BlockProcess(GSC_Params &_params)
    {
        params = _params;
        no_copy = 0;
        // prev_pos = 0;
        no_samples_index = 0;
        for (int i = 0; i < 8; ++i)
            perm_lut8[i] = 1 << (7 - i);
        for (int i = 0; i < 64; ++i)
            perm_lut64[i] = 1ull << i;
    }
    ~BlockProcess()
    {
      
    }
    void SetCurBlock(uint64_t _cur_no_vec, uint8_t *cur_data);
    void ProcessSquareBlock(vector<uint32_t> &perm,vector<bool> &zeros, vector<bool> &copies, vector<uint32_t> &origin_of_copy, vector<uint8_t> &samples_indexes, bool permute = true);
    void ProcessLastBlock(vector<bool> &zeros, vector<bool> &copies, vector<uint32_t> &origin_of_copy, vector<uint8_t> &samples_indexes);
    void ProcessVariant(vector<uint32_t> &perm,vector<variant_desc_t> &v_vcf_data_io);
    // void ProcessLastPerm(vector<uint32_t> &perm,vector<vector<uint8_t>> &_vint_last_perm);
    void addSortFieldBlock(fixed_field_block &_fixed_field_block_io,vector<bool> &_all_zeros,vector<bool> &_all_copies,vector<uint32_t> &_comp_pos_copy,
    vector<bool> &_zeros_only, vector<bool> &_copies, vector<uint32_t> &_origin_of_copy,vector<uint8_t> &_samples_indexes,vector<variant_desc_t> &_v_vcf_data_io,int64_t &prev_pos);
    // void addSortFieldBlock(sort_field_block &sort_fixed_field_block_io,vector<bool> &_all_zeros,vector<bool> &_all_copies,vector<uint32_t> &_comp_pos_copy,
    // vector<bool> &_zeros_only, vector<bool> &_copies, vector<uint32_t> &_origin_of_copy,vector<uint8_t> &_samples_indexes,FieldsPackage &fields_pck,int64_t &prev_pos);
};

