#pragma once

#include <iostream>
#include <vector>
#include "gsc_params.h"
#include "defs.h"
#include <sdsl/bit_vectors.hpp>
#include "bit_memory.h"
#include "file_handle.h"
#include "compression_reader.h"
#include "variant.h"
#include "block_processing.h"
#include "queues.h"
#include <tuple>
#include "bsc.h"
#include "zstd_compress.h"
// #include <filesystem>
using namespace std;

class Compressor
{
    GSC_Params params;

    // CompOtherFields<int,uint8_t,uint8_t> comp_other_fields;
    // CompVarBlockQueue<fixed_field_block> comp_var_block_queue;
    CompVarBlockQueue<fixed_field_block> comp_sort_block_queue;
    vector<thread> part_compress_thread;
    vector<thread> block_process_thread;

    FILE *comp = nullptr;
    bool is_stdout = false;
    FILE *temp_file = nullptr;
    uint64_t sdsl_offset;
    uint64_t other_fields_offset;
    bool mode_type;
    string fname ;
    char *temp_file1_fname = nullptr;
    string temp_file2_fname;
    File_Handle_2 * file_handle2 = nullptr;


    sdsl::bit_vector zeros_only_bit_vector[2];
    sdsl::bit_vector copy_bit_vector[2];
    sdsl::bit_vector unique;
    sdsl::rank_support_v5<> rank_unique;
    sdsl::rank_support_v5<> rank_zeros_only_vector[2];
    sdsl::rank_support_v5<> rank_copy_bit_vector[2];


    uint64_t copy_no = 0;
    uint64_t  unique_no = 0;
    vector<uint32_t> comp_pos_copy;
    vector<bool> all_zeros;
    vector<bool> all_copies;
    uint32_t no_blocks = 0;
    // uint64_t start = 0;
    int64_t prev_pos = 0;
    // uint32_t sort_field_block_id = 0;
    uint32_t fixed_field_block_id = 0;
    // sort_field_block sort_field_block_io,sort_field_block_compress,comp_sort_field_block;
    // fixed_field_block fixed_field_block_compress,fixed_field_block_io;
    fixed_field_block fixed_field_block_io,fixed_field_block_compress;


    
    bool end_of_processing = false;
    uint32_t no_curr_chrom_block =  0;
    vector<int64_t> chunks_min_pos;
    uint32_t cur_chunk_actual_pos = 0;
    map<uint32_t,vector<uint8_t>> vint_last_perm;

    map<int, chunk_stream> chunks_streams;   

    CBitMemory bm_comp_copy_orgl_id;
    uint32_t used_bits_cp;
    vector<pair<std::string, uint32_t>> where_chrom;

    uint64_t no_vec;
    size_t block_size;
    mutex mtx_gt_block;
	condition_variable cv_gt_block;
    int cur_block_id = 0;


    //压缩Meta
    vector<uint8_t> all_v_header, comp_v_header;
	vector<uint8_t> all_v_samples, comp_v_samples;

    //统计压缩后各字段的大小（单位BYTE）
    uint64_t Meta_comp_size = 0;
    uint64_t CHORM_comp_size = 0;
    uint64_t POS_comp_size = 0;
    uint64_t ID_comp_size = 0;
    uint64_t REF_comp_size = 0;
    uint64_t ALT_comp_size = 0;
    uint64_t QUAL_comp_size = 0;
    uint64_t GT_comp_size = 0;

    mutex mtx_v_coder;
	condition_variable cv_v_coder;
	vector<uint32_t> v_coder_part_ids;
    vector<CBSCWrapper*> v_bsc_size;
    vector<CBSCWrapper*> v_bsc_data;

    
    bool input_pos;
    size_t toal_all_size = 0;
    uint32_t no_keys;
    vector<key_desc> keys;
    int key_gt_id;

    

    

    bool OpenForWriting(const string &out_file_name);
    char bits_used(unsigned int n);
    void compressReplicatedRow();;
    bool writeCompressFlie();

    void compress_other_fileds(SPackage& pck, vector<uint8_t>& v_compressed, vector<uint8_t>& v_tmp);
    void compress_INT_fileds(SPackage& pck, vector<uint8_t>& v_compressed, vector<uint8_t>& v_tmp);
    void lock_coder_compressor(SPackage& pck);
    bool check_coder_compressor(SPackage& pck);
    void unlock_coder_compressor(SPackage& pck);
    void lock_gt_block_process(int &_block_id);
    bool check_gt_block_process(int &_block_id);
    void unlock_gt_block_process();
    void Encoder(vector<uint8_t>& v_data, vector<uint8_t>& v_tmp);
    bool compress_meta(vector<string> v_samples,const string& v_header);
    void InitCompressParams();
    bool compressFixedFields(fixed_field_block &field_block_io);
    bool OpenTempFile(const string &out_file_name);
    bool writeTempFlie(fixed_field_block &fixed_field_block_io);
public:
    ~Compressor()
    {
        for (auto p : v_bsc_size)
		    if (p)
			    delete p;

	    for (auto p : v_bsc_data)
		    if (p)
			    delete p;

	    if (file_handle2)
		    delete file_handle2;
        // if(fname)
        //     free(fname);
        if(temp_file1_fname)
            free(temp_file1_fname);
        // if(comp)
        //     fclose(comp);
           
    }

    Compressor(GSC_Params &_params)
    {
        params = _params;
        // curr_vec_id = 0;
        copy_no = 0;
        unique_no = 0;
        no_blocks = 0;
        input_pos = true;
        all_zeros.reserve(no_variants_in_buf);
        all_copies.reserve(no_variants_in_buf);
        comp_pos_copy.reserve(no_variants_in_buf);

        chunks_min_pos.reserve(no_variants_in_buf);
        
    }
    bool CompressProcess();

};


