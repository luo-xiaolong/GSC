
#pragma once
#include <stdio.h>
#include <iostream>
#include <string>

#include <regex>
#include <sdsl/bit_vectors.hpp>
#include "bit_memory.h"
#include <stdio.h>
#include<algorithm>
#include "gsc_params.h"
#include "decompression_reader.h"

#include "samples.h"

#include "variant.h"

// #include <deque>

#include "nmmintrin.h"

#include "htslib/vcf.h"

#include "htslib/hts.h"

// #include <condition_variable>


#include "defs.h"

#include "utils.h"

#include "vint_code.h"
#include <unordered_map>
#include <array>
#define MMAP

#ifdef MMAP
#include <cpp-mmf/memory_mapped_file.hpp>
#endif
#define MAX 0x7FFFFFFFFFFFFFFF




class Decompressor {

    
    DecompressionReader decompression_reader;

    GSC_Params params;

    vector<string> v_samples;
    string header;
    
    // size_t num_chunks;
    vector<block_t> fixed_variants_chunk,fixed_variants_chunk_io;
	vector<vector<uint32_t>> sort_perm,sort_perm_io;
    vector<uint8_t> decompress_gt_indexes,decompress_gt_indexes_io;
    vector<uint32_t> sparse_matrix_cols;



    Samples smpl;

    // string samples_to_decompress;

    // Range

    string range;
    string in_file_name, out_file_name;
    char compression_level;
    file_type out_type;
    bool out_genotypes;
    bool decompression_mode_type;

    uint32_t start_chunk_id = 0;
    uint32_t end_chunk_id = 0;
    uint32_t cur_chunk_id = 0;
    int64_t range_1,range_2;

    size_t chunk_size;


    // Output and output settings
    htsFile *out_file = nullptr;
    vector<htsFile*> split_files; 

    bcf_hdr_t * out_hdr = nullptr;
    bcf1_t * rec;

    COutFile out_fam;
    COutFile out_bim;
    COutFile out_bed;

    string cur_chrom = "";
    int cur_file = -1;


    uint32_t records_to_process;

    // uint32_t max_MB_memory;

    uint32_t minAC,maxAC;

    float minAF, maxAF;

    float min_qual, max_qual; 

    // vector<vector<int>> s_perm;

    kstring_t str = {0,0,0};

    // vector<block_t> fixed_variants_chunk_io;

    vector<int> pos_to_block;
    
    // vector<pair<std::string,int>> d_where_chrom;

    // vector<int64_t> d_last_pos;

    // uint16_t* decomp_samples_index = nullptr;

    // vector<uint8_t> decomp_no_index;

    // vector<int> index_pos;

    vector<uint32_t> id_pos;

    
    // uint32_t cur_first_row;

    uint32_t cur_vec_id = 0;

    uint32_t end_chunk_actual_pos = 0;

    uint32_t start_chunk_actual_pos = 0;

    uint32_t * sampleIDs = nullptr;
    
    size_t standard_block_size;
    
    uint8_t *decomp_data = nullptr;

    uint8_t *decomp_data_perm = nullptr;

    int full_byte_count,trailing_bits;

    long long *tmp_arr = nullptr;

    long long * lookup_table_ptr = nullptr;

    int fields_pos;
    // string genotype;
    vector<uint8_t> genotype;
    string record;
    variant_desc_t temp_desc;
    vector<field_desc> temp_fields;
    int count = 0;
    vector<vector<field_desc>> all_fields;
    vector<vector<field_desc>> all_fields_io;

    // uint32_t no_fields = 0;
    // vector<key_desc> decomp_keys;
    // size_t no_other_fields;
    // int no_keys;
   
    // bool skip_processing;
    std::unordered_map<uint64_t, uint8_t *> done_unique;

    std::unordered_map<uint64_t, uint8_t *>::const_iterator got_it;

    alignas(8) uint8_t gt_lookup_table[256][256][8];
    // alignas(8) uint8_t lut1[256][256][8];
    uint8_t map_t256[256];
    std::unordered_map<char, std::unordered_map<char, std::array<uint8_t, 2>>> bits_lut;

    // vector<uint8_t> pt;
    string genetic_distance = "0";
    void initialLut();
    // void initialLut2();
    void initialXORLut();
    void getRangeGT(uint8_t *a, size_t no_rec, string &str);
    void getRangeGT(uint8_t *a,const vector<uint32_t> &rev_perm,size_t no_rec, vector<uint8_t> &str);
    void getSamplesGT(string &str, uint8_t res1, uint8_t res2, int p,int cur_blo ,size_t &standard_block_size);

    static bool comp(const variant_desc_t &v1, const variant_desc_t &v2);
    
    static bool comp1(const variant_desc_t &v1, const variant_desc_t &v2);

    bool initalIndex();


    void decoded_vector_row(uint64_t vec_id, uint64_t offset, uint64_t v_offset,uint64_t length, int pos, uint8_t *decomp_data);
    
    int decompressAll();

    int BedFormatDecompress();

    int decompressRange(const string & range);

    int decompressSampleSmart(const string & range);


    uint8_t extract_partial_bytes(uint64_t vec_id, std::vector< std::pair<uint32_t, uint32_t> > & whichByte_whereInRes, uint8_t * resUnique, bool & is_uniqe_id, uint64_t & curr_zeros, uint64_t & curr_copy,uint64_t curr_non_copy_vec_id_offset, uint8_t * resAll, uint32_t unique_pos_first_in_block,  bool full_decode = true); //uint8_t extract_partial_bytes(uint64_t vec_id, std::vector< std::pair<uint32_t, uint8_t> > & whichByte_whereInRes, uint8_t * resUnique, bool & is_uniqe_id, uint64_t & curr_zeros, uint64_t & curr_copy, uint8_t * resAll);

    int decompressRangeSample(const string & range);

    void inline decode_perm_rev(int vec2_start, const vector<uint32_t> &rev_perm, uint8_t *decomp_data_perm, uint8_t *decomp_data);


    void inline reverse_perm(const vector<uint32_t> &perm, vector<uint32_t> &rev_perm, int no_haplotypes);


    uint32_t perm_lut8[8];


    uint8_t *zeros_only_vector = nullptr;


    // Create a deque containing ids of remembered vectors

    std::deque<uint64_t> stored_unique;

    uint64_t max_stored_unique = 0;

    sdsl::rrr_vector<>::rank_1_type rrr_rank_zeros_bit_vector[2];

    sdsl::rrr_vector<>::rank_1_type rrr_rank_copy_bit_vector[2];

    sdsl::bit_vector zeros_bit_vector[2];

    // sdsl::rank_support_v5<> rank_zeros_only_vector[2];

    sdsl::bit_vector copy_bit_vector[2];



    bool SetVariantToRec(variant_desc_t& desc, vector<field_desc>& fields, vector<key_desc> &keys, vector<uint8_t> &_my_str, size_t _standard_block_size);
    void appendVCFToRec(variant_desc_t &_desc, vector<uint8_t> &_genotype, size_t _standard_block_size, vector<field_desc> &_fields, vector<key_desc> &_keys);
    void appendVCF(variant_desc_t &_desc,vector<uint8_t> &_my_str,size_t _no_haplotypes);
    // void appendVCFToRec(variant_desc_t &_desc, string_view _genotype, size_t _no_haplotypes,vector<field_desc>& _fields, vector<key_desc> &_keys);
    // bool SetVariantToRec(variant_desc_t& desc, vector<field_desc>& fields, vector<key_desc> &keys, string_view _my_str, size_t _standard_block_size);
    bool SetVariant(variant_desc_t &desc, vector<uint8_t> &_my_str, size_t _standard_block_size);
    // void decompress();
    bool splitFileWriting(int file_num);
    bool OpenForWriting();
    
    void setMemoryUsage();
    int initOutSplitFile();
    int initOut();
    void WriteBEDMagicNumbers();
    int analyzeInputSamples(vector<string>& v_samples);   
    bool initDecompression(DecompressionReader &decompression_reader);
    bool Close();
    // bool decompress_meta(vector<uint8_t> &compress_v_samples,vector<uint8_t>& compress_v_header);
    bool analyzeInputRange(uint32_t &start_chunk_id,uint32_t &end_chunks);
    void calculate_start_position(int &_start_block,int &_start_position);
    void calculate_end_position(int &_end_block,int &_end_position);


public:

    //htsFile *out;

    Decompressor()

    {
       
        // skip_processing = true;
        // // compression_level = '1';

        // out_type = file_type::VCF_File;

        // out_file_name = "";

        // samples_to_decompress = "";

        // max_MB_memory = 0;

        // MB_memory = true;

        // range = "";

        records_to_process = UINT32_MAX;
        decompression_mode_type = false;
        // out_genotypes = true;
        
        // out_ohter_fields = false;

        // minAC = 0;

        // maxAC = INT32_MAX;

        // minAF = 0;

        // maxAF = 1;

        // pos_hash.reserve(50000000);

    }

    Decompressor(GSC_Params & _params)
    {

        params = _params;
        
        // skip_processing = true;

        in_file_name = params.in_file_name;

        compression_level = params.compression_level;

        out_type = params.out_type;

        out_file_name = params.out_file_name;
        decompression_mode_type = false;
        // // out_samples_name = params.out_samples_name;
        
        // // out_samples_file_name = params.out_samples_file_name;

        // // out_file_flag = params.out_file_flag;

        // // out_ohter_fields = params.out_ohter_fields;

        // // max_block_num = params.max_block_num;

        // samples_to_decompress = params.samples;

        // max_MB_memory = params.max_MB_memory;

        // MB_memory = params.MB_memory;

        range = params.range;

        out_genotypes = params.out_genotypes;

        // records_to_process = params.records_to_process;

        minAF = params.minAF;

        maxAF = params.maxAF;

        minAC = params.minAC;

        maxAC = params.maxAC;

        min_qual = params.min_qual;
        
        max_qual = params.max_qual;

        for (int i = 0; i < 8; ++i)
            perm_lut8[i] = 1 << (7 - i);
        fixed_variants_chunk.reserve(no_variants_in_buf);
	    sort_perm.reserve(no_variants_in_buf);
        decompress_gt_indexes.reserve(no_variants_in_buf);
        fixed_variants_chunk_io.reserve(no_variants_in_buf);
	    sort_perm_io.reserve(no_variants_in_buf);
        decompress_gt_indexes_io.reserve(no_variants_in_buf);
        all_fields.reserve(no_variants_in_buf);
        {
            bits_lut['0']['0'] = {1, 1};
            bits_lut['0']['1'] = {0, 1};
            bits_lut['0']['2'] = {1, 1};
            bits_lut['0']['.'] = {1, 0};
            bits_lut['1']['0'] = {0, 1};
            bits_lut['1']['1'] = {0, 0};
            bits_lut['1']['2'] = {0, 1};
            bits_lut['1']['.'] = {1, 0};
            bits_lut['2']['0'] = {1, 1};
            bits_lut['2']['1'] = {0, 1};
            bits_lut['2']['2'] = {1, 1};
            bits_lut['2']['.'] = {1, 0};
            bits_lut['.']['0'] = {1, 0};
            bits_lut['.']['1'] = {1, 0};
            bits_lut['.']['2'] = {1, 0};
            bits_lut['.']['.'] = {1, 0};
        }

    }

    ~Decompressor()

    {
        // // test_out.close();
        // if(!s_perm.empty())

        //     s_perm.clear();    
        if(decomp_data){
            delete [] decomp_data;
            decomp_data = nullptr;
        }
        if(decomp_data_perm){
            delete [] decomp_data_perm;
            decomp_data_perm = nullptr;
        }
        if(zeros_only_vector)

            delete [] zeros_only_vector;

        if(sampleIDs)

            delete [] sampleIDs;

        
        
        

    }
    // bool getChrom();

    // // bool loadFile();

    // void setMemoryUsage();

    // int initOut();

    // // int loadBCF();

    // int initSample(vector<string>& v_samples);               //  2022年8月3日添加 vector<string>& v_samples

    // void clearHash();//

    // void setVariant(const vector<block_t> &de_q_block,uint32_t &_cur_first_row);//

    // void getPosHash();//

    // void setperm(vector<vector<int>> perm);

    // void getperm(int block_id, uint32_t * perm);

    

    // void setsamples(string str);


    // void setoutname(string str);

	// void setHeader(const string str);

    

    // bool DecompressVarinat();
    bool decompressProcess();
};



