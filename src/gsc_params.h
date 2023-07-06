#pragma once

#include <string>
#include <stdlib.h>
#include "defs.h"
#include "bsc.h"
enum class task_mode_t {none, mcompress, mdecompress};
enum class file_type {VCF_File, BCF_File};


const bsc_params_t p_bsc_fixed_fields = { 25, 16, 64,LIBBSC_BLOCKSORTER_BWT, LIBBSC_CODER_QLFC_ADAPTIVE };
	
const bsc_params_t p_bsc_meta = { 25, 16, 64,LIBBSC_BLOCKSORTER_BWT, LIBBSC_CODER_QLFC_ADAPTIVE };

const bsc_params_t p_bsc_size = { 25, 19, 128,LIBBSC_BLOCKSORTER_BWT, LIBBSC_CODER_QLFC_ADAPTIVE };
const bsc_params_t p_bsc_flag = { 25, 19, 64,LIBBSC_BLOCKSORTER_BWT, LIBBSC_CODER_QLFC_ADAPTIVE };
const bsc_params_t p_bsc_text = { 25, 19, 128,LIBBSC_BLOCKSORTER_BWT, LIBBSC_CODER_QLFC_ADAPTIVE };
const bsc_params_t p_bsc_int = { 25, 19, 128,LIBBSC_BLOCKSORTER_BWT, LIBBSC_CODER_QLFC_ADAPTIVE };
const bsc_params_t p_bsc_real = { 25, 19, 64,LIBBSC_BLOCKSORTER_BWT, LIBBSC_CODER_QLFC_ADAPTIVE };


struct GSC_Params{
    
    task_mode_t task_mode;
    file_type in_type, out_type;
    uint32_t max_replication_depth, max_MB_memory;
    std::string in_file_name;
    std::string out_file_name;
    std::string out_samples_file_name;
    std::string samples;
    std::string range;


    uint32_t no_threads;
    uint32_t no_gt_threads;
    uint32_t var_in_block;
    uint32_t ploidy;
    uint64_t vec_len;
    uint32_t records_to_process;
    uint32_t n_samples;

    char compression_level, mode;
    bool MB_memory;
    size_t no_blocks ;

    bool compress_all;
    bool out_AC_AN;
    bool out_genotypes;
    bool out_samples_name;
    bool out_file_flag;
    bool out_ohter_fields;
    uint32_t minAC, maxAC;
    float minAF, maxAF;
    
    GSC_Params()
    {
        task_mode = task_mode_t::none;
        in_type = file_type::VCF_File;
        out_type = file_type::VCF_File;
        max_replication_depth = 100;
        ploidy = 2;
        no_threads = 5;
        no_gt_threads = 1;
        var_in_block = 0;     
        vec_len =  0;
        n_samples = 0;
        no_blocks = 0;
        in_file_name = "";
        out_file_name = "";
        out_samples_file_name = "";
        range = "";
        samples = "";
        MB_memory = true;  // Remember some of the decoded vectors
        max_MB_memory = 0;  // 0 means no limit (if MB_memory == true)
        compression_level = '1';
        compress_all = false;
        out_samples_name = false;
        out_AC_AN = false;
        out_genotypes = true;
        out_file_flag = false;
        out_ohter_fields = false;
        records_to_process = UINT32_MAX;
        mode = '\0';
        minAC = 0;
        maxAC = INT32_MAX;
        minAF = 0;
        maxAF = 1;
    }

};






