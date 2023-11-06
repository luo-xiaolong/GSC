
#pragma once

#include <iostream>
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "defs.h"
#include "gsc_params.h"
#include "bit_memory.h"
#include "queues.h"
#include "variant.h"
#include <string>
#include <vector>
#include "utils.h"
#include "file_handle.h"
#include <unordered_set>
#include <unordered_map>
#include <stack>
class CompressionReader {

    
	htsFile* in_file = nullptr;
    vector<htsFile*> merge_files;
    bcf_hdr_t * vcf_hdr= nullptr;
    // vector<bcf_hdr_t*> vcf_hdrs;
	bcf1_t * vcf_record;
    file_type in_type;
    std::string in_file_name;
    bool in_open;
    bool vcf_hdr_read;
    compress_mode_t compress_mode;
    bool merge_flag;
    bool merge_failure_flag;

    uint32_t no_samples;
    uint32_t ploidy;
    uint64_t vec_len;   
    uint64_t no_vec;
    vector<string> samples_list;
    vector<variant_desc_t> v_vcf_data_compress;
    int *cur_g_data = nullptr;
    int ncur_g_data;
    int *gt_data = nullptr;
    int temp;
    int64_t cur_pos;
    int32_t tmpi;
    size_t chunk_size;  

    CBitMemory bv;
    fixed_field_block fixed_field_block_buf;
    int64_t block_max_size;
    uint32_t no_vec_in_block, vec_read_in_block, block_id;
    uint32_t no_fixed_fields,vec_read_fixed_fields,fixed_fields_id;
    GtBlockQueue * Gt_queue = nullptr;
    VarBlockQueue<fixed_field_block>* Var_queue = nullptr;


    void *dst_int = nullptr;
    void *dst_real = nullptr;
    void *dst_str = nullptr;
    void *dst_flag = nullptr;
    int  ndst_int = 0;
    int  ndst_real = 0;
    int  ndst_str = 0;
    int  ndst_flag = 0;
    std::vector<int> FilterIdToFieldId;
    std::vector<int> InfoIdToFieldId;
    std::vector<int> FormatIdToFieldId;
    int no_flt_keys, no_info_keys, no_fmt_keys;
    uint32_t no_keys;
    vector<key_desc> keys;
    int key_gt_id;
    uint32_t no_actual_variants;
    vector<uint32_t> actual_variants;
    vector<uint32_t> v_size;
	vector<uint8_t> v_data;	
    vector<CBuffer> v_o_buf;
    vector<int> v_buf_ids_size;
	vector<int> v_buf_ids_data;
    File_Handle_2 * file_handle2 = nullptr;
    PartQueue<SPackage> * part_queue = nullptr;
    unordered_map<int, unordered_set<int>> field_order_graph;
    bool field_order_flag;
    // unordered_map<int, unordered_set<int>> graph;
    unordered_map<int, int> inDegree;
    vector<int> order;
    size_t no_chrom_num;
    string no_chrom;
    vector<pair<std::string,uint32_t>> where_chrom;
    string cur_chrom;
    vector<int64_t> chunks_min_pos;
    bool start_flag;




    // int temp_count = 0;
    #ifdef LOG_INFO
	unordered_map<int, unordered_set<int>> distinct_values;
    #endif
    bool ReadFile();
    bool setBitVector();
	void addVariant(int * gt_data, int ngt_data,variant_desc_t &desc);
    bool GetVariantFromRec(bcf1_t* rec, vector<field_desc>& fields);
    bool GetFilterInfoFormatKeys(int &no_flt_keys, int &no_info_keys, int &no_fmt_keys, vector<key_desc> &keys);
    void ProcessFixedVariants(bcf1_t *vcf_record, variant_desc_t &desc);
	bool SetVariantOtherFields(vector<field_desc> &fields);
    vector<int> topoSort(unordered_map<int, unordered_set<int>>& graph);
    vector<int> topo_sort(unordered_map<int, unordered_set<int>> &graph,unordered_map<int, int> inDegree);


    
public:

    CompressionReader() {
        in_open = false;
        vcf_hdr_read = false;
        no_samples = 0;
        ploidy = 0;
        no_chrom_num = 0;
        no_vec = 0;
        start_flag = true;
        field_order_flag = false;
    }
    CompressionReader(const GSC_Params & params) {
        in_open = false;
        vcf_hdr_read = false;
        no_samples = 0;
        ploidy = params.ploidy;
        no_chrom_num = 0;
        in_file_name = params.in_file_name;
        in_type = params.in_type;
        compress_mode = params.compress_mode;
        merge_flag = params.merge_file_flag;
        merge_failure_flag = false;
        v_vcf_data_compress.reserve(no_variants_in_buf);
        actual_variants.reserve(no_variants_in_buf);
        no_vec = 0;
        start_flag = true;
        field_order_flag = false;
    }
    
  
    
    ~CompressionReader() {

    if(in_open)
    {
        bcf_hdr_destroy(vcf_hdr);
        vcf_hdr = nullptr;
        hts_close(in_file);
        in_file = nullptr;
        in_open = false;
        // bcf_destroy1(vcf_record);
    }
        // test.close();
    //    var_out.close();
    }
          
    void setQueue(GtBlockQueue * _queue)
    {
        Gt_queue = _queue;
    };
    void setPartQueue(PartQueue<SPackage> *_part_queue){
        part_queue = _part_queue;
    }
    uint64_t getNoVec()
    {
        return no_vec;
    }

    bool OpenForReading(string & file_name);
    uint32_t GetSamples(vector<string> &s_list);
    bool GetHeader(string &v_header);
    void InitVarinats(File_Handle_2 *_file_handle2);
	bool ProcessInVCF();
	uint32_t setNoVecBlock(GSC_Params & params);
	void GetWhereChrom(vector<pair<std::string,uint32_t>> &_where_chrom,vector<int64_t> &chunks_min_pos);
    uint32_t GetOtherFieldsBlockSum();
    void GetOtherField(vector<key_desc> &_keys,uint32_t &_no_keys,int &_key_gt_id);
    vector<uint32_t> GetActualVariants();
    void UpdateKeys(vector<key_desc> &_keys);
    void CloseFiles();
};
