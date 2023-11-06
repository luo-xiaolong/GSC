#include "decompressor.h"
#include <bitset>
#include <chrono>

using namespace std::chrono;



bool Decompressor::analyzeInputRange(uint32_t & start_chunk_id,uint32_t & end_chunk_id){

    size_t prev_chrom_size;
    bool query_flag = false;

    // cout << "Begin to extract the genotype sparse matrix" << endl;
    if (range == "")
    {
        
        for (size_t i = 0; i < decompression_reader.d_where_chrom.size(); i++)
        {
            prev_chrom_size = i ? decompression_reader.d_where_chrom[i - 1].second : 0;
            end_chunk_id += (decompression_reader.d_where_chrom[i].second - prev_chrom_size + params.no_blocks-1) / params.no_blocks;
            // no_block += (d_where_chrom[i].second - prev_chr_size) / no_comp_bolck + (((d_where_chrom[i].second - prev_chr_size) % no_comp_bolck) ? 1 : 0);
            
        }

    }
    else
    {
        std::regex pattern(R"(^(\w+)(?::(-?\d+))?(?:,(-?\d+))?:?,?$)");
        std::smatch matches;
        std::string cur_query_chrom;
        if (std::regex_match(range, matches, pattern)) {
            cur_query_chrom = matches[1].str();
            if(matches[2].matched){
                range_1 = std::stoi(matches[2].str());
                if(matches[3].matched)
                    range_2 = std::stoi(matches[3].str());
                else
                    range_2 = MAX;
            }else{
                if(matches[3].matched){
                    std::cout << "Invalid input format." << std::endl;
                    return 0;
                }else{
                    range_1 = 0;
                    range_2 = MAX;
                }
            }
            // range_1 = (matches[2].matched) ? std::stoi(matches[2].str()) : 0 ;
            // range_2 = (matches[3].matched) ? std::stoi(matches[3].str()) : MAX;

            // std::cout << "Chromosome: " << cur_query_chrom << std::endl;
            // std::cout << "Start position: " << range_1 << std::endl;
            // std::cout << "End position: " << range_2 << std::endl;
        } else {
            std::cout << "Invalid input format." << std::endl;

        }

        // string cur_query_chrom = range.substr(0, range.find(':'));
        // auto curr_pos = range.find(':');
        // string_view cur_query_chrom(range.c_str(),curr_pos);
        
        // try {
        //     string cur_range = range.substr(curr_pos + 1);
            
        //     curr_pos = cur_range.find(',');
        //     if(curr_pos == string::npos){
                
        //         if(cur_range == ""){
        //             range_1 = 0;
        //             range_2 = MAX;
        //         }
        //         else{
        //             range_1 = stoll(cur_range);
        //             range_2 = MAX;
                    
        //         }
        //     }
        //     else{
        //         range_1 = stoll(cur_range.substr(0, curr_pos));
        //         cur_range = cur_range.substr(curr_pos + 1);
        //         if (cur_range == "")
        //             range_2 = MAX;    
        //         else
        //         {
        //             range_2 = stoll(cur_range);
        //         }    
                           
        //     }
        // } 
        // catch (const std::invalid_argument& e) {
        //     std::cerr << "Invalid argument: " << e.what() << std::endl;
        // }
        //  catch (const std::out_of_range& e) {
        //     std::cerr << "Out of range: " << e.what() << std::endl;
        // }
        if (range_2 < range_1)
        {
            cout << "error range!" << endl;
            return 0;
        }

        for (size_t i = 0; i < decompression_reader.d_where_chrom.size(); i++)
        {
            prev_chrom_size = i ? decompression_reader.d_where_chrom[i - 1].second : 0;
            
            if (cur_query_chrom == decompression_reader.d_where_chrom[i].first)
            {
                query_flag = true;
                start_chunk_id = end_chunk_id;
                end_chunk_id += (decompression_reader.d_where_chrom[i].second - prev_chrom_size + params.no_blocks-1) / params.no_blocks;

                break;
            }
            end_chunk_id += (decompression_reader.d_where_chrom[i].second - prev_chrom_size + params.no_blocks-1) / params.no_blocks;
            
        }

        if (!query_flag)
        {
            cout << "The specified chromosome wass not found!!!\n";
            return 0;
        }

        int64_t start_chunk = lower_bound(decompression_reader.chunks_min_pos.begin() + start_chunk_id, decompression_reader.chunks_min_pos.begin() + end_chunk_id, range_1) - decompression_reader.chunks_min_pos.begin() - 1;
        int64_t end_chunk = upper_bound(decompression_reader.chunks_min_pos.begin() + start_chunk_id, decompression_reader.chunks_min_pos.begin() + end_chunk_id, range_2) - decompression_reader.chunks_min_pos.begin() - 1;

        if (end_chunk < 0)
        {
            cout << "no var int the range!" << endl;
            return 0;
        }
        start_chunk_id = start_chunk > (int)start_chunk_id ? start_chunk : start_chunk_id;
        end_chunk_id = end_chunk + 1;

        
    }
    return 1;
}
bool Decompressor::initDecompression(DecompressionReader &decompression_reader){

    standard_block_size = decompression_reader.n_samples *(uint32_t)decompression_reader.ploidy;
    max_stored_unique = standard_block_size*2;

    if(standard_block_size < 1024)
    {
        chunk_size = CHUNK_SIZE1;
    }
    else if(standard_block_size < 4096 )
    {
        chunk_size = CHUNK_SIZE2;
    }
    else if(standard_block_size < 8192)
    {
        chunk_size = CHUNK_SIZE3;
    }
    else
    {
        chunk_size = standard_block_size;
    }
    params.no_blocks = chunk_size/standard_block_size;
    rrr_rank_zeros_bit_vector[0] = sdsl::rrr_vector<>::rank_1_type(&decompression_reader.rrr_zeros_bit_vector[0]);
    rrr_rank_zeros_bit_vector[1] = sdsl::rrr_vector<>::rank_1_type(&decompression_reader.rrr_zeros_bit_vector[1]);
    rrr_rank_copy_bit_vector[0] = sdsl::rrr_vector<>::rank_1_type(&decompression_reader.rrr_copy_bit_vector[0]);
    rrr_rank_copy_bit_vector[1] = sdsl::rrr_vector<>::rank_1_type(&decompression_reader.rrr_copy_bit_vector[1]);
    if (params.samples == ""){
        decomp_data = new uint8_t[decompression_reader.vec_len * 2];
        decomp_data_perm = new uint8_t[decompression_reader.vec_len * 2];
        zeros_only_vector = new uint8_t[decompression_reader.vec_len]();
    }
    int no_haplotypes = standard_block_size;
    if (params.samples != ""){
        no_haplotypes = smpl.no_samples * decompression_reader.ploidy;
        uint64_t bv_size = decompression_reader.rrr_zeros_bit_vector[0].size();

        zeros_bit_vector[0] = sdsl::bit_vector(bv_size);
        zeros_bit_vector[1] = sdsl::bit_vector(bv_size);
        copy_bit_vector[0] = sdsl::bit_vector(bv_size);
        copy_bit_vector[1] = sdsl::bit_vector(bv_size);
        uint64_t v_pos;

        for (v_pos = 0; v_pos + 64 < bv_size; v_pos += 64)
        {
            zeros_bit_vector[0].set_int(v_pos, decompression_reader.rrr_zeros_bit_vector[0].get_int(v_pos, 64), 64);
            zeros_bit_vector[1].set_int(v_pos, decompression_reader.rrr_zeros_bit_vector[1].get_int(v_pos, 64), 64);
            copy_bit_vector[0].set_int(v_pos, decompression_reader.rrr_copy_bit_vector[0].get_int(v_pos, 64), 64);
            copy_bit_vector[1].set_int(v_pos, decompression_reader.rrr_copy_bit_vector[1].get_int(v_pos, 64), 64);
        }

        uint64_t tail_len = bv_size - v_pos;
        zeros_bit_vector[0].set_int(v_pos, decompression_reader.rrr_zeros_bit_vector[0].get_int(v_pos, tail_len), tail_len);
        zeros_bit_vector[1].set_int(v_pos, decompression_reader.rrr_zeros_bit_vector[1].get_int(v_pos, tail_len), tail_len);
        copy_bit_vector[0].set_int(v_pos, decompression_reader.rrr_copy_bit_vector[0].get_int(v_pos, tail_len), tail_len);
        copy_bit_vector[1].set_int(v_pos, decompression_reader.rrr_copy_bit_vector[1].get_int(v_pos, tail_len), tail_len);
    } 
    else
    {
        if(standard_block_size & 7)//%8)
        {
            full_byte_count = decompression_reader.vec_len - 1;
            trailing_bits = standard_block_size & 7;//%8;
        }
        else
        {
            full_byte_count = decompression_reader.vec_len;
            trailing_bits = 0;
        }
        tmp_arr = new long long[full_byte_count];  
        
    }
    initialXORLut();
    initialLut();
    if(out_type != file_type::BED_File){
        int fmt_id = bcf_hdr_id2int(out_hdr,BCF_DT_ID,"GT");
        bcf_enc_int1(&str, fmt_id);
        bcf_enc_size(&str, decompression_reader.ploidy, BCF_BT_INT8);
        char *tmp;
        str.m = str.l + no_haplotypes + 1;
            
        kroundup32(str.m);

        if ((tmp = (char*)realloc(str.s, str.m)))
            str.s = tmp;
        else
            exit(8);
        
        str.l = 3;
    }
    
    return true;
}
//*************************************************************************************************************************************
// Official Decompress Program Entry
bool Decompressor::decompressProcess()
{
    MyBarrier  my_barrier(3);
    
    decompression_reader.SetNoThreads(params.no_threads);
    // unique_ptr<CompressedFileLoading> cfile(new CompressedFileLoading());
    if(!decompression_reader.OpenReading(in_file_name))
        return false;
    if(params.compress_mode == compress_mode_t::lossless_mode){
        if(!decompression_reader.OpenReadingPart2(in_file_name))
            return false;
    }
    
    decompression_reader.decompress_meta(v_samples, header);
    
    if(analyzeInputSamples(v_samples)) // Retrieving sample name.
        return false;


    if(params.split_flag){

        if(!splitFileWriting(static_cast<int>(decompression_reader.d_where_chrom.size())))
            return false;
        initOutSplitFile();
    }else{
        if(!OpenForWriting())
            return false;
        initOut();
    }

    initDecompression(decompression_reader);
    
    if(!analyzeInputRange(start_chunk_id,end_chunk_id)){
        return false;
    }

    if(!decompression_reader.setStartChunk(start_chunk_id)){
        return false;
    } 
    cur_chunk_id = start_chunk_id;
    unique_ptr<thread> decompress_thread(new thread([&]{
        
        while(cur_chunk_id < end_chunk_id){
            my_barrier.count_down_and_wait();
            my_barrier.count_down_and_wait();
            initalIndex();
            if(out_type == file_type::BED_File){
                BedFormatDecompress();
                
            }else{
                if(params.compress_mode == compress_mode_t::lossless_mode){
                    
                    uint32_t no_actual_varians =  decompression_reader.actual_varians[cur_chunk_id-1];

                    decompressAll();

                    for(uint32_t i = 0; i < no_actual_varians; ++i)
                        for(size_t j = 0; j < decompression_reader.keys.size(); ++j){
                            if(all_fields_io[i][j].data_size)
                            {
                                delete[] all_fields_io[i][j].data;
                                all_fields_io[i][j].data = nullptr;
                                all_fields_io[i][j].data_size = 0;
                            }
                            else
                                all_fields_io[i][j].data = nullptr;

                        }
                    all_fields_io.clear();
                    
                }
                else{
                    if (params.samples == "")
                        decompressRange(range);
                    else
                        decompressSampleSmart(range);
                }
            }
            fixed_variants_chunk_io.clear();
            sort_perm_io.clear();
            decompress_gt_indexes_io.clear();
            
        }
    }));

    //用智能指针创建线程
    unique_ptr<thread> process_thread(new thread([&]{
        while (cur_chunk_id < end_chunk_id)
        {
            if(!decompression_reader.readFixedFields()){
                
                break;
            }
            if(params.compress_mode == compress_mode_t::lossless_mode){

                uint32_t no_actual_varians =  decompression_reader.actual_varians[cur_chunk_id];
                
                while (no_actual_varians--)
                {
                    all_fields.emplace_back(vector<field_desc>(decompression_reader.keys.size()));
                    decompression_reader.GetVariants(all_fields.back());
                    
                }
            }
            decompression_reader.Decoder(fixed_variants_chunk,sort_perm,decompress_gt_indexes,cur_chunk_id);
            
            my_barrier.count_down_and_wait();
            my_barrier.count_down_and_wait();
            
        }
        
    })); 
    while (cur_chunk_id < end_chunk_id)
	{
        my_barrier.count_down_and_wait();
		swap(fixed_variants_chunk, fixed_variants_chunk_io);
        swap(sort_perm, sort_perm_io);
        swap(decompress_gt_indexes, decompress_gt_indexes_io);
        if(params.compress_mode == compress_mode_t::lossless_mode){
            swap(all_fields, all_fields_io);
        }
        start_chunk_actual_pos = decompression_reader.getActualPos(cur_chunk_id);
        end_chunk_actual_pos = decompression_reader.getActualPos(cur_chunk_id+1);

        cur_chunk_id++;
		my_barrier.count_down_and_wait();
          
    }

    
    process_thread->join();
    decompress_thread->join();
    if(params.compress_mode == compress_mode_t::lossless_mode)
        decompression_reader.close();
    Close();
    return true;
}
bool Decompressor::Close(){

    if(out_type == file_type::BED_File){
        out_fam.Close();
        out_bed.Close();
        out_bim.Close();
    }else{
        if(out_hdr){
            bcf_hdr_destroy(out_hdr);
            out_hdr = nullptr;
        }
        if (out_file)
        {
            hts_close(out_file);
            out_file = nullptr;
        }
        if(params.split_flag){
            if(split_files[cur_file]){
                hts_close(split_files[cur_file]);
                split_files[cur_file] = nullptr;            
            }
        }

        if(rec){
            bcf_destroy1(rec);
        }
        if(str.m)
            free(str.s);
        if(tmp_arr)
            delete[] tmp_arr;
    }
    return true;
}
// *****************************************************************************************************************
// bool Decompressor::createSamplesNameFile(vector<string> &v_samples)
// {
//     std::ofstream sn_file(out_samples_file_name + ".sn");
//     if (sn_file)
//     {
//         for (uint32_t i = 0; i < v_samples.size(); i++)
//         {
//             sn_file << v_samples[i] << std::endl;
//         }
//         std::cout << "File with list of samples (" << out_samples_file_name + ".sn"
//                   << ") created." << std::endl;

//         sn_file.close();
//         return true;
//     }
//     else
//     {
//         std::cout << "Could not open " << out_samples_file_name + ".sn"
//                   << "file with list of samples." << std::endl;
//         return false;
//     }
//     return true;
// }

// 建立异或表
// *****************************************************************************************************************
void Decompressor::initialXORLut(){

    uint8_t cur_xor_result = 0;
    uint8_t temp = 0;
    for (int i = 0; i < 256; ++i){
        cur_xor_result = i & perm_lut8[0];
        for(int j = 7; j > 0; --j){
            temp = ((i >> j) ^ (i >> (j-1))) & 1;
            cur_xor_result += temp * perm_lut8[8 - j];
        }
        map_t256[cur_xor_result] = i;
    }

}
// *****************************************************************************************************************
// void Decompressor::initialLut1()
// {
    
//     uint8_t mask;
//     for (int c = 0; c < 8; c++)
//     {
//         mask = 0x80 >> c;
//         for (int i = 0; i < 256; i++)
//         {
//             if (i & mask)
//             {
//                 for (int j = 0; j < 256; j++)
//                 {
//                     if (j & mask) // 11
//                     {
//                         gt_lookup_table[i][j][c] = '0';
//                     }
//                     else // 10
//                     {
//                         gt_lookup_table[i][j][c] = '.';
//                     }
//                 }
//             }
//             else
//             {
//                 for (int j = 0; j < 256; j++)
//                 {
//                     if (j & mask) // 01
//                     {
//                         gt_lookup_table[i][j][c] = '1';
//                     }
//                     else // 00
//                     {
//                         gt_lookup_table[i][j][c] = '0';
//                     }
//                 }
//             }
//         }
//     }
// }
// *****************************************************************************************************************
void Decompressor::initialLut()
{
    // xor_map_table();
    uint8_t mask;
    for (int c = 0; c < 8; c++)
    {
        mask = 0x80 >> c;
        for (int i = 0; i < 256; i++)
        {
            if (i & mask)
            {
                for (int j = 0; j < 256; j++)
                {
                    if (j & mask) // 11
                    {
                        gt_lookup_table[i][j][c] = '2';
                    }
                    else // 10
                    {
                        gt_lookup_table[i][j][c] = '.';
                    }
                }
            }
            else
            {
                for (int j = 0; j < 256; j++)
                {
                    if (j & mask) // 01
                    {
                        gt_lookup_table[i][j][c] = '1';
                    }
                    else // 00
                    {
                        gt_lookup_table[i][j][c] = '0';
                    }
                }
            }
        }
    }
}
void Decompressor::getRangeGT(uint8_t *a,const vector<uint32_t> &rev_perm,size_t no_rec, vector<uint8_t> &str){

    for(size_t i = 0;i<no_rec;i++){
        int cur_byte = rev_perm[i] >> 3;
        int cur_pos = rev_perm[i] % 8;
        int r1 = a[cur_byte];
        int r2 = a[decompression_reader.vec_len + cur_byte];
        str[i]= gt_lookup_table[r1][r2][cur_pos];
    }   
    
}
// // *****************************************************************************************************************
void Decompressor::getRangeGT(uint8_t *a, size_t no_rec, string &str)
{
    // int end,g;
    // if(no_rec & 7)//%8)
    // {
    //     end = decompression_reader.vec_len - 1;
    //     g = no_rec & 7;//%8;
    // }
    // else
    // {
    //     end = decompression_reader.vec_len;
    //     g = 0;
    // }
    // long long * lookup_table_ptr = nullptr;
    // long long *tmp_arr = new long long[end];
    // int vec1_start = 0,vec2_start=decompression_reader.vec_len;
    // for (vec1_start = 0; vec1_start < end; ++vec1_start)
    // {
    //     //memcpy(pt + (vec1_start << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]], 8);
    //     lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
    //     //*(tmp_arr+ vec1_start) = *lookup_table_ptr;
    //     tmp_arr[vec1_start] = *lookup_table_ptr;
    // }
    // memcpy(str.data(), tmp_arr, end << 3);
    // if(g)
    // {
    //     memcpy(str.data() + (end << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], g);
        
    // }   
    // if(tmp_arr)
    //     delete [] tmp_arr;  
}
// void Decompressor::getSamplesGT(vector<uint8_t> &str, uint8_t res1, uint8_t res2, int p)
// {

//     str += gt_lookup_table[res1][res2][p];

// }
// *****************************************************************************************************************
inline bool Decompressor::comp(const variant_desc_t &v1, const variant_desc_t &v2)
{
    return v1.pos <= v2.pos;
}
// *****************************************************************************************************************
inline bool Decompressor::comp1(const variant_desc_t &v1, const variant_desc_t &v2)
{
    return v1.pos < v2.pos;
}

// *****************************************************************************************************************
bool Decompressor::initalIndex()
{
    id_pos.clear();
    id_pos.reserve(no_variants_in_buf);
    for(size_t i = 0; i < decompress_gt_indexes_io.size(); ++i)
    {
        if(decompress_gt_indexes_io[i] == '\0')
            id_pos.emplace_back(i);
    }
    id_pos.shrink_to_fit();
    return true;
}

// // *****************************************************************************************************************
void Decompressor::appendVCF(variant_desc_t &_desc, vector<uint8_t> &_my_str, size_t _no_haplotypes)
{

    bcf_clear(rec);
    string record;
    record = _desc.chrom + "\t0\t" + _desc.id + "\t" + _desc.ref + "\t" + _desc.alt + "\t" + _desc.qual + "\t" + _desc.filter + "\t" + _desc.info;
    kstring_t s;
    s.s = (char *)record.c_str();
    s.m = record.length();
    s.l = 0;
    vcf_parse(&s, out_hdr, rec);
    rec->pos = (int32_t)(_desc.pos - 1);
  
    if (out_genotypes)
    {
        if(count <= INT8_MAX && count > INT8_MIN+1){
            vector<int8_t> gt_arr(_no_haplotypes);

            int t = 0;
            for (auto cur_genotype : _my_str)
            {
                if (cur_genotype == '.')
                    gt_arr[t++] = bcf_gt_missing;
                else
                    gt_arr[t++] = bcf_gt_phased(int(cur_genotype - '0'));

            }
            
            str.l = 3;
            memcpy(str.s + str.l, gt_arr.data(), _no_haplotypes * sizeof(int8_t));
            str.l += _no_haplotypes;
            str.s[str.l] = 0;
            // for(int i=0;i<str.l;i++)
            //     cout<<(int)str.s[i]<;
    
            // bcf_update_genotypes(out_hdr, rec, gt_arr.data(), _no_haplotypes);
            bcf_update_genotypes_fast(out_hdr, rec,str);
        }
        else{

            vector<int> gt_arr(_no_haplotypes);

            int t = 0;
            for (auto cur_genotype : _my_str)
            {
                if (cur_genotype == '.')
                    gt_arr[t++] = bcf_gt_missing;
                else
                    gt_arr[t++] = bcf_gt_phased(int(cur_genotype - '0'));

            }
            bcf_update_genotypes(out_hdr, rec, gt_arr.data(), _no_haplotypes);
        }
    }
    if(params.split_flag){

        if(_desc.chrom != cur_chrom){
            if(cur_file != -1 ){
                if (split_files[cur_file])
                {
                    hts_close(split_files[cur_file]);
                    split_files[cur_file] = nullptr;
                }
            }
            cur_chrom = _desc.chrom;
            cur_file++;
            
        }
        
        bcf_write1(split_files[cur_file], out_hdr, rec);
    }else
        bcf_write1(out_file, out_hdr, rec);
    // bcf_write(out_file, out_hdr, rec);

}
// // *****************************************************************************************************************

void Decompressor::appendVCFToRec(variant_desc_t &_desc, vector<uint8_t> &_genotype, size_t _standard_block_size, vector<field_desc> &_fields, vector<key_desc> &_keys)
{
    bcf_clear(rec);
    record = _desc.chrom + "\t0\t" + _desc.id + "\t" + _desc.ref + "\t" + _desc.alt + "\t" + _desc.qual + "\t" + "." + "\t" + ".";
    kstring_t s;
    s.s = (char *)record.c_str();
    s.m = record.length();
    s.l = 0;
    vcf_parse(&s, out_hdr, rec);
    rec->pos = (int32_t)(_desc.pos - 1);
    int curr_size = 0;
    for (size_t i = 0; i < _fields.size(); i++)
    {
        uint32_t id = _keys[i].actual_field_id;
        if (_keys[id].keys_type == key_type_t::flt)
        {
            
            if (_fields[id].present)
                bcf_add_filter(out_hdr, rec, _keys[id].key_id);
        }
    }
    
    for (size_t i = 0; i < _keys.size(); i++)
    {
        uint32_t id = _keys[i].actual_field_id;
        if (_keys[id].keys_type == key_type_t::info)
        {
            
            if (_fields[id].present)
                switch (_keys[id].type)
		        {
		        case BCF_HT_INT:
                    curr_size = _fields[id].data_size >> 2;
                    // cout<<curr_size<<endl;
                    bcf_update_info_int32(out_hdr,rec,bcf_hdr_int2id(out_hdr, BCF_DT_ID, _keys[id].key_id), _fields[id].data, curr_size);
			        break;
		        case BCF_HT_REAL:
                    curr_size = _fields[id].data_size >> 2;
                    bcf_update_info_float(out_hdr,rec,bcf_hdr_int2id(out_hdr, BCF_DT_ID, _keys[id].key_id), _fields[id].data, curr_size);
			        break;
		        case BCF_HT_STR:
                    curr_size = _fields[id].data_size;
                    bcf_update_info(out_hdr, rec, bcf_hdr_int2id(out_hdr, BCF_DT_ID, _keys[id].key_id), _fields[id].data, curr_size, BCF_HT_STR);
			        break;
		        case BCF_HT_FLAG:
                    curr_size = 1;
                    bcf_update_info_flag(out_hdr,rec,bcf_hdr_int2id(out_hdr, BCF_DT_ID, _keys[id].key_id), _fields[id].data, curr_size);
			        break;
		        }

                
        }
    }

    // // FORMAT

    // field_desc gt_phased;
    for (size_t i = 0; i < _keys.size(); i++)
    {
        int id = _keys[i].actual_field_id;

        if(id == decompression_reader.key_gt_id)
        {
        //     // gt_phased.resize(_fields[id].data_size);
        // //    memcpy(gt_phased.data(), _fields[id].data, _fields[id].data_size);
        //     // gt_phased = move(_fields[id]);
            bool GT_NULL_flag = false;
            for(uint32_t i = 0; i < decompression_reader.n_samples; i++)
                if(_genotype[i] == '.'){
                    GT_NULL_flag = true;
                    break;
                }      
            if(count <= INT8_MAX && count > INT8_MIN+1&&!GT_NULL_flag){
                
                vector<uint8_t> gt_arr(_standard_block_size);
                uint32_t t = 0;
                int cur_gt = 0;
                for(uint32_t i = 0; i < decompression_reader.n_samples; i++){
                    cur_gt = i*decompression_reader.ploidy;
                    if(_genotype[cur_gt] == '.')
                        gt_arr[cur_gt] = bcf_gt_missing;
                    else 
                        gt_arr[cur_gt] = bcf_gt_unphased((int)(_genotype[cur_gt]-'0'));
                    for(uint32_t j = 1;j< (uint32_t)decompression_reader.ploidy;j++){
                        cur_gt++;
                        if(_fields[id].data[t] == '/'){
                            if(_genotype[cur_gt] == '.')
                                    gt_arr[cur_gt] = bcf_gt_missing;
                                else 
                                    gt_arr[cur_gt] = bcf_gt_unphased((int)(_genotype[cur_gt]-'0'));

                        } 
                        else if(_fields[id].data[t] == '|'){
                            if(_genotype[cur_gt] == '.')
                                gt_arr[cur_gt] = bcf_next_gt_missing;
                            else 
                                gt_arr[cur_gt] = bcf_gt_phased((int)(_genotype[cur_gt]-'0'));
                        }
                        t++;

                    }
                }
                
                str.l = 3;
                
                memcpy(str.s + str.l, gt_arr.data(), _standard_block_size * sizeof(uint8_t));
                str.l += _standard_block_size;
                str.s[str.l] = 0;

                bcf_update_genotypes_fast(out_hdr, rec,str);
            }
            else{
                
                vector<int> gt_arr(_standard_block_size);
                uint32_t t = 0;
                int cur_gt = 0;
                for(uint32_t i = 0; i < decompression_reader.n_samples; i++){
                    cur_gt = i*decompression_reader.ploidy;
                    if(_genotype[cur_gt] == '.')
                        gt_arr[cur_gt] = bcf_gt_missing;
                    else 
                        gt_arr[cur_gt] = bcf_gt_unphased(_genotype[cur_gt]-'0');
                    for(uint32_t j = 1;j< (uint32_t)decompression_reader.ploidy;j++){
                        cur_gt++;
                        if(_fields[id].data[t] == '/'){
                        if(_genotype[cur_gt] == '.')
                                gt_arr[cur_gt] = bcf_gt_missing;
                            else 
                                gt_arr[cur_gt] = bcf_gt_unphased(_genotype[cur_gt]-'0');

                        } 
                        else if(_fields[id].data[t] == '|'){
                            if(_genotype[cur_gt] == '.')
                                gt_arr[cur_gt] = bcf_next_gt_missing;
                            else 
                                gt_arr[cur_gt] = bcf_gt_phased(_genotype[cur_gt]-'0');
                        }
                        else{
                            gt_arr[cur_gt] =  GT_NOT_CALL;  
                        }
                        t++;
                    }
                }
                bcf_update_genotypes(out_hdr, rec, gt_arr.data(), _standard_block_size);
            }
            continue;
        }
        if (_keys[id].keys_type == key_type_t::fmt)
        {
            
            if (_fields[id].present)
            {
                switch (_keys[id].type)
		        {
		        case BCF_HT_INT:
                    curr_size = _fields[id].data_size >> 2;
                    bcf_update_format(out_hdr, rec, bcf_hdr_int2id(out_hdr, BCF_DT_ID, _keys[id].key_id), _fields[id].data, curr_size, BCF_HT_INT);
			        break;
		        case BCF_HT_REAL:
                    curr_size = _fields[id].data_size >> 2;
                    bcf_update_format(out_hdr, rec, bcf_hdr_int2id(out_hdr, BCF_DT_ID, _keys[id].key_id), _fields[id].data, curr_size, BCF_HT_REAL);
			        break;
		        case BCF_HT_STR:
                    curr_size = _fields[id].data_size;
                    bcf_update_format(out_hdr, rec, bcf_hdr_int2id(out_hdr, BCF_DT_ID, _keys[id].key_id), _fields[id].data, curr_size, BCF_HT_STR);
			        break;
		        case BCF_HT_FLAG:
                    assert(0);
			        break;
		        }

                
            }
        }
    }
    if(params.split_flag){

        if(_desc.chrom != cur_chrom){
            if(cur_file != -1 ){
                if (split_files[cur_file])
                {
                    hts_close(split_files[cur_file]);
                    split_files[cur_file] = nullptr;
                }
            }
            cur_chrom = _desc.chrom;
            cur_file++;
            
        }
        
        bcf_write1(split_files[cur_file], out_hdr, rec);
    }else
        bcf_write1(out_file, out_hdr, rec);
}
// 
bool Decompressor::SetVariantToRec(variant_desc_t &desc, vector<field_desc> &fields, vector<key_desc> &keys, vector<uint8_t> &_my_str, size_t _standard_block_size)
{
    
    if(desc.alt.find("<M>") == string::npos){

        if(desc.alt.find("<N>") == string::npos){
            if(count){
                
                appendVCFToRec(temp_desc, genotype, _standard_block_size, temp_fields, keys);  
            }
            
            appendVCFToRec(desc, _my_str, _standard_block_size, fields, keys);
            count = 0;
        
        }
        else{
            if(count){
                appendVCFToRec(temp_desc, genotype, _standard_block_size, temp_fields, keys);  
            }
            
            genotype = _my_str;
            // temp_fields.clear();
            // cout<<"start move"<<endl;
            temp_fields.resize(fields.size());
            // temp_fields = std::move(fields);
            for(size_t i = 0;i<fields.size();i++)
                temp_fields[i] = std::move(fields[i]);
            // cout<<"end move"<<endl;
            temp_desc = desc;
            temp_desc.alt = desc.alt.substr(0, desc.alt.find_first_of(','));
            count = 1;

        }
        fields_pos++;
    }
    else{
        if(count == 0){
            appendVCFToRec(desc, _my_str, _standard_block_size, fields, keys);
            count = 0;
            fields_pos++;
        }else{
            temp_desc.alt.append(",");
            temp_desc.alt += desc.alt.substr(0, desc.alt.find_first_of(','));
            uint8_t target = uint8_t('2' + count-1);
            for (size_t i = 0; i < _my_str.size(); i++)
            {
                if (genotype[i] == target && _my_str[i] == '2')
                    genotype[i]++;
            }
            count++;
        }
       
    }
    
    return true;
}
// // // **********************************************************************************************************************************************************

bool Decompressor::SetVariant(variant_desc_t &desc, vector<uint8_t> &_my_str, size_t _standard_block_size)
{
    if(desc.alt.find("<M>") == string::npos){
        if(desc.alt.find("<N>") == string::npos){
            if(count){
                appendVCF(temp_desc, genotype, _standard_block_size);
                
            }
            appendVCF(desc, _my_str, _standard_block_size);
            count = 0;
        }
        else{
            if(count)
                appendVCF(temp_desc, genotype, _standard_block_size);
            
            genotype = _my_str;
            temp_desc = desc;
            temp_desc.alt = desc.alt.substr(0, desc.alt.find_first_of(','));
            count = 1;
        }
    }
    else{
        if(count == 0){
            appendVCF(desc, _my_str, _standard_block_size);
            count = 0;
        }else{
            temp_desc.alt.append(",");
            temp_desc.alt += desc.alt.substr(0, desc.alt.find_first_of(','));
            char target = char('2' + count-1);
            for (size_t i = 0; i < _my_str.size(); i++)
            {

                if (genotype[i] == target && _my_str[i] == '2')
                    genotype[i]++;
            }
            count++;
        }
    }
    
    return true;
}
//*****************************************************************************************************************
int Decompressor::BedFormatDecompress(){

    done_unique.clear();
    stored_unique.clear();
    uint32_t cur_block_id = 0;
    uint32_t c_out_line = 0;
    uint32_t no_var = 0;
    uint32_t start_var = 0;
    int vec1_start,vec2_start = decompression_reader.vec_len;
    // fields_pos  = 0;
    vector<uint8_t> my_str(standard_block_size);
    vector<uint32_t> rev_perm(standard_block_size);
    
    uint64_t curr_non_copy_vec_id_offset = start_chunk_actual_pos * 2 - rrr_rank_zeros_bit_vector[0](start_chunk_actual_pos) - 
                rrr_rank_zeros_bit_vector[1](start_chunk_actual_pos) -rrr_rank_copy_bit_vector[0](start_chunk_actual_pos) - 
                rrr_rank_copy_bit_vector[1](start_chunk_actual_pos);

    no_var = end_chunk_actual_pos - start_chunk_actual_pos;

    size_t cur_var;

    for (cur_var = start_var; cur_var + standard_block_size <= no_var; cur_var += standard_block_size )
    {
        
        cur_block_id = cur_var / standard_block_size;
        for(size_t i = 0; i < standard_block_size; i++){
            
            vec2_start = decompression_reader.vec_len;
            fill_n(decomp_data, decompression_reader.vec_len*2, 0);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, 0, decomp_data_perm);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, vec2_start, decomp_data_perm);

            decode_perm_rev(vec2_start, sort_perm_io[cur_block_id], decomp_data_perm, decomp_data);
            
            
            for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
            {
                lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                tmp_arr[vec1_start] = *lookup_table_ptr;
            }
            memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
            if(trailing_bits)
            {
                memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);
                
            }   



            // 优化1: 避免重复搜索
            // size_t comma_pos = desc.alt.find(",");
            // std::string desc_alt = (comma_pos == std::string::npos) ? desc.alt : desc.alt.substr(0, comma_pos);


            variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[i];
            string desc_pos = to_string(desc.pos);
            size_t comma_pos = desc.alt.find(",");
            string desc_alt = comma_pos == string::npos ? desc.alt : desc.alt.substr(0, comma_pos);
            string bim_line = desc.chrom + "\t" + desc.id + "\t" + genetic_distance + "\t" + desc_pos + "\t" + desc_alt + "\t" + desc.ref + "\n";
            
            out_bim.Write(bim_line.c_str(), bim_line.size());
            for(int i = 0; i < (int)standard_block_size; i += 2){
                char first = my_str[i];
                char second = my_str[i + 1];
                out_bed.PutBit(bits_lut[first][second][0]);
                out_bed.PutBit(bits_lut[first][second][1]);
                // if(first == '0' &&  second == '0'){
                    
                //     out_bed.PutBit(1);
                //     out_bed.PutBit(1);
                // }
                // else if((first == '0' &&  second == '1') || (first == '1' &&  second == '0')){
                //     out_bed.PutBit(0);
                //     out_bed.PutBit(1);
                // }
                // else if(first == '1' &&  second == '1'){
                //     out_bed.PutBit(0);
                //     out_bed.PutBit(0);
                // }
                // else if(first == '.' ||  second == '.'){
                    
                //     out_bed.PutBit(1);
                    
                //     out_bed.PutBit(0);
                // }
            } 
            out_bed.FlushPartialByteBuffer();
          
            
        
        }
    }
    if(no_var % standard_block_size)
    {
        cur_block_id = cur_var / standard_block_size;            
        reverse_perm(sort_perm_io[cur_block_id], rev_perm, standard_block_size);
        
        for(;cur_var < no_var;++cur_var){
            vec2_start = decompression_reader.vec_len;
            fill_n(decomp_data, decompression_reader.vec_len*2, 0);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, 0, decomp_data_perm);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, vec2_start, decomp_data_perm);
            // if(!(no_var%standard_block_size))
            //     decode_perm_rev(vec2_start, sort_perm_io[cur_block_id], decomp_data_perm, decomp_data);
            // else
            //     memcpy(decomp_data, decomp_data_perm, decompression_reader.vec_len*2);

            // for(int i = 0;i<rev_perm.size();i++)
            //     cout<<rev_perm[i]<<" ";
            // cout<<endl;
            decode_perm_rev(vec2_start, rev_perm, decomp_data_perm, decomp_data);
            for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
            {
                lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                tmp_arr[vec1_start] = *lookup_table_ptr;
            }
            memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
            if(trailing_bits)
            {
                memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);
                
            } 
                
            c_out_line = cur_var % standard_block_size;

            variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[c_out_line];

            string desc_pos = to_string(desc.pos);
            string desc_alt = desc.alt.find(",") == string::npos ? desc.alt : desc.alt.substr(0,desc.alt.find(","));
            string bim_line = desc.chrom + "\t" + desc.id + "\t" + genetic_distance + "\t" + desc_pos + "\t" + desc_alt + "\t" + desc.ref + "\n";
            out_bim.Write(bim_line.c_str(), bim_line.size());
            for(int i = 0; i < (int)standard_block_size; i += 2){
                char first = my_str[i];
                char second = my_str[i + 1];
                out_bed.PutBit(bits_lut[first][second][0]);
                out_bed.PutBit(bits_lut[first][second][1]);
                // if(first == '0' &&  second == '0'){
                //     out_bed.PutBit(1);
                //     out_bed.PutBit(1);
                // }
                // else if(first == '0' &&  second == '1' || first == '1' &&  second == '0'){
                //     out_bed.PutBit(0);
                //     out_bed.PutBit(1);
                // }
                // else if(first == '1' &&  second == '1'){
                //     out_bed.PutBit(0);
                //     out_bed.PutBit(0);
                // }
                // else if(first == '.' ||  second == '.'){
                //     out_bed.PutBit(1);
                //     out_bed.PutBit(0);
                // }
            } 
            out_bed.FlushPartialByteBuffer();
            // SetVariantToRec(desc, all_fields_io[fields_pos], decompression_reader.keys, my_str, standard_block_size);
        }
    }
    // if(cur_chunk_id == end_chunk_id && count){
      
    //     appendVCFToRec(temp_desc, genotype, static_cast<uint32_t>(standard_block_size), temp_fields, decompression_reader.keys);
    // }

    cout<< cur_chunk_id << "\r";
    fflush(stdout);

    for (auto &it : done_unique)
        delete[] it.second;

    done_unique.clear();

    return 0;
    
}

//*****************************************************************************************************************
int Decompressor::decompressAll(){

    done_unique.clear();
    stored_unique.clear();
    uint32_t cur_block_id = 0;
    uint32_t c_out_line = 0;
    uint32_t no_var = 0;
    uint32_t start_var = 0;
    int vec1_start,vec2_start=decompression_reader.vec_len;
    fields_pos  = 0;
    vector<uint8_t> my_str(standard_block_size);
    vector<uint32_t> rev_perm(standard_block_size);
    uint64_t curr_non_copy_vec_id_offset = start_chunk_actual_pos * 2 - rrr_rank_zeros_bit_vector[0](start_chunk_actual_pos) - 
                rrr_rank_zeros_bit_vector[1](start_chunk_actual_pos) -rrr_rank_copy_bit_vector[0](start_chunk_actual_pos) - 
                rrr_rank_copy_bit_vector[1](start_chunk_actual_pos);

    no_var = end_chunk_actual_pos - start_chunk_actual_pos;

    size_t cur_var;

    for (cur_var = start_var; cur_var + standard_block_size <= no_var; cur_var += standard_block_size )
    {
        
        cur_block_id = cur_var / standard_block_size;
        for(size_t i = 0; i < standard_block_size; i++){
            
            vec2_start = decompression_reader.vec_len;
            fill_n(decomp_data, decompression_reader.vec_len*2, 0);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, 0, decomp_data_perm);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, vec2_start, decomp_data_perm);

            decode_perm_rev(vec2_start, sort_perm_io[cur_block_id], decomp_data_perm, decomp_data);
            
            for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
            {
                lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                tmp_arr[vec1_start] = *lookup_table_ptr;
            }
            memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
            if(trailing_bits)
            {
                memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);
                
            }   
            variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[i];

            SetVariantToRec(desc, all_fields_io[fields_pos], decompression_reader.keys, my_str,standard_block_size);
            
        
        }
    }
    if(no_var % standard_block_size)
    {
        cur_block_id = cur_var / standard_block_size;            
        reverse_perm(sort_perm_io[cur_block_id], rev_perm, standard_block_size);
        
        for(;cur_var < no_var;++cur_var){
            vec2_start = decompression_reader.vec_len;
            fill_n(decomp_data, decompression_reader.vec_len*2, 0);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, 0, decomp_data_perm);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, vec2_start, decomp_data_perm);
            // if(!(no_var%standard_block_size))
            //     decode_perm_rev(vec2_start, sort_perm_io[cur_block_id], decomp_data_perm, decomp_data);
            // else
            //     memcpy(decomp_data, decomp_data_perm, decompression_reader.vec_len*2);

            // for(int i = 0;i<rev_perm.size();i++)
            //     cout<<rev_perm[i]<<" ";
            // cout<<endl;
            decode_perm_rev(vec2_start, rev_perm, decomp_data_perm, decomp_data);
            for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
            {
                lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                tmp_arr[vec1_start] = *lookup_table_ptr;
            }
            memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
            if(trailing_bits)
            {
                memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);
                
            } 
                
            c_out_line = cur_var % standard_block_size;

            variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[c_out_line];
            
            SetVariantToRec(desc, all_fields_io[fields_pos], decompression_reader.keys, my_str, standard_block_size);
        }
    }
    if(cur_chunk_id == end_chunk_id && count){
      
        appendVCFToRec(temp_desc, genotype, static_cast<uint32_t>(standard_block_size), temp_fields, decompression_reader.keys);
    }

    cout<< cur_chunk_id << "\r";
    fflush(stdout);

    for (auto &it : done_unique)
        delete[] it.second;

    done_unique.clear();

    return 0;
    
}
// // *****************************************************************************************************************
// Decompress by range
int Decompressor::decompressRange(const string &range)
{

    done_unique.clear();
    stored_unique.clear();
    uint32_t cur_block_id = 0;
    uint32_t c_out_line = 0;
    uint32_t no_var = 0;
    uint32_t start_var = 0;
    int vec1_start,vec2_start;
    bool skip_processing = false;
    vector<uint32_t> rev_perm(standard_block_size);
    vector<uint8_t> my_str(standard_block_size);

    uint64_t curr_non_copy_vec_id_offset = start_chunk_actual_pos * 2 - rrr_rank_zeros_bit_vector[0](start_chunk_actual_pos) - 
                rrr_rank_zeros_bit_vector[1](start_chunk_actual_pos) -rrr_rank_copy_bit_vector[0](start_chunk_actual_pos) - 
                rrr_rank_copy_bit_vector[1](start_chunk_actual_pos);
    
    if (range == "")
    {
        no_var = end_chunk_actual_pos - start_chunk_actual_pos;

        size_t cur_var;
        for (cur_var = start_var; cur_var + standard_block_size <= no_var; cur_var += standard_block_size )
        {
            cur_block_id = cur_var / standard_block_size;
            // reverse_perm(sort_perm_io[cur_block_id], rev_perm, standard_block_size);
            for(size_t i = 0; i < standard_block_size; i++){
                vec2_start = decompression_reader.vec_len;
                fill_n(decomp_data, decompression_reader.vec_len*2, 0);
                decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, 0, decomp_data_perm);
                decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, vec2_start, decomp_data_perm);
                decode_perm_rev(vec2_start, sort_perm_io[cur_block_id], decomp_data_perm, decomp_data);
            
                for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
                {
                    lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);

                    tmp_arr[vec1_start] = *lookup_table_ptr;
                }
                memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
                if(trailing_bits)
                {
                    memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);
                    
                }
            
                // string my_str = "";
                // getRangeGT(decomp_data, standard_block_size, my_str);
                // for(int i = 0;i<(int)decompression_reader.vec_len*2;i++)
                //     decomp_data_perm[i] = map_t256[decomp_data_perm[i]];
                // getRangeGT(decomp_data_perm,rev_perm,standard_block_size, my_str);    
                // c_out_line = i;

                variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[i];
                
                if(params.out_AC_AN){
                    
                    uint32_t AN = standard_block_size;
                    uint32_t AC = std::count_if(my_str.begin(), my_str.end(),
                                        [](char c)
                                        { return c == '1' || c == '2' || c == '.'; });
                    
                    skip_processing = ((float)AC / AN > maxAF || (float)AC / AN < minAF || AC > maxAC || AC < minAC);

                    if (skip_processing)
                        continue;

                    desc.info = "AN=" + std::to_string(AN) + ";AC=" + std::to_string(AC);

                    
                } 
                skip_processing = (atoi(desc.qual.c_str()) > max_qual || atoi(desc.qual.c_str()) < min_qual || !(params.out_id == desc.id || params.out_id.empty()));

                if (skip_processing)
                    continue;  

                SetVariant(desc, my_str, standard_block_size);
            
            }
        }
        if(no_var % standard_block_size)
        {
            cur_block_id = cur_var / standard_block_size;
                           
            reverse_perm(sort_perm_io[cur_block_id], rev_perm, standard_block_size);
            for(;cur_var < no_var;++cur_var){
                vec2_start = decompression_reader.vec_len;
                fill_n(decomp_data, decompression_reader.vec_len*2, 0);
                decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, 0, decomp_data_perm);
                decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, vec2_start, decomp_data_perm);
                // if(!(no_var%standard_block_size))
                //     decode_perm_rev(vec2_start, sort_perm_io[cur_block_id], decomp_data_perm, decomp_data);
                // else
                //     memcpy(decomp_data, decomp_data_perm, decompression_reader.vec_len*2);

                decode_perm_rev(vec2_start, rev_perm, decomp_data_perm, decomp_data);
                for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
                {
                    lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                    tmp_arr[vec1_start] = *lookup_table_ptr;
                }
                memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
                if(trailing_bits)
                {
                    memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);
                    
                } 
                // if(!(no_var%standard_block_size))
                //     for(int i = 0;i< (int)decompression_reader.vec_len*2;i++)
                //         decomp_data_perm[i] = map_t256[decomp_data_perm[i]];
                        
                // // string my_str = "";
                // // vector<uint8_t> my_str(standard_block_size);
                // // getRangeGT(decomp_data_perm, standard_block_size, my_str);
                // getRangeGT(decomp_data_perm,rev_perm,standard_block_size, my_str);
                    
                c_out_line = cur_var % standard_block_size;
                variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[c_out_line];
                if(params.out_AC_AN){
                    
                    uint32_t AN = standard_block_size;
                    uint32_t AC = std::count_if(my_str.begin(), my_str.end(),
                                        [](char c)
                                        { return c == '1' || c == '2' || c == '.'; });
                    
                    skip_processing = ((float)AC / AN > maxAF || (float)AC / AN < minAF || AC > maxAC || AC < minAC);

                    if (skip_processing)
                        continue;

                    desc.info = "AN=" + std::to_string(AN) + ";AC=" + std::to_string(AC);

                    
                }   
                
                skip_processing = (atoi(desc.qual.c_str()) > max_qual || atoi(desc.qual.c_str()) < min_qual || !(params.out_id == desc.id || params.out_id.empty()));

                if (skip_processing)
                    continue;
                // for(int i = 0; i < standard_block_size; i++){
                //     cout << my_str[i] << " ";
                // }
                // cout << endl; 
                SetVariant(desc, my_str, standard_block_size);
            }
        }

        if(cur_chunk_id == end_chunk_id && count)
            appendVCF(temp_desc, genotype, standard_block_size);
        for (auto &it : done_unique)
            delete[] it.second;

        done_unique.clear();
    }
    else
    {
       
        int start_block, end_block;
        int start_position, end_position;
        
        if(cur_chunk_id-1 == start_chunk_id){
            
            calculate_start_position(start_block,start_position);
            start_var = uint32_t(start_block * standard_block_size + start_position);
            
        }else{
            start_var = 0;
        }
        if(cur_chunk_id == end_chunk_id){
            
            calculate_end_position(end_block,end_position);
            
            no_var = uint32_t(end_block * standard_block_size + end_position);
        }else{

            no_var = end_chunk_actual_pos - start_chunk_actual_pos;
        }       

        cur_vec_id = (start_chunk_actual_pos + start_var) * 2;
        bool  tail_flag = true;
        for (size_t cur_var = start_var;cur_var < no_var; cur_var++)
        {

            cur_block_id = cur_var / standard_block_size;
            if(fixed_variants_chunk_io[cur_block_id].data_compress.size() != standard_block_size && tail_flag){
                reverse_perm(sort_perm_io[cur_block_id], rev_perm, standard_block_size);
                sort_perm_io[cur_block_id] = rev_perm;
                tail_flag = false;
            }
            vec2_start = decompression_reader.vec_len;
            fill_n(decomp_data, decompression_reader.vec_len*2, 0);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, 0, decomp_data_perm);
            decoded_vector_row(cur_vec_id++, 0, curr_non_copy_vec_id_offset, decompression_reader.vec_len, vec2_start, decomp_data_perm);
            // if(fixed_variants_chunk_io[cur_block_id].data_compress.size() == standard_block_size)
            //     
            // else
            //     memcpy(decomp_data, decomp_data_perm, decompression_reader.vec_len*2);
            decode_perm_rev(vec2_start, sort_perm_io[cur_block_id], decomp_data_perm, decomp_data);
            for (vec1_start = 0; vec1_start < full_byte_count; ++vec1_start)
            {
                lookup_table_ptr = (long long *)(gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start++]]);
                tmp_arr[vec1_start] = *lookup_table_ptr;
            }
            memcpy(my_str.data(), tmp_arr, full_byte_count << 3);
            if(trailing_bits)
            {
                memcpy(my_str.data() + (full_byte_count << 3), gt_lookup_table[decomp_data[vec1_start]][decomp_data[vec2_start]], trailing_bits);
                
            } 
                
            // if(fixed_variants_chunk_io[cur_block_id].data_compress.size() == standard_block_size){
            //     for(int i = 0;i< (int)decompression_reader.vec_len*2;i++)
            //         decomp_data_perm[i] = map_t256[decomp_data_perm[i]];
            // }
            // getRangeGT(decomp_data_perm,rev_perm,standard_block_size, my_str);

            // // string my_str = "";

            // // getRangeGT(decomp_data, standard_block_size, my_str);

            c_out_line = cur_var % standard_block_size;
            variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[c_out_line];
            if(params.out_AC_AN){
                uint32_t AN = standard_block_size;

                uint32_t AC = std::count_if(my_str.begin(), my_str.end(),
                                    [](char c)
                                    { return c == '1' || c == '2' || c == '.'; });
                skip_processing = ((float)AC / AN > maxAF || (float)AC / AN < minAF || AC > maxAC || AC < minAC);
                if (skip_processing)
                    continue;
                desc.info = "AN=" + std::to_string(AN) + ";AC=" + std::to_string(AC);
            }

            skip_processing = (atoi(desc.qual.c_str()) > max_qual || atoi(desc.qual.c_str()) < min_qual || !(params.out_id == desc.id || params.out_id.empty()));

            if (skip_processing)
                continue;
            SetVariant(desc, my_str, standard_block_size);
            
        }
        if(cur_chunk_id == end_chunk_id && count)
            appendVCF(temp_desc, genotype, standard_block_size);
        for (auto &it : done_unique)
            delete[] it.second;
        done_unique.clear();
    }
    cout<< cur_chunk_id << "\r";
    fflush(stdout);
    // if (decomp_data)
    //     delete[] decomp_data;
    // if (decomp_data_perm)
    //     delete[] decomp_data_perm;
    // if (new_decomp_data_perm)
    //     delete[] new_decomp_data_perm;

    return 0;
}
// // *****************************************************************************************************************
int Decompressor::decompressSampleSmart(const string &range)
{
    
    
    uint32_t no_var;

    size_t cur_block_id;

    size_t c_out_line;

    if (smpl.no_samples > NO_SAMPLE_THRESHOLD) // || range != "")
    {
        return decompressRangeSample(range);
    }
    else if (range != "")
    {
        ;
    }
    uint32_t no_haplotypes = smpl.no_samples * decompression_reader.ploidy;

    std::vector<std::pair<uint32_t, uint32_t>> whichByte_whereInRes(no_haplotypes + 1); // +1 for guard

    bool is_unique = false;
    
    bool skip_processing = false;

    uint32_t unique_pos = 0, unique_pos_first_in_block = 0;

    uint64_t curr_zeros = 0, curr_copy = 0;

    vector<uint32_t> rev_perm(standard_block_size);

    vector<uint8_t> gt_variant_data(no_haplotypes);

    uint32_t  prev_block_id = 0xFFFF;

    uint32_t where; // uchar where;

    uint32_t ind_id_orig;

    unique_pos = 0;

    curr_zeros = 0, curr_copy = 0;

    uint32_t written_records = 0;

    uint8_t *resUnique = nullptr;

    resUnique = new uint8_t[standard_block_size*2 * no_haplotypes];

    uint8_t *resAll = nullptr;

    resAll = new uint8_t[2 * no_haplotypes];

    uint32_t first_vec_in_block = 0;

    // int var_block_num=fixed_variants_chunk_io.size();
    uint64_t prev_chr_zeros_copy = rrr_rank_zeros_bit_vector[0](start_chunk_actual_pos) + rrr_rank_zeros_bit_vector[1](start_chunk_actual_pos) +
                                       rrr_rank_copy_bit_vector[0](start_chunk_actual_pos) + rrr_rank_copy_bit_vector[1](start_chunk_actual_pos);
    uint64_t curr_non_copy_vec_id_offset = start_chunk_actual_pos*2 - prev_chr_zeros_copy;
    
    if (range != "")
    {

        int start_block, end_block;
        int start_position, end_position;
        uint32_t start_var;
        if(cur_chunk_id-1 == start_chunk_id){
            calculate_start_position(start_block,start_position);
            start_var = uint32_t(start_block * standard_block_size + start_position);
            
        }else{
            start_var = 0;
        }
        if(cur_chunk_id == end_chunk_id){
            calculate_end_position(end_block,end_position);
            no_var = uint32_t(end_block * standard_block_size + end_position);
        }else{
            no_var = end_chunk_actual_pos - start_chunk_actual_pos;
        }
        
        cur_vec_id = (start_var + start_chunk_actual_pos) * 2;
  
        for (size_t cur_var = start_var;cur_var < no_var ; cur_var++)
        {
            cur_block_id = cur_var / standard_block_size;

            if (cur_block_id != prev_block_id) // Get perm and find out which bytes need decoding
            {
                first_vec_in_block = start_chunk_actual_pos*2 + cur_block_id * standard_block_size*2;
                unique_pos = 0;
                if (prev_block_id == 0xFFFF) // to set curr_zeros, curr_copy (first processed block)
                {
                    uint8_t parity = first_vec_in_block & 1;
                    uint32_t id = first_vec_in_block >> 1;
                    curr_zeros = rrr_rank_zeros_bit_vector[0](id + ((parity))) + rrr_rank_zeros_bit_vector[1](id);
                    curr_copy = rrr_rank_copy_bit_vector[0](id + ((parity))) + rrr_rank_copy_bit_vector[1](id);
                }
                unique_pos_first_in_block = first_vec_in_block - curr_zeros - curr_copy - curr_non_copy_vec_id_offset;
                if(fixed_variants_chunk_io[cur_block_id].data_compress.size() != standard_block_size)
                    rev_perm = sort_perm_io[cur_block_id];
                else
                    reverse_perm(sort_perm_io[cur_block_id], rev_perm, standard_block_size);

            //                 for(int i = 0;i<rev_perm.size();i++)
            //     cout<<rev_perm[i]<<" ";
            // cout<<endl;
                for (uint32_t s = 0; s < smpl.no_samples; s++)
                {
                    for (uint32_t p = 0; p < decompression_reader.ploidy; p++)
                    {
                        ind_id_orig = sampleIDs[s] * decompression_reader.ploidy + p;
                        where = s * decompression_reader.ploidy + p;
                        // cout<<perm[ind_id_orig]<<endl;
                        whichByte_whereInRes[where] = std::make_pair(rev_perm[ind_id_orig] >> 3, where); // now original_id_only
                    }
                }
                whichByte_whereInRes[no_haplotypes] = std::make_pair(0xFFFFFFFF, no_haplotypes); // guard
                // Sort by byte_no
                std::sort(whichByte_whereInRes.begin(), whichByte_whereInRes.end());
                prev_block_id = cur_block_id;

                // Get vectors from all block
                // copies only within the same block, so only part of resUnique (within the same block and using the same perm) will be used - no need to fix previous unique bytes (these are in different blocks, with diff perm)
                
                for (uint64_t i = first_vec_in_block; i < cur_vec_id; i++)
                {

                    extract_partial_bytes(i, whichByte_whereInRes, resUnique, is_unique, curr_zeros, curr_copy, curr_non_copy_vec_id_offset, resAll, unique_pos_first_in_block, false);
                    if (is_unique)
                    {
                        memcpy(resUnique + unique_pos * no_haplotypes, resAll + (i & 1) * no_haplotypes, no_haplotypes);
                        unique_pos++;
                    }
                }
            }
        //     //   else
            {

        //         // previous vectors in block are decompressed already
        
                for (uint64_t i = cur_vec_id; i <= cur_vec_id + 1; i++)
                {
                    extract_partial_bytes(i, whichByte_whereInRes, resUnique, is_unique, curr_zeros, curr_copy,curr_non_copy_vec_id_offset, resAll, unique_pos_first_in_block);

                    if (is_unique)
                    {
                        
                        //  for(uint32_t idx = 0; idx < no_haplotypes; idx++)
                        {
                            memcpy(resUnique + unique_pos * no_haplotypes, resAll + (i & 1) * no_haplotypes, no_haplotypes);
                        }
                        unique_pos++;
      
                    }
                }
            }
            written_records++;
            prev_block_id = cur_block_id;
            cur_vec_id += 2;
            int r_i = 0;
            // string my_str;

            size_t c_out_line = cur_var % standard_block_size;

            // if (fixed_variants_chunk_io[cur_block_id].data_compress.size() != standard_block_size){
            //     for (int g = 0; g < (int)smpl.no_samples; g++)
            //     {
            //         int samples_start = sampleIDs[g] * decompression_reader.ploidy;
            //         for (size_t p = 0; p < decompression_reader.ploidy; p++, r_i++)
            //         {
            //             uint32_t gt_pos = rev_perm[samples_start+p] % 8;
            //             gt_variant_data[r_i] = gt_lookup_table[resAll[r_i]][resAll[r_i + no_haplotypes]][gt_pos];
            //             // getSamplesGT(gt_variant_data[r_i], resAll[r_i], resAll[r_i + no_haplotypes],gt_pos,standard_block_size); // 2022年9月15日新增

            //         }
            //     }
            // }
            // else{
            for (int g = 0; g < (int)smpl.no_samples; g++)
            {
                int samples_start = sampleIDs[g] * decompression_reader.ploidy;
                for (size_t p = 0; p < decompression_reader.ploidy; p++, r_i++)
                {
                    uint32_t gt_pos = rev_perm[samples_start+ p] % 8;
                    resAll[r_i] = map_t256[resAll[r_i]];
                    resAll[r_i + no_haplotypes] = map_t256[resAll[r_i + no_haplotypes]];
                    gt_variant_data[r_i] = gt_lookup_table[resAll[r_i]][resAll[r_i + no_haplotypes]][gt_pos];
                    // getSamplesGT(gt_variant_data[r_i], resAll[r_i], resAll[r_i + no_haplotypes],gt_pos, standard_block_size); // 2022年9月15日新增

                }
            }
            // }
            variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[c_out_line];
            if(params.out_AC_AN){
                uint32_t AN = no_haplotypes;

                uint32_t AC = std::count_if(gt_variant_data.begin(), gt_variant_data.end(),
                                    [](char c)
                                    { return c == '1' || c == '2' || c == '.'; });
                skip_processing = ((float)AC / AN > maxAF || (float)AC / AN < minAF || AC > maxAC || AC < minAC);
                if (skip_processing)
                    continue;

                desc.info = "AN=" + std::to_string(AN) + ";AC=" + std::to_string(AC);
            }
            skip_processing = (atoi(desc.qual.c_str()) > max_qual || atoi(desc.qual.c_str()) < min_qual || !(params.out_id == desc.id || params.out_id.empty()));

            if (skip_processing)
                continue;          
            SetVariant(desc, gt_variant_data, static_cast<size_t> (no_haplotypes));
            
        }
        if(cur_chunk_id == end_chunk_id && count)
            appendVCF(temp_desc, genotype, static_cast<size_t> (no_haplotypes));

    }
    else
    {
        no_var = end_chunk_actual_pos - start_chunk_actual_pos;

        for (size_t cur_var = 0; cur_var < no_var; cur_var++)
        {
        
            // Permutations
            cur_block_id = cur_var  / standard_block_size;

            if (cur_block_id != prev_block_id) // get perm and find out which bytes need decoding
            {
                first_vec_in_block = start_chunk_actual_pos*2 +cur_block_id * standard_block_size*2;
                unique_pos = 0;
                if (prev_block_id == 0xFFFF) // to set curr_zeros, curr_copy (first processed block)
                {
                    uint8_t parity = first_vec_in_block & 1;
                    uint32_t id = first_vec_in_block >> 1;
                    curr_zeros = rrr_rank_zeros_bit_vector[0](id + ((parity))) + rrr_rank_zeros_bit_vector[1](id);
                    curr_copy = rrr_rank_copy_bit_vector[0](id + ((parity))) + rrr_rank_copy_bit_vector[1](id);
                }

                unique_pos_first_in_block = first_vec_in_block - curr_zeros - curr_copy - curr_non_copy_vec_id_offset;
                if(fixed_variants_chunk_io[cur_block_id].data_compress.size() != standard_block_size)
                    rev_perm = sort_perm_io[cur_block_id];
                else
                    reverse_perm(sort_perm_io[cur_block_id], rev_perm, standard_block_size);

                for (uint32_t s = 0; s < smpl.no_samples; s++)
                {
                    for (uint32_t p = 0; p < decompression_reader.ploidy; p++)
                    {
                        ind_id_orig = sampleIDs[s] * decompression_reader.ploidy + p;
                        where = s * decompression_reader.ploidy + p;

                        whichByte_whereInRes[where] = std::make_pair(rev_perm[ind_id_orig] >> 3, where); // now original_id_only
                    }
                }
                whichByte_whereInRes[no_haplotypes] = std::make_pair(0xFFFFFFFF, no_haplotypes); // guard

                std::sort(whichByte_whereInRes.begin(), whichByte_whereInRes.end());
                // cout<<"prev_block_id:"<<prev_block_id<<endl;
                // cout<<"start_chunk_actual_pos:"<<start_chunk_actual_pos<<endl;
                prev_block_id = cur_block_id;
                // Get vectors from all block
                
                // for (uint64_t i = cur_vec_id; i <= cur_vec_id + 1; i++)
                // {
                //     extract_partial_bytes(i, whichByte_whereInRes, resUnique, is_unique, curr_zeros, curr_copy,curr_non_copy_vec_id_offset, resAll, unique_pos_first_in_block);

                //     if (is_unique)
                //     {
                //         memcpy(resUnique + unique_pos * no_haplotypes, resAll + (i & 1) * no_haplotypes, no_haplotypes);
                //         unique_pos++;
                //     }
                // }
            }
            // else
            {
                for (uint64_t i = cur_vec_id; i <= cur_vec_id + 1; i++)
                {

                    extract_partial_bytes(i, whichByte_whereInRes, resUnique, is_unique, curr_zeros, curr_copy,curr_non_copy_vec_id_offset, resAll, unique_pos_first_in_block);

                    if (is_unique)
                    {
                        memcpy(resUnique + unique_pos * no_haplotypes, resAll + (i & 1) * no_haplotypes, no_haplotypes);
                        unique_pos++;
                        
                    }
                }
            }
            cur_vec_id += 2;
            prev_block_id = cur_block_id;
            c_out_line = cur_var % standard_block_size;
            int r_i = 0;
            // if (fixed_variants_chunk_io[cur_block_id].data_compress.size() != standard_block_size){
            //     for (int g = 0; g < (int)smpl.no_samples; g++)
            //     {
            //         int samples_start = sampleIDs[g] * decompression_reader.ploidy;
            //         for (size_t p = 0; p < decompression_reader.ploidy; p++, r_i++)
            //         {
            //             uint32_t gt_pos = rev_perm[samples_start+ p] % 8;
            //             gt_variant_data[r_i] = gt_lookup_table[resAll[r_i]][resAll[r_i + no_haplotypes]][gt_pos];
            //             // getSamplesGT(gt_variant_data[r_i], resAll[r_i], resAll[r_i + no_haplotypes],gt_pos,standard_block_size); // 2022年9月15日新增

            //         }
            //     }
            // }
            // else{
            for (int g = 0; g < (int)smpl.no_samples; g++)
            {
                int samples_start = sampleIDs[g] * decompression_reader.ploidy;
                for (size_t p = 0; p < decompression_reader.ploidy; p++, r_i++)
                {
                    uint32_t gt_pos = rev_perm[samples_start+ p] % 8;
                    resAll[r_i] = map_t256[resAll[r_i]];
                    resAll[r_i + no_haplotypes] = map_t256[resAll[r_i + no_haplotypes]];
                    gt_variant_data[r_i] = gt_lookup_table[resAll[r_i]][resAll[r_i + no_haplotypes]][gt_pos];
                    // getSamplesGT(gt_variant_data[r_i], resAll[r_i], resAll[r_i + no_haplotypes],gt_pos, standard_block_size); // 2022年9月15日新增

                }
            }
            // }

            variant_desc_t desc = fixed_variants_chunk_io[cur_block_id].data_compress[c_out_line];
            
            if(params.out_AC_AN){

                uint32_t AN = no_haplotypes;

                uint32_t AC = std::count_if(gt_variant_data.begin(), gt_variant_data.end(),
                                    [](char c)
                                    { return c == '1' || c == '2' || c == '.'; });
                skip_processing = ((float)AC / AN > maxAF || (float)AC / AN < minAF || AC > maxAC || AC < minAC);
                if (skip_processing)
                    continue;

                desc.info = "AN=" + std::to_string(AN) + ";AC=" + std::to_string(AC);
            }
            skip_processing = (atoi(desc.qual.c_str()) > max_qual || atoi(desc.qual.c_str()) < min_qual || !(params.out_id == desc.id || params.out_id.empty()));

            if (skip_processing)
                continue;         
            SetVariant(desc, gt_variant_data, static_cast<size_t> (no_haplotypes));

        }
        if(cur_chunk_id == end_chunk_id && count)
            appendVCF(temp_desc, genotype, static_cast<size_t> (no_haplotypes));
    }
    cout<< cur_chunk_id << "\r";
    fflush(stdout);
    if (resUnique)
        delete[] resUnique;
    if (resAll)
        delete[] resAll;

    return 0;
}
int Decompressor::decompressRangeSample(const string &range)
{
    return 0;
}
// // ***************************************************************************************************************************************************************
// // Get the byte corresponding to the sample genotype in the matrix
uint8_t Decompressor::extract_partial_bytes(uint64_t vec_id, std::vector<std::pair<uint32_t, uint32_t>> &whichByte_whereInRes, uint8_t *resUnique, bool &is_uniqe_id, uint64_t &curr_zeros, uint64_t &curr_copy, uint64_t curr_non_copy_vec_id_offset, uint8_t *resAll, uint32_t unique_pos_first_in_block, bool full_decode) // uint8_t Decompressor::extract_partial_bytes(uint64_t vec_id, std::vector< std::pair<uint32_t, uint8_t> > & whichByte_whereInRes, uint8_t * resUnique, bool & is_uniqe_id, uint64_t & curr_zeros, uint64_t & curr_copy, uint8_t * resAll)
{
    uint32_t next_haplotype = 0;
    uint32_t last_byte = whichByte_whereInRes.back().first; // Last byte to decode
    uint32_t tmp = 0;
    uint8_t parity = vec_id & 1;
    uint64_t vector = vec_id >> 1, curr_non_copy_vec_id;
    uint64_t no_haplotypes = whichByte_whereInRes.size() - 1;
    uint64_t vec_start = parity * no_haplotypes; //(vec_id - first_vec_in_block) * no_haplotypes;
    // size_t index_pos_size;
    size_t id_pos_size;
    uint64_t index_start = 0;
    size_t cur_i = 0;
    size_t count = 0;
    // uint8_t prev_index1 = 0, prev_index2 = 0;
    uint32_t cur_index;
    size_t prev_i = 999999;
    uint32_t cur_id = 0;
    int cur_pos = 0;
    std::vector<uint8_t> curr_gt_indexes;
    if (copy_bit_vector[parity][vector]) // If vector is a copy od other vector (certainly it is placed within the same block)
    {
        if (!full_decode)
        {
            is_uniqe_id = false;
            curr_copy++;
            return 0;
        }
        unsigned long long bit_pos = (curr_copy)*decompression_reader.used_bits_cp;
        decompression_reader.bm_comp_copy_orgl_id.SetPos(bit_pos >> 3);
        decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, bit_pos & 7); // %8
        decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, decompression_reader.used_bits_cp);

        // Here curr_non_copy_vec_id is a fake curr_non_copy_vec_id (it is ud of the next non_copy_vec_id)

        curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy - curr_non_copy_vec_id_offset;

        curr_non_copy_vec_id = curr_non_copy_vec_id - tmp - 1;
        // cout<<"unique_pos_first_in_block:"<<curr_non_copy_vec_id<<" "<<unique_pos_first_in_block<<endl;
        is_uniqe_id = false;
        curr_copy++;
        memcpy(resAll + vec_start, resUnique + (curr_non_copy_vec_id - unique_pos_first_in_block) * no_haplotypes, no_haplotypes);
        

        return 0;
    }
    else if (zeros_bit_vector[parity][vector]) 
    {
        is_uniqe_id = false;
        curr_zeros++;
        if (!full_decode)
            return 0;
        fill_n(resAll + vec_start, no_haplotypes, 0);
        return 0;
    }
    else
    {
        curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy - curr_non_copy_vec_id_offset;
        is_uniqe_id = true;
    }
    fill_n(resAll + no_haplotypes*parity, no_haplotypes, 0);

    if (curr_non_copy_vec_id == 0)
        index_start = 0;
    else
        index_start = id_pos[curr_non_copy_vec_id - 1] + 1;
    curr_gt_indexes.assign(decompress_gt_indexes_io.begin()+ index_start, decompress_gt_indexes_io.begin() + id_pos[curr_non_copy_vec_id]);
    sparse_matrix_cols = vint_code::DecodeArray(curr_gt_indexes);
    id_pos_size = sparse_matrix_cols.size();
    int prev = -1;
    while (whichByte_whereInRes[next_haplotype].first < last_byte && next_haplotype < whichByte_whereInRes.back().second)
    {
        if (cur_i < id_pos_size)
        {
            if (cur_i != prev_i)
            {
                cur_index = sparse_matrix_cols[cur_i]+prev;
                prev = cur_index;
                cur_pos = cur_index % 8;
                cur_id = cur_index / 8;
            }
            prev_i = cur_i;
            if (cur_id < whichByte_whereInRes[next_haplotype].first)
                cur_i += 1;
            else if (cur_id == whichByte_whereInRes[next_haplotype].first)
            {
                resAll[vec_start + whichByte_whereInRes[next_haplotype].second] += perm_lut8[cur_pos];
                cur_i += 1;
                count++;
            }
            else
            {
                next_haplotype++;
                if (count)
                {
                    if (next_haplotype < whichByte_whereInRes.back().second)
                        resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = whichByte_whereInRes[next_haplotype].first == whichByte_whereInRes[next_haplotype - 1].first ? resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] : 0;
                }
                else
                    resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] = 0;
            }
        }
        else
        {
            next_haplotype++;
            if (count)
            {
                if (next_haplotype < whichByte_whereInRes.back().second)
                    resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = whichByte_whereInRes[next_haplotype].first == whichByte_whereInRes[next_haplotype - 1].first ? resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] : 0;
            }
            else
                resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] = 0;
        }
    }
    // next_haplotype = 0;



    // if (parity)
    // {
    //     if (copy_bit_vector[parity][vector]) // If vector is a copy od other vector (certainly it is placed within the same block)
    //     {
    //         if (!full_decode)
    //         {
    //             is_uniqe_id = false;
    //             curr_copy++;
    //             return 0;
    //         }
    //         unsigned long long bit_pos = (curr_copy)*decompression_reader.used_bits_cp;
    //         decompression_reader.bm_comp_copy_orgl_id.SetPos(bit_pos >> 3);
    //         decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, bit_pos & 7); // %8
    //         decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, decompression_reader.used_bits_cp);

    //         // Here curr_non_copy_vec_id is a fake curr_non_copy_vec_id (it is ud of the next non_copy_vec_id)

    //         curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy - curr_non_copy_vec_id_offset;

    //         curr_non_copy_vec_id = curr_non_copy_vec_id - tmp - 1;
    //         // cout<<"unique_pos_first_in_block:"<<curr_non_copy_vec_id<<" "<<unique_pos_first_in_block<<endl;
    //         is_uniqe_id = false;
    //         curr_copy++;
    //         memcpy(resAll + vec_start, resUnique + (curr_non_copy_vec_id - unique_pos_first_in_block) * no_haplotypes, no_haplotypes);
            

    //         return 0;
    //     }
    //     else if (zeros_bit_vector[1][vector])
    //     {
    //         cout<<zeros_bit_vector[1][vector]<<":"<<zeros_bit_vector[0][vector]<<":"<<vector<<":"<<vec_id<<endl;
    //         is_uniqe_id = false;
    //         curr_zeros++;
    //         if (!full_decode)
    //             return 0;
    //         fill_n(resAll + vec_start, no_haplotypes, 0);
    //         return 0;
    //     }
    //     else
    //     {
    //         curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy - curr_non_copy_vec_id_offset;
    //         is_uniqe_id = true;
    //     }
    //     fill_n(resAll + no_haplotypes, no_haplotypes, 0);

    //     if (curr_non_copy_vec_id == 0)
    //         index_start = 0;
    //     else
    //         index_start = id_pos[curr_non_copy_vec_id - 1] + 1;
    //     curr_gt_indexes.assign(decompress_gt_indexes_io.begin()+ index_start, decompress_gt_indexes_io.begin() + id_pos[curr_non_copy_vec_id]);
    //     sparse_matrix_cols = vint_code::DecodeArray(curr_gt_indexes);
    //     id_pos_size = sparse_matrix_cols.size();
    //     int prev = -1;
    //     while (whichByte_whereInRes[next_haplotype].first < last_byte && next_haplotype < whichByte_whereInRes.back().second)
    //     {
    //         if (cur_i < id_pos_size)
    //         {
    //             if (cur_i != prev_i)
    //             {
    //                 cur_index = sparse_matrix_cols[cur_i]+prev;
    //                 prev = cur_index;
    //                 cur_pos = cur_index % 8;
    //                 cur_id = cur_index / 8;
    //             }
    //             prev_i = cur_i;
    //             if (cur_id < whichByte_whereInRes[next_haplotype].first)
    //                 cur_i += 1;
    //             else if (cur_id == whichByte_whereInRes[next_haplotype].first)
    //             {
    //                 resAll[vec_start + whichByte_whereInRes[next_haplotype].second] += perm_lut8[cur_pos];
    //                 cur_i += 1;
    //                 count++;
    //             }
    //             else
    //             {
    //                 next_haplotype++;
    //                 if (count)
    //                 {
    //                     if (next_haplotype < whichByte_whereInRes.back().second)
    //                         resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = whichByte_whereInRes[next_haplotype].first == whichByte_whereInRes[next_haplotype - 1].first ? resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] : 0;
    //                 }
    //                 else
    //                     resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] = 0;
    //             }
    //         }
    //         else
    //         {
    //             next_haplotype++;
    //             if (count)
    //             {
    //                 if (next_haplotype < whichByte_whereInRes.back().second)
    //                     resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = whichByte_whereInRes[next_haplotype].first == whichByte_whereInRes[next_haplotype - 1].first ? resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] : 0;
    //             }
    //             else
    //                 resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] = 0;
    //         }
    //     }
    //     next_haplotype = 0;
    // }
    // else if (!zeros_bit_vector[0][vector])
    // {
    //     is_uniqe_id = false;
    //     curr_zeros++;
    //     if (!full_decode)
    //         return 0;
        
    //     fill_n(resAll + vec_start, no_haplotypes, 0);
    //     return 0;
    // }
    // else
    // {
    //     if (copy_bit_vector[parity][vector])
    //     {
    //         if (!full_decode)
    //         {
    //             is_uniqe_id = false;
    //             curr_copy++;
    //             return 0;
    //         }
    //         unsigned long long bit_pos = (curr_copy)*decompression_reader.used_bits_cp;
    //         decompression_reader.bm_comp_copy_orgl_id.SetPos(bit_pos >> 3);
    //         decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, bit_pos & 7); // %8
    //         decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, decompression_reader.used_bits_cp);
    //         // Here curr_non_copy_vec_id is a fake curr_non_copy_vec_id (it is ud of the next non_copy_vec_id)
    //         curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy - curr_non_copy_vec_id_offset;
    //         curr_non_copy_vec_id = curr_non_copy_vec_id - tmp - 1;

    //         is_uniqe_id = false;
    //         curr_copy++;
    //         memcpy(resAll + vec_start, resUnique + (curr_non_copy_vec_id - unique_pos_first_in_block) * no_haplotypes, no_haplotypes);
            
    //         return 0;
    //     }
    //     else
    //     {
    //         curr_non_copy_vec_id = vec_id - curr_zeros - curr_copy - curr_non_copy_vec_id_offset;
    //         is_uniqe_id = true;
    //     }

    //     fill_n(resAll, no_haplotypes, 0);

    //     if (curr_non_copy_vec_id == 0)
    //         index_start = 0;
    //     else
    //         index_start = id_pos[curr_non_copy_vec_id - 1] + 1;
    //     curr_gt_indexes.assign(decompress_gt_indexes_io.begin()+ index_start, decompress_gt_indexes_io.begin() + id_pos[curr_non_copy_vec_id]);
    //     sparse_matrix_cols = vint_code::DecodeArray(curr_gt_indexes);
    //     int prev = -1;
    //     id_pos_size = sparse_matrix_cols.size();
    //     while (whichByte_whereInRes[next_haplotype].first < last_byte && next_haplotype < whichByte_whereInRes.back().second)
    //     {
    //         if (cur_i < id_pos_size)
    //         {
    //             if (cur_i != prev_i)
    //             {
    //                 cur_index = sparse_matrix_cols[cur_i]+prev;
    //                 prev = cur_index;
    //                 cur_pos = cur_index % 8;
    //                 cur_id = cur_index / 8;
    //             }
    //             prev_i = cur_i;
    //             if (cur_id < whichByte_whereInRes[next_haplotype].first)
    //                 cur_i += 1;
    //             else if (cur_id == whichByte_whereInRes[next_haplotype].first)
    //             {
    //                 resAll[vec_start + whichByte_whereInRes[next_haplotype].second] += perm_lut8[cur_pos];
    //                 cur_i += 1;
    //                 count++;
    //             }
    //             else
    //             {
    //                 next_haplotype++;
    //                 if (count)
    //                 {
    //                     if (next_haplotype < whichByte_whereInRes.back().second)
    //                         resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = whichByte_whereInRes[next_haplotype].first == whichByte_whereInRes[next_haplotype - 1].first ? resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] : 0;
    //                 }
    //                 else
    //                     resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] = 0;
    //             }
    //         }
    //         else
    //         {
    //             next_haplotype++;
    //             if (count)
    //             {
    //                 if (next_haplotype < whichByte_whereInRes.back().second)
    //                     resAll[vec_start + whichByte_whereInRes[next_haplotype].second] = whichByte_whereInRes[next_haplotype].first == whichByte_whereInRes[next_haplotype - 1].first ? resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] : 0;
    //             }
    //             else
    //                 resAll[vec_start + whichByte_whereInRes[next_haplotype - 1].second] = 0;
    //         }
    //     }  
    // }
    return 0;
}

// // *****************************************************************************************************************
void Decompressor::calculate_start_position(int &_start_block,int &_start_position)
{
    pos_to_block.resize(fixed_variants_chunk_io.size(), 0);

    for (size_t i = 0; i < fixed_variants_chunk_io.size(); i++)
    {

        pos_to_block[i] = fixed_variants_chunk_io[i].data_compress[0].pos;
        
    }
    variant_desc_t cur_desc;
    _start_block = lower_bound(pos_to_block.begin(), pos_to_block.end(), range_1) - pos_to_block.begin() - 1;
    
    cur_desc.pos = range_1;

    if (_start_block < 0)
    {

        _start_block = 0;
        _start_position = 0;
    }
    else
    {
        if (cur_desc.pos > fixed_variants_chunk_io[_start_block].data_compress[fixed_variants_chunk_io[_start_block].data_compress.size() - 1].pos)
        {
            _start_position = 0;
            _start_block++;
        }
        else
        {
            _start_position = upper_bound(fixed_variants_chunk_io[_start_block].data_compress.begin(), 
                fixed_variants_chunk_io[_start_block].data_compress.end(), cur_desc, comp) - fixed_variants_chunk_io[_start_block].data_compress.begin();
            _start_position = _start_position > 0 ? _start_position : 0;
        }
    }    
    
}
void Decompressor::calculate_end_position(int &_end_block,int &_end_position){


    pos_to_block.resize(fixed_variants_chunk_io.size(), 0);

    for (size_t i = 0; i < fixed_variants_chunk_io.size(); i++)
    {

        pos_to_block[i] = fixed_variants_chunk_io[i].data_compress[0].pos;
    }
    variant_desc_t cur_desc;
    _end_block = upper_bound(pos_to_block.begin(), pos_to_block.end(), range_2) - pos_to_block.begin() - 1;
    
    cur_desc.pos = range_2;
    _end_position = upper_bound(fixed_variants_chunk_io[_end_block].data_compress.begin(), 
        fixed_variants_chunk_io[_end_block].data_compress.end(), cur_desc, comp1) - fixed_variants_chunk_io[_end_block].data_compress.begin();
    _end_position = _end_position > 0 ? _end_position : 0;

}
bool Decompressor::splitFileWriting(int file_num){

    split_files.resize(file_num);
    for(int i = 0; i < file_num; i++){
        char write_mode[5] = "wb-";
        if (out_type == file_type::BCF_File)
        {
            write_mode[3] = compression_level;
            write_mode[4] = '\0';
        }
        if (out_file_name != "")
        {
            char *gz_fname = (char *)malloc(strlen((out_file_name + decompression_reader.d_where_chrom[i].first).c_str()) + 5);

            if (out_type == file_type::VCF_File)
            {
                snprintf(gz_fname, strlen((out_file_name + decompression_reader.d_where_chrom[i].first).c_str()) + 5, "%s.vcf", (out_file_name + decompression_reader.d_where_chrom[i].first).c_str());
                split_files[i] = hts_open(gz_fname, "w");
            }
            else
            {

                snprintf(gz_fname, strlen((out_file_name + decompression_reader.d_where_chrom[i].first).c_str()) + 5, "%s.bcf", (out_file_name + decompression_reader.d_where_chrom[i].first).c_str());
                split_files[i] = hts_open(gz_fname, write_mode);
            }

            free(gz_fname);
            if (!split_files[i]){
                std::cout << "could not open " << split_files[i] << " file" << std::endl;
                return false;
            }
            else
            {
                hts_set_opt(split_files[i], HTS_OPT_CACHE_SIZE, 32000000);
                rec = bcf_init();
            }
        }
    }
    return true;
}
// *****************************************************************************************************************
bool Decompressor::OpenForWriting()
{
    
    char write_mode[5] = "wb-";
    if (out_type == file_type::BCF_File)
    {
        write_mode[3] = compression_level;
        write_mode[4] = '\0';
    }


    if (out_file_name != "")
    {
        char *gz_fname = (char *)malloc(strlen(out_file_name.c_str()) + 5);

        if (out_type == file_type::VCF_File)
        {
            snprintf(gz_fname, strlen(out_file_name.c_str()) + 5, "%s.vcf", out_file_name.c_str());
            out_file = hts_open(gz_fname, "w");
        }
        else if(out_type == file_type::BCF_File)
        {

            snprintf(gz_fname, strlen(out_file_name.c_str()) + 5, "%s.bcf", out_file_name.c_str());
            out_file = hts_open(gz_fname, write_mode);
        }
        else{
            
            if(!out_fam.Open(out_file_name + ".fam" , "w"))
	        {
		        cerr << "Cannot open " << out_file_name << ".fam file\n";
		        return false;
	        }
            if(!out_bim.Open(out_file_name + ".bim" , "w"))
	        {
		        cerr << "Cannot open " << out_file_name << ".bim file\n";
		        return false;
	        }
            if(!out_bed.Open(out_file_name + ".bed" , "wb"))
	        {
		        cerr << "Cannot open " << out_file_name << ".bed file\n";
		        return false;
	        }
        }
        free(gz_fname);
    }
    else
    {
        if (out_type == file_type::VCF_File)
            out_file = hts_open("-", "w");

        else if(out_type == file_type::BCF_File)
            out_file = hts_open("-", write_mode);
        else{
            if(!out_fam.Open("-" , "w"))
	        {

		        return false;
	        }
            if(!out_bim.Open("-" , "w"))
	        {
		        
		        return false;
	        }
            if(!out_bed.Open("-" , "wb"))
	        {
		        
		        return false;
	        }
        } 

    }
    if(out_type != file_type::BED_File){
        if (!out_file){
            std::cout << "could not open " << out_file << " file" << std::endl;
            return false;
        }
        else
        {
            hts_set_opt(out_file, HTS_OPT_CACHE_SIZE, 32000000);
            rec = bcf_init();
        }
    }
    return true;
}
// *****************************************************************************************************************
int Decompressor::analyzeInputSamples(vector<string> &v_samples)
{
    
    if (smpl.loadSamples(v_samples))
        return 1;
    
    return 0;
}
int Decompressor::initOutSplitFile(){

    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t";
    string str = "";
    smpl.get_all_samples(str);
    out_hdr = bcf_hdr_init("r");
    bcf_hdr_parse(out_hdr, (char *)header.c_str());
    bcf_hdr_sync(out_hdr);

    char delim = '\t';
    std::stringstream ss(str);
    std::string item;
    while (getline(ss, item, delim))
    {
        bcf_hdr_add_sample(out_hdr, item.c_str());
    }
    bcf_hdr_sync(out_hdr);
    
    for(size_t i = 0; i < split_files.size(); i++){
        if (bcf_hdr_write(split_files[i], out_hdr) < 0)
            return 1;
    }
    return 0;
}
// *****************************************************************************************************************
void Decompressor::WriteBEDMagicNumbers() {
    // 写入前两个固定的魔术数字字节
    out_bed.PutByte(0x6C);  // 01101100
    out_bed.PutByte(0x1B);  // 00011011
    
    // 写入第三个字节，这里我们假设以样本为主的存储格式
    out_bed.PutByte(0x01);  // 00000001

    // 如果你需要写入以SNP为主的存储格式，那么你会使用
    // PutByte(0x00);  // 00000000
}
int Decompressor::initOut()
{
    if(out_type == file_type::BED_File){
        string str = "";
        smpl.get_all_samples(str);
        char delim = '\t';
        std::stringstream ss(str);
        std::string item;
        string fam_line;
        while (getline(ss, item, delim))
        {
            fam_line = "0\t"+ item +"\t0\t0\t0\t-9\n";
            // cout<< fam_line << endl;
            out_fam.Write(fam_line.c_str(), fam_line.size());
        }
        WriteBEDMagicNumbers();
    }
    else{
        header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"; 
        string str = "";
        if (params.samples == "")
        {
            smpl.get_all_samples(str);
        }
        else
        {
            
            sampleIDs = smpl.setSamples(params.samples.c_str(), str);
        }
           
        if (!out_file)
        {
            return 1;
        }
        out_hdr = bcf_hdr_init("r");
        if (out_genotypes)
            header += "FORMAT";
        bcf_hdr_parse(out_hdr, (char *)header.c_str());
        if(params.compress_mode == compress_mode_t::lossly_mode){
            bcf_hdr_remove(out_hdr,BCF_HL_INFO, NULL);
            bcf_hdr_remove(out_hdr,BCF_HL_FMT, NULL);
            bcf_hdr_append(out_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        }
        bcf_hdr_sync(out_hdr);
        if(params.out_AC_AN){
            bcf_hdr_append(out_hdr, "##INFO=<ID=AC,Number=A,Type=String,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">");
            bcf_hdr_append(out_hdr, "##INFO=<ID=AN,Number=A,Type=String,Description=\"Total number of alleles in called genotypes\">");
            bcf_hdr_sync(out_hdr);
        }

        if (out_genotypes)
        {
            char delim = '\t';
            std::stringstream ss(str);
            std::string item;
            while (getline(ss, item, delim))
            {
                if (out_genotypes)
                    bcf_hdr_add_sample(out_hdr, item.c_str());
            }
            bcf_hdr_sync(out_hdr);
        }

        if (bcf_hdr_write(out_file, out_hdr) < 0)
            return 1;
    // }
    }
    return 0;
}
// *****************************************************************************************************************
// void Decompressor::setMemoryUsage()
// {
//     uint32_t no_haplotypes = smpl.no_samples * decompression_reader.ploidy;
//     // Set memory usage
//     if (params.MB_memory)
//     {
//         if (params.max_MB_memory)
//             max_stored_unique = (params.max_MB_memory * 1000000) / decompression_reader.vec_len;
//         else
//             max_stored_unique = decompression_reader.no_vec;
//         if (BLOCK_SIZE < max_stored_unique)
//             max_stored_unique = no_haplotypes*2;
//     }
//     else
//         max_stored_unique = 0;

// }

// // *****************************************************************************************************************
void Decompressor::decoded_vector_row(uint64_t vec_id, uint64_t offset, uint64_t _curr_non_copy_vec_id_offset, uint64_t length, int pos, uint8_t *decomp_data)
{

    uint64_t id, curr_non_copy_vec_id, toDelete;
    uint32_t tmp = 0;
    size_t id_pos_size;
    uint64_t index_start;
    int cur_id;
    int cur_pos;
    uint32_t cur_index;

    id = vec_id >> 1; // vec_id/2
    uint8_t parity = vec_id & 1;
    if (decompression_reader.rrr_zeros_bit_vector[parity][id])
    {

        memcpy(decomp_data + pos, zeros_only_vector, decompression_reader.vec_len);

        // pos += decompression_reader.vec_len;
        
        return;
    }
    else if (decompression_reader.rrr_copy_bit_vector[parity][id])
    {

        unsigned long long bit_pos = (rrr_rank_copy_bit_vector[0](id + ((parity))) + rrr_rank_copy_bit_vector[1](id)) * decompression_reader.used_bits_cp;
        decompression_reader.bm_comp_copy_orgl_id.SetPos(bit_pos >> 3);      // /8
        decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, bit_pos & 7); // %8
        decompression_reader.bm_comp_copy_orgl_id.GetBits(tmp, decompression_reader.used_bits_cp);

        // Here curr_non_copy_vec_id is a fake curr_non_copy_vec_id (it is ud of the next non_copy_vec_id)
        // prev_chr_zeros_copy:Count the total number of all 1s and copied rows before the start position
        curr_non_copy_vec_id = vec_id - rrr_rank_zeros_bit_vector[0](id + ((parity))) - rrr_rank_zeros_bit_vector[1](id) -
                                rrr_rank_copy_bit_vector[0](id + ((parity))) - rrr_rank_copy_bit_vector[1](id) -_curr_non_copy_vec_id_offset;
        curr_non_copy_vec_id = curr_non_copy_vec_id - tmp - 1;

        got_it = done_unique.find(curr_non_copy_vec_id);
        
        if (got_it != done_unique.end())
        {
            memcpy(decomp_data + pos, got_it->second, decompression_reader.vec_len);
            // pos += decompression_reader.vec_len;
            return;
        }
    }
    else // unique and not a copy - no need to check if it was previously decompressed (got_it)
    {

        curr_non_copy_vec_id = vec_id - rrr_rank_zeros_bit_vector[0](id + ((parity))) - rrr_rank_zeros_bit_vector[1](id) -
                                rrr_rank_copy_bit_vector[0](id + ((parity))) - rrr_rank_copy_bit_vector[1](id) -_curr_non_copy_vec_id_offset;

        
    }
    uint8_t *cur_decomp_data = new uint8_t[decompression_reader.vec_len];
    fill_n(cur_decomp_data, decompression_reader.vec_len, 0);
  
    if (curr_non_copy_vec_id == 0)
        index_start = 0;
    else
        index_start = id_pos[curr_non_copy_vec_id - 1] + 1;

    // id_pos_size = id_pos[curr_non_copy_vec_id] - index_start;
    vector<uint8_t> curr_gt_indexes;
    curr_gt_indexes.assign(decompress_gt_indexes_io.begin()+ index_start, decompress_gt_indexes_io.begin() + id_pos[curr_non_copy_vec_id]);
    sparse_matrix_cols = vint_code::DecodeArray(curr_gt_indexes);
    id_pos_size = sparse_matrix_cols.size();
    int prev = -1;
    for(size_t i = 0; i < id_pos_size; ++i)
    {
        cur_index = sparse_matrix_cols[i]+prev;
        prev = cur_index;
        cur_pos = cur_index % 8;
        cur_id = cur_index / 8;
        cur_decomp_data[cur_id] += perm_lut8[cur_pos];
    }

    memcpy(decomp_data + pos, cur_decomp_data, decompression_reader.vec_len);
    // pos += decompression_reader.vec_len;

    if (max_stored_unique)
    {

        uint8_t *vector;
        if (done_unique.size() > max_stored_unique)
        {
            toDelete = stored_unique.back();
            stored_unique.pop_back();
            vector = done_unique[toDelete];
            done_unique.erase(toDelete);
        }
        else
        {
            vector = new uint8_t[decompression_reader.vec_len];
        }
        memcpy(vector, decomp_data + pos ,decompression_reader.vec_len);

        done_unique[curr_non_copy_vec_id] = vector;
        stored_unique.push_front(curr_non_copy_vec_id);

    }
    delete[] cur_decomp_data;
}
//获取原始每行数据的字节
// *****************************************************************************************************************
void Decompressor::decode_perm_rev(int vec2_start, const vector<uint32_t> &rev_perm, uint8_t *decomp_data_perm, uint8_t *decomp_data) 
{
    
    for(int i = 0;i<(int)decompression_reader.vec_len*2;i++)
        decomp_data_perm[i] = map_t256[decomp_data_perm[i]];
    size_t x;
    for (x = 0; x + 8 < standard_block_size;)
    {
        int x8 = x / 8;
        int d_x8 = decomp_data_perm[x8];
        int d2_x8 = decomp_data_perm[vec2_start + x8];
        if (!d_x8 && !d2_x8)
        {
            x += 8;
            continue;
        }

        int x_p = 1 << 7;   
        for (int i = 0; i < 8; ++i, ++x, x_p >>= 1)
        {
            auto j = rev_perm[x];
            uint8_t j_p = perm_lut8[j % 8];
            int j8 = j / 8;
            if (d_x8 & x_p)
                decomp_data[j8] += j_p;
            if (d2_x8 & x_p)
                decomp_data[vec2_start + j8] += j_p;
        }
    }
    int x8 = x / 8;
    uint8_t d_x8 = decomp_data_perm[x8];
    uint8_t d2_x8 = decomp_data_perm[vec2_start + x8];
    
    for (; x < standard_block_size; ++x)
    {
        auto j = rev_perm[x];
        uint8_t x_p = perm_lut8[x % 8];
        uint8_t j_p = perm_lut8[j % 8];
        int j8 = j / 8;
        if (d_x8 & x_p)
            decomp_data[j8] += j_p;
        if (d2_x8 & x_p)
            decomp_data[vec2_start + j8] += j_p;
    }
    
}
//获取原始的perm顺序
// *****************************************************************************************************************
void inline Decompressor::reverse_perm(const vector<uint32_t> &perm, vector<uint32_t> &rev_perm, int no_haplotypes)
{
    for (int i = 0; i < no_haplotypes; ++i)
        rev_perm[perm[i]] = i;
}
