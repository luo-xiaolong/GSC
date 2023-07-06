#include "compression_reader.h"

// ***************************************************************************************************************************************
// ***************************************************************************************************************************************
bool CompressionReader::OpenForReading(string & file_name)
{

    if(in_file)
        hts_close(in_file);
    
    if (in_type == file_type::VCF_File)
    {
        in_file = hts_open(file_name.c_str(), "r");
    }
    else
    {
        in_file = hts_open(file_name.c_str(), "rb");
    }
    if(!in_file){
        std::cout << "could not open " << in_file_name << " file" << std::endl;
        return false;
    }
    hts_set_opt(in_file, HTS_OPT_CACHE_SIZE, 32 << 20);
    hts_set_opt(in_file, HTS_OPT_BLOCK_SIZE, 32 << 20);
    if(vcf_hdr)
        bcf_hdr_destroy(vcf_hdr);
    vcf_hdr = bcf_hdr_read(in_file);
    vcf_record = bcf_init();
    return true;
}

// ***************************************************************************************************************************************
// Get number of samples in VCF
bool CompressionReader::ReadFile()
{
    if(!in_file || !vcf_hdr)
        return -1;
        
    no_samples = bcf_hdr_nsamples(vcf_hdr);
    
    if (!vcf_hdr_read)
    {
        
        for (size_t i = 0; i < no_samples; i++)
            samples_list.emplace_back(vcf_hdr->samples[i]);

        
        vcf_hdr_read = true;
    }
    return true;
}
// Get number of samples in VCF
// ***************************************************************************************************************************************
uint32_t CompressionReader::GetSamples(vector<string> &s_list)
{

    if (!vcf_hdr_read)
    {
        ReadFile();
    }
    s_list = samples_list;

    return no_samples;
}

// ***************************************************************************************************************************************
bool CompressionReader::GetHeader(string &v_header)
{

    if (!vcf_hdr)
        return false;
    kstring_t str = {0, 0, nullptr};
    // bcf_hdr_remove(vcf_hdr,BCF_HL_INFO, NULL);
    bcf_hdr_format(vcf_hdr, 0, &str);
    char *ptr = strstr(str.s, "#CHROM");
    v_header.assign(str.s, ptr - str.s);
    if(str.m)
        free(str.s);
    return true;
}


// ***************************************************************************************************************************************
void CompressionReader::GetWhereChrom(vector<pair<std::string, uint32_t>> &_where_chrom,vector<int64_t> &_chunks_min_pos)
{
    _where_chrom = where_chrom;
    _chunks_min_pos = chunks_min_pos;
}
// ***************************************************************************************************************************************
vector<uint32_t> CompressionReader::GetActualVariants()
{
    return actual_variants;
}
// ***************************************************************************************************************************************
bool CompressionReader::setBitVector()
{
    if (no_samples == 0)
        return false;
    vec_len = (no_samples * ploidy) / 8 + (((no_samples * ploidy) % 8) ? 1 : 0);
    block_max_size = vec_len * no_vec_in_block + 1;
    // phased_block_max_size = ((no_samples *  (ploidy-1)) / 8 + (((no_samples *  (ploidy-1)) % 8)?1:0))*no_vec_in_block/2;
    bv.Create(block_max_size);
    
    vec_read_fixed_fields = 0;
    fixed_fields_id = 0;
    block_id = 0;
    vec_read_in_block = 0;
    return true;
}
// ***************************************************************************************************************************************
void CompressionReader::GetOtherField(vector<key_desc> &_keys, uint32_t &_no_keys,int &_key_gt_id)
{
    _keys = keys;
    _no_keys = no_keys;
    _key_gt_id = key_gt_id;

}
void CompressionReader::UpdateKeys(vector<key_desc> &_keys)
{
    // vector<key_desc>  tmp_keys = _keys;
    size_t k = 0;
    if(order.size() < no_keys - no_flt_keys){
        for(size_t i = 0; i < _keys.size();i++)
        {
            if(_keys[i].keys_type == key_type_t::info)
            {

                auto it = std::find(order.begin(), order.end(), _keys[i].actual_field_id);
                if (it != order.end()) {
                    _keys[i].actual_field_id = order[k++];
                }
            }
            else if(_keys[i].keys_type == key_type_t::fmt)
            {
                auto it = std::find(order.begin(), order.end(), _keys[i].actual_field_id);
                if (it != order.end()) {
                    _keys[i].actual_field_id = order[k++];
                }
            }
        }
    }else{
        for(size_t i = 0; i < _keys.size();i++)
        {
            if(_keys[i].keys_type == key_type_t::info)
            {

                _keys[i].actual_field_id = order[k++];
            }
            else if(_keys[i].keys_type == key_type_t::fmt)
            {
                _keys[i].actual_field_id = order[k++];

            }
        }        
    }

}
// ***************************************************************************************************************************************
void CompressionReader::InitVarinats(File_Handle_2 *_file_handle2)
{
    file_handle2 = _file_handle2;
    GetFilterInfoFormatKeys(no_flt_keys, no_info_keys, no_fmt_keys, keys);
    no_keys = no_flt_keys + no_info_keys + no_fmt_keys;
    v_o_buf.resize(no_keys);
    for (uint32_t i = 0; i < no_keys; i++)
        v_o_buf[i].SetMaxSize(max_buffer_size, 0);
    v_buf_ids_size.resize(no_keys, -1);
	v_buf_ids_data.resize(no_keys, -1);
    // cout << "no_keys:" << no_keys << endl;
    for (uint32_t i = 0; i < no_keys; i++)
    {
        switch (keys[i].keys_type)
        {
        case key_type_t::flt:
            if (keys[i].key_id >= FilterIdToFieldId.size())
                FilterIdToFieldId.resize(keys[i].key_id + 1, -1);
            FilterIdToFieldId[keys[i].key_id] = i;
            break;
        case key_type_t::info:
            // cout << "no_samples:" <<keys[i].key_id << endl;
            if (keys[i].key_id >= InfoIdToFieldId.size())
                InfoIdToFieldId.resize(keys[i].key_id + 1, -1);
            InfoIdToFieldId[keys[i].key_id] = i;
            break;
        case key_type_t::fmt:
            if (keys[i].key_id >= FormatIdToFieldId.size())
                FormatIdToFieldId.resize(keys[i].key_id + 1, -1);
            FormatIdToFieldId[keys[i].key_id] = i;
            // cout<<"InfoIdToFieldId:"<<FormatIdToFieldId[keys[i].key_id]<<endl;
            break;
        }
    }
    
    // file_handle2->RegisterStream("start");
    for (uint32_t i = 0; i < no_keys; i++)
    {
        v_buf_ids_size[i] = file_handle2->RegisterStream("key_" + to_string(i) + "_size");
        // cout<<v_buf_ids_size[i]<<endl;
    }
    for (uint32_t i = 0; i < no_keys; i++)
    {
        v_buf_ids_data[i] = file_handle2->RegisterStream("key_" + to_string(i) + "_data");
        // cout<<v_buf_ids_data[i]<<endl;

    }    
    

}
// ***************************************************************************************************************************************
bool CompressionReader::GetFilterInfoFormatKeys(int &no_flt_keys, int &no_info_keys, int &no_fmt_keys, vector<key_desc> &keys)
{
    if (!vcf_hdr)
        return false;

    no_flt_keys = 0;
    no_info_keys = 0;
    no_fmt_keys = 0;

    key_desc new_key;
    // cout<<vcf_hdr->nhrec<<endl;
    for (int i = 0; i < vcf_hdr->nhrec; i++)
    {

        if (vcf_hdr->hrec[i]->type == BCF_HL_FLT || vcf_hdr->hrec[i]->type == BCF_HL_INFO || vcf_hdr->hrec[i]->type == BCF_HL_FMT)
        {
            // checking if it is a duplicate; if so, curr_dict_id different (and not increased)
            int id = bcf_hdr_id2int(vcf_hdr, BCF_DT_ID, vcf_hdr->hrec[i]->vals[0]);

            new_key.key_id = id;
            if (vcf_hdr->hrec[i]->type == BCF_HL_FLT)
            {
                no_flt_keys++;
                new_key.keys_type = key_type_t::flt;
                new_key.type = BCF_HT_FLAG;
            }
            else if (vcf_hdr->hrec[i]->type == BCF_HL_INFO || vcf_hdr->hrec[i]->type == BCF_HL_FMT)
            {
                if (vcf_hdr->hrec[i]->type == BCF_HL_FMT) // FORMAT
                {
                    no_fmt_keys++;
                    new_key.keys_type = key_type_t::fmt;
                    new_key.type = bcf_hdr_id2type(vcf_hdr, BCF_HL_FMT, id);
                    if (strcmp(vcf_hdr->id[BCF_DT_ID][id].key, "GT") == 0)
                    {
                        key_gt_id = (int) keys.size();
                    }
                   
                }
                else // INFO
                {
                    no_info_keys++;
                    new_key.keys_type = key_type_t::info;
                    new_key.type = bcf_hdr_id2type(vcf_hdr, BCF_HL_INFO, id);
                }
            }
            new_key.actual_field_id = (uint32_t) keys.size();
            // cout<<"new_key.actual_field_id:"<<new_key.key_id<<":"<<new_key.actual_field_id<<endl;
            keys.emplace_back(new_key);
        }
    }
    keys.shrink_to_fit();
    return true;
}
// ************************************************************************************
bool CompressionReader::GetVariantFromRec(bcf1_t *rec, vector<field_desc> &fields)
{
    vector<int> field_order;
    field_order.reserve(no_keys);
    if (rec->d.n_flt)
    {
        
        for (int i = 0; i < rec->d.n_flt; ++i)
        {
            fields[FilterIdToFieldId[rec->d.flt[i]]].present = true; // DEV, TMP
        }
        
    }

    int curr_size = 0;
                 
    // INFO
    if (rec->n_info)
    {
        for (int i = 0; i < rec->n_info; ++i)
        {

            bcf_info_t *z = &rec->d.info[i];
            // if(order.size() < no_keys - no_flt_keys || field_order_flag)
                field_order.emplace_back(InfoIdToFieldId[z->key]);
            auto &cur_field = fields[InfoIdToFieldId[z->key]];
            cur_field.present = true;
            if (!z->vptr)
                continue;
            if (z->key >= vcf_hdr->n[BCF_DT_ID])
            {
                hts_log_error("Invalid BCF, the INFO index is too large");
                return -1;
            }
            if (z->len <= 0)
                continue;

            int type = bcf_hdr_id2type(vcf_hdr, BCF_HL_INFO, z->key);
            int curr_size = 0;

            switch (type)
            {
            case BCF_HT_INT:

                curr_size = bcf_get_info_values(vcf_hdr, rec, vcf_hdr->id[BCF_DT_ID][z->key].key, &dst_int, &ndst_int, type);
                break;
            case BCF_HT_REAL:

                curr_size = bcf_get_info_values(vcf_hdr, rec, vcf_hdr->id[BCF_DT_ID][z->key].key, &dst_real, &ndst_real, type);
                break;
            case BCF_HT_STR:

                curr_size = bcf_get_info_values(vcf_hdr, rec, vcf_hdr->id[BCF_DT_ID][z->key].key, &dst_str, &ndst_str, type);
                break;
            case BCF_HT_FLAG:
                curr_size = bcf_get_info_values(vcf_hdr, rec, vcf_hdr->id[BCF_DT_ID][z->key].key, &dst_flag, &ndst_flag, type);
                break;
            }

            if (curr_size)
            {
                switch (type)
                {
                case BCF_HT_INT:
                    cur_field.data_size = (uint32_t)curr_size * 4;
                    cur_field.data = new char[cur_field.data_size];
                    memcpy(cur_field.data, (char *)dst_int, cur_field.data_size);
                    break;
                case BCF_HT_REAL:
                    cur_field.data_size = (uint32_t)curr_size * 4;
                    cur_field.data = new char[cur_field.data_size];
                    memcpy(cur_field.data, (char *)dst_real, cur_field.data_size);
                    break;
                case BCF_HT_STR:
                    cur_field.data_size = curr_size;
                    cur_field.data = new char[cur_field.data_size];
                    memcpy(cur_field.data, (char *)dst_str, cur_field.data_size);
                    break;
                case BCF_HT_FLAG:
                    cur_field.data = nullptr;
                    cur_field.data_size = 0;
                    break;
                }
            }
            else
            {
                cout << "Error getting variant" << endl;
                return false;
            }
        }
    }
    
    // FORMAT and individual information
    if (rec->n_sample)
    {
        
        int i; // , j;
        if (rec->n_fmt)
        {
            
            bcf_fmt_t *fmt = rec->d.fmt;
            for (i = 0; i < (int)rec->n_fmt; ++i)
            {
                
                if (!fmt[i].p)
                    continue;
                if (fmt[i].id < 0) //! bcf_hdr_idinfo_exists(h,BCF_HL_FMT,fmt[i].id) )
                {
                    hts_log_error("Invalid BCF, the FORMAT tag id=%d not present in the header", fmt[i].id);
                    abort();
                }

                // if(order.size() < no_keys - no_flt_keys || field_order_flag)
                
                    field_order.emplace_back(FormatIdToFieldId[fmt[i].id]);
                auto &cur_field = fields[FormatIdToFieldId[fmt[i].id]];
                cur_field.present = true;
                
                int bcf_ht_type = bcf_hdr_id2type(vcf_hdr, BCF_HL_FMT, fmt[i].id); // BCF_HT_INT or BCF_HT_REAL or BCF_HT_FLAG or BCF_HT_STR
                auto vcf_hdr_key = vcf_hdr->id[BCF_DT_ID][fmt[i].id].key;
                // cout<<vcf_hdr_key<<endl;
                if (strcmp(vcf_hdr_key, "GT") == 0)
                {
                    int *gt_arr = NULL, ngt_arr = 0;
                    curr_size = bcf_get_genotypes(vcf_hdr, rec, &gt_arr, &ngt_arr);
                    cur_field.data_size = curr_size - rec->n_sample;
                    cur_field.data = new char[cur_field.data_size];
                    int cur_phased_pos = 0;
                    for(uint32_t j = 0; j < rec->n_sample; ++j)
                    {
                        for(uint32_t k = 1; k < ploidy; ++k)
                        {   
                            if(bcf_gt_is_phased(gt_arr[j*ploidy+k]))
                                if(gt_arr[j*ploidy+k] == GT_NOT_CALL)
                                    cur_field.data[cur_phased_pos++] = '.';
                                else
                                    cur_field.data[cur_phased_pos++] = '|';
                            else
                                cur_field.data[cur_phased_pos++] = '/';
                        }

                    }
                    free(gt_arr);
                    continue;

                }
                else
                {
                    switch (bcf_ht_type)
                    {
                    case BCF_HT_INT:
                        curr_size = bcf_get_format_values(vcf_hdr, rec, vcf_hdr_key, &dst_int, &ndst_int, bcf_ht_type);
                        break;
                    case BCF_HT_REAL:
                        curr_size = bcf_get_format_values(vcf_hdr, rec, vcf_hdr_key, &dst_real, &ndst_real, bcf_ht_type);
                        break;
                    case BCF_HT_STR:
                        curr_size = bcf_get_format_values(vcf_hdr, rec, vcf_hdr_key, &dst_str, &ndst_str, bcf_ht_type);
                        break;
                    case BCF_HT_FLAG:
                        curr_size = bcf_get_format_values(vcf_hdr, rec, vcf_hdr_key, &dst_flag, &ndst_flag, bcf_ht_type);
                        break;
                    }
                }
                
                if (curr_size)
                {
                    // cur_field.data_size = (uint32_t)curr_size;

                    if (bcf_ht_type == BCF_HT_INT || bcf_ht_type == BCF_HT_REAL)
                    {
                        cur_field.data_size = (uint32_t)curr_size * 4;
                        cur_field.data = new char[cur_field.data_size];

                        if (bcf_ht_type == BCF_HT_INT)
                        {
                            memcpy(cur_field.data, dst_int, cur_field.data_size);
                        }
                        else if (bcf_ht_type == BCF_HT_REAL)
                            memcpy(cur_field.data, dst_real, cur_field.data_size);
                        // else if (bcf_ht_type == BCF_HT_STR)
                        //     memcpy(cur_field.data, dst_int, curr_size * 4); // GTs are ints!
                    }
                    else if (bcf_ht_type == BCF_HT_STR)
                    {
                        cur_field.data_size = (uint32_t)curr_size;
                        cur_field.data = new char[cur_field.data_size];
                        memcpy(cur_field.data, dst_str, cur_field.data_size);
                    }
                    else
                        assert(0);
                }
            }
        }
    }
    // if(order.size() < no_keys - no_flt_keys || field_order_flag){
    //     field_order_flag =  false;
        for (size_t i = 0; i < field_order.size() - 1; ++i) {
            // cout<<field_order[i]<<" ";
            if (field_order_graph[field_order[i]].find(field_order[i+1]) == field_order_graph[field_order[i]].end()) {
                field_order_graph[field_order[i]].insert(field_order[i+1]);
                inDegree[field_order[i+1]]++;
            }
            if (inDegree.find(field_order[i]) == inDegree.end()) {
               
                inDegree[field_order[i]] = 0;
            }    
        }
        // cout<<endl;
        // order = topo_sort(field_order_graph,inDegree);

    // }

    return true;
}

// *******************************************************************************************************************************
bool CompressionReader::SetVariantOtherFields(vector<field_desc> &fields)
{
    for (uint32_t i = 0; i < no_keys; i++)
    {

        switch (keys[i].type)
        {
        case BCF_HT_INT:

            v_o_buf[i].WriteInt(fields[i].data, fields[i].present ? fields[i].data_size : 0);

#ifdef LOG_INFO
            {
                uint32_t *pp = (uint32_t *)fields[i].data;
                for (int j = 0; j < fields[i].data_size; ++j)
                    distinct_values[i].insert(pp[j]);
            }
#endif

            break;
        case BCF_HT_REAL:
            v_o_buf[i].WriteReal(fields[i].data, fields[i].present ? fields[i].data_size : 0);

#ifdef LOG_INFO
            {
                uint32_t *pp = (uint32_t *)fields[i].data;

                for (int j = 0; j < fields[i].data_size; ++j)
                    distinct_values[i].insert(pp[j]);
            }
#endif

            break;
        case BCF_HT_STR:
            v_o_buf[i].WriteText(fields[i].data, fields[i].present ? fields[i].data_size : 0);
            break;
        case BCF_HT_FLAG:

            v_o_buf[i].WriteFlag(fields[i].present);
            break;
        }
        
        if (v_o_buf[i].IsFull())
		{
			// cout<<"v_buf_ids_size[i]:"<<v_buf_ids_size[i]<<endl;
			// cout<<"v_buf_ids_size[i]:"<<v_buf_ids_size[i]<<endl;
			auto part_id = file_handle2->AddPartPrepare(v_buf_ids_size[i]);
			file_handle2->AddPartPrepare(v_buf_ids_data[i]);
            // cout<<part_id<<endl;
			vector<uint32_t> v_size;
			vector<uint8_t> v_data;
			
			v_o_buf[i].GetBuffer(v_size, v_data);
            
			
			SPackage pck(i, v_buf_ids_size[i], v_buf_ids_data[i], part_id, v_size, v_data);

			part_queue->Push(pck);
		}
    }
    
    return true;
}
// ***************************************************************************************************************************************
// Splits multiple alleles sites, reads genotypes, creates blocks of bytes to process, fills out [archive_name].bcf file
bool CompressionReader::ProcessInVCF()
{
    
    if (!vcf_hdr_read)
        ReadFile();
    setBitVector();
    no_vec = 0;
    tmpi = 0;
    temp = 0;
    cur_pos = 0;
    gt_data = new int[no_samples*ploidy];
    no_actual_variants = 0;
    while (bcf_read1(in_file, vcf_hdr, vcf_record) >= 0)
    {
        variant_desc_t desc;
        if (vcf_record->errcode)
        {
            std::cout << "Repair VCF file\n";
            exit(9);
        }
        bcf_unpack(vcf_record, BCF_UN_ALL);
        if (vcf_record->d.fmt->n != (int)ploidy)
        {
            std::cout << "Wrong ploidy (not equal to " << ploidy << ") for record at position " << vcf_record->pos << ".\n";
            std::cout << "Repair VCF file OR set correct ploidy using -p option\n";
            exit(9);
        }
        if (tmpi % 100000 == 0){
            cout << tmpi << "\r";
            fflush(stdout);
        }
        if(compress_all){

            std::vector<field_desc> curr_field(keys.size());

            GetVariantFromRec(vcf_record, curr_field);

            SetVariantOtherFields(curr_field);
        

            for(size_t j = 0; j < keys.size(); ++j)
            {
                if(curr_field[j].data_size > 0)
                {
				    if(curr_field[j].data)
					    delete[] curr_field[j].data;
                        curr_field[j].data = nullptr;
                        curr_field[j].data_size = 0;
                }
            } 
            curr_field.clear();
        }
            
        ++no_actual_variants;
        ProcessFixedVariants(vcf_record, desc);

        tmpi++;
    }
    if(compress_all)
        order = topo_sort(field_order_graph,inDegree);
    if (cur_g_data)
    {
        delete[] cur_g_data;
    }
    if (gt_data)
    {
        delete[] gt_data;
    }

    std::cout << "Read all the variants and gentotypes" << endl;
    // Last pack (may be smaller than block size）
    
        CloseFiles();
    return true;
}
// void CompressionReader::addPhased(int * gt_data, int ngt_data){

//     for(int i = 0;i<ngt_data;i++)
//         if(bcf_gt_is_phased(gt_data[i&1]) )
//         {

//         }
//         else {
//             cout<<cur_g_data[i]<<" ";
//         }

// }
// ***************************************************************************************************************************************
void CompressionReader::ProcessFixedVariants(bcf1_t *vcf_record, variant_desc_t &desc){
    
    if (vcf_record->n_allele <= 2)
    {
        desc.ref.erase(); // REF
        if (vcf_record->n_allele > 0)
        {
            if (vcf_record->pos + 1 == cur_pos)
            {

                desc.ref = to_string(++temp) + vcf_record->d.allele[0];

            }
            else
            {
                temp = 0;
                desc.ref = vcf_record->d.allele[0];
            }
        }
        else
            desc.ref = '.';
        if (vcf_record->n_allele > 1)
        {
            desc.alt.erase(); // ALT
            desc.alt += vcf_record->d.allele[1];
        }

        bcf_get_genotypes(vcf_hdr, vcf_record, &cur_g_data, &ncur_g_data);

        for (int i = 0; i < ncur_g_data; i++)
        {
            gt_data[i] = bcf_gt_allele(cur_g_data[i]);
            
        }
        
        addVariant(gt_data, ncur_g_data, desc);
        bcf_clear(vcf_record);
    }
    else
    {
        if (vcf_record->n_allele == 3 && (strcmp(vcf_record->d.allele[2], "<M>") == 0 || strcmp(vcf_record->d.allele[2], "<N>") == 0))
        {

            desc.ref.erase(); // REF
            if (vcf_record->n_allele > 0)
            {
                if (vcf_record->pos + 1 == cur_pos)
                {
                    ;
                }
                else
                {
                    temp = 0;
                }
                desc.ref = to_string(++temp) + vcf_record->d.allele[0];
                // snprintf(temp_str, sizeof(temp_str), "%d", ++temp);
                // desc.ref = temp_str;
                // desc.ref += vcf_record->d.allele[0];
            }
            else
                desc.ref = '.';
            desc.alt.erase(); // ALT
            desc.alt += vcf_record->d.allele[1];
            desc.alt += ',';
            desc.alt += vcf_record->d.allele[2];
            bcf_get_genotypes(vcf_hdr, vcf_record, &cur_g_data, &ncur_g_data);
            // if (gt_data)
            //     delete[] gt_data;
            // gt_data = new int[ncur_g_data];
            for (int i = 0; i < ncur_g_data; i++)
            {
                gt_data[i] = bcf_gt_allele(cur_g_data[i]);
            }

            addVariant(gt_data, ncur_g_data, desc);
            bcf_clear(vcf_record);
        }
        else
        {
            
            if (vcf_record->pos + 1 != cur_pos)
                temp = 0;
            for (int n = 1; n < vcf_record->n_allele; n++) // create one line for each single allele
            {


                // size_t len_ref = strlen(vcf_record->d.allele[0]), len_alt = strlen(vcf_record->d.allele[n]);
                // if (len_ref > 1 && len_alt > 1)
                // {
                //     while ((vcf_record->d.allele[0][len_ref - 1] == vcf_record->d.allele[n][len_alt - 1]) && len_ref > 1 && len_alt > 1)
                //     {
                //         len_ref--;
                //         len_alt--;
                //     }
                // }

                desc.ref.erase(); // REF

                if (vcf_record->n_allele > 0)
                {
                    desc.ref += vcf_record->d.allele[0];
                    // desc.ref = to_string(++temp) + desc.ref.substr(0, len_ref);
                    desc.ref = to_string(++temp) + desc.ref;

                }
                else
                    desc.ref = '.';
                desc.alt.erase(); // ALT
                desc.alt += vcf_record->d.allele[n];
                // desc.alt = desc.alt.substr(0, len_alt);
                if(n == 1)
                    desc.alt += ",<N>";
                else
                    desc.alt += ",<M>";
                
                bcf_get_genotypes(vcf_hdr, vcf_record, &cur_g_data, &ncur_g_data);
                // if (gt_data)
                //     delete[] gt_data;
                // gt_data = new int[ncur_g_data];
                // fill_n(gt_data, ncur_g_data, 0);
                for (int i = 0; i < ncur_g_data; i++)
                { // gt_arr needed to create bit vectors
                    // cout<<bcf_gt_allele(cur_g_data[i])<<" ";
                    if (bcf_gt_allele(cur_g_data[i]) != 0)
                    {
                        if (bcf_gt_allele(cur_g_data[i]) == n)
                            cur_g_data[i] = bcf_gt_phased(1);
                        else if (bcf_gt_is_missing(cur_g_data[i]) || cur_g_data[i] == GT_NOT_CALL)
                            cur_g_data[i] = bcf_gt_missing;
                        else
                            cur_g_data[i] = bcf_gt_phased(2);
                    }

                    gt_data[i] = bcf_gt_allele(cur_g_data[i]);
                }
                addVariant(gt_data, ncur_g_data, desc);

            }
            tmpi += vcf_record->n_allele - 1;
                // temp_count--;
            bcf_clear(vcf_record);
        }
    }
    cur_pos = desc.pos;
}
// ***************************************************************************************************************************************
void CompressionReader::addVariant(int *gt_data, int ngt_data, variant_desc_t &desc)
{
    // for (int i = 0; i < ngt_data; i++)
    //     cout<<gt_data[i]<<" ";
    // cout<<endl;
    // for (int i = 0; i < ngt_data; i++)
    //     if (gt_data[i] == 0 || gt_data[i] == 1)
    //     {
    //         cout<<"0"<<" ";
    //     }
    //     else // if(bcf_gt_is_missing(gt_arr[i]) || bcf_gt_allele(gt_arr[i]) == 2)
    //     {
    //         cout<<"1"<<" ";
    //     }
    // cout<<endl;
    // for (int i = 0; i < ngt_data; i++)
    // {

    //     if (gt_data[i] == 1 || gt_data[i] == 2)
    //     {
    //         cout<<"1"<<" ";
    //     }
    //     else // 0
    //     {
    //         cout<<"0"<<" ";
    //     }
    // }
    // cout<<endl;
    desc.chrom = vcf_hdr->id[BCF_DT_CTG][vcf_record->rid].key; // CHROM

    desc.pos = vcf_record->pos + 1;                      // POS
    desc.id = vcf_record->d.id ? vcf_record->d.id : "."; // ID
                                                      //    desc.qual.erase(); // QUAL
    if (bcf_float_is_missing(vcf_record->qual))
        desc.qual = ".";
    else
    {
        kstring_t s = {0, 0, 0};
        kputd(vcf_record->qual, &s);
        //      desc.qual += s.s;
        desc.qual = s.s;
        free(s.s);
    }
    if(!vec_read_fixed_fields)
        chunks_min_pos.emplace_back(desc.pos);
    
    if (start_flag)
    {
        cur_chrom = desc.chrom;
        start_flag = false;
    }

    if (desc.chrom == cur_chrom)
    {
        for (int i = 0; i < ngt_data; i++)
        {
            if (gt_data[i] == 0 || gt_data[i] == 1)
            {
                bv.PutBit(0);
            }
            else // if(bcf_gt_is_missing(gt_arr[i]) || bcf_gt_allele(gt_arr[i]) == 2)
            {
                bv.PutBit(1);
            }
        }

        bv.FlushPartialWordBuffer();

        // Set vector with less significant bits of dibits
        for (int i = 0; i < ngt_data; i++)
        {

            if (gt_data[i] == 1 || gt_data[i] == 2)
            {
                bv.PutBit(1);
            }
            else // 0
            {
                bv.PutBit(0);
            }
        }
        bv.FlushPartialWordBuffer();

        v_vcf_data_compress.emplace_back(desc);
        vec_read_fixed_fields++;
        if( vec_read_fixed_fields == no_fixed_fields)
        {
            actual_variants.emplace_back(no_actual_variants);
            no_actual_variants = 0;
            vec_read_fixed_fields = 0;
        }
        vec_read_in_block += 2; // Two vectors added

        if (vec_read_in_block == no_vec_in_block) // Insert complete block into queue of blocks
        {

            bv.TakeOwnership();

            Gt_queue->Push(block_id, bv.mem_buffer, vec_read_in_block, v_vcf_data_compress);
            v_vcf_data_compress.clear();
            no_chrom_num++;
            no_vec = no_vec + vec_read_in_block;
            block_id++;
            bv.Close();
            bv.Create(block_max_size);
            vec_read_in_block = 0;
        }
        no_chrom = cur_chrom;
    }
    else
    {
        if( vec_read_fixed_fields)
        {
            actual_variants.emplace_back(no_actual_variants-1);
            no_actual_variants = 1;
            vec_read_fixed_fields = 0;
        }
        if(!vec_read_fixed_fields)
            chunks_min_pos.emplace_back(desc.pos);
        vec_read_fixed_fields ++;    
        cur_chrom = desc.chrom;
        if (vec_read_in_block)
        {
            
            bv.TakeOwnership();
            Gt_queue->Push(block_id, bv.mem_buffer, vec_read_in_block, v_vcf_data_compress);
            v_vcf_data_compress.clear();
            no_chrom_num++;
            no_vec = no_vec + vec_read_in_block;
            block_id++;
            bv.Close();
            bv.Create(block_max_size);
            vec_read_in_block = 0;
        }
        for (int i = 0; i < ngt_data; i++)
        {

            if (gt_data[i] == 0 || gt_data[i] == 1)
            {
                bv.PutBit(0);
            }
            else // if(bcf_gt_is_missing(gt_arr[i]) || bcf_gt_allele(gt_arr[i]) == 2)
            {
                bv.PutBit(1);
            }
        }
        bv.FlushPartialWordBuffer();
        // Set vector with less significant bits of dibits
        for (int i = 0; i < ngt_data; i++)
        {
            if (gt_data[i] == 1 || gt_data[i] == 2)
            {
                bv.PutBit(1);
            }
            else // 0
            {
                bv.PutBit(0);
            }
        }
        bv.FlushPartialWordBuffer();
        vec_read_in_block += 2; // Two vectors added
        v_vcf_data_compress.emplace_back(desc);
        where_chrom.emplace_back(make_pair(no_chrom, no_chrom_num));
    }
}


// ***************************************************************************************************************************************
uint32_t CompressionReader::setNoVecBlock(GSC_Params &params)
{
    params.var_in_block = no_samples * params.ploidy;

    if(params.var_in_block < 1024)
    {
        chunk_size = CHUNK_SIZE1;
    }
    else if(params.var_in_block < 4096 )
    {
        chunk_size = CHUNK_SIZE2;
    }
    else if(params.var_in_block < 8192)
    {
        chunk_size = CHUNK_SIZE3;
    }
    else
    {
        chunk_size = params.var_in_block;
    }
    params.no_blocks = chunk_size/params.var_in_block;

    no_fixed_fields = params.no_blocks*params.var_in_block;

    if (params.task_mode == task_mode_t::mcompress)
    {
        no_vec_in_block = params.var_in_block * 2;
        params.vec_len = params.var_in_block / 8 + ((params.var_in_block % 8) ? 1 : 0);
        params.n_samples = no_samples;
    }
    return 0;
}
void CompressionReader::CloseFiles(){
    // cout<<"Closing files"<<endl;
    if( vec_read_fixed_fields)
    {
        actual_variants.emplace_back(no_actual_variants);
        no_actual_variants = 0;
        vec_read_fixed_fields = 0;
    } 
    if (vec_read_in_block)
    {
        bv.TakeOwnership();
        // Gt_queue->Push(block_id, bv.mem_buffer, vec_read_in_block,v_vcf_data_compress,chrom_flag::none);
        Gt_queue->Push(block_id, bv.mem_buffer, vec_read_in_block, v_vcf_data_compress);
        // v_vcf_data_compress.clear();
        no_chrom_num++;
        block_id++;
        no_vec = no_vec + vec_read_in_block;
        vec_read_in_block = 0;
        bv.Close();
    }
    
    v_vcf_data_compress.emplace_back(variant_desc_t());
    Gt_queue->Push(block_id, nullptr, 0 ,v_vcf_data_compress);
    
    block_id = 0;
    Gt_queue->Complete();
    

    where_chrom.emplace_back(make_pair(no_chrom, no_chrom_num));

    no_chrom_num = 0;
    if(compress_all){
	    for (uint32_t i = 0; i < no_keys; ++i)
	    {
        
		    v_o_buf[i].GetBuffer(v_size, v_data);
        
		    auto part_id = file_handle2->AddPartPrepare(v_buf_ids_size[i]);
		    file_handle2->AddPartPrepare(v_buf_ids_data[i]);

		    SPackage pck(i, v_buf_ids_size[i], v_buf_ids_data[i], part_id, v_size, v_data);

		    part_queue->Push(pck);

	    }

	    part_queue->Complete();
    }
	// file_handle2->Close();
	
    
}
// ***************************************************************************************************************************************
vector<int> CompressionReader::topo_sort(unordered_map<int, unordered_set<int>> &graph,unordered_map<int, int> inDegree) {

    vector<int> result;
    queue<int> q;
    for (const auto& p : inDegree) {
        if (p.second == 0) {
            q.push(p.first);
        }
    }
    while (!q.empty()) {
        if(q.size() > 1)
            field_order_flag = true;
        int cur = q.front();
        q.pop();
        result.push_back(cur);
        for (int next : graph[cur]) {
            
            if (--inDegree[next] == 0) {
                q.push(next);
            }
        }
    }

    // Step 3: Check if there's a cycle in the graph
    if (result.size() != graph.size()) {
        return vector<int>();
    }

    return result;
}