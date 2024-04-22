#include "compressor.h"
#include <iostream>


// **************************************************************************************************
//write the compressed file
bool Compressor::writeCompressFlie()
{
    std::cerr << "The BSC algorithm was used to compress the genotype part" << endl;
    // cout << "Gt_index_original_size:" << toal_all_size<< endl;

    uint32_t curr_no_blocks = 0;
    uint32_t where_chrom_size = static_cast<uint32_t>(where_chrom.size());
    uint32_t chunks_streams_size = static_cast<uint32_t>(chunks_streams.size());
    fwrite(&chunks_streams_size, sizeof(uint32_t), 1, comp);

    for(auto cur_chunk:chunks_streams)
    {
        fwrite(&cur_chunk.second.cur_chunk_actual_pos, sizeof(uint32_t), 1, comp);
        fwrite(&cur_chunk.second.offset, sizeof(size_t), 1, comp);
    }
    fwrite(&params.ploidy, sizeof(uint8_t), 1, comp);

    fwrite(&params.vec_len, sizeof(params.vec_len), 1, comp);

    fwrite(&no_vec, sizeof(no_vec), 1, comp);

    fwrite(&copy_no, sizeof(copy_no), 1, comp);
    fwrite(&used_bits_cp, sizeof(char), 1, comp);
    fwrite(&bm_comp_copy_orgl_id.mem_buffer_pos, sizeof(int32_t), 1, comp);
    fwrite(bm_comp_copy_orgl_id.mem_buffer, 1, bm_comp_copy_orgl_id.mem_buffer_pos, comp);
    // fwrite(&no_blocks, sizeof(no_blocks), 1, comp);
    // fwrite(&params.max_no_vec_in_block, sizeof(params.max_no_vec_in_block), 1, comp);
    fwrite(&params.n_samples, sizeof(params.n_samples), 1, comp);

    uint32_t chunks_min_pos_size = static_cast<uint32_t>(chunks_min_pos.size());
    fwrite(&chunks_min_pos_size,sizeof(uint32_t) , 1, comp);
    fwrite(&chunks_min_pos[0], sizeof(int64_t), chunks_min_pos_size, comp);

    fwrite(&where_chrom_size, sizeof(uint32_t), 1, comp);
    // fwrite(&all_indexes_pos[0], all_indexes_pos.size() * sizeof(uint64_t), 1, comp);
    for (const auto &elem : where_chrom)
    {
        size_t chrom_size = elem.first.size();
        fwrite(&chrom_size, sizeof(size_t), 1, comp);
        fwrite(elem.first.data(), sizeof(char), chrom_size, comp);
        fwrite(&elem.second, sizeof(int), 1, comp);
    }
    uint32_t vint_last_perm_size = static_cast<uint32_t>(vint_last_perm.size());
    fwrite(&vint_last_perm_size, sizeof(uint32_t), 1, comp);
    for (const auto &data : vint_last_perm)
    {

        fwrite(&data.first, sizeof(uint32_t), 1, comp);
        uint32_t data_size = static_cast<uint32_t>(data.second.size());
        fwrite(&data_size, sizeof(uint32_t), 1, comp);
        fwrite(data.second.data(), sizeof(uint8_t), data_size, comp);
    }

    bm_comp_copy_orgl_id.Close();
    Meta_comp_size += comp_v_header.size();
    uint32_t  comp_size = static_cast<uint32_t>(comp_v_header.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    fwrite(comp_v_header.data(), 1, comp_v_header.size(), comp);

    Meta_comp_size += comp_v_samples.size();
    comp_size = static_cast<uint32_t>(comp_v_samples.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    fwrite(comp_v_samples.data(), 1, comp_v_samples.size(), comp);
    
    uint64_t size;
    uint8_t *temp_buffer = NULL;
    while (fread(&size, sizeof(size), 1, temp_file) == 1) { 
        chunks_streams[curr_no_blocks+1].offset = ftell(comp);
        
        temp_buffer = (uint8_t*)realloc(temp_buffer, size);
        if (!temp_buffer) {
            perror("Memory allocation failed");
            fclose(temp_file);
            return 1;
        }

        if (fread(temp_buffer, 1, size, temp_file) != size) { 
            perror("Error reading data block");
            free(temp_buffer);
            fclose(temp_file);
            return 1;
        }
        fwrite(temp_buffer, 1, size, comp);
        curr_no_blocks++;
    }

    // while (comp_sort_block_queue.Pop(fixed_field_block_id,fixed_field_block_io))
    // {
    //     // std::cerr<<"fixed_fixed_field_block_id: "<<fixed_fixed_field_block_id<<":"<<sort_fixed_field_block_id<<endl;
    //     assert(fixed_field_block_id == curr_no_blocks);
    //     chunks_streams[curr_no_blocks+1].offset = ftell(comp);
    //     curr_no_blocks++;
        
    //     fwrite(&fixed_field_block_io.no_variants, sizeof(uint32_t), 1, comp);
    //     CHORM_comp_size += fixed_field_block_io.chrom.size();
    //     comp_size = static_cast<uint32_t>(fixed_field_block_io.chrom.size());
    //     fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    //     fwrite(fixed_field_block_io.chrom.data(), 1, fixed_field_block_io.chrom.size(), comp);
        
    //     POS_comp_size += fixed_field_block_io.pos.size();
    //     comp_size = static_cast<uint32_t>(fixed_field_block_io.pos.size());
    //     fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    //     fwrite(fixed_field_block_io.pos.data(), 1, fixed_field_block_io.pos.size(), comp);

    //     ID_comp_size += fixed_field_block_io.id.size();
    //     comp_size = static_cast<uint32_t>(fixed_field_block_io.id.size());
    //     fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    //     fwrite(fixed_field_block_io.id.data(), 1, fixed_field_block_io.id.size(), comp);


    //     REF_comp_size += fixed_field_block_io.ref.size();
    //     comp_size = static_cast<uint32_t>(fixed_field_block_io.ref.size());
    //     fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    //     fwrite(fixed_field_block_io.ref.data(), 1, fixed_field_block_io.ref.size(), comp);

        
    //     ALT_comp_size += fixed_field_block_io.alt.size();
    //     comp_size = static_cast<uint32_t>(fixed_field_block_io.alt.size());
    //     fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    //     fwrite(fixed_field_block_io.alt.data(), 1, fixed_field_block_io.alt.size(), comp);

        
    //     QUAL_comp_size += fixed_field_block_io.qual.size();
    //     comp_size = static_cast<uint32_t>(fixed_field_block_io.qual.size());
    //     fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    //     fwrite(fixed_field_block_io.qual.data(), 1, fixed_field_block_io.qual.size(), comp);

    //     GT_comp_size += fixed_field_block_io.gt_block.size();
    //     comp_size = static_cast<uint32_t>(fixed_field_block_io.gt_block.size());
    //     fwrite(&comp_size, sizeof(uint32_t), 1, comp);
    //     fwrite(fixed_field_block_io.gt_block.data(), 1, fixed_field_block_io.gt_block.size(), comp);
    // }
    // std::cerr << "Meta_comp_size:    " << Meta_comp_size <<"\tByte"<<endl;
    // std::cerr << "CHORM_comp_size:   " << CHORM_comp_size <<"\tByte"<< endl;
    // std::cerr << "POS_comp_size:     " << POS_comp_size <<"\tByte"<< endl;
    // std::cerr << "ID_comp_size:      " << ID_comp_size <<"\tByte"<< endl;
    // std::cerr << "REF_comp_size:     " << REF_comp_size <<"\tByte"<< endl;
    // std::cerr << "ALT_comp_size:     " << ALT_comp_size <<"\tByte"<< endl;
    // std::cerr << "QUAL_comp_size:    " << QUAL_comp_size <<"\tByte"<< endl;
    // cout << "GT_index_comp_size:      " << GT_comp_size <<"\tByte"<< endl;
    free(temp_buffer);
    fclose(temp_file);
    if (remove(temp_file1_fname) != 0) {
        perror("Error deleting temp1 file");
        return EXIT_FAILURE;
    }
    // remove(temp_file1_fname);
    other_fields_offset = ftell(comp);

	FILE *other_f = fopen(temp_file2_fname.c_str(),"rb"); 
    
    if (other_f) {
     
        const size_t buffer_size = 1024; 
        char buffer[buffer_size];
        size_t bytes_read;
        while ((bytes_read = fread(buffer, 1, buffer_size, other_f)) > 0) {
            fwrite(buffer, 1, bytes_read, comp);
        }
        fclose(other_f);
        other_f = nullptr;
        if (remove(temp_file2_fname.c_str()) != 0) {
            perror("Error deleting temp2 file");
            return EXIT_FAILURE;
        }

        mode_type = true;
    }   
    
    sdsl_offset = ftell(comp);
    
    // sdsl_offset = end_offset;
    fseek(comp, 0, SEEK_SET);
    fwrite(&mode_type, sizeof(mode_type), 1, comp);
    fwrite(&other_fields_offset, sizeof(other_fields_offset), 1, comp);
    fwrite(&sdsl_offset, sizeof(sdsl_offset), 1, comp);
    std::cerr<<int(mode_type) <<" " << other_fields_offset <<" " << sdsl_offset <<endl;
    fseek(comp, 21, SEEK_SET);
    for(auto cur_chunk:chunks_streams)
    {
        fwrite(&cur_chunk.second.cur_chunk_actual_pos, sizeof(uint32_t), 1, comp);
        fwrite(&cur_chunk.second.offset, sizeof(size_t), 1, comp);
    }
    fseek(comp, 0, SEEK_END);

    if (comp && comp != stdout)
    {
        fclose(comp);
        comp = nullptr;
    }
    if (is_stdout) {
     
        for (int v = 0; v < 2; v++) {
            sdsl::rrr_vector<> rrr_bit_vector(zeros_only_bit_vector[v]);
            rrr_bit_vector.serialize(std::cout);
            sdsl::util::clear(rrr_bit_vector);
        }
        for (int v = 0; v < 2; v++) {
            sdsl::rrr_vector<> rrr_bit_vector(copy_bit_vector[v]);
            rrr_bit_vector.serialize(std::cout);
            sdsl::util::clear(rrr_bit_vector);
        }
    } else {
   
        sdsl::osfstream out(fname, std::ios::binary | std::ios::app);
        if (!out) {
            if (sdsl::util::verbose) {
                std::cerr << "ERROR: store_to_file not successful for: `" << fname << "`" << std::endl;
            }
            exit(1);
        }
        sdsl::rrr_vector<> rrr_bit_vector[5];
        for (int v = 0; v < 2; v++) {
            sdsl::rrr_vector<> rrr_bit_vector(zeros_only_bit_vector[v]);
            rrr_bit_vector.serialize(out);
            sdsl::util::clear(rrr_bit_vector);
        }
        for (int v = 0; v < 2; v++) {
            rrr_bit_vector[v] = sdsl::rrr_vector<>(copy_bit_vector[v]);
            rrr_bit_vector[v].serialize(out);
            sdsl::util::clear(rrr_bit_vector[v]);
        }

        out.close();
    }

    if (sdsl::util::verbose) {
        std::cerr << "INFO: store_to_file: `" << fname << "`" << std::endl;
    }
    // sdsl::osfstream out(fname, std::ios::binary | std::ios::app);
    // if (!out)
    // {
    //     if (sdsl::util::verbose)
    //     {
    //     std::cerr << "ERROR: store_to_file not successful for: `" << fname << "`" << std::endl;
    //     }
    //     exit(1);
    // }
    // sdsl::rrr_vector<> rrr_bit_vector[5];
    
    // for (int v = 0; v < 2; v++)
    // {
    //     rrr_bit_vector[v] = sdsl::rrr_vector<>(zeros_only_bit_vector[v]);
    //     rrr_bit_vector[v].serialize(out);
    //     sdsl::util::clear(rrr_bit_vector[v]);
    // }
    // for (int v = 0; v < 2; v++)
    // {
    //     rrr_bit_vector[v + 2] = sdsl::rrr_vector<>(copy_bit_vector[v]);
    //     rrr_bit_vector[v + 2].serialize(out);
    //     sdsl::util::clear(rrr_bit_vector[v + 2]);
    // }

    // out.close();
    
    // if (sdsl::util::verbose)
    // {
    // std::cerr << "INFO: store_to_file: `" << fname << "`" << std::endl;
    // }

    
    // std::cerr << "genotype compress file (" << params.out_file_name + ".gsc"<< ") created." << std::endl;    


    return 0;
}
//open file for writing
// *******************************************************************************************************************
bool Compressor::OpenForWriting(const string &out_file_name)
{
    if(out_file_name != "-"){
        fname = out_file_name;
        // fname = (char *)malloc(strlen(out_file_name.c_str()) + 5);
        // snprintf(fname, strlen(out_file_name.c_str()) + 5, "%s.gsc", out_file_name.c_str());
        comp = fopen(fname.c_str(), "wb");
        if (!comp)
        {

            std::cerr << "ERROR: storing archive not successful for: `" << fname << "`" << std::endl;
            exit(1);
        }
    }else {
        
        comp = stdout;
        is_stdout = true;
        
    }
    if (setvbuf(comp, nullptr, _IOFBF, 64 << 20) != 0) {
        std::cerr << "ERROR: Buffer setup failed for: `" << fname << "`" << std::endl;
        // if (fname != nullptr) {
        //     free(fname);
        // }
        exit(1);
    }


    mode_type = false;

    fwrite(&mode_type, sizeof(mode_type), 1, comp);
    
    other_fields_offset = ftell(comp) + sizeof(uint64_t);
    fwrite(&other_fields_offset, sizeof(other_fields_offset), 1, comp);

    sdsl_offset = ftell(comp) + sizeof(uint64_t);
    
    if (fwrite(&sdsl_offset, sizeof(sdsl_offset), 1, comp) != 1) {
        std::cerr << "ERROR: Write operation failed for: `" << fname << "`" << std::endl;
    }

    
    // setvbuf(comp, nullptr, _IOFBF, 64 << 20);
    // sdsl_offset = ftell(comp) + sizeof(uint64_t);
    
    // fwrite(&sdsl_offset, sizeof(sdsl_offset), 1, comp);

    int id = (int) chunks_streams.size();
                                                
	chunks_streams[id] = chunk_stream(0,0);

    
    return true;
}
//*******************************************************************************************************************
bool Compressor::OpenTempFile(const string &out_file_name)
{
    temp_file1_fname = (char *)malloc(strlen(out_file_name.c_str()) + 5);

    snprintf(temp_file1_fname, strlen(out_file_name.c_str()) + 5, "%s.temp", out_file_name.c_str());

    temp_file = fopen(temp_file1_fname, "wb");
    if (!temp_file)
    {

        std::cerr << "ERROR: storing archive not successful for: `" << temp_file1_fname << "`" << std::endl;
        exit(1);
    }
    
    setvbuf(temp_file, nullptr, _IOFBF, 64 << 20);


    if(params.compress_mode == compress_mode_t::lossless_mode){
        if (file_handle2)
		    delete file_handle2;
	    file_handle2 = new File_Handle_2(false);

        temp_file2_fname = out_file_name + ".com_tmp_gsc";

        if (!file_handle2->Open(temp_file2_fname))
	    {
		    cerr << "Cannot open " << temp_file2_fname << "\n";
		    return false;
	    }
    }
    return true;
}
// *******************************************************************************************************************
//compess pragma entry
bool Compressor::CompressProcess()
{
    //it is important to initialize the library before using it
    CBSCWrapper::InitLibrary(p_bsc_features);
    
    MyBarrier  my_barrier(3);

    unique_ptr<CompressionReader> compression_reader(new CompressionReader(params));

    if (!compression_reader->OpenForReading(params.in_file_name))
	{
		cerr << "Cannot open: " << params.in_file_name << endl;
		return false;
	}
    if (!OpenForWriting(params.out_file_name))
		return false;
    
    if (!OpenTempFile(params.out_file_name))
		return false;
        
    string header;
    
    vector<string> v_samples;

    compression_reader->GetHeader(header);
    
    uint32_t no_samples = compression_reader->GetSamples(v_samples);
    // for(int i =0;i<v_samples.size();i++)
    //     std::cerr<<v_samples[i]<<endl;
    std::cerr << "no_samples:" << no_samples << endl;

    if (!no_samples)
    {
        std::cerr << "The number of genotype samples is zero and cannot be compressed!\n";
        return false;
    }
    
    //compress meta data
    compress_meta(v_samples , header);


    compression_reader->setNoVecBlock(params);

    std::cerr<<"no_gt_threads:"<<params.no_gt_threads<<endl;
    GtBlockQueue inGtBlockQueue(max((int)(params.no_blocks*params.no_gt_threads),8));

    // VarBlockQueue<fixed_fixed_field_block> inVarBlockQueue(max((int)params.no_threads * 2, 8));
    VarBlockQueue<fixed_field_block>  sortVarBlockQueue(max((int)params.no_threads * 2, 8));
    compression_reader->setQueue(&inGtBlockQueue);

    PartQueue<SPackage> part_queue(max((int)params.no_threads * 2, 8)); 

    if(params.compress_mode == compress_mode_t::lossless_mode){

        compression_reader->setPartQueue(&part_queue); 

        compression_reader->InitVarinats(file_handle2);
    
        compression_reader->GetOtherField(keys,no_keys,key_gt_id);
        
        InitCompressParams();
        
        part_compress_thread.reserve(params.no_threads); 
        
        for (uint32_t i = 0; i < params.no_threads; ++i)
        {
            part_compress_thread.emplace_back(thread([&]() {

                SPackage pck;
		        vector<uint8_t> v_compressed;
		        vector<uint8_t> v_tmp;
		
		        auto fo = [this](SPackage& pck)->bool {return check_coder_compressor(pck); };

		        while (true)
		        {
                    
			        if (!part_queue.Pop<SPackage>(pck, fo))
				        break;
                    compress_other_fileds(pck, v_compressed, v_tmp);
                
		        }      


            }));
        }
    }    
    block_size = no_samples * params.ploidy * 2;
    
    unique_ptr<thread> compress_thread(new thread([&] {
        fixed_field_block fixed_field_block_process;
        while (true)
        {

            if (!sortVarBlockQueue.Pop(fixed_field_block_id,fixed_field_block_process)){
                break;
            }
            
            compressFixedFields(fixed_field_block_process);

        }

    }));

    // create multiple threads to handle individual blocks
    block_process_thread.reserve(params.no_gt_threads);
    string prev_chrom = "";
    int chunk_id = 0;
    for(uint32_t i = 0; i < params.no_gt_threads; ++i)
        block_process_thread.emplace_back(thread([&]() {
            
            int block_id = 0;
            unsigned long num_rows;
            unsigned char *data = nullptr;
            vector<variant_desc_t> v_vcf_data_io;             
            vector<uint32_t> origin_of_copy;
            origin_of_copy.reserve(no_variants_in_buf);
            vector<uint8_t> samples_indexes;   //Index of the location where the block storing 1 is stored.
            vector<uint32_t> perm;
            perm.clear();
            perm.resize(no_samples * params.ploidy, 0);
            for (size_t i_p = 0; i_p < perm.size(); i_p++)
                perm[i_p] = i_p;
            BlockProcess block_process(params);
                                    
            while (true)
            {
                if (!inGtBlockQueue.Pop(block_id, data, num_rows, v_vcf_data_io)){
                    break;
                }

                vector<bool> zeros_only(num_rows, false);
                vector<bool> copies(num_rows, false);                           
                block_process.SetCurBlock(num_rows, data);
                //Gets the sparse encoding for each block        
                          
                if (num_rows){    
                    // if(num_rows == block_size)  
                        block_process.ProcessSquareBlock(perm, zeros_only, copies, origin_of_copy,samples_indexes,true);
                    // else
                    //     block_process.ProcessLastBlock(zeros_only, copies, origin_of_copy,samples_indexes);

                    if(num_rows == block_size)

                        block_process.ProcessVariant(perm , v_vcf_data_io);

                }
                
                if(data != nullptr)
                    delete[] data;
                //Gets the sparse encoding for each block
                lock_gt_block_process(block_id);
                {
                    // if(!cur_block_id)
                    //     prev_chrom = v_vcf_data_io[0].chrom;
                    // if(num_rows != block_size)                    
                    //     block_process.ProcessLastPerm(perm,vint_last_perm); 

                   
                    
                    if(prev_chrom != v_vcf_data_io[0].chrom){

                        prev_chrom = v_vcf_data_io[0].chrom;
                        if(no_curr_chrom_block){

                            sortVarBlockQueue.Push(chunk_id,fixed_field_block_io);
                            // compressFixedFields(fixed_field_block_io);
                            toal_all_size += fixed_field_block_io.gt_block.size();
                            int id = (int) chunks_streams.size();                
	                        chunks_streams[id] = chunk_stream(cur_chunk_actual_pos,0);
                            
                            no_curr_chrom_block = 0;
                            prev_pos = 0;
                            chunk_id++;
                            fixed_field_block_io.Clear();
                        }
                        if(num_rows){

                            cur_chunk_actual_pos += (uint32_t)v_vcf_data_io.size();
                            block_process.addSortFieldBlock(fixed_field_block_io,all_zeros,all_copies,comp_pos_copy,zeros_only, copies, origin_of_copy,samples_indexes,v_vcf_data_io,prev_pos);
                            no_curr_chrom_block++;
                            if(num_rows % block_size){
                                vint_last_perm.emplace(chunk_id,vint_code::EncodeArray(perm));
                            }

                            if(no_curr_chrom_block == params.no_blocks){
                                sortVarBlockQueue.Push(chunk_id,fixed_field_block_io);
                                // compressFixedFields(fixed_field_block_io);
                                toal_all_size += fixed_field_block_io.gt_block.size();
                                int id = (int) chunks_streams.size();               
                                chunks_streams[id] = chunk_stream(cur_chunk_actual_pos,0);
                                no_curr_chrom_block = 0;
                                prev_pos = 0;
                                chunk_id++;
                                fixed_field_block_io.Clear();
                            }
                        }
                        else{

                            sortVarBlockQueue.Complete();
                        }
                            
                    }
                    else{

                        cur_chunk_actual_pos += (uint32_t)v_vcf_data_io.size();
                        block_process.addSortFieldBlock(fixed_field_block_io,all_zeros,all_copies,comp_pos_copy,zeros_only, copies, origin_of_copy,samples_indexes,v_vcf_data_io,prev_pos);
                        no_curr_chrom_block++;

                        if(num_rows % block_size){
                            vint_last_perm.emplace(chunk_id,vint_code::EncodeArray(perm));
                        }

                        if(no_curr_chrom_block == params.no_blocks){
                            sortVarBlockQueue.Push(chunk_id,fixed_field_block_io);
                            // compressFixedFields(fixed_field_block_io);
                            toal_all_size += fixed_field_block_io.gt_block.size();
                            int id = (int) chunks_streams.size();               
	                        chunks_streams[id] = chunk_stream(cur_chunk_actual_pos,0);
                            no_curr_chrom_block = 0;
                            prev_pos = 0;
                            chunk_id++;
                            fixed_field_block_io.Clear();
                        }
                        

                    }   
                    

           
                }         
                unlock_gt_block_process();                   
                
            }                            
                                    
        })); 

    
    if (!compression_reader->ProcessInVCF())
        return false;

    no_vec = compression_reader->getNoVec();

    // obtain the variation description

    compression_reader->GetWhereChrom(where_chrom,chunks_min_pos);
    
    if(params.compress_mode == compress_mode_t::lossless_mode){
        
        for (uint32_t i = 0; i < params.no_threads; ++i)
		    part_compress_thread[i].join();
        
        auto stream_id = file_handle2->RegisterStream("part2_params");
        vector<uint8_t> v_desc;
        vector<uint32_t> actual_variants = compression_reader->GetActualVariants();
        compression_reader->UpdateKeys(keys);
        append(v_desc, static_cast<uint32_t>(actual_variants.size()));
        
        for (uint32_t i = 0; i < actual_variants.size(); ++i)
        {
            append(v_desc, actual_variants[i]);
            
        }
        append(v_desc, no_keys);
        append(v_desc, key_gt_id);
	    for (uint32_t i = 0; i < no_keys; ++i)
	    {
            // std::cerr<<"key_id:"<<keys[i].key_id<<":"<<keys[i].actual_field_id<<endl;
		    append(v_desc, keys[i].key_id);
            append(v_desc, keys[i].actual_field_id);
		    append(v_desc, static_cast<uint64_t>(keys[i].keys_type));
		    append(v_desc, keys[i].type);
	    }
        vector<uint8_t> v_desc_compressed;
        CBSCWrapper bsc;
		bsc.InitCompress(p_bsc_fixed_fields);
		bsc.Compress(v_desc, v_desc_compressed);
        // zstd::zstd_compress(v_desc, v_desc_compressed);
        file_handle2->AddParamsPart(stream_id,v_desc_compressed);
        file_handle2->Close();
    }

     for(uint32_t i = 0; i < params.no_gt_threads; ++i)
        block_process_thread[i].join();
    compress_thread->join();

    if(temp_file)
        fclose(temp_file);
    temp_file = fopen(temp_file1_fname, "rb");
    std::cerr << "Complete and process the chunking." << std::endl;


    compressReplicatedRow();

    writeCompressFlie();
    
    return true;
}
//process the all_zeros and all_copies
// **************************************************************************************************
void Compressor::compressReplicatedRow()
{
    zeros_only_bit_vector[0] = sdsl::bit_vector(no_vec / 2 + no_vec % 2, 0);
    zeros_only_bit_vector[1] = sdsl::bit_vector(no_vec / 2 + no_vec % 2, 0);
    copy_bit_vector[0] = sdsl::bit_vector(no_vec / 2 + no_vec % 2, 0);
    copy_bit_vector[1] = sdsl::bit_vector(no_vec / 2 + no_vec % 2, 0);
   
    // comp_pos_copy = new uint32_t[all_copy_num]();
 
    unique = sdsl::bit_vector(no_vec, 0);
    
    for (uint64_t i = 0; i < all_zeros.size(); i++)
    {
        zeros_only_bit_vector[i % 2][i / 2] = all_zeros[i];
        copy_bit_vector[i % 2][i / 2] = all_copies[i];
    }
    comp_pos_copy.shrink_to_fit();
    copy_no = comp_pos_copy.size();
    
    // std::cerr << "no_vec:" << no_vec << endl;
    // rank_copy_bit_vector[0] = sdsl::rank_support_v5<>(&copy_bit_vector[0]);
    // rank_copy_bit_vector[1] = sdsl::rank_support_v5<>(&copy_bit_vector[1]);
    // rank_zeros_only_vector[0] = sdsl::rank_support_v5<>(&zeros_only_bit_vector[0]);
    // rank_zeros_only_vector[1] = sdsl::rank_support_v5<>(&zeros_only_bit_vector[1]);

    uint32_t n_copies = 0;

    uint64_t i = 0, zeros_no = 0, uqq = 0;

    for (i = 0; i < no_vec; i++)
    {
        if (zeros_only_bit_vector[i % 2][i / 2])
        {
            zeros_no++;
            continue;
        }
        else if (!copy_bit_vector[i % 2][i / 2])
        {
            unique[i] = 1;
            uqq++;
        }
        else
            n_copies++;
    }
    rank_unique = sdsl::rank_support_v5<>(&unique);
    uint64_t curr_non_copy_vec_id = 0;
    uint64_t copy_id = 0, origin_unique_id;
    uint32_t max_diff_copy = 0;
    
    for (uint64_t i = 0; i < no_vec; i++)
    {
        if (zeros_only_bit_vector[i % 2][i / 2])
            continue;
        if (!copy_bit_vector[i % 2][i / 2])
        {
            curr_non_copy_vec_id++;
            
            continue;
        }
        // Here curr_non_copy_vec_id is a fake curr_non_copy_vec_id (it is ud of the next non_copy_vec_id)
        
        origin_unique_id = rank_unique(comp_pos_copy[copy_id]);
        
        // Store difference -1, to not waste one value
        comp_pos_copy[copy_id] = curr_non_copy_vec_id - origin_unique_id - 1;
        
        if (comp_pos_copy[copy_id] > max_diff_copy)
        {
            max_diff_copy = comp_pos_copy[copy_id];
        }
        copy_id++;
    }
    used_bits_cp = bits_used(max_diff_copy);
    bm_comp_copy_orgl_id.Create(copy_no * 4);

    for (i = 0; i < copy_no; i++)
    {
        // std::cerr<<comp_pos_copy[i]<<endl;
        bm_comp_copy_orgl_id.PutBits(comp_pos_copy[i], (int32_t)used_bits_cp);
    }
    bm_comp_copy_orgl_id.FlushPartialWordBuffer();

    sdsl::util::clear(rank_unique);
}

//get the bits of the number
// **************************************************************************************************
char Compressor::bits_used(unsigned int n)
{
    char bits = 0;
    while (n)
    {
        n = n >> 1;
        bits++;
    }
    return bits;
}
// ************************************************************************************
void Compressor::lock_coder_compressor(SPackage& pck)
{
	unique_lock<mutex> lck(mtx_v_coder);
	cv_v_coder.wait(lck, [&, this] {
		int sid = pck.key_id;

		return (int) v_coder_part_ids[sid] == pck.part_id;
		});
}

// ************************************************************************************
bool Compressor::check_coder_compressor(SPackage& pck)
{
	unique_lock<mutex> lck(mtx_v_coder);
	int sid = pck.key_id;

	
	return (int) v_coder_part_ids[sid] == pck.part_id;
}

// ************************************************************************************
void Compressor::unlock_coder_compressor(SPackage& pck)
{
	lock_guard<mutex> lck(mtx_v_coder);
	int sid = pck.key_id;

	++v_coder_part_ids[sid];
	cv_v_coder.notify_all();
}
void Compressor::compress_other_fileds(SPackage& pck, vector<uint8_t>& v_compressed, vector<uint8_t>& v_tmp)
{
    
    lock_coder_compressor(pck);
    CBSCWrapper *cbsc_size = v_bsc_size[pck.key_id];
	if (pck.v_data.size())
	{
        
		v_compressed.clear();
        // lzma2::lzma2_compress(pck.v_data, v_compressed, 10, 3);
          
        CBSCWrapper *cbsc = v_bsc_data[pck.key_id];
        if(keys[pck.key_id].type != BCF_HT_INT){
            cbsc->Compress(pck.v_data, v_compressed);
            // zstd::zstd_compress(pck.v_data, v_compressed);
        }
        else{
            Encoder(pck.v_data, v_tmp);
            cbsc->Compress(v_tmp, v_compressed);
            // zstd::zstd_compress(v_tmp, v_compressed);
        }
            
        // cbsc->Compress(pck.v_data, v_compressed);
        
		file_handle2->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed);
	}
	else
	{
		v_compressed.clear();

		file_handle2->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed);
        
	}
    
    v_compressed.clear();
    
    // lzma2::lzma2_compress(pck.v_size, v_compressed, 10, 1);
    v_tmp.resize(pck.v_size.size() * 4);

	copy_n((uint8_t*)pck.v_size.data(), v_tmp.size(), v_tmp.data());
    
    cbsc_size->Compress(v_tmp, v_compressed);
    // zstd::zstd_compress(v_tmp, v_compressed);
    
	file_handle2->AddPartComplete(pck.stream_id_size, pck.part_id, v_compressed);

	unlock_coder_compressor(pck);
}
// void Compressor::compress_INT_fileds(SPackage& pck, vector<uint8_t>& v_compressed, vector<uint8_t>& v_tmp)
// {
    
//     vector<uint32_t> vec_temp1;
//     vector<uint32_t> vec_temp2;
//     lock_coder_compressor(pck);
// 	if (pck.v_data.size())
// 	{
// 		vec_temp1.resize(pck.v_data.size()/4);
//         v_compressed.clear();
//         for (size_t i = 0; i < vec_temp1.size(); i++) {
//             vec_temp1[i] = ((uint32_t)pck.v_data[i*4] ) |
//                         ((uint32_t)pck.v_data[i*4+1] << 8) |
//                         ((uint32_t)pck.v_data[i*4+2] << 16) |
//                         ((uint32_t)pck.v_data[i*4+3]<< 24);
            
//         }
        
//         // lzma2::lzma2_compress(pck.v_data, v_compressed, 10, 3);
        
//         FastPForCompress::FastPFor_Compress(vec_temp1, vec_temp2);
        

//         v_compressed.resize(vec_temp2.size()*4);

//         copy_n((uint8_t*)vec_temp2.data(),vec_temp2.size()*4,v_compressed.data());
// 		file_handle2->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed);
// 	}
// 	else
// 	{
// 		v_compressed.clear();

// 		file_handle2->AddPartComplete(pck.stream_id_data, pck.part_id, v_compressed);
// 	}

//     v_compressed.clear();
//     lzma2::lzma2_compress(pck.v_size, v_compressed, 10, 1);
//     // FastPForCompress::FastPFor_Compress(pck.v_size, vec_temp2);
//     // v_compressed.resize(vec_temp2.size()*4);
// 	// copy_n((uint8_t*)vec_temp2.data(),vec_temp2.size()*4,v_compressed.data());
    
// 	file_handle2->AddPartComplete(pck.stream_id_size, pck.part_id, v_compressed);

// 	unlock_coder_compressor(pck);
// }

void Compressor::Encoder(vector<uint8_t>& v_data, vector<uint8_t>& v_tmp)
{
    v_tmp.resize(v_data.size());
    size_t size = v_data.size()/4;
    for (size_t i = 0; i < v_data.size()/4; i++) {
        v_tmp[i] = v_data[i*4];
        v_tmp[i+size] = v_data[i*4+1];
        v_tmp[i+size*2]  = v_data[i*4+2];
        v_tmp[i+size*3]  = v_data[i*4+3];
    }
    // v_data = tmp;
    // uint8_t count = 1;
    // for (size_t i = 1; i < v_data.size(); ++i) {
    //     if(v_data[i] == v_data[i-1] && count < 255)
    //     {
    //         count++;
    //     }
    //     else
    //     {

    //         v_tmp.emplace_back(count);
    //         v_tmp.emplace_back(v_data[i-1]);
    //         count = 1;
    //     }
        
    // }
    // v_tmp.emplace_back(count);
    // v_tmp.emplace_back(v_data[v_data.size() - 1]);

    // uint8_t count = 0;
    // for (size_t i = 0; i < v_data.size(); ++i) {
    //     if (v_data[i] == 0) {
    //         count++;
    //         if(count == 255)
    //         {
    //             v_tmp.emplace_back(0);
    //             v_tmp.emplace_back(count);
    //             count = 0;
    //         }
    //     }
    //     else {
    //         if(count)
    //         {
    //             v_tmp.emplace_back(0);
    //             v_tmp.emplace_back(count);
    //             count = 0;
    //         }
    //         v_tmp.emplace_back(v_data[i]);
    //     }
    // }
    // if(count)
    // {
    //     v_tmp.emplace_back(0);
    //     v_tmp.emplace_back(count);
    // }

    v_tmp.shrink_to_fit();
    // for (size_t i = 0; i < v_data.size(); i++)
    // {
    //     std::cerr<<(int)v_data[i]<<" ";
    // }
    // std::cerr<<endl;
    // std::cerr<<endl;
    // for(size_t i = 0; i < v_tmp.size(); i++)
    // {
    //     std::cerr<<(int)v_tmp[i]<<" ";
    // }
    // std::cerr<<endl;
    // std::cerr<<endl;
}
//compressor meta data
// *******************************************************************************************************************************************
bool Compressor::compress_meta(vector<string> v_samples,const string& v_header)
{

	append_str(all_v_header, v_header);
    
	for (auto &x : v_samples)
		append_str(all_v_samples, x);
        
	for (auto data : {
			make_tuple(ref(all_v_header), ref(comp_v_header), "header"),
			make_tuple(ref(all_v_samples), ref(comp_v_samples), "samples"),
		 })
	{
        CBSCWrapper bsc;

		bsc.InitCompress(p_bsc_meta);
       
		bsc.Compress(get<0>(data), get<1>(data));


        // zstd::zstd_compress(get<0>(data), get<1>(data));

        // LZMACompress::Compress(get<0>(data), get<1>(data), 9);

        // lz4:: lz4_compress(get<0>(data), get<1>(data), 12);
        // BrotliUtils::compressData(get<0>(data), get<1>(data));
        // std::cerr<<get<0>(data).size()<<":"<<get<1>(data).size()<<endl;    
		// lzma2::lzma2_compress(get<0>(data), get<1>(data), get<2>(data), 10);
		// std::cerr << get<2>(data) << " size: " << get<1>(data).size() << endl;

		// fh.WriteUInt(get<1>(data).size(), 4, f);
		// fh.Write(get<1>(data).data(), get<1>(data).size(), f);
	}

	return true;
}
//init bsc compress params
//******************************************************************************************************************************************
void Compressor::InitCompressParams(){
    
    v_coder_part_ids.resize(no_keys, 0);
    v_bsc_data.resize(no_keys);
    v_bsc_size.resize(no_keys);
   
    for (uint32_t i = 0; i < no_keys; ++i)
    {   
        v_bsc_data[i] = new CBSCWrapper();
        v_bsc_size[i] = new CBSCWrapper();
        v_bsc_size[i]->InitCompress(p_bsc_size);
        
        switch (keys[i].type)
		{
		case BCF_HT_FLAG:
			v_bsc_data[i]->InitCompress(p_bsc_flag);
			break;
		case BCF_HT_INT:
            v_bsc_data[i]->InitCompress(p_bsc_int);
			break;
		case BCF_HT_REAL:
			v_bsc_data[i]->InitCompress(p_bsc_real);
			break;
		case BCF_HT_STR:
			v_bsc_data[i]->InitCompress(p_bsc_text);
			break;
		}
    }

 
}
void Compressor::lock_gt_block_process(int &_block_id)
{
	unique_lock<mutex> lck(mtx_gt_block);
	cv_gt_block.wait(lck, [&, this] {return cur_block_id == _block_id;});
}

// ************************************************************************************
bool Compressor::check_gt_block_process(int &_block_id)
{
	unique_lock<mutex> lck(mtx_gt_block);

	return (int) cur_block_id == _block_id;
}

// ************************************************************************************
void Compressor::unlock_gt_block_process()
{
	lock_guard<mutex> lck(mtx_gt_block);
    ++cur_block_id;
	cv_gt_block.notify_all();
}
bool Compressor::compressFixedFields(fixed_field_block &fixed_field_block_io){

    fixed_field_block_compress.no_variants = fixed_field_block_io.no_variants;
    for (auto data : {
        make_tuple(ref(fixed_field_block_io.chrom), ref(fixed_field_block_compress.chrom),  "chrom"),
        make_tuple(ref(fixed_field_block_io.id), ref(fixed_field_block_compress.id),  "id"),
        make_tuple(ref(fixed_field_block_io.alt), ref(fixed_field_block_compress.alt), "alt"),
        make_tuple(ref(fixed_field_block_io.qual), ref(fixed_field_block_compress.qual), "qual"),          
        make_tuple(ref(fixed_field_block_io.pos), ref(fixed_field_block_compress.pos),  "pos"),
        make_tuple(ref(fixed_field_block_io.ref), ref(fixed_field_block_compress.ref), "ref"),
        // make_tuple(ref(fixed_field_block_io.gt_block), ref(fixed_field_block_compress.gt_block), "GT"),
    })
    {         

        //BSC

        CBSCWrapper cbsc;
        cbsc.InitCompress(p_bsc_fixed_fields);
        cbsc.Compress(get<0>(data), get<1>(data));   

        //ZSTD  

        // zstd::zstd_compress(get<0>(data), get<1>(data));

        //LZMA

        // LZMACompress::Compress(get<0>(data), get<1>(data), 9);

        // lz4:: lz4_compress(get<0>(data), get<1>(data), 12);
        // BrotliUtils::compressData(get<0>(data), get<1>(data));

        // std::cerr<<get<0>(data).size()<<":"<<get<1>(data).size()<<endl;    
    }
    if(fixed_field_block_io.gt_block.size() < (2<<20)){
        CBSCWrapper cbsc;
        cbsc.InitCompress(p_bsc_fixed_fields);
        cbsc.Compress(fixed_field_block_io.gt_block, fixed_field_block_compress.gt_block);  
        fixed_field_block_compress.gt_block.emplace_back(0);                  
    }else{
                    // std::cerr<<"zstd"<<endl;
        zstd::zstd_compress(fixed_field_block_io.gt_block,fixed_field_block_compress.gt_block);
        // std::cerr<<fixed_field_block_compress.no_variants<<":"<<fixed_field_block_io.gt_block.size()<<":"<<fixed_field_block_compress.gt_block.size()<<endl;
        fixed_field_block_compress.gt_block.emplace_back(1);
    }
    writeTempFlie(fixed_field_block_compress);
    // comp_sort_block_queue.Push(fixed_field_block_id,fixed_field_block_compress);
    fixed_field_block_compress.Clear();    
    return true;
}

bool Compressor::writeTempFlie(fixed_field_block &fixed_field_block_io){

    uint64_t offset = ftell(temp_file) + sizeof(uint64_t);
    
    uint64_t start_offset = ftell(temp_file);
    size_t comp_size = 0;
    fwrite(&offset, sizeof(offset), 1, temp_file);
    
    fwrite(&fixed_field_block_io.no_variants, sizeof(uint32_t), 1, temp_file);
    comp_size = static_cast<uint32_t>(fixed_field_block_io.chrom.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, temp_file);
    fwrite(fixed_field_block_io.chrom.data(), 1, fixed_field_block_io.chrom.size(), temp_file);
    
    comp_size = static_cast<uint32_t>(fixed_field_block_io.pos.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, temp_file);
    fwrite(fixed_field_block_io.pos.data(), 1, fixed_field_block_io.pos.size(), temp_file);

    comp_size = static_cast<uint32_t>(fixed_field_block_io.id.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, temp_file);
    fwrite(fixed_field_block_io.id.data(), 1, fixed_field_block_io.id.size(), temp_file);

    comp_size = static_cast<uint32_t>(fixed_field_block_io.ref.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, temp_file);
    fwrite(fixed_field_block_io.ref.data(), 1, fixed_field_block_io.ref.size(), temp_file);

    comp_size = static_cast<uint32_t>(fixed_field_block_io.alt.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, temp_file);
    fwrite(fixed_field_block_io.alt.data(), 1, fixed_field_block_io.alt.size(), temp_file);

    comp_size = static_cast<uint32_t>(fixed_field_block_io.qual.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, temp_file);
    fwrite(fixed_field_block_io.qual.data(), 1, fixed_field_block_io.qual.size(), temp_file);

    comp_size = static_cast<uint32_t>(fixed_field_block_io.gt_block.size());
    fwrite(&comp_size, sizeof(uint32_t), 1, temp_file);
    fwrite(fixed_field_block_io.gt_block.data(), 1, fixed_field_block_io.gt_block.size(), temp_file);
    
    offset = ftell(temp_file) - offset;

    fseek(temp_file, start_offset, SEEK_SET);
  
    fwrite(&offset, sizeof(offset), 1, temp_file);
    fseek(temp_file, 0, SEEK_END);

    return  true;

}