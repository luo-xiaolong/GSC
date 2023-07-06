#include "decompression_reader.h"
#include <chrono>

using namespace std::chrono;
using namespace std;

// *******************************************************************************************************************************
void DecompressionReader::out_perm(vector<uint32_t> &perm, vector<variant_desc_t> &v_vcf_data_io)
{

	vector<sblock> temp;
	for (size_t i = 0; i < v_vcf_data_io.size(); i++)
	{
		sblock sb(v_vcf_data_io[i].pos, i, v_vcf_data_io[i].ref);
		temp.emplace_back(sb);
	}
	my_merge_sort(temp, temp.size());
	// my_merge_sort(temp);
	for (size_t i = 0; i < v_vcf_data_io.size(); i++)
	{
		perm[i] = static_cast<uint32_t>(temp[i].p);
	}
	for (size_t i = 0; i < v_vcf_data_io.size(); i++)
	{

		if (atoi(temp[i].s_ref.c_str()))
		{

			temp[i].s_ref = temp[i].s_ref.substr(to_string(atoi(temp[i].s_ref.c_str())).length());
		}
		v_vcf_data_io[i].pos = temp[i].val;
		v_vcf_data_io[i].ref = temp[i].s_ref;
	}
	
}

// *******************************************************************************************************************************
bool DecompressionReader::OpenReading(const string &in_file_name){
	
	CBSCWrapper::InitLibrary(p_bsc_features);
	string fname = in_file_name + ".gti";
	// sdsl vectors
    sdsl::isfstream in(fname, std::ios::binary | std::ios::in);
    if (!in)
    {
        // if (sdsl::util::verbose)
        {
            std::cerr << "Could not load file `" << fname << "`" << std::endl;
        }
        exit(1);
    }
	uint64_t sdsl_offset;
	in.read((char *)&sdsl_offset, sizeof(uint64_t));
	
	in.seekg(sdsl_offset, std::ios::beg);

    rrr_zeros_bit_vector[0].load(in);
    rrr_zeros_bit_vector[1].load(in);
    rrr_copy_bit_vector[0].load(in);
    rrr_copy_bit_vector[1].load(in);
	in.seekg(sizeof(sdsl_offset), std::ios::beg);
    uint64_t FileStartPosition = in.tellg();
	
    in.close();
    if (sdsl::util::verbose)
    {
        std::cerr << "Load file `" << fname << "`" << std::endl;
    }

    // rest of archive
    buf_pos = 0;
    uint64_t arch_size = 0;
#ifndef MMAP
    {
        FILE *comp = fopen(fname.c_str(), "rb");

        if (!comp)
        {
            cout << "Input file (" << fname << ")error\n";
            exit(1);
        }
        arch_size = sdsl_offset - FileStartPosition;
        fseek(comp, FileStartPosition, SEEK_SET);

        buf = new uchar[arch_size];
        fread(buf, sizeof(uchar) * arch_size, 1, comp);

        fclose(comp);

        cout << FileStartPosition << endl;
        cout << arch_size << endl;
    }
#else // if BOOST
    /*    const boost::interprocess::mode_t mode = boost::interprocess::read_only;
        fm = new boost::interprocess::file_mapping(fname.c_str(), mode);
        region = new boost::interprocess::mapped_region(*fm, mode, 0, 0);
        buf = reinterpret_cast<unsigned char*>(region->get_address()) + FileStartPosition;
        arch_size = region->get_size() - FileStartPosition;*/
    fm = new memory_mapped_file::read_only_mmf(fname.c_str());
    if (!fm->is_open())
    {
        cerr << "No file: " << fname << endl;
        exit(1);
    }
    buf = (uint8_t *)fm->data() + FileStartPosition;
    arch_size = sdsl_offset - FileStartPosition;

#endif
	uint32_t chunks_streams_size;
	memcpy(&chunks_streams_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	// cout<<"chunks_streams_size: "<<chunks_streams_size<<endl;
	for (uint32_t i = 0; i < chunks_streams_size; i++){
		size_t offset;
		uint32_t cur_chunk_actual_pos;

		memcpy(&cur_chunk_actual_pos,buf + buf_pos, sizeof(uint32_t));
		buf_pos = buf_pos + sizeof(uint32_t);

		memcpy(&offset, buf + buf_pos, sizeof(size_t));
		buf_pos = buf_pos + sizeof(size_t);

		chunks_streams[i].cur_chunk_actual_pos = cur_chunk_actual_pos;
		chunks_streams[i].offset = offset;
		
	}
    memcpy(&ploidy, buf + buf_pos, sizeof(uint8_t));
    buf_pos = buf_pos + sizeof(uint8_t);
	
    memcpy(&vec_len, buf + buf_pos, sizeof(uint64_t));
    buf_pos = buf_pos + sizeof(uint64_t);
	
    memcpy(&no_vec, buf + buf_pos, sizeof(uint64_t));
    buf_pos = buf_pos + sizeof(uint64_t);
	
    memcpy(&no_copy, buf + buf_pos, sizeof(uint64_t));
    buf_pos = buf_pos + sizeof(uint64_t);
	
    memcpy(&used_bits_cp, buf + buf_pos, sizeof(char));
    buf_pos = buf_pos + sizeof(char);

    memcpy(&bm_comp_cp_size, buf + buf_pos, sizeof(int));
    buf_pos = buf_pos + sizeof(int);

    bm_comp_copy_orgl_id.Open(buf + buf_pos, bm_comp_cp_size);
    buf_pos = buf_pos + sizeof(uint8_t) * bm_comp_cp_size;
	
    // memcpy(&max_no_vec_in_block, buf + buf_pos, sizeof(uint32_t));
    // buf_pos = buf_pos + sizeof(uint32_t);

    memcpy(&n_samples, buf + buf_pos, sizeof(uint32_t));
    buf_pos = buf_pos + sizeof(uint32_t);
	
	uint32_t chunks_min_pos_size;
	memcpy(&chunks_min_pos_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	// cout<<"chunks_min_pos_size: "<<chunks_min_pos_size<<endl;
	chunks_min_pos.resize(chunks_min_pos_size);
	memcpy(&chunks_min_pos[0], buf + buf_pos, chunks_min_pos_size * sizeof(int64_t));
	buf_pos = buf_pos + chunks_min_pos_size * sizeof(int64_t);
	// for(uint32_t i = 0; i < chunks_min_pos_size; i++){
	// 	cout<<chunks_min_pos[i]<<endl;
	// }
	uint32_t where_chrom_size;	
	memcpy(&where_chrom_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);

	d_where_chrom.resize(where_chrom_size);
	// cout<<where_chrom_size<<endl;
	for (size_t i = 0; i < where_chrom_size; ++i)
	{
    	size_t chrom_size;
		memcpy(&chrom_size, buf + buf_pos, sizeof(size_t));
		buf_pos = buf_pos + sizeof(size_t);
		// cout<<chrom_size<<endl;
    	d_where_chrom[i].first.resize(chrom_size);
    	memcpy(&d_where_chrom[i].first[0], buf + buf_pos, chrom_size*sizeof(char));
		buf_pos = buf_pos + chrom_size*sizeof(char);
		
		memcpy(&d_where_chrom[i].second, buf + buf_pos, sizeof(uint32_t));
		buf_pos = buf_pos + sizeof(uint32_t);
		
	}
	uint32_t vint_last_perm_size;
	memcpy(&vint_last_perm_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	// cout<<"vint_last_perm_size: "<<vint_last_perm_size<<endl;
	for (uint32_t i = 0; i < vint_last_perm_size; ++i)
    {
		uint32_t data_size;
		uint32_t first;
		
		memcpy(&first, buf + buf_pos, sizeof(uint32_t));
        buf_pos = buf_pos + sizeof(uint32_t);
        memcpy(&data_size, buf + buf_pos, sizeof(uint32_t));
		buf_pos = buf_pos + sizeof(uint32_t);
		// cout<<first<<":"<<data_size<<endl;
		vector<uint8_t> data(data_size);
		
        memcpy(&data[0], buf + buf_pos, data_size*sizeof(uint8_t));
		buf_pos = buf_pos + data_size*sizeof(uint8_t);

		vint_last_perm.emplace(first,data);
    }
	// cout<<"vint_last_perm_size: "<<vint_last_perm_size<<endl;
	uint32_t  comp_size;
	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);

	comp_v_header.resize(comp_size);
	memcpy(&comp_v_header[0], buf + buf_pos, comp_size*sizeof(uint8_t));
	buf_pos = buf_pos + comp_size*sizeof(uint8_t);

	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);

	comp_v_samples.resize(comp_size);
	memcpy(&comp_v_samples[0], buf + buf_pos, comp_size*sizeof(uint8_t));
	buf_pos = buf_pos + comp_size*sizeof(uint8_t);
	
	buf = buf-FileStartPosition;

	return true;
}
void DecompressionReader::SetNoThreads(uint32_t _no_threads){
	no_threads = _no_threads;
}
bool DecompressionReader::OpenReadingPart2(const string &in_file_name){

	vector<uint8_t> part2_params_data;
	vector<uint8_t> v_desc;
	size_t pos_part2_params = 0;
	uint64_t temp = 0;
	if (file_handle2)
		delete file_handle2;
	file_handle2 = new File_Handle_2(true);
    if (!file_handle2->Open(in_file_name))
	{
		cerr << "Cannot open " << in_file_name << "\n";
		return false;
	}
	int stream_id = file_handle2->GetStreamId("part2_params");
	
	if (stream_id < 0)
	{
		std::cerr << "Corrupted part2!\n";
		exit(0);
	}
	if(!file_handle2->GetPart(stream_id,part2_params_data))
	{
		std::cerr << "Corrupted part2!\n";
		exit(0);
	}
	CBSCWrapper bsc;
	bsc.InitDecompress();
	bsc.Decompress(part2_params_data, v_desc);
	uint32_t actual_varians_size = 0;
	read(v_desc, pos_part2_params, actual_varians_size);
	actual_varians.resize(actual_varians_size);
	
	for(uint32_t i = 0; i < actual_varians_size; i++){
		read(v_desc, pos_part2_params, actual_varians[i]);
		
	}
	read(v_desc, pos_part2_params, no_keys);
	read(v_desc, pos_part2_params, key_gt_id);
	keys.resize(no_keys);
	
	for(uint32_t i = 0; i < no_keys; i++){
		read(v_desc, pos_part2_params, keys[i].key_id);
		read(v_desc, pos_part2_params, keys[i].actual_field_id);
		read(v_desc, pos_part2_params, temp);
		keys[i].keys_type = (key_type_t)temp;
		read(v_desc, pos_part2_params, keys[i].type);
		// cout<<keys[i].key_id<<" "<<keys[i].actual_field_id<<endl;
	}

	InitDecompressParams();
	for (uint32_t i = 0; i < no_keys; ++i){
		decomp_part_queue->PushQueue(i);
		
	}
	part_decompress_thread.reserve(no_threads); 
    
    for (uint32_t i = 0; i < no_threads; ++i)
    {
        part_decompress_thread.emplace_back(thread([&]() {

            
			
		    while (!decomp_part_queue->IsComplete())
		    {
				uint32_t p_id;
				SPackage* pck = new SPackage;
			    if (!decomp_part_queue->PopQueue(p_id)){
					delete pck;
					break;
				}
				    
				
				pck->key_id = p_id;
				pck->stream_id_size = file_handle2->GetStreamId("key_" + to_string(p_id) + "_size");
				pck->stream_id_data = file_handle2->GetStreamId("key_" + to_string(p_id) + "_data");
				
                if(decompress_other_fileds(pck)){
					
					lock_guard<mutex> lck(m_packages);
				
					v_packages[pck->key_id] = pck;
					
				}
				else{
					
					pck->v_size.clear();
					pck->v_data.clear();

					lock_guard<mutex> lck(m_packages);
					v_packages[pck->key_id] = pck;
				}
				cv_packages.notify_all();
			
		    }      


        }));
    }

	return true;
}
//**********************************************************************************************************************
void DecompressionReader::close(){
	decomp_part_queue->Complete();

		for (uint32_t i = 0; i < no_threads; ++i)
			part_decompress_thread[i].join();
			
		// file_handle2->Close();
}
//**********************************************************************************************************************
bool DecompressionReader::decompress_other_fileds(SPackage* pck)
{
    vector<uint8_t> v_compressed;
	vector<uint8_t> v_tmp;
	CBSCWrapper *cbsc_size = v_bsc_size[pck->key_id];

	CBSCWrapper *cbsc = v_bsc_data[pck->key_id];
	
	if(!file_handle2->GetPart(pck->stream_id_size,v_compressed))
		return false;
	
	cbsc_size->Decompress(v_compressed, v_tmp);

	pck->v_size.resize(v_tmp.size()/4);

	copy_n(v_tmp.data(), v_tmp.size(), (uint8_t*)pck->v_size.data());
	
	v_compressed.clear();

	v_tmp.clear();

	file_handle2->GetPart(pck->stream_id_data,v_compressed);
	

	if(v_compressed.size()){

		if(keys[pck->key_id].type == BCF_HT_INT)
		{
			
            cbsc->Decompress(v_compressed, v_tmp);
			
			pck->v_data.resize(v_tmp.size());
			Decoder(v_tmp, pck->v_data);
			// cout<<pck->key_id<<":"<<v_compressed.size()<<":"<<v_tmp.size()<<":"<<pck->v_data.size()<<endl;
			// // cout<<pck->v_data.size()<<endl;
		}
		else{
			cbsc->Decompress(v_compressed, pck->v_data);
			// cout<<pck->key_id<<":"<<v_compressed.size()<<":"<<pck->v_data.size()<<endl;
		}
	}
	else{
		pck->v_data.clear();
		// cout<<pck->key_id<<":"<<v_compressed.size()<<":"<<pck->v_data.size()<<endl;
	}
	
	return true;
}
void DecompressionReader::Decoder(vector<uint8_t>& v_tmp, vector<uint8_t>& v_data)
{
	
    size_t size = v_tmp.size()/4;
    for (size_t i = 0; i < v_data.size()/4; i++) {
		v_data[i*4] = v_tmp[i];
		v_data[i*4+1] = v_tmp[i+size];
		v_data[i*4+2] = v_tmp[i+size*2];
		v_data[i*4+3] = v_tmp[i+size*3];

    }
	

}
//**********************************************************************************************************************
void DecompressionReader::InitDecompressParams(){

    
	v_packages.resize(no_keys,nullptr);
	decomp_part_queue = new DecompressPartQueue<uint32_t>(1);
    v_coder_part_ids.resize(no_keys, 0);
    v_bsc_data.resize(no_keys);
    v_bsc_size.resize(no_keys);
	v_i_buf.resize(no_keys);
    for (uint32_t i = 0; i < no_keys; ++i)
    {   
        v_bsc_data[i] = new CBSCWrapper();
        v_bsc_size[i] = new CBSCWrapper();
        v_bsc_size[i]->InitDecompress();
        
        switch (keys[i].type)
		{
		case BCF_HT_FLAG:
			v_bsc_data[i]->InitDecompress();
			break;
		case BCF_HT_INT:
            v_bsc_data[i]->InitDecompress();
			break;
		case BCF_HT_REAL:
			v_bsc_data[i]->InitDecompress();
			break;
		case BCF_HT_STR:
			v_bsc_data[i]->InitDecompress();
			break;
		}
    }

 
}
//**********************************************************************************************************************
void DecompressionReader::GetVariants(vector<field_desc> &fields){
	
    // Load and set fields
    for(uint32_t i = 0; i < no_keys; i++)
    {
		if (v_i_buf[i].IsEmpty())
		{
			
			unique_lock<mutex> lck(m_packages);

			cv_packages.wait(lck, [&, this] {return v_packages[i] != nullptr; });

			v_i_buf[i].SetBuffer(v_packages[i]->v_size, v_packages[i]->v_data);
			
			delete v_packages[i];

			v_packages[i] = nullptr;
			
			decomp_part_queue->PushQueue(i);
			
		}
		
		switch (keys[i].type)
		{
		case BCF_HT_INT:
			v_i_buf[i].ReadInt(fields[i].data, fields[i].data_size);
			fields[i].present = fields[i].data != nullptr;
			break;
		case BCF_HT_REAL:
			v_i_buf[i].ReadReal(fields[i].data, fields[i].data_size);
			fields[i].present = fields[i].data != nullptr;
			break;
		case BCF_HT_STR:
			v_i_buf[i].ReadText(fields[i].data, fields[i].data_size);
			fields[i].present = fields[i].data != nullptr;
			break;
		case BCF_HT_FLAG:
			uint8_t temp;
			v_i_buf[i].ReadFlag(temp);
			fields[i].present = (bool)temp;
			fields[i].data_size = 0;
			break;
		}
    }
    

}
//**********************************************************************************************************************
void DecompressionReader::decompress_meta(vector<string> &v_samples, string &header){
    vector<uint8_t> all_v_header;
    vector<uint8_t> all_v_samples;
    size_t p_header;
	size_t p_samples;
    
    for (auto data : {
		make_tuple(ref(all_v_header), ref(comp_v_header), ref(p_header), "header"),
		make_tuple(ref(all_v_samples), ref(comp_v_samples), ref(p_samples), "samples"),
		})
	{
        CBSCWrapper bsc;

		bsc.InitDecompress();
       
		bsc.Decompress(get<1>(data), get<0>(data));
		get<2>(data) = 0;
	}
	header.clear();
	read_str(all_v_header, p_header, header);
	v_samples.clear();
	string sample;
	for (uint32_t i = 0; i < n_samples; ++i)
	{
		read_str(all_v_samples, p_samples, sample);
		v_samples.emplace_back(sample);
        
	}
	cout<<"Decompress header and samples done!"<<endl;
    
}
bool DecompressionReader::setStartChunk(uint32_t start_chunk_id){
	
	auto it = chunks_streams.find(start_chunk_id + 1);
	
	if(it == chunks_streams.end())
		return false;
	buf_pos = it->second.offset;
	
	return true;
}
//get current chunk actual position
uint32_t DecompressionReader::getActualPos(uint32_t chunk_id){
	auto it = chunks_streams.find(chunk_id);
	if(it == chunks_streams.end())
		return 0;
	
	return it->second.cur_chunk_actual_pos;
}
//**********************************************************************************************************************
//read current chunk fixed fields compressed data
bool DecompressionReader::readFixedFields(){
	
	// auto it = chunks_streams.find(start_chunk_id);
	// buf_pos = it->second - buf_pos;
	memcpy(&fixed_field_block_compress.no_variants, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	uint32_t comp_size ;

	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	fixed_field_block_compress.chrom.resize(comp_size);
	memcpy(&fixed_field_block_compress.chrom[0], buf + buf_pos, comp_size*sizeof(uint8_t));
	buf_pos = buf_pos + comp_size*sizeof(uint8_t);

	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	fixed_field_block_compress.pos.resize(comp_size);
	memcpy(&fixed_field_block_compress.pos[0], buf + buf_pos, comp_size*sizeof(uint8_t));
	buf_pos = buf_pos + comp_size*sizeof(uint8_t);

	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	fixed_field_block_compress.id.resize(comp_size);
	memcpy(&fixed_field_block_compress.id[0], buf + buf_pos, comp_size*sizeof(uint8_t));
	buf_pos = buf_pos + comp_size*sizeof(uint8_t);

	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	fixed_field_block_compress.ref.resize(comp_size);
	memcpy(&fixed_field_block_compress.ref[0], buf + buf_pos, comp_size*sizeof(uint8_t));
	buf_pos = buf_pos + comp_size*sizeof(uint8_t);

	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	fixed_field_block_compress.alt.resize(comp_size);
	memcpy(&fixed_field_block_compress.alt[0], buf + buf_pos, comp_size*sizeof(uint8_t));
	buf_pos = buf_pos + comp_size*sizeof(uint8_t);

	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	fixed_field_block_compress.qual.resize(comp_size);
	memcpy(&fixed_field_block_compress.qual[0], buf + buf_pos, comp_size*sizeof(uint8_t));
	buf_pos = buf_pos + comp_size*sizeof(uint8_t);

	memcpy(&comp_size, buf + buf_pos, sizeof(uint32_t));
	buf_pos = buf_pos + sizeof(uint32_t);
	fixed_field_block_compress.gt_block.resize(comp_size);
	memcpy(&fixed_field_block_compress.gt_block[0], buf + buf_pos, comp_size*sizeof(uint8_t));
	buf_pos = buf_pos + comp_size*sizeof(uint8_t);
	
	return true;
}
void DecompressionReader::initDecoderParams(){
	p_chrom = 0;
	p_pos = 0;
	p_id = 0;
	p_ref = 0;
	p_alt = 0;
	p_gt = 0;
	p_qual = 0;
	prev_pos = 0;
	fixed_field_block_io.Clear();
	fixed_field_block_io.Initalize();
	fixed_field_block_io.Clear();
	fixed_field_block_io.Initalize();
}
//*******************************************************************************************************************************
bool DecompressionReader::Decoder(vector<block_t> &v_blocks,vector<vector<uint32_t>> &s_perm,vector<uint8_t> &gt_index,uint32_t cur_chunk_id){
	
	
	initDecoderParams();
	no_variants =  fixed_field_block_compress.no_variants;

	// cout<<"no_variants: "<<no_variants<<endl;
	vector<variant_desc_t> v_vcf_fixed_data_io;
	vector<variant_desc_t> v_vcf_sort_data_io;
	uint32_t block_size = n_samples * uint32_t(ploidy);
	v_vcf_fixed_data_io.reserve(no_variants_in_buf);
	v_vcf_sort_data_io.reserve(no_variants_in_buf);
	vector<uint32_t> perm(block_size,0);
	unique_ptr<thread> decomp_fixed_field_thread(new thread([&]{
		for (auto data : {
        	make_tuple(ref(fixed_field_block_io.chrom), ref(fixed_field_block_compress.chrom), p_chrom, "chrom"),
        	make_tuple(ref(fixed_field_block_io.id), ref(fixed_field_block_compress.id), p_id, "id"),
        	make_tuple(ref(fixed_field_block_io.alt), ref(fixed_field_block_compress.alt), p_alt, "alt"),
    		make_tuple(ref(fixed_field_block_io.qual), ref(fixed_field_block_compress.qual), p_qual, "qual"),


    	})
		{
			CBSCWrapper bsc;
			bsc.InitDecompress();
			bsc.Decompress(get<1>(data), get<0>(data));
			// cout<<get<3>(data)<<":"<<get<1>(data).size()<<":"<<get<0>(data).size()<<endl;

		}
		uint32_t i_variant;
		variant_desc_t desc;
		for(i_variant = 0; i_variant < no_variants; i_variant++) 	
		{
			// Load variant description
			read_str(fixed_field_block_io.chrom, p_chrom, desc.chrom);
			read_str(fixed_field_block_io.id, p_id, desc.id);
			read_str(fixed_field_block_io.alt, p_alt, desc.alt);
			read_str(fixed_field_block_io.qual, p_qual, desc.qual);
			v_vcf_fixed_data_io.emplace_back(desc);
		}
	
	}));
	unique_ptr<thread> decomp_sort_field__thread(new thread([&]{

		for (auto data : {
			
        	make_tuple(ref(fixed_field_block_io.pos), ref(fixed_field_block_compress.pos), p_pos, "pos"),
        	make_tuple(ref(fixed_field_block_io.ref), ref(fixed_field_block_compress.ref), p_ref, "ref"),
        	// make_tuple(ref(fixed_field_block_io.gt_block), ref(fixed_field_block_compress.gt_block), p_gt, "GT")

    	})
		{
			CBSCWrapper bsc;
			bsc.InitDecompress();
			bsc.Decompress(get<1>(data), get<0>(data));
			// cout<<get<3>(data)<<":"<<get<1>(data).size()<<":"<<get<0>(data).size()<<endl;

		}
		if(fixed_field_block_compress.gt_block.back() == 0)
		{
			fixed_field_block_compress.gt_block.pop_back();
			CBSCWrapper bsc;
			bsc.InitDecompress();
			bsc.Decompress(fixed_field_block_compress.gt_block, fixed_field_block_io.gt_block);			
			
		}else{
			fixed_field_block_compress.gt_block.pop_back();	
			zstd::zstd_decompress(fixed_field_block_compress.gt_block, fixed_field_block_io.gt_block);
			// cout<<"fixed_field_block_io:"<<fixed_field_block_io.gt_block.size()<<endl;		
		}
		uint32_t i_variant;
		variant_desc_t desc;
		for(i_variant = 0; i_variant < no_variants; i_variant++) 	
		{
			int64_t pos = 0;
			// Load variant description
			
			read(fixed_field_block_io.pos, p_pos, pos);
			
			pos += prev_pos;
			prev_pos = pos;
			desc.pos = pos;
			read_str(fixed_field_block_io.ref, p_ref, desc.ref);
			desc.filter = ".";
			desc.info = ".";
			v_vcf_sort_data_io.emplace_back(desc);
			if ((i_variant+1) % block_size == 0)
			{
				out_perm(perm, v_vcf_sort_data_io);
				s_perm.emplace_back(perm);
				v_blocks.emplace_back(v_vcf_sort_data_io);
				v_vcf_sort_data_io.clear();
			}
		}
		if (i_variant % block_size)
		{
			auto it = vint_last_perm.find(cur_chunk_id);
			// cout<<cur_chunk_id<<endl;
			// cout<<vint_last_perm.size()<<":"<<it->first<<endl;
			// if(it == vint_last_perm.end())
			// 	return false;
			// for (size_t i_p = 0; i_p < perm.size(); i_p++)
			// {
			// 	perm[i_p] = i_p;
			// }

			perm = vint_code::DecodeArray(it->second);
			s_perm.emplace_back(perm);

			for (size_t i_p = 0; i_p < i_variant % block_size; i_p++)
			{
				if (atoi(v_vcf_sort_data_io[i_p].ref.c_str()))
				{
					v_vcf_sort_data_io[i_p].ref = v_vcf_sort_data_io[i_p].ref.substr(to_string(atoi(v_vcf_sort_data_io[i_p].ref.c_str())).length());
				}
			}

			v_blocks.emplace_back(v_vcf_sort_data_io);
			v_vcf_sort_data_io.clear();
		}
	}));
	decomp_fixed_field_thread->join();
	decomp_sort_field__thread->join();
	uint32_t i_variant = 0;

	for(size_t i = 0; i < v_blocks.size(); i++) 	
	{
		for (size_t j = 0; j < v_blocks[i].data_compress.size(); j++)
		{
			v_blocks[i].data_compress[j].chrom = std::move(v_vcf_fixed_data_io[i_variant].chrom);
			v_blocks[i].data_compress[j].id = std::move(v_vcf_fixed_data_io[i_variant].id);
			v_blocks[i].data_compress[j].alt = std::move(v_vcf_fixed_data_io[i_variant].alt);
			v_blocks[i].data_compress[j].qual = std::move(v_vcf_fixed_data_io[i_variant].qual);

			i_variant++;
		}

	}
	// v_vcf_fixed_data_io.clear();
	gt_index = move(fixed_field_block_io.gt_block);
	return true;
}

// // *******************************************************************************************************************************
// void DecompressionReader::SetNoSamples(uint32_t _no_samples)
// {
// 	no_samples = _no_samples;
// }

// // *******************************************************************************************************************************
// bool DecompressionReader::GetHeader(string &_v_header)
// {
// 	_v_header = v_header;
// 	return true;
// }
// // *******************************************************************************************************************************
// bool DecompressionReader::SetHeader(string &_v_header)
// {
// 	v_header = _v_header;
// 	return true;
// }
// // *******************************************************************************************************************************
// bool DecompressionReader::SetSamples(vector<string> &_v_samples)
// {
// 	v_samples = _v_samples;
// 	return true;
// }
// // *******************************************************************************************************************************
// bool DecompressionReader::GetSamples(vector<string> &_v_samples)
// {
// 	_v_samples = v_samples;
// 	return true;
// }


