#pragma once
#include <cstring>
#include <vector>
#include "defs.h"
using namespace std;
// ************************************************************************************
typedef struct variant_desc_tag {
    string chrom;
    uint64_t pos;
    string id;
    string ref;
    string alt;
    string qual;
    string filter;
    string info;
    string format;
    int alt_num;
    variant_desc_tag():pos(0),ref(""),alt_num(0){}
    bool operator==(const struct variant_desc_tag &x)
    {
        if (chrom == "" && x.chrom == "")
            return true;

        return chrom == x.chrom &&
               pos == x.pos;
        //		&&
        //			id == x.id;
        //		&&
        //			ref == x.ref &&
        //			alt == x.alt &&
        //			qual == x.qual &&
        //			filter == x.filter;
        //			info == x.info;
    }

    bool operator !=(const struct variant_desc_tag &x)
    {
        return !operator==(x);
    }

    bool operator<(const struct variant_desc_tag &x)
    {
        if (chrom != x.chrom)
        {
            if (chrom.empty())
                return false;
            if (x.chrom.empty())
                return true;

            return chrom < x.chrom;
        }
        if (pos != x.pos)
            return pos < x.pos;
        return alt < x.alt;
    }
} variant_desc_t;
typedef struct block_tag
{
    
    std::vector<variant_desc_t> data_compress;
    block_tag(std::vector<variant_desc_t> v_vcf_data) :
            data_compress(v_vcf_data)
    {}
} block_t;
// typedef struct sortblock {
// 	int64_t val;
// 	int p;
// 	string s_ref;
// 	// string s_alt;
// 	// sortblock() : val(0), p(0), s_ref("") {}
//     sortblock(int64_t a, int b, string str_ref) : val(a), p(b), s_ref(str_ref){}
//     //  sortblock(int64_t a, int b, string str_ref,string s_alt) : val(a), p(b), s_ref(str_ref),s_alt(s_alt){}

// } sblock;
// ************************************************************************************
enum  class key_type_t {flt, info, fmt};  // FILTER / INFO / FORMAT

// ************************************************************************************
typedef struct key_desc_tag {
    uint32_t key_id;
    uint32_t actual_field_id;
    key_type_t keys_type;
    int8_t type; //one of: BCF_HT_FLAG 0 / BCF_HT_INT  1 / BCF_HT_REAL 2 / BCF_HT_STR 3
} key_desc;


// ************************************************************************************
typedef struct field_desc_tag {
    bool present = false;  // true if present in description
    char *data = nullptr;
    uint32_t data_size = 0; //current size of allocated memory

    //构造函数
    field_desc_tag() = default;
    // //拷贝构造函数
    field_desc_tag(field_desc_tag& other) {
        // cout<<"copy"<<endl;
        if (data_size) {
            delete[] data;
            data = nullptr;
        }
        present = other.present;
        data_size = other.data_size;
        
        if(data_size){
            data = new char[data_size+1];
            std::memcpy(data, other.data, data_size);
            data[data_size] = '\0';
        }

    }
    //移动构造函数
    field_desc_tag(field_desc_tag&& other) noexcept
    {
        // cout<<"Move"<<endl;
        present = other.present;
        data_size = other.data_size;
        data = other.data;
        other.present = false;
        other.data_size = 0;
        other.data = nullptr;
    }
    field_desc_tag& operator=(field_desc_tag&& other) noexcept{
        // cout<<"Move="<<endl;
        if(this != &other){
            delete [] data;
            data = nullptr;
            present = other.present;
            data_size = other.data_size;
            data = other.data;
            other.present = false;
            other.data_size = 0;
            other.data = nullptr;
        }
        return *this;
    }
    //析构函数
    ~field_desc_tag() {
        // cout<<"~field_desc_tag"<<endl;
        if(data != nullptr){
            // cout<<"delete"<<endl;
            delete[] data;
            data = nullptr;
        }
    }
    // //写一个赋值深度拷贝的函数
    // field_desc_tag& operator=(const field_desc_tag& other) {
    //     if (this != &other) {
    //         if (data_size) {
     
    //             delete[] data;
    //             data = nullptr;
    //         }
    //         present = other.present;
    //         data_size = other.data_size;
    //         if(data_size){
    //             data = new char[data_size];

    //             std::memcpy(data, other.data, data_size);
    //         }
    //         else
    //             data = nullptr;
            
    //     }
    //     return *this;
    // }
    
} field_desc;
// ************************************************************************************
// struct sort_field_block {
   

//     std::vector<uint8_t> pos;
//     std::vector<uint8_t> ref;
//     std::vector<uint8_t> gt_block;

    
//     sort_field_block(){ 
//        pos.reserve(no_variants_in_buf); 
//        ref.reserve(no_variants_in_buf); 
//        gt_block.reserve(no_variants_in_buf); 
//     }
//     sort_field_block(  std::vector<uint8_t> _pos, std::vector<uint8_t> _ref, std::vector<uint8_t> _gt_block)
//         : pos(_pos)
//         , ref(_ref)
//         , gt_block(_gt_block)
//         {}

//     void Clear(){
//         pos.clear();
//         ref.clear();
//         gt_block.clear();
//     }
//     int empty(){
//         if(pos.empty() && ref.empty() && gt_block.empty())
//             return 1;
//     }
//     void Initalize(){
//         pos.reserve(no_variants_in_buf); 
//         ref.reserve(no_variants_in_buf); 
//         gt_block.reserve(no_variants_in_buf); 
//     }

// };
// struct fixed_field_block {

//     std::vector<uint8_t> chrom;
//     std::vector<uint8_t> id;
//     std::vector<uint8_t> alt;
//     std::vector<uint8_t> qual;
//     uint32_t no_variants;

    
//     fixed_field_block(){
//        chrom.reserve(no_variants_in_buf); 
//        id.reserve(no_variants_in_buf); 
//        alt.reserve(no_variants_in_buf); 
//        qual.reserve(no_variants_in_buf); 
//        no_variants = 0;
//     }
//     fixed_field_block( std::vector<uint8_t> _chrom, std::vector<uint8_t> _id, std::vector<uint8_t> _alt, std::vector<uint8_t> _qual, uint32_t _no_variants)
//         : chrom(_chrom)
//         , id(_id)
//         , alt(_alt)
//         , qual(_qual)
//         ,no_variants(_no_variants)
//         {}
//     void Initalize(){
//         chrom.reserve(no_variants_in_buf); 
//         id.reserve(no_variants_in_buf); 
//         alt.reserve(no_variants_in_buf); 
//         qual.reserve(no_variants_in_buf);
//         no_variants = 0; 
//     }
//     void Clear(){
//         chrom.clear();
//         id.clear();
//         alt.clear();
//         qual.clear();
//         no_variants = 0;
//     }


// };
struct fixed_field_block {

    std::vector<uint8_t> chrom;
    std::vector<uint8_t> id;
    std::vector<uint8_t> alt;
    std::vector<uint8_t> qual;
    uint32_t no_variants;
    std::vector<uint8_t> pos;
    std::vector<uint8_t> ref;
    std::vector<uint8_t> gt_block;
    
    fixed_field_block(){
       chrom.reserve(no_variants_in_buf); 
       id.reserve(no_variants_in_buf); 
       alt.reserve(no_variants_in_buf); 
       qual.reserve(no_variants_in_buf); 
       no_variants = 0;
       pos.reserve(no_variants_in_buf); 
       ref.reserve(no_variants_in_buf); 
       gt_block.reserve(no_variants_in_buf); 
    }
    fixed_field_block( std::vector<uint8_t> _chrom, std::vector<uint8_t> _id, std::vector<uint8_t> _alt, std::vector<uint8_t> _qual, uint32_t _no_variants,std::vector<uint8_t> _pos, 
    std::vector<uint8_t> _ref, std::vector<uint8_t> _gt_block)
        : chrom(_chrom)
        , id(_id)
        , alt(_alt)
        , qual(_qual)
        ,no_variants(_no_variants)
        , pos(_pos)
        , ref(_ref)
        , gt_block(_gt_block)
        {}
    void Initalize(){
        chrom.reserve(no_variants_in_buf); 
        id.reserve(no_variants_in_buf); 
        alt.reserve(no_variants_in_buf); 
        qual.reserve(no_variants_in_buf);
        no_variants = 0; 
        pos.reserve(no_variants_in_buf); 
        ref.reserve(no_variants_in_buf); 
        gt_block.reserve(no_variants_in_buf); 

    }
    void Clear(){
        chrom.clear();
        id.clear();
        alt.clear();
        qual.clear();
        no_variants = 0;
        pos.clear();
        ref.clear();
        gt_block.clear();
    }


};
struct SPackage {

	int key_id;
	uint32_t stream_id_size;
	uint32_t stream_id_data;
	int part_id;
	vector<uint32_t> v_size;
	vector<uint8_t> v_data;
	int stream_id_src;
	bool is_func;
	SPackage()
	{

		key_id = -1;

		stream_id_size = 0;
		stream_id_data = 0;
		part_id = -1;
		stream_id_src = -1;
		is_func = false;
	}

	SPackage(int _key_id, uint32_t _stream_id_size, uint32_t _stream_id_data, int _part_id, vector<uint32_t>& _v_size, vector<uint8_t>& _v_data)
	{
		key_id = _key_id;
		stream_id_size = _stream_id_size;
		stream_id_data = _stream_id_data;
		stream_id_src = -1;
		part_id = _part_id;
		v_size = move(_v_size);
		v_data = move(_v_data);
		is_func = false;

		_v_size.clear();
		_v_data.clear();
	}

	SPackage( int _key_id, int _db_id, uint32_t _stream_id_size, uint32_t _stream_id_data, int _part_id, vector<uint32_t>& _v_size, int _stream_id_src)
	{
		key_id = _key_id;
		stream_id_size = _stream_id_size;
		stream_id_data = _stream_id_data;
		stream_id_src = _stream_id_src;
		part_id = _part_id;
		is_func = true;


	}
};
// struct FieldsPackage {
//     uint32_t block_id;
//     int64_t data_size;
//     uint8_t *data = nullptr;
//     uint32_t num_rows;
//     vector<variant_desc_t> v_vcf_data_compress;
//     FieldsPackage() = default;
//     //拷贝构造函数
//     FieldsPackage(uint32_t _block_id, int64_t _data_size, uint8_t* _data, uint32_t _num_rows, const vector<variant_desc_t>& _v_vcf_data_compress) 
//     : block_id(_block_id), data_size(_data_size), num_rows(_num_rows), v_vcf_data_compress(_v_vcf_data_compress)
//     {
//         if(data != nullptr){
//             cout<<"FieldsPackage destructor"<<endl;
//             delete[] data;
//             data = nullptr;
//         }
//         data = new unsigned char[data_size];
//         std::memcpy(data, _data, data_size);
//     }

//     // FieldsPackage(const FieldsPackage& other) {
//     //     block_id = other.block_id;
//     //     data_size = other.data_size;
//     //     data = new uint8_t[data_size];
//     //     std::memcpy(data, other.data, data_size);
//     //     num_rows = other.num_rows;
//     //     v_vcf_data_compress = other.v_vcf_data_compress;
//     // }
//     //析构函数
//     ~FieldsPackage() {
       

//     }
// };