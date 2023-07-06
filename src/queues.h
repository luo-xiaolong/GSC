#pragma once

#include <thread>
#include <mutex>
#include <condition_variable>
#include <numeric>
#include <queue>
#include <map>
#include <list>
#include <stack>
#include <tuple>
#include <set>
#include <functional>
#include "variant.h"
#include "defs.h"

using namespace std;

// ********************************************************************************
#include <atomic>
#include <vector>



class GtBlockQueue
{  
public:
    GtBlockQueue() : flag(false), capacity(0)
    {}
    
    GtBlockQueue(size_t _capacity) : flag(false), capacity(_capacity)
    {}
    
    ~GtBlockQueue()
    {}

    void Push(int id_block, unsigned char *data, size_t num_rows,vector<variant_desc_t> &v_vcf_data_compress)
    {
        unique_lock<std::mutex> lck(m_mutex);
        // cout<<"PUSH"<<endl;
        cv_push.wait(lck, [this] {return g_blocks.size() < capacity;});

        g_blocks.push_back(genotype_block_t(id_block, data, num_rows,std::move(v_vcf_data_compress)));

        cv_pop.notify_all();
    }
    
    bool Pop(int &id_block, unsigned char *&data, size_t &num_rows,vector<variant_desc_t> &v_vcf_data_compress)
    {
        unique_lock<std::mutex> lck(m_mutex);
        
        cv_pop.wait(lck, [this] {return !g_blocks.empty() || flag; });
        // cout<<"Pop"<<endl;
        if (flag && g_blocks.empty())
            return false;
        auto block = move(g_blocks.front());
        g_blocks.pop_front();
        id_block = block.block_id;
        data = block.data;
        num_rows   = block.num_rows;
        v_vcf_data_compress = std::move(block.v_vcf_data_compress);
        
        cv_push.notify_all();
        
        return true;
    }
    
    void Complete()
    {
        unique_lock<std::mutex> lck(m_mutex);
        
        flag = true;
        
        cv_pop.notify_all();
    }
private:
    typedef struct genotype_block
    {
        int block_id;
        unsigned char *data;
        size_t num_rows;
        vector<variant_desc_t> v_vcf_data_compress;
        genotype_block(int _block_id, unsigned char* _data, size_t _num_rows,vector<variant_desc_t> &&_v_vcf_data_compress):
        block_id(_block_id), data(_data), num_rows(_num_rows), v_vcf_data_compress(move(_v_vcf_data_compress))
        {}
    } genotype_block_t;
    
    list<genotype_block_t> g_blocks;
    
    bool flag;
    size_t capacity;
    
    mutex m_mutex;
    condition_variable cv_pop, cv_push;
};
// class GtBlockQueue
// {  
// public:
//     GtBlockQueue() : flag(false), capacity(0)
//     {}
    
//     GtBlockQueue(size_t _capacity) : flag(false), capacity(_capacity)
//     {}
    
//     ~GtBlockQueue()
//     {}

//     void Push(int id_block, unsigned char *data, size_t num_rows,vector<variant_desc_t> &v_vcf_data_compress,chrom_flag chrom_end_flag)
//     {
//         unique_lock<std::mutex> lck(m_mutex);
//         cv_push.wait(lck, [this] {return g_blocks.size() < capacity;});
        
//         g_blocks.emplace(genotype_block_t(id_block, data, num_rows,v_vcf_data_compress,chrom_end_flag));
        
//         cv_pop.notify_all();
//     }
    
//     bool Pop(int &id_block, unsigned char *&data, size_t &num_rows,vector<variant_desc_t> &v_vcf_data_compress,chrom_flag &chrom_end_flag)
//     {
//         unique_lock<std::mutex> lck(m_mutex);
//         cv_pop.wait(lck, [this] {return !g_blocks.empty() || flag; });
        
//         if (flag && g_blocks.empty())
//             return false;
//         auto block = move(g_blocks.front());
//         g_blocks.pop();
//         id_block = block.block_id;
//         data = block.data;
//         num_rows   = block.num_rows;
//         v_vcf_data_compress = std::move(block.v_vcf_data_compress);
//         chrom_end_flag = block.chrom_end_flag;
        
//         cv_push.notify_all();
        
//         return true;
//     }
    
//     void Complete()
//     {
//         unique_lock<std::mutex> lck(m_mutex);
        
//         flag = true;
        
//         cv_pop.notify_all();
//     }
// private:
//     typedef struct genotype_block
//     {
//         int block_id;
//         unsigned char *data;
//         size_t num_rows;
//         vector<variant_desc_t> v_vcf_data_compress;
//         chrom_flag chrom_end_flag;
//         genotype_block(int _block_id, unsigned char* _data, size_t _num_rows,vector<variant_desc_t> &_v_vcf_data_compress,chrom_flag chrom_end_flag):
//         block_id(_block_id), data(_data), num_rows(_num_rows), v_vcf_data_compress(_v_vcf_data_compress),chrom_end_flag(chrom_end_flag)
//         {}
//     } genotype_block_t;
    
//     queue<genotype_block_t> g_blocks;
    
//     bool flag;
//     size_t capacity;
    
//     mutex m_mutex;
//     condition_variable cv_pop, cv_push;
// };
// ********************************************************************************
// class ProcessingBlockQueue
// {
// private:
//     struct processing_block
//     {
//         int block_id;
//         vector<bool> zeros;
//         vector<bool> copies;
//         const uint32_t* origin_of_copy;
//         uint32_t copy_num;

//         processing_block(int _block_id, vector<bool>&& _zeros, vector<bool>&&_copies, const uint32_t* _origin_of_copy, uint32_t _copy_num) :
//             block_id(_block_id), zeros(std::move(_zeros)), copies(std::move(_copies)), origin_of_copy(_origin_of_copy), copy_num(_copy_num)
//         {}

//         bool operator<(const processing_block& other) const
//         {
//             return block_id < other.block_id;
//         }
//     };

//     set<processing_block> p_blocks;

//     mutex m_mutex;

// public:
//     ProcessingBlockQueue() {}

//     ~ProcessingBlockQueue() {}

//     void Push(int id_block, vector<bool>&& zeros, vector<bool>&& copies, const uint32_t* origin_of_copy, uint32_t copy_num)
//     {
//         lock_guard<std::mutex> lck(m_mutex);
//         p_blocks.emplace(id_block, std::move(zeros), std::move(copies), origin_of_copy, copy_num);
//     }

//     bool Pop(int& id_block, vector<bool>& zeros, vector<bool>& copies, uint32_t*& origin_of_copy, uint32_t& copy_num)
//     {
//         unique_lock<std::mutex> lck(m_mutex);

//         if (p_blocks.empty())
//             return false;

//         auto p = p_blocks.begin();

//         id_block = p->block_id;
//         zeros = std::move(p->zeros);
//         copies = std::move(p->copies);
//         origin_of_copy = const_cast<uint32_t*>(p->origin_of_copy);
//         copy_num = p->copy_num;
//         p_blocks.erase(p_blocks.begin());

//         return true;
//     }

//     void get_all_copy_num(uint64_t& all_copy_num) const
//     {
//         all_copy_num = accumulate(p_blocks.begin(), p_blocks.end(), 0ULL, [](uint64_t sum, const processing_block& pb) { return sum + pb.copy_num; });
//     }
// };

// template<typename BlockIDType, typename DataType>
// class BlockVariants
// {  
// public:
//     BlockVariants() : flag(false), capacity(0)
//     {}
    
//     BlockVariants(size_t _capacity) : flag(false), capacity(_capacity)
//     {}
    
//     ~BlockVariants()
//     {}
    
//     void Push(BlockIDType _block_id,DataType &v_block)
//     {
//         unique_lock<std::mutex> lck(m_mutex);
//         cv_push.wait(lck, [this] {return (int) comp_blocks.size() < capacity;});
        
//         comp_blocks.emplace(block_t(_block_id,v_block));
        
//         cv_pop.notify_all();
//     }

//     bool Pop(BlockIDType &block_id,DataType &v_block)
//     {
//         unique_lock<std::mutex> lck(m_mutex);
//         cv_pop.wait(lck, [this] {return !comp_blocks.empty() || flag; });
        
//         if (flag && comp_blocks.empty())
//             return false;

//         block_id = comp_blocks.front().block_id;
//         v_block = comp_blocks.front().v_block;


//         comp_blocks.pop();
        
//         cv_push.notify_all();
        
//         return true;
//     }
    
//     void Complete()
//     {
//         unique_lock<std::mutex> lck(m_mutex);
        
//         flag = true;
        
//         cv_pop.notify_all();
//     }
// private:
//     typedef struct block
//     {   
//         BlockIDType block_id;
//         DataType v_block;
        
//         block(BlockIDType _block_id,DataType &v_block):
//             block_id(_block_id),v_block(v_block){}
//     } block_t;
//     queue<block_t> comp_blocks;

//     bool flag;
//     int capacity;
    
//     mutex m_mutex;
//     condition_variable cv_pop, cv_push;
// };

// class DecompBlockVariants
// {
// private:
//      typedef struct decomp_block_variants
//     {
//         int block_id;
//         vector<block_t> q_blocks;
//         vector<vector<int>> s_perm;
//         vector<uint8_t> samples_indexes;
//         decomp_block_variants(int _block_id,vector<block_t> &_q_blocks,vector<vector<int>> &_s_perm,vector<uint8_t> &_samples_indexes):block_id(_block_id),q_blocks(_q_blocks),s_perm(_s_perm),samples_indexes(_samples_indexes){}
//         friend bool operator<(const decomp_block_variants &x, const decomp_block_variants &y)
//         {
//             return x.block_id < y.block_id;
//         }
//     } decomp_block_variants_t;
    
//     set<decomp_block_variants_t> decomp_blocks;

//     mutex m_mutex;

// public:
//     DecompBlockVariants()
//     {}
    
//     ~DecompBlockVariants()
//     {}
    
//     void Push(int _block_id,vector<block_t> &_q_blocks,vector<vector<int>> &_s_perm,vector<uint8_t> &_samples_indexes)
//     {
//         lock_guard<std::mutex> lck(m_mutex);

//         decomp_blocks.insert(decomp_block_variants_t(_block_id,_q_blocks,_s_perm,_samples_indexes));
        
//     }
    
//     bool Pop(int &block_id,vector<block_t> &q_blocks,vector<vector<int>> &s_perm,vector<uint8_t> &samples_indexes)
//     {
//         unique_lock<std::mutex> lck(m_mutex);
        
//         if (decomp_blocks.empty())
//             return false;
        
//         auto de_block = decomp_blocks.begin();
//         block_id = de_block->block_id;
//         q_blocks = de_block->q_blocks;
//         s_perm = de_block->s_perm;
//         samples_indexes=de_block->samples_indexes;
//         decomp_blocks.erase(decomp_blocks.begin());
        
//         return true;
//     }
   
// };
// class CBlockVariants
// {  
// public:
//     CBlockVariants() 
//     {}
    
//     ~CBlockVariants()
//     {}
    
//     void Push(int _block_id,VBlock &_v_block,uint32_t no_variants)
//     {
//         lock_guard<std::mutex> lck(m_mutex);

        
//         comp_blocks.insert(cblock_variants_t(_block_id,_v_block,no_variants));

//     }
    
//     bool Pop(int &block_id,VBlock &v_block,uint32_t &no_variants)
//     {

//         unique_lock<std::mutex> lck(m_mutex);
//         if (comp_blocks.empty())
//             return false;
//         auto c_block = comp_blocks.begin();
//         block_id = c_block->block_id;
//         v_block = c_block->v_block;
//         no_variants = c_block->no_variants;
//         comp_blocks.erase(comp_blocks.begin());
        
//         return true;
//     }

// private:
//     typedef struct cblock_variants
//     {   
//         int block_id;
//         VBlock v_block;
//         uint32_t no_variants;
        
//         cblock_variants(int _block_id,VBlock &_v_block,uint32_t _no_variants):
//             block_id(_block_id),v_block(_v_block),no_variants(_no_variants){}
//         friend bool operator<(const cblock_variants &x, const cblock_variants &y)
//         {
//             return x.block_id < y.block_id;
//         }
//     } cblock_variants_t;
    
//     set<cblock_variants_t> comp_blocks;
    
//     mutex m_mutex;

// };


// template<typename BlockIDType, typename SizeType, typename DataType>
// class BlockingQueue
// {
// public:
//     BlockingQueue() : flag(false), capacity(0) {}

//     BlockingQueue(size_t _capacity) : flag(false), capacity(_capacity) {}

//     ~BlockingQueue() {}

//     void Push(BlockIDType block_id, vector<SizeType>&& v_size, vector<DataType>&& v_data)
//     {
//         unique_lock<std::mutex> lck(m_mutex);
//         cv_push.wait(lck, [this] { return blocks.size() < capacity; });

//         blocks.emplace(Blocks_t(block_id, move(v_size), move(v_data)));

//         cv_pop.notify_all();
//     }

//     bool Pop(BlockIDType& block_id, vector<SizeType>& v_size, vector<DataType>& v_data)
//     {
//         unique_lock<std::mutex> lck(m_mutex);
//         cv_pop.wait(lck, [this] { return !blocks.empty() || flag; });

//         if (flag && blocks.empty())
//             return false;

//         auto temp = move(blocks.front());
//         blocks.pop();

//         block_id = temp.block_id;
//         v_size = move(temp.v_size);
//         v_data = move(temp.v_data);

//         cv_push.notify_all();

//         return true;
//     }

//     void Complete()
//     {
//         unique_lock<std::mutex> lck(m_mutex);

//         flag = true;

//         cv_pop.notify_all();
//     }

// private:
//     typedef struct Blocks
//     {
//         BlockIDType block_id;
//         vector<SizeType> v_size;
//         vector<DataType> v_data;
//         Blocks(){}
//         Blocks(BlockIDType _block_id,vector<SizeType> &&_v_size,vector<DataType> &&_v_data):block_id(_block_id),v_size(std::move(_v_size)),v_data(std::move(_v_data)){}
     
//     } Blocks_t;
//     queue<Blocks_t> blocks;

//     bool flag;
//     size_t capacity;

//     mutex m_mutex;
//     condition_variable cv_pop, cv_push;
// };

// template<typename BlockIDType, typename SizeType, typename DataType>
// class CompOtherFields
// {
// private:
//     typedef struct comp_other_fields
//     {
//         BlockIDType block_id;
//         vector<SizeType> comp_v_size;
//         vector<DataType> comp_v_data;
//         comp_other_fields(BlockIDType _block_id,vector<SizeType> &_comp_v_size,vector<DataType> &_comp_v_data):block_id(_block_id),comp_v_size(_comp_v_size),comp_v_data(_comp_v_data){}
//         friend bool operator<(const comp_other_fields &x, const comp_other_fields &y)
//         {
//             return x.block_id < y.block_id;
//         }
//     } comp_other_fields_t;
    
//     set<comp_other_fields_t> other_fields;

//     mutex m_mutex;

// public:
//     CompOtherFields()
//     {}
    
//     ~CompOtherFields()
//     {}
    
//     void Push(BlockIDType _block_id,vector<SizeType> &_comp_v_size,vector<DataType> &_comp_v_data)
//     {
//         lock_guard<std::mutex> lck(m_mutex);

//         other_fields.insert(comp_other_fields_t(_block_id,_comp_v_size,_comp_v_data));
        
//     }
    
//     bool Pop(BlockIDType &block_id,vector<SizeType> &comp_v_size,vector<DataType> &comp_v_data)
//     {
        
//         unique_lock<std::mutex> lck(m_mutex);
        
//         if (other_fields.empty())
//             return false;
        
//         auto o_fields = other_fields.begin();
//         block_id = o_fields->block_id;
//         comp_v_size = o_fields->comp_v_size;
//         comp_v_data = o_fields->comp_v_data;
       
//         other_fields.erase(other_fields.begin());
        
//         return true;
//     }
   
// };
template<typename DataType>
class VarBlockQueue
{  
public:
    VarBlockQueue() : flag(false), capacity(0)
    {}
    
    VarBlockQueue(size_t _capacity) : flag(false), capacity(_capacity)
    {}
    
    ~VarBlockQueue()
    {}
    
    void Push(uint32_t id_block, DataType &data)
    {
        unique_lock<std::mutex> lck(m_mutex);
        cv_push.wait(lck, [this] {return var_blocks.size() < capacity;});

        var_blocks.emplace(variant_block_t(id_block,data));
        
        cv_pop.notify_all();
    }
    
    bool Pop(uint32_t &id_block, DataType &data)
    {
        unique_lock<std::mutex> lck(m_mutex);
        cv_pop.wait(lck, [this] {return !var_blocks.empty() || flag; });
        if (flag && var_blocks.empty())
            return false;
        auto block = std::move(var_blocks.front());
        var_blocks.pop();
        id_block = block.block_id;

        data = std::move(block.data);

        
        cv_push.notify_all();
        
        return true;
    }
    
    void Complete()
    {
        unique_lock<std::mutex> lck(m_mutex);
        
        flag = true;
        
        cv_pop.notify_all();
    }
private:
    typedef struct variant_block
    {
        uint32_t block_id;
        DataType data;

        variant_block(uint32_t _block_id,DataType &_data):
        block_id(_block_id), data(_data)
        {}
    } variant_block_t;
    
    queue<variant_block_t> var_blocks;
    
    bool flag;
    size_t capacity;
    
    mutex m_mutex;
    condition_variable cv_pop, cv_push;
};

template<typename DataType>

class CompVarBlockQueue
{
private:
    typedef struct comp_fixed_fields
    {
        uint32_t block_id;
        DataType comp_v_data;
        comp_fixed_fields(uint32_t _block_id,DataType &_comp_v_data):block_id(_block_id),comp_v_data(_comp_v_data){}
        friend bool operator<(const comp_fixed_fields &x, const comp_fixed_fields &y)
        {
            return x.block_id < y.block_id;
        }
    } comp_fixed_fields_t;
    
    set<comp_fixed_fields_t> fixed_fields;

    mutex m_mutex;

public:
    CompVarBlockQueue()
    {}
    
    ~CompVarBlockQueue()
    {}
    
    void Push(uint32_t _block_id,DataType &_comp_v_data)
    {
        lock_guard<std::mutex> lck(m_mutex);

        fixed_fields.insert(comp_fixed_fields_t(_block_id,_comp_v_data));
        
    }
    
    bool Pop(uint32_t &block_id,DataType &comp_v_data)
    {
        
        unique_lock<std::mutex> lck(m_mutex);
        
        if (fixed_fields.empty())
            return false;
        
        auto f_fields = fixed_fields.begin();
        block_id = f_fields->block_id;
        comp_v_data = move(f_fields->comp_v_data);
       
        fixed_fields.erase(fixed_fields.begin());
        
        return true;
    }
   
};

template<typename PartType>
class PartQueue
{  
public:
    PartQueue() : flag(false), capacity(0)
    {}
    
    PartQueue(size_t _capacity) : flag(false), capacity(_capacity)
    {}
    
    ~PartQueue()
    {}
    
    void Push(PartType &data)
    {
        unique_lock<std::mutex> lck(m_mutex);
        cv_push.wait(lck, [this] {return part_queue.size() < capacity;});

        part_queue.emplace_back(move(data));;
        
        cv_pop.notify_all();
    }

    // void PushQueue(PartType &data)
    // {
    //     unique_lock<std::mutex> lck(m_mutex);
	// 	bool was_empty = n_elements == 0;
	// 	q.push_back(data);
	// 	++n_elements;

	// 	if(was_empty)
	// 		cv_queue_empty.notify_all();
    // }

    // bool PopQueue(PartType &data)
	// {
	// 	unique_lock<std::mutex> lck(m_mutex);
	// 	cv_pop.wait(lck, [this]{return !q.empty() || flag;}); 

    //     if (flag && q.empty())
    //         return false;
        
	// 	data = std::move(q.front());
	// 	q.pop();
    //     if(q.empty())
	// 	    cv_push.notify_all();
	// 	return true;
	// }

    template<typename S>
    bool Pop(PartType &data,const std::function<bool(S &item)> fo)
    {
        unique_lock<std::mutex> lck(m_mutex);
        cv_pop.wait(lck, [this] {return !part_queue.empty() || flag; });
        if (flag && part_queue.empty())
            return false;
        auto cur_part = part_queue.begin();
        for (; cur_part != part_queue.end(); ++cur_part)
			if (fo(*cur_part))
				break;
        if (cur_part == part_queue.end())
			cur_part = part_queue.begin();
        data = move(*cur_part);
		part_queue.erase(cur_part);
        // auto cur_part = std::move(part_queue.front());
        // part_queue.pop();

        // data = std::move(cur_part.data);

        cv_push.notify_all();
        
        return true;
    }
    
    void Complete()
    {
        unique_lock<std::mutex> lck(m_mutex);
        
        flag = true;
        
        cv_pop.notify_all();
    }

private:

    list<PartType> part_queue;
    queue<PartType> q;
    bool flag;
    size_t capacity;
    
    mutex m_mutex;
    condition_variable cv_pop, cv_push;
};
template<typename PartType>
class DecompressPartQueue
{  
public:
    DecompressPartQueue() : flag(false), n_producers(0), n_elements(0)
    {}
    
    DecompressPartQueue(int _n_producers) :flag(false),n_producers(_n_producers),n_elements(0)
    {}
    
    ~DecompressPartQueue()
    {}
    
    void PushQueue(PartType data)
    {
        unique_lock<std::mutex> lck(m_mutex);
        
		bool was_empty = n_elements == 0;
		
        part_queue.push_back(data);

		++n_elements;

		if(was_empty)
			cv_queue_empty.notify_all();
    }

    bool PopQueue(PartType &data)
	{
        unique_lock<std::mutex> lck(m_mutex);
		cv_queue_empty.wait(lck, [this]{return !this->part_queue.empty() || !this->n_producers;}); 

		if(n_elements == 0)
			return false;

		data = move(part_queue.front());
//		q.pop();
		part_queue.pop_front();
		--n_elements;
		if(n_elements == 0)
			cv_queue_empty.notify_all();

		return true;
	}

    bool IsComplete()
	{
		lock_guard<mutex> lck(m_mutex);

		return n_elements == 0 && n_producers == 0;
	}
    void Complete()
    {
        unique_lock<std::mutex> lck(m_mutex);
        
		n_producers--;

		if(!n_producers)
			cv_queue_empty.notify_all();
    }

private:

    list<PartType> part_queue;
    bool flag;
    int n_producers;
	uint32_t n_elements;
    mutex m_mutex;
    condition_variable cv_queue_empty;
};
