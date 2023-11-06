#pragma once

//#define DEVELOPMENT_MODE
#define MAX_NUMBER_OF_GROUP 32

#define TWO_ZEROS_VECTOR

#define FORCED_BV_TYPE RRR

#define OUR_STRTOL

#ifdef WIN32

#define my_fopen fopen

#define my_fseek _fseeki64

#define my_ftell _ftelli64


#else
#define my_fopen fopen
#define my_fseek fseek
#define my_ftell ftell
#define _TCHAR char

#define _tmain main

#define MAX 0x7FFFFFFFFFFFFFFF
#define GT_NOT_CALL -0x7FFFFFFF
#define bcf_next_gt_missing 1

#define CHUNK_SIZE1 1<<16
#define CHUNK_SIZE2 1<<15
#define CHUNK_SIZE3 1<<14
// typedef unsigned short int uint16;

// typedef int int32;

// typedef short int int16;

// typedef unsigned int uint32;

// typedef long long int64;

// typedef unsigned long long uint64;
// typedef unsigned char uchar;
// typedef unsigned char uchar_t;
// typedef unsigned int uint32_t;
#include <stdio.h>

#include <stdlib.h>

#include <math.h>

#include <limits>

#include <string>

#include <cstddef>
#endif 
// const uint32 MC_TRIALS = 100000;

// const uint32 MAX_HUF_LUT_LEN = 16;
enum class chrom_flag {none, last_chrom_flag, end_chrom_flag};

// const size_t BLOCK_SIZE = 1<<16;

const uint32_t NO_SAMPLE_THRESHOLD = 4000;

const size_t COMP_BLOCK_SIZE = 1 << 30;

const int max_part_size = 4;

const size_t no_variants_in_buf = 8192u;

const uint32_t p_bsc_features = 1u;

const uint32_t max_buffer_size = 8 << 20;

const float min_float = std::numeric_limits<float>::lowest();

const float max_float = std::numeric_limits<float>::max();
// variable "where" and "whichByte_whereInRes" changed to uint32_t from uchar

// const uint32_t MATCH_BITS_HUF = 8;
// 8-10

// const uint32_t END_FLAG = 0xFFFFFFFFu;

// const uint32_t REMOVED = 0xFFFFFFFEu;

// const uint32_t ZEROS = 0xFFFFFFFDu;

// const uint32_t NO_POS = ~((uint32_t)0);

// const uint32_t HT_INIT_SIZE = 1 << 10; // Power of 2, cannot be smaller than PART_SIZE*2

// const uint32_t HASH_KEY_LEN1 = 32;

// const uint32_t HASH_KEY_LEN2 = 5;

// const uint32_t MIN_ZERO_RUN_LEN = 2;

// const uint32_t MIN_ONES_RUN_LEN = 2;

// const uint32_t FULL_POS_STEP = 1025; // So there are 1024 * bits_used between full positions

// const uint32_t BITS_IN_BYTE = 8;

// const uint32_t MAX_LITERAL_RUN = 252; // max 252

// const uint32_t MIN_LITERAL_RUN = 20; // min 2

// const uint32_t PART_SIZE = 12; // Number of variants processed in block (number of vectors processed in block: 2*PART_SIZE)

// const uint32_t PART_TRIALS = 10 * PART_SIZE * 2 / 20; 

// const uint32_t MC_ARRAY_SIZE = (PART_TRIALS + 63) / 64;

