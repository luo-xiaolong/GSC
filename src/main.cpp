
#include <iostream>

#include <string>

#include "gsc_params.h"

#include "bit_memory.h"

#include "defs.h"

#include "queues.h"

#include "compression_reader.h"

#include "block_processing.h"

#include "compressor.h"

#include "decompressor.h"

// #include <algorithm>

#include <cstdio>

#include <vector>

// #include <list>

// #include <unordered_map>

// #include <nmmintrin.h>

#include <chrono>

#include <time.h>


using namespace std;
using namespace std::chrono;

int usage();

int usage_compress();

int usage_decompress();

int compress_entry();

int decompress_entry();

int params_options(int argc, const char *argv[]);

void decom(Decompressor &decompressor);
//--------------------------------------------------------------------------------
GSC_Params params;
// Show execution options

int usage()

{

    cerr << "Usage: gsc [option] [arguments] " << endl;

    cerr << "Available options: " << endl;

    cerr << "\tcompress - compress VCF/BCF file" << endl;

    cerr << "\tdecompress     - query and decompress to VCF/BCF file " << endl;
    
    exit(0);
}

int usage_compress()

{

    cerr << "Compress usage: " << endl;
    
    cerr << "\tgsc compress <options> [out_file] [in_file]   \n" << endl;

    cerr << "Input: " << endl;

    cerr << "\t[in_file]\t path to input file (a VCF or VCF.GZ file by default)"<< endl;
    
    cerr << "\t--bcf    \t- input is a BCF file (input is a VCF or VCF.GZ file by default)\t" << endl;

    cerr << "\t--polidy [x]\t- set ploidy of samples in input VCF to [x] (number >= 1; 2 by default)" << endl;

    cerr << "Output: " << endl;

    cerr << "\t--out [name]\t- output to a file and set output name to [name] \t" << endl;
         
    cerr << "\t[out_file]\t- path to output file \n"<< endl;

    cerr << "Parameters options : " << endl;

    cerr << "\t--all\t- set to fully compressed mode(compress the entire VCF file,by default only genotype and VCF fixed fields are compressed)" << endl;

    cerr << "\t--threads [x]\t- set number of threads to [x] (number >= 1; 2 by default)" << endl;

    cerr << "\t--depth [x]\t- set the maximum replication depth to [x] (number >= 0; 0 means no matches; 100 by default)" << endl;

    exit(0);
}

int usage_decompress()

{

    cerr << "Decompress and Query usage:" << endl;

    cerr << "\tgsc decompress <options> [out_file] [in_file]" << endl;

    cerr << "Input: " << endl;

    cerr << "\t[in_file]\t path to input file (prefix of the file name to be decompressed)"<< endl;

    cerr << "Output: " << endl;

    cerr << "\t--out [name]\t- output to a file and set output name to [name] \t" << endl;

    // cerr << "\t-S [name]\t- output a file with all the samples and set output samples file name to [name] \t" << endl;
    cerr << "\t--all\t- set to fully decompressed mode(decompress the entire VCF file,by default only genotype and VCF fixed fields are decompressed)" << endl;

    cerr << "\t--bcf\t- output a BCF file and please use it together with param '-o' (output is a VCF file by default)\t" << endl;

    cerr<<"\t[out_file]\t- you need to enter the output file path "<< endl;

    cerr << "\t-G \t- don't output sample genotypes (only #CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO columns)" << endl;

    cout << "\t-C \t- write AC/AN to the INFO field (always set when using -minAC, -maxAC, -minAF or -maxAF)"<< endl;

    cerr << "Query: " << endl;
    
    cerr << "\tWhen querying, you must switch to the default mode (that is, remove the '--all' in the parameter), and only need to use the parameter '--all' when you need to decompress the entire VCF file" << endl;
    
    cerr << "\t--range\t- range in format [start],[end] (for example: -r 4999756,4999852). By default all variants are decompressed." << endl;

    cerr << "\t-samples-name\t- sample name(s), separated by comms (for example: -s HG03861,NA18639) OR '@' sign followed by the name of a file with sample name(s) separated by whitespaces (for exaple: -s @file_with_IDs.txt). By default all samples/individuals are decompressed" << endl;

    // cerr << "\t-n X \t- process at most X records (by default: all from the selected range)" << endl;

    cerr << "Settings: " << endl;

    cerr << "\t--minAC X \t- report only sites with count of alternate alleles among selected samples smaller than or equal to X (default: no limit)" << endl;

    cerr << "\t--maxAC X \t- report only sites with count of alternate alleles among selected samples greater than or equal to X" << endl;

    cerr << "\t--minAF X \t- report only sites with allele frequency among selected samples greather than or equal to X (X - number between 0 and 1; default: 0)" << endl;

    cerr << "\t--maxAF X \t- report only sites with allele frequency among selected samples smaller than or equal to X (X - number between 0 and 1; default: 1)" << endl;

    // cerr << "\t-m X\t- limit maximum memory usage to remember previous vectors to X MB (no limit by default)\t" << endl;

    cerr << endl;

    exit(0);
}
//Main program entry
int main(int argc, const char *argv[])
{
    
    high_resolution_clock::time_point start = high_resolution_clock::now();

    int result = 0;

    if (!params_options(argc, argv))
		return 1;

    if(params.task_mode == task_mode_t::mcompress){
        result = compress_entry();
        if (result)
		    std::cout << "Compression error!!!\n";
    }

    else if(params.task_mode == task_mode_t::mdecompress){
        result = decompress_entry();
        if (result)
		    std::cout << " Query error!!!\n";
    }

    high_resolution_clock::time_point end = high_resolution_clock::now();

	duration<double> time_duration = duration_cast<duration<double>>(end - start);

	std::cout << "Total processing time: " << time_duration.count() << " seconds.\n";

    return result;
}

// Parse the parameters
int params_options(int argc, const char *argv[]){

    if (argc < 2)
	{
		return usage();
	}

    if (string(argv[1]) == "compress")
		params.task_mode = task_mode_t::mcompress;

	else if (string(argv[1]) == "decompress")
		params.task_mode = task_mode_t::mdecompress;
    
    else
        return usage();

    if(params.task_mode == task_mode_t::mcompress){
         
        if (argc < 4)
            return usage_compress();
        
        int i;
        int temp;
        
        for(i = 2; i < argc - 1; ++i){

            if (argv[i][0] != '-'){
                return usage_decompress();
                break;
            }
            if (strncmp(argv[i], "--out", 5) == 0){

                params.out_file_flag = true;

                i++;

                if (i >= argc)
                    return usage_decompress();

                params.out_file_name = string(argv[i]);
            }
            else if (strncmp(argv[i], "--bcf", 5) == 0)
                params.in_type = file_type::BCF_File;

            else if (strncmp(argv[i], "--all", 5) == 0)
                params.compress_all = true;
            
            else if (strncmp(argv[i], "--ploidy", 8) == 0){

                i++;

                if (i >= argc)
                    return usage_compress();

                temp = atoi(argv[i]);

                if (temp < 1)
                    usage_compress();

                params.ploidy = temp;
            }
            else if (strncmp(argv[i], "--depth", 7) == 0){
                
                i++;

                if (i >= argc)
                    return usage_compress();

                temp = atoi(argv[i]);

                if (temp < 0)
                    usage_compress();

                params.max_replication_depth = temp;
            }

            else if (strncmp(argv[i], "--threads", 2) == 0){

                i++;

                if (i >= argc)
                    return usage_compress();

                temp = atoi(argv[i]);

                if (temp < 1)
                    usage_compress();

                params.no_threads = temp;
            }
        }
        if (i >= argc)
            return usage_compress();
        if(params.out_file_flag == false)
            params.out_file_name = string(argv[argc - 1]);
        params.in_file_name = string(argv[argc - 1]);

    }

    else if(params.task_mode == task_mode_t::mdecompress){
        
        if (argc < 3)
            return usage_decompress();

        int i,temp;
        float temp_f;

        for (i = 2; i < argc - 1; ++i){
            
            if (argv[i][0] != '-'){
                usage_decompress();
                break;
            }
            if (strncmp(argv[i], "--out", 5) == 0){

                params.out_file_flag = true;

                i++;

                if (i >= argc)
                    return usage_decompress();

                params.out_file_name = string(argv[i]);
            }
            else if (strncmp(argv[i], "--all", 5) == 0)
                params.compress_all = true;
            
            else if (strncmp(argv[i], "--bcf", 5) == 0)

                params.out_type = file_type::BCF_File;
            

            else if (strncmp(argv[i], "--samples-name", 14) == 0){

                i++;

                params.samples = string(argv[i]);
            }

            else if (strncmp(argv[i], "--range", 7) == 0){

                i++;

                params.range = string(argv[i]);
            }

            else if (strncmp(argv[i], "-n", 2) == 0){

                i++;

                if (i >= argc)
                    return usage_decompress();

                temp = atoi(argv[i]);

                if (temp < 0)
                    usage_decompress();
                
                else
                    params.records_to_process = temp;
            
            }
            else if (strncmp(argv[i], "-S", 2) == 0){
                
                params.out_samples_name = true;
                
                i++;

                if (i >= argc)
                    return usage_decompress();

                params.out_samples_file_name = string(argv[i]);

            }

            else if (strncmp(argv[i], "-C", 2) == 0)

                params.out_AC_AN = true;

            else if (strncmp(argv[i], "-G", 2) == 0)

                params.out_genotypes = false;

            else if(strncmp(argv[i], "-V", 2) == 0)
                params.out_ohter_fields = true;

            else if (strncmp(argv[i], "--minAC", 7) == 0){

                i++;

                if (i >= argc)
                    return usage_decompress();

                temp = atoi(argv[i]);

                if (temp < 0)
                    usage_decompress();

                params.minAC = temp;

                params.out_AC_AN = true;
            }

            else if (strncmp(argv[i], "--maxAC", 7) == 0){

                i++;

                if (i >= argc)
                    return usage_decompress();

                temp = atoi(argv[i]);

                if (temp < 0)
                    usage_decompress();

                params.maxAC = temp;

                params.out_AC_AN = true;
            }

            else if (strncmp(argv[i], "--minAF", 7) == 0){

                i++;

                if (i >= argc)
                    return usage_decompress();

                temp_f = atof(argv[i]);

                if (temp_f < 0.0 || temp_f > 1.0)
                    usage_decompress();

                params.minAF = temp_f;

                params.out_AC_AN = true;
            }

            else if (strncmp(argv[i], "--maxAF", 7) == 0){

                i++;

                if (i >= argc)
                return usage_decompress();

                temp_f = atof(argv[i]);

                if (temp_f < 0.0 || temp_f > 1.0)
                    usage_decompress();
                
                params.maxAF = temp_f;

                params.out_AC_AN = true;
            }
      
        }

        if (i >= argc)
            return usage_decompress();

        params.in_file_name = string(argv[i]);
    
    } 
    return 1;  
}
 //**********************************************************************************************************************************

//  Program compression inlet
int compress_entry()

{

    Compressor compressor(params);  //Passing compression parameters.
    if(!compressor.CompressProcess())
        return 1;


    return 0;
}    
// *********************************************************************************************************************
//  Program decompression inlet
int decompress_entry(){

    // bool result = true;    
    if(params.out_type == file_type::BCF_File && params.out_file_name =="")

        return usage_decompress();

    Decompressor decompressor(params);    // Load settings and data

    // decompressor.getChrom();              //Obtaining chromosome information.

       
    if(!decompressor.decompressProcess())
        return 1;
        

    return 0;
}
