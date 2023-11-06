
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

    cerr << endl;

    cerr << "Mode options: " << endl;
    
    cerr << "\t-M,  --mode_lossly\t\t choose lossy compression mode (lossless compression mode by default)\t" << endl;

    cerr << endl;

    cerr << "Input\\Output options: " << endl;

    cerr << "\t[in_file]\t\t\t path to input file (a VCF or VCF.GZ file by default)"<< endl;
    
    cerr << "\t-b,  --bcf\t\t\t input is a BCF file (input is a VCF or VCF.GZ file by default)\t" << endl;

    cerr << "\t-p,  --polidy [X]\t\t set ploidy of samples in input VCF to [X] (number >= 1; 2 by default)" << endl;

    cerr << "\t-o,  --out [out_file]\t\t output to a file and set output out_file to [out_file] \t" << endl;
         
    cerr << "\t[out_file]\t\t\t path to output file \n"<< endl;

    cerr << endl;

    cerr << "Parameters options : " << endl;

    cerr << "\t-t,  --threads [X]\t\t set number of threads to [X] (number >= 1; 2 by default)" << endl;

    cerr << "\t-d,  --depth [X]\t\t set the maximum replication depth to [X] (number >= 0; 0 means no matches; 100 by default)" << endl;
   
    cerr << "\t-m,  --merge [X]\t\t [X] separated by comms (for example: -m chr1.vcf,chr2.vcf) OR '@' sign followed by the name of a file with VCF file path separated by whitespaces (for exaple: -m @file_with_IDs.txt). By default all VCF flies are compressed" << endl;
     
    exit(0);
}

int usage_decompress()

{

    cerr << "Decompress and Query usage:" << endl;

    cerr << "\tgsc decompress <options> [out_file] [in_file]" << endl;

    cerr << endl;

    cerr << "Mode options: " << endl;
    
    cerr << "\t-M,  --mode_lossly\t\tchoose lossy compression mode (lossless compression mode by default)\t" << endl;

    cerr << endl;

    cerr << "Input\\Output options: " << endl;

    cerr << "\t[in_file]\t\t\tpath to input file (prefix of the file name to be decompressed)"<< endl;

    cerr << "\t-b,  --bcf\t\t\toutput a BCF file and please use it together with param '-o' (output is a VCF file by default)\t" << endl;
    
    cerr << "\t-o,  --out [out_file]\t\toutput to a file and set output out_file to [out_file] \t" << endl;
    
    cerr << "\t[out_file]\t\t\tyou need to enter the output file path "<< endl;

    cerr << "\th\\H, --header-only\\--no-header\tonly output the header\\don't output the header (only genotypes)" << endl;

    cerr << "\t-G,  --no-genotype\t\tdon't output sample genotypes (only #CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO columns)" << endl;

    cerr << "\t-C,  --out-ac-an\t\twrite AC/AN to the INFO field (always set when using -minAC, -maxAC, -minAF or -maxAF)"<< endl;

    cerr << "\t-S,  --split\t\t\tsplit output into multiple files (one per chromosome)" << endl;

    cerr << endl;

    cerr << "Filter options:: " << endl;
    
    cerr << "\t-r,  --range [X]\t\trange in format [start],[end] (for example: -r 4999756,4999852). By default all variants are decompressed." << endl;

    cerr << "\t-s,  --samples [X]\t\tsamples separated by comms (for example: -s HG03861,NA18639) OR '@' sign followed by the name of a file with sample name(s) separated by whitespaces (for exaple: -s @file_with_IDs.txt). By default all samples/individuals are decompressed" << endl;

    // cerr << "\t-n X \t- process at most X records (by default: all from the selected range)" << endl;

    cerr << "\t--minAC [X] \t\t\treport only sites with count of alternate alleles among selected samples smaller than or equal to X (default: no limit)" << endl;

    cerr << "\t--maxAC [X] \t\t\treport only sites with count of alternate alleles among selected samples greater than or equal to X" << endl;

    cerr << "\t--minAF [X] \t\t\treport only sites with allele frequency among selected samples greather than or equal to X (X - number between 0 and 1; default: 0)" << endl;

    cerr << "\t--maxAF [X] \t\t\treport only sites with allele frequency among selected samples smaller than or equal to X (X - number between 0 and 1; default: 1)" << endl;

    cerr << "\t--min-qual [X] \t\t\treport only sites with QUAL greater than or equal to X (default: 0)" << endl;

    cerr << "\t--max-qual [X] \t\t\treport only sites with QUAL smaller than or equal to X (default: 1000000)" << endl;

    cerr << "\t-i [ID=^] \t\t\treport only sites with ID equal to ID(for example: -i \"ID=rs6040355\")(default: all)" << endl;

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

            if (strcmp(argv[i], "--mode_lossly") == 0 || strcmp(argv[i], "-M") == 0)

                params.compress_mode = compress_mode_t::lossly_mode;

            else if (strcmp(argv[i], "--out") == 0 || strcmp(argv[i], "-o") == 0){

                params.out_file_flag = true;

                i++;

                if (i >= argc)
                    return usage_decompress();

                params.out_file_name = string(argv[i]);
            }
            else if (strcmp(argv[i], "--bcf") == 0 || strcmp(argv[i], "-b") == 0)

                params.in_type = file_type::BCF_File;

            else if (strcmp(argv[i], "--merge") == 0 || strcmp(argv[i], "-m") == 0){

                params.merge_file_flag = true;

                i++;

                if (i >= argc)
                    return usage_compress();     

                params.in_file_name = string(argv[i]);
                cout<<params.in_file_name<<endl;
            }
            else if (strcmp(argv[i], "--ploidy") == 0 || strcmp(argv[i], "-p") == 0){

                i++;

                if (i >= argc)
                    return usage_compress();

                temp = atoi(argv[i]);

                if (temp < 1)
                    usage_compress();

                params.ploidy = temp;
            }
            else if (strcmp(argv[i], "--depth") == 0 || strcmp(argv[i], "-d") == 0){
                
                i++;

                if (i >= argc)
                    return usage_compress();

                temp = atoi(argv[i]);

                if (temp < 0)
                    usage_compress();

                params.max_replication_depth = temp;
            }
            else if (strcmp(argv[i], "--threads") == 0 || strcmp(argv[i], "-t") == 0){

                i++;

                if (i >= argc)
                    return usage_compress();

                temp = atoi(argv[i]);

                if (temp < 1) 
                    usage_compress();

                params.no_threads = temp;
            }
        }
        if (i > argc)
            return usage_compress();

        if(params.out_file_flag == false)
            params.out_file_name = string(argv[argc - 1]);
        
        if(!params.merge_file_flag)
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
            if (strcmp(argv[i], "--mode_lossly") == 0 || strcmp(argv[i], "-M") == 0)

                params.compress_mode = compress_mode_t::lossly_mode;

            else if (strcmp(argv[i], "--out") == 0 || strcmp(argv[i], "-o") == 0){

                params.out_file_flag = true;

                i++;

                if (i >= argc)
                    return usage_decompress();

                params.out_file_name = string(argv[i]);
            }

            else if (strcmp(argv[i], "--bcf") == 0 || strcmp(argv[i], "-b") == 0)

                params.out_type = file_type::BCF_File;

            else if (strcmp(argv[i], "--make-bed") == 0)

                params.out_type = file_type::BED_File;    

            else if (strcmp(argv[i], "--threads") == 0 || strcmp(argv[i], "-t") == 0){

                i++;

                if (i >= argc)
                    return usage_decompress();

                temp = atoi(argv[i]);

                if (temp < 1) 
                    usage_decompress();

                params.no_threads = temp;
            }
            else if (strcmp(argv[i], "--level") == 0 || strcmp(argv[i], "-l") == 0){
                
                i++;
                if(i >= argc)
                    return usage_decompress();
                temp = atoi(argv[i]);
                if(temp < 0 || temp > 9)
                    return usage_decompress();
                else
                {
                    if(temp)
                        params.compression_level = argv[i][0];
                    else
                        params.compression_level = 'u';
                }
            }
         
            else if (strcmp(argv[i], "--split") == 0 || strcmp(argv[i], "-S") == 0)

                params.split_flag = true ;

            else if (strcmp(argv[i], "--samples") == 0 || strcmp(argv[i], "-s") == 0){

                i++;

                params.samples = string(argv[i]);
            }

            else if (strcmp(argv[i], "--range") == 0 || strcmp(argv[i], "-r") == 0){

                i++;

                params.range = string(argv[i]);
            }

            else if (strcmp(argv[i], "-n") == 0){

                i++;

                if (i >= argc)
                    return usage_decompress();

                temp = atoi(argv[i]);

                if (temp < 0)
                    usage_decompress();
                
                else
                    params.records_to_process = temp;
            
            }
            else if (strcmp(argv[i], "-O") == 0){
                
                params.out_samples_name = true;
                
                i++;

                if (i >= argc)
                    return usage_decompress();

                params.out_samples_file_name = string(argv[i]);

            }

            else if (strcmp(argv[i], "--out-ac-an") == 0 || strcmp(argv[i], "-C") == 0)

                params.out_AC_AN = true;

            else if (strcmp(argv[i], "--no-genotype") == 0 || strcmp(argv[i], "-G") == 0)

                params.out_genotypes = false;
            
            else if (strcmp(argv[i], "--no-header") == 0 || strcmp(argv[i], "-H") == 0)

                params.out_header_flag = false;

            else if (strcmp(argv[i], "--header-only") == 0 || strcmp(argv[i], "-h") == 0)

                params.out_header_flag = true;

            else if (strcmp(argv[i], "-i") == 0){

                i++;

                if (i >= argc)
                    return usage_decompress();
                if(strncmp(argv[i], "ID=",3) != 0)
                    return usage_decompress();
                else
                    params.out_id = argv[i]+3;
            }

            else if (strcmp(argv[i], "--min-qual") == 0){

                i++;

                if (i >= argc)
                    return usage_decompress();

                temp_f = atoi(argv[i]);

                params.min_qual = temp_f;

            }
            else if (strcmp(argv[i], "--max-qual") == 0){

                i++;

                if (i >= argc)
                    return usage_decompress();

                temp_f = atoi(argv[i]);

                params.max_qual = temp_f;

            }

            else if (strcmp(argv[i], "--minAC") == 0){

                i++;

                if (i >= argc)
                    return usage_decompress();

                temp = atoi(argv[i]);

                if (temp < 0)
                    usage_decompress();

                params.minAC = temp;

                params.out_AC_AN = true;
            }

            else if (strcmp(argv[i], "--maxAC") == 0){

                i++;

                if (i >= argc)
                    return usage_decompress();

                temp = atoi(argv[i]);

                if (temp < 0)
                    usage_decompress();

                params.maxAC = temp;

                params.out_AC_AN = true;
            }

            else if (strcmp(argv[i], "--minAF") == 0){

                i++;

                if (i >= argc)
                    return usage_decompress();

                temp_f = atof(argv[i]);

                if (temp_f < 0.0 || temp_f > 1.0)
                    usage_decompress();

                params.minAF = temp_f;

                params.out_AC_AN = true;
            }

            else if (strcmp(argv[i], "--maxAF") == 0){

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
