
#include <iostream>

#include <unistd.h>

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

std::cerr << "Usage of gsc compress:\n\n"
            << "\tgsc compress [options] [--in [in_file]] [--out [out_file]]\n\n"
            << "Where:\n\n"
            << "\t[options]              Optional flags and parameters for compression.\n"
            << "\t-i,  --in [in_file]    Specify the input file (default: VCF or VCF.GZ). If omitted, input is taken from standard input (stdin).\n"
            << "\t-o,  --out [out_file]  Specify the output file. If omitted, output is sent to standard output (stdout).\n\n"
              
            << "Options:\n\n"
            << "\t-M,  --mode_lossly     Choose lossy compression mode (lossless by default).\n"
            << "\t-b,  --bcf             Input is a BCF file (default: VCF or VCF.GZ).\n"
            << "\t-p,  --ploidy [X]      Set ploidy of samples in input VCF to [X] (default: 2).\n"
            << "\t-t,  --threads [X]     Set number of threads to [X] (default: 1).\n"
            << "\t-d,  --depth [X]       Set maximum replication depth to [X] (default: 100, 0 means no matches).\n"
            << "\t-m,  --merge [X]       Specify files to merge, separated by commas (e.g., -m chr1.vcf,chr2.vcf), or '@' followed by a file containing a list of VCF files (e.g., -m @file_with_IDs.txt). By default, all VCF files are compressed.\n"
            << std::endl;

    exit(0);
}

int usage_decompress()

{
std::cerr << "Usage of gsc decompress and query:\n\n"
          << "\tgsc decompress [options] --in [in_file] --out [out_file]\n\n"

          << "Where:\n"
          << "\t[options]              Optional flags and parameters for compression.\n"
          << "\t-i,  --in [in_file]    Specify the input file . If omitted, input is taken from standard input (stdin).\n"
          << "\t-o,  --out [out_file]  Specify the output file (default: VCF). If omitted, output is sent to standard output (stdout).\n\n"

          << "Options:\n\n"
          << "    General Options:\n\n"
          << "\t-M,  --mode_lossly\tChoose lossy compression mode (default: lossless).\n"
          << "\t-b,  --bcf\t\tOutput a BCF file (default: VCF).\n\n"

          << "    Filter options (applicable in lossy compression mode only): \n\n"
          << "\t-r,  --range [X]\tSpecify range in format [start],[end] (e.g., -r 4999756,4999852).\n"
          << "\t-s,  --samples [X]\tSamples separated by comms (e.g., -s HG03861,NA18639) OR '@' sign followed by the name of a file with sample name(s) separated by whitespaces (for exaple: -s @file_with_IDs.txt). By default all samples/individuals are decompressed. \n"
          << "\t--header-only\t\tOutput only the header of the VCF/BCF.\n"
          << "\t--no-header\t\tOutput without the VCF/BCF header (only genotypes).\n"
          << "\t-G,  --no-genotype\tDon't output sample genotypes (only #CHROM, POS, ID, REF, ALT, QUAL, FILTER, and INFO columns).\n"
          << "\t-C,  --out-ac-an\tWrite AC/AN to the INFO field.\n"
          << "\t-S,  --split\t\tSplit output into multiple files (one per chromosome).\n"
          << "\t-I, [ID=^]\t\tInclude only sites with specified ID (e.g., -I \"ID=rs6040355\").\n"
          << "\t--minAC [X]\t\tInclude only sites with AC <= X.\n"
          << "\t--maxAC [X]\t\tInclude only sites with AC >= X.\n"
          << "\t--minAF [X]\t\tInclude only sites with AF >= X (X: 0 to 1).\n"
          << "\t--maxAF [X]\t\tInclude only sites with AF <= X (X: 0 to 1).\n"
          << "\t--min-qual [X]\t\tInclude only sites with QUAL >= X.\n"
          << "\t--max-qual [X]\t\tInclude only sites with QUAL <= X.\n"
          << std::endl;

    // cerr << "Decompress and Query usage:" << endl;

    // cerr << "\tgsc decompress <options> [out_file] [in_file]" << endl;

    // cerr << endl;

    // cerr << "Mode options: " << endl;
    
    // cerr << "\t-M,  --mode_lossly\t\tchoose lossy compression mode (lossless compression mode by default)\t" << endl;

    // cerr << endl;

    // cerr << "Input\\Output options: " << endl;

    // cerr << "\t[in_file]\t\t\tpath to input file (prefix of the file name to be decompressed)"<< endl;

    // cerr << "\t-b,  --bcf\t\t\toutput a BCF file and please use it together with param '-o' (output is a VCF file by default)\t" << endl;
    
    // cerr << "\t-o,  --out [out_file]\t\toutput to a file and set output out_file to [out_file] \t" << endl;
    
    // cerr << "\t[out_file]\t\t\tyou need to enter the output file path "<< endl;

    // cerr << "\th\\H, --header-only\\--no-header\tonly output the header\\don't output the header (only genotypes)" << endl;

    // cerr << "\t-G,  --no-genotype\t\tdon't output sample genotypes (only #CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO columns)" << endl;

    // cerr << "\t-C,  --out-ac-an\t\twrite AC/AN to the INFO field (always set when using -minAC, -maxAC, -minAF or -maxAF)"<< endl;

    // cerr << "\t-S,  --split\t\t\tsplit output into multiple files (one per chromosome)" << endl;

    // cerr << endl;

    // cerr << "Filter options:: " << endl;
    
    // cerr << "\t-r,  --range [X]\t\trange in format [start],[end] (for example: -r 4999756,4999852). By default all variants are decompressed." << endl;

    // cerr << "\t-s,  --samples [X]\t\tsamples separated by comms (for example: -s HG03861,NA18639) OR '@' sign followed by the name of a file with sample name(s) separated by whitespaces (for exaple: -s @file_with_IDs.txt). By default all samples/individuals are decompressed" << endl;

    // // cerr << "\t-n X \t- process at most X records (by default: all from the selected range)" << endl;

    // cerr << "\t--minAC [X] \t\t\treport only sites with count of alternate alleles among selected samples smaller than or equal to X (default: no limit)" << endl;

    // cerr << "\t--maxAC [X] \t\t\treport only sites with count of alternate alleles among selected samples greater than or equal to X" << endl;

    // cerr << "\t--minAF [X] \t\t\treport only sites with allele frequency among selected samples greather than or equal to X (X - number between 0 and 1; default: 0)" << endl;

    // cerr << "\t--maxAF [X] \t\t\treport only sites with allele frequency among selected samples smaller than or equal to X (X - number between 0 and 1; default: 1)" << endl;

    // cerr << "\t--min-qual [X] \t\t\treport only sites with QUAL greater than or equal to X (default: 0)" << endl;

    // cerr << "\t--max-qual [X] \t\t\treport only sites with QUAL smaller than or equal to X (default: 1000000)" << endl;

    // cerr << "\t-I [ID=^] \t\t\treport only sites with ID equal to ID(for example: -I \"ID=rs6040355\")(default: all)" << endl;

    // // cerr << "\t-m X\t- limit maximum memory usage to remember previous vectors to X MB (no limit by default)\t" << endl;

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
		    std::cerr << "Compression error!!!\n";
    }

    else if(params.task_mode == task_mode_t::mdecompress){
        result = decompress_entry();
        if (result)
		    std::cerr << " Decompression error!!!\n";
    }

    high_resolution_clock::time_point end = high_resolution_clock::now();

	duration<double> time_duration = duration_cast<duration<double>>(end - start);

	std::cerr << "Total processing time: " << time_duration.count() << " seconds.\n";

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
             
        int temp;
        int i = 2;
           
        for(; i < argc; ++i){

            if (argv[i][0] != '-'){
                return usage_compress();
                break;
            }

            if (strcmp(argv[i], "--mode_lossly") == 0 || strcmp(argv[i], "-M") == 0)

                params.compress_mode = compress_mode_t::lossly_mode;
            
            else if (strcmp(argv[i], "--in") == 0 || strcmp(argv[i], "-i") == 0){
                
                if(!isatty(STDIN_FILENO)){

                    std::cerr << "Error: Conflicting inputs - both filename and stdin data detected." << std::endl;
                    return usage_compress();
                }
                i++;

                if (i >= argc)
                    return usage_compress();

                params.in_file_name = string(argv[i]);
            }                

            else if (strcmp(argv[i], "--out") == 0 || strcmp(argv[i], "-o") == 0){
                
                if(!isatty(STDOUT_FILENO))
                    return usage_compress();
                // params.out_file_flag = true;

                i++;

                if (i >= argc)
                    return usage_compress();

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
        if(isatty(STDIN_FILENO) && params.in_file_name == "-"){

            std::cerr << "Error: No input file specified and no data provided via stdin!" << std::endl;
            return usage_compress();
        }
        if(isatty(STDOUT_FILENO) && params.out_file_name == "-"){

            std::cerr << "Warning: No output file specified and no data provided via stdout!" << std::endl;

        }
        

    }
    else if(params.task_mode == task_mode_t::mdecompress){
        
        


        int temp, i = 2;;
        float temp_f;
        for (; i < argc; ++i){
            
            if (argv[i][0] != '-'){
                usage_decompress();
                break;
            }
            if (strcmp(argv[i], "--mode_lossly") == 0 || strcmp(argv[i], "-M") == 0)

                params.compress_mode = compress_mode_t::lossly_mode;

            else if (strcmp(argv[i], "--in") == 0 || strcmp(argv[i], "-i") == 0){
                
                if(!isatty(STDIN_FILENO)){

                    std::cerr << "Error: Conflicting inputs - both filename and stdin data detected!" << std::endl;
                    return usage_decompress();
                }
                i++;

                if (i >= argc)
                    return usage_decompress();

                params.in_file_name = string(argv[i]);
            }

            else if (strcmp(argv[i], "--out") == 0 || strcmp(argv[i], "-o") == 0){

                // params.out_file_flag = true;
                if(!isatty(STDOUT_FILENO))
                    return usage_decompress();

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

                if (i >= argc)
                    return usage_decompress();               


                params.samples = string(argv[i]);
            }

            else if (strcmp(argv[i], "--range") == 0 || strcmp(argv[i], "-r") == 0){

                i++;

                if (i >= argc)
                    return usage_decompress();

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

            else if (strcmp(argv[i], "-I") == 0){

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

        if(isatty(STDIN_FILENO) && params.in_file_name == "-"){

            std::cerr << "Error: No input file specified and no data provided via stdin!" << std::endl;
            return usage_decompress();
        }
        if(isatty(STDOUT_FILENO) && params.out_file_name == "-"){

            std::cerr << "Warning: No output file specified and no data provided via stdout!" << std::endl;

        }        
        
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
