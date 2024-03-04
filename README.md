# GSC (Genotype Sparse Compression)
Genotype Sparse Compression (GSC) is an advanced tool for lossless compression of VCF files, designed to efficiently store and manage VCF files in a compressed format. It accepts VCF/BCF files as input and utilizes advanced compression techniques to significantly reduce storage requirements while ensuring fast query capabilities. In our study, we successfully compressed the VCF files from the 1000 Genomes Project (1000Gpip3), consisting of 2504 samples and 80 million variants, from an uncompressed VCF file of 803.70GB to approximately 1GB.
## Requirements 
### GSC requires:

- **Compiler Compatibility**: GSC requires a modern C++14-ready compiler, such as:
  - g++ version 10.1.0 or higher

- **Build System**: Make build system is necessary for compiling GSC.

- **Operating System**: GSC supports 64-bit operating systems, including:
  - Linux (Ubuntu)
  
## Installation
To download, build and install GSC use the following commands.
```bash
git clone https://github.com/luo-xiaolong/GSC.git
cd GSC
make
```
To clean the GSC build use:
```bash
make clean
```
## Usage
```bash
Usage: gsc [option] [arguments] 
Available options: 
        compress - compress VCF/BCF file
        decompress     - query and decompress to VCF/BCF file
```
- Compress the input VCF/BCF file
```bash
Compress usage: 
        gsc compress <options> [out_file] [in_file]   


Mode options: 
        -M,  --mode_lossly               choose lossy compression mode (lossless compression mode by default)

Input\Output options: 
        [in_file]                        path to input file (a VCF or VCF.GZ file by default)
        -b,  --bcf                       input is a BCF file (input is a VCF or VCF.GZ file by default)
        -p,  --polidy [X]                set ploidy of samples in input VCF to [X] (number >= 1; 2 by default)
        -o,  --out [out_file]            output to a file and set output out_file to [out_file] 
        [out_file]                       path to output file 


Parameters options : 
        -t,  --threads [X]               set number of threads to [X] (number >= 1; 2 by default)
        -d,  --depth [X]                 set the maximum replication depth to [X] (number >= 0; 0 means no matches; 100 by default)
        -m,  --merge [X]                 [X] separated by comms (for example: -m chr1.vcf,chr2.vcf) OR '@' sign followed by the name of a file with VCF file path separated by whitespaces (for exaple: -m @file_with_IDs.txt). By default all VCF flies are compressed
```
- Decompress / Query
```bash
Decompress and Query usage:
        gsc decompress <options> [out_file] [in_file]

Mode options: 
        -M,  --mode_lossly              choose lossy compression mode (lossless compression mode by default)

Input\Output options: 
        [in_file]                       path to input file (prefix of the file name to be decompressed)
        -b,  --bcf                      output a BCF file and please use it together with param '-o' (output is a VCF file by default)
        -o,  --out [out_file]           output to a file and set output out_file to [out_file] 
        [out_file]                      you need to enter the output file path 
        h\H, --header-only\--no-header  only output the header\don't output the header (only genotypes)
        -G,  --no-genotype              don't output sample genotypes (only #CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO columns)
        -C,  --out-ac-an                write AC/AN to the INFO field (always set when using -minAC, -maxAC, -minAF or -maxAF)
        -S,  --split                    split output into multiple files (one per chromosome)

Filter options:: 
        -r,  --range [X]                range in format [start],[end] (for example: -r 4999756,4999852). By default all variants are decompressed.
        -s,  --samples [X]              samples separated by comms (for example: -s HG03861,NA18639) OR '@' sign followed by the name of a file with sample name(s) separated by whitespaces (for exaple: -s @file_with_IDs.txt). By default all samples/individuals are decompressed
        --minAC [X]                     report only sites with count of alternate alleles among selected samples smaller than or equal to X (default: no limit)
        --maxAC [X]                     report only sites with count of alternate alleles among selected samples greater than or equal to X
        --minAF [X]                     report only sites with allele frequency among selected samples greather than or equal to X (X - number between 0 and 1; default: 0)
        --maxAF [X]                     report only sites with allele frequency among selected samples smaller than or equal to X (X - number between 0 and 1; default: 1)
        --min-qual [X]                  report only sites with QUAL greater than or equal to X (default: 0)
        --max-qual [X]                  report only sites with QUAL smaller than or equal to X (default: 1000000)
        -i [ID=^]                       report only sites with ID equal to ID(for example: -i "ID=rs6040355")(default: all)
```
## Example
There is an example VCF/VCF.gz/BCF file, toy.vcf/toy.vcf.gz/toy.bcf, in the toy_ex folder, which can be used to test GSC
### compress

lossless compression:
```bash
The input file format is VCF:
./gsc compress -o toy/toy_compress_result toy/toy.vcf
```
lossly compression:
```bash
The input file format is VCF:
./gsc compress -M -o toy/toy_compress_result toy/toy.vcf

```
### Decompress
lossless decompression:
```bash
./gsc decompress -o toy/toy_decompress_result toy/toy_compress_result
```
lossly decompression:
```bash
./gsc decompress -M -o toy/toy_decompress_result toy/toy_compress_result
```
