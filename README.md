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
Usage of gsc compress:

        gsc compress [options] [--in [in_file]] [--out [out_file]]

Where:

        [options]              Optional flags and parameters for compression.
        -i,  --in [in_file]    Specify the input file (default: VCF or VCF.GZ). If omitted, input is taken from standard input (stdin).
        -o,  --out [out_file]  Specify the output file. If omitted, output is sent to standard output (stdout).

Options:

        -M,  --mode_lossly     Choose lossy compression mode (lossless by default).
        -b,  --bcf             Input is a BCF file (default: VCF or VCF.GZ).
        -p,  --ploidy [X]      Set ploidy of samples in input VCF to [X] (default: 2).
        -t,  --threads [X]     Set number of threads to [X] (default: 1).
        -d,  --depth [X]       Set maximum replication depth to [X] (default: 100, 0 means no matches).
        -m,  --merge [X]       Specify files to merge, separated by commas (e.g., -m chr1.vcf,chr2.vcf), or '@' followed by a file containing a list of VCF files (e.g., -m @file_with_IDs.txt). By default, all VCF files are compressed.
```
- Decompress / Query
```bash
Usage of gsc decompress and query:

        gsc decompress [options] --in [in_file] --out [out_file]

Where:
        [options]              Optional flags and parameters for compression.
        -i,  --in [in_file]    Specify the input file . If omitted, input is taken from standard input (stdin).
        -o,  --out [out_file]  Specify the output file (default: VCF). If omitted, output is sent to standard output (stdout).

Options:

    General Options:

        -M,  --mode_lossly      Choose lossy compression mode (default: lossless).
        -b,  --bcf              Output a BCF file (default: VCF).

    Filter options (applicable in lossy compression mode only): 

        -r,  --range [X]        Specify range in format [start],[end] (e.g., -r 4999756,4999852).
        -s,  --samples [X]      Samples separated by comms (e.g., -s HG03861,NA18639) OR '@' sign followed by the name of a file with sample name(s) separated by whitespaces (for exaple: -s @file_with_IDs.txt). By default all samples/individuals are decompressed. 
        --header-only           Output only the header of the VCF/BCF.
        --no-header             Output without the VCF/BCF header (only genotypes).
        -G,  --no-genotype      Don't output sample genotypes (only #CHROM, POS, ID, REF, ALT, QUAL, FILTER, and INFO columns).
        -C,  --out-ac-an        Write AC/AN to the INFO field.
        -S,  --split            Split output into multiple files (one per chromosome).
        -I, [ID=^]              Include only sites with specified ID (e.g., -I "ID=rs6040355").
        --minAC [X]             Include only sites with AC <= X.
        --maxAC [X]             Include only sites with AC >= X.
        --minAF [X]             Include only sites with AF >= X (X: 0 to 1).
        --maxAF [X]             Include only sites with AF <= X (X: 0 to 1).
        --min-qual [X]          Include only sites with QUAL >= X.
        --max-qual [X]          Include only sites with QUAL <= X.
```
## Example
There is an example VCF/VCF.gz/BCF file, `toy.vcf`/`toy.vcf.gz`/`toy.bcf`, in the toy folder, which can be used to test GSC
### compress

lossless compression:
```bash
The input file format is VCF:
./gsc compress -o toy/toy_lossless toy/toy.vcf
```
This will create a file:
* `toy_lossless.gsc` - The compressed archive of the entire VCF file.

lossy compression:
```bash
The input file format is VCF:
./gsc compress -M -o toy/toy_lossy toy/toy.vcf
```
This will create a file:
* `toy_lossy.gsc` - The compressed archive of the entire VCF file.
### Decompress
lossless decompression:

To decompress the compressed toy_lossless (which includes: toy_lossless.gti and toy_lossless.dbs) into a VCF file named toy_lossless_decomp.vcf:
```bash
./gsc decompress -o toy/toy_lossless_decomp toy/toy_lossless
```
lossy decompression:

To decompress the compressed toy_lossy (which include: toy_lossy.gti) into a VCF file named toy_lossy_decomp.vcf:
```bash
./gsc decompress -M -o toy/toy_lossy_decomp toy/toy_lossy
```
