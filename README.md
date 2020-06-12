# GT-Pro - the GenoTyper for PROkaryotes  

lightweight and rapid tool for prokaryotes genotyping in metagenomes

Rationale:
Large public databases of genome and metagenome sequences contain a wealth of information on the population structure of many microbes. For some well-studied environments, like the human gut, these data should be sufficient to build high resolution maps of genomic variation for common species. Maps of genome variation can in turn be used to identify core-genome regions and common SNPs within the core-genome.

Based on these maps of genomic variation, we design a lightweight method to scan new short-read datasets, and genotype the populations.   


## What GT-Pro does

GT-Pro utilizes an exact matching algorithm to perform ultra-rapid and accurate genotyping of known SNPs from metagenomes.

It takes as inputs metagenomic sequencing read files in FASTQ format and uses a concise table-shaped format for its output, in which every row represents a bi-allelic SNP site. Each row has exactly 8 fields: species, SNP ID, contig, contig position, allele 1, allele 2 and coverage of allele 1 and coverage of allele 2.    

## Dependencies

Once compiled sucessfully, GT-Pro does not require any hard dependencies to run. The dependencies listed are only optional and used for specific cases.

If you have input sequenncing data in gzip, bzip2 and lz4 format, the following dependencies are required to help GT-Pro decode these files.
* pigz (A parallel implementation of gzip for modern multi-processor, multi-core machines; https://zlib.net/pigz/)
* lbzip2 (A free, multi-threaded compression utility with support for bzip2 compressed file format; http://lbzip2.org/)
* lz4 (Extremely Fast Compression algorithm; http://www.lz4.org)

## Installation

`git clone https://github.com/zjshi/gt-pro.git`
`cd /path/to/gt-pro/`  
`make`  

Two binary files should be found in the same directory as /path/to/gt-pro/, they are sckmerdb_build and gt_pro. The programs can be put under your favorite system path or directly referenced through full path.  

<b>Notes for C++ compiler</b>

gt-pro requires C++ compiler to work properly. The compiler should be compatible with C++ 11 standards. All the tests have been done and passed with clang-900.0.38, but it should be compatible for GNU C Compiler (newer than 5.4.0).

## Simple Tutorial (for quick reference)

<b>if you are start a trial of GT-Pro, please directly refer to "Step-by-step usage" section</b>

### Necessary downloads:

#### Bit encoded sckmers and sc-spans  

`aws s3 cp s3://microbiome-bdimitrov/gt-pro2.0/databases/20190723_881species/20190723_881species_optimized_db_kmer_index.bin ./`  
`aws s3 cp s3://microbiome-bdimitrov/gt-pro2.0/databases/20190723_881species/20190723_881species_optimized_db_snps.bin ./`  

#### Dictionary for parsing gt-pro raw output  

`aws s3 cp s3://jason.shi-bucket/public/variants_main.covered.hq.snp_dict.tsv ./`  

### Quick usage:  

#### create an optimized version of the database for your machine  

`/path/to/gt_pro -d /path/to/database_prefix </dev/null`  

#### metagenotyping in sequencing samples  

`/path/to/gt_pro -d /path/to/database_prefix -C /path/to/my_inputs 1.fastq[.lz4, .gz, .bz2], 2.fastq....`  

#### parse gt-pro raw output  

`python3 ./script/gtp_parse.py --dict /path/to/snp_dict.tsv --in </path/to/gt_pro_output> --v2`  

## Step-by-step Tutorial (for starter):

### Step 1: clone GT-Pro repo to your favorite location and change current directory to it

`git clone https://github.com/zjshi/gt-pro.git`  
`cd /path/to/gt-pro`  
`make`  

### Step 2: download raw material. The raw material consists of two binary files stored in a AWS S3 bucket and can be downloaded with the following command

`aws s3 cp s3://jason.shi-bucket/public/20190723_881species_optimized_db_kmer_index.bin ./`  
`aws s3 cp s3://jason.shi-bucket/public/20190723_881species_optimized_db_snps.bin ./`  

This will emit the following files to the same directory  

`20190723_881species_optimized_db_kmer_index.bin`  
`20190723_881species_optimized_db_snps.bin`  

which currently add up to 13GB and represent the "optimized" or compressed database.  

### Step 3: create an optimized version of the DB that was produced in its raw ".bin" format by step 2. This step can be performed using the following sample command line:

`/path/to/gt_pro -d 20190723_881species </dev/null`  

and it will take a while (several minutes on AWS c5.18xlarge instance) but will emit the following files to the same directory.  

That exact same command will also, at the same time, emit files that are named something like:  

`20190723_881species_optimized_db_mmer_bloom_xx.bin`  
`20190723_881species_optimized_db_lmer_index_yy.bin`  

xx and yy should be a double digit number representing size of l and m filter in bits, respectively. These are the suffix array and bloom filter constructed for the ideal L and M parameter values for your machine's RAM size, according to our benchmark information.   

These two values are important to the best computing performance and RAM dependent, and it is automatically determined by GT-Pro along with the configuration of the right number of threads, the ideal index size, etc. GT-Pro just does it automatically for you.  

### Step 4. Having completed the above steps, you may now simply process any number of files you wish, with a command like

`/path/to/gt_pro -d /path/to/20190723_881species -C /path/to/my_inputs 1.fastq[.lz4, .gz, .bz2], 2.fastq....`  

Here -d flag specifies a complete path to the prefix of database and -C flag specifies the location GT-Pro will look for input files given next to the location.  

The outputs will be deposited in the current directory and automatically named so that each output name will contain a portion of the input name and a portion of the database name.  

If you prefer the previous style of numbered outputs, you may obtain that via the flag -f -o out.%{n}. That is a powerful flag, documented in the help text for the gtpro executable.  

For more flags and advanced usage, simply type in  

`/path/to/gt_pro`

### Step 5: GT-Pro raw output is hard encoded, we provide a small script to parse the raw output  

First, you may need to dowload the parsing dictionary using the following command line  

`aws s3 cp s3://jason.shi-bucket/public/variants_main.covered.hq.snp_dict.tsv ./`  

Next, use the following the command line for parsing  

`python3 ./scripts/gtp_parse.py --dict variants_main.covered.hq.snp_dict.tsv --in /path/to/gt_pro_raw_output --v2`  

gtp_parse.py is a parser written in Python script, please see its helper text for more detailed usage.  

### A Mini test case:

To verify GT-Pro is good to go at this stage, simply try the following command lines

`./gt_pro -d ./20190723_881species -C test/ SRR413665_2.fastq.gz`  
`gunzip ./SRR413665_2__gtpro__20190723_881species.tsv.gz`  
`python3 ./scripts/gtp_parse.py --dict variants_main.covered.hq.snp_dict.tsv --in ./SRR413665_2__gtpro__20190723_881species.tsv --v2`  

If everything runs to finishing without error, parsed GT-Pro will be printed to standard output.

## Other downloading information

### In general, AWS tools are recommnded for downloading files with large volume.

#### Bit encoded sckmers and sc-spans  
`aws s3 cp s3://jason.shi-bucket/public/20190723_881species_optimized_db_kmer_index.bin ./`  
`aws s3 cp s3://jason.shi-bucket/public/20190723_881species_optimized_db_snps.bin ./`  

#### species taxonomy metadata  
`aws s3 cp s3://jason.shi-bucket/public/species_taxonomy_ext.tsv`  

#### gt-pro raw output parsing dictionary (updated to include chrom and local position fields)  
`aws s3 cp s3://jason.shi-bucket/public/variants_main.covered.hq.snp_dict.tsv`  

### Alternatively, these files can be retrieved with more standard tools, e.g. wget, through http addresses.

You might consume less network bandwidth but a 10X speed drop is expected.

#### Bit encoded sckmers and sc-spans  
`wget http://jason.shi-bucket.s3.amazonaws.com/public/20190723_881species_optimized_db_kmer_index.bin`  
`wget http://jason.shi-bucket.s3.amazonaws.com/public/20190723_881species_optimized_db_snps.bin`  

#### species taxonomy metadata  
`wget http://jason.shi-bucket.s3.amazonaws.com/public/species_taxonomy_ext.tsv`  

#### gt-pro raw output parsing dictionary (updated to include chrom and local position fields)  
`wget http://jason.shi-bucket.s3.amazonaws.com/public/variants_main.covered.hq.snp_dict.tsv`  

