# GT-Pro - the GenoTyper for PROkaryotes  

ultra-rapid genotyping of prokaryotes in shotgun metagenomes

Rationale:
Large public databases of genome and metagenome sequences contain a wealth of information on the within-species genetic variation and population structure of many microbes, especially those from well-studied environments, like the human gut. Micobiome single nucleotide polymorohisms (SNPs) can be detected in shotgun metagenomics sequencing via alignment to databases of genes or genomes. However, this approach is computationally intensive (requiring high-performance computing) and is only accurate for abundant species (with at least 5x sequencing coverage). 

Gut Database: 
We built a database of common SNPs within the core genomes of species that are prevalent in human stool. Leveraging these SNPs, we designed a "virtual genotyping array". This highly compressed data structure is comprised of 31-basepair (bp) "probes" which contain each allele of 56 million bi-allelic SNPs (from 881 gut species) with 15-bp flanking sequence on either side and which have no known exact match to any other species.   

## What GT-Pro does

GT-Pro utilizes an exact matching algorithm to perform ultra-rapid and accurate genotyping of known SNPs from metagenomes.

To genotype a microbiome, GT-Pro takes as input one or more shotgun metagenomics sequencing libraries in FASTQ format. It returns counts of reads exactly matching each allele of each SNP in a concise table-shaped format for its output, with one row for each bi-allelic SNP site that has exactly 8 fields: species, SNP ID, contig, contig position, allele 1, allele 2 and coverage of allele 1 and coverage of allele 2. The k-mer exact match based genotyping algorithm is optimized for machine specificiations, and it can run on a personal computer.

## Dependencies

Once compiled sucessfully, GT-Pro does not require any hard dependencies to run. The dependencies listed are only optional and used for specific cases.

If you have input sequenncing data in gzip, bzip2 and lz4 format, the following dependencies are required to help GT-Pro decode these files.
* pigz (A parallel implementation of gzip for modern multi-processor, multi-core machines; https://zlib.net/pigz/)
* lbzip2 (A free, multi-threaded compression utility with support for bzip2 compressed file format; http://lbzip2.org/)
* lz4 (Extremely Fast Compression algorithm; http://www.lz4.org)

## Installation

First, retrieve a copy of GT-Pro to your local computing environment  

`git clone https://github.com/zjshi/gt-pro.git`  

Change your current working directory into where you put GT-Pro  
`cd /path/to/gt-pro/`  

Type in the command line to compile the source code of GT-Pro  
`make`  

Two binary files should be found in the same directory as /path/to/gt-pro/, they are sckmerdb_build and gt_pro. The programs can be put under your favorite system path or directly referenced through full path.  

<b>Notes for C++ compiler</b>

gt-pro requires C++ compiler to work properly. The compiler should be compatible with C++ 11 standards. All the tests have been done and passed with clang-900.0.38, but it should be compatible for GNU C Compiler (newer than 5.4.0).

## Database

### Download using wget
`wget http://jason.shi-bucket.s3.amazonaws.com/public/20190723_881species_optimized_db_kmer_index.bin`  
`wget http://jason.shi-bucket.s3.amazonaws.com/public/20190723_881species_optimized_db_snps.bin`  
`wget http://jason.shi-bucket.s3.amazonaws.com/public/species_taxonomy_ext.tsv`  
`wget http://jason.shi-bucket.s3.amazonaws.com/public/variants_main.covered.hq.snp_dict.tsv`  

### Alternatively, download using aws tools (faster)
`aws s3 cp s3://microbiome-bdimitrov/gt-pro2.0/databases/20190723_881species/20190723_881species_optimized_db_kmer_index.bin ./`  
`aws s3 cp s3://microbiome-bdimitrov/gt-pro2.0/databases/20190723_881species/20190723_881species_optimized_db_snps.bin ./`  
`aws s3 cp s3://jason.shi-bucket/public/variants_main.covered.hq.snp_dict.tsv ./` 
`aws s3 cp s3://jason.shi-bucket/public/species_taxonomy_ext.tsv`  

These downloads currently add up to 13GB and represent the compressed database of 56 million common, bi-allelic gut microbiome SNPs. 

### Finally, create an optimized version of the database for your machine  

`/path/to/gt_pro -d /path/to/database_prefix </dev/null`  

xx and yy should be a double digit number representing size of l and m filter in bits, respectively. These are the suffix array and bloom filter constructed for the ideal L and M parameter values for your machine's RAM size, according to our benchmark information.   

These two values are important to ensure the best computing performance and RAM usage. They are automatically determined by GT-Pro along with the configuration of the right number of threads, the ideal index size, and other parameters. GT-Pro does this automatically for you. 


## Quick usage:  

#### metagenotyping in sequencing samples  

`/path/to/gt_pro -d /path/to/database_prefix -C /path/to/my_inputs 1.fastq[.lz4, .gz, .bz2], 2.fastq....`  

Here -d flag specifies a complete path to the prefix of database and -C flag specifies the location GT-Pro will look for input files given next to the location.  

The outputs will be deposited in the current directory and automatically named so that each output name will contain a portion of the input name and a portion of the database name.  

If you prefer the style of numbered outputs, you may obtain that via the flag -f -o out.%{n}. That is a powerful flag, documented in the help text for the gtpro executable.  

For more flags and advanced usage, simply type in  

`/path/to/gt_pro`

#### decompress output file

`gunzip ./path/to/gt_pro_raw_output` 

#### parse gt-pro raw output  

`python3 ./script/gtp_parse.py --dict /path/to/snp_dict.tsv --in </path/to/gt_pro_output> --v2`  


 
### Example: A Mini test case:

To verify GT-Pro is ready to run, simply try the following command lines

`./gt_pro -d ./20190723_881species -C test/ SRR413665_2.fastq.gz`  
`gunzip ./SRR413665_2__gtpro__20190723_881species.tsv.gz`  
`python3 ./scripts/gtp_parse.py --dict variants_main.covered.hq.snp_dict.tsv --in ./SRR413665_2__gtpro__20190723_881species.tsv --v2`  

If everything runs to finishing without error, parsed GT-Pro will be printed to standard output.

## How to cite

If you find GT-Pro helpful, please consider citing our Biorxiv paper for now: Zhou Jason Shi, Boris Dimitrov, Chunyu Zhao, Stephen Nayfach, Katherine S. Pollard; "Ultra-rapid metagenotyping of the human gut microbiome"; bioRxiv 2020.06.12.149336; doi: https://doi.org/10.1101/2020.06.12.149336
