# gt-pro2.0 

lightweight and rapid tool for prokaryotes genotyping in metagenomes

Rationale:
Large public databases of genome and metagenome sequences contain a wealth of information on the population structure of many microbes. For some well-studied environments, like the human gut, these data should be sufficient to build high resolution maps of genomic variation for common species. Maps of genome variation can in turn be used to identify core-genome regions and common SNPs within the core-genome.

Based on these maps of genomic variation, it should be possible to design a lightweight method to scan new short-read datasets, and genotype the populations.   

## Installation

`cd /path/to/gt-pro2.0/`  
`make`  

Three binary files should be found in the same directory as /path/to/gt-pro2.0/, they are sckmerdb_build, sckmerdb_inspect and gt_pro. The programs should be put under your favorite system path or directly referenced through full path.  

<b>Notes for C++ compiler</b>

gt-pro2.0 requires C++ compiler to work properly. The compiler should be compatible with C++ 11 standards. All the tests have been done and passed with clang-900.0.38, but it should be compatible for GNU C Compiler (newer than 5.4.0).

## Tutorial

program usage:  

build your build sckmerdb from k-mer profile files
`sckmerdb_build fpath [fpath ...]`  

genotyping in metagenome samples  
`gt_pro -d <sckmerdb_path: string> -r <read_len; int; default 90> -t <n_threads; int; default 1> [-h] input1 [input2 ...]`  

inspect sckmerdb  
`sckmerdb_inspect fpath`  

test case:  
`./sckmerdb_build ./test/275577.sckmers.db.tsv ./test/276044.sckmers.db.tsv > ./test/sckmer_2sps.bin`  
`./sckmerdb_inspect ./test/sckmer_2sps.bin | head -n 10`  
`./gt_pro -d ./test/sckmer_2sps.bin <(gzip -dc ./test/SRR413665_2.fastq.gz)`  

## Download sckmerdb

k-mer profile files for all 974 species  
`wget http://jason.shi-bucket.s3.amazonaws.com/sckmerdb/sckmer_profiles.tar.bz2`  

raw k-mer database files for all 974 species  
`wget http://jason.shi-bucket.s3.amazonaws.com/sckmerdb/sckmer_dbs.tar.bz2`

all-in-one encoded k-mer database file  
`wget http://jason.shi-bucket.s3.amazonaws.com/sckmerdb/sckmerdb_sp974.bin`

species taxonomy metadata  
`wget http://jason.shi-bucket.s3.amazonaws.com/sckmerdb/gut_species_taxonomy.tsv`  
