# GT-Pro - the GenoTyper for PROkaryotes  

ultra-rapid genotyping of prokaryotes in shotgun metagenomes

Rationale:
Large public databases of genome and metagenome sequences contain a wealth of information on the within-species genetic variation and population structure of many microbes, especially those from well-studied environments, like the human gut. Micobiome single nucleotide polymorohisms (SNPs) can be detected in shotgun metagenomics sequencing via alignment to databases of genes or genomes. However, this approach is computationally intensive (requiring high-performance computing) and is only accurate for abundant species (with at least 5x sequencing coverage). 

Gut Database: 
We built a database of common SNPs within the core genomes of species that are prevalent in human stool. Leveraging these SNPs, we designed a "virtual genotyping array". This highly compressed data structure is comprised of 31-basepair (bp) "probes" which contain each allele of 56 million bi-allelic SNPs (from 881 gut species) with 15-bp flanking sequence on either side and which have no known exact match to any other species.   

## What GT-Pro does

GT-Pro utilizes an exact matching algorithm to perform ultra-rapid and accurate genotyping of known SNPs from metagenomes.

To genotype a microbiome, GT-Pro takes as input one or more shotgun metagenomics sequencing libraries in FASTQ format. It returns counts of reads exactly matching each allele of each SNP in a concise table-shaped format for its output, with one row for each bi-allelic SNP site that has exactly 8 fields: species, global position, contig, contig position, allele 1, allele 2 and coverage of allele 1 and coverage of allele 2. The k-mer exact match based genotyping algorithm is optimized for machine specificiations, and it can run on a personal computer.

## Step-by-step tutorial

If you are a first-time user, we recommend to use this [step-by-step tutorial](ExampleTutorial.md) to get started with GT-Pro. 

Feel free to skip this section if you are looking for quick usage or examples.

## Dependencies

* Python3  
* pigz (Optional; A parallel implementation of gzip for modern multi-processor, multi-core machines; https://zlib.net/pigz/)
* lbzip2 (Optional; A free, multi-threaded compression utility with support for bzip2 compressed file format; http://lbzip2.org/)
* lz4 (Optional; Extremely Fast Compression algorithm; http://www.lz4.org)

Note:

* Once compiled sucessfully, GT-Pro relies on Python3 for easy access and interface display.
* If you have input sequencing data in gzip, bzip2 and lz4 format, pigz, lbzip2 and lz4 are required to help GT-Pro better decode these files. GT-Pro will fall back to system default decompressor (e.g. gzip) if these depedencies not detected.

## Installation

First, retrieve a copy of GT-Pro to your local computing environment  

`git clone https://github.com/zjshi/gt-pro.git`  

Change your current working directory into where you put GT-Pro  
`cd /path/to/gt-pro/`  

Type in the command line to compile the source code of GT-Pro  
`make`  

Type in the command line to make GT-Pro ready to execute  
`chmod 755 GT_Pro`  

The main program (`GT_Pro`) should be found in the same directory as `/path/to/gt-pro/`. The GT-Pro can be added to the system path so that the main program can be accessed from anywhere. Reference through full path is also allowed.  

<b>Notes for C++ compiler</b>

GT-Pro requires a C++ compiler that is compatible with C++ 11 standards to work properly. All the tests have been done and passed with clang-900.0.38, but it should be compatible for GNU C Compiler (newer than 5.4.0). We have not tested GT-Pro with older compilers, but we expect it to run similiarly as long as it compiles successfully.

## Dowloading default database 

### Main species-specific k-mer database

`wget --content-disposition https://fileshare.czbiohub.org/s/daK4Wj3N7EwSSrd/download`  
`wget --content-disposition https://fileshare.czbiohub.org/s/XNCzWziB4JydWFH/download`  

Note: 

* Upon the completion, two files can be found `20190723_881species_optimized_db_kmer_index.bin` and `20190723_881species_optimized_db_snps.bin`, which togethers represent the compressed database of species-specific k-mers targeting 56 million common, bi-allelic gut microbiome SNPs. 

* The name of this default database is "20190723_881species".

### SNP dictionary for parsing

`wget --content-disposition https://fileshare.czbiohub.org/s/cYDwCAB539kE5yt/download`  

Note:  

* Upon the completion, one file named "variants_main.covered.hq.snp_dict.tsv" can be found in your current directory.

### Species taxonomy metadata

`wget --content-disposition https://fileshare.czbiohub.org/s/XBzwFpfJpmJpBSQ/download`

Note:  

* Upon the completion, one file named "species_taxonomy_ext.tsv" can be found in your current directory.

## Quick usage:  

### Optimize GT-Pro database before metagenotyping  

`/path/to/GT_Pro optimize -d /path/to/database_name -i /path/to/GT-Pro/test/SRR413665_2.fastq.gz`  

Here -d flag again specifies a complete path to the prefix of database, 

The -i flag specifies the location of a testing input file. Any file in FASTQ format, compressed (.lz4, .gz and .bz2 accepted) or not, should work. In case you don't have such a file at hand, we included a small one in the test directory which is used in this example.

Note: we recommend run `/path/to/GT_Pro optimize` in a new environment or before metagenotyping a large number of samples.

### Metagenotyping samples/metagenomes 

`/path/to/GT_Pro genotype -d /path/to/database_name /path/to/1.fastq[.lz4, .gz, .bz2] /path/to/2.fastq[.lz4, .gz, .bz2] ...`  

Here -d flag again specifies a complete path to the prefix of database, which is similiar to the optimize example above.

`GT_Pro genotype` accepts more than one input in a commmand line. The inputs files should be always in FASTQ format and can be supplied as compressed or uncompressed files. GT-Pro can handle multiple mainstream compression algorithms including LZ4 (.lz4), gzip (.gz) and bzip2 (.bz2) and accepts inputs with mixed file types. 
  
For more flags and advanced usage, simply type in  

`/path/to/GT_Pro genotype`   

or 

`/path/to/GT_Pro genotype -h`

See more use examples below.

### Parsing GT-Pro raw output  

`/path/to/GT_Pro parse --dict /path/to/snp_dict.tsv --in /path/to/GT_Pro/raw/output`  

`/path/to/GT_Pro parse` is a simple utility which converts GT_Pro raw output into a more human-friendly format.

## Output format

### GT-Pro raw output

The GT_Pro raw output file is text-based and in the format of "tab-separated values" (tsv), which may be viewed using a text editor.

It has two fields, including an encoded genotype field and a genotype count field. An example of such looks like the following:

| Genotype         | Count       |
| :---             |    :----:   |
| 10001308280      | 1           |
| 10001318280      | 10          |
| 10001316909      | 2           |
| 1036811486637    | 4           |
| 10370212571198   | 1           |
| 1005620676855    | 5           |
| 1024781625558    | 6           |
| ...              | ...         |

The genotype code consists of digits only. Its first 6 digits represent a species ID, the next digit (7th) represents it is a reference or alternative allele and starting from the 8th digit to the end represent the global pos of a SNP site in a species.

We note the raw output file itself is ready for further analysis involving genotypes within or across samples, which may or may not involves linking the decoded genotypes to two external files:

1. species taxonomy metadata (species_taxonomy_ext.tsv)
2. SNP dictionary (variants_main.covered.hq.snp_dict.tsv)

see the section of "Dowloading default database" about how to retrieve them.

### GT-Pro parsed output

The GT_Pro parsed output file is also text-based and in the format of "tab-separated values" (tsv), which may be viewed using a text editor.

It has eight fields as the following:

1. Species ID: six digit ID which specifies a species 
2. Global Pos: up to seven digits which specifies the global position of a SNP in a species
3. Contig: string type with arbitary length which specifies the contig of a representative genome where a SNP is from
4. Local Pos: up to seven digits which specifies the local position of a SNP on a contig
5. Allele 1: single character, A, C, G or T, which specifies allele 1 of a SNP
6. Allele 2: similiar as Ref allele but specifies allele 2 of a SNP
7. Allele 1 Cnt: an integer specifying the count of detected allele 1 in a metagenome
8. Allele 2 Cnt: an integer specifying the count of detected allele 2 in a metagenome

An example of such looks like the following:

| Species ID       | Global Pos    | Contig         | Local Pos      | Allele 1       | Allele 2       | Allele 1 Cnt   | Allele 2 Cnt   |
| :---             |    :----:     |    :----:      |     :----:     |    :----:      |    :----:      |    :----:      |    :----:      |
| 100099           | 349759        | Ga0310508_101  | 349759         | C              | T              | 7              | 0              |
| 100099           | 472876        | Ga0310508_102  | 20713          | C              | T              | 8              | 0              |
| 101349           | 140977        | NZ_GG770218.1  | 131457         | A              | T              | 0              | 11             |
| 101349           | 150940        | NZ_GG770217.1  | 4457           | G              | A              | 11             | 0              |
| 102506           | 1937345       | QFSG01000067   | 4553           | T              | C              | 2              | 0              |
| 103681           | 982849        | AM36-9BH.Scaf4 | 151893         | A              | G              | 5              | 5              |
| 103702           | 1756265       | .14207_7_45.8  | 101338         | C              | A              | 0              | 20             |
| ...              | ...           | ...            | ...            |  ...           |  ...           |  ...           |  ...           |

## More regular examples

### The simplest use:
`GT_Pro genotype -d /path/to/dbname path/to/input.fastq.gz`

Note: 

* By default the output file can be found in your current directory and named to path_to_input__gtpro__dbname.tsv.lz4

* A long name like path_to_input__gtpro__dbname.tsv.lz4 is a purposeful design to help avoid accidental file overwriting.

* The file type .tsv.lz4 reveals that the output is in the format of tab-separated values and compressed with lz4 algorithm.

Question: how to find the dbname?  
A: In a directory hosting a GT-Pro database, highly likily in your current directory, you may see two or all of four following files:
* `[dbname]_optimized_db_kmer_index.bin`
* `[dbname]_optimized_db_snps.bin`
* `[dbname]_optimized_db_lmer_index_xx.bin`, xx is a two-digit number, e.g. 30
* `[dbname]_optimized_db_mmer_bloom_yy.bin`, yy is also a two-digit number, e.g. 35

dbname is all the characters in the brackets while brakets itself not included.

Question: how many CPUs does GT-Pro use by default?  
A: By defauly, GT-Pro automatically detects the number of available CPUs in a computing environment and then uses all of them

Question: why does GT-Pro mentioning 'Skipped x input files due to pre-existing results.'  
A: GT-Pro automatically detects whether there is a conflict on the location it is about to write an output file. There a conflict exists, GT-Pro will by default cautiously skip output writing. You may have a quick check and manually resolve the confict, or use the -f flag to overwrite pre-existing files.

### The simplest use + use a certain number of CPUs:
`GT_Pro genotype -d /path/to/dbname -t 8 path/to/input.fastq.gz`

Note:

* You may consult this example if you do not want to let GT-Pro use all of CPUs.

* The flag of -t can be used for designating a maximum number of CPUs GT-Pro uses.

* For performance reasons, we recommend never supplying -t with a number more than the total number of CPUs in a computing environment

### Genotype more than one sample/metagenome:
`GT_Pro genotype -d /path/to/dbname path/to/input1.fastq.gz path/to/input2.fastq.gz path/to/input3.fastq.gz ...`

Note: 
* by default all the output files can be found in your current directory and each is named similiarly as the single input use.

### Genotype more than one sample/metagenome with mixed file types (e.g. .gz, .lz4, .bz2):
`GT_Pro genotype -d /path/to/dbname path/to/input1.fastq.gz path/to/input2.fastq.lz4 path/to/input3.fq.bz2 ...`

### GT_Pro does not care if input files are mixed with the compressed or uncompressed as long as they are in FASTQ format:
`GT_Pro genotype -d /path/to/dbname path/to/input1.fastq.gz path/to/input2.fastq path/to/input3.fq ...`

### Genotype more than one sample/metagenome in the same directory:
`GT_Pro genotype -d /path/to/dbname -C /path/to/input/directory test576/r1.fastq.lz4 test576/r2.fq.bz2`

Note: 
* you might want to use -C to avoid super long commandline

### Genotype a lot of samples with the same file type (e.g. gzipped fastq) in the same directory:
`GT_Pro genotype -d /path/to/dbname -C /path/to/input/directory *.fastq.gz`

### Genotype a lot of samples with mixed file type (e.g. .gz, .lz4, .bz2) in the same directory:
`GT_Pro genotype -d /path/to/dbname -C /path/to/input/directory *.fastq.gz *.fastq.bz2 *.fq.lz4`

### Genotype everything in one directory:
`GT_Pro genotype -d /path/to/dbname -C /path/to/input/directory *`

Note: 
* you might want to make sure all files in the directory are compatible with GT-Pro

### Designate an full output path for a sample/metagenome :
`GT_Pro genotype -d /path/to/dbname -o /path/to/output/name path/to/input.fastq.gz`

Note: 
* Please consult this example if you dislike the default way how GT-Pro writes an output file. In the example, the output can be found at /path/to/output/name.tsv.lz4

### Designate a different output directory than the current directory for a sample/metagenome :
`GT_Pro genotype -d /path/to/dbname -o /path/to/output/directory/%{in}__gtpro__%{db} path/to/input.fastq.gz`

Note: 

* Suppose that you are fine the way how GT-Pro names an output file but want to keep the current directory clean. In this example, the output can be found at in the designated output directory with a file name like path_to_input__gtpro__dbname.tsv.lz4

* You may also realized that the example of the simplest use is equivalent as `GT_Pro genotype -d /path/to/dbname -o ./%{in}__gtpro__%{db} path/to/input.fastq.gz`

### Designate a different output directory than the current directory when genotyping more than one sample/metagenome :
`GT_Pro genotype -d /path/to/dbname -o /path/to/output/directory/%{in}__gtpro__%{db} path/to/input1.fastq.gz path/to/input2.fastq.gz path/to/input3.fastq.gz`
    
Note: 

* This is similiar as designating different output directory for a single input.

* The output files can be found at in the designated output directory, with names like path_to_input1__gtpro__dbname.tsv.lz4, path_to_input2__gtpro__dbname.tsv.lz4 and path_to_input3__gtpro__dbname.tsv.lz4

## Build customized database

### Additional dependencies

* KMC3  
KMC3 and its installation guidelines can be found at [here](https://github.com/refresh-bio/KMC)

* [Biopython](https://biopython.org/)  
May be installed with `pip install biopython`  

* [Numpy](https://numpy.org/)  
May be installed with `pip install numpy`

### Building usage 

GT-Pro also allows users to build their own database using the `build` subcommand. An example can be found as the following:

* extract test species directories as inputs  
`tar xzvf test/100035.tar.gz && tar xzvf test/101747.tar.gz && tar xzvf test/102779.tar.gz`

* build a three-species GT-Pro database   
`./GT_Pro build --in test/build.list --out ./test/my_db --dbname tri_db --threads 4`

The database building of GT-Pro has a simple interface, which only requires a input file (specified by --in) containing a list of species directories and a output directory (specified by --out) for hosting intermediate and end files related to database building.

In the the example above, GT-Pro reads a input file which contains paths to three species directories and generates a three-species sck-mers database using four threads. Both `test/build.list` and three species directories are included in the test directory.  

Each species directory should contains the following items:
1. a subdirectory named 'genomes' containing whole genomes of strains (n >= 10) in the species
2. a reference genome file named 'reference.fna' 
3. a multiple sequence alignment (MSA) named 'msa.fa', which should be generated by mapping the genomes in 1 to the reference genome in 2
4. a VCF file named 'core_snps.vcf' describing SNP sites in the coordinate of the reference genome in 2 and genotypes across genomes in 1
5. (Optional) a tsv file named 'coords.tsv' describing core genome regions that GT-Pro will use for k-mer extraction.

Both `msa.fa` and `core_snps.vcf` are in standard format. In our paper, we mainly used [MUMmer4](https://mummer4.github.io/) to generate these files. Note that MUMmer4 does not include both files as its standard output and so please expect some addition effort converting MUMmer4 output. If you want to skip that, we recommend using [CallM](https://github.com/zjshi/CallM) which takes 1 and 2 and generates 3, 4 and 5 for you.

The `coords.tsv` has three fields including chrom, start and end and each row specifies a core genome region on a contig (chrom) starting from start to the end. We generated the file again mainly using MUMmer4. This argument is optional, if not supplied, GT-Pro will use all genomic regions when possible.

The --dbname can used for assigning a database name. It is optional and when not supplied, GT-Pro will use a default name of 'sckmer_db'. The --threads flag tells GT-Pro how many threads to use to build the database. 

For more descriptions on other arguments, simply type in

`/path/to/GT_Pro build -h`

Note: building a database for a large number of species may consume a lot of computing resource as well as take long time. We do not recommend to build large databases on a personal computer.

GT-Pro chooses a default k=31 and both database building and metagenotyping code is heavily optimized over this choice. It does not take other k for database building.

## How to cite

If you find GT-Pro helpful, please consider citing our Biorxiv paper for now: Zhou Jason Shi, Boris Dimitrov, Chunyu Zhao, Stephen Nayfach, Katherine S. Pollard; "Ultra-rapid metagenotyping of the human gut microbiome"; bioRxiv 2020.06.12.149336; doi: https://doi.org/10.1101/2020.06.12.149336
