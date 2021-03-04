# Step-by-step Tutorial

## Purpose of this tutorial

Provide a step-by-step tutorial, which helps starter to get familiar with GT-Pro.

After following through this tutorial, a user should have a version of GT-Pro that is installed, optimized according to local computing environment and ready to run.

## Downloading notes

This tutorial involve multiple downloads using aws cli (https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-welcome.html) which is fast and recommended. However, if you don't have access to aws account and cli, the downloads can be achieved with more standard tools, e.g. wget, through http addresses.

## Step 1: clone GT-Pro repo to your favorite location and change current directory to it

<b>if you already followed through installation section, please skip this step and jump to step 2</b>

`git clone https://github.com/zjshi/gt-pro.git`  
`cd /path/to/gt-pro && make`  
`chmod 755 ./GT_Pro` 

Upon the successful completion of the steps above, a main program GT_Pro can be found in the current directory.

Feel free to add /path/to/gt-pro to your system path so that GT_Pro can be reference from any location, but this is not required. GT-Pro can be used normally by referencing its main program GT_Pro by full ## Step 2: download raw material. The raw material consists of two binary files stored in a AWS S3 bucket and can be downloaded with the following command

## Step 2: download raw material. The raw material consists of two binary files stored in a nextcloud location and can be downloaded with the following command

`wget --content-disposition https://fileshare.czbiohub.org/s/daK4Wj3N7EwSSrd/download`  
`wget --content-disposition https://fileshare.czbiohub.org/s/XNCzWziB4JydWFH/download`  

This will emit the following files to the same directory  

`20190723_881species_optimized_db_kmer_index.bin`  
`20190723_881species_optimized_db_snps.bin`  

which currently add up to 13GB and represent the compressed database of 56 million common, bi-allelic gut microbiome SNPs.  

Upon the completion of this step, note that you successfully downloaded a GT-Pro database with a name "20190723_881species", which is ready for the next step. 

## Step 3: create an optimized version of the DB that was produced in its raw ".bin" format by step 2. This step can be performed using the following sample command line:

`./GT_Pro optimize --db 20190723_881species --in test/SRR413665_2.fastq.gz --out test/ `  

and it will take a while (several minutes on AWS c5.18xlarge instance) but will emit the following two files to the same directory that are named something like:  

`20190723_881species_optimized_db_mmer_bloom_xx.bin`  
`20190723_881species_optimized_db_lmer_index_yy.bin`  

xx and yy should be a double digit number (e.g. 32) representing size of l and m in bits for the L-index and M-filter, respectively. 

These are two important data structures that have a non-trivial impact on GT-Pro's performance. 

They are automatically determined by GT-Pro along with the configuration of the right number of threads, the ideal index size, and other parameters. 

If you are interested, please consult our paper for more details. Not understanding these two values does not affect the regular use of GT-Pro

## Step 4. Having completed the above steps, you may now process a test file in gzipped FASTQ format, with a command like

`./GT_Pro genotype -d ./20190723_881species -C ./test SRR413665_2.fastq.gz`  

Here -d flag specifies a complete path to the prefix of database and -C flag specifies the location GT-Pro will look for input files given next to the location.  
For more flags, examples and advanced usage, simply type in  

`./GT_Pro genotype -h`

The outputs will be deposited in the current directory and automatically named so that each output name will contain a portion of the input name and a portion of the database name. 

In this example, you may find an output file named SRR413665_2__gtpro__20190723_881species.tsv.gz in your current directory. 

Note that the raw output file could be in compressed format (.gz), so decompress it using the following command if necessary  

`gunzip ./SRR413665_2__gtpro__20190723_881species.tsv.gz`  

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

can be retrieved using 

`wget --content-disposition https://fileshare.czbiohub.org/s/XBzwFpfJpmJpBSQ/download` 

2. SNP dictionary (variants_main.covered.hq.snp_dict.tsv)

can be retrieved using

`wget --content-disposition https://fileshare.czbiohub.org/s/cYDwCAB539kE5yt/download`

## Step 5: (Optional) GT-Pro raw output is encoded, we provide a small utility to parse the raw output  

First, you may have already dowloaded the parsing dictionary (variants_main.covered.hq.snp_dict.tsv) 

If not, it can be retrieved using the following command line  

`wget --content-disposition https://fileshare.czbiohub.org/s/cYDwCAB539kE5yt/download`

Note that the raw output file could be in compressed format (.gz) by default, so decompress it using the following command if necessary  

`gunzip ./SRR413665_2__gtpro__20190723_881species.tsv.gz`  

Next, parse the uncompressed output file with the following command   

`./GT_Pro parse --dict variants_main.covered.hq.snp_dict.tsv --in ./SRR413665_2__gtpro__20190723_881species.tsv`  

./GT_Pro parse is a simple utility which converts GT_Pro raw output into a more human-friendly format.  

If everything runs to finishing without error, parsed GT-Pro output will be printed to standard output.

The GT_Pro parsed output file is also text-based and in the format of "tab-separated values" (tsv), which may be viewed using a text editor.

It has two fields, including an encoded genotype field and a genotype count field. An example of such looks like the following:

| Species ID       | Global Pos    | Contig         | Local Pos      | Ref Allele     | Alt Allele     | Ref Allele Cnt | Alt Allele Cnt |
| :---             |    :----:     |    :----:      |     :----:     |    :----:      |    :----:      |    :----:      |    :----:      |
| 100099           | 349759        | Ga0310508_101  | 349759         | C              | T              | 7              | 0              |
| 100099           | 472876        | Ga0310508_102  | 20713          | C              | T              | 8              | 0              |
| 101349           | 140977        | NZ_GG770218.1  | 131457         | A              | T              | 0              | 11             |
| 101349           | 150940        | NZ_GG770217.1  | 4457           | G              | A              | 11             | 0              |
| 102506           | 1937345       | QFSG01000067   | 4553           | T              | C              | 2              | 0              |
| 103681           | 982849        | AM36-9BH.Scaf4 | 151893         | A              | G              | 5              | 5              |
| 103702           | 1756265       | .14207_7_45.8  | 101338         | C              | A              | 0              | 20             |
| ...              | ...           | ...            | ...            |  ...           |  ...           |  ...           |  ...           |
