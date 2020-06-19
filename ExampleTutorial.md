# Step-by-step Tutorial

## Purpose of this tutorial

Provide a step-by-step tutorial, which helps starter to get familiar with GT-Pro.

After following through this tutorial, a user should have a version of GT-Pro that is installed, optimized according to local computing environment and ready to run.

## Downloading notes

This tutorial involve multiple downloads using aws cli (https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-welcome.html) which is fast and recommended. However, if you don't have access to aws account and cli, the downloads can be achieved with more standard tools, e.g. wget, through http addresses.

## Step 1: clone GT-Pro repo to your favorite location and change current directory to it

<b>if you already followed through installation section, please skip this step and jump to step 2</b>

`git clone https://github.com/zjshi/gt-pro.git`  
`cd /path/to/gt-pro`  
`make`  

## Step 2: download raw material. The raw material consists of two binary files stored in a AWS S3 bucket and can be downloaded with the following command

<b>Before that make sure (use pwd to find out) you are already in the gt-pro root directory. Change your current directory if necessary. </b>

`wget http://jason.shi-bucket.s3.amazonaws.com/public/20190723_881species_optimized_db_kmer_index.bin`  
`wget http://jason.shi-bucket.s3.amazonaws.com/public/20190723_881species_optimized_db_snps.bin`  

Or (faster but require aws cli)

`aws s3 cp s3://jason.shi-bucket/public/20190723_881species_optimized_db_kmer_index.bin ./`  
`aws s3 cp s3://jason.shi-bucket/public/20190723_881species_optimized_db_snps.bin ./`  

This will emit the following files to the same directory  

`20190723_881species_optimized_db_kmer_index.bin`  
`20190723_881species_optimized_db_snps.bin`  

which currently add up to 13GB and represent the compressed database of 56 million common, bi-allelic gut microbiome SNPs.  

## Step 3: create an optimized version of the DB that was produced in its raw ".bin" format by step 2. This step can be performed using the following sample command line:

`./gt_pro -d 20190723_881species </dev/null`  

and it will take a while (several minutes on AWS c5.18xlarge instance) but will emit the following two files to the same directory that are named something like:  

`20190723_881species_optimized_db_mmer_bloom_xx.bin`  
`20190723_881species_optimized_db_lmer_index_yy.bin`  

xx and yy should be a double digit number representing size of l and m filter in bits, respectively. These are the suffix array and bloom filter constructed for the ideal L and M parameter values for your machine's RAM size, according to our benchmark information.   

These two values are important to ensure the best computing performance and RAM usage. They are automatically determined by GT-Pro along with the configuration of the right number of threads, the ideal index size, and other parameters. GT-Pro does this automatically for you. 

## Step 4. Having completed the above steps, you may now process a number of files you wish, with a command like

`./gt_pro -d ./20190723_881species -C ./test SRR413665_2.fastq.gz`  

Here -d flag specifies a complete path to the prefix of database and -C flag specifies the location GT-Pro will look for input files given next to the location.  

The outputs will be deposited in the current directory and automatically named so that each output name will contain a portion of the input name and a portion of the database name.  

For more flags and advanced usage, simply type in  

`./gt_pro`

## Step 5: GT-Pro raw output is hard encoded, we provide a small script to parse the raw output  

First, you may need to dowload the parsing dictionary using the following command line  

`wget http://jason.shi-bucket.s3.amazonaws.com/public/variants_main.covered.hq.snp_dict.tsv`

Or (faster)

`aws s3 cp s3://jason.shi-bucket/public/variants_main.covered.hq.snp_dict.tsv ./`  

Note that the raw output file could be in compressed format (.gz), so decompress it using the following command if necessary  

`gunzip ./SRR413665_2__gtpro__20190723_881species.tsv.gz`  

Next, parse the uncompressed output file with the following command   

`python3 ./scripts/gtp_parse.py --dict variants_main.covered.hq.snp_dict.tsv --in ./SRR413665_2__gtpro__20190723_881species.tsv --v2`  

gtp_parse.py is a parser written in Python script, please see its helper text for more detailed usage.  

If everything runs to finishing without error, parsed GT-Pro will be printed to standard output.

