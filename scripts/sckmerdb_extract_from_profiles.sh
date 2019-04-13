#/bin/bash
#
# This converts kmer profile files into text files as an intermediate step in database generation.
#
# INPUT:  A directory with *.sckmers.tsv files, one per species, with content like:
#
# ...
# 343 25  TGGATACCACGGCGCAAGAGCACGCACAGAA TGGATACCACGGCGCAAGAGCACGCGCAGAA TTCTGTGCGTGCTCTTGCGCCGTGGTATCCA TTCTGCGCGTGCTCTTGCGCCGTGGTATCCA 1   11  259564  8   3
# 343 24  GGATACCACGGCGCAAGAGCACGCACAGAAT GGATACCACGGCGCAAGAGCACGCGCAGAAT ATTCTGTGCGTGCTCTTGCGCCGTGGTATCC ATTCTGCGCGTGCTCTTGCGCCGTGGTATCC 1   11  259564  8   3
# ...
#
# where each line contains, in this order:
#
#      $1 - genomic position of SNP
#      $2 - zero-based offset of SNP in forward 31-mer
#      $3 - forward 31-mer covering the SNP's major allele
#      $4 - forward 31-mer covering the SNP's minor allele
#      $5 - reverse complement of $3
#      $6 - reverse complement of $4
#      $7 - unused
#      $8 - 6 decimal digit species ID
#      $9 - unused
#
# OUTPUT:  A *.db.tsv output file for each corresponding *.sckmer.tsv input file, with content like:
#
# ...
# CGCGTGCTCTTGCGCCGTGGTATCCATAGGC 0   2595641343
# GCCTATGGATACCACGGCGCAAGAGCACGCA 30  2595640343
# GCCTATGGATACCACGGCGCAAGAGCACGCG 30  2595641343
# TGCGTGCTCTTGCGCCGTGGTATCCATAGGC 0   2595640343
# ...
#
# That's 4 lines for each input line; one kmer per line.  The columns are
#
#       $1 - kmer
#       $2 - zero-based offset of SNP within kmer (= input column $2 for the forward kmers, 30 - $2 for RCs).
#       $3 - decimal literal "snp coordinate" obtained by concatenating the following character strings:
#               * the 6 decimal digit species id (input column $8),
#               * a single digit allele type:
#                      (0 = major allele, matching input $3 and $5,
#                       1 = minor allele, matching input $4 and $6),
#               * genomic position of SNP (input column $1)
#
# The output is sorted so that identical kmers (covering different SNPs) are consecutive.
#
# Those DB files can be merged and packed into binary format through the C++ program sckmer_build.cpp, for use in the gtpro query tool.
#

# If we can't figure it out
DEFAULT_NUMCPUS=8

# This works on Linux and Mac OS;  and otherwise will fall back on the default value.
let NUMCPUS=`getconf _NPROCESSORS_ONLN || echo ${DEFAULT_NUMCPUS}`

# ** CAUTION ** The "30" below is assuming 31-mers.  Should be one less than kmer's K.
ls | grep sckmers.tsv | cut -d'.' -f1 | xargs -I{} -n1 -P${NUMCPUS} bash -c 'awk "$0" {}.sckmers.tsv | LC_ALL=C sort -k1 > {}.db.tsv' '{print $3 "\t" $2 "\t" $10 "0" $1 "\n" $4 "\t" $2 "\t" $10 "1" $1 "\n" $5 "\t" (30-$2) "\t" $10 "0" $1 "\n" $6 "\t" (30-$2) "\t" $10 "1" $1}'
