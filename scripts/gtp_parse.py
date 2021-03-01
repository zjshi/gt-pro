# This script, gtp_parse.py, is a parser written in Python which converts the raw output of GT-Pro to a more human readable style.
# The raw output of GT-Pro is text-based and minimalized, and ready to link to external files, e.g. a SNP dictionary.
# Ths script essentially joins (left) a raw output file of GT-Pro and a complete set of SNPs with site-specific discription. 
# Running this script is not mandatory but can provide an example on how a typical genotype profile looks like.
# If everything runs to finishing without error, parsed GT-Pro will be printed to standard output on default.

import os, sys, argparse

def parse_args():
    """ Return dictionary of command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        usage=argparse.SUPPRESS)
	# An example of snp_dict can be found at https://github.com/zjshi/gt-pro
    parser.add_argument('--dict', type=str, dest='snp_dict', required=True,
        help="""Path to disctionary of SNPs.""")
	# The main input file should be the raw output fron GT-Pro.
    parser.add_argument('--in', type=str, dest='in_path', required=True,
        help="""Path to input file, should be the direct output from gt-pro or gt-pro2.0""")
    parser.add_argument('--out', type=str, dest='out', default='/dev/stdout',
        help="""Path to output file""")

    return vars(parser.parse_args())

# Read genotypes from GT-Pro raw output and sort them by species ID, SNP ID and allele type.
# The raw output is text-based and coded to include the minimal information required for linking external files, e.g. a SNP dictionary.
# The output has two fields: genotype ID and genotype count.
# The genotype ID is a string no longer than 14 chars and consisting of digits only.
# First six digits represents a species ID.
# The seventh digit is either 0 or 1, which identifies the allele type.
# From the eighth digit to the end represents a SNP site ID.
def read_gtp_out(fpath):
	gtp_array = []
	with open(fpath, 'r') as fh:
		for line in fh:
			items = line.rstrip().split()
			if line[0] == '#':
				pass
			elif items[0] == 'snp_id':
				pass
			else:
				sp_id = int(items[0][0:6])
				snp_type = int(items[0][6])
				# SNP ID is the same as the global position of the SNP, 
				# together with species ID, can uniquely identifies a SNP site across hundreds species.  
				pos = int(items[0][7:])
				gtp_array.append([sp_id, pos, snp_type, items[1]])

	
	return sorted(gtp_array, key= lambda x: (x[0], x[1], x[2]))

# Join (left) detected genotypes and a SNP dictionary by species and SNP ID, 
# then organize genotype by SNP sites.
# The SNP dictionary is a regular table,
# which contains milliions of rows and six fields per row. Each row represents a SNP site.
# The fields are species ID, SNP ID, contig ID, local position, reference allele and alternative allele.
# The SNP dictionary is sorted by species and SNP ID numerically upon the retrieval.
# The parsed GT-Pro output file will have the same six decriptive fields in the SNP dictionary 
# and two additional fields: reference allele count and alternative allele count 
def snp_dict_lookup(fpath, gtp_array, out='/dev/stdout'):
	snp_dict = []
	parsed_results = []

	gtp_i = 0

	with open(fpath, 'r') as fh:
		for line in fh:
			if line[0] == '#':
				continue
				
			items = line.rstrip().split()
			sp_id = int(items[0])
			pos = int(items[1])
			
			if sp_id < gtp_array[gtp_i][0]:
				continue
			
			if sp_id == gtp_array[gtp_i][0] and pos < gtp_array[gtp_i][1]:
				continue

			while sp_id > gtp_array[gtp_i][0] or (sp_id == gtp_array[gtp_i][0] and pos > gtp_array[gtp_i][1]):
				gtp_i = gtp_i + 1	
				if gtp_i >= len(gtp_array):
					break

			if gtp_i >= len(gtp_array):
				break

			if sp_id == gtp_array[gtp_i][0] and pos == gtp_array[gtp_i][1]:
				# merge two genotypes by each bi-allelic site
				ct_pair = ["-1", "-1"]

				oppo_type = 1 - gtp_array[gtp_i][2]

				ct_pair[gtp_array[gtp_i][2]] = gtp_array[gtp_i][3]

				if gtp_i == len(gtp_array) - 1:
					ct_pair[oppo_type] = "0"
				else:
					if gtp_array[gtp_i][0] == gtp_array[gtp_i+1][0] and gtp_array[gtp_i][1] == gtp_array[gtp_i+1][1]:
						ct_pair[oppo_type] = gtp_array[gtp_i+1][3]
					else:
						ct_pair[oppo_type] = "0"

				parsed_results.append(items + ct_pair)

	
	with open(out, 'w') as fw:
		fw.write('species\tglobal_pos\tcontig\tlocal_pos\tref_allele\talt_allele\tref_count\talt_count\n')
		for result in parsed_results:
			fw.write("{}\n".format("\t".join(result)))

def main():
	args = parse_args()	

	gtp_fpath = args['in_path']
	snp_dict = args['snp_dict']	

	gtp_array = read_gtp_out(gtp_fpath)
	snp_dict_lookup(snp_dict, gtp_array, args['out'])

if __name__ == "__main__":
	main()
		
