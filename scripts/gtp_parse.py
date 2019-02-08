import os, sys, argparse

def parse_args():
    """ Return dictionary of command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        usage=argparse.SUPPRESS)
    parser.add_argument('--dict', type=str, dest='snp_dict', required=True,
        help="""Path to disctionary of SNPs""")
    parser.add_argument('--in', type=str, dest='in_path', required=True,
        help="""Path to input file, should be the direct output from gt-pro or gt-pro2.0""")
    parser.add_argument('--v2', dest='ver_2', action='store_true', 
		help="""Toggle for gt-pro2.0 style input""")
    parser.add_argument('--out', type=str, dest='out', default='/dev/stdout',
        help="""Path to output file""")

    return vars(parser.parse_args())

def read_gtp_out(fpath, v2):
	gtp_array = []
	with open(fpath, 'r') as fh:
		for line in fh:
			items = line.rstrip().split()
			if line[0] == '#':
				pass
			elif items[0] == 'snp_id':
				pass
			else:
				if v2 is False:
					sp_id = int(items[0][0:6])
					pos = int(items[0][6:])
					snp_type = int(items[1])
					gtp_array.append([sp_id, pos, snp_type, items[2]])
				else:
					sp_id = int(items[0][0:6])
					snp_type = int(items[0][6])
					pos = int(items[0][7:])
					gtp_array.append([sp_id, pos, snp_type, items[1]])

	
	return sorted(gtp_array, key= lambda x: (x[0], x[1], x[2]))

def snp_dict_lookup(fpath, gtp_array, v2, out='/dev/stdout'):
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
				ct_pair = ["-1", "-1"]

				if v2 is False:
					ct_pair[gtp_array[gtp_i][2]] = gtp_array[gtp_i][3]
					ct_pair[gtp_array[gtp_i+1][2]] = gtp_array[gtp_i+1][3]
				else:
					oppo_type = 1 - gtp_array[gtp_i][2]

					ct_pair[gtp_array[gtp_i][2]] = gtp_array[gtp_i][3]
					ct_pair[oppo_type] = "0"

				parsed_results.append(items + ct_pair)

	
	with open(out, 'w') as fw:
		for result in parsed_results:
			fw.write("{}\n".format("\t".join(result)))

def main():
	args = parse_args()	

	gtp_fpath = args['in_path']
	snp_dict = args['snp_dict']	

	gtp_array = read_gtp_out(gtp_fpath, args['ver_2'])
		
	snp_dict_lookup(snp_dict, gtp_array, args['ver_2'], args['out'])

if __name__ == "__main__":
	main()
		
