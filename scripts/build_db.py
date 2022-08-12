#!/usr/bin/env python3
import sys, os, argparse, shutil
import signal, time
import multiprocessing as mp
import extract_kmers 

# build_db.py implements a pipeline for building GT-Pro customized database
# it requires a list of paths to species folders as the main input
# details and guidelines about preparing such a species folder can be found
# in GT-Pro main docs as well as the helper text of this program
def parse_args():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""Usage: GT_Pro build --in <input file> --out <directory of ouput> [--dbname <name of the database>] [--threads <no. of threads>] [--overwrite <overwrite output directory>] [--keep-tmp <do not purge temporary files>]""")

	parser.add_argument('--in', dest='in', type=str, required=True,
		help="""A list of paths, each path leads to a species directory containing the following, names should be matched exactly: 
		a directory hosting genome sequences: genomes/;
		a reference genome: reference.fna;
		a vcf file hosting SNPs in the coordinate of reference genome: core_snps.vcf;
		a multiple sequences alignment of genomic sequences in genomes/ to the reference genome: msa.fa 
		and a coords file (optional) for targeting specific regions on the reference genome: coords.tsv""")
	parser.add_argument('--out', type=str, dest='out', required=True,
		help="""A directory to output files""")
	parser.add_argument('--dbname', type=str, dest='dbname', default="sckmer_db",
		help="""Name for the output sckmer database (default=sckmer_db)""")
	parser.add_argument('--overwrite', dest='overwrite', action='store_true', help="""Overwrite existing output files""")
	parser.add_argument('--keep-tmp', dest='keep_tmp', action='store_true', help="""Do not purging all itermediate output files at the end of database building""")
	parser.add_argument('--threads', dest='n_threads', default=1, type=int,
		help="""No. of threads to use (default=1) for database building""")

	args = vars(parser.parse_args())

	return args

def run_command(cmd, env=None):
	import subprocess as sp
	if env:
		p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, env=env)
	else:
		p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
	out, err = p.communicate()
	if p.returncode != 0:
		err_msg =  "\nError: the following returned non-zero status: '%s':\n" % cmd
		err_msg += "\n%s" % err
		sys.exit(err_msg)
	else:
		return out, err

def parallel(function, argument_list, threads):
	def init_worker():
		signal.signal(signal.SIGINT, signal.SIG_IGN)

	pool = mp.Pool(int(threads), init_worker)

	try:
		results = []
		for arguments in argument_list:
			p = pool.apply_async(function, args=arguments)
			results.append(p)
		pool.close()

		while True:
			if all(r.ready() for r in results):
				return [r.get() for r in results]
			time.sleep(1)

	except KeyboardInterrupt:
		pool.terminate()
		pool.join()
		sys.exit("\nKeyboardInterrupt")

def read_input_list(in_path):
	input_path_array = []
	with open(in_path, 'r') as fh:
		for line in fh:
			input_path_array.append(line.rstrip())
			sys.stderr.write("read {}\n".format(line.rstrip()))

	sys.stderr.write("[OK] a total of {} species found.\n".format(len(input_path_array)))

	return input_path_array

# this function validates each species folder
# it particularly looks for file following files or directory in a species folder
# 1. a directory named "genomes" containing whole genome sequences in FASTA format
# 2. a reference genome named "reference.fna"
# 3. a multiple sequence alignment file named "msa.fa" which hosting the alignment of whole genomes in "genomes" to "reference.fna"
# 4. a VCF file named "core_snps.vcf" describing core snps in the genomic coordinate of "reference.fna" 
# 5. (optional) a tab separated values file named "coords.tsv" 
# It also makes an attempt to assign species lable to input valid species
def validate_input_paths(input_path_array):
	path_objs = []
	dir_target = 'genomes'
	val_targets = ['reference.fna', 'msa.fa', 'core_snps.vcf']
	opt_target = 'coords.tsv'
	for i, inpath in enumerate(input_path_array):
		path_map = dict()
		sys.stderr.write("check {} for required files\n".format(inpath))
		# assign species label to each input species
		# species label must be a six digit string, and the first digit should be from 1-9
		species_label = str(100001 + i)
		path_map['species_lab'] = species_label
		path_map['species_dir'] = inpath.rstrip('/') + '/'
		path_map['genome_dir'] = path_map['species_dir'] + dir_target + '/'
		path_map['genome_paths'] = locate_genomes(path_map['genome_dir'])
		assert len(path_map['genome_paths']) > 10
		
		vpaths = []
		for vtarget in val_targets:
			assert os.path.isfile(path_map['species_dir'] + vtarget)
			sys.stderr.write("[OK] {} found.\n".format(path_map['species_dir'] + vtarget))
			vpaths.append(path_map['species_dir'] + vtarget)
		path_map['vtarget_paths'] = vpaths

		if os.path.isfile(path_map['species_dir'] + opt_target):
			path_map['otarget_paths'] = path_map['species_dir'] + opt_target

		path_objs.append(path_map)
	
	return path_objs 

# locate all genome FASTA files (.fna or .fa) in a directory
def locate_genomes(in_dir):
	assert os.path.isdir(in_dir) 
	fpaths = []
	for f in os.listdir(in_dir):
		fpath = in_dir.rstrip('/')+'/'+f
		if os.path.isfile(fpath) and fpath[-4:] == ".fna" or fpath[-3:] == ".fa":
			fpaths.append(fpath)
			sys.stderr.write("\tgenome path found: {}\n".format(fpath))

	sys.stderr.write("{} genomes sequences will be used for database building\n".format(len(fpaths)))

	return fpaths

# call extract_kmers.py submodule which extracts snp-covering k-mers
def run_extract_kmers(path_objs, output_dir, n_threads=1):
	temp_dir = output_dir.rstrip('/')+'/temp/extract/'
	if not os.path.isdir(temp_dir):
		os.makedirs(temp_dir)

	arg_list = []
	for path_obj in path_objs:
		path_obj['kmer_stage1'] = temp_dir + path_obj['species_lab'] + "-snp-kmer.tsv"
		arg_list.append([path_obj, path_obj['kmer_stage1']])
	parallel(extract_kmers.run, arg_list, n_threads)
	
# call db_val to validate snp-covering k-mers (sck-mers) within species
def intraspec_eval(path_obj, output_dir, n_threads=1):
	assert 'kmer_stage1' in path_obj
	
	temp_dir = output_dir.rstrip('/')+'/temp/inspec_eval/'
	if not os.path.isdir(temp_dir):
		os.makedirs(temp_dir)
	
	genome_list = temp_dir + path_obj['species_lab'] + '.list'
	with open(genome_list, 'w') as fh:
		for genome_path in path_obj['genome_paths']:
			fh.write("{}\n".format(genome_path))

	command = "db_val "
	command += "-d {} ".format(path_obj['kmer_stage1'])
	command += "-n {} ".format(path_obj['species_lab'])
	command += "-t {} ".format(n_threads) 
	command += "-L {} ".format(genome_list)
	command += "1> {}{}_kmer_profiles.tsv ".format(temp_dir, path_obj['species_lab']) 
	command += "2> {}{}_kmer_profiles.log".format(temp_dir, path_obj['species_lab'])
	environ = os.environ.copy()
	run_command(command, environ)

	path_obj['genome_list'] = genome_list
	raw_profile = "{}{}_kmer_profiles.tsv".format(temp_dir, path_obj['species_lab'])
	sckmer_profile = "{}{}.sckmer_profiles.tsv".format(temp_dir, path_obj['species_lab'])
	
	code_kmers = "{}{}.kmers.db.tsv".format(temp_dir, path_obj['species_lab'])
	with open(code_kmers, 'w') as fw:
		with open(sckmer_profile, 'w') as sfw:
			with open(raw_profile, 'r') as fh:
				for line in fh:
					items = line.rstrip().split('\t')
					if int(items[8]) > 0:
						continue
					if float(items[7]) / (float(items[6]) + float(items[7])) < 0.5:
						continue
					if int(items[10]) == 0 or int(items[11]) == 0:
						continue
					sfw.write(line)
					fw.write("{}\t{}0{}\n".format(items[2], items[9], items[0]))
					fw.write("{}\t{}1{}\n".format(items[3], items[9], items[0]))
					fw.write("{}\t{}0{}\n".format(items[4], items[9], items[0]))
					fw.write("{}\t{}1{}\n".format(items[5], items[9], items[0]))
	path_obj['kmer_stage2'] = "{}{}.sckmer_profiles.tsv".format(temp_dir, path_obj['species_lab'])
	path_obj['kmer_stage3'] = code_kmers
	
# extract and pool all possible k-mers from whole genome sequences per species 
# this step call kmc, kmc_dump and mk_pool
def mk_kpool(path_obj, output_dir, n_threads=1):
	assert 'genome_list' in path_obj
	temp_dir = output_dir.rstrip('/')+'/temp/kpools/'
	if not os.path.isdir(temp_dir):
		os.makedirs(temp_dir)

	command = "mkdir {}/{} && ".format(temp_dir, path_obj['species_lab'])
	command += "kmc "
	command += "-k31 -m32 -ci1 "
	command += "-t{} ".format(n_threads)
	command += "-fm @{} ".format(path_obj['genome_list']) 
	command += "{}{}_res ".format(temp_dir, path_obj['species_lab'])
	command += "{}{} ".format(temp_dir, path_obj['species_lab'])
	command += "&> {}{}_res.log".format(temp_dir, path_obj['species_lab'])
	environ = os.environ.copy()
	run_command(command, environ)	

	command = "kmc_dump "
	command += "-ci1 "
	command += "{}{}_res ".format(temp_dir, path_obj['species_lab'])
	command += "{}{}_res.list ".format(temp_dir, path_obj['species_lab'])
	environ = os.environ.copy()
	run_command(command, environ)	

	command = "mk_pool "
	command += "{}{}_res.list ".format(temp_dir, path_obj['species_lab'])
	command += "1> {}{}_pool_wtrc.bin ".format(temp_dir, path_obj['species_lab'])
	command += "2> {}{}_pool_wtrc.log".format(temp_dir, path_obj['species_lab'])
	environ = os.environ.copy()
	run_command(command, environ)	
	path_obj['kmer_pool'] = "{}{}_pool_wtrc.bin".format(temp_dir, path_obj['species_lab'])

# apply n-1 species filter to validate sck-mers
# those who pass the filter are identified as species-specific sck-mers
def ss_screen(path_objs, output_dir, n_threads=1):
	for i, path_obj in enumerate(path_objs):
		assert 'kmer_pool' in path_obj
		kmer_pools = []
		for j, obj in enumerate(path_objs):
			if i == j:
				continue
			else:
				kmer_pools.append(obj['kmer_pool'])
		path_obj['n_minus_1_filter'] = kmer_pools
				
	temp_dir = output_dir.rstrip('/')+'/temp/ss_kmers/'
	if not os.path.isdir(temp_dir):
		os.makedirs(temp_dir)

	for path_obj in path_objs:
		filter_list = temp_dir + path_obj['species_lab'] + '.filter.list'
		with open(filter_list, 'w') as fh:
			for filter_path in path_obj['n_minus_1_filter']:
				fh.write("{}\n".format(filter_path))
		path_obj['n_minus_1_filter_path'] = filter_list
	
	arg_list = []
	for path_obj in path_objs:
		assert 'kmer_stage3' in path_obj 
		ss_bin = "{}{}.ss.hq.bin".format(temp_dir, path_obj['species_lab'])
		sckmer_allowed = "{}{}.sckmer_allowed.tsv".format(temp_dir, path_obj['species_lab'])
		ss_log = "{}{}.ss.hq.log".format(temp_dir, path_obj['species_lab'])
		arg_list.append([path_obj['kmer_stage3'], ss_bin, path_obj['n_minus_1_filter_path'], sckmer_allowed, ss_log])
		path_obj['kmer_stage4'] = "{}{}.sckmer_allowed.tsv".format(temp_dir, path_obj['species_lab'])

	parallel(ss_screen_single, arg_list, n_threads)

def ss_screen_single(kmer_stage3, ss_bin, filter_list, sckmer_allowed, ss_log):
	command = "db_uniq "
	command += "-d {} ".format(kmer_stage3)
	command += "-o {} ".format(ss_bin)
	command += "-L {}".format(filter_list)
	environ = os.environ.copy()
	out, err = run_command(command, environ)

	command = "db_dump "
	command += "{} ".format(ss_bin)
	command += "1> {} ".format(sckmer_allowed)
	command += "2> {}".format(ss_log)
	environ = os.environ.copy()
	run_command(command, environ)

# this step pools all species-specific sck-mers and charaterize snp centered spans (sc-spans) for compression
# it also purges all temporary files unless otherwise specified
def merge_build(path_objs, output_dir, dbname):
	temp_dir = output_dir.rstrip('/')+'/temp/final_build/'
	if not os.path.isdir(temp_dir):
		os.makedirs(temp_dir)

	for path_obj in path_objs:	
		assert 'kmer_stage2' in path_obj
		assert 'kmer_stage4' in path_obj
		
		sckmer_profiles = "{}{}.sckmer_profiles.tsv".format(temp_dir, path_obj['species_lab'])
		sckmer_allowed = "{}{}.sckmer_allowed.tsv".format(temp_dir, path_obj['species_lab'])

		shutil.move(path_obj['kmer_stage2'], sckmer_profiles)
		shutil.move(path_obj['kmer_stage4'], sckmer_allowed)
		
		path_obj['kmer_stage2'] = sckmer_profiles
		path_obj['kmer_stage4'] = sckmer_allowed
	
	command = "db_build "
	command += "`ls {}*allowed*tsv | sort` ".format(temp_dir)
	command += "> {}{}.bin".format(output_dir.rstrip('/')+'/', dbname)
	environ = os.environ.copy()
	run_command(command, environ)


# finally generates two accessory files
# 1. dbname.snp_dict.tsv which can be used for parsing raw GT-Pro genotype output into a more human readable format
# 2. dbname.species_map.tsv which provides a one-to-one mapping relationship between species label and input specie folder
def generate_accessories(path_objs, output_dir, dbname, keep_tmp=False):
	snp_dict = output_dir.rstrip('/')+'/'+dbname+'.snp_dict.tsv'
	spec_lab_map = output_dir.rstrip('/')+'/'+dbname+'.species_map.tsv'

	with open(spec_lab_map, 'w') as fw:
		fw.write('species\tinput_dir\n')
		for path_obj in path_objs:
			assert 'species_lab' in path_obj and 'species_dir' in path_obj
			fw.write('{}\t{}\n'.format(path_obj['species_lab'], path_obj['species_dir']))
	
	with open(snp_dict+"_tmp", 'w') as fw:
		fw.write('#species\tglobal_pos\tcontig\tlocal_pos\tref_allele\talt_allele\n')
		for path_obj in path_objs:
			assert 'species_lab' in path_obj and 'vtarget_paths' in path_obj and 'kmer_stage4' in path_obj
			assert os.path.isfile(path_obj['vtarget_paths'][2])
			assert os.path.isfile(path_obj['kmer_stage4'])

			allowed_kmers = dict()
			with open(path_obj['kmer_stage4'], 'r') as fh:
				for line in fh:
					items = line.rstrip().split('\t')
					gb_pos = items[1][7:]
					allowed_kmers[gb_pos] = 1
			
			with open(path_obj['vtarget_paths'][2], 'r') as fh:
				for line in fh:
					if line[0] == '#':
						continue
					items = line.rstrip().split('\t')
					contig = items[0]
					loc_pos = items[1]
					gb_pos = items[2]
					ref = items[3]
					alt = items[4]

					if gb_pos not in allowed_kmers:
						continue
					fw.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(path_obj['species_lab'], gb_pos, contig, loc_pos, ref, alt))

	command = "(head -n 1 {}_tmp && ".format(snp_dict)
	command += "sed '1d' {}_tmp | ".format(snp_dict)
	command += "sort -k1,1 -k2,2n) "
	command += "> {} && ".format(snp_dict)
	command += "rm {}_tmp".format(snp_dict)
	environ = os.environ.copy()
	run_command(command, environ)

	if keep_tmp is not True:
		shutil.rmtree(output_dir.rstrip('/')+'/temp/')

# main function presents a pipeline which calls the submodules orderly 
# read input -> validate neccessary files -> extract sck-mers ->
# validate sck-mers within species -> validate sck-mers with n-1 species filter ->
# compress species-specific sck-mers and generate sc-spans ->
# generate accessory files
def main():
	args = parse_args()

	input_list = args['in']
	output_dir = args['out']

	overwrite = args['overwrite']
	keep_tmp = args['keep_tmp']
	dbname = args['dbname']
	n_threads = args['n_threads']
	
	if os.path.isdir(output_dir):
		if overwrite is False:
			sys.stderr.write("{} exists; can not use it as output directory\n".format(output_dir))
			sys.stderr.write("consider using --overwrite or choosing another location\n")
			sys.exit(1)
		else:
			sys.stderr.write("{} exists and will be overwritten\n".format(output_dir))
			shutil.rmtree(output_dir)
	os.makedirs(output_dir)

	input_array = read_input_list(input_list)
	path_objs = validate_input_paths(input_array)

	run_extract_kmers(path_objs, output_dir, n_threads)

	for path_obj in path_objs:
		intraspec_eval(path_obj, output_dir, n_threads)
		mk_kpool(path_obj, output_dir, n_threads)
	
	ss_screen(path_objs, output_dir, n_threads)
	merge_build(path_objs, output_dir, dbname)
	generate_accessories(path_objs, output_dir, dbname, keep_tmp)


if __name__ == "__main__":
	main()
