#!/usr/bin/env python3
import sys, os, argparse, shutil
import extract_kmers 

def parse_args():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="""Usage: optimize --db <GT-Pro database name> --in <Fastq file path> --out <directory of ouput>""")

	parser.add_argument('-d', '--db', dest='db', type=str, required=True,
		help="""GT-Pro database to optimize""")
	parser.add_argument('-i', '--in', dest='in', type=str, required=True,
		help="""A path to a test fastq data. Consider using the one (SRR413665_2.fastq.gz) included in the test directory for convenience""") 
	parser.add_argument('--tmp', type=str, dest='tmp', default="./",
		help="""A directory to host temporary files""")

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

def optimize(dbname):
	sys.stderr.write("[OK] start initial optimization\n")
	raw_bin = dbname + ".bin"
	kmer_index_bin = dbname + "_optimized_db_kmer_index.bin"
	db_snps_bin = dbname + "_optimized_db_snps.bin"
	if os.path.isfile(kmer_index_bin) and os.path.isfile(db_snps_bin):
		sys.stderr.write("[OK] database found\n")
		sys.stderr.write("[OK] optimize from {} and {}\n".format(kmer_index_bin, db_snps_bin))
		command = "gt_pro "
		command += "-d {} ".format(dbname)
		command += "</dev/null"
		environ = os.environ.copy()
		run_command(command, environ)
	elif os.path.isfile(raw_bin):
		sys.stderr.write("[OK] database found\n")
		sys.stderr.write("[OK] optimize from {}\n".format(raw_bin))
		command = "gt_pro "
		command += "-d {} ".format(raw_bin)
		command += "</dev/null"
		environ = os.environ.copy()
		run_command(command, environ)
	else:
		assert False 
	sys.stderr.write("[OK] initial optimization done\n")

def breakin_test(dbname, inpath, tmp_dir):
	sys.stderr.write("[OK] finalize optimization with a break-in test\n")
	assert os.path.isfile(inpath)
	
	tmp_dir = tmp_dir.rstrip('/')+'/'
	if not os.path.isdir(tmp_dir):
		os.makedirs(tmp_dir)

	command = "gt_pro "
	command += "-d {} ".format(dbname)
	command += "{} ".format(inpath)
	command += "-o {} ".format(tmp_dir + '/breakin')
	environ = os.environ.copy()
	run_command(command, environ)

	shutil.rmtree(tmp_dir)

	sys.stderr.write("[OK] optimization done\n")

	
def main():
	args = parse_args()

	dbname = args['db']
	inpath = args['in']
	tmp_dir = args['tmp']

	optimize(dbname)
	breakin_test(dbname, inpath, tmp_dir)


if __name__ == "__main__":
	main()
