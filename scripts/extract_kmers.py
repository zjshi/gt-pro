import sys, os, argparse, copy, signal
import numpy as np
import multiprocessing as mp

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def parallel(function, argument_list, threads):
	""" Based on: https://gist.github.com/admackin/003dd646e5fadee8b8d6 """

	def init_worker():
		signal.signal(signal.SIGINT, signal.SIG_IGN)

	sys.stderr.write("initializing the process pool\n")
	pool = mp.Pool(int(threads), init_worker)
	sys.stderr.write("Done initializing, the pool now has {} workers\n".format(threads))

	try:
		results = []
		for i, arguments in enumerate(argument_list):
			p = pool.apply_async(function, args=arguments)
			results.append(p)
			sys.stderr.write("worker {} committed to work\n".format(i))

		pool.close()
		sys.stderr.write("multlpleprocessing pool is now closed\n")

		while True:
			if all(r.ready() for r in results):
				return [r.get() for r in results]
			sleep(1)

	except KeyboardInterrupt:
		pool.terminate()
		pool.join()
		sys.exit("\nKeyboardInterrupt")

"""
This function fetch all possible kmers of size specified by kmer_size argument,
all fetched kmers will serve as the database used for verify the uniqueness of kmers that actually snp
build a dictionary for all possible kmers is quite time consuming, one possible way to speed up is to pre-compute it and store it in form of binary file that is directly loadable, thus avoid hash computation again.

Upadate: 07/01/18
Using the following format for the sake of simplicity
kmer_seq: ATGC
"""
def fetch_all_kmers(genome_seq, kmer_size, coords=None):
	kmers = []

	if len(genome_seq) < kmer_size:
		return kmers
	else:
		for i in range(len(genome_seq)-kmer_size+1):
			kmer = genome_seq[i:(i+kmer_size)]
			kmers.append(kemr)
			rc_kmer = revcomp(kmer)
		return kmers

def build_kmer_db(genome_seq, kmer_size):
	kmers = dict()

	if len(genome_seq) < kmer_size:
		return kmers
	else:
		for i in range(len(genome_seq)-kmer_size+1):
			kmer = genome_seq[i:(i+kmer_size)]
			if kmer not in kmers:
				kmers[kmer] = 1
			else:
				kmers[kmer] = kmers[kmer] + 1

		#for kmer, count in kmers.iteritems():
		#	sys.stderr.write("{}\t{}\n".format(kmer, count))

		return kmers

"""
The module takes kmer size argument and allows to build kmer database of different size.
It also takes a kmer_type argument for allowing two approaches to search kmers that cover target snps:
1) fetch all eligible kmers; fetch_all_snp_kmers()
2) fetch kmers whose target snp was at the center; fetch_center_snp_kmers().

Upadate: 07/01/18
Using the following format for the sake of simplicity
id/glob_pos: 11111
allele_pos_on_kmer: 3
kmer_seq(REF/+): ATGC
kmer_seq(ALT/+): ATTC
kmer_seq(REF/-): GCAT
kmer_seq(ALT/-): GAAT
"""
def fetch_snp_kmers(genome_seq, snp_pos, snp_alleles, kmer_size, kmer_type, coords=None):
	if kmer_type == 'all':
		return fetch_all_snp_kmers(genome_seq, snp_pos, snp_alleles, kmer_size, coords)
	elif kmer_type == 'center':
		return fetch_center_snp_kmers(genome_seq, snp_pos, snp_alleles, kmer_size, coords)
	else:
		sys.exit("the specified kmer_type value was not recognized by the program: {}".format(kmer_type))

def fetch_all_snp_kmers(genome_seq, snp_pos, snp_alleles, kmer_size, coords=None):
	sys.stderr.write("[searching] start to search {}-mers\n".format(kmer_size))

	inds_map = None
	if coords is not None:
		inds_map = [None for i in range(len(genome_seq))]

		for i, coord in enumerate(coords):
			for j in range(int(coord[1]), int(coord[2])+1):
				inds_map[j] = i

	kmers = []
	for ri, pos in enumerate(snp_pos):
		kmer_start = int(pos)-kmer_size+1
		kmer_end = int(pos)+kmer_size-1

		if inds_map is not None:
			if inds_map[int(pos)] is None:
				continue

			cur_coord = coords[inds_map[int(pos)]]
			coord_start, coord_end = int(cur_coord[1]), int(cur_coord[2])
			kmer_start, kmer_end = max(coord_start, kmer_start), min(coord_end, kmer_end)

		if kmer_end - kmer_start + 1 >= kmer_size:
			subseq = genome_seq[kmer_start:(kmer_end+1)]

			for i in range(len(subseq)-kmer_size+1):
				kmer = subseq[i:(i+kmer_size)]

				var_pos = kmer_size-i-1
				
				kmer = kmer[:var_pos]+snp_alleles[ri][0]+kmer[var_pos+1:]
				akmer = kmer[:var_pos]+snp_alleles[ri][1]+kmer[var_pos+1:]

				rc_kmer = revcomp(kmer)
				rc_akmer = revcomp(akmer)

				kmers.append([pos, var_pos, kmer, akmer, rc_kmer, rc_akmer])
	sys.stderr.write("	a total of {} kmers was found\n".format(len(kmers)))
	return kmers

def load_msa(msa_path):
	genome_msa = dict()

	with open(msa_path, 'r') as fh:
		for line in fh:
			if line[0] == '>':
				working_id = line.split(' ')[0][1:]
			elif line[0] == '=':
				pass
			else:
				if working_id not in genome_msa:
					genome_msa[working_id] = ""
				
				genome_msa[working_id] = genome_msa[working_id] + line.rstrip()
	
	genome_seqs = [genome_msa[key] for key in genome_msa.keys()]

	return genome_seqs


def fetch_all_from_msa(genome_seqs, ref_seq, snp_pos, snp_alleles, kmer_size, coords=None):
	sys.stderr.write("[searching] start to search {}-mers\n".format(kmer_size))

	inds_map = None
	if coords is not None:
		inds_map = [None for i in range(len(ref_seq))]

		for i, coord in enumerate(coords):
			for j in range(int(coord[1]), int(coord[2])+1):
				inds_map[j] = i

	kmer_records = []
	for ri, pos in enumerate(snp_pos):
		kmer_start = int(pos)-kmer_size+1
		kmer_end = int(pos)+kmer_size-1

		if inds_map is not None:
			if inds_map[int(pos)] is None:
				continue

			cur_coord = coords[inds_map[int(pos)]]
			coord_start, coord_end = int(cur_coord[1]), int(cur_coord[2])
			kmer_start, kmer_end = max(coord_start, kmer_start), min(coord_end, kmer_end)

		if kmer_end - kmer_start + 1 >= kmer_size:
			subseqs = [genome_seq[kmer_start:(kmer_end+1)] for genome_seq in genome_seqs]

			for i in range(len(subseqs[0])-kmer_size+1):
				raw_kmers = [subseq[i:(i+kmer_size)] for subseq in subseqs]

				kmers = []
				for rk in raw_kmers:
					if '-' not in rk and 'N' not in rk:
						kmers.append(rk)

				ukmers, counts = np.unique(kmers, return_counts=True)
				uk_inds = np.argsort(counts)[::-1]

				var_pos = kmer_size-i-1

				kmer = ""
				akmer = ""
				kflag = False
				akflag = False
				for ukmer in ukmers[uk_inds]:
					if kflag is False:
						if ukmer[var_pos] == snp_alleles[ri][0]:	
							kmer = ukmer
							kflag = True
					
					if akflag is False:
						if ukmer[var_pos] == snp_alleles[ri][1]:
							akmer = ukmer
							akflag = True

					if kflag is True and akflag is True:
						break
				
				if len(kmer) != 31 or len(akmer) != 31:
					continue

				rc_kmer = revcomp(kmer)
				rc_akmer = revcomp(akmer)

				kmer_records.append([pos, var_pos, kmer, akmer, rc_kmer, rc_akmer])

	sys.stderr.write("	a total of {} kmer records was found\n".format(len(kmer_records)))
	return kmer_records

def fetch_center_snp_kmers(genome_seq, snp_pos, snp_alleles, kmer_size, coords=None):
	sys.stderr.write("[searching] start to search {}-mers\n".format(kmer_size))

	inds_map = None
	if coords is not None:
		inds_map = [None for i in range(len(genome_seq))]

		for i, coord in enumerate(coords):
			for j in range(int(coord[1]), int(coord[2])+1):
				inds_map[j] = i

	is_even = (kmer_size % 2 == 0)

	kmers = []
	for ri, pos in enumerate(snp_pos):
		kmer_start, kmer_end, var_pos = 0, 0, 0

		if is_even:
			var_pos = int(kmer_size/2)
			kmer_start = int(pos)-int(kmer_size/2)+1
			kmer_end = int(pos)+int(kmer_size/2)
		else:
			var_pos = int(kmer_size/2)+1
			kmer_start = int(pos)-int(kmer_size/2)
			kmer_end = int(pos)+int(kmer_size/2)

		if inds_map is not None:
			if inds_map[int(pos)] is None:
				continue

			cur_coord = coords[inds_map[int(pos)]]

			coord_start = int(cur_coord[1])
			coord_end = int(cur_coord[2])

			if kmer_start < coord_start or kmer_end > coord_end:
				continue

		kmer = genome_seq[kmer_start:(kmer_end+1)]
		akmer = kmer[:var_pos]+snp_alleles[ri][0]+kmer[var_pos+1:]

		akmer = kmer
		akmer = kmer[:var_pos]+snp_alleles[ri][1]+kmer[var_pos+1:]

		rc_kmer = revcomp(kmer)
		rc_akmer = revcomp(akmer)

		kmers.append([pos, var_pos, kmer, akmer, rc_kmer, rc_akmer])
	sys.stderr.write("	a total of {} kmers was found\n".format(len(kmers)))
	return kmers

def revcomp(seq):
	""" Reverse complement sequence

	Args:
		seq:	string from alphabet {A,T,C,G,N}

	Returns:
		reverse complement of seq
	"""
	complement = {
		'A':'T',
		'T':'A',
		'G':'C',
		'C':'G',
		'N':'N',
		'R':'N',
		'Y':'N',
		'K':'N',
		'M':'N',
		'S':'N',
		'W':'N',
		'B':'N',
		'D':'N',
		'H':'N',
		'V':'N'
	}
	return ''.join([complement[_] for _ in seq[::-1]])

def dump_tsv(kmers, output):
	snp_kmer_path = output.rstrip('/') + "-snp-kmer.tsv"
	with open(snp_kmer_path, 'w') as fh:
		for kmer in kmers:
			if len(kmer) == 6:
				fh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(*kmer))
			elif len(kmer) == 9:
				fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(*kmer))
			else:
				assert False

# a mini function to load all coordinates in to memory
def read_coords(fpath):
	sys.stderr.write("[load] loading key coordinates on core-genome from {}\n".format(fpath))

	coords = []
	with open(fpath, "r") as fh:
		fh.readline()
		for line in fh:
			coords.append(line.rstrip('\n').split('\t'))

	sys.stderr.write("	a total of {} divisions was found\n".format(str(len(coords))))
	return coords

def open_vcf_file(fpath):
	"""
	* ``Record.CHROM``; string
	* ``Record.POS``; int
	* ``Record.ID``; None
	* ``Record.REF``; string
	* ``Record.ALT``; list
	* ``Record.QUAL``; None
	* ``Record.FILTER``; list
	* ``Record.INFO``; dictionary

	additional attributes:
	* ``Record.FORMAT``; string
	* ``Record.samples``; list
	* ``Record.genotype``; object
	"""

	sys.stderr.write("[load] loading core snps from {}\n".format(fpath))
	snp_gb_pos, snp_alleles = [], []

	with open(fpath, 'r') as fh:
		for line in fh:
			if line[0] == '#':
				continue
			items = line.rstrip().split('\t')[0:5]
			if len(items[4]) > 1:
				continue
			else:
				snp_gb_pos.append(int(items[2]))
				snp_alleles.append([items[3], items[4]])

	sys.stderr.write("	a total of {} core bi-allelic snps was found\n".format(len(snp_gb_pos)))

	return snp_gb_pos, snp_alleles 

def open_genome_seq(genome_path):
	sys.stderr.write("[load] loading core-genome consensus sequence from {}\n".format(genome_path))

	records = list(SeqIO.parse(genome_path, "fasta"))
	main_genome = ""
	for record in records:
		main_genome = main_genome + str(record.seq).upper()

	sys.stderr.write("	the loaded core-genome has a consensus sequence of {} bases\n".format(str(len(main_genome))))

	return main_genome

def run(path_obj, output_dir):
	spec_lab = path_obj['species_lab']
	ref_path = path_obj['vtarget_paths'][0] 
	msa_path = path_obj['vtarget_paths'][1] 
	vcf_path = path_obj['vtarget_paths'][2]

	coords_path = None
	if 'otarget_paths' in path_obj:
		coords_path = path_obj['otarget_paths']

	k_size = 31
	k_type = 'all'
	
	ref_seq = open_genome_seq(ref_path)
	snp_gb_pos, snp_alleles = open_vcf_file(vcf_path)

	coords = None
	if coords_path is not None:
		coords = read_coords(coords_path)

	genome_seqs = load_msa(msa_path)
	snp_kmers = fetch_all_from_msa(genome_seqs, ref_seq, snp_gb_pos, snp_alleles, k_size, coords)

	outname = output_dir + spec_lab 
	dump_tsv(snp_kmers, outname)
	path_obj['kmer_stage1'] = outname + "-snp-kmer.tsv"
