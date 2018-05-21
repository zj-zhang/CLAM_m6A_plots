## count peak intervals across different datasets
## Zijun Zhang
## 4.4.2018

from collections import defaultdict
import os
from pygr import seqdb
import gzip
import re


def read_file_peak(peak_fn, fn_id, peak_dict):
	with open(peak_fn, 'r') as f:
		for line in f:
			ele = line.strip().split()
			peak_id = ':'.join([ele[0], ele[1], ele[2], ele[5] ])
			peak_dict[peak_id][fn_id] = float(ele[6])
	return peak_dict


def count_peak():
	par_dir = '/u/flashscratch/flashscratch1/f/frankwoe/m6A_CLAM/Snakemake_m6A_pipeline_2/projects'
	peak_dict = defaultdict(dict)
	GSE_list = [ x for x in os.listdir(par_dir) if x.startswith('GSE') ] 
	for GSE in GSE_list:
		try:
			print GSE
			peak_list = [x for x in os.listdir( os.path.join(par_dir, GSE, 'projects', GSE, 'clam') ) if x.startswith('peaks-') ]
			for peak in peak_list:
					fn_id = peak.split('-')[1]+'.'+GSE
					peak_fn = os.path.join(par_dir, GSE, 'projects', GSE, 'clam', peak, 'narrow_peak.combined.bed')
					peak_dict = read_file_peak(peak_fn, fn_id, peak_dict)
		except:
			pass
	return peak_dict


def _fetch_seq(genome, chr, start, end, strand):
		"""Fetch a genmoic sequence upon calling `pygr`
		`genome` is initialized and stored in self.genome
		Args:
			chr (str): chromosome
			start (int): start locus
			end (int): end locuse
			strand (str): either '+' or '-'
		Returns:
			seq (str): genomic sequence in upper cases
		"""
		start = int(start)
		end = int(end)
		try:
			seq = genome[chr][start:end]
			if strand == "-":
				seq = -seq
		except:
			raise Exception('pygr cannot fetch sequences')
		return str(seq).upper()


def _match_motif(seq, motif_list):
	match_list = []
	for motif in motif_list:
		if not motif.startswith("("): motif = "(" + motif
		if not motif.endswith(")"): motif = motif + ")"
		match_list.append( str( len( re.findall(motif, seq) ) ) )
	return match_list


def parse_motif():
	peak_dict = count_peak()
	genome = seqdb.SequenceFileDB('/u/home/f/frankwoe/nobackup/hg19/hg19.noRand.fa')
	
	# motifs include: GGACT, RRACT, DRACH
	motif_list = ["(GGACT)", "([GA][GA]ACT)", "([AGT][AG]AC[ACT])"]
	
	with gzip.GzipFile('peak_sum.tsv.gz', 'wb') as f:
		# header line
		f.write("\t".join( [ 'peak', 'exp_num', 'exp', 'occurence_num', 'occurence', 'fc', 'GGACT,RRACT,DRACH', 'seq' ] ) + '\n')
		for peak in peak_dict:
			chrom, start, end, strand = peak.split(':')
			peak_seq = _fetch_seq(genome, chrom, start, end, strand)
			occurence = peak_dict[peak].keys()
			occured_exp = list(set([x.split('.')[1] for x in occurence ]))
			fc = [ str(peak_dict[peak][x]) for x in occurence ]
			motif_count = _match_motif(peak_seq, motif_list)
			
			f.write('\t'.join( [
				peak, 
				str(len(occured_exp)),
				','.join(occured_exp), 
				str(len(occurence)),
				','.join(occurence), 
				','.join(fc), 
				','.join(motif_count), 
				peak_seq ]) + '\n' )
