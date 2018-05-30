'''
Perform permutation test to find enriched Transposable Elements
Zijun Zhang
5.28.2018
Workflow:
a) compile a list of all 100bp bins in all genes, and their corresponding transcript location and repeats; need location->bin
b) shuffle the target set of peaks by randomly choosing another bin with the same location
c) compute the permutation p-value
'''

import os
import sys
import tempfile
import pybedtools
import re
import pandas as pd
import numpy as np
from collections import Counter, defaultdict
from multiprocessing import Pool

import parse_peak_location as pl
import parse_peak_repeats as pr

def compile_background_bins(gtf_fn='/u/nobackup/yxing/NOBACKUP/frankwoe/hg19/gencode.v19.annotation.gtf', binsize=100):
	BED_formatter = '{chrom}\t{start}\t{end}\t{peak_id}\t.\t{strand}\t{gene}\n'
	fo = tempfile.NamedTemporaryFile()
	with open(gtf_fn, 'r') as fi:
		for line in fi:
			if line.startswith('#'):
				continue
			ele = line.strip().split('\t')
			if ele[2]!="gene":
				continue
			chrom = ele[0]
			start = int(ele[3])
			end = int(ele[4])
			strand = ele[6]
			gene = re.findall(r"(\w+)", ele[-1])[1].split('.')[0]
			for i in range(start, end, binsize):
				peak_id = ':'.join([chrom, str(i), str(i+binsize), strand])
				fo.write(BED_formatter.format(
					chrom=chrom, 
					start=i,
					end=i+binsize,
					peak_id=peak_id,
					gene=gene,
					strand=strand))
	fo.flush()
	return fo


def get_background_dataframe():
	tmp_fo = compile_background_bins()
	bg_loc_dict = pl.parser(tmp_fo.name)
	bg_repeat_dict = pr.parser(tmp_fo.name)
	return (bg_loc_dict, bg_repeat_dict)


def read_allpeak(peak_fn ='/u/flashscratch/flashscratch1/f/frankwoe/m6A_CLAM/src/annotate_peaks/annot_peak.df.gz'):
	return pd.read_table(peak_fn)


def get_single_permutation(bg_loc_dict, bg_repeat_dict, feature_counter):
	bg_repeat_counter = Counter()
	#bg_repeat_dict = {}
	for feature in feature_counter:
		bg_loc = np.random.choice(bg_loc_dict[feature], feature_counter[feature], replace=False)
		bg_repeat_counter += Counter([x for peak in bg_loc for x in bg_repeat_dict[peak] if x!="."])
		#bg_repeat_dict[feature] = Counter([x for peak in bg_loc for x in bg_repeat_dict[peak] if x!="."])
	return bg_repeat_counter

def permute_TE():
	
	# read in files
	print("read in peak file")
	peak_fn ='/u/flashscratch/flashscratch1/f/frankwoe/m6A_CLAM/src/annotate_peaks/annot_peak.df.gz'
	peak_df = read_allpeak(peak_fn)
	usable_peak_df = peak_df.iloc[[i for i in range(peak_df.shape[0]) 
		if peak_df.occurence_num[i]>1 and 
		peak_df.intron[i]!=1 and 
		peak_df.CDS_only[i]+peak_df.UTR_only[i]+peak_df.CDS_and_UTR[i]>0 ] ]
	print("prepare background")
	bg_loc_dict, bg_repeat_dict = get_background_dataframe()
	
	# stats of the observed pattern
	feature_counter = Counter(
		{x: sum(usable_peak_df[x]) for x in ['CDS_only', 'UTR_only', 'CDS_and_UTR'] }
		)

	repeat_counter = Counter([x for x in usable_peak_df['repeat'] if x!="."])
	
	# basically, fix feature_counter distribution, compute p-value for repeat_counter
	# using `get_single_permutation` to derive the background
	B = 9999
	#pool = Pool()
	#bg = pool.map(get_single_permutation, [(bg_loc_dict, bg_repeat_dict, feature_counter, b) for b in ])
	
	bg_dist = defaultdict( lambda : 1 )
	for b in range(B):
		if not b%50:
			print("permuted %i"%b)
		bg_repeat_counter = get_single_permutation(bg_loc_dict, bg_repeat_dict, feature_counter)
		for repeat in bg_repeat_counter:
			if bg_repeat_counter[repeat] >= repeat_counter[repeat]:
				bg_dist[repeat] += 1
				
	
	# write out results
	with open('permutation.out.txt', 'w') as fo:
		for repeat in bg_dist:
			pval = bg_dist[repeat]/(B+1.)
			obs = repeat_counter[repeat]
			if obs>10:
				fo.write("{repeat}\t{obs}\t{pval}\n".format(repeat=repeat, obs=obs, pval=pval))
	
	return 0

if __name__ == '__main__':
	permute_TE()