# read in the sequencing depth for each project
# and compare that with the number of peaks called
# Zijun Zhang
# 5.16.2018

import os
from collections import defaultdict
import pandas as pd

def read_file_peak(peak_fn):
	i = 0
	with open(peak_fn, 'r') as f:
		for line in f:
			i += 1
			#ele = line.strip().split()
			#peak_id = ':'.join([ele[0], ele[1], ele[2], ele[5] ])
			#peak_dict[peak_id][fn_id] = float(ele[6])
	return i

def read_stats_mapping(ms_fn):
	stats = {}
	with open(ms_fn, 'r') as f:
		for line in f:
			ele = line.strip().split('\t')
			if ele[1].endswith('%'):
				stats[ele[0]] = float(ele[1].rstrip('%'))
			else:
				stats[ele[0]] = int(ele[1])
	return stats


def peak_vs_depth():
	par_dir = '/u/flashscratch/flashscratch1/f/frankwoe/m6A_CLAM/Snakemake_m6A_pipeline_2/projects'
	peak_dict = defaultdict(dict)
	GSE_list = [ x for x in os.listdir(par_dir) if x.startswith('GSE') ] 
	for GSE in GSE_list:
		print GSE
		if not os.path.isdir(os.path.join(par_dir, GSE, 'projects', GSE, 'clam')):
			continue
		if not os.path.isdir(os.path.join(par_dir, GSE, 'projects', GSE, 'star')):
			continue
		peak_list = [x for x in os.listdir( os.path.join(par_dir, GSE, 'projects', GSE, 'clam') ) if x.startswith('peaks-') ]
		for peak in peak_list:
			fn_id = peak.split('-')[1]+'.'+GSE
			peak_dict[fn_id]['unique_peak_num'] = read_file_peak(os.path.join(par_dir, GSE, 'projects', GSE, 'clam', peak, 'narrow_peak.unique.bed'))
			try:
				ip_stats = read_stats_mapping(os.path.join(par_dir, GSE, 'projects', GSE, 'star', peak.split('-')[1], 'mapping_stats.txt'))
			except:
				ip_stats = read_stats_mapping(os.path.join(par_dir, GSE, 'projects', GSE, 'star', 'ip', peak.split('-')[1], 'mapping_stats.txt'))
			peak_dict[fn_id]['ip_depth'] = ip_stats['Number of input reads']
			peak_dict[fn_id]['ip_uniq'] = ip_stats['Uniquely mapped reads %']
			peak_dict[fn_id]['ip_multi'] = ip_stats['% of reads mapped to multiple loci']
			try:
				con_stats = read_stats_mapping(os.path.join(par_dir, GSE, 'projects', GSE, 'star', peak.split('-')[2], 'mapping_stats.txt'))
			except:
				con_stats = read_stats_mapping(os.path.join(par_dir, GSE, 'projects', GSE, 'star', 'con', peak.split('-')[2], 'mapping_stats.txt'))
			peak_dict[fn_id]['con_depth'] = con_stats['Number of input reads']
			peak_dict[fn_id]['con_uniq'] = con_stats['Uniquely mapped reads %']
			peak_dict[fn_id]['con_multi'] = con_stats['% of reads mapped to multiple loci']
			try:
				peak_dict[fn_id]['rescued_peak_num'] = read_file_peak(os.path.join(par_dir, GSE, 'projects', GSE, 'clam', peak, 'narrow_peak.rescue.bed'))
			except:
				pass

	peak_df = pd.DataFrame.from_dict(peak_dict, orient='index')
	peak_df.to_csv('peak_vs_depth.csv', sep='\t')