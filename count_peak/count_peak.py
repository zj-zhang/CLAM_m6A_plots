## count peak intervals across different datasets
## Zijun Zhang
## 4.4.2018

from collections import defaultdict
import os


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
			