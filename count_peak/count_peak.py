## count peak intervals across different datasets
## Zijun Zhang
## 4.4.2018
## revised 5.23.2018: add gene name

from collections import defaultdict
import os
import math

def read_file_peak(peak_fn, fn_id, peak_dict, peak_to_gene):
	with open(peak_fn, 'r') as f:
		for line in f:
			ele = line.strip().split()
			peak_id = ':'.join([ele[0], ele[1], ele[2], ele[5] ])
			peak_dict[peak_id][fn_id] = float(ele[6])
			peak_to_gene[peak_id] = ele[3].split('.')[0]
	return peak_dict, peak_to_gene


def count_peak():
	par_dir = '/u/flashscratch/flashscratch1/f/frankwoe/m6A_CLAM/Snakemake_m6A_pipeline_2/projects'
	peak_dict = defaultdict(dict)
	peak_to_gene = {}
	GSE_list = [ x for x in os.listdir(par_dir) if x.startswith('GSE') ] 
	for GSE in GSE_list:
		try:
			print GSE
			peak_list = [x for x in os.listdir( os.path.join(par_dir, GSE, 'projects', GSE, 'clam') ) if x.startswith('peaks-') ]
			for peak in peak_list:
					fn_id = peak.split('-')[1]+'.'+GSE
					peak_fn = os.path.join(par_dir, GSE, 'projects', GSE, 'clam', peak, 'narrow_peak.combined.bed')
					peak_dict, peak_to_gene = read_file_peak(peak_fn, fn_id, peak_dict, peak_to_gene)
		except:
			pass
	return peak_dict, peak_to_gene

def write_to_BED():
	peak_dict, peak_to_gene = count_peak()
	BED_formatter = "{chrom}\t{start}\t{end}\t{name}\t{avg_score}\t{strand}\t{exp}\t{rep}\t{score}\t{gene}\n"
	with open('all_peak.bed', 'w') as f:
		for peak in peak_dict:
			chrom, start, end, strand = peak.split(':')
			avg_score = int(math.exp( sum(peak_dict[peak].values())/len(peak_dict[peak]) ) )
			avg_score = 1000 if avg_score >= 1000 else avg_score
			rep = peak_dict[peak].keys()
			exp = ','.join(list(set([ x.split('.')[1] for x in rep ])))
			score = ','.join([ str(peak_dict[peak][x]) for x in rep ])
			f.write(BED_formatter.format( chrom=chrom, start=start, end=end, name=peak,
				strand=strand, avg_score=avg_score, exp=exp, rep=','.join(rep), score=score,
				gene=peak_to_gene[peak]))
	
	os.system('sort -k1,1 -k2,2n all_peak.bed > all_peak.sorted.bed')
	#os.system('bedToBigBed -type=bed6+4 all_peak.sorted.bed ~/nobackup/programs/UCSC/hg19.chrom.sizes all_peak.sorted.bb')
	os.remove('all_peak.bed')
	#os.remove('all_peak.sorted.bed')


if __name__ == '__main__':
	write_to_BED()