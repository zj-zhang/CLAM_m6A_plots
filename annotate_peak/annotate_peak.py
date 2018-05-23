## annotate all peaks in BED format
## Zijun Zhang
## 5.21.2018

import pyBigWig
from collections import defaultdict


def parse_bigBed((chr, start, end), bb):
	a = bb.entries(chr, start, end)
	if a is None:
		a = []
	return a

def parse_bigWig((chr, start, end), bw):
	return bw.stats(chr, start, end)


def parse_tissue(peak_exp_dict):
	tissue_dict = {}
	with open('replicate_tissue_info.txt', 'r') as f:
		for line in f:
			if line.startswith('#'):
				continue
			ele = line.strip().split('\t')
			tissue_dict[ele[0]] = ele[1::]
	peak_tissue_dict = {}
	for peak in peak_exp_dict:
		peak_tissue_dict[peak] = [tissue_dict[x] for x in peak_exp_dict[peak]]
		
	return peak_tissue_dict


def parse_repeats(peak_exp_dict):
	peak_repeat_dict = defaultdict(list)
	with open('all_peak.with_repeats.bed', 'r') as f:
		for line in f:
			ele = line.strip().split()
			name = ele[3]
			repeat = ele[12]
			if repeat!='.':
				repeat = repeat if ele[5]==ele[14] else 'anti.'+repeat
			peak_repeat_dict[name].append(repeat)
			
	return peak_repeat_dict



def annotate_peak():
	snp_bb = pyBigWig.open('/u/flashscratch/flashscratch1/f/frankwoe/m6A_CLAM/src/SNP/hg19_snp150Common.bb')
	cons_bw = pyBigWig.open('/u/nobackup/yxing/NOBACKUP/frankwoe/hg19/hg19.100way.phyloP100way.bw')
	
	peak_exp_dict = {}
	peak_gene_dict = {}
	peak_cons_dict = {}  # peak -> [phylop_score, snp_num]
	with open('all_peak.sorted.bed', 'r') as f:
		for line in f:
			try:
				ele = line.strip().split()
				gene = ele[-1]
				chrom, start, end = ele[0:3]
				start = int(start)
				end = int(end)
				peak = ele[3]
				exp_list = ele[7].split(',')
				#if len(exp_list)==1:
				#	continue
				peak_exp_dict[peak] = exp_list
				peak_gene_dict[peak] = gene
				
				phylop_score = parse_bigWig((chrom, start, end), cons_bw)[0]
				snp_num = len(parse_bigBed((chrom, start, end), snp_bb))
				
				if phylop_score is None:
					phylop_score = -1000
				
				peak_cons_dict[peak] = [phylop_score, snp_num]
			except:
				pass
	
	peak_tissue_dict = parse_tissue(peak_exp_dict)
	peak_repeat_dict = parse_repeats(peak_exp_dict)
	
	with open('all_peak.annotate.txt', 'w') as fo:
		fo.write('peak\tgene\tcell_line\ttissue\ttreatment\trepeat\tavg_cons\tsnp_num\n')
		for peak in peak_exp_dict:
			if not peak in peak_cons_dict:
				continue
			fo.write('{peak}\t{gene}\t{cell_line}\t{tissue}\t{treatment}\t{repeat}\t{avg_cons}\t{snp_num}\n'.format(
				peak=peak,
				gene=peak_gene_dict[peak],
				cell_line=','.join([x[0] for x in peak_tissue_dict[peak] ]),
				tissue=','.join([x[2] for x in peak_tissue_dict[peak] ]),
				treatment=','.join([x[3] for x in peak_tissue_dict[peak] ]),
				repeat=','.join(peak_repeat_dict[peak]),
				avg_cons=str(peak_cons_dict[peak][0]),
				snp_num=str(peak_cons_dict[peak][1])
				))
		
if __name__ == '__main__':
	annotate_peak()