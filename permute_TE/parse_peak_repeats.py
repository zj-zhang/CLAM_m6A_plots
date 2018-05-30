'''
parse peak repetitive elements
Zijun Zhang
5.28.2018
'''

import pybedtools
import pandas as pd
from collections import defaultdict

def	intersect_with_repMask(fn, repMasker='/u/flashscratch/flashscratch1/f/frankwoe/hg19_gtf_bed/hg19_repeatMasker.sorted.bed'):
	out = pybedtools.BedTool(fn).intersect(repMasker, wao=True, f=0.5)
	repeat_dict = defaultdict(list)
	with open(out.fn, 'r') as fi:
		for line in fi:
			ele = line.strip().split('\t')
			if ele[10]!='.':
				repeat = ele[10] if ele[5]==ele[12] else 'anti.'+ele[10]
			else:
				repeat = '.'
			repeat_dict[ele[3]].append(repeat)
	return repeat_dict

def parser(fn, repMasker='/u/flashscratch/flashscratch1/f/frankwoe/hg19_gtf_bed/hg19_repeatMasker.sorted.bed'):
	return intersect_with_repMask(fn, repMasker)