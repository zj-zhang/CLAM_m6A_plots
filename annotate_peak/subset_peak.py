## subset peaks based on their features
## e.g. intronic vs exonic vs UTR peaks
## repeat-derived peaks vs non-repeat peaks
## peaks with motif vs peaks w/o motif
## peaks on different cell lines vs. tissue-specific peaks
## conserved peaks vs SNP density
## Zijun Zhang
## 5.21.2018

import pandas as pd
import numpy as np
import maptlotlib.pyplot as plt

def read_annot():
	# read data and do some clean-up
	df_annot = pd.read_table('all_peak.annotate.txt.gz')
	df_annot.loc[df_annot['avg_cons']==-1000, 'avg_cons'] = np.nan
	df_motif = pd.read_table('../count_peak_and_motif/peak_sum.tsv.gz')
	motif = np.asarray([[int(i) for i in x.split(',')] for x in df_motif['GGACT,RRACT,DRACH']])
	df_motif['motif.GGACT'] = pd.Series(motif[:,0], index=df_motif.index)
	df_motif['motif.RRACT'] = pd.Series(motif[:,1], index=df_motif.index)
	df_motif['motif.DRACH'] = pd.Series(motif[:,2], index=df_motif.index)
	
	# merge
	df_loc = pd.read_table('peak_localization_df.csv.gz')
	df_loc.columns = ['peak',  'CDS_only',  'UTR_only',  'CDS_and_UTR',  'intron']
	tmp = pd.merge(df_annot, df_motif, on='peak')
	df = pd.merge(tmp, df_loc, on='peak')
	
	# plot
	df.boxplot(column=['avg_cons', 'snp_num'], by='CDS_and_UTR')
	plt.show()
	
	# save
	df.to_csv("annot_peak.df.gz", sep='\t', compression='gzip', index=False)
	
	return df

