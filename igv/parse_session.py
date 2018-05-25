## parse an IGV session for CLAM_m6A Snakemake pipeline
## Zijun Zhang
## 5.24.2018

import sys
import os
from collections import defaultdict

path_list, out_fp, locus = sys.argv[1], sys.argv[3], sys.argv[2]
path_list = path_list.split(',')

locus_parts = locus.split(':')
if len(locus_parts)>2:
	locus = locus_parts[0] + ':' + locus_parts[1] + '-' +locus_parts[2]

gse_to_rep = defaultdict(list)
for x in path_list:
	gse_to_rep[x.split('.')[1]].append(x.split('.')[0])

s = ''
s += '''<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'''
s += '''<Session genome="hg19" hasGeneTrack="true" hasSequenceTrack="true" locus="{locus}" path="{fp}" version="8">\n'''.format(
	locus=locus,
	fp=out_fp
	)

avail_path = {}
par_dir = '/u/home/f/frankwoe/scratch/m6A_CLAM/Snakemake_m6A_pipeline_2/projects/'
s += "\t<Resources>\n"
s += '''\t\t<Resource path="/u/home/f/frankwoe/scratch/hg19_gtf_bed/hg19_repeatMasker.sorted.bed"/>\n'''
s += '''\t\t<Resource path="/u/flashscratch/flashscratch1/f/frankwoe/m6A_CLAM/src/igv/hg19/gencode.v19.annotation.filtered.gtf"/>\n'''

for gse in gse_to_rep:
	bw_dir = os.path.join(par_dir, gse, 'projects', gse, 'bigwig')
	if not os.path.isdir(bw_dir):
		continue
	ip_to_full = {x.split('-')[0]:x for x in os.listdir(bw_dir)}
	for name in gse_to_rep[gse]:
		full_name = ip_to_full[name]
		this_ip_combined_fp = os.path.join(par_dir, gse, 'projects', gse, 'bigwig', full_name, full_name.split('-')[0], "combined_pos.bw")
		this_ip_unique_fp = os.path.join(par_dir, gse, 'projects', gse, 'bigwig', full_name, full_name.split('-')[0], "unique_pos.bw")
		this_inp_combined_fp = os.path.join(par_dir, gse, 'projects', gse, 'bigwig', full_name, full_name.split('-')[1], "combined_pos.bw")
		if not (os.path.isfile(this_ip_combined_fp) and \
			os.path.isfile(this_ip_unique_fp) and \
			os.path.isfile(this_inp_combined_fp)):
			continue

		s += '''\t\t<Resource path="{path}"/>\n'''.format(
			path=this_ip_combined_fp
			)
		s += '''\t\t<Resource path="{path}"/>\n'''.format(
			path=this_ip_unique_fp
			)
		s += '''\t\t<Resource path="{path}"/>\n'''.format(
			path=this_inp_combined_fp
			)
		avail_path[name] = [this_ip_combined_fp, this_inp_combined_fp, this_ip_unique_fp]


s+= '\t</Resources>\n'
s+= '\t<Panel height="400" name="DataPanel" width="1200">\n'

for name in avail_path:
	# add IP combined
	s += '''\t\t<Track altColor="0,0,178" autoScale="true" clazz="org.broad.igv.track.DataSourceTrack" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="{id}" name="{name}" normalize="false" renderer="BAR_CHART" sortable="true" visible="true" windowFunction="mean">\n'''.format(
		id=avail_path[name][0],
		name=name+'.IP'		)
	s += '''\t\t\t<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="200.0" minimum="0.0" type="LINEAR"/>\n'''
	s +='\t\t</Track>\n'
	# add IP unique
	s += '''\t\t<Track altColor="0,0,178" autoScale="true" clazz="org.broad.igv.track.DataSourceTrack" color="0,178,0" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="{id}" name="{name}" normalize="false" renderer="BAR_CHART" sortable="true" visible="true" windowFunction="mean">\n'''.format(
		id=avail_path[name][2],
		name=name+'.IP-uniq')
	s += '''\t\t\t<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="200.0" minimum="0.0" type="LINEAR"/>\n'''
	s +='\t\t</Track>\n'

	# add Input combined
	s += '''\t\t<Track altColor="0,0,178" autoScale="true" clazz="org.broad.igv.track.DataSourceTrack" color="178,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="{id}" name="{name}" normalize="false" renderer="BAR_CHART" sortable="true" visible="true" windowFunction="mean">\n'''.format(
		id=avail_path[name][1],
		name=name+'.Inp'	)
	s += '''\t\t\t<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="200.0" minimum="0.0" type="LINEAR"/>\n'''
	s+='\t\t</Track>\n'

s += '\t</Panel>\n'

s += '\t<Panel height="90" name="FeaturePanel" width="1200">\n'
# add default Reference sequence and RefSeq
s += '''\t\t<Track altColor="0,0,178" autoScale="false" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="Reference sequence" name="Reference sequence" sortable="false" visible="true"/>\n'''
s += '''\t\t<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" colorScale="ContinuousColorScale;0.0;308.0;255,255,255;0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="35" id="hg19_genes" name="RefSeq Genes" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count">\n'''
s += '''\t\t\t<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="308.0" minimum="0.0" type="LINEAR"/>\n'''
s += '\t\t</Track>\n'

# add Gencode
s += '''<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="153,51,0" displayMode="COLLAPSED" featureVisibilityWindow="1000000" fontSize="10" id="/u/flashscratch/flashscratch1/f/frankwoe/m6A_CLAM/src/igv/hg19/gencode.v19.annotation.filtered.gtf" name="Gencode.V19" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/>\n'''

# add RepeatMasker BED
s += '''\t\t<Track altColor="102,102,102" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="102,102,102" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="10" id="/u/home/f/frankwoe/scratch/hg19_gtf_bed/hg19_repeatMasker.sorted.bed" name="RepeatMasker" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/>\n'''

# add RRACU Forward
s += '''\t\t<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="10000000" fontSize="10" id="[GA][GA]ACT" name="Forward RRACU" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count">\n'''
s += '''\t\t\t<SequenceMatchSource featureWindowSize="10000000" genome="hg19" pattern="[GA][GA]ACT" strand="POSITIVE"/>\n'''
s += '''\t\t</Track>\n'''
# add RRACU Reverse
s += '''\t\t<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="255,0,0" displayMode="COLLAPSED" featureVisibilityWindow="10000000" fontSize="10" id="[GA][GA]ACT Negative" name="Reverse RRACU" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count">\n'''
s += '''\t\t\t<SequenceMatchSource featureWindowSize="10000000" genome="hg19" pattern="[GA][GA]ACT" strand="NEGATIVE"/>\n'''
s += '''\t\t</Track>\n'''

s += '\t</Panel>\n'


s += '\t<PanelLayout dividerFractions="0.7"/>\n'
s += '\t<HiddenAttributes>\n'
s += '\t\t<Attribute name="NAME"/>\n'
s += '\t\t\t<Attribute name="DATA FILE"/>\n'
s += '\t\t\t<Attribute name="DATA TYPE"/>\n'
s += '\t</HiddenAttributes>\n'
s += '</Session>'


with open(out_fp, 'w') as fo:
	fo.write(s)
