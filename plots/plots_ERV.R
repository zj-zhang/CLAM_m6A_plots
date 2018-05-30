## make plots for the final data of peaks
## Zijun Zhang
## 5.23.2018

library(ggplot2)

data = read.table(gzfile('annot_peak.df.gz'), sep='\t', header=T, stringsAsFactors=F)

data_to_use = data[data$occurence_num>1,]
data_not_use = data[data$occurence_num==1,]


length(grep('ERV', data_to_use[,"repeat."]))
length(Reduce(intersect, list(
	grep('ERV', data_to_use[,"repeat."]), 
	which(data_to_use$intron!=1),
	which(data_to_use$CDS_only+data_to_use$UTR_only+data_to_use$CDS_and_UTR >0) ) ))
idx = Reduce(intersect, list(
	grep('ERV', data_to_use[,"repeat."]), 
	which(data_to_use$intron!=1),
	which(data_to_use$CDS_only+data_to_use$UTR_only+data_to_use$CDS_and_UTR >0),
	which(data_to_use$motif.RRACT>0)
	) 
	)
length(idx)


length(grep('ERV', data_not_use[,"repeat."]))
length(Reduce(intersect, list(
	grep('ERV', data_not_use[,"repeat."]), 
	which(data_not_use$intron!=1),
	which(data_not_use$CDS_only+data_not_use$UTR_only+data_not_use$CDS_and_UTR >0) ) ))

idx2 = Reduce(intersect, list(
	grep('ERV', data_not_use[,"repeat."]), 
	which(data_not_use$intron!=1),
	which(data_not_use$CDS_only+data_not_use$UTR_only+data_not_use$CDS_and_UTR >0),
	which(data_not_use$motif.RRACT>0)
	) 
	)
length(idx2)

## get examples

get_commandline = function(this, name) cat('python parse_session.py', this$occurence, this$peak, name, sep=' ', file='cmd.txt')

alu_df = data_to_use[idx,]

## YY1, in 6 exp 12 rep
this = alu_df[alu_df$gene=="ENSG00000100811",]
get_commandline(this, 'YY1.xml')

## PHLPP2, 1 exp 2 rep
this = alu_df[alu_df$gene=="ENSG00000040199",]
get_commandline(this, 'PHLPP2.xml')

## IRF3, in 2 exp 2 rep
this = alu_df[alu_df$gene=="ENSG00000126456",] 
get_commandline(this, 'IRF3.xml')

## SRSF5, in 1 exp 2 rep
this = alu_df[alu_df$gene=="ENSG00000100650",]  
get_commandline(this, 'SRSF5.xml')

## BRCA1, in 5 exp 15 rep
this = alu_df[alu_df$gene=="ENSG00000012048",]  
get_commandline(this, 'BRCA1.xml')

## FGFR1, in 4 exp 8 rep
this = alu_df[alu_df$gene=="ENSG00000077782",]  
get_commandline(this, 'FGFR1.xml')

## CPM, in 4 exp 6 rep
this = alu_df[alu_df$gene=="ENSG00000135678",]  
get_commandline(this, 'CPM.xml')

alu_df[alu_df$gene=="ENSG00000244038",]  ## DDOST, in 3 exp 5 rep




## get gene sets

geneset = unique(data_to_use$gene[idx])
geneset = as.data.frame(geneset, stringsAsFactors=F)
write.table(geneset, 'Alu_genes.txt', quote=F, row.names=F, col.names=F)
	

geneset2 = unique(data_not_use$gene[idx2])
geneset2 = as.data.frame(geneset2, stringsAsFactors=F)
write.table(geneset2, 'Alu_negative_genes.txt', quote=F, row.names=F, col.names=F)


## plot conservation
pdf('cons_ecdf.pdf', width=8, height=4)
par(mfrow=c(1,2))
plot(ecdf(data_to_use$avg_cons[idx]), main='Alu-m6A peaks', xlim=c(-2, 7), 
	do.points=F, xlab='Avg.PhyloP')
lines(ecdf(data_not_use$avg_cons[idx2]), col='red', do.points=F)
#legend('bottomright', legend=c('Positive', 'Negative'), col=c('black', 'red'), lwd=1)

plot(ecdf(data_to_use$avg_cons), main='All peaks', xlim=c(-2, 7), 
	do.points=F, xlab='Avg.PhyloP')
lines(ecdf(data_not_use$avg_cons), col='red', do.points=F)
legend('bottomright', legend=c('Peak occurrence>1', 'Peak occurrence=1'), col=c('black', 'red'), lwd=1)

dev.off()

## plot snp_num
pdf('snp_ecdf.pdf', width=8, height=4)
par(mfrow=c(1,2))
plot(ecdf(data_to_use$snp_num[idx]), main='Alu-m6A peaks', xlim=c(-2, 7), 
	do.points=F, xlab='SNP num')
lines(ecdf(data_not_use$snp_num[idx2]), col='red', do.points=F)
#legend('bottomright', legend=c('Positive', 'Negative'), col=c('black', 'red'), lwd=1)

plot(ecdf(data_to_use$snp_num), main='All peaks', xlim=c(-2, 7), 
	do.points=F, xlab='SNP num')
lines(ecdf(data_not_use$snp_num), col='red', do.points=F)
legend('bottomright', legend=c('Peak occurrence>1', 'Peak occurrence=1'), col=c('black', 'red'), lwd=1)

dev.off()


