## plot the ratio of motif over peak confidence
## Zijun Zhang
## 5.17.2018

library(ggplot2)

data = read.table(gzfile('peak_sum.tsv.gz'), header=T, sep='\t', stringsAsFactors=F)

total_exp = unique( unlist( sapply(data$exp, function(x) strsplit(as.character(x), ",")[[1]] ) ) )

N = nrow(data)
motif_df = data.frame(peak=rep(NA,N), exp_num=rep(0,N), motif1=rep(0,N), motif2=rep(0,N),
	motif3=rep(0,N), stringsAsFactors=F)
exp_df = matrix(0, nrow=length(total_exp), ncol=length(total_exp))
rownames(exp_df) = total_exp
colnames(exp_df) = total_exp

for(i in 1:nrow(data))
{
	if(! i%%10000) print(i)
	this = data[i,]
	exps = strsplit(as.character(this$exp), ",")[[1]]
	motifs = as.numeric(strsplit(as.character(this[5]), ",")[[1]])
	
	#motif_list = list(peak=this[1], exp_num=length(exps), motif1=motifs[1], motif2=motifs[2], motif3=motifs[3])
	
	#motif_df = rbind.data.frame(motif_df, motif_list, stringsAsFactors=F)
	motif_df[i,] = c(peak=this[1], exp_num=length(exps), motif1=motifs[1], motif2=motifs[2], motif3=motifs[3])
	idx = which(total_exp %in% exps)
	if(length(idx)==1) {
		exp_df[idx, idx] = exp_df[idx, idx] + 1
	} else{
		for(m in 1:length(idx)) {
			for(n in m:length(idx)) {
				exp_df[idx[m], idx[n]] = exp_df[idx[m], idx[n]] + 1
				exp_df[idx[n], idx[m]] = exp_df[idx[n], idx[m]] + 1
			}
		}
	}
}


