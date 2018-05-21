## plot the ratio of motif over peak confidence
## Zijun Zhang
## 5.17.2018

library(ggplot2)


prepare_data = function() {
	data = read.table(gzfile('peak_sum.tsv.gz'), header=T, sep='\t', stringsAsFactors=F)
	total_exp = unique( unlist( sapply(data$exp, function(x) strsplit(as.character(x), ",")[[1]] ) ) )

	N = nrow(data)

	motif_df = data.frame(peak=rep(NA,N), exp_num=rep(0,N), motif1=rep(0,N), motif2=rep(0,N),
		motif3=rep(0,N), stringsAsFactors=F)
	exp_df = matrix(0, nrow=length(total_exp), ncol=length(total_exp))
	rownames(exp_df) = total_exp
	colnames(exp_df) = total_exp

	## takes a huge amount of time; don't run
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


	write.table(exp_df, 'exp_df.tsv', sep='\t', col.names=T, quote=F)
	write.table(motif_df, 'motif_df.tsv', sep='\t', quote=F)
}

plot_motif_dist = function() {	
	motif_df = read.table(gzfile('motif_df.tsv.gz', 'r'), header=T, sep='\t', stringsAsFactors=F)
	max_n = 35
	plot_df = data.frame(
		exp_num=seq(1,max_n),
		peak_num = rep(0,max_n),
		mean_motif1=rep(0,max_n), sd_motif1=rep(0,max_n),
		prop_motif1=rep(0,max_n), 
		mean_motif2=rep(0,max_n), sd_motif2=rep(0,max_n),
		prop_motif2=rep(0,max_n), 
		mean_motif3=rep(0,max_n), sd_motif3=rep(0,max_n),
		prop_motif3=rep(0,max_n) )
	for(i in 1:max_n)
	{
		if(i!=max_n) {
			this = motif_df[motif_df$occurence_num==i,]
		} else {
			this = motif_df[motif_df$occurence_num>=i,]
		}
		plot_df$peak_num[i] = nrow(this)
		
		plot_df$prop_motif1[i] = mean(this$motif1>0)
		plot_df$mean_motif1[i] = mean(this$motif1[this$motif1>0])
		plot_df$sd_motif1[i] = sd(this$motif1[this$motif1>0])
		
		plot_df$prop_motif2[i] = mean(this$motif2>0)
		plot_df$mean_motif2[i] = mean(this$motif2[this$motif2>0])
		plot_df$sd_motif2[i] = sd(this$motif2[this$motif2>0])
		
		plot_df$prop_motif3[i] = mean(this$motif3>0)
		plot_df$mean_motif3[i] = mean(this$motif3[this$motif3>0])
		plot_df$sd_motif3[i] = sd(this$motif3[this$motif3>0])
	}
	
	
	pdf('motif_analysis.pdf', width=12, height=6)
	par(mfrow=c(1,2))
	## plot the proportion of peaks with motif
	plot(-plot_df$exp_num, plot_df$prop_motif1, type='b', ylim=c(0,1),
		xlab='No. of occurrences (replicates)', ylab='% peaks with motif',
		main='Percentage of peaks with Motifs')
	lines(-plot_df$exp_num, plot_df$prop_motif2, type='b', col='blue')
	lines(-plot_df$exp_num, plot_df$prop_motif3, type='b', col='red')
	legend('topright', legend=c("GGACU","RRACU","DRACH"), col=c('black','blue','red'), lwd=1)
	abline(h=0.080, col='black', lty=2)
	abline(h=0.338, col='blue', lty=2)
	abline(h=0.816, col='red', lty=2)
	
	## plot the average motif count for peaks with motif
	plot(-plot_df$exp_num, plot_df$mean_motif1, type='b', ylim=c(1,3),
		xlab='No. of occurrences (replicates)', ylab='Avg. No. of motifs per peak',
		main='Average number of Motifs in peak with motif')
	lines(-plot_df$exp_num, plot_df$mean_motif2, type='b', col='blue')
	lines(-plot_df$exp_num, plot_df$mean_motif3, type='b', col='red')
	legend('topright', legend=c("GGACU","RRACU","DRACH"), col=c('black','blue','red'), lwd=1)
	dev.off()
	
	pdf('peak_num_cdf.pdf', width=4, height=4)
	plot(ecdf(motif_df$occurence_num), verticals=T, do.points=T, col.01line=NULL, main='CDF of peak occurrence number',
		xlab='No. of Occurrence (replicates)')
	dev.off()
	
	#gg_df = melt(plot_df, id.vars='exp_num')
	#p = ggplot(gg_df, aes(x=exp_num, y=value, group=variable))
	
}