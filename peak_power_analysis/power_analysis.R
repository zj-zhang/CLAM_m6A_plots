## some plots to show the relationship between peaks and depth
## Zijun Zhang
## 5.16.2018

library(LSD)
library(ggplot2)
library(Darts)


heat_plot = function(x, y, xlab='', ylab='', xlim=c(0,1), ylim=c(0,1))
{

	library(MASS)
	DF <- data.frame(x,y)
	dens <- kde2d(x,y)
	gr <- data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z))
	names(gr) <- c("xgr", "ygr", "zgr")

	# Fit a model
	mod <- loess(zgr~xgr*ygr, data=gr)

	# Apply the model to the original data to estimate density at that point
	
	DF$density <- predict(mod, newdata=data.frame(xgr=x, ygr=y))
	
	# Draw plot
	library(ggplot2)
	library(RColorBrewer)
	
	myPalette = colorRampPalette(c("darkgrey", "darkblue", "red", "orange"))
	sc = scale_colour_gradientn(colours = myPalette(100))
	#lab = paste0("Spearman~rho==", round(cor(x,y, method='sp'),2))
	lab = paste0("Spearman's Rho=", round(cor(x,y, method='sp'),2))
	
	p = ggplot(DF, aes(x=x,y=y, color=density)) + geom_point() + 
		sc + 
		geom_abline(intercept=0, slope=1, linetype='dashed')+
		myTheme+
		no_border +
		coord_cartesian(xlim=xlim, ylim=ylim) +
		xlab(xlab) +
		ylab(ylab) +
		annotate('text', x=max(DF$x),y=min(DF$y),hjust=1, vjust=0, label=lab) +
		theme(legend.position=c(1,0), legend.justification=c(1,0)) + 
		guides(color=F)
	return(p)
}



data=read.table('peak_vs_depth.csv')

#png('Depth_v_Upeak.png', width=100, height=100)
pdf('Depth_v_Upeak.pdf', width=6, height=6)
heatscatter(log10(data$ip_depth), log10(data$unique_peak_num), cor=T, xlab='log10(IP coverage)', ylab='log10(#unique peaks)',
	main='Depth vs Unique peaks', xlim=c(5,10), ylim=c(0,7))
dev.off()


pdf('Mread_depth.pdf', width=6, height=6)
rescued_coverage = log10(data$ip_depth*data$ip_multi/100)
hist(rescued_coverage,  prob=T, breaks=10, xlab='log10(Rescued IP coverage)', main='Histogram of multi-mapped reads', col='grey')
lines(density(rescued_coverage), col='blue', lwd=2)
lines(density(rescued_coverage, adjust=2), lty="dotted", col="darkgreen", lwd=2)
dev.off()


pdf('RescuedDepth_v_Mpeak.pdf', width=6, height=6)
heatscatter(rescued_coverage, log10(data$rescued_peak_num), cor=T, xlab='log10(rescued IP coverage)', ylab='log10(#rescued peaks)', xlim=c(4,8), ylim=c(0,6),
	main='Rescued Depth vs. rescued peaks')
dev.off()


