#!/usr/bin/env Rscript

source('scripts/_include_options.R')

e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')

mash.results = readRDS('checkpoints/cellprop_mashr_results.rds')
emma.results = readRDS('checkpoints/cellprop_emma_results.rds')
keep.genes = readRDS('checkpoints/keep_genes.rds')

emma.beta = emma.results[,paste('beta',predictor,sep='.'),]
emma.pval = emma.results[,paste('pval',predictor,sep='.'),]
emma.qval = apply(emma.pval,2,function(x) p.adjust(x,'fdr'))

library(mashr)
mash.beta = get_pm(mash.results)
mash.lfsr = get_lfsr(mash.results)

# Significant genes per region
apply(emma.qval,2,function(x) sum(x < fsr.cutoff,na.rm=TRUE))
apply(mash.lfsr,2,function(x) sum(x < fsr.cutoff))

# How many genes are significant in n regions
table(apply(mash.lfsr,1,function(x) sum(x < fsr.cutoff)))

# Genes that are significant in just one region
apply(mash.lfsr[apply(mash.lfsr,1,function(x) sum(x < fsr.cutoff)) == 1,],2,function(x) sum(x < fsr.cutoff))

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

# Histogram of p values for each region
p = ggplot(melt(emma.pval),aes(value,fill=Var2)) +
	geom_histogram(binwidth=0.005) +
	facet_wrap(~Var2,nrow=region.rows) +
	scale_x_continuous(
		limits = c(0,1),
		breaks = seq(0,1,0.1),
		labels = c('0.0',rep('',4),'0.5',rep('',4),'1.0')
	) +
	scale_fill_manual(values=region.colors) +
	xlab('p value') + ylab('Count') + theme_classic() +
	theme(legend.position='none')
ggsave(p,file=paste0('figures/cellprop_model_results_pval_enrichment_emma_',tolower(predictor.label),'_histogram.pdf'),useDingbats=FALSE,height=5)

# QQ-plot of p values for each region
# p = ggplot(melt(emma.pval),aes(sample=value,color=Var2)) +
# 	stat_qq(distribution=stats::qunif) +
# 	geom_abline(slope=1,col='black',size=0.2) +
# 	facet_wrap(~Var2,nrow=region.rows) +
# 	coord_fixed() +
# 	scale_color_manual(values=region.colors) +
# 	scale_x_continuous(breaks=seq(0,1,0.1),labels=c(0,rep('',9),1)) +
# 	scale_y_continuous(breaks=seq(0,1,0.1),labels=c(0,rep('',9),1)) +
# 	xlab('Expected') + ylab('Observed') + theme_classic() +
# 	theme(legend.position='none')
# ggsave(p,file=paste0('figures/cellprop_model_results_pval_enrichment_emma_',tolower(predictor.label),'_qqplot.pdf'),useDingbats=FALSE,height=5)

# QQ plot
p = ggplot(
		do.call(rbind,lapply(split(melt(emma.pval),melt(emma.pval)$Var2),function(x) {
			data.frame(
				region=unique(x$Var2),
				expected=-log10(seq(1/sum(complete.cases(x)),1,1/sum(complete.cases(x)))),
				observed=-log10(quantile(x$value,seq(1/sum(complete.cases(x)),1,1/sum(complete.cases(x))),na.rm=TRUE))
			)
		})),
		aes(expected,observed,color=region)
	) +
	geom_smooth(method=lm,se=FALSE,color='black',size=0.2,linetype=2) +
	geom_point(size=0.5) +
	geom_abline(slope=1,col='black',size=0.2) +
	facet_wrap(~region,nrow=region.rows) +
#	coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
	scale_color_manual(values=region.colors) +
 	scale_x_continuous(breaks=seq(0,ceiling(-log10(min(1/(max(table(melt(emma.pval)$Var2),na.rm=TRUE)-1)))),1)) +
 	scale_y_continuous(breaks=seq(0,ceiling(max(-log10(melt(emma.pval)$value),na.rm=TRUE)),1)) +
	xlab('Expected') + ylab('Observed') + theme_classic() +
	theme(legend.position='none')
ggsave(p,file=paste0('figures/cellprop_model_results_pval_enrichment_emma_',tolower(predictor.label),'_qqplot.pdf'),useDingbats=FALSE,height=5)


# Put together longform data frame of emma results
b.emma = data.frame(
	expand.grid(dimnames(emma.results)[[1]],names(keep.genes)),
	pval=as.numeric(emma.results[,paste('pval',predictor,sep='.'),]),
	qval=as.numeric(apply(emma.results[,paste('pval',predictor,sep='.'),],2,p.adjust,method='fdr')),
	beta=as.numeric(emma.results[,paste('beta',predictor,sep='.'),])
)

p = ggplot(subset(b.emma,qval < fsr.cutoff),aes(beta,color=Var2)) +
	scale_color_manual(name='Region',values=region.colors) +
	geom_density() +
	theme_classic(base_size=12) +
	guides(color = guide_legend(ncol = 2)) +
	xlab(expression(italic(beta))) +
	ylab('Density')
ggsave(p,file=paste0('figures/cellprop_model_results_beta_density_emma_',tolower(predictor.label),'.pdf'),width=7,height=3,useDingbats=FALSE)

p = ggplot(melt(tapply(b.emma$qval,b.emma$Var2,function(x) sum(x < fsr.cutoff,na.rm=TRUE))),aes(Var1,value,fill=Var1)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Region',values=region.colors) +
	theme_classic(base_size=12) +
	theme(legend.position='none',axis.text.x=element_text(angle=-45,hjust=0,vjust=1)) +
	xlab('Regions') +
	ylab('Number of genes')
ggsave(p,file=paste0('figures/cellprop_model_results_beta_count_emma_',tolower(predictor.label),'.pdf'),width=7,height=3,useDingbats=FALSE)

# Put together longform data frame of mashr results
b.mash = data.frame(
	expand.grid(rownames(mash.beta),colnames(mash.beta)),
	qval=as.numeric(mash.lfsr),
	beta=as.numeric(mash.beta),
	stringsAsFactors=FALSE
)

p = ggplot(subset(b.mash,qval < fsr.cutoff),aes(beta,color=Var2)) +
	scale_color_manual(name='Region',values=region.colors) +
	geom_density() +
	theme_classic(base_size=12) +
	guides(color = guide_legend(ncol = 2)) +
	xlab(expression(italic(beta))) +
	ylab('Density')
ggsave(p,file=paste0('figures/cellprop_model_results_beta_density_mash_',tolower(predictor.label),'.pdf'),width=7,height=3,useDingbats=FALSE)

p = ggplot(melt(tapply(b.mash$qval,b.mash$Var2,function(x) sum(x < fsr.cutoff,na.rm=TRUE))),aes(Var1,value,fill=Var1)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Region',values=region.colors) +
	theme_classic(base_size=12) +
	theme(legend.position='none',axis.text.x=element_text(angle=-45,hjust=0,vjust=1)) +
	xlab('Regions') +
	ylab('Number of genes')
ggsave(p,file=paste0('figures/cellprop_model_results_beta_count_mash_',tolower(predictor.label),'.pdf'),width=7,height=3,useDingbats=FALSE)

b.emma$qval.signed = with(b.emma,ifelse(beta>0,1,-1) * qval)
b.mash$qval.signed = with(b.mash,ifelse(beta>0,1,-1) * qval)

# Put together longform data frame with gene counts (up vs. downregulated, split by method)
model.counts.combined = rbind(
	within(melt(tapply(b.emma$qval.signed,b.emma$Var2,function(x) -sum(abs(x) < fsr.cutoff & x < 0,na.rm=TRUE))), {method='EMMA'; direction='down'} ),
	within(melt(tapply(b.emma$qval.signed,b.emma$Var2,function(x) sum(abs(x) < fsr.cutoff & x >= 0,na.rm=TRUE))), {method='EMMA'; direction='up'} ),
	within(melt(tapply(b.mash$qval.signed,b.mash$Var2,function(x) -sum(abs(x) < fsr.cutoff & x < 0,na.rm=TRUE))), {method='MASH'; direction='down'} ),
	within(melt(tapply(b.mash$qval.signed,b.mash$Var2,function(x) sum(abs(x) < fsr.cutoff & x >= 0,na.rm=TRUE))), {method='MASH'; direction='up'} )
)

ylimit = ceiling(with(model.counts.combined,max(abs(value)))/500) * 500

levels(model.counts.combined$Var1) = gsub('ACCg','ACC',levels(model.counts.combined$Var1))

# Note that the code below will only format correctly with ggplot >= 3.3.0
if (packageVersion('ggplot2') < 3.3) warning('Some plots may not display correctly with ggplot2 version < 3.3.0')

p = ggplot(model.counts.combined,aes(Var1,value,fill=Var1,alpha=direction)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Region',values=region.colors) +
	scale_alpha_manual(values=c(0.75,1)) +
	scale_y_continuous(
		limits = c(-ylimit,ylimit),
		breaks = c(-ylimit,-ylimit*0.5,0,ylimit*0.5,ylimit),
		labels = c(formatC(ylimit,width=5,flag=' '),'Decrease',formatC(0,width=5,flag=' '),'Increase',formatC(ylimit,width=5,flag=' '))
	) +
	theme_classic(base_size=12) +
	theme(
		axis.ticks.y = element_line(linetype=c(1,0,1,0,1)),
		axis.title.x = element_blank(),
		axis.title.y = element_text(),
		axis.text.x = element_text(
			angle = -45, hjust = 0, vjust = 1
		),
		axis.text.y=element_text(
			face = c('plain','bold','plain','bold','plain'),
#			size = axis.text.size * c(1,2,1,2,1),
			angle = c(0,90,0,90,0), hjust=0.5
		)
	) +
	facet_wrap(~method,nrow=2) +
	theme(legend.position='none') +
	ylab('Number of genes')
ggsave(p,file=paste0('figures/cellprop_model_results_beta_count_comparison_',tolower(predictor.label),'.pdf'),width=7,height=5,useDingbats=FALSE)

p = ggplot(droplevels(subset(model.counts.combined,method == 'MASH')),aes(Var1,value,fill=Var1,alpha=direction)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Region',values=region.colors) +
	scale_alpha_manual(values=c(0.75,1)) +
	scale_y_continuous(
		limits = c(-ylimit,ylimit),
		breaks = c(-ylimit,-ylimit*0.5,0,ylimit*0.5,ylimit),
		labels = c(formatC(ylimit,width=5,flag=' '),'Decrease',formatC(0,width=5,flag=' '),'Increase',formatC(ylimit,width=5,flag=' '))
	) +
	theme_classic(base_size=24) +
	theme(
		axis.ticks.y = element_line(linetype=c(1,0,1,0,1)),
		axis.title.x = element_blank(),
		axis.title.y = element_text(),
		axis.text.x = element_text(
			angle = -45, hjust = 0, vjust = 1
		),
		axis.text.y=element_text(
			face = c('plain','bold','plain','bold','plain'),
#			size = axis.text.size * c(1,2,1,2,1),
			angle = c(0,90,0,90,0), hjust=0.5
		)
	) +
	theme(legend.position='none') +
	ylab('Number of genes')
ggsave(p,file=paste0('figures/cellprop_model_results_beta_count_mash_',tolower(predictor.label),'.pdf'),width=7,height=5,useDingbats=FALSE)



# pair.share = get_pairwise_sharing(mash.results, factor = 0.5)

# Metadata counting of significant effects

b.mash.split.genes = split(b.mash,b.mash$Var1)

# Tally up sig region-combos, as well as up-regulated and down-regulated
region.combinations = table(unlist(lapply(b.mash.split.genes,function(x) {
	paste(subset(x,qval < fsr.cutoff & beta != 0)$Var2,collapse='-')
})))

# Get rid of blanks
region.combinations = region.combinations[as.logical(nchar(names(region.combinations)))]

# Sort in same order (by total)
region.combinations = region.combinations[order(region.combinations,decreasing=TRUE)]

region.combinations.inc = region.combinations.dec = integer(length(region.combinations))
names(region.combinations.inc) = names(region.combinations.dec) = names(region.combinations)

# Ensure all have the same order
region.combinations.inc[names(region.combinations)] = table(unlist(lapply(b.mash.split.genes,function(x) {
	paste(subset(x,qval < fsr.cutoff & beta > 0)$Var2,collapse='-')
})))[names(region.combinations)]

region.combinations.dec[names(region.combinations)] = table(unlist(lapply(b.mash.split.genes,function(x) {
	paste(subset(x,qval < fsr.cutoff & beta < 0)$Var2,collapse='-')
})))[names(region.combinations)]

# Set NAs to 0
region.combinations.inc[is.na(region.combinations.inc)] = 0
region.combinations.dec[is.na(region.combinations.dec)] = 0

# The sum of up-regulated genes and down-regulated genes do not add up to the sum of significant genes
# This is because genes can have different directions, causing their regions to be tallied differently.
# Thus, report based on the sum of up- and down-regulated region combos
region.combinations.sum = region.combinations.inc + region.combinations.dec

region.combinations.sum = sort(region.combinations.sum,decreasing=TRUE)
region.combinations = region.combinations[names(region.combinations.sum)]
region.combinations.inc = region.combinations.inc[names(region.combinations.sum)]
region.combinations.dec = region.combinations.dec[names(region.combinations.sum)]

# Set NAs to 0
region.combinations.inc[is.na(region.combinations.inc)] = 0
region.combinations.dec[is.na(region.combinations.dec)] = 0

width.of.bars = 0.8
region.combinations.results = do.call(rbind,lapply(1:length(region.combinations),function(i) {
	x = names(region.combinations)[i]
	n = unlist(lapply(strsplit(names(region.combinations),'-'),length))[i]
	count.all = as.integer(region.combinations[x])
	count.inc = as.integer(region.combinations.inc[x])
	count.dec = as.integer(region.combinations.dec[x])
	out = integer(length(keep.genes))
	names(out) = c(names(keep.genes))
	out[unlist(strsplit(x,split='-'))] = 1
	out = rbind(
		data.frame(
			combination=i,
			n_regions = n,
			share_region = n > fraction.shared.cutoff * length(keep.genes),
			region=factor(region.levels,levels=region.levels),
			region_sig = factor(ifelse(as.logical(out),names(out),NA),levels=region.levels),
			value = 1,
			xmin = seq(1,length(keep.genes)) - (width.of.bars)/2,
			xmax = seq(1,length(keep.genes)) + (width.of.bars)/2,
			ymin = i - (width.of.bars)/2,
			ymax = i + (width.of.bars)/2,
			chart = 'meta'
		),
		data.frame(
			combination=i,
			n_regions = n,
			share_region = n > fraction.shared.cutoff * length(keep.genes),
			region=NA,
			region_sig=NA,
			value=count.all,
			xmin = NA,
			xmax = NA,
			ymin = NA,
			ymax = NA,
			chart='count_all'
		),
		data.frame(
			combination=i,
			n_regions = n,
			share_region = n > fraction.shared.cutoff * length(keep.genes),
			region=NA,
			region_sig=NA,
			value=count.inc,
			xmin = NA,
			xmax = NA,
			ymin = NA,
			ymax = NA,
			chart='count_increase'
		),
		data.frame(
			combination=i,
			n_regions = n,
			share_region = n > fraction.shared.cutoff * length(keep.genes),
			region=NA,
			region_sig=NA,
			value=count.dec,
			xmin = NA,
			xmax = NA,
			ymin = NA,
			ymax = NA,
			chart='count_decrease'
		)
	)
	out
}))

max.plot = length(keep.genes) * 2
base.color = '#000000'
fade.color = '#000000'
highlight.color = '#ff0000'
base.size = 14

meta.combinations.results.plot = subset(region.combinations.results,combination < max.plot)

p1 = ggplot(subset(meta.combinations.results.plot,chart=='meta')) +
	geom_blank(aes(x = region,y=combination)) +
	geom_rect(
		aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=region,alpha=region_sig,color=region_sig),
		linetype=1
	) +
	scale_fill_manual(values=region.colors) +
	scale_alpha_manual(values=rep(1,length(region.colors)),na.value=0) +
	scale_color_manual(values=rep('black',length(region.colors)),na.value=0) +
	scale_y_continuous(trans='reverse',expand=c(0,0)) +
	scale_x_discrete(position='top',expand=c(0,0)) +
	coord_equal() +
	theme_classic(base_size=base.size) + 
	theme(
		legend.position='none',
		axis.line=element_blank(),
		axis.title=element_blank(),
		axis.text.x=element_text(angle=45,vjust=0,hjust=0,margin=margin(t=-1)),
		axis.text.y=element_blank(),
		axis.ticks=element_blank()
	)
p2 = ggplot(subset(meta.combinations.results.plot,chart=='count_increase')) +
	geom_bar(aes(combination,value),stat='identity',fill=base.color,width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0),limits=c(0,ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100)) +
	coord_flip() +
	ylab('Upregulated') +
	theme_classic(base_size=base.size) +
	theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
p3 = ggplot(subset(meta.combinations.results.plot,chart=='count_decrease')) +
	geom_bar(aes(combination,value),stat='identity',fill=base.color,width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(trans='reverse',expand=c(0,0),limits=c(ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100,0)) +
	coord_flip() +
	ylab('Downregulated') +
	theme_classic(base_size=base.size) +
	theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())

library(egg)

pdf(file=paste0('figures/cellprop_model_results_beta_count_comparison_',tolower(predictor.label),'_by_meta.pdf'),useDingbats=FALSE,height=7,width=11)
	ggarrange(p3,p1,p2,ncol=3,nrow=1,widths=c(1,1,1),newpage=FALSE)
dev.off()

p4 = ggplot(subset(meta.combinations.results.plot,chart=='count_increase')) +
	geom_bar(aes(combination,value,fill=share_region),stat='identity',width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0),limits=c(0,ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100)) +
	scale_fill_manual(values=c(fade.color,highlight.color)) +
	coord_flip() +
	ylab('Upregulated') +
	theme_classic(base_size=base.size) +
	theme(legend.position='none',axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
p5 = ggplot(subset(meta.combinations.results.plot,chart=='count_decrease')) +
	geom_bar(aes(combination,value,fill=share_region),stat='identity',width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(trans='reverse',expand=c(0,0),limits=c(ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100,0)) +
	scale_fill_manual(values=c(fade.color,highlight.color)) +
	coord_flip() +
	ylab('Downregulated') +
	theme_classic(base_size=base.size) +
	theme(legend.position='none',axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())

pdf(file=paste0('figures/cellprop_model_results_beta_count_comparison_',tolower(predictor.label),'_by_meta_shared.pdf'),useDingbats=FALSE,height=7,width=11)
	ggarrange(p5,p1,p4,ncol=3,nrow=1,widths=c(1,1,1),newpage=FALSE)
dev.off()

p6 = ggplot(subset(meta.combinations.results.plot,chart=='count_increase')) +
	geom_bar(aes(combination,value,fill=n_regions==1),stat='identity',width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0),limits=c(0,ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100)) +
	scale_fill_manual(values=c(fade.color,highlight.color)) +
	coord_flip() +
	ylab('Upregulated') +
	theme_classic(base_size=base.size) +
	theme(legend.position='none',axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
p7 = ggplot(subset(meta.combinations.results.plot,chart=='count_decrease')) +
	geom_bar(aes(combination,value,fill=n_regions==1),stat='identity',width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(trans='reverse',expand=c(0,0),limits=c(ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100,0)) +
	scale_fill_manual(values=c(fade.color,highlight.color)) +
	coord_flip() +
	ylab('Downregulated') +
	theme_classic(base_size=base.size) +
	theme(legend.position='none',axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())

pdf(file=paste0('figures/cellprop_model_results_beta_count_comparison_',tolower(predictor.label),'_by_meta_single.pdf'),useDingbats=FALSE,height=7,width=11)
	ggarrange(p7,p1,p6,ncol=3,nrow=1,widths=c(1,1,1),newpage=FALSE)
dev.off()

save(list=ls(),file='checkpoints/checkpoint_cellprop_model_results_visualization.RData')

# Pull in pre-deconvolution MASH results
m.results = readRDS('checkpoints/mashr_results.rds')

m.beta = get_pm(m.results)
m.lfsr = get_lfsr(m.results)
m.sbet = m.beta / get_psd(m.results)

mash.sbet = mash.beta / get_psd(mash.results)

m1 = mutate(reshape2::melt(m.sbet),method='bulk')
m2 = mutate(reshape2::melt(mash.sbet),method='deconvoluted')

m1$deconvoluted = m2$value

levels(m1$Var2) = gsub('ACCg','ACC',levels(m1$Var2))

library(egg)
library(ggrastr)

limit = with(m1,max(abs(c(value,deconvoluted))))

p = ggplot(m1,aes(value,deconvoluted,color=Var2)) +
	geom_smooth(se=FALSE,method=lm,color='black',size=0.5) +
	geom_hline(yintercept=0,linetype=2,size=0.25) +
	geom_vline(xintercept=0,linetype=2,size=0.25) +
	geom_abline(intercept=0,slope=1,linetype=2,size=0.25) +
	geom_point_rast(size=0.1,alpha=0.02,shape=21,fill=NA) +
	facet_wrap(~Var2,nrow=region.rows) +
	scale_x_continuous(limits=c(-limit,limit),breaks=seq(-5,5,5)) +
	scale_y_continuous(limits=c(-limit,limit),breaks=seq(-5,5,5)) +
	theme_article(base_size=16) +
	coord_fixed() +
	scale_color_manual(values=region.colors) +
	theme(legend.position='none') +
	xlab(expression('Standardized '*beta['BULK'])) + ylab(expression('Standardized '*beta['DECON.']))
ggsave(p,file=paste0('figures/cellprop_model_standardized_beta_comparison_',tolower(predictor.label),'.pdf'),width=7,height=5,useDingbats=FALSE)


# wbaDEGs

wbadegs = c(
	names(which(unlist(lapply(mashr.genes,function(x) {
		(sum(m.lfsr[x,] < fsr.cutoff) >= fraction.shared.cutoff * length(keep.genes)) && sum(m.beta[x,][m.lfsr[x,] < fsr.cutoff] > 0) >= fraction.shared.cutoff * length(keep.genes)
	})))),
	names(which(unlist(lapply(mashr.genes,function(x) {
		(sum(m.lfsr[x,] < fsr.cutoff) >= fraction.shared.cutoff * length(keep.genes)) && sum(m.beta[x,][m.lfsr[x,] < fsr.cutoff] < 0) >= fraction.shared.cutoff * length(keep.genes)
	}))))
)

p = ggplot(subset(m1,Var1 %in% wbadegs),aes(value,deconvoluted,color=Var2)) +
	geom_smooth(se=FALSE,method=lm,color='black',size=0.5) +
	geom_hline(yintercept=0,linetype=2,size=0.25) +
	geom_vline(xintercept=0,linetype=2,size=0.25) +
	geom_abline(intercept=0,slope=1,linetype=2,size=0.25) +
	geom_point_rast(size=0.1,alpha=0.1,shape=21,fill=NA) +
	facet_wrap(~Var2,nrow=region.rows) +
	scale_x_continuous(limits=c(-limit,limit),breaks=seq(-5,5,5)) +
	scale_y_continuous(limits=c(-limit,limit),breaks=seq(-5,5,5)) +
	theme_article(base_size=16) +
	coord_fixed() +
	scale_color_manual(values=region.colors) +
	theme(legend.position='none') +
	xlab(expression('Standardized '*beta['BULK'])) + ylab(expression('Standardized '*beta['DECON.']))
ggsave(p,file=paste0('figures/cellprop_model_standardized_beta_comparison_',tolower(predictor.label),'_wbadegs.pdf'),width=7,height=5,useDingbats=FALSE)
