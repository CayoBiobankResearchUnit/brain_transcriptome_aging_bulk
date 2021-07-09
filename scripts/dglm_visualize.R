#!/usr/bin/env Rscript

source('scripts/_include_options.R')

dglm.results = readRDS('checkpoints/dglm_results.rds')
mash.results = readRDS('checkpoints/dglm_mashr_results.rds')
keep.genes = readRDS('checkpoints/keep_genes.rds')
e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')

dglm.beta = dglm.results[,'beta',]
dglm.pval = dglm.results[,'pval',]
dglm.qval = apply(dglm.pval,2,function(x) p.adjust(x,'fdr'))
dglm.sbet = dglm.results[,'beta',] / sqrt(dglm.results[,'bvar',])

library(mashr)
mash.beta = get_pm(mash.results)
mash.lfsr = get_lfsr(mash.results)
mash.sbet = mash.beta / get_psd(mash.results)

# Significant genes per region
apply(dglm.qval,2,function(x) sum(x < fsr.cutoff,na.rm=TRUE))
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
p = ggplot(melt(dglm.pval),aes(value,fill=Var2)) +
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
ggsave(p,file=paste0('figures/model_results_dglm_pval_enrichment_dglm_',tolower(predictor.label),'_histogram.pdf'),useDingbats=FALSE,height=5)

# QQ-plot of p values for each region
# p = ggplot(melt(dglm.pval),aes(sample=value,color=Var2)) +
# 	stat_qq(distribution=stats::qunif) +
# 	geom_abline(slope=1,col='black',size=0.2) +
# 	facet_wrap(~Var2,nrow=region.rows) +
# 	coord_fixed() +
# 	scale_color_manual(values=region.colors) +
# 	scale_x_continuous(breaks=seq(0,1,0.1),labels=c(0,rep('',9),1)) +
# 	scale_y_continuous(breaks=seq(0,1,0.1),labels=c(0,rep('',9),1)) +
# 	xlab('Expected') + ylab('Observed') + theme_classic() +
# 	theme(legend.position='none')
# ggsave(p,file=paste0('figures/model_results_dglm_pval_enrichment_dglm_',tolower(predictor.label),'_qqplot.pdf'),useDingbats=FALSE,height=5)

library(ggrastr)

# QQ plot
p = ggplot(
		do.call(rbind,lapply(split(melt(dglm.pval),melt(dglm.pval)$Var2),function(x) {
			data.frame(
				region=factor(unique(x$Var2),levels=region.levels,labels=gsub('ACCg','ACC',region.levels)),
				expected=-log10(seq(1/sum(complete.cases(x)),1,1/sum(complete.cases(x)))),
				observed=-log10(quantile(x$value,seq(1/sum(complete.cases(x)),1,1/sum(complete.cases(x))),na.rm=TRUE))
			)
		})),
		aes(expected,observed,color=region)
	) +
	geom_smooth(method=lm,se=FALSE,color='black',size=0.2,linetype=2) +
	geom_point_rast(size=0.5) +
	geom_abline(slope=1,col='black',size=0.2) +
	facet_wrap(~region,nrow=region.rows) +
	coord_cartesian(xlim=c(0,4),ylim=c(0,10)) +
#	coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
	scale_color_manual(values=region.colors) +
 	scale_x_continuous(breaks=seq(0,ceiling(-log10(min(1/(max(table(melt(dglm.pval)$Var2),na.rm=TRUE)-1)))),1)) +
 	scale_y_continuous(breaks=seq(0,ceiling(max(-log10(melt(dglm.pval)$value),na.rm=TRUE)),1)) +
	xlab('Expected') + ylab('Observed') + theme_classic() +
	theme(legend.position='none')
ggsave(p,file=paste0('figures/model_results_dglm_pval_enrichment_dglm_',tolower(predictor.label),'_qqplot.pdf'),useDingbats=FALSE,height=5)

library(egg)
# QQ plot
p = ggplot(
		do.call(rbind,lapply(split(melt(dglm.pval),melt(dglm.pval)$Var2),function(x) {
			data.frame(
				region=factor(unique(x$Var2),levels=region.levels,labels=gsub('ACCg','ACC',region.levels)),
				expected=-log10(seq(1/sum(complete.cases(x)),1,1/sum(complete.cases(x)))),
				observed=-log10(quantile(x$value,seq(1/sum(complete.cases(x)),1,1/sum(complete.cases(x))),na.rm=TRUE))
			)
		})),
		aes(expected,observed,color=region)
	) +
#	geom_smooth(method=lm,se=FALSE,color='black',size=0.2,linetype=2) +
	geom_point_rast(size=0.5) +
	geom_abline(slope=1,col='black',size=0.2) +
	facet_wrap(~region,nrow=region.rows) +
	coord_cartesian(xlim=c(0,4),ylim=c(0,10)) +
#	coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
	scale_color_manual(values=region.colors) +
 	scale_x_continuous(breaks=seq(0,ceiling(-log10(min(1/(max(table(melt(dglm.pval)$Var2),na.rm=TRUE)-1)))),1)) +
 	scale_y_continuous(breaks=seq(0,ceiling(max(-log10(melt(dglm.pval)$value),na.rm=TRUE)),1)) +
	xlab(expression(-log[10] * ('Expected'~italic(p)))) +
	ylab(expression(-log[10] * ('Observed'~italic(p)))) +
	theme_article() +
	theme(legend.position='none')
ggsave(p,file=paste0('figures/model_results_dglm_pval_enrichment_dglm_',tolower(predictor.label),'_qqplot_article.pdf'),useDingbats=FALSE,height=5)


# Put together longform data frame of dglm results
b.dglm = data.frame(
	expand.grid(dimnames(dglm.results)[[1]],names(keep.genes)),
	pval=as.numeric(dglm.results[,'pval',]),
	qval=as.numeric(apply(dglm.results[,'pval',],2,p.adjust,method='fdr')),
	beta=as.numeric(dglm.results[,'beta',]),
	sbet=as.numeric(dglm.results[,'beta',] / sqrt(dglm.results[,'bvar',]))
)
levels(b.dglm$Var2) = gsub('ACCg','ACC',levels(b.dglm$Var2))

p = ggplot(subset(b.dglm,qval < fsr.cutoff),aes(beta,color=Var2)) +
	scale_color_manual(name='Region',values=region.colors) +
	geom_density() +
	theme_classic(base_size=12) +
	guides(color = guide_legend(ncol = 2)) +
	xlab(expression(italic(beta))) +
	ylab('Density')
ggsave(p,file=paste0('figures/model_results_dglm_beta_density_dglm_',tolower(predictor.label),'.pdf'),width=7,height=3,useDingbats=FALSE)

p = ggplot(melt(tapply(b.dglm$qval,b.dglm$Var2,function(x) sum(x < fsr.cutoff,na.rm=TRUE))),aes(Var1,value,fill=Var1)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Region',values=region.colors) +
	theme_classic(base_size=12) +
	theme(legend.position='none',axis.text.x=element_text(angle=-45,hjust=0,vjust=1)) +
	xlab('Regions') +
	ylab('Number of genes')
ggsave(p,file=paste0('figures/model_results_dglm_beta_count_dglm_',tolower(predictor.label),'.pdf'),width=7,height=3,useDingbats=FALSE)

# Put together longform data frame of mashr results
b.mash = data.frame(
	expand.grid(rownames(mash.beta),colnames(mash.beta)),
	qval=as.numeric(mash.lfsr),
	beta=as.numeric(mash.beta),
	sbet=as.numeric(mash.sbet),
	stringsAsFactors=FALSE
)
levels(b.mash$Var2) = gsub('ACCg','ACC',levels(b.mash$Var2))

p = ggplot(b.mash,aes(sbet,color=Var2)) +
	geom_density() +
	geom_vline(aes(xintercept=value,color=Var2),data=within(reshape2::melt(with(b.mash,tapply(sbet,Var2,median))),{Var2=factor(Var1,levels=levels(b.mash$Var2))}),show.legend=FALSE) +
	facet_wrap(~Var2,ncol=1) +
	scale_x_continuous(limits=c(with(b.mash,c(-1,1)*max(abs(sbet))))) +
	scale_color_manual(name='Region',values=region.colors) +
	theme_classic(base_size=12) +
	theme(
		strip.text=element_blank(),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank()
	) +
	guides(color = guide_legend(ncol = 1)) +
	xlab(expression('Standardized'~italic(beta))) +
	ylab('Density')
ggsave(p,file=paste0('figures/model_results_dglm_sbeta_density_mash_',tolower(predictor.label),'.pdf'),width=7,height=5,useDingbats=FALSE)

p = ggplot(subset(b.mash,qval < fsr.cutoff),aes(beta,color=Var2)) +
	scale_color_manual(name='Region',values=region.colors) +
	geom_density() +
	theme_classic(base_size=12) +
	guides(color = guide_legend(ncol = 2)) +
	xlab(expression(italic(beta))) +
	ylab('Density')
ggsave(p,file=paste0('figures/model_results_dglm_beta_density_mash_',tolower(predictor.label),'.pdf'),width=7,height=3,useDingbats=FALSE)

p = ggplot(melt(tapply(b.mash$qval,b.mash$Var2,function(x) sum(x < fsr.cutoff,na.rm=TRUE))),aes(Var1,value,fill=Var1)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Region',values=region.colors) +
	theme_classic(base_size=12) +
	theme(legend.position='none',axis.text.x=element_text(angle=-45,hjust=0,vjust=1)) +
	xlab('Regions') +
	ylab('Number of genes')
ggsave(p,file=paste0('figures/model_results_dglm_beta_count_mash_',tolower(predictor.label),'.pdf'),width=7,height=3,useDingbats=FALSE)

b.dglm$qval.signed = with(b.dglm,ifelse(beta>0,1,-1) * qval)
b.mash$qval.signed = with(b.mash,ifelse(beta>0,1,-1) * qval)

# Put together longform data frame with gene counts (up vs. Decreased variance, split by method)
model.counts.combined = rbind(
	within(melt(tapply(b.dglm$qval.signed,b.dglm$Var2,function(x) -sum(abs(x) < fsr.cutoff & x < 0,na.rm=TRUE))), {method='DGLM'; direction='down'} ),
	within(melt(tapply(b.dglm$qval.signed,b.dglm$Var2,function(x) sum(abs(x) < fsr.cutoff & x >= 0,na.rm=TRUE))), {method='DGLM'; direction='up'} ),
	within(melt(tapply(b.mash$qval.signed,b.mash$Var2,function(x) -sum(abs(x) < fsr.cutoff & x < 0,na.rm=TRUE))), {method='MASH'; direction='down'} ),
	within(melt(tapply(b.mash$qval.signed,b.mash$Var2,function(x) sum(abs(x) < fsr.cutoff & x >= 0,na.rm=TRUE))), {method='MASH'; direction='up'} )
)

ylimit = ceiling(with(model.counts.combined,max(abs(value)))/500) * 500

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
	theme_article(base_size=12) +
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
ggsave(p,file=paste0('figures/model_results_dglm_beta_count_comparison_',tolower(predictor.label),'.pdf'),width=7,height=5,useDingbats=FALSE)


ylimit = ceiling(with(subset(model.counts.combined,method=='MASH'),max(abs(value)))/500) * 500

p = ggplot(subset(model.counts.combined,method=='MASH'),aes(Var1,value,fill=Var1,alpha=direction)) +
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
ggsave(p,file=paste0('figures/model_results_dglm_beta_count_mash_',tolower(predictor.label),'.pdf'),width=7,height=5,useDingbats=FALSE)


# pair.share = get_pairwise_sharing(mash.results, factor = 0.5)

# Metadata counting of significant effects

b.mash.split.genes = split(b.mash,b.mash$Var1)

# # Tally up sig region-combos, as well as up-regulated and down-regulated
# region.combinations = table(unlist(lapply(b.mash.split.genes,function(x) {
# 	paste(as.character(subset(x,qval < fsr.cutoff & beta != 0)$Var2),collapse='-')
# })))
# 
# # Get rid of blanks
# region.combinations = region.combinations[as.logical(nchar(names(region.combinations)))]
# 
# # Sort in same order (by total)
# region.combinations = region.combinations[order(region.combinations,decreasing=TRUE)]
# 
# region.combinations.inc = region.combinations.dec = integer(length(region.combinations))
# names(region.combinations.inc) = names(region.combinations.dec) = names(region.combinations)

# Ensure all have the same order
region.combinations.inc = table(unlist(lapply(b.mash.split.genes,function(x) {
	paste(as.character(subset(x,qval < fsr.cutoff & beta > 0)$Var2),collapse='-')
})))

region.combinations.dec = table(unlist(lapply(b.mash.split.genes,function(x) {
	paste(subset(x,qval < fsr.cutoff & beta < 0)$Var2,collapse='-')
})))

region.combinations.labels = union(names(region.combinations.inc),names(region.combinations.dec))
region.combinations.labels = region.combinations.labels[region.combinations.labels != '']

region.combinations.inc = region.combinations.inc[region.combinations.labels]
region.combinations.dec = region.combinations.dec[region.combinations.labels]
names(region.combinations.inc) = names(region.combinations.dec) = region.combinations.labels
region.combinations.inc[is.na(region.combinations.inc)] = 0
region.combinations.dec[is.na(region.combinations.dec)] = 0

region.combinations = region.combinations.inc + region.combinations.dec

region.combinations.order = names(region.combinations[order(region.combinations,decreasing=TRUE)])

region.combinations.inc = region.combinations.inc[region.combinations.order]
region.combinations.dec = region.combinations.dec[region.combinations.order]

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

fraction.shared.cutoff=1/3

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
			share_region = n >= fraction.shared.cutoff * length(keep.genes),
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
			share_region = n >= fraction.shared.cutoff * length(keep.genes),
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
			share_region = n >= fraction.shared.cutoff * length(keep.genes),
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
			share_region = n >= fraction.shared.cutoff * length(keep.genes),
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
	ylab('Increased variance') +
	theme_classic(base_size=base.size) +
	theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
p3 = ggplot(subset(meta.combinations.results.plot,chart=='count_decrease')) +
	geom_bar(aes(combination,value),stat='identity',fill=base.color,width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(trans='reverse',expand=c(0,0),limits=c(ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100,0)) +
	coord_flip() +
	ylab('Decreased variance') +
	theme_classic(base_size=base.size) +
	theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())

pdf(file=paste0('figures/model_results_dglm_beta_count_comparison_',tolower(predictor.label),'_by_meta.pdf'),useDingbats=FALSE,height=7,width=11)
	ggarrange(p3,p1,p2,ncol=3,nrow=1,widths=c(1,1,1),newpage=FALSE)
dev.off()

p4 = ggplot(subset(meta.combinations.results.plot,chart=='count_increase')) +
	geom_bar(aes(combination,value,fill=share_region),stat='identity',width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0),limits=c(0,ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100)) +
	scale_fill_manual(values=c(fade.color,highlight.color)) +
	coord_flip() +
	ylab('Increased variance') +
	theme_classic(base_size=base.size) +
	theme(legend.position='none',axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
p5 = ggplot(subset(meta.combinations.results.plot,chart=='count_decrease')) +
	geom_bar(aes(combination,value,fill=share_region),stat='identity',width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(trans='reverse',expand=c(0,0),limits=c(ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100,0)) +
	scale_fill_manual(values=c(fade.color,highlight.color)) +
	coord_flip() +
	ylab('Decreased variance') +
	theme_classic(base_size=base.size) +
	theme(legend.position='none',axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())

pdf(file=paste0('figures/model_results_dglm_beta_count_comparison_',tolower(predictor.label),'_by_meta_shared.pdf'),useDingbats=FALSE,height=7,width=11)
	ggarrange(p5,p1,p4,ncol=3,nrow=1,widths=c(1,1,1),newpage=FALSE)
dev.off()

p6 = ggplot(subset(meta.combinations.results.plot,chart=='count_increase')) +
	geom_bar(aes(combination,value,fill=n_regions==1),stat='identity',width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0),limits=c(0,ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100)) +
	scale_fill_manual(values=c(fade.color,highlight.color)) +
	coord_flip() +
	ylab('Increased variance') +
	theme_classic(base_size=base.size) +
	theme(legend.position='none',axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
p7 = ggplot(subset(meta.combinations.results.plot,chart=='count_decrease')) +
	geom_bar(aes(combination,value,fill=n_regions==1),stat='identity',width=width.of.bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(trans='reverse',expand=c(0,0),limits=c(ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100,0)) +
	scale_fill_manual(values=c(fade.color,highlight.color)) +
	coord_flip() +
	ylab('Decreased variance') +
	theme_classic(base_size=base.size) +
	theme(legend.position='none',axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())

pdf(file=paste0('figures/model_results_dglm_beta_count_comparison_',tolower(predictor.label),'_by_meta_single.pdf'),useDingbats=FALSE,height=7,width=11)
	ggarrange(p7,p1,p6,ncol=3,nrow=1,widths=c(1,1,1),newpage=FALSE)
dev.off()


# Inputs
# mash.lfsr: matrix with regions as columns, genes as rows, and mashr LFSRs as values
# mash.beta: matrix with regions as columns, genes as rows, and mashr betas as values
# region.levels: vector with regions in preferred order
# region.colors: vector with color choices in hex RGB.
threshold.range = seq(0.01,0.2,0.01)
# Loop through a series of cutoffs
which.region.tally = do.call(rbind,lapply(threshold.range,function(threshold) {
	# raw counts
	# unique counts
	# Create matrices on only those genes that are significant in that direction in only one region
	this.inc.lfsr = mash.lfsr[rowSums(mash.lfsr < threshold & mash.beta > 0) == 1,]
	this.inc.beta = mash.beta[rowSums(mash.lfsr < threshold & mash.beta > 0) == 1,]
	this.dec.lfsr = mash.lfsr[rowSums(mash.lfsr < threshold & mash.beta < 0) == 1,]
	this.dec.beta = mash.beta[rowSums(mash.lfsr < threshold & mash.beta < 0) == 1,]
	out = rbind(
		# Counts/proportions of total significant genes with positive betas in each region
		data.frame(
			threshold,
			region = colnames(mash.lfsr),
			count = as.integer(colSums(mash.lfsr < threshold & mash.beta > 0)),
			proportion = as.numeric(
				colSums(mash.lfsr < threshold & mash.beta > 0) /
					sum(mash.lfsr < threshold & mash.beta > 0)
			),
			direction = 'Increase',
			type = 'total'
		),
		# Counts/proportions of total significant genes with negative betas in each region
		data.frame(
			threshold,
			region = colnames(mash.lfsr),
			count = as.integer(colSums(mash.lfsr < threshold & mash.beta < 0)),
			proportion = as.numeric(
				colSums(mash.lfsr < threshold & mash.beta < 0) /
					sum(mash.lfsr < threshold & mash.beta < 0)
			),
			direction = 'Decrease',
			type = 'total'
		),
		# Counts/proportions of unique significant genes with positive betas in each region
		data.frame(
			threshold,
			region = colnames(mash.lfsr),
			count = as.integer(colSums(this.inc.lfsr < threshold & this.inc.beta > 0)),
			proportion = as.numeric(
				colSums(this.inc.lfsr < threshold & this.inc.beta > 0) /
					sum(this.inc.lfsr < threshold & this.inc.beta > 0)
			),
			direction = 'Increase',
			type = 'unique'
		),
		# Counts/proportions of unique significant genes with negative betas in each region
		data.frame(
			threshold,
			region = colnames(mash.lfsr),
			count = as.integer(colSums(this.inc.lfsr < threshold & this.inc.beta < 0)),
			proportion = as.numeric(
				colSums(this.inc.lfsr < threshold & this.inc.beta < 0) /
					sum(this.inc.lfsr < threshold & this.inc.beta < 0)
			),
			direction = 'Decrease',
			type = 'unique'
		)
	)
	out$direction = factor(out$direction,levels=c('Increase','Decrease'))
	out$region = factor(out$region,levels=region.levels)
	out
}))

# The y-axis below can be switched between proportion and count below. The rank order will stay the same
# Plot for total number of significant genes per region
p = ggplot(subset(which.region.tally,type=='total'),aes(threshold,count,color=region)) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (count)')
ggsave(p,file=paste0('figures/model_results_dglm_sig_counts_total_',tolower(predictor.label),'_region_tally.pdf'),useDingbats=FALSE,height=5)

# Plot of unique significant genes per region
p = ggplot(subset(which.region.tally,type=='unique'),aes(threshold,count,color=region)) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (count)')
ggsave(p,file=paste0('figures/model_results_dglm_sig_counts_unique_',tolower(predictor.label),'_region_tally.pdf'),useDingbats=FALSE,height=5)

# Plot for total number of significant genes per region
p = ggplot(subset(which.region.tally,type=='total'),aes(threshold,proportion,color=region)) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (proportion)')
ggsave(p,file=paste0('figures/model_results_dglm_sig_proportions_total_',tolower(predictor.label),'_region_tally.pdf'),useDingbats=FALSE,height=5)

# Plot of unique significant genes per region
p = ggplot(subset(which.region.tally,type=='unique'),aes(threshold,proportion,color=region)) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (proportion)')
ggsave(p,file=paste0('figures/model_results_dglm_sig_proportions_unique_',tolower(predictor.label),'_region_tally.pdf'),useDingbats=FALSE,height=5)

# Number of genes in n regions
n.region.tally = do.call(rbind,lapply(threshold.range,function(threshold) {
	# raw counts
	table(rowSums(mash.lfsr < threshold & mash.beta > 0))
	table(rowSums(mash.lfsr < threshold & mash.beta < 0))
	out = rbind(
		data.frame(
			threshold,
			n.regions = as.integer(names(table(rowSums(mash.lfsr < threshold & mash.beta > 0)))),
			count = as.integer(table(rowSums(mash.lfsr < threshold & mash.beta > 0))),
			proportion = as.numeric(
				table(rowSums(mash.lfsr < threshold & mash.beta > 0)) / nrow(mash.lfsr)
			),
			direction = 'Increase'
		),
		data.frame(
			threshold,
			n.regions = as.integer(names(table(rowSums(mash.lfsr < threshold & mash.beta < 0)))),
			count = as.integer(table(rowSums(mash.lfsr < threshold & mash.beta < 0))),
			proportion = as.numeric(
				table(rowSums(mash.lfsr < threshold & mash.beta < 0)) / nrow(mash.lfsr)
			),
			direction = 'Decrease'
		)
	)
	out$direction = factor(out$direction,levels=c('Increase','Decrease'))
	out
}))
p = ggplot(droplevels(subset(n.region.tally,n.regions>0)),aes(threshold,count,color=factor(n.regions))) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
#	scale_color_manual(values=c('#000000',region.colors)) +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (count)')
ggsave(p,file=paste0('figures/model_results_dglm_sig_counts_',tolower(predictor.label),'_n_regions_tally.pdf'),useDingbats=FALSE,height=5)

p = ggplot(droplevels(subset(n.region.tally,n.regions>0)),aes(threshold,proportion,color=factor(n.regions))) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
#	scale_color_manual(values=c('#000000',region.colors)) +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (proportion)')
ggsave(p,file=paste0('figures/model_results_dglm_sig_proportions_',tolower(predictor.label),'_n_regions_tally.pdf'),useDingbats=FALSE,height=5)







plot.gene = function(ensembl_gene_id,scale=TRUE) {
	do.call(rbind,lapply(region.levels,function(x) {
		out = subset(meta,Region %in% x)
		if (scale) {
			out$gene = scale(matrix(e.keep[ensembl_gene_id,rownames(out)],ncol=1))
		} else {
			out$gene = matrix(e.keep[ensembl_gene_id,rownames(out)],ncol=1)
		}
		out
	}))
}

mash.sig = unlist(lapply(1:nrow(mash.lfsr),function(i) sum(mash.beta[i,] > 0 & mash.lfsr[i,] < fsr.cutoff) >= 5))

best.median = names(which.max(apply(dglm.sbet[mash.sig,],1,median)))
best.mean = names(which.max(apply(dglm.sbet[mash.sig,],1,mean)))
best.max = names(which.max(apply(dglm.sbet[mash.sig,],1,max)))

moo = data.frame(
	count = sapply(rownames(dglm.beta[mash.sig,]),function(x) sum(dglm.qval[x,] < fsr.cutoff & dglm.beta[x,] > 0)),
	median = sapply(rownames(dglm.beta[mash.sig,]),function(x) median(dglm.sbet[x,]))
)
best.count = rownames(moo[with(moo,order(count,median,decreasing=TRUE)),])[1]

moo = data.frame(
	count = sapply(rownames(dglm.beta[mash.sig,]),function(x) sum(dglm.sbet[x,] > 0)),
	median = sapply(rownames(dglm.beta[mash.sig,]),function(x) median(dglm.sbet[x,]))
)
best.mode = rownames(moo[with(moo,order(count,median,decreasing=TRUE)),])[1]

best.mash = names(which.max(rowMeans(mash.sbet[mash.sig,])))

example.genes = c(best.median,best.mean,best.max,best.count,best.mode)

for (i in example.genes) {
	foo = plot.gene(i,scale=TRUE)
	levels(foo$Region) = gsub('ACCg','ACC',levels(foo$Region))
	
	foo$predictor_group = factor(foo[[predictor]] <= min(foo[[predictor]]) + diff(range(foo[[predictor]]))/2,levels=c('TRUE','FALSE'))
	levels(foo$predictor_group) = c(g1.label,g2.label)

	p = ggplot(foo,aes_string(predictor,'gene')) +
		geom_smooth(aes_string(fill='Region'),method='lm',size=0.5,linetype=2,alpha=0.75,color='#ffffff') +
		geom_point(color='#000000',size=1,shape=21) +
		facet_wrap(~Region,nrow=region.rows) +
		scale_fill_manual(values=region.colors) +
#		scale_color_manual(values=region.colors) +
		scale_y_continuous(limits=c(-max(abs(foo$gene)),max(abs(foo$gene)))) +
		xlab(predictor.label) + ylab(expression(log[2]('Expression'))) +
		theme_classic() +
		theme(legend.position='none',axis.text.y=element_blank(),axis.ticks.y=element_blank())
	ggsave(p,file=paste0('figures/dglm_variance_example_',predictor,'_by_gene_',i,'.pdf'),useDingbats=FALSE,height=5)

	foo.bar = do.call(rbind,lapply(split(foo,with(foo,list(predictor_group,Region))),function(x) data.frame( predictor_group = factor(unique(x$predictor_group),levels=c(g1.label,g2.label)), Region = factor(unique(x$Region),levels=levels(x$Region)), gene = mean(x$gene), e.sd = sd(x$gene), e.var = with(x,mean((gene - mean(gene))^2)), e.ub = mean(x$gene) + with(x,mean((gene - mean(gene))^2)), e.lb = mean(x$gene) - with(x,mean((gene - mean(gene))^2)))))

	p = ggplot(foo,aes_string('predictor_group','gene',color='Region')) +
		geom_point(aes_string('predictor_group','gene',color='Region'),data=foo.bar,size=2) +
		geom_errorbar(aes_string('predictor_group',ymin='e.lb',ymax='e.ub',color='Region'),data=foo.bar,size=0.5,width=0.25) +
		geom_jitter(width=0.25,height=0,size=0.5) +
		facet_wrap(~Region,nrow=region.rows) +
		scale_color_manual(values=region.colors) +
		scale_y_continuous(limits=c(-max(abs(foo$gene)),max(abs(foo$gene)))) +
		xlab(predictor.label) + ylab(expression(log[2]('Expression'))) +
		theme_classic() +
		theme(legend.position='none',axis.text.y=element_blank(),axis.ticks.y=element_blank())
	ggsave(p,file=paste0('figures/dglm_variance_example_',predictor,'_by_gene_',i,'_variance.pdf'),useDingbats=FALSE,height=5)

	p = ggplot(foo,aes_string('predictor_group','gene',color='Region')) +
		geom_point(aes_string('predictor_group','gene',color='Region'),data=foo.bar,size=2) +
		geom_errorbar(aes_string('predictor_group',ymin='e.lb',ymax='e.ub',color='Region'),data=foo.bar,size=0.5,width=0.25) +
		geom_jitter(width=0.25,height=0,size=0.5) +
		facet_wrap(~Region,nrow=region.rows) +
		scale_color_manual(values=region.colors) +
		scale_y_continuous(limits=c(-max(abs(foo$gene)),max(abs(foo$gene))),breaks=seq(-4,4,2)) +
		xlab(paste0(predictor.label,' (',tolower(predictor.units),')')) + ylab(expression(italic(Z)*'-score')) +
		theme_article(base_size=16) +
		theme(legend.position='none')
	ggsave(p,file=paste0('figures/dglm_variance_example_',predictor,'_by_gene_',i,'_variance_article.pdf'),useDingbats=FALSE,height=5)

}

go.enrichment.results = readRDS('checkpoints/dglm_topgo_results.rds')

library(ggrastr)
for (i in c('FET','KST')) {

	go.enrichment.results = subset(go.enrichment.results.all,grepl(paste0('^',toupper(substr(i,1,2))),test))

	p = ggplot(
		do.call(rbind,lapply(split(subset(go.enrichment.results,region=='all'),subset(go.enrichment.results,region=='all')$direction),function(x) {
			within(data.frame(unique(x[c('direction','region')]),expected=-log10(seq(1/nrow(x),1,1/nrow(x))),observed=-log10(quantile(x$pval,seq(1/nrow(x),1,1/nrow(x)),na.rm=TRUE))),{
				direction=factor(direction,levels=c('increase','decrease'),labels=c('Increase','Decrease'))
			})
		})),
		aes(expected,observed,color=direction)) +
#		geom_smooth(method=lm,se=FALSE,color='black',size=0.2,linetype=2) +
		geom_point_rast(size=0.5) +
		geom_abline(slope=1,col='black',size=0.2) +
		scale_x_continuous(breaks=seq(0,ceiling(-log10(min(1/(max(table(subset(go.enrichment.results,region=='all')$direction,subset(go.enrichment.results,region=='all')$region))-1)))),1)) +
		scale_y_continuous(breaks=seq(0,ceiling(max(-log10(subset(go.enrichment.results,region=='all')$pval))),1)) +
		scale_color_manual(name='Direction',values=c('#386cb0','#f0027f')) +
	#	coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
		xlab(expression(-log[10] * ('Expected'~italic(p)))) +
		ylab(expression(-log[10] * ('Observed'~italic(p)))) +
		theme_classic() +
		theme() +
		guides(color=guide_legend(override.aes=list(size=3)))
	ggsave(p,file=paste0('figures/dglm_variance_pval_enrichment_topgo_',i,'_',tolower(predictor.label),'_qqplot_union.pdf'),useDingbats=FALSE)

	p = ggplot(
		subset(do.call(rbind,lapply(split(subset(go.enrichment.results,region=='all'),subset(go.enrichment.results,region=='all')$direction),function(x) {
			within(data.frame(unique(x[c('direction','region')]),expected=-log10(seq(1/nrow(x),1,1/nrow(x))),observed=-log10(quantile(x$pval,seq(1/nrow(x),1,1/nrow(x)),na.rm=TRUE))),{
				direction=factor(direction,levels=c('increase','decrease'),labels=c('Increase','Decrease'))
			})
		})),direction=='Increase'),
		aes(expected,observed,color=direction)) +
#		geom_smooth(method=lm,se=FALSE,color='black',size=0.2,linetype=2) +
		geom_point_rast(size=0.5) +
		geom_abline(slope=1,col='black',size=0.2) +
		scale_x_continuous(breaks=seq(0,ceiling(-log10(min(1/(max(table(subset(go.enrichment.results,region=='all')$direction,subset(go.enrichment.results,region=='all')$region))-1)))),1)) +
		scale_y_continuous(breaks=seq(0,ceiling(max(-log10(subset(go.enrichment.results,region=='all')$pval))),1)) +
		scale_color_manual(name='Direction',values=c('#386cb0','#f0027f')) +
	#	coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
		xlab(expression(-log[10] * ('Expected'~italic(p)))) +
		ylab(expression(-log[10] * ('Observed'~italic(p)))) +
		theme_classic(base_size=24) +
		theme(legend.position='none') +
		guides(color=guide_legend(override.aes=list(size=3)))
	ggsave(p,file=paste0('figures/dglm_variance_pval_enrichment_topgo_',i,'_',tolower(predictor.label),'_qqplot_union_presentation.pdf'),useDingbats=FALSE)
}


save(list=ls(),file='checkpoints/checkpoint_model_results_dglm_visualization.RData')

