#!/usr/bin/env Rscript

source('scripts/_include_options.R')

e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')

mash.results = readRDS('checkpoints/mashr_results.rds')
emma.results = readRDS('checkpoints/emma_results.rds')
keep.genes = readRDS('checkpoints/keep_genes.rds')

emma.beta = emma.results[,paste('beta',predictor,sep='.'),]
emma.pval = emma.results[,paste('pval',predictor,sep='.'),]
emma.qval = apply(emma.pval,2,function(x) p.adjust(x,'fdr'))

library(mashr)
mash.beta = get_pm(mash.results)
mash.lfsr = get_lfsr(mash.results)

mash.sbet = get_pm(mash.results) / get_psd(mash.results)

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
	xlab(expression(italic(P)~'value')) + ylab('Count') + theme_classic() +
	theme(legend.position='none')
ggsave(p,file=paste0('figures/model_results_pval_enrichment_emma_',tolower(predictor.label),'_histogram.pdf'),useDingbats=FALSE,height=5)

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
# ggsave(p,file=paste0('figures/model_results_pval_enrichment_emma_',tolower(predictor.label),'_qqplot.pdf'),useDingbats=FALSE,height=5)

library(ggrastr)
library(egg)
# QQ plot
p = ggplot(
		within(do.call(rbind,lapply(split(melt(emma.pval),melt(emma.pval)$Var2),function(x) {
			data.frame(
				region=unique(x$Var2),
				expected=-log10(seq(1/sum(complete.cases(x)),1,1/sum(complete.cases(x)))),
				observed=-log10(quantile(x$value,seq(1/sum(complete.cases(x)),1,1/sum(complete.cases(x))),na.rm=TRUE))
			)
		})),{levels(region) = gsub('ACCg','ACC',levels(region))}),
		aes(expected,observed,color=region)
	) +
	geom_smooth(method=lm,se=FALSE,color='black',size=0.2,linetype=2) +
	geom_point_rast(size=0.5) +
	geom_abline(slope=1,col='black',size=0.2) +
	facet_wrap(~region,nrow=region.rows) +
#	coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
	scale_color_manual(values=region.colors) +
 	scale_x_continuous(breaks=seq(0,ceiling(-log10(min(1/(max(table(melt(emma.pval)$Var2),na.rm=TRUE)-1)))),1)) +
 	scale_y_continuous(breaks=seq(0,ceiling(max(-log10(melt(emma.pval)$value),na.rm=TRUE)),1)) +
	xlab(expression(-log[10] * ('Expected'~italic(p)))) +
	ylab(expression(-log[10] * ('Observed'~italic(p)))) +
	theme_classic() +
	theme(legend.position='none')
ggsave(p,file=paste0('figures/model_results_pval_enrichment_emma_',tolower(predictor.label),'_qqplot.pdf'),useDingbats=FALSE,height=5)

p = ggplot(
		within(do.call(rbind,lapply(split(melt(emma.pval),melt(emma.pval)$Var2),function(x) {
			data.frame(
				region=unique(x$Var2),
				expected=-log10(seq(1/sum(complete.cases(x)),1,1/sum(complete.cases(x)))),
				observed=-log10(quantile(x$value,seq(1/sum(complete.cases(x)),1,1/sum(complete.cases(x))),na.rm=TRUE))
			)
		})),{levels(region) = gsub('ACCg','ACC',levels(region))}),
		aes(expected,observed,color=region)
	) +
#	geom_smooth(method=lm,se=FALSE,color='black',size=0.2,linetype=2) +
	geom_point_rast(size=0.5) +
	geom_abline(slope=1,col='black',size=0.2) +
	facet_wrap(~region,nrow=region.rows) +
#	coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
	scale_color_manual(values=region.colors) +
 	scale_x_continuous(breaks=seq(0,ceiling(-log10(min(1/(max(table(melt(emma.pval)$Var2),na.rm=TRUE)-1)))),1)) +
 	scale_y_continuous(breaks=seq(0,ceiling(max(-log10(melt(emma.pval)$value),na.rm=TRUE)),1)) +
	xlab(expression(-log[10] * ('Expected'~italic(p)))) +
	ylab(expression(-log[10] * ('Observed'~italic(p)))) +
	theme_article() +
	theme(legend.position='none')
ggsave(p,file=paste0('figures/model_results_pval_enrichment_emma_',tolower(predictor.label),'_qqplot_article.pdf'),useDingbats=FALSE,height=5)


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
ggsave(p,file=paste0('figures/model_results_beta_density_emma_',tolower(predictor.label),'.pdf'),width=7,height=3,useDingbats=FALSE)

p = ggplot(melt(tapply(b.emma$qval,b.emma$Var2,function(x) sum(x < fsr.cutoff,na.rm=TRUE))),aes(Var1,value,fill=Var1)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Region',values=region.colors) +
	theme_classic(base_size=12) +
	theme(legend.position='none',axis.text.x=element_text(angle=-45,hjust=0,vjust=1)) +
	xlab('Regions') +
	ylab('Number of genes')
ggsave(p,file=paste0('figures/model_results_beta_count_emma_',tolower(predictor.label),'.pdf'),width=7,height=3,useDingbats=FALSE)

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
ggsave(p,file=paste0('figures/model_results_beta_density_mash_',tolower(predictor.label),'.pdf'),width=7,height=3,useDingbats=FALSE)

p = ggplot(melt(tapply(b.mash$qval,b.mash$Var2,function(x) sum(x < fsr.cutoff,na.rm=TRUE))),aes(Var1,value,fill=Var1)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Region',values=region.colors) +
	theme_classic(base_size=12) +
	theme(legend.position='none',axis.text.x=element_text(angle=-45,hjust=0,vjust=1)) +
	xlab('Regions') +
	ylab('Number of genes')
ggsave(p,file=paste0('figures/model_results_beta_count_mash_',tolower(predictor.label),'.pdf'),width=7,height=3,useDingbats=FALSE)

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
	facet_wrap(~method,nrow=2) +
	theme(legend.position='none') +
	ylab('Number of genes')
ggsave(p,file=paste0('figures/model_results_beta_count_comparison_',tolower(predictor.label),'.pdf'),width=7,height=5,useDingbats=FALSE)

levels(model.counts.combined$Var1) = gsub('ACCg','ACC',levels(model.counts.combined$Var1))
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
ggsave(p,file=paste0('figures/model_results_beta_count_comparison_',tolower(predictor.label),'_article.pdf'),width=7,height=5,useDingbats=FALSE)


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
ggsave(p,file=paste0('figures/model_results_beta_count_mash_',tolower(predictor.label),'.pdf'),width=7,height=5,useDingbats=FALSE)


# pair.share = get_pairwise_sharing(mash.results, factor = 0.5)

# Metadata counting of significant effects

# First, simply count how many genes are significant and in one direction in only one region
unique.genes = unlist(mclapply(rownames(mash.beta),function(i) {
	if ((sum(mash.lfsr[i,] < fsr.cutoff & mash.beta[i,] > 0) == 1) || (sum(mash.lfsr[i,] < fsr.cutoff & mash.beta[i,] < 0) == 1)) {
		i
	} else {
		character(0)
	}
},mc.cores=n.cores))

print(length(unique.genes))

# First, simply count how many genes are significant and in one direction in only one region
conflicting.genes = unlist(mclapply(rownames(mash.beta),function(i) {
	if ((sum(mash.lfsr[i,] < fsr.cutoff & mash.beta[i,] > 0) == 1) && (sum(mash.lfsr[i,] < fsr.cutoff & mash.beta[i,] < 0) == 1)) {
		i
	} else {
		character(0)
	}
},mc.cores=n.cores))

print(length(conflicting.genes))

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

meta.combinations.results.plot = subset(region.combinations.results,combination <= max.plot)

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

pdf(file=paste0('figures/model_results_beta_count_comparison_',tolower(predictor.label),'_by_meta.pdf'),useDingbats=FALSE,height=7,width=11)
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

pdf(file=paste0('figures/model_results_beta_count_comparison_',tolower(predictor.label),'_by_meta_shared.pdf'),useDingbats=FALSE,height=7,width=11)
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

pdf(file=paste0('figures/model_results_beta_count_comparison_',tolower(predictor.label),'_by_meta_single.pdf'),useDingbats=FALSE,height=7,width=11)
	ggarrange(p7,p1,p6,ncol=3,nrow=1,widths=c(1,1,1),newpage=FALSE)
dev.off()



pair.effects = data.frame(reshape2::melt(mash.lfsr),beta=reshape2::melt(mash.beta)$value)
names(pair.effects) = c('ensembl_gene_id','region','lfsr','beta')
levels(pair.effects$region) = gsub('ACCg','ACC',levels(pair.effects$region))

#region.pairs = subset(expand.grid(levels(pair.effects$region),levels(pair.effects$region)),Var1 != Var2)

region.pairs = gsub('ACCg','ACC',matrix(region.levels[t(combn(length(region.levels),2))],ncol=2))

mash.genes = rownames(mash.beta)

library(parallel)

gene.share = do.call(rbind,mclapply(1:nrow(region.pairs),function(i) {
	this.pair = region.pairs[i,]
	
	e1 = subset(pair.effects,region == this.pair[1])
	e2 = subset(pair.effects,region == this.pair[2])
	rownames(e1) = e1$ensembl_gene_id
	rownames(e2) = e2$ensembl_gene_id
	e1 = e1[mash.genes,]
	e2 = e2[mash.genes,]
	if (!identical(rownames(e1),rownames(e2))) stop('error')

	colnames(e1)[2:ncol(e1)] = paste0(colnames(e1)[2:ncol(e1)],'.1')
	colnames(e2)[2:ncol(e2)] = paste0(colnames(e2)[2:ncol(e2)],'.2')

	data.frame(e1,e2[2:ncol(e2)])

},mc.cores=n.cores))

gene.share = within(gene.share,{
	pass.magn = ifelse(beta.1 * beta.2 > 0,apply(cbind(beta.1,beta.2),1,function(x) max(abs(x)) < min(abs(x)) * 2),FALSE)
	pass.sign = beta.1 * beta.2 > 0
	pass.lfsr = lfsr.1 < fsr.cutoff | lfsr.2 < fsr.cutoff
})

region.share.mag = do.call(rbind,mclapply(1:nrow(region.pairs),function(i) {
	this.pair = region.pairs[i,]
	
	this = subset(gene.share,region.1 == this.pair[1] & region.2 == this.pair[2] & pass.lfsr)
	out = unique(subset(this,select=c('region.1','region.2')))
	out$shared = mean(this$pass.magn)
	rownames(out) = NULL
	out
},mc.cores=n.cores))

region.share.sgn = do.call(rbind,mclapply(1:nrow(region.pairs),function(i) {
	this.pair = region.pairs[i,]
	
	this = subset(gene.share,region.1 == this.pair[1] & region.2 == this.pair[2] & pass.lfsr)
	out = unique(subset(this,select=c('region.1','region.2')))
	out$shared = mean(this$pass.sign)
	rownames(out) = NULL
	out
},mc.cores=n.cores))

p = ggplot() +
	geom_text(
		data=data.frame(label=gsub('ACCg','ACC',region.levels),x=(1:15)-0.15,y=(15:1)-0.45),
		aes(x,y,label=label),
		hjust=1, vjust=0, angle=45
	) +
	geom_point(
		data=data.frame(label=factor(gsub('ACCg','ACC',region.levels),levels=gsub('ACCg','ACC',region.levels)),x=(1:15),y=(15:1)),
		aes(x,y,color=label,shape=label),
		size=2,
		show.legend=FALSE
	) +
	geom_point(
		data=subset(data.frame(label=factor(gsub('ACCg','ACC',region.levels),levels=gsub('ACCg','ACC',region.levels)),x=(1:15),y=16),x>1),
		aes(x,y,color=label,shape=label),
		size=2,
		show.legend=FALSE
	) +
	geom_point(
		data=subset(data.frame(label=factor(gsub('ACCg','ACC',region.levels),levels=gsub('ACCg','ACC',region.levels)),x=16,y=(15:1)),y>1),
		aes(x,y,color=label,shape=label),
		size=2,
		show.legend=FALSE
	) +
	geom_tile(data=mutate(region.share.mag,region.1=as.integer(factor(region.1,levels=gsub('ACCg','ACC',region.levels[length(region.levels):1]))),region.2=as.integer(factor(region.2,levels=gsub('ACCg','ACC',region.levels)))),aes(x=region.2,y=region.1,fill=shared)) +
	scale_fill_viridis(option='viridis',limits=c(0,1),breaks=seq(0,1,0.25)) +
	scale_color_manual(values=region.colors) +
	scale_shape_manual(values=region.shapes) +
	scale_x_continuous(limits=c(0,length(region.levels)+1)) +
	scale_y_continuous(limits=c(0,length(region.levels)+1)) +
	coord_fixed() +
	guides(fill=guide_colorbar(barheight=15)) +
	theme_classic(base_size=16) +
	theme(
		axis.text=element_blank(),
		axis.ticks=element_blank(),
		axis.line=element_blank(),
		axis.title=element_blank(),
		legend.title=element_blank()
	)
ggsave(p,file=paste0('figures/model_results_shared_magnitude_',tolower(predictor.label),'.pdf'),width=7,height=5,useDingbats=FALSE)


p = ggplot() +
	geom_text(
		data=data.frame(label=gsub('ACCg','ACC',region.levels),x=(1:15)-0.15,y=(15:1)-0.45),
		aes(x,y,label=label),
		hjust=1, vjust=0, angle=45
	) +
	geom_point(
		data=data.frame(label=factor(gsub('ACCg','ACC',region.levels),levels=gsub('ACCg','ACC',region.levels)),x=(1:15),y=(15:1)),
		aes(x,y,color=label,shape=label),
		size=2,
		show.legend=FALSE
	) +
	geom_point(
		data=subset(data.frame(label=factor(gsub('ACCg','ACC',region.levels),levels=gsub('ACCg','ACC',region.levels)),x=(1:15),y=16),x>1),
		aes(x,y,color=label,shape=label),
		size=2,
		show.legend=FALSE
	) +
	geom_point(
		data=subset(data.frame(label=factor(gsub('ACCg','ACC',region.levels),levels=gsub('ACCg','ACC',region.levels)),x=16,y=(15:1)),y>1),
		aes(x,y,color=label,shape=label),
		size=2,
		show.legend=FALSE
	) +
	geom_tile(data=mutate(region.share.sgn,region.1=as.integer(factor(region.1,levels=gsub('ACCg','ACC',region.levels[length(region.levels):1]))),region.2=as.integer(factor(region.2,levels=gsub('ACCg','ACC',region.levels)))),aes(x=region.2,y=region.1,fill=shared)) +
	scale_fill_viridis(option='viridis',limits=c(0,1),breaks=seq(0,1,0.25)) +
	scale_color_manual(values=region.colors) +
	scale_shape_manual(values=region.shapes) +
	scale_x_continuous(limits=c(0,length(region.levels)+1)) +
	scale_y_continuous(limits=c(0,length(region.levels)+1)) +
	coord_fixed() +
	guides(fill=guide_colorbar(barheight=15)) +
	theme_classic(base_size=16) +
	theme(
		axis.text=element_blank(),
		axis.ticks=element_blank(),
		axis.line=element_blank(),
		axis.title=element_blank(),
		legend.title=element_blank()
	)
ggsave(p,file=paste0('figures/model_results_shared_sign_',tolower(predictor.label),'.pdf'),width=7,height=5,useDingbats=FALSE)


# Repeat but calculate sharing across all regions
all.share = do.call(rbind,mclapply(1:length(mash.genes),function(i) {
	these.lfsr = mash.lfsr[i,]
	these.beta = mash.beta[i,]
	these.sbet = mash.sbet[i,]
	if (any(these.lfsr < fsr.cutoff)) {
		# j = which.max(abs(these.sbet))
		j = which.min(these.lfsr)
		if (length(j) > 1) stop('error')
		# Code: 0: directions not shared; 1: directions shared but not magnitude; 2: directions and magnitude shared
		# 
		out = ifelse(
			these.beta[j] * these.beta[-j] > 0,
			as.integer(abs(these.beta[-j]) > (abs(these.beta[j])/2))+1,
			0
		)
		data.frame(ensembl_gene_id=mash.genes[i],n.shared.direction=sum(out>0)+1,n.shared.magnitude=sum(out>1)+1)
	} else {
		data.frame(ensembl_gene_id=mash.genes[i],n.shared.direction=0,n.shared.magnitude=0)
	}
},mc.cores=n.cores))

all.share.proportions = do.call(rbind,lapply(names(all.share)[2:3],function(i) {
	deg = subset(all.share,n.shared.direction > 0)
	table(deg[[i]])/nrow(deg)
	out = reshape2::melt(table(deg[[i]])/nrow(deg))
	names(out) = c('n.regions','value')
	out$variable = factor(i,levels=c('n.shared.direction','n.shared.magnitude'),labels=c('Shared by direction','Shared by magnitude'))
	out
}))

p = ggplot(all.share.proportions,aes(n.regions,value)) +
	geom_bar(stat='identity') +
	facet_wrap(~variable,ncol=1,scales='fixed') +
	scale_x_continuous(breaks=1:15) +
	theme_classic(base_size=16) +
	theme(strip.background=element_blank(),strip.text=element_text(hjust=0)) +
	xlab('Number of regions') + ylab(paste0('Proportion of ',substr(tolower(predictor.label),1,1),'DEGs'))
ggsave(p,file=paste0('figures/model_results_n_shared_',tolower(predictor.label),'.pdf'),width=7,height=5,useDingbats=FALSE)

# Repeat but calculate sharing across set1 regions
set1.share = do.call(rbind,mclapply(1:length(mash.genes),function(i) {
	these.lfsr = mash.lfsr[i,region.set1]
	these.beta = mash.beta[i,region.set1]
	these.sbet = mash.sbet[i,region.set1]
	if (any(these.lfsr < fsr.cutoff)) {
		# j = which.max(abs(these.sbet))
		j = which.min(these.lfsr)
		if (length(j) > 1) stop('error')
		# Code: 0: directions not shared; 1: directions shared but not magnitude; 2: directions and magnitude shared
		# 
		out = ifelse(
			these.beta[j] * these.beta[-j] > 0,
			as.integer(abs(these.beta[-j]) > (abs(these.beta[j])/2))+1,
			0
		)
		data.frame(ensembl_gene_id=mash.genes[i],n.shared.direction=sum(out>0)+1,n.shared.magnitude=sum(out>1)+1)
	} else {
		data.frame(ensembl_gene_id=mash.genes[i],n.shared.direction=0,n.shared.magnitude=0)
	}
},mc.cores=n.cores))

set2.share = do.call(rbind,mclapply(1:length(mash.genes),function(i) {
	these.lfsr = mash.lfsr[i,region.set2]
	these.beta = mash.beta[i,region.set2]
	these.sbet = mash.sbet[i,region.set2]
	if (any(these.lfsr < fsr.cutoff)) {
		# j = which.max(abs(these.sbet))
		j = which.min(these.lfsr)
		if (length(j) > 1) stop('error')
		# Code: 0: directions not shared; 1: directions shared but not magnitude; 2: directions and magnitude shared
		# 
		out = ifelse(
			these.beta[j] * these.beta[-j] > 0,
			as.integer(abs(these.beta[-j]) > (abs(these.beta[j])/2))+1,
			0
		)
		data.frame(ensembl_gene_id=mash.genes[i],n.shared.direction=sum(out>0)+1,n.shared.magnitude=sum(out>1)+1)
	} else {
		data.frame(ensembl_gene_id=mash.genes[i],n.shared.direction=0,n.shared.magnitude=0)
	}
},mc.cores=n.cores))

set1.share.proportions = do.call(rbind,lapply(names(set1.share)[2:3],function(i) {
	deg = subset(set1.share,n.shared.direction > 0)
	table(deg[[i]])/nrow(deg)
	out = reshape2::melt(table(deg[[i]])/nrow(deg))
	names(out) = c('n.regions','value')
	out$variable = factor(i,levels=c('n.shared.direction','n.shared.magnitude'),labels=c('Shared by direction','Shared by magnitude'))
	out
}))

set2.share.proportions = do.call(rbind,lapply(names(set2.share)[2:3],function(i) {
	deg = subset(set2.share,n.shared.direction > 0)
	table(deg[[i]])/nrow(deg)
	out = reshape2::melt(table(deg[[i]])/nrow(deg))
	names(out) = c('n.regions','value')
	out$variable = factor(i,levels=c('n.shared.direction','n.shared.magnitude'),labels=c('Shared by direction','Shared by magnitude'))
	out
}))

all.set.share.proportions = rbind(
	mutate(all.share.proportions,set='all regions'),
	mutate(set1.share.proportions,set='cortical regions'),
	mutate(set2.share.proportions,set='subcortical regions')
)

p = ggplot(subset(all.set.share.proportions,variable=='Shared by direction'),aes(n.regions,value)) +
	geom_bar(stat='identity') +
	facet_wrap(~set,nrow=1,scales='free_x',drop=TRUE) +
	scale_x_continuous(breaks=1:15) +
	theme_classic(base_size=16) +
	theme(strip.background=element_blank(),strip.text=element_text(hjust=0)) +
	xlab('Number of regions') + ylab(paste0('Proportion of ',substr(tolower(predictor.label),1,1),'DEGs'))
ggsave(p,file=paste0('figures/model_results_n_shared_',tolower(predictor.label),'_direction_by_region_group.pdf'),width=15,height=5,useDingbats=FALSE)

p = ggplot(subset(all.set.share.proportions,variable=='Shared by magnitude'),aes(n.regions,value)) +
	geom_bar(stat='identity') +
	facet_wrap(~set,nrow=1,scales='free_x',drop=TRUE) +
	scale_x_continuous(breaks=1:15) +
	theme_classic(base_size=16) +
	theme(strip.background=element_blank(),strip.text=element_text(hjust=0)) +
	xlab('Number of regions') + ylab(paste0('Proportion of ',substr(tolower(predictor.label),1,1),'DEGs'))
ggsave(p,file=paste0('figures/model_results_n_shared_',tolower(predictor.label),'_magnitude_by_region_group.pdf'),width=15,height=5,useDingbats=FALSE)


p = ggplot(subset(all.set.share.proportions,TRUE),aes(n.regions,value)) +
	geom_bar(stat='identity') +
	facet_grid(cols=vars(set),rows=vars(variable),scales='free_x',space='fixed') +
	scale_x_continuous(breaks=1:15) +
	theme_classic(base_size=16) +
	theme(strip.background=element_blank(),strip.text=element_text(hjust=0),strip.background.y=element_blank(),strip.text.y=element_blank()) +
	xlab('Number of regions') + ylab(paste0('Proportion of ',substr(tolower(predictor.label),1,1),'DEGs'))
ggsave(p,file=paste0('figures/model_results_n_shared_',tolower(predictor.label),'_all_by_region_group.pdf'),width=12,height=7,useDingbats=FALSE)





# Tree comparison

library(dendextend)
library(ape)
library(parallel)


# mash.beta.strong = mash.beta[apply(mash.lfsr,1,function(x) any(x < fsr.cutoff)),]

# m.s = matrix(as.numeric(scale(as.numeric(t(mash.beta.strong)))),nrow=ncol(mash.beta.strong),ncol=nrow(mash.beta.strong),dimnames=list(colnames(mash.beta.strong),rownames(mash.beta.strong)))

for (method in c('emma','mash')) {

	qval = if (method == 'emma') 'qval' else if (method == 'mash') 'lfsr'
	beta.matrix = get(paste0(method,'.beta'))
	qval.matrix = get(paste(method,qval,sep='.'))

	# Scale beta matrix
 	m.s = matrix(scale(t(beta.matrix)),nrow=ncol(beta.matrix),dimnames=list(colnames(beta.matrix),rownames(beta.matrix)))
#	m.s = matrix(as.numeric(scale(as.numeric(t(beta.matrix)))),nrow=ncol(beta.matrix),dimnames=list(colnames(beta.matrix),rownames(beta.matrix)))

	m.cluster = hclust(dist(m.s))
	m.cluster.tree = as.dendrogram(m.cluster)

	tree.colors = region.colors
	names(tree.colors) = names(keep.genes)

	colReg = function(n) {
		if (is.leaf(n)) {
			a = attributes(n)
			labCol = as.character(tree.colors[a$label])
			attr(n, 'nodePar') = list(lab.cex = 1, lab.col = labCol, col = labCol, pch=19, cex=0.5)
			attr(n, 'edgePar') = list(col = labCol, lwd=2.5)
		}
		n
	}

	edgeReg = list(lwd=2.5)

	m.cluster.tree = dendrapply(m.cluster.tree, colReg)

	pdf(file=paste0('figures/model_results_',method,'_hclust_dendrogram_regions.pdf'),useDingbats=FALSE,height=10)
		plot(m.cluster.tree,edgePar=edgeReg,horiz=TRUE,yaxt='n')
	dev.off()

	# Bootstrap libraries

	m.cluster.bootstrap = mclapply(1:1000,function(i) {
		m.x = m.s[,apply(qval.matrix,1,function(x) any(x < rnorm(1,mean=fsr.cutoff,sd=fsr.cutoff/4)))]
		m.x = m.x[,sample(1:ncol(m.x),replace=TRUE)]
	#	m.x = m.s[,sample(1:ncol(m.s),replace=TRUE)]
		as.phylo(hclust(dist(m.x)))
	},mc.cores=n.cores)

	# Visualize uncertainty
	library(phangorn)

	m.cluster.bootstrap.multiphylo = do.call(c,m.cluster.bootstrap)

	r.cluster.tree = readRDS('checkpoints/dendrogram_gene_expression_best_tree.rds')

	pdf(file=paste0('figures/model_results_',method,'_hclust_dendrogram_regions_bootstrap.pdf'),useDingbats=FALSE,height=10)
		densiTree(
			m.cluster.bootstrap.multiphylo,
			alpha=0.005,
			consensus=as.phylo(r.cluster.tree),
			direction='leftwards',
			scaleX=TRUE,
			col='#000000',
			width=1, lty=1, cex=1, font=1,
			tip.color=rapply(r.cluster.tree,function(n) {
				if (is.leaf(n)) attr(n, 'edgePar')[['col']]
			}), srt=0, adj=0,
			label.offset=0,
			scale.bar = FALSE,
			jitter=list(amount = 0.02, random = TRUE)
		)
		# Add the consensus manually
		consensus = as.phylo(m.cluster.tree)
		xy = ape::plotPhyloCoor(consensus, tip.height = 1:length(keep.genes), direction = 'leftwards')
		xx = xy[,1]; yy = xy[,2]; xx = xx/max(xx); xx = xx + (1 - max(xx))
	#	ape::cladogram.plot(consensus$edge, xx, yy, edge.color = '#0000ff', edge.width = 1, edge.lty = 1)
	dev.off()

	b.cluster.tree = dendrapply(as.dendrogram(hclust(dist(m.s[,apply(qval.matrix,1,function(x) any(x < fsr.cutoff))]))),colReg)
#	b.cluster.tree = consensus(m.cluster.bootstrap.multiphylo)

	pdf(file=paste0('figures/model_results_',method,'_hclust_dendrogram_regions_bootstrap_reordered.pdf'),useDingbats=FALSE,height=10)
		densiTree(
			m.cluster.bootstrap.multiphylo,
			alpha=0.005,
			consensus=as.phylo(b.cluster.tree),
			direction='leftwards',
			scaleX=TRUE,
			col='#000000',
			width=1, lty=1, cex=1, font=1,
			tip.color=rapply(b.cluster.tree,function(n) {
				if (is.leaf(n)) attr(n, 'edgePar')[['col']]
			}), srt=0, adj=0,
			label.offset=0,
			scale.bar = FALSE,
			jitter=list(amount = 0.02, random = TRUE)
		)
		# Add the consensus manually
		consensus = as.phylo(b.cluster.tree)
		xy = ape::plotPhyloCoor(consensus, tip.height = 1:length(keep.genes), direction = 'leftwards')
		xx = xy[,1]; yy = xy[,2]; xx = xx/max(xx); xx = xx + (1 - max(xx))
#		ape::cladogram.plot(consensus$edge, xx, yy, edge.color = '#0000ff', edge.width = 1, edge.lty = 1)
	dev.off()

	r.cluster.bootstrap.multiphylo = readRDS('checkpoints/dendrogram_gene_expression_bootstrap_tree.rds')

	pdf(file=paste0('figures/model_results_',method,'_hclust_dendrogram_regions_bootstrap_expression_reordered.pdf'),useDingbats=FALSE,height=10)
		densiTree(
			r.cluster.bootstrap.multiphylo,
			alpha=0.005,
			consensus=as.phylo(b.cluster.tree),
			direction='rightwards',
			scaleX=TRUE,
			col='#000000',
			width=1, lty=1, cex=1, font=1,
			tip.color=rapply(b.cluster.tree,function(n) {
				if (is.leaf(n)) attr(n, 'edgePar')[['col']]
			}), srt=0, adj=0,
			label.offset=0,
			scale.bar = FALSE,
			jitter=list(amount = 0.02, random = TRUE)
		)
		# Add the consensus manually
		consensus = as.phylo(r.cluster.tree)
		xy = ape::plotPhyloCoor(consensus, tip.height = 1:length(keep.genes), direction = 'rightwards')
		xx = xy[,1]; yy = xy[,2]; xx = xx/max(xx); xx = xx + (1 - max(xx))
#		ape::cladogram.plot(consensus$edge, xx, yy, edge.color = '#0000ff', edge.width = 1, edge.lty = 1)
	dev.off()
}

go.enrichment.results.all = readRDS('checkpoints/topgo_results.rds')
disease.enrichment.results.all = readRDS('checkpoints/disease_enrichment_results.rds')

library(egg)

for (i in c('fet','kst')) {

	go.enrichment.results = subset(go.enrichment.results.all,grepl(paste0('^',toupper(substr(i,1,2))),test))

	# Visualize GO results
	p = ggplot(
		do.call(rbind,lapply(split(go.enrichment.results,list(go.enrichment.results$direction, go.enrichment.results$region)),function(x) {
			within(data.frame(unique(x[c('direction','region')]),expected=-log10(seq(1/nrow(x),1,1/nrow(x))),observed=-log10(quantile(x$pval,seq(1/nrow(x),1,1/nrow(x)),na.rm=TRUE))),{
				region=factor(region,levels=c('all',region.levels),labels=gsub('ACCg','ACC',c('wbaDEGs',region.levels)))
				direction=factor(direction,levels=c('increase','decrease'),labels=c('Increase','Decrease'))
			})
		})),
		aes(expected,observed,color=direction)) +
#		geom_smooth(method=lm,se=FALSE,color='black',size=0.2,linetype=2) +
		geom_point_rast(size=0.1) +
		geom_abline(slope=1,col='black',size=0.2) +
		facet_wrap(~region,nrow=region.rows + 1) +
		scale_x_continuous(breaks=seq(0,ceiling(-log10(min(1/(max(table(go.enrichment.results$direction,go.enrichment.results$region))-1)))),1)) +
		scale_y_continuous(breaks=seq(0,ceiling(max(-log10(go.enrichment.results$pval))),5)) +
		scale_color_manual(name='Direction',values=c('#386cb0','#f0027f')) +
	#	coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
		xlab(expression(-log[10] * ('Expected'~italic(p)))) +
		ylab(expression(-log[10] * ('Observed'~italic(p)))) +
		theme_article(base_size=16) +
		theme() +
		guides(color=guide_legend(override.aes=list(size=3)))
	ggsave(p,file=paste0('figures/model_results_pval_enrichment_topgo_',i,'_',tolower(predictor.label),'_qqplot.pdf'),useDingbats=FALSE)

	p = ggplot(
		do.call(rbind,lapply(split(subset(go.enrichment.results,region=='all'),subset(go.enrichment.results,region=='all')$direction),function(x) {
			within(data.frame(unique(x[c('direction','region')]),expected=-log10(seq(1/nrow(x),1,1/nrow(x))),observed=-log10(quantile(x$pval,seq(1/nrow(x),1,1/nrow(x)),na.rm=TRUE))),{
				direction=factor(direction,levels=c('increase','decrease'),labels=c('Increase','Decrease'))
			})
		})),
		aes(expected,observed,color=direction)) +
		geom_smooth(method=lm,se=FALSE,color='black',size=0.2,linetype=2) +
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
	ggsave(p,file=paste0('figures/model_results_pval_enrichment_topgo_',i,'_',tolower(predictor.label),'_qqplot_union.pdf'),useDingbats=FALSE)

	disease.enrichment.results = subset(disease.enrichment.results.all,grepl(paste0('^',toupper(substr(i,1,2))),test) & dataset == 'DISEASES')

	# Visualize GO results
	p = ggplot(
		do.call(rbind,lapply(split(disease.enrichment.results,list(disease.enrichment.results$direction, disease.enrichment.results$region)),function(x) {
			within(data.frame(unique(x[c('direction','region')]),expected=-log10(seq(1/nrow(x),1,1/nrow(x))),observed=-log10(quantile(x$pval,seq(1/nrow(x),1,1/nrow(x)),na.rm=TRUE))),{
				region=factor(region,levels=c('all',region.levels),labels=gsub('ACCg','ACC',c('wbaDEGs',region.levels)))
				direction=factor(direction,levels=c('increase','decrease'),labels=c('Increase','Decrease'))
			})
		})),
		aes(expected,observed,color=direction)) +
#		geom_smooth(method=lm,se=FALSE,color='black',size=0.2,linetype=2) +
		geom_point_rast(size=0.1) +
		geom_abline(slope=1,col='black',size=0.2) +
		facet_wrap(~region,nrow=region.rows + 1) +
		scale_x_continuous(breaks=seq(0,ceiling(-log10(min(1/(max(table(disease.enrichment.results$direction,disease.enrichment.results$region))-1)))),1)) +
		scale_y_continuous(breaks=seq(0,ceiling(max(-log10(disease.enrichment.results$pval))),5)) +
		scale_color_manual(name='Direction',values=c('#386cb0','#f0027f')) +
	#	coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
		xlab(expression(-log[10] * ('Expected'~italic(p)))) +
		ylab(expression(-log[10] * ('Observed'~italic(p)))) +
		theme_article(base_size=16) +
		theme() +
		guides(color=guide_legend(override.aes=list(size=3)))
	ggsave(p,file=paste0('figures/model_results_pval_enrichment_disease_',i,'_',tolower(predictor.label),'_qqplot.pdf'),useDingbats=FALSE)

	p = ggplot(
		do.call(rbind,lapply(split(subset(disease.enrichment.results,region=='all'),subset(disease.enrichment.results,region=='all')$direction),function(x) {
			within(data.frame(unique(x[c('direction','region')]),expected=-log10(seq(1/nrow(x),1,1/nrow(x))),observed=-log10(quantile(x$pval,seq(1/nrow(x),1,1/nrow(x)),na.rm=TRUE))),{
				direction=factor(direction,levels=c('increase','decrease'),labels=c('Increase','Decrease'))
			})
		})),
		aes(expected,observed,color=direction)) +
		geom_smooth(method=lm,se=FALSE,color='black',size=0.2,linetype=2) +
		geom_point_rast(size=0.5) +
		geom_abline(slope=1,col='black',size=0.2) +
		scale_x_continuous(breaks=seq(0,ceiling(-log10(min(1/(max(table(subset(disease.enrichment.results,region=='all')$direction,subset(disease.enrichment.results,region=='all')$region))-1)))),1)) +
		scale_y_continuous(breaks=seq(0,ceiling(max(-log10(subset(disease.enrichment.results,region=='all')$pval))),1)) +
		scale_color_manual(name='Direction',values=c('#386cb0','#f0027f')) +
	#	coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
		xlab(expression(-log[10] * ('Expected'~italic(p)))) +
		ylab(expression(-log[10] * ('Observed'~italic(p)))) +
		theme_classic() +
		theme() +
		guides(color=guide_legend(override.aes=list(size=3)))
	ggsave(p,file=paste0('figures/model_results_pval_enrichment_disease_',i,'_',tolower(predictor.label),'_qqplot_union.pdf'),useDingbats=FALSE)

}


# QQ plot for HOMER
homer.1 = readRDS('checkpoints/homer_results_1.rds')
homer.2 = readRDS('checkpoints/homer_results_2.rds')

homer.1$direction = 'increase'
homer.2$direction = 'decrease'

homer = rbind(homer.1,homer.2)

p = ggplot(
	do.call(rbind,lapply(split(homer,homer$direction),function(x) {
		within(data.frame(direction=unique(x$direction),expected=-log10(seq(1/nrow(x),1,1/nrow(x))),observed=-log10(quantile(x$p.value,seq(1/nrow(x),1,1/nrow(x)),na.rm=TRUE))),{
			direction=factor(direction,levels=c('increase','decrease'),labels=c('Increase','Decrease'))
		})
	})),
	aes(expected,observed,color=direction)) +
#	geom_smooth(method=lm,se=FALSE,color='black',size=0.2,linetype=2) +
	geom_point_rast(size=0.5) +
	geom_abline(slope=1,col='black',size=0.2) +
	scale_x_continuous(breaks=seq(0,ceiling(-log10(min(1/(max(table(homer$direction))-1)))),1)) +
	scale_y_continuous(breaks=seq(0,ceiling(max(-log10(homer$pval))),1)) +
	scale_color_manual(name='Direction',values=c('#386cb0','#f0027f')) +
#	coord_fixed(xlim=c(0,1),ylim=c(0,1)) +
	xlab(expression(-log[10] * ('Expected'~italic(p)))) +
	ylab(expression(-log[10] * ('Observed'~italic(p)))) +
	theme_classic(base_size=24) +
	theme() +
	guides(color=guide_legend(override.aes=list(size=3)))
ggsave(p,file=paste0('figures/model_results_pval_enrichment_homer_',tolower(predictor.label),'_qqplot_union.pdf'),useDingbats=FALSE)


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
			count = as.integer(colSums(mash.lfsr < threshold)),
			proportion = as.numeric(
				colSums(mash.lfsr < threshold) /
					sum(mash.lfsr < threshold)
			),
			direction = 'Increase or decrease',
			type = 'total'
		),
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
			count = as.integer(colSums(this.inc.lfsr < threshold)),
			proportion = as.numeric(
				colSums(this.inc.lfsr < threshold) /
					sum(this.inc.lfsr < threshold)
			),
			direction = 'Increase or decrease',
			type = 'unique'
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
	out$direction = factor(out$direction,levels=c('Increase or decrease','Increase','Decrease'))
	out$region = factor(out$region,levels=region.levels)
	out
}))

# The y-axis below can be switched between proportion and count below. The rank order will stay the same
# Plot for total number of significant genes per region
p = ggplot(droplevels(subset(which.region.tally,direction %in% c('Increase','Decrease') & type=='total')),aes(threshold,count,color=region)) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (count)')
ggsave(p,file=paste0('figures/model_results_sig_counts_total_',tolower(predictor.label),'_region_tally.pdf'),useDingbats=FALSE,height=5)

# Plot of unique significant genes per region
p = ggplot(droplevels(subset(which.region.tally,direction %in% c('Increase','Decrease') & type=='unique')),aes(threshold,count,color=region)) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (count)')
ggsave(p,file=paste0('figures/model_results_sig_counts_unique_',tolower(predictor.label),'_region_tally.pdf'),useDingbats=FALSE,height=5)

# Plot for total number of significant genes per region
p = ggplot(droplevels(subset(which.region.tally,direction %in% c('Increase','Decrease') & type=='total')),aes(threshold,proportion,color=region)) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (proportion)')
ggsave(p,file=paste0('figures/model_results_sig_proportions_total_',tolower(predictor.label),'_region_tally.pdf'),useDingbats=FALSE,height=5)

# Plot of unique significant genes per region
p = ggplot(droplevels(subset(which.region.tally,direction %in% c('Increase','Decrease') & type=='unique')),aes(threshold,proportion,color=region)) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (proportion)')
ggsave(p,file=paste0('figures/model_results_sig_proportions_unique_',tolower(predictor.label),'_region_tally.pdf'),useDingbats=FALSE,height=5)




# The y-axis below can be switched between proportion and count below. The rank order will stay the same
# Plot for total number of significant genes per region
p = ggplot(droplevels(subset(which.region.tally,direction %in% c('Increase or decrease','Increase','Decrease') & type=='total')),aes(threshold,count,color=region)) +
	geom_line() +
	facet_wrap(~direction,scales='free_y',ncol=1) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank(),strip.background=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (count)')
ggsave(p,file=paste0('figures/model_results_sig_counts_total_',tolower(predictor.label),'_region_tally_withtotal.pdf'),useDingbats=FALSE,height=5)

# Plot of unique significant genes per region
p = ggplot(droplevels(subset(which.region.tally,direction %in% c('Increase or decrease','Increase','Decrease') & type=='unique')),aes(threshold,count,color=region)) +
	geom_line() +
	facet_wrap(~direction,scales='free_y',ncol=1) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (count)')
ggsave(p,file=paste0('figures/model_results_sig_counts_unique_',tolower(predictor.label),'_region_tally_withtotal.pdf'),useDingbats=FALSE,height=5)

# Plot for total number of significant genes per region
p = ggplot(droplevels(subset(which.region.tally,direction %in% c('Increase or decrease','Increase','Decrease') & type=='total')),aes(threshold,proportion,color=region)) +
	geom_line() +
	facet_wrap(~direction,scales='free_y',ncol=1) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank(),strip.background=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (proportion)')
ggsave(p,file=paste0('figures/model_results_sig_proportions_total_',tolower(predictor.label),'_region_tally_withtotal.pdf'),useDingbats=FALSE,height=5)

# Plot of unique significant genes per region
p = ggplot(droplevels(subset(which.region.tally,direction %in% c('Increase or decrease','Increase','Decrease') & type=='unique')),aes(threshold,proportion,color=region)) +
	geom_line() +
	facet_wrap(~direction,scales='free_y',ncol=1) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=16) +
	theme(legend.title=element_blank()) +
	guides(color=guide_legend(override.aes=list(linetype=1,size=1))) +
	xlab('LFSR threshold') +
	ylab('Significant genes (proportion)')
ggsave(p,file=paste0('figures/model_results_sig_proportions_unique_',tolower(predictor.label),'_region_tally_withtotal.pdf'),useDingbats=FALSE,height=5)





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
ggsave(p,file=paste0('figures/model_results_sig_counts_',tolower(predictor.label),'_n_regions_tally.pdf'),useDingbats=FALSE,height=5)

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
ggsave(p,file=paste0('figures/model_results_sig_proportions_',tolower(predictor.label),'_n_regions_tally.pdf'),useDingbats=FALSE,height=5)







save(list=ls(),file='checkpoints/checkpoint_model_results_visualization.RData')

