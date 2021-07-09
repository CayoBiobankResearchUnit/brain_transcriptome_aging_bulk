#!/usr/bin/env Rscript

source('scripts/_include_options.R')

meta.animal = readRDS('checkpoints/cayo_bulkbrain_animal_metadata.rds')
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')

library(ggplot2)
library(egg)

# Histogram of age

p = ggplot(within(meta.animal,{sex=factor(sex,labels=c('Female','Male'))}),aes_string(age.variable,fill=sex.variable)) +
	geom_histogram(binwidth=0.5) +
	facet_wrap(as.formula(paste0('~',sex.variable)),nrow=2) +
	scale_fill_manual(values=sex.colors) +
	scale_x_continuous(limits=c(0,ceiling(max(meta[[age.variable]])/5)*5+5)) +
	xlab('Age (years)') + ylab('Count') +
	theme_classic(base_size=16) +
	theme(legend.position='none',strip.background=element_blank(),strip.text=element_text())
ggsave(p,file=paste0('figures/metadata_age_histogram.pdf'),useDingbats=FALSE,height=4)

p = ggplot(within(meta.animal,{sex=factor(sex,labels=c('Female','Male'))}),aes_string(age.variable,fill=sex.variable)) +
	geom_histogram(binwidth=0.5) +
	facet_wrap(as.formula(paste0('~',sex.variable)),nrow=2) +
	scale_fill_manual(values=sex.colors) +
	scale_x_continuous(limits=c(0,30)) +
	xlab('Age (years)') + ylab('Count') +
	theme_classic(base_size=24) +
	theme(legend.position='none',strip.background=element_blank(),strip.text=element_text())
ggsave(p,file=paste0('figures/metadata_age_histogram_article.pdf'),useDingbats=FALSE,height=5)

# Histogram of rank

p = ggplot(within(meta.animal,{sex=factor(sex,labels=c('Female','Male'))}),aes_string(rank.variable,fill=sex.variable)) +
	geom_histogram(binwidth=0.02) +
	facet_wrap(as.formula(paste0('~',sex.variable)),nrow=2) +
	scale_fill_manual(values=sex.colors) +
	scale_x_continuous(limits=c(0,1)) +
	scale_y_continuous(limits=c(0,2),breaks=seq(0,2,1)) +
	xlab('Rank (scaled)') + ylab('Count') +
	theme_classic() +
	theme(legend.position='none')
ggsave(p,file=paste0('figures/metadata_rank_histogram.pdf'),useDingbats=FALSE,height=4)

# Dotplot of rank

p = ggplot(within(meta.animal,{sex=factor(sex,labels=c('Female','Male'))}),aes_string(rank.variable,'1',color=sex.variable)) +
	geom_point() +
	facet_wrap(as.formula(paste0('~',sex.variable)),nrow=2) +
	scale_color_manual(values=sex.colors) +
	scale_x_continuous(limits=c(0,1)) +
	scale_y_continuous(limits=c(0,2),breaks=seq(0,2,1)) +
	xlab('Rank (scaled)') +
	theme_classic() +
	theme(legend.position='none',axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
ggsave(p,file=paste0('figures/metadata_rank_dotplot.pdf'),useDingbats=FALSE,height=4)

# Age by rank

p = ggplot(within(meta.animal,{sex=factor(sex,labels=c('Female','Male'))}),aes_string(age.variable,rank.variable,color=sex.variable)) +
	geom_point() +
	geom_smooth(method=lm) +
	scale_color_manual(name='Sex',values=sex.colors) +
	scale_x_continuous() +
	scale_y_continuous(limits=c(0,1)) +
	xlab('Age (years)') + ylab('Rank (scaled)') +
	theme_classic() +
	theme() +
	guides(color = guide_legend(override.aes = list(size = 3,linetype=0,fill=NA)))
ggsave(p,file=paste0('figures/metadata_age_by_rank.pdf'),useDingbats=FALSE,height=4)

# Histogram of age per region

region.labeller = as.character(apply(matrix(1:length(region.levels),nrow=(length(region.levels)/region.rows),ncol=region.rows),2,function(x) rep(region.levels[x],2)))
names(region.labeller) = as.character(apply(matrix(1:length(region.levels),nrow=(length(region.levels)/region.rows),ncol=region.rows),2,function(x) c(paste(region.levels[x],'f',sep='.'),paste(region.levels[x],'m',sep='.'))))

region.labeller = gsub('ACCg','ACC',region.labeller)

p = ggplot(
		within(meta,{
			sex.labeled=factor(sex,labels=c('Female','Male'))
			facet=factor(paste(Region,sex,sep='.'),levels=names(region.labeller))
		}),
		aes_string(age.variable,fill='sex.labeled')
	) +
	geom_histogram(binwidth=0.5) +
	facet_wrap(
		~facet,
		nrow=region.rows * 2,
		labeller=labeller(facet=region.labeller)
	) +
	scale_fill_manual(name='Sex',values=sex.colors) +
	scale_x_continuous(limits=c(0,ceiling(max(meta[[age.variable]])/5)*5+5)) +
	xlab('Age (years)') + ylab('Count') +
	theme_classic() +
	theme(strip.background=element_blank())

# ggplotGrob draws to device. Disable this by sending the output to /dev/null (not ideal but I give up on better solutions)
pdf(file='/dev/null')
	g = ggplotGrob(p)
dev.off()

library(grid)

# Remove strip text
pdf(file=paste0('figures/metadata_age_histogram_facet_regions.pdf'),useDingbats=FALSE,height=4)
	strips = subset(g$layout,(as.logical(as.integer(gsub('strip-t-[0-9]+-','',g$layout$name)) %% 2) %in% FALSE))
	grid.newpage()
	grid.draw(g[-unique(strips$b),])
dev.off()

# Sex by rank

# RIN histogram

p = ggplot(within(meta,{sex=factor(sex,labels=c('Female','Male'))}),aes_string('RIN',fill=sex.variable)) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(as.formula(paste0('~',sex.variable)),nrow=2) +
	scale_fill_manual(values=sex.colors) +
	scale_x_continuous(limits=c(0,max(meta$RIN)),breaks=seq(0,10,1)) +
	xlab('RNA integrity') + ylab('Count') +
	theme_classic() +
	theme(legend.position='none')
ggsave(p,file=paste0('figures/metadata_rin_histogram.pdf'),useDingbats=FALSE,height=4)

# Histogram of RIN per region

p = ggplot(
		within(meta,{
			sex.labeled=factor(sex,labels=c('Female','Male'))
			facet=factor(paste(Region,sex,sep='.'),levels=names(region.labeller))
		}),
		aes_string('RIN',fill='sex.labeled')
	) +
	geom_histogram(binwidth=0.2) +
	facet_wrap(
		~facet,
		nrow=region.rows * 2,
		labeller=labeller(facet=region.labeller)
	) +
	scale_fill_manual(name='Sex',values=sex.colors) +
	scale_x_continuous(limits=c(0,max(meta$RIN)),breaks=seq(0,10,2)) +
	scale_y_continuous(limits=c(0,6),breaks=seq(0,6,3)) +
	xlab('RNA integrity') + ylab('Count') +
	theme_classic() +
	theme(strip.background=element_blank())

# ggplotGrob draws to device. Disable this by sending the output to /dev/null (not ideal but I give up on better solutions)
pdf(file='/dev/null')
	g = ggplotGrob(p)
dev.off()

# Remove strip text
pdf(file=paste0('figures/metadata_rin_histogram_facet_regions.pdf'),useDingbats=FALSE,height=4)
	strips = subset(g$layout,(as.logical(as.integer(gsub('strip-t-[0-9]+-','',g$layout$name)) %% 2) %in% FALSE))
	grid.newpage()
	grid.draw(g[-unique(strips$b),])
dev.off()

# RIN by age

p = ggplot(within(meta,{sex=factor(sex,labels=c('Female','Male'))}),aes_string(age.variable,'RIN',color=sex.variable)) +
	geom_point() +
	geom_smooth(method=lm,se=FALSE) +
	scale_color_manual(name='Sex',values=sex.colors) +
	scale_x_continuous(limits=c(0,ceiling(max(meta[[age.variable]])/5)*5+5)) +
	scale_y_continuous(limits=c(0,max(meta$RIN)),breaks=seq(0,10,1)) +
	coord_fixed(ratio=(ceiling(max(meta[[age.variable]])/5)*5+5)/10) +
	xlab('Age (years)') + ylab('RNA integrity') +
	theme_classic(base_size=16) +
	theme() +
	guides(color = guide_legend(override.aes = list(size = 3,linetype=0,fill=NA)))
ggsave(p,file=paste0('figures/metadata_rin_by_age.pdf'),useDingbats=FALSE)

p = ggplot(within(meta,{sex=factor(sex,labels=c('Female','Male'))}),aes_string(age.variable,'RIN',color=sex.variable)) +
	geom_point() +
#	geom_smooth(method=lm,se=FALSE) +
	facet_wrap(~Region,nrow=region.rows) +
	scale_color_manual(name='Sex',values=sex.colors) +
	scale_x_continuous(limits=c(0,ceiling(max(meta[[age.variable]])/5)*5+5)) +
	scale_y_continuous(limits=c(0,max(meta$RIN)),breaks=seq(0,10,2)) +
	coord_fixed(ratio=(ceiling(max(meta[[age.variable]])/5)*5+5)/10) +
	xlab('Age (years)') + ylab('RNA integrity') +
	theme_classic() +
	theme() +
	guides(color = guide_legend(override.aes = list(size = 3,linetype=0,fill=NA)))
ggsave(p,file=paste0('figures/metadata_rin_by_age_facet_regions.pdf'),useDingbats=FALSE,height=5)


levels(meta$Region) = gsub('ACCg','ACC',levels(meta$Region))

p = ggplot(within(meta,{sex=factor(sex,labels=c('Female','Male'))}),aes_string(age.variable,'RIN',color=sex.variable)) +
	geom_point() +
#	geom_smooth(method=lm,se=FALSE) +
	facet_wrap(~Region,nrow=region.rows) +
	scale_color_manual(name='Sex',values=sex.colors) +
	scale_x_continuous(limits=c(0,ceiling(max(meta[[age.variable]])/5)*5+5)) +
	scale_y_continuous(limits=c(0,max(meta$RIN)),breaks=seq(0,10,2)) +
	coord_fixed(ratio=(ceiling(max(meta[[age.variable]])/5)*5+5)/10) +
	xlab('Age (years)') + ylab('RNA integrity') +
	theme_classic() +
	theme(strip.background=element_blank()) +
	guides(color = guide_legend(override.aes = list(size = 3,linetype=0,fill=NA)))
ggsave(p,file=paste0('figures/metadata_rin_by_age_facet_regions_article.pdf'),useDingbats=FALSE,height=5)

# Social metadata

for (i in 1:length(social.metrics)) {
	sm = social.metrics[i]
	sl = social.labels[i]

	p = ggplot(within(meta.animal,{sex=factor(sex,labels=c('Female','Male'))}),aes_string(sm,fill=sex.variable)) +
		geom_histogram(binwidth=0.1) +
		facet_wrap(as.formula(paste0('~',sex.variable)),nrow=2) +
		scale_fill_manual(values=sex.colors) +
		scale_x_continuous(breaks=seq(0,0.5,1)) +
		xlab(sl) + ylab('Count') +
		theme_classic() +
		theme(legend.position='none')
	ggsave(p,file=paste0('figures/metadata_',tolower(sm),'_histogram.pdf'),useDingbats=FALSE,height=4)

}

# Holding time

p = ggplot(within(meta.animal,{sex=factor(sex,labels=c('Female','Male'))}),aes_string('days_between',fill=sex.variable)) +
	geom_histogram(binwidth=1) +
	facet_wrap(as.formula(paste0('~',sex.variable)),nrow=2) +
	scale_fill_manual(values=sex.colors) +
	scale_x_continuous(breaks=seq_len(max(meta$days_between))) +
	xlab('Holding time (days)') + ylab('Count') +
	theme_classic() +
	theme(legend.position='none')
ggsave(p,file=paste0('figures/metadata_','daysbetween','_histogram.pdf'),useDingbats=FALSE,height=4)

p = ggplot(within(meta.animal,{sex=factor(sex,labels=c('Female','Male'))}),aes_string('GroomOUT','GroomIN',color=sex.variable)) +
	geom_point() +
	coord_equal() +
	scale_color_manual(values=sex.colors) +
	scale_x_continuous(breaks=seq(0,1,0.5)) +
	scale_y_continuous(breaks=seq(0,1,0.5)) +
	xlab('Active grooming index') + ylab('Passive grooming index') +
	theme_classic() +
	theme(legend.position='none')
ggsave(p,file=paste0('figures/metadata_','groom_inVsout','_histogram.pdf'),useDingbats=FALSE,width=4,height=4)
