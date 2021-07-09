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

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(egg)
library(ggrepel)

go.enrichment.results.all = readRDS('checkpoints/topgo_results.rds')
disease.enrichment.results.all = readRDS('checkpoints/disease_enrichment_results.rds')
dglm.go.enrichment.results.all = readRDS('checkpoints/dglm_topgo_results.rds')
dglm.disease.enrichment.results.all = readRDS('checkpoints/dglm_disease_results.rds')


go.results = subset(go.enrichment.results.all,set == 'union' & region == 'all' & test == 'FET')

go.results = do.call(rbind,lapply(split(go.results,go.results$direction),function(x) {
	x = x[order(x$qval,x$pval),]
	x$i = 1:nrow(x)
	x
}))
go.results$direction = factor(go.results$direction,levels=c('increase','decrease'),labels=c('Increase','Decrease'))

go.sig = structure(list(go_id = structure(c(5L, 2L, 8L, 3L, 4L, 9L, 6L, 
7L, 1L, 10L), .Label = c("GO:0006629", "GO:0007154", "GO:0007165", 
"GO:0007610", "GO:0023052", "GO:0040011", "GO:0044283", "GO:0050807", 
"GO:0050896", "GO:0055114"), class = "factor"), go_label = structure(c(9L, 
2L, 6L, 8L, 1L, 7L, 4L, 10L, 3L, 5L), .Label = c("behavior", 
"cell communication", "lipid metabolic process", "locomotion", 
"oxidation-reduction process", "regulation of synapse organization", 
"response to stimulus", "signal transduction", "signaling", "small molecule biosynthetic process"
), class = "factor"), direction = structure(c(2L, 2L, 2L, 2L, 
2L, 2L, 2L, 1L, 1L, 1L), .Label = c("Increase", "Decrease"), class = "factor")), row.names = c(NA, 
-10L), class = "data.frame")

go.sig$sig = TRUE

go.sig = merge(go.results,go.sig,by=c('go_id','direction'),all.x=TRUE,all.y=TRUE)
go.sig$sig[is.na(go.sig$sig)] = FALSE

go.sig$display = factor(go.sig$sig,levels=c('TRUE','FALSE'))

disease.results = subset(disease.enrichment.results.all,set == 'union' & region == 'all' & test == 'FET' & dataset == 'DISEASES')

disease.results = do.call(rbind,lapply(split(disease.results,disease.results$direction),function(x) {
	x = x[order(x$qval,x$pval),]
	x$i = 1:nrow(x)
	x
}))
disease.results$direction = factor(disease.results$direction,levels=c('increase','decrease'),labels=c('Increase','Decrease'))

disease.sig = structure(list(do_id = structure(c(7L, 3L, 5L, 2L, 1L, 4L, 8L, 
9L, 6L), .Label = c("DOID:1059", "DOID:1094", "DOID:12849", "DOID:2030", 
"DOID:3312", "DOID:3613", "DOID:5419", "DOID:8927", "DOID:9351"
), class = "factor"), do_label = structure(c(9L, 3L, 4L, 2L, 
7L, 1L, 8L, 6L, 5L), .Label = c("Anxiety disorder", "ADHD", 
"Autistic disorder", "Bipolar disorder", "Canavan disease", "Diabetes mellitus", 
"Intellectual disability", "Learning disability", "Schizophrenia"
), class = "factor"), direction = structure(c(2L, 2L, 2L, 2L, 
2L, 2L, 2L, 1L, 1L), .Label = c("Increase", "Decrease"), class = "factor")), row.names = c(NA, 
-9L), class = "data.frame")

disease.sig$sig = TRUE

disease.sig = merge(disease.results,disease.sig,by=c('do_id','direction'),all.x=TRUE,all.y=TRUE)
disease.sig$sig[is.na(disease.sig$sig)] = FALSE

disease.sig$display = factor(disease.sig$sig,levels=c('TRUE','FALSE'))


p1 = ggplot(
		data=subset(go.sig,pval < 0.001),
		aes(
			x=-i,y=-log10(pval),
			label=go_label
		)
	) +
	geom_point(aes(color=display,shape=display,size=display)) +
	geom_point(data=subset(go.sig,display=='TRUE'),aes(x=-i,y=-log10(pval)),color='black',shape=19,size=1) +
	geom_label_repel(
		size=5,
		max.time=10,
		max.iter=1e5,
		nudge_x=-12,
		nudge_y=2,
		force=5e2,
		na.rm=TRUE
	) +
	facet_wrap(~direction,nrow=2,scales='free_x') +
	scale_y_continuous(limits=c(0,max(c(-log10(go.results$pval),-log10(disease.results$pval))))) +
	scale_color_manual(values=c('red','black')) +
	scale_shape_manual(values=c(21,19)) +
	scale_size_manual(values=c(3,1)) +
	theme_article(base_size=24) +
	theme(
		strip.background=element_blank(),
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		legend.position='none',
		panel.spacing=unit(0,'lines')
	) +
	xlab('Biological process') + ylab(expression(-log[10]*'('*italic(P)*')'))

p2 = ggplot(
		data=subset(disease.sig,pval < 0.001),
		aes(
			x=-i,y=-log10(pval),
			label=do_label
		)
	) +
	geom_point(aes(color=display,shape=display,size=display)) +
	geom_point(data=subset(disease.sig,display=='TRUE'),aes(x=-i,y=-log10(pval)),color='black',shape=19,size=1) +
	geom_label_repel(
		size=5,
		max.time=10,
		max.iter=1e5,
		nudge_x=-12,
		nudge_y=2,
		force=5e2,
		na.rm=TRUE
	) +
	facet_wrap(~direction,nrow=2,scales='free_x') +
	scale_y_continuous(limits=c(0,max(c(-log10(go.results$pval),-log10(disease.results$pval))))) +
	scale_color_manual(values=c('red','black')) +
	scale_shape_manual(values=c(21,19)) +
	scale_size_manual(values=c(3,1)) +
	theme_article(base_size=24) +
	theme(
		strip.background=element_blank(),
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		legend.position='none',
		axis.title.y=element_blank(),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		panel.spacing=unit(0,'lines')
	) +
	xlab('Disease')

pdf(file=paste0('figures/enrichment_results_mash_go_diseases_presentation.pdf'),useDingbats=FALSE,height=7,width=9)
	ggarrange(p1,p2,ncol=2,nrow=1,widths=c(1,1),newpage=FALSE,padding=0)
dev.off()

