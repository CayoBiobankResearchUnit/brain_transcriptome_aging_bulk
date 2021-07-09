#!/usr/bin/env Rscript

source('scripts/_include_options.R')

e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')
keep.genes = readRDS('checkpoints/keep_genes.rds')
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')

library(parallel)
library(doParallel)
library(ggplot2)
library(egg)

remove.terms = function(model,covariates) {
	remove.term = function(model,term) {
		# Find term in model matrix/output
		search.string = gsub('\\.','\\\\.',paste0('^',term))
	
		# Number of classes of that term
		n.classes = sum(grepl(search.string,names(coef(model))))

		# Extract betas
		betas = coef(model)[grep(search.string,names(coef(model)))]
	
		# Extract model terms
		mat = model.matrix(model)[,grep(search.string,colnames(model.matrix(model)))]
		if (n.classes > 1) {
			names(betas) = colnames(mat) = gsub(search.string,'',names(betas))
			regress.term = rowSums(t(betas * t(mat)))
		} else {
			regress.term = betas * mat
		}
		as.numeric(replace(regress.term,which(is.na(regress.term)),0))
	}
	Reduce(`+`,lapply(covariates,function(x) {
		remove.term(model,x)
	}))
}

clus = makeCluster(n.cores)
registerDoParallel(cores=n.cores)
clusterExport(clus,varlist=c('e.keep','meta','remove.terms','full.covariates','batch.covariates'),envir=environment())

e.regressed = t(parApply(clus,e.keep,1,function(y) {
	
	# "this" is a data frame with the metadata and the expression for a single gene
	this = meta
	this$e = y


	model = lm(as.formula(paste('e',paste(full.covariates,collapse=' + '),sep=' ~ ')),data=this)

	partial.resid = y - remove.terms(model,batch.covariates)

	return(partial.resid)
	# return(resid(model)) # Return residuals of model NOT including biological terms of interest
}))

stopCluster(clus)

if (FALSE) {

# Variance partitioning

library(variancePartition)
# m = meta

m = meta
m = within(m,{
	Age = exact_age_years
	Batch = Sequencing.batch
	Rank = factor(as.character(ordinal.rank),levels=c('L','M','H'))
	Sex = sex
	Reads = reads_mapped
})

# m$sex = as.integer(factor(m$sex))
f = ~ Age + (1|Rank) + (1|Sex) + RIN + log(Reads) + (1|Region) + (1|Batch) + (1|Individual)

var.part = fitExtractVarPartModel(e.keep,f,m)
vp = sortCols(var.part)

p = plotVarPart( vp ) + theme_classic()

f = ~ Age + (1|Rank) + (1|Sex) + (1|Region) + (1|Individual)

var.part2 = fitExtractVarPartModel(e.regressed,f,m)
vp2 = sortCols(var.part2)
p = plotVarPart( vp2 ) + theme_classic()

out = vector(length=15)
names(out) = region.levels

out = lapply(region.levels,function(i) {
	message(i)
	m.this = subset(m,Region == i)
	e.this = e.regressed[keep.genes[[i]],rownames(m.this)]
	f = ~ Age + (1|Rank) + (1|Sex)
	var.part.this = fitExtractVarPartModel(e.this,f,m.this)
	var.part.this
})

names(out) = region.levels

saveRDS(var.part,file='checkpoints/variance_partitioning_with_batch.rds')
saveRDS(var.part2,file='checkpoints/variance_partitioning_overall.rds')

saveRDS(out,file='checkpoints/variance_partitioning_by_region.rds')

vp.all = do.call(rbind,lapply(region.levels,function(i) {
	x = out[[i]]
	y = reshape2::melt(as.data.frame(x))
	y$Region = factor(i,levels=region.levels)
	y
}))

levels(vp.all$Region) = gsub('ACCg','ACC',levels(vp.all$Region))

vp.all$variable = factor(vp.all$variable,levels=intersect(colnames(vp),unique(vp.all$variable)))

v.colors = c(brewer.pal(ncol(vp)-1,'Set3'),'grey85')
names(v.colors) = colnames(vp)
names(v.colors)[names(v.colors) == 'log(Reads)'] = 'Reads'

library(ggrastr)

p = ggplot() +
	geom_violin(data=vp.all,aes(variable,value * 100,fill=variable),scale='width') +
	geom_boxplot_jitter(data=vp.all,aes(variable,value * 100),width=0.075,fill='#fcfcfc',outlier.colour='black',outlier.size=0.1,outlier.jitter.width=0,outlier.jitter.height=0) +
	facet_wrap(~Region,nrow=region.rows) +
	scale_fill_manual(values=v.colors[levels(vp.all$variable)]) +
	scale_y_continuous(limit=c(0,100),breaks=seq(0,100,25)) +
	theme_article(base_size=16) +
	theme(legend.position='none',axis.text.x=element_text(angle=45,hjust=1,vjust=1),axis.title.x=element_blank()) +
	ylab('Variance explained (%)')
ggsave(p,file=paste0('figures/data_visualization_variance_partition_region.pdf'),useDingbats=FALSE,height=6,width=8)

vp.batch = reshape2::melt(as.data.frame(var.part))
vp.batch$variable = factor(vp.batch$variable,levels=colnames(vp))

levels(vp.batch$variable)[levels(vp.batch$variable) == 'log(Reads)'] = 'Reads'

p = ggplot() +
	geom_violin(data=vp.batch,aes(variable,value * 100,fill=variable),scale='width') +
	geom_boxplot_jitter(data=vp.batch,aes(variable,value * 100),width=0.075,fill='#fcfcfc',outlier.colour='black',outlier.size=0.1,outlier.jitter.width=0,outlier.jitter.height=0) +
	scale_fill_manual(values=v.colors[levels(vp.batch$variable)]) +
	scale_y_continuous(limit=c(0,100),breaks=seq(0,100,25)) +
	theme_classic(base_size=24) +
	theme(legend.position='none',axis.text.x=element_text(angle=45,hjust=1,vjust=1),axis.title.x=element_blank()) +
	ylab('Variance explained (%)')
ggsave(p,file=paste0('figures/data_visualization_variance_partition_batch.pdf'),useDingbats=FALSE,height=7,width=7)

vp.resid = reshape2::melt(as.data.frame(var.part2))
vp.resid$variable = factor(vp.resid$variable,levels=colnames(vp2))

p = ggplot() +
	geom_violin(data=vp.resid,aes(variable,value * 100,fill=variable),scale='width') +
	geom_boxplot_jitter(data=vp.resid,aes(variable,value * 100),width=0.075,fill='#fcfcfc',outlier.colour='black',outlier.size=0.1,outlier.jitter.width=0,outlier.jitter.height=0) +
	scale_fill_manual(values=v.colors[levels(vp.resid$variable)]) +
	scale_y_continuous(limit=c(0,100),breaks=seq(0,100,25)) +
	theme_classic(base_size=24) +
	theme(legend.position='none',axis.text.x=element_text(angle=45,hjust=1,vjust=1),axis.title.x=element_blank()) +
	ylab('Variance explained (%)')
ggsave(p,file=paste0('figures/data_visualization_variance_partition_resid.pdf'),useDingbats=FALSE,height=7,width=7)

}

library(tidyverse)
library(ggplot2)
library(umap)
library(RColorBrewer)

set.seed(3)

a = umap(t(e.regressed), n_neighbors = 50, min_dist = 0.5)
a.umap = data.frame(as.data.frame(a$layout),meta[,c(full.covariates,'Region')])

p = ggplot(a.umap,aes_string('V1','V2',color='Region',shape=batch.variable)) +
	geom_point(size=1.25) +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=18) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file='figures/data_visualization_dimension_reduction_umap_batched.pdf',useDingbats=FALSE)

p = ggplot(a.umap,aes_string('V1','V2',color='Region')) +
	geom_point(size=1.25) +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=18) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file='figures/data_visualization_dimension_reduction_umap.pdf',useDingbats=FALSE)

p = ggplot(a.umap,aes_string('V1','V2',color='Region',shape='Region')) +
	geom_point(size=1.25) +
	scale_color_manual(values=region.colors) +
	scale_shape_manual(values=region.shapes) +
	theme_classic(base_size=18) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	guides(shape = guide_legend(override.aes = list(size = 3))) +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file='figures/data_visualization_dimension_reduction_umap_shapes.pdf',useDingbats=FALSE)

p = ggplot(a.umap,aes_string('V1','V2',color='Region',shape='Region')) +
	geom_point(size=1.5) +
	scale_color_manual(values=region.colors) +
	scale_shape_manual(values=region.shapes) +
	theme_classic(base_size=24) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	guides(shape = guide_legend(override.aes = list(size = 3))) +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file='figures/data_visualization_dimension_reduction_umap_shapes_presentation.pdf',useDingbats=FALSE,width=12,height=7)
ggsave(p,file='figures/data_visualization_dimension_reduction_umap_shapes_presentation1.pdf',useDingbats=FALSE,width=12,height=7)

p = ggplot(a.umap,aes_string('V2','V1',color='Region',shape='Region')) +
	geom_point(size=1.5) +
	scale_color_manual(values=region.colors) +
	scale_shape_manual(values=region.shapes) +
	theme_classic(base_size=24) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	guides(shape = guide_legend(override.aes = list(size = 3))) +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file='figures/data_visualization_dimension_reduction_umap_shapes_presentation2.pdf',useDingbats=FALSE,width=12,height=7)

p = ggplot(within(a.umap,{V1=-V1}),aes_string('V1','V2',color='Region',shape='Region')) +
	geom_point(size=1.5) +
	scale_color_manual(values=region.colors) +
	scale_shape_manual(values=region.shapes) +
	theme_classic(base_size=24) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	guides(shape = guide_legend(override.aes = list(size = 3))) +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file='figures/data_visualization_dimension_reduction_umap_shapes_presentation3.pdf',useDingbats=FALSE,width=12,height=7)

p = ggplot(within(a.umap,{V2=-V2}),aes_string('V2','V1',color='Region',shape='Region')) +
	geom_point(size=1.5) +
	scale_color_manual(values=region.colors) +
	scale_shape_manual(values=region.shapes) +
	theme_classic(base_size=24) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	guides(shape = guide_legend(override.aes = list(size = 3))) +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file='figures/data_visualization_dimension_reduction_umap_shapes_presentation4.pdf',useDingbats=FALSE,width=12,height=7)

library(viridis)
p = ggplot(a.umap,aes_string('V1','V2',color=predictor)) +
	geom_point(size=1.5) +
	scale_color_viridis(name=predictor.label,option='C') +
	theme_classic(base_size=24) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
#	guides(color = guide_legend(override.aes = list(size = 3))) +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_umap_',tolower(predictor.label),'_presentation.pdf'),useDingbats=FALSE,width=12,height=7)

p = ggplot(a.umap,aes_string('V1','V2',color=rank.ordinal.variable)) +
	geom_point(size=1.5) +
	scale_color_manual(name='Rank',values=c('#cbc9e2','#9e9ac8','#6a51a3')) +
	theme_classic(base_size=24) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	guides(color = guide_legend(override.aes = list(size = 3))) +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_umap_',tolower('rank_ordinal'),'_presentation.pdf'),useDingbats=FALSE,width=12,height=7)

p = ggplot(a.umap,aes_string('V1','V2',color=sex.variable)) +
	geom_point(size=1.5) +
	scale_color_manual(name='Sex',values=sex.colors) +
	theme_classic(base_size=24) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	guides(color = guide_legend(override.aes = list(size = 3))) +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_umap_',tolower('sex'),'_presentation.pdf'),useDingbats=FALSE,width=12,height=7)

p = ggplot(a.umap,aes_string('V1','V2',color='RIN')) +
	geom_point(size=1.5) +
	scale_color_viridis(name='RIN') +
	theme_classic(base_size=24) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_umap_',tolower('rin'),'_presentation.pdf'),useDingbats=FALSE,width=12,height=7)

a.umap$mapped_reads = log10(a.umap$reads_mapped)

p = ggplot(a.umap,aes_string('V1','V2',color='mapped_reads')) +
	geom_point(size=1.5) +
	scale_color_viridis(name=expression(log[10]('mapped reads'))) +
	theme_classic(base_size=24) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_umap_',tolower('mapped_reads'),'_presentation.pdf'),useDingbats=FALSE,width=12,height=7)

p = ggplot(a.umap,aes_string('V1','V2',color='Sequencing.batch')) +
	geom_point(size=1.5) +
	scale_color_manual(name='Batch',values=c('#7fc97f','#beaed4')) +
	theme_classic(base_size=24) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	guides(color = guide_legend(override.aes = list(size = 3))) +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_umap_',tolower('batch'),'_presentation.pdf'),useDingbats=FALSE,width=12,height=7)


# p = ggplot(a.umap,aes_string('V1','V2',color=rank.variable)) +
# 	geom_point(size=1.5) +
# 	scale_color_viridis(name='Rank') +
# 	theme_classic(base_size=24) +
# 	xlab('UMAP 1') +
# 	ylab('UMAP 2') +
# 	coord_fixed() +
# #	guides(shape = guide_legend(override.aes = list(size = 3))) +
# 	theme(axis.ticks=element_blank(),axis.text=element_blank())
# ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_umap_',tolower('rank_scaled'),'_presentation.pdf'),useDingbats=FALSE)


plot.gene = function(ensembl_gene_id,normalize=TRUE) {
	do.call(rbind,lapply(region.levels,function(x) {
		out = subset(meta,Region %in% x)
		if (normalize) {
			out$gene = scale(matrix(e.regressed[ensembl_gene_id,rownames(out)],ncol=1))
		} else {
			out$gene = e.keep[ensembl_gene_id,rownames(out)]
		}
		out
	}))
}

library(egg)

for (i in plot.ensembl.genes) {
	foo = plot.gene(i)
	levels(foo$Region) = gsub('ACCg','ACC',levels(foo$Region))
	
	p = ggplot(foo,aes_string(predictor,'gene',color='Region',fill='Region')) +
		geom_point() + geom_smooth(method=loess,se=TRUE,span=1) +
		facet_wrap(~Region,nrow=region.rows) +
		scale_color_manual(values=region.colors) +
		scale_fill_manual(values=region.colors) +
		scale_x_continuous(limits=c(0,20)) +
		scale_y_continuous(limits=c(-max(abs(foo$gene)),max(abs(foo$gene)))) +
		coord_cartesian(ylim=0.75 * c(-max(abs(foo$gene)),max(abs(foo$gene)))) +
		xlab(paste0(predictor.label,' (',tolower(predictor.units),')')) + ylab(expression(italic(Z)*'-score')) +
		theme_classic() +
		theme(legend.position='none',strip.background=element_blank())
	ggsave(p,file=paste0('figures/data_visualization_',predictor,'_by_gene_loess_',i,'.pdf'),useDingbats=FALSE,height=5)

	p = ggplot(foo,aes_string(predictor,'gene',color='Region',fill='Region')) +
		geom_point() + geom_smooth(method=lm,se=TRUE) +
		facet_wrap(~Region,nrow=region.rows) +
		scale_color_manual(values=region.colors) +
		scale_fill_manual(values=region.colors) +
		scale_x_continuous(limits=c(0,20)) +
		scale_y_continuous(limits=c(-max(abs(foo$gene)),max(abs(foo$gene)))) +
		coord_cartesian(ylim=0.75 * c(-max(abs(foo$gene)),max(abs(foo$gene)))) +
		xlab(paste0(predictor.label,' (',tolower(predictor.units),')')) + ylab(expression(italic(Z)*'-score')) +
		theme_classic() +
		theme(legend.position='none',strip.background=element_blank())
	ggsave(p,file=paste0('figures/data_visualization_',predictor,'_by_gene_',i,'.pdf'),useDingbats=FALSE,height=5)

	p = ggplot(foo,aes_string(predictor,'gene',color='Region',fill='Region')) +
		geom_point() + geom_smooth(method=lm,se=TRUE) +
		facet_wrap(~Region,nrow=region.rows) +
		scale_color_manual(values=region.colors) +
		scale_fill_manual(values=region.colors) +
		scale_x_continuous(limits=c(0,20)) +
		scale_y_continuous(limits=c(-max(abs(foo$gene)),max(abs(foo$gene)))) +
		coord_cartesian(ylim=0.75 * c(-max(abs(foo$gene)),max(abs(foo$gene)))) +
		xlab(paste0(predictor.label,' (',tolower(predictor.units),')')) + ylab(expression(italic(Z)*'-score')) +
		theme_article(base_size=16) +
		theme(legend.position='none',strip.background=element_blank())
	ggsave(p,file=paste0('figures/data_visualization_',predictor,'_by_gene_',i,'_article.pdf'),useDingbats=FALSE,height=5)

	foo = plot.gene(i,normalize=FALSE)
	levels(foo$Region) = gsub('ACCg','ACC',levels(foo$Region))
	
	p = ggplot(foo,aes_string(predictor,'gene',color='Region',fill='Region')) +
		geom_point() + geom_smooth(method=lm,se=TRUE) +
		facet_wrap(~Region,nrow=region.rows) +
		scale_color_manual(values=region.colors) +
		scale_fill_manual(values=region.colors) +
		scale_x_continuous(limits=c(0,20)) +
#		scale_y_continuous(limits=c(-max(abs(foo$gene)),max(abs(foo$gene)))) +
#		coord_cartesian(ylim=0.75 * c(-max(abs(foo$gene)),max(abs(foo$gene)))) +
		xlab(paste0(predictor.label,' (',tolower(predictor.units),')')) + ylab(expression(log[2]('TPM'))) +
		theme_classic() +
		theme(legend.position='none',strip.background=element_blank())
	ggsave(p,file=paste0('figures/data_visualization_',predictor,'_by_gene_',i,'_raw.pdf'),useDingbats=FALSE,height=5)

#	p = ggplot(subset(foo,Region %in% region.levels[1]),aes_string(predictor,'gene',color='Region',fill='Region')) +
#		geom_point(data=foo,aes_string(predictor,'gene'),alpha=0) + geom_point() + geom_smooth(method=lm,se=TRUE) +
#		facet_wrap(~Region,nrow=region.rows) +
#		scale_color_manual(values=region.colors) +
#		scale_fill_manual(values=region.colors) +
#		scale_x_continuous(limits=c(0,ceiling(foo[[predictor]]/5) * 5)) +
#		scale_y_continuous(limits=c(-max(abs(foo$gene)),max(abs(foo$gene)))) +
#		coord_cartesian(ylim=0.75 * c(-max(abs(foo$gene)),max(abs(foo$gene)))) +
#		xlab(predictor.label) + ylab(expression(log[2]('TPM'))) +
#		theme_classic() +
#		theme(legend.position='none',axis.text.y=element_blank(),axis.ticks.y=element_blank(),strip.background=element_blank())
#	ggsave(p,file=paste0('figures/data_visualization_',predictor,'_by_gene_',i,'_example.pdf'),useDingbats=FALSE,height=5)
}

# 3D umap
# a = umap(t(e.regressed), n_neighbors = 50, min_dist = 0.5, n_components = 3)

b = prcomp(cor(e.regressed))
b.pca = data.frame(as.data.frame(b$rotation)[,1:which(summary(b)$importance['Cumulative Proportion',] > 0.99)[1]],meta[,c(full.covariates,'Region')])

b.importance = gsub(' ','',paste0(format(summary(b)$importance['Proportion of Variance',1:which(summary(b)$importance['Cumulative Proportion',] > 0.99)[1]] * 100,digits=2),'%'))

p = ggplot(b.pca,aes_string('PC1','PC2',color='Region')) +
	geom_point(size=1.25) +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=18) +
	xlab(paste0('PC1 (',b.importance[1],')')) +
	ylab(paste0('PC2 (',b.importance[2],')')) +
#	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file='figures/data_visualization_dimension_reduction_pca.pdf',useDingbats=FALSE)

p = ggplot(b.pca,aes_string('PC1','PC2',color=predictor)) +
	geom_point(size=1.25) +
	scale_color_viridis(name=predictor.label,option='C') +
	theme_classic(base_size=18) +
	xlab(paste0('PC1 (',b.importance[1],')')) +
	ylab(paste0('PC2 (',b.importance[2],')')) +
#	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_pca_',tolower(predictor.label),'.pdf'),useDingbats=FALSE)

p = ggplot(b.pca,aes_string('PC1','PC2',color='ordinal.rank')) +
	geom_point(size=1.25) +
	scale_color_brewer(name='Ordinal rank',palette='Purples') +
	theme_classic(base_size=18) +
	xlab(paste0('PC1 (',b.importance[1],')')) +
	ylab(paste0('PC2 (',b.importance[2],')')) +
#	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_pca_','rank','.pdf'),useDingbats=FALSE)

for (i in 1:which(summary(b)$importance['Cumulative Proportion',] > 0.99)[1]) {
	p = ggplot(b.pca,aes_string(predictor,paste0('PC',i),color='Region')) +
		geom_point(size=1.25) +
		scale_color_manual(values=region.colors) +
		theme_classic(base_size=18) +
		scale_x_continuous(trans='log2') +
		xlab(paste0(predictor.label,' (',tolower(predictor.units),')')) +
		ylab(paste0('PC',i,' (',b.importance[i],')')) +
	#	coord_fixed() +
		theme()
	ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_',tolower(predictor.label),'_pc',i,'.pdf'),useDingbats=FALSE)
}

b.pca.long = tidyr::pivot_longer(b.pca,cols=paste0('PC',1:which(summary(b)$importance['Cumulative Proportion',] > 0.99)[1]))
b.pca.long = within(b.pca.long,{
	name=factor(name,levels=paste0('PC',1:which(summary(b)$importance['Cumulative Proportion',] > 0.99)[1]))
	pc.label=name
	levels(pc.label) = paste0(
		paste0('PC',1:which(summary(b)$importance['Cumulative Proportion',] > 0.99)[1]),
		' (',b.importance,')'
	)
})
levels(b.pca.long$Region) = gsub('ACCg','ACC',levels(b.pca.long$Region))

p = ggplot(b.pca.long) +
	geom_point(aes_string(predictor,'value',color='Region'),size=0.5) +
	geom_smooth(aes_string(predictor,'value'),method=lm,se=FALSE,color='#000000',size=0.5) +
	scale_color_manual(values=region.colors) +
	facet_wrap(~pc.label,ncol=1,strip.position='right') +
	theme_classic(base_size=18) +
	scale_x_continuous(trans='log2') +
	xlab(paste0(predictor.label,' (',tolower(predictor.units),')')) +
	ylab('PC loading') +
	coord_fixed(ratio=2) +
	theme(axis.text.y=element_blank(),strip.text.y.right = element_text(angle = 0))
ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_',tolower(predictor.label),'_pc_all.pdf'),useDingbats=FALSE)




p = ggplot(b.pca,aes_string('PC1','PC2',color='Region',shape='Region')) +
	geom_point(size=1.5) +
	scale_color_manual(values=region.colors) +
	scale_shape_manual(values=region.shapes) +
	theme_classic(base_size=24) +
	xlab(paste0('PC1 (',b.importance[1],')')) +
	ylab(paste0('PC2 (',b.importance[2],')')) +
	coord_fixed() +
	guides(shape = guide_legend(override.aes = list(size = 3))) +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file='figures/data_visualization_dimension_reduction_PC_shapes_presentation.pdf',useDingbats=FALSE,width=12,height=7)

p = ggplot(b.pca,aes_string('PC1','PC2',color=predictor)) +
	geom_point(size=1.5) +
	scale_color_viridis(name=predictor.label,option='C') +
	theme_classic(base_size=24) +
	xlab(paste0('PC1 (',b.importance[1],')')) +
	ylab(paste0('PC2 (',b.importance[2],')')) +
	coord_fixed() +
#	guides(color = guide_legend(override.aes = list(size = 3))) +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_PC_',tolower(predictor.label),'_presentation.pdf'),useDingbats=FALSE,width=12,height=7)

p = ggplot(b.pca,aes_string('PC1','PC2',color=rank.ordinal.variable)) +
	geom_point(size=1.5) +
	scale_color_manual(name='Rank',values=c('#cbc9e2','#9e9ac8','#6a51a3')) +
	theme_classic(base_size=24) +
	xlab(paste0('PC1 (',b.importance[1],')')) +
	ylab(paste0('PC2 (',b.importance[2],')')) +
	coord_fixed() +
	guides(color = guide_legend(override.aes = list(size = 3))) +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_PC_',tolower('rank_ordinal'),'_presentation.pdf'),useDingbats=FALSE,width=12,height=7)

p = ggplot(b.pca,aes_string('PC1','PC2',color=sex.variable)) +
	geom_point(size=1.5) +
	scale_color_manual(name='Sex',values=sex.colors) +
	theme_classic(base_size=24) +
	xlab(paste0('PC1 (',b.importance[1],')')) +
	ylab(paste0('PC2 (',b.importance[2],')')) +
	coord_fixed() +
	guides(color = guide_legend(override.aes = list(size = 3))) +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_PC_',tolower('sex'),'_presentation.pdf'),useDingbats=FALSE,width=12,height=7)

p = ggplot(b.pca,aes_string('PC1','PC2',color='RIN')) +
	geom_point(size=1.5) +
	scale_color_viridis(name='RIN') +
	theme_classic(base_size=24) +
	xlab(paste0('PC1 (',b.importance[1],')')) +
	ylab(paste0('PC2 (',b.importance[2],')')) +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_PC_',tolower('rin'),'_presentation.pdf'),useDingbats=FALSE,width=12,height=7)

b.pca$mapped_reads = log10(b.pca$reads_mapped)

p = ggplot(b.pca,aes_string('PC1','PC2',color='mapped_reads')) +
	geom_point(size=1.5) +
	scale_color_viridis(name=expression(log[10]('mapped reads'))) +
	theme_classic(base_size=24) +
	xlab(paste0('PC1 (',b.importance[1],')')) +
	ylab(paste0('PC2 (',b.importance[2],')')) +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_PC_',tolower('mapped_reads'),'_presentation.pdf'),useDingbats=FALSE,width=12,height=7)

p = ggplot(b.pca,aes_string('PC1','PC2',color='Sequencing.batch')) +
	geom_point(size=1.5) +
	scale_color_manual(name='Batch',values=c('#7fc97f','#beaed4')) +
	theme_classic(base_size=24) +
	xlab(paste0('PC1 (',b.importance[1],')')) +
	ylab(paste0('PC2 (',b.importance[2],')')) +
	coord_fixed() +
	guides(color = guide_legend(override.aes = list(size = 3))) +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file=paste0('figures/data_visualization_dimension_reduction_PC_',tolower('batch'),'_presentation.pdf'),useDingbats=FALSE,width=12,height=7)











# PCA for each region
pca = do.call(rbind,lapply(names(keep.genes),function(i) {
	genes = keep.genes[[i]]
	ids = rownames(subset(meta,Region == i))
	
	b = prcomp(cor(e.regressed[genes,ids]))
	b.pca = do.call(rbind,lapply(1:(length(keep.genes)*2),function(x) {
		data.frame(
			meta[ids,full.covariates],
			pc = x,
			pc.label = paste0('PC',formatC(x,width=2,flag=0)),
			loading = b$rotation[,paste0('PC',x)],
			region = i,
			prop.variance = summary(b)$importance['Proportion of Variance',paste0('PC',x)],
			cum.variance = summary(b)$importance['Cumulative Proportion',paste0('PC',x)],
			label.region = paste0(
				i, ' (',
				format(round(summary(b)$importance['Proportion of Variance',paste0('PC',x)] * 100,2),nsmall=2),
				'%)'
			),
			label.pc = paste0(
				paste0('PC',formatC(x,width=2,flag=0)), ' (',
				format(round(summary(b)$importance['Proportion of Variance',paste0('PC',x)] * 100,2),nsmall=2),
				'%)'
			),
			stringsAsFactors=FALSE
		)
	}))
	b.pca
}))

# Just print the first half
for (i in 1:length(keep.genes)) {
	p = ggplot(subset(pca,pc == i),aes_string(predictor,'loading')) +
		geom_point() +
		geom_smooth(method=lm,se=TRUE) +
		facet_wrap(~label.region,scales='free_y',nrow=region.rows) +
		theme_classic() +
		theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) + 
		xlab(predictor.label) +
		ylab(paste0('PC',i))
	ggsave(p,file=paste0('figures/data_visualization_',tolower(predictor.label),'_vs_pc',formatC(i,width=2,flag=0),'.pdf'),useDingbats=FALSE)
}
for (r in names(keep.genes)) {
	p = ggplot(subset(pca,region == r & pc <= 16),aes_string(predictor,'loading')) +
		geom_point() +
		geom_smooth(method=lm,se=TRUE) +
		facet_wrap(~label.pc,scales='free_y',nrow=region.rows) +
		theme_classic() +
		theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()) + 
		xlab(predictor.label) +
		ylab(paste0('PC'))
	ggsave(p,file=paste0('figures/data_visualization_',tolower(predictor.label),'_vs_pcs_region_',r,'.pdf'),useDingbats=FALSE)
}

# Look at age effect on each PC:

pca.split = split(pca,list(pca$pc,pca$region))

pca.predictor = do.call(rbind,lapply(pca.split,function(x) {
	model = lm(as.formula(paste0('loading~',paste(setdiff(full.covariates,batch.covariates),collapse='+'))),data=x)
	data.frame(unique(x[c('pc','pc.label','region','prop.variance','cum.variance','label.region','label.pc')]),
		coefficient = coef(model)[predictor],
		pval = coef(summary(model))[predictor,'Pr(>|t|)']
	)
}))

# Now split by region

pca.split = split(pca.predictor,pca.predictor$region)

pca.predictor = do.call(rbind,lapply(pca.split,function(x) {
	x[1:which(x$cum.variance >= 0.99)[1],]
}))

pca.predictor = within(pca.predictor,{
	Region = factor(region,levels=region.levels)
	pc.label = factor(pc.label,levels=sort(unique(pc.label)))
})

levels(pca.predictor$Region) = gsub('ACCg','ACC',levels(pca.predictor$Region))
p = ggplot(pca.predictor,aes(pc,-log10(pval),color=Region,size=prop.variance)) +
	geom_point() +
	geom_hline(data=reshape2::melt(0.05 / tapply(pca.predictor$pc,pca.predictor$Region,max)) %>% mutate(Region=Var1,pval=value),aes(yintercept=-log10(pval)),linetype=3,size=0.2) +
	geom_hline(aes(yintercept=-log10(0.05)),size=0.5,linetype=2) +
	facet_wrap(~Region,ncol=region.rows) +
	scale_size_continuous(range=c(0,3),guide=FALSE) +
	scale_color_manual(values=region.colors) +
	scale_x_continuous(breaks=1:max(pca.predictor$pc),labels=ifelse(1:max(pca.predictor$pc) %% 5 == 0 | 1:max(pca.predictor$pc) == 1,1:max(pca.predictor$pc),'')) +
	theme_classic() +
	theme(strip.text=element_blank()) +
	guides(color = guide_legend(override.aes = list(size = 3))) +
	xlab('PC') + ylab(parse(text=paste0('-log[10]*"("*italic(P)[',tolower(predictor.label),']*")"')))
gt = ggplot_gtable(ggplot_build(p))
gt$layout$clip[grepl('panel',gt$layout$name)] = 'off'

pdf(file=paste0('figures/data_visualization_',tolower(predictor.label),'_vs_pc_pvals.pdf'),height=5,useDingbats=FALSE)
	grid::grid.draw(gt)
dev.off()
# ggsave(p,file=paste0('figures/data_visualization_',tolower(predictor.label),'_vs_pc_pvals.pdf'),height=5,useDingbats=FALSE)

# K-means clustering

library(cluster)

# Scale across entire matrix (but not per-gene to preserve relative scales between genes)
e.k = matrix(as.numeric(scale(as.numeric(t(e.regressed)))),nrow=ncol(e.regressed),ncol=nrow(e.regressed),dimnames=list(colnames(e.regressed),rownames(e.regressed)))

# Scale per-gene
# e.k = scale(t(e.regressed))

# Set number of iterations
niter = length(keep.genes)

# Initialize results vectors
wss = numeric(niter)
wss[1] = (nrow(e.k) - 1) * sum(apply(e.k,2,var))

avg.sil.values = numeric(niter - 1)

km.results = vector('list',niter)

for (i in 2:niter) {
	message('Now running k = ',i)
	km.results[[i]] = kmeans(e.k,centers=i,nstart=25)
	km.result = km.results[[i]]

	# Save elbow scores (within-cluster sum of squares)
	wss[i] = sum(km.result$withinss)

	# Save silhouette score
	avg.sil.values[i-1] = mean(silhouette(km.result$cluster, dist(e.k))[,3])

	km.clusters = data.frame(meta,cluster=km.result$cluster)
	p = ggplot(data.frame(a.umap,cluster=km.clusters$cluster),aes(V1,V2,color=factor(cluster))) +
		geom_point(size=2) +
		theme_classic(base_size=18) +
		xlab('UMAP 1') +
		ylab('UMAP 2') +
		coord_fixed() +
		theme(axis.ticks=element_blank(),axis.text=element_blank())
	if (i > 8) {
		p = p + scale_color_manual(name='Cluster',values=colorRampPalette(brewer.pal(8, 'Set1'))(i))
	} else {
		p = p + scale_color_brewer(name='Cluster',palette='Set1')
	}
	ggsave(p,file=paste0('figures/data_visualization_umap_kmeans_samples_clusters_k_',formatC(i,width=2,flag=0),'.pdf'),useDingbats=FALSE)
}

best.sil = (2:niter)[which.max(avg.sil.values)]
best.elb = which.min(unlist(lapply(1:niter,function(x) {
	d = diff(wss[(x-1):(x+1)])
	d[2]/d[1]
})))

pdf(file='figures/data_visualization_kmeans_samples_elbow.pdf',useDingbats=FALSE)
	plot(1:niter, wss, type='b', pch = 19, frame = FALSE, xlab='Number of clusters', ylab='Within groups sum of squares',xlim=c(1,niter),xaxt = 'n')
	points(wss[best.elb] ~ best.elb, col='red', cex=2)
	axis(side = 1, at=1:niter)
dev.off()

pdf(file='figures/data_visualization_kmeans_samples_silhouette.pdf',useDingbats=FALSE)
	plot(2:niter, avg.sil.values, type = 'b', pch = 19, frame = FALSE,  xlab = 'Number of clusters', ylab = 'Average silhoettes', xlim=c(1,niter),xaxt = 'n')
	points(avg.sil.values[best.sil-1] ~ best.sil, col='red', cex=2)
	axis(side = 1, at=1:niter)
dev.off()

# Mapping percentage
p = ggplot(km.clusters,aes(factor(cluster),perc_mapped)) +
	geom_boxplot(outlier.shape=NA) +
	geom_jitter() +
	xlab('Cluster') +
	ylab('Proportion mapped') +
	theme_classic()
ggsave(p,file='figures/data_visualization_sequencing_stats_cluster_vs_map_perc.pdf',useDingbats=FALSE)

# Mapped reads
p = ggplot(km.clusters,aes(factor(cluster),reads_mapped)) +
	geom_boxplot(outlier.shape=NA) +
	geom_jitter() +
	xlab('Cluster') +
	ylab('Reads mapped') +
	theme_classic()
ggsave(p,file='figures/data_visualization_sequencing_stats_cluster_vs_reads_mapped.pdf',useDingbats=FALSE)

# Total reads
p = ggplot(km.clusters,aes(factor(cluster),reads)) +
	geom_boxplot(outlier.shape=NA) +
	geom_jitter() +
	xlab('Cluster') +
	ylab('Total reads') +
	theme_classic()
ggsave(p,file='figures/data_visualization_sequencing_stats_cluster_vs_reads_total.pdf',useDingbats=FALSE)

# Duplicate rate
p = ggplot(km.clusters,aes(factor(cluster),perc_duplicate)) +
	geom_boxplot(outlier.shape=NA) +
	geom_jitter() +
	xlab('Cluster') +
	ylab('Duplication rate') +
	theme_classic()
ggsave(p,file='figures/data_visualization_sequencing_stats_cluster_vs_dup_perc.pdf',useDingbats=FALSE)

# Quality mean
p = ggplot(km.clusters,aes(factor(cluster),qual_mean)) +
	geom_boxplot(outlier.shape=NA) +
	geom_jitter() +
	xlab('Cluster') +
	ylab('Mean quality') +
	theme_classic()
ggsave(p,file='figures/data_visualization_sequencing_stats_cluster_vs_qual_mean.pdf',useDingbats=FALSE)

# Dendrogram

library(dendextend)

h.cluster = hclust(dist(e.k))
h.cluster.tree = as.dendrogram(h.cluster)

library2region = meta$Region
names(library2region) = meta$Library

tree.colors = region.colors
names(tree.colors) = names(keep.genes)

colLab = function(n) {
    if (is.leaf(n)) {
        a = attributes(n)
        labCol = as.character(tree.colors[library2region[a$label]])
        attr(n, 'nodePar') = list(lab.cex = 0.1, lab.col = labCol, col = labCol, pch=19, cex=0.05)
        attr(n, 'edgePar') = list(col = labCol, lwd=0.5)
    }
    n
}

h.cluster.tree = dendrapply(h.cluster.tree, colLab)

edgePar = list(lwd=0.5)

pdf(file='figures/data_visualization_hclust_dendrogram_all.pdf',useDingbats=FALSE,height=10)
	plot(h.cluster.tree,edgePar=edgePar,horiz=TRUE,leaflab = 'none', yaxt='n')
dev.off()

pdf(file='figures/data_visualization_hclust_dendrogram_all_labels.pdf',useDingbats=FALSE,height=10)
	plot(h.cluster.tree,edgePar=edgePar,horiz=TRUE, yaxt='n')
dev.off()

# Repeat on mean expression across regions
region.libraries = with(meta,split(Library,Region))

e.r = do.call(rbind,lapply(names(keep.genes),function(x) {
	matrix(colMeans(e.k[region.libraries[[x]],]),nrow=1,dimnames=list(x,colnames(e.k)))
}))

r.cluster = hclust(dist(e.r))
r.cluster.tree = as.dendrogram(r.cluster)

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

r.cluster.tree = dendrapply(r.cluster.tree, colReg)

pdf(file='figures/data_visualization_hclust_dendrogram_regions.pdf',useDingbats=FALSE,height=10)
	plot(r.cluster.tree,edgePar=edgeReg,horiz=TRUE,yaxt='n')
dev.off()

# Bootstrap libraries

library(ape)

r.cluster.bootstrap = mclapply(1:1000,function(i) {
	e.x = do.call(rbind,lapply(names(keep.genes),function(x) {
		matrix(colMeans(e.k[sample(region.libraries[[x]],replace=TRUE),]),nrow=1,dimnames=list(x,colnames(e.k)))
	}))
	as.phylo(hclust(dist(e.x)))
},mc.cores=n.cores)

# Visualize uncertainty
library(phangorn)

r.cluster.bootstrap.multiphylo = do.call(c,r.cluster.bootstrap)

pdf(file='figures/data_visualization_hclust_dendrogram_regions_bootstrap.pdf',useDingbats=FALSE,height=10)
	densiTree(
		r.cluster.bootstrap.multiphylo,
		alpha=0.005,
		consensus=as.phylo(r.cluster.tree),
		direction='rightwards',
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
	consensus = as.phylo(r.cluster.tree)
	xy = ape::plotPhyloCoor(consensus, tip.height = 1:length(keep.genes), direction = 'rightwards')
	xx = xy[,1]; yy = xy[,2]; xx = xx/max(xx); xx = xx + (1 - max(xx))
	ape::cladogram.plot(consensus$edge, xx, yy, edge.color = '#0000ff', edge.width = 1, edge.lty = 1)
dev.off()

pdf(file='figures/data_visualization_hclust_dendrogram_regions_bootstrap_noconsensus.pdf',useDingbats=FALSE,height=10)
	densiTree(
		r.cluster.bootstrap.multiphylo,
		alpha=0.005,
		consensus=as.phylo(r.cluster.tree),
		direction='rightwards',
		scaleX=TRUE,
		col='#000000',
		width=1, lty=1, cex=2, font=1,
		tip.color=rapply(r.cluster.tree,function(n) {
			if (is.leaf(n)) attr(n, 'edgePar')[['col']]
		}), srt=0, adj=0,
		label.offset=0,
		scale.bar = FALSE,
		jitter=list(amount = 0.02, random = TRUE)
	)
	# Add the consensus manually
#	consensus = as.phylo(r.cluster.tree)
#	xy = ape::plotPhyloCoor(consensus, tip.height = 1:length(keep.genes), direction = 'rightwards')
#	xx = xy[,1]; yy = xy[,2]; xx = xx/max(xx); xx = xx + (1 - max(xx))
#	ape::cladogram.plot(consensus$edge, xx, yy, edge.color = '#0000ff', edge.width = 1, edge.lty = 1)
dev.off()

pdf(file='figures/data_visualization_hclust_dendrogram_regions_bootstrap_withconsensus.pdf',useDingbats=FALSE,height=10)
	densiTree(
		r.cluster.bootstrap.multiphylo,
		alpha=0.005,
		consensus=as.phylo(r.cluster.tree),
		direction='rightwards',
		scaleX=TRUE,
		col='#000000',
		width=1, lty=1, cex=2, font=1,
		tip.color=rapply(r.cluster.tree,function(n) {
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
	ape::cladogram.plot(consensus$edge, xx, yy, edge.color = '#0000ff', edge.width = 1, edge.lty = 1)
dev.off()

# # Can also use the java program DensiTree
# 
# if (!dir.exists('results')) dir.create('results')
# write.nexus(r.cluster.bootstrap,file='results/gene_expression_clustering_region_bootstraps.nex')
# 
# write.tree(as.phylo(r.cluster.tree),file='results/gene_expression_clustering_region_consensus.nwk')
# 
# # Create binaries folder if doesn't exist
# if (!file.exists('bin/DensiTree.jar')) {
# 	if (!dir.exists('bin')) dir.create('bin')
# 	download.file('https://github.com/rbouckaert/DensiTree/releases/download/v2.2.3/DensiTree.v2.2.7.jar',destfile='bin/DensiTree.jar')
# }

# Save best tree and multiTrees for later
saveRDS(as.phylo(r.cluster.tree),'checkpoints/dendrogram_gene_expression_best_tree.rds')

saveRDS(r.cluster.bootstrap.multiphylo,file='checkpoints/dendrogram_gene_expression_bootstrap_tree.rds')

saveRDS(e.regressed,file='checkpoints/regressed_expression_matrix.rds')
message('Successfully wrote tree nexus file!\n\nEdit it by running the following:\njava -jar bin/DensiTree.jar results/gene_expression_clustering_region_bootstraps.nex')



save(list=ls(),file='checkpoints/checkpoint_data_visualization.RData')