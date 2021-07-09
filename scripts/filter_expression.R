#!/usr/bin/env Rscript

source('scripts/_include_options.R')

meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')

txi.gene = readRDS('checkpoints/kallisto_genes.rds')

# Do not include libraries that have been filtered out of animal metadata
abundance = txi.gene$abundance[,rownames(meta)]
counts = txi.gene$counts[,rownames(meta)]

# Filter to only autosomal genes + X
if (ignore.checkpoints || !file.exists('checkpoints/rnaseq_gene_metadata.rds')) {
	library(biomaRt)

	mmul = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
		dataset='mmulatta_gene_ensembl')

	mmul.genes = getBM(
		attributes=c('ensembl_gene_id','chromosome_name','gene_biotype'),
		filters = 'ensembl_gene_id',
		values = rownames(counts),
		mart = mmul)
	saveRDS(mmul.genes,file='checkpoints/rnaseq_gene_metadata.rds')
} else {
	message('Checkpoint found!\nLoading gene metadata from file.')

	mmul.genes = readRDS('checkpoints/rnaseq_gene_metadata.rds')
}

autosome.x.genes = subset(mmul.genes,chromosome_name %in% c(1:20,'X'))$ensembl_gene_id

abundance = abundance[rownames(abundance) %in% autosome.x.genes,]
counts = counts[rownames(counts) %in% autosome.x.genes,]

mmul.genes = subset(mmul.genes,ensembl_gene_id %in% autosome.x.genes)

# Check biotype distribution by sample

biotype.values = list(
	'IG_gene' = c('IG_C_gene','IG_D_gene','IG_J_gene','IG_V_gene'),
	'lncRNA' = c('antisense_RNA','lincRNA','lncRNA'),
	'ncRNA' = c('miRNA','misc_RNA','Mt_rRNA','Mt_tRNA','ncRNA','piRNA','pre_miRNA','ribozyme','rRNA','scaRNA','scRNA','siRNA','snRNA','snoRNA','sRNA','tRNA','vaultRNA'),
	'Pseudogene' = c('IG_C_pseudogene','IG_J_pseudogene','IG_pseudogene','IG_V_pseudogene','polymorphic_pseudogene','processed_pseudogene','pseudogene','rRNA_pseudogene','TR_J_pseudogene','TR_V_pseudogene','transcribed_processed_pseudogene','transcribed_unitary_pseudogene','transcribed_unprocessed_pseudogene','translated_processed_pseudogene','translated_unprocessed_pseudogene','unitary_pseudogene','unprocessed_pseudogene'),
	'Protein_coding' = c('protein_coding'),
	'TEC' = c('TEC'),
	'TR_gene' = c('TR_C_gene','TR_D_gene','TR_J_gene','TR_V_gene')
)

biotype.proportions = apply(do.call(rbind,lapply(with(mmul.genes,split(ensembl_gene_id,gene_biotype)),function(x) {
	colSums(abundance[x,])
})),2,function(x) x/sum(x))

# Merge like-categories
bioclass.proportions = do.call(rbind,lapply(biotype.values,function(i) {
	j = i[i %in% rownames(biotype.proportions)]
	out = colSums(matrix(biotype.proportions[j,],ncol=ncol(biotype.proportions),dimnames=list(j,colnames(biotype.proportions))),na.rm=TRUE)
	if (!sum(out)) {
		matrix(0,nrow=0,ncol(biotype.proportions),dimnames=list(character(0),colnames(biotype.proportions)))
	} else {
		out
	}
}))

library(reshape2)
biotype.proportions.df = melt(biotype.proportions)
biotype.proportions.df = merge(biotype.proportions.df,meta,by.x='Var2',by.y='Library',all.x=TRUE)
bioclass.proportions.df = melt(bioclass.proportions)
bioclass.proportions.df = merge(bioclass.proportions.df,meta,by.x='Var2',by.y='Library',all.x=TRUE)

library(ggplot2)
library(RColorBrewer)
p = ggplot(biotype.proportions.df,aes_string('Var2','value',fill=batch.variable)) +
	geom_bar(stat='identity') +
	scale_fill_brewer(palette='Dark2') +
	facet_wrap(as.formula(paste('~ Var1',sex.variable,sep=' + ')),ncol=2,scales='free') +
	xlab('Library') +
	ylab('Proportion') +
	theme_classic() +
	theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5,size=2))
ggsave(p,file='figures/data_qc_gene_biotypes.pdf',height=20,width=20)

p = ggplot(bioclass.proportions.df,aes_string('Var2','value',fill=batch.variable)) +
	geom_bar(stat='identity') +
	scale_fill_brewer(palette='Dark2') +
	facet_wrap(as.formula(paste('~ Var1',sex.variable,sep=' + ')),ncol=2,scales='free') +
	xlab('Library') +
	ylab('Proportion') +
	theme_classic() +
	theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5,size=2))
ggsave(p,file='figures/data_qc_gene_bioclasses.pdf',height=10,width=20)

# Filter to genes with mean >= TPM cutoff in any region

regions = sort(unique(meta$Region))

samples.by.region = split(meta$Library,meta$Region)

keep.genes = lapply(samples.by.region,function(x) {
	names(which(rowMeans(abundance[,x]) >= tpm.cutoff))
})

keep = Reduce(union,keep.genes)
# keep = names(which(rowMeans(abundance) >= 2))

library(limma)

v = voom(counts)

e = v$E[keep,rownames(meta)]

saveRDS(e,file='checkpoints/filtered_expression_matrix.rds')
saveRDS(keep.genes,file='checkpoints/keep_genes.rds')
