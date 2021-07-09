#!/usr/bin/env Rscript

library(tximport)
library(rhdf5)

samples = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')

files = file.path('kallisto',samples,'abundance.h5')
names(files) = samples

## Import hd5 kallisto-mapped files
txi.kallisto = tximport(files, type = 'kallisto', txOut = TRUE)

if (ignore.checkpoints || !file.exists('checkpoints/rnaseq_transcripts_genes.rds')) {
	library(biomaRt)
	mmul = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='mmulatta_gene_ensembl')

	# Fetch transcript-to-gene join table
	tx2gene = getBM(
		attributes = c('ensembl_transcript_id_version','ensembl_gene_id'),
		filters = 'ensembl_transcript_id_version',
		values = rownames(txi.kallisto$abundance),
		mart = mmul)

	saveRDS(tx2gene,file='checkpoints/rnaseq_transcripts_genes.rds')
} else {
	message('Checkpoint found!\nLoading transcript metadata from file.')

	tx2gene = readRDS('checkpoints/rnaseq_transcripts_genes.rds')
}

# Summarize mapped data to gene level
txi.gene = summarizeToGene(txi.kallisto, tx2gene)

saveRDS(txi.gene,file='checkpoints/kallisto_genes.rds')