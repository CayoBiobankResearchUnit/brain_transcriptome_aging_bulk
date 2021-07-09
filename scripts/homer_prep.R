#!/usr/bin/env Rscript

mash.results = readRDS('checkpoints/mashr_results.rds')

library(mashr)

mash.lfsr = get_lfsr(mash.results)
mash.beta = get_pm(mash.results)

up.genes = names(which(unlist(lapply(rownames(mash.beta),function(i) {
	out = any(mash.lfsr[i,] < 0.05 & mash.beta[i,] > 0)
	names(out) = i
	out
}))))

down.genes = names(which(unlist(lapply(rownames(mash.beta),function(i) {
	out = any(mash.lfsr[i,] < 0.05 & mash.beta[i,] < 0)
	names(out) = i
	out
}))))

background.up = setdiff(rownames(mash.beta),up.genes)
background.down = setdiff(rownames(mash.beta),down.genes)

# keep.genes = readRDS('checkpoints/keep_genes.rds')
# 
# all.genes = Reduce(union,keep.genes)
# 
# best.genes = readRDS('checkpoints/predictions_default_genes.rds')
# 
# background.genes = setdiff(all.genes,best.genes)

dir.create('results',showWarnings=FALSE)

write(up.genes,sep='\n',file='results/mmul10_best_genes_1.txt')
write(down.genes,sep='\n',file='results/mmul10_best_genes_2.txt')
write(background.up,sep='\n',file='results/mmul10_background_genes_1.txt')
write(background.down,sep='\n',file='results/mmul10_background_genes_2.txt')
