#!/usr/bin/env Rscript

source('scripts/_include_options.R')

cell.proportions = readRDS('checkpoints/cayo_bulkbrain_cell_proportions.rds')

e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')
keep.genes = readRDS('checkpoints/keep_genes.rds')
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')

pc = prcomp(cell.proportions)

# Keep smallest number of PCs that together explain >90% of the variation
pc.keep = 1:which(summary(pc)$importance['Cumulative Proportion',] > 0.9)[1]

cell.proportions.pc = pc$x[,pc.keep]

model.covariates = c(model.covariates,paste0('PC',pc.keep))

meta = data.frame(meta,cell.proportions.pc[rownames(meta),])

# Put together k and z matrices

animals = sort(unique(meta$Individual))

# Calculate k matrix (pairwise relatedness)

k = matrix(0,nrow=length(animals),ncol=length(animals),dimnames=list(animals,animals))

if (length(list.files(path='data',pattern='cayo_brain_lcmlkin_results\\.txt'))) {
	# If kinship output is found
	kinship = read.delim('data/cayo_brain_lcmlkin_results.txt',stringsAsFactors=FALSE)
	kinship = kinship[kinship$Ind1 %in% animals & kinship$Ind2 %in% animals,]
	i = as.matrix(kinship[,c('Ind1','Ind2')])
	k[i[,1:2]] = kinship$pi_HAT
	k[i[,2:1]] = kinship$pi_HAT
	diag(k) = 1
} else {
	# If kinship output is not found, substitute an identity matrix for testing purposes
	diag(k) = 1
}

# Calculate z matrix (matching libraries to genotype)
z = matrix(0,nrow=nrow(meta),ncol=nrow(k))
rownames(z) = meta$Library
colnames(z) = rownames(k)
i = as.matrix(meta[,c('Library','Individual')])
z[i] = 1

library(parallel)
library(doParallel)

# Within region

regions = names(keep.genes)

# Initialize output
out = vector('list',length(regions))
names(out) = regions

# Would like to add rank.scaled and reads_mapped but right now that is leading to a computational singularity problem
# Design model covariates

for (i in 1:length(regions)) {
	message('Now analyzing region ',regions[i])
	m = droplevels(subset(meta,Region %in% regions[i]))
	z.this = z[rownames(z) %in% m$Library,]
	# Subset expression matrix to only genes passing filters (in keep.genes) and including only libraries of this region
	e.this = e.keep[keep.genes[[regions[i]]],colnames(e.keep) %in% m$Library]

	# Drop covariates if they are uniform across dataset
	c.this = model.covariates[apply(m[,model.covariates],2,function(x) length(unique(x))) > 1]

	design = model.matrix(as.formula(paste('~',paste(c.this,collapse=' + '))), data=m)

	clus = makeCluster(n.cores)
	registerDoParallel(cores=n.cores)  
	clusterExport(clus,varlist=c('e.this','k','z.this','design'),envir=environment())

	out[[regions[i]]] = t(parApply(clus,e.this,1,function(y) {
		require(EMMREML)

		emma=emmreml(y = y,X = design,Z = z.this,K = k,varbetahat = T,varuhat = T,PEVuhat = T,test = T)

		p = emma$pvalbeta
		varb = emma$varbetahat
		b = emma$betahat

		c(b,varb,p[,"none"])
	}))

	stopCluster(clus)

	colnames(out[[regions[i]]])[(ncol(design) * 0 + 1):(ncol(design) * 1)] = paste('beta',colnames(design),sep='.')
	colnames(out[[regions[i]]])[(ncol(design) * 1 + 1):(ncol(design) * 2)] = paste('bvar',colnames(design),sep='.')
	colnames(out[[regions[i]]])[(ncol(design) * 2 + 1):(ncol(design) * 3)] = paste('pval',colnames(design),sep='.')

}

regions.dimnames = list(genes = Reduce(union,lapply(out,rownames)), outputs = Reduce(union,lapply(out,colnames)), regions = regions)
regions.dim = unlist(lapply(regions.dimnames,length))
regions.numeric = numeric(Reduce(`*`,regions.dim))
regions.numeric[!regions.numeric] = NA
regions.array = array(unname(regions.numeric),dim=unname(regions.dim),dimnames=unname(regions.dimnames))

for (i in 1:length(out)) {
	foo = reshape2::melt(out[[i]])
	foo$Var3 = names(out)[i]
	j = as.matrix(foo[,paste0('Var',1:3)])
	regions.array[j] = foo$value
	rm(foo)
}

saveRDS(regions.array,file='checkpoints/cellprop_emma_results.rds')
