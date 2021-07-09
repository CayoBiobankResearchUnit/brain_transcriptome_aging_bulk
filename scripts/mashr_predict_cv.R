#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

if (!length(arguments)) {
	stop('At least 1 argument must be provided')
} else if (length(arguments) > 4) {
	stop('Only 4 arguments supported, ',length(arguments),' provided')
} else if (length(arguments) == 2) {
	this.i = as.integer(arguments[1])
	kf = as.integer(arguments[2])
	do.social = FALSE
	which.social = NULL
	message('No argument social provided, assuming social = FALSE')
} else if (length(arguments) == 4) {
	this.i = as.integer(arguments[1])
	kf = as.integer(arguments[2])
	do.social = as.logical(as.integer(arguments[3]))
	which.social = as.integer(arguments[3])
	this.cutoff = if (do.social) as.numeric(arguments[4]) else NULL
} else {
	this.i = as.integer(arguments[1])
	kf = 0
	do.social = FALSE
	which.social = NULL
	message('No argument k provided, assuming leave-one-out cross validation')
}

source('scripts/_include_options.R')

social.prefix = if (do.social) paste0('social_',if (which.social > 1) 'universal_' else '',formatC(this.cutoff),'_') else ''

e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')
keep.genes = readRDS('checkpoints/keep_genes.rds')
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')
mash.best.genes = readRDS(paste0('checkpoints/',social.prefix,'mashr_predictions_best_genes.rds'))

# Put together k and z matrices

animals = sort(unique(meta$Individual))

# Perform k-fold partitioning

# Set seed is essential so the partitioning is consistent across replicates
set.seed(seed)

# If leave-one out, k is equal to the number of animals
if (!kf) kf = length(animals)
kf.partitions = data.frame(animals,i=sample(rep(1:kf,times=ceiling(length(animals)/kf)),length(animals)),stringsAsFactors=FALSE)

training.animals = with(kf.partitions,animals[i != this.i])
evaluation.animals = with(kf.partitions,animals[i == this.i])

message('Training on ',length(training.animals), ' animals and evaluating on ',length(evaluation.animals),' animals.')

# Back up animals and meta objects
animals.backup = animals
meta.backup = meta

# Now filter to only training animals
animals = training.animals
meta = subset(meta,Individual %in% training.animals)

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

# Initialize output
out = vector('list',length(region.levels))
names(out) = region.levels

for (i in 1:length(region.levels)) {
	message('Now analyzing region ',region.levels[i])
	m = droplevels(subset(meta,Region %in% region.levels[i]))
	z.this = z[rownames(z) %in% m$Library,]
	# Subset expression matrix to only genes passing filters (in keep.genes) and including only libraries of this region
	e.this = e.keep[keep.genes[[region.levels[i]]],colnames(e.keep) %in% m$Library]

	# Drop covariates if they are uniform across dataset
	c.this = model.covariates[apply(m[,model.covariates],2,function(x) length(unique(x))) > 1]

	design = model.matrix(as.formula(paste('~',paste(c.this,collapse=' + '))), data=m)

	clus = makeForkCluster(n.cores)
	registerDoParallel(cores=n.cores)  
	clusterExport(clus,varlist=c('e.this','k','z.this','design'),envir=environment())

	out[[region.levels[i]]] = t(parApply(clus,e.this,1,function(y) {
		require(EMMREML)

		emma=emmreml(y = y,X = design,Z = z.this,K = k,varbetahat = T,varuhat = T,PEVuhat = T,test = T)

		p = emma$pvalbeta
		varb = emma$varbetahat
		b = emma$betahat

		c(b,varb,p[,"none"])
	}))

	stopCluster(clus)

	colnames(out[[region.levels[i]]])[(ncol(design) * 0 + 1):(ncol(design) * 1)] = paste('beta',colnames(design),sep='.')
	colnames(out[[region.levels[i]]])[(ncol(design) * 1 + 1):(ncol(design) * 2)] = paste('bvar',colnames(design),sep='.')
	colnames(out[[region.levels[i]]])[(ncol(design) * 2 + 1):(ncol(design) * 3)] = paste('pval',colnames(design),sep='.')

}

regions.dimnames = list(genes = Reduce(union,lapply(out,rownames)), outputs = Reduce(union,lapply(out,colnames)), regions = region.levels)
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

emma.results = regions.array

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# Run adaptive shrinkage
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #

Bhat = emma.results[,paste('beta',predictor,sep='.'),]
Shat = sqrt(emma.results[,paste('bvar',predictor,sep='.'),])

# For missing data, set beta to 0 and standard error to 100
# See https://github.com/stephenslab/mashr/issues/17#issuecomment-330628338
Bhat[is.na(Bhat)] = 0
Shat[is.na(Shat)] = 1000

library(mashr)

# Create the mashr data object
mash.data = mash_set_data(Bhat,Shat)

# Compute canonical covariance matrices
U.c = cov_canonical(mash.data)  

m.1by1 = mash_1by1(mash.data)
strong.subset = get_significant_results(m.1by1, thresh = 0.05)

# Get random subset (randomly choose half of all genes)
set.seed(seed)
random.subset = sample(1:nrow(Bhat),ceiling(nrow(Bhat)/2))

# Set temporary objects in order to estimate null correlation structure
temp = mash_set_data(Bhat[random.subset,],Shat[random.subset,])
temp.U.c = cov_canonical(temp)
Vhat = estimate_null_correlation(temp,temp.U.c)
rm(list=c('temp','temp.U.c'))

mash.random = mash_set_data(Bhat[random.subset,],Shat[random.subset,],V=Vhat)
mash.strong = mash_set_data(Bhat[strong.subset,],Shat[strong.subset,], V=Vhat)

# Perform PCA and extreme deconvolution to obtain data-driven covariances
U.pca = cov_pca(mash.strong,5)
U.ed = cov_ed(mash.strong, U.pca)

# Fit mash model
U.c = cov_canonical(mash.random)

m.r = mash(mash.random, Ulist = c(U.ed,U.c), outputlevel = 1)

m = mash(mash.data, g=get_fitted_g(m.r), fixg=TRUE)

mash.results = m

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# Predict ages
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #

# Restore meta
meta = meta.backup

mash.beta = get_pm(mash.results)
mash.lfsr = get_lfsr(mash.results)
mash.sbet = mash.beta / get_psd(mash.results)

mash.gene.beta = mash.beta[mash.best.genes,]
mash.gene.sbet = mash.sbet[mash.best.genes,]

e.scaled = t(apply(e.keep,1,scale))
colnames(e.scaled) = colnames(e.keep)

e.scaled = e.scaled[,with(meta,Library[Individual %in% evaluation.animals])]

mash.gene.predictions.merged  = unlist(lapply(colnames(e.scaled),function(x) sum(e.scaled[mash.best.genes,x] * rowMeans(mash.gene.beta))))
mash.gene.predictions.single  = unlist(lapply(colnames(e.scaled),function(x) sum(e.scaled[mash.best.genes,x] * mash.gene.beta[,with(meta,as.character(Region)[Library==x])])))

mash.gene.predictions.combined = rbind(data.frame(
	x = colnames(e.scaled),
	p = mash.gene.predictions.merged,
	i = this.i,
	k = kf,
	type = 'region-agnostic',
	stringsAsFactors=FALSE
),
data.frame(
	x = colnames(e.scaled),
	p = mash.gene.predictions.single,
	i = this.i,
	k = kf,
	type = 'region-specific',
	stringsAsFactors=FALSE
))

saveRDS(mash.gene.predictions.combined,file=paste0('checkpoints/',social.prefix,'predictions_cv_k_',formatC(kf,width=2,flag=0),'_i_',formatC(this.i,width=3,flag=0),'.rds'))
