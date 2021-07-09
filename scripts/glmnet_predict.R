#!/usr/bin/env Rscript

source('scripts/_include_options.R')

meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')

arguments = commandArgs(trailing=TRUE)

if (length(arguments) == 2) {
	this.i = as.integer(arguments[1])
	this.index = as.integer(arguments[2])
	do.social = FALSE
	which.social = NULL
	this.cutoff = NULL
} else if (length(arguments == 4)) {
	this.i = as.integer(arguments[1])
	this.index = as.integer(arguments[2])
	do.social = as.logical(as.integer(arguments[3]))
	which.social = as.integer(arguments[3])
	this.cutoff = as.numeric(arguments[4])
}

this.region = if (this.i == 0) 'all' else region.levels[this.i]

social.prefix = if (do.social) paste0('social_',if (which.social > 1) 'universal_' else '',formatC(this.cutoff),'_') else ''

meta.region = subset(meta,Region == this.region | this.region == 'all')
this.animal = sort(unique(meta$Individual))[this.index]

# Do nothing if the integer surpasses the number of individuals
if (this.index <= nrow(meta.region)) {

this.ind = rownames(subset(meta.region,Individual == this.animal))

e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')

library(glmnet)
library(parallel)

# Scale gene expression matrix
if (do.social) {
	social.genes = if (which.social == 1) {
		readRDS(paste0('checkpoints/social_',formatC(this.cutoff),'_predictions_best_genes.rds'))
	} else if (which.social > 1) {
		readRDS(paste0('checkpoints/social_universal_',formatC(this.cutoff),'_predictions_best_genes.rds'))
	}
	e.this = e.keep[social.genes,rownames(subset(meta,Region == this.region | this.region == 'all'))]
} else {
	e.this = e.keep[,rownames(subset(meta,Region == this.region | this.region == 'all'))]
}
e.scaled.region = t(apply(e.this,1,scale))
colnames(e.scaled.region) = colnames(e.this)

e.train = e.scaled.region[,setdiff(rownames(meta.region),this.ind)]
e.test = matrix(e.scaled.region[,this.ind],ncol=length(this.ind),dimnames=list(rownames(e.train),this.ind))

meta.train = subset(meta.region,!Library %in% this.ind)
meta.test = subset(meta.region,Library %in% this.ind)

model = cv.glmnet(
			t(e.train),
			y=meta.train[[predictor]],
			family='gaussian',
			alpha=0.5,
			nfolds=10,
			standardize=FALSE)

out = data.frame(meta[this.ind,c('Library','Individual','Region',predictor)],
		prediction = as.numeric(predict(model,t(e.test),s=model$lambda.min)),
		lambda = model$lambda.min)

betas = coef(model, s=model$lambda.min)[-1,]
betas = betas[betas > 0]

if (this.region == 'all') this.ind = this.animal

saveRDS(betas,file=paste0('checkpoints/',social.prefix,'glmnet_betas_',this.ind,'_',this.region,'.rds'))
saveRDS(out,file=paste0('checkpoints/',social.prefix,'glmnet_predictions_',this.ind,'_',this.region,'.rds'))

}
