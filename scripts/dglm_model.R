#!/usr/bin/env Rscript

source('scripts/_include_options.R')

meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')
e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')
keep.genes = readRDS('checkpoints/keep_genes.rds')

library(parallel)
library(doParallel)

regions = names(keep.genes)

# Initialize output
out = vector('list',length(regions))
names(out) = regions

for (i in 1:length(region.levels)) {
	message('Now analyzing region ',region.levels[i])
	m = droplevels(subset(meta,Region %in% region.levels[i]))
	# Subset expression matrix to only genes passing filters (in keep.genes) and including only libraries of this region
	e.this = e.keep[keep.genes[[region.levels[i]]],colnames(e.keep) %in% m$Library]

#	# Mean-center and scale
#	e.this = t(scale(t(e.this)))

	id.join = rownames(m)
	names(id.join) = m$Individual

	# Drop covariates if they are uniform across dataset
	c.this = model.covariates[apply(m[,model.covariates],2,function(x) length(unique(x))) > 1]

	m.this = m[,c.this]

	clus = makeCluster(n.cores)
	registerDoParallel(cores=n.cores)
	clusterExport(clus,varlist=c('e.this','m.this','c.this','predictor'),envir=environment())

	out[[regions[i]]] = t(parApply(clus,e.this,1,function(y) {
		require(dglm)

		d = m.this
		d$e = y
		results = try(dglm(
			as.formula(paste('e ~',paste(c.this,collapse=' + '))),
			as.formula(paste('~',predictor)),
			family=gaussian(),
			dlink='log',
			data=d
		))
		if ('try-error' %in% class(results)) {
			c(NA,NA,NA)
		} else {
			coef(summary(results)$dispersion.summary)[predictor,c(1,2,4)]
		}
	}))
	stopCluster(clus)

	colnames(out[[regions[i]]]) = c('beta','bvar','pval')

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

saveRDS(regions.array,file='checkpoints/dglm_results.rds')
