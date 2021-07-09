#!/usr/bin/env Rscript

source('scripts/_include_options.R')

arguments = commandArgs(trailing=TRUE)

predictor = arguments[1]

emma.results = readRDS('checkpoints/emma_results.rds')

message('Running predictor: ',predictor)

Bhat = emma.results[,paste('beta',predictor,sep='.'),]
Shat = sqrt(emma.results[,paste('bvar',predictor,sep='.'),])

# For missing data, set beta to 0 and standard error to 1000
# See https://github.com/stephenslab/mashr/issues/17#issuecomment-330628338
# See also pmid:31320509

Bhat[is.na(Bhat)] = 0
Shat[is.na(Shat)] = 1000

library(mashr)

# Create the mashr data object
mash.data = mash_set_data(Bhat,Shat)

# Compute canonical covariance matrices
U.c = cov_canonical(mash.data)  

# # Apply multivariate adaptive shrinkage method in initial model ("naive" run)
# now = Sys.time()
# m.c = mash(mash.data, U.c)
# print(Sys.time() - now) # print time elapsed cause this sometimes takes awhile
# 
# message('Log likelihood (canonical covariance matrix): ',format(get_loglik(m.c),digits=10))
# 
# # Save pairwise sharing matrix
# mash.shared = get_pairwise_sharing(m.c, factor = 0.5)

m.1by1 = mash_1by1(mash.data)
strong.subset = get_significant_results(m.1by1, thresh = 0.05)

# # Get significant results
# strong.subset = get_significant_results(m.c, thresh = 0.05, sig_fn = ashr::get_lfdr)

# # The code below is the preferred method (computes covariance matrix on a null dataset to learn correlation structure among null tests)

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

now = Sys.time()
m.r = mash(mash.random, Ulist = c(U.ed,U.c), outputlevel = 1)
print(Sys.time() - now)

now = Sys.time()
m = mash(mash.data, g=get_fitted_g(m.r), fixg=TRUE)
print(Sys.time() - now)

message('Log likelihood (final model): ',format(get_loglik(m),digits=10))


# # The code below is an alternative method (does not compute covariance matrix on a null dataset)
# # Perform PCA and extreme deconvolution
# U.pca = cov_pca(mash.data,5,subset=strong.subset)
# U.ed = cov_ed(mash.data, U.pca,subset=strong.subset)
# 
# # Run with data-driven covariance matrix
# now = Sys.time()
# m.ed = mash(mash.data, U.ed)
# print(Sys.time() - now)
# 
# message('Log likelihood (data-driven covariance matrix): ',format(get_loglik(m.ed),digits=10))
# 
# # Finally, run with both data-driven and canonical covariance matrices
# now = Sys.time()
# m = mash(mash.data, c(U.c,U.ed))
# print(Sys.time() - now)
# 
# message('Log likelihood (canonical + data-driven covariance matrices): ',format(get_loglik(m),digits=10))
# 
# # Get number of significant genes in each region
# apply(get_lfsr(m),2,function(x) sum(x < fsr.cutoff))

social.lfsr = ashr::get_lfsr(m)
social.beta = ashr::get_pm(m)

coarse.thresholds = round(10^-seq(1/3,4,1/3) * 10^(floor(seq(1/3,4,1/3)-1/3)+1)) / (10^(floor(seq(1/3,4,1/3)-1/3)+1))

fraction.social.regions = fraction.shared.cutoff
for (i in coarse.thresholds) {
	social.genes = names(which(apply(social.lfsr,1,function(x) sum(x < i)) > 0))
	message(length(social.genes),' genes passing LFSR < ',formatC(i))
	saveRDS(
		social.genes,
		file=paste0('checkpoints/social_',formatC(i),'_predictions_best_genes.rds')
	)
	social.genes.universal = rownames(social.beta)[unlist(lapply(rownames(social.beta),function(x) {
		(sum(social.lfsr[x,] < i) >= fraction.social.regions * length(region.levels)) && sum(social.beta[x,][social.lfsr[x,] < i] > 0) >= fraction.social.regions * length(region.levels)
	})) | unlist(lapply(rownames(social.beta),function(x) {
		(sum(social.lfsr[x,] < i) >= fraction.social.regions * length(region.levels)) && sum(social.beta[x,][social.lfsr[x,] < i] < 0) >= fraction.social.regions * length(region.levels)
	}))]
	message(length(social.genes.universal),' genes passing LFSR < ',formatC(i),' in nearly all (',fraction.social.regions*100,'%) regions')
	saveRDS(
		social.genes.universal,
		file=paste0('checkpoints/social_universal_',formatC(i),'_predictions_best_genes.rds')
	)
}

saveRDS(social.lfsr,file='checkpoints/social_genes_lfsr.rds')

saveRDS(m,file=paste0('checkpoints/mashr_results_social_',predictor,'.rds'))
