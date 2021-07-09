#!/usr/bin/env Rscript

source('scripts/_include_options.R')

arguments = commandArgs(trailing=TRUE)

if (length(arguments)) {
	do.social = as.logical(as.integer(arguments[1]))
	which.social = as.integer(arguments[1])
	this.cutoff = if (do.social) as.numeric(arguments[2]) else NULL
} else {
	do.social = FALSE
	which.social = NULL
	this.cutoff = NULL
}

social.prefix = if (do.social) paste0('social_',if (which.social > 1) 'universal_' else '',formatC(this.cutoff),'_') else ''

mash.results = readRDS('checkpoints/mashr_results.rds')

e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')

meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')

library(mashr)

mash.beta = get_pm(mash.results)
mash.lfsr = get_lfsr(mash.results)

coarse.thresholds = round(10^-seq(1/3,4,1/3) * 10^(floor(seq(1/3,4,1/3)-1/3)+1)) / (10^(floor(seq(1/3,4,1/3)-1/3)+1))

if (do.social) {
	social.genes = if (which.social == 1) {
		readRDS(paste0('checkpoints/social_',formatC(this.cutoff),'_predictions_best_genes.rds'))
	} else if (which.social > 1) {
		readRDS(paste0('checkpoints/social_universal_',formatC(this.cutoff),'_predictions_best_genes.rds'))
	}
	mash.beta = mash.beta[social.genes,]
	mash.lfsr = mash.lfsr[social.genes,]
}

get.prediction.genes = function(threshold) {
	rownames(mash.beta)[unlist(lapply(rownames(mash.beta),function(x) {
		(sum(mash.lfsr[x,] < threshold) >= fraction.shared.cutoff * length(region.levels)) && sum(mash.beta[x,][mash.lfsr[x,] < threshold] > 0) >= fraction.shared.cutoff * length(region.levels)
	})) | unlist(lapply(rownames(mash.beta),function(x) {
		(sum(mash.lfsr[x,] < threshold) >= fraction.shared.cutoff * length(region.levels)) && sum(mash.beta[x,][mash.lfsr[x,] < threshold] < 0) >= fraction.shared.cutoff * length(region.levels)
	}))]
}

if (!do.social) saveRDS(get.prediction.genes(fsr.cutoff),file='checkpoints/predictions_default_genes.rds')

run.predictions = function(threshold,merge.regions=TRUE,genes=NULL) {
	# Identify genes showing robust effects across tissues
	prediction.genes = if (is.null(genes)) get.prediction.genes(threshold) else genes

	message('Using ',length(prediction.genes),' most robustly ',tolower(predictor.label),'-predictive genes.')

	# Filter to just those genes
	mash.gene.beta = matrix(mash.beta[prediction.genes,],ncol=ncol(mash.beta),dimnames=list(prediction.genes,colnames(mash.beta)))
	mash.gene.lfsr = matrix(mash.lfsr[prediction.genes,],ncol=ncol(mash.lfsr),dimnames=list(prediction.genes,colnames(mash.lfsr)))

	if (merge.regions) {
		# Construct single multi-region predictor vector
		mash.gene.predictors = rowMeans(mash.gene.beta)
	} else {
		# Construct separate predictors for each region
		mash.gene.predictors = mash.gene.beta
	}

#	# Weight predictors by their standard deviation
#	mash.gene.sd = apply(mash.gene.beta,1,sd)
#	mash.gene.predictors.weighted = mash.gene.predictors / mash.gene.sd

	# Scale gene expression matrix
	e.scaled = t(apply(e.keep,1,scale))
	colnames(e.scaled) = colnames(e.keep)

	# Calculate predictions
	if (merge.regions) {
		# If merging region predictors, multiply 1 x n vector of scaled expression by the generic 1 x n vector of predictors
		mash.gene.predictions = unlist(lapply(colnames(e.scaled),function(x) sum(e.scaled[names(mash.gene.predictors),x] * mash.gene.predictors)))
	} else {
		# If not merging region predictors, multiply 1 x n vector of scaled expression by 1 x n vector for the matching region
		mash.gene.predictions = unlist(lapply(colnames(e.scaled),function(x) sum(e.scaled[rownames(mash.gene.predictors),x] * mash.gene.predictors[,with(meta,as.character(Region)[Library==x])])))
	}
	names(mash.gene.predictions) = colnames(e.scaled)

	# Scale to match real distribution
	mash.standardized.predictions = mean(meta[[predictor]]) + (mash.gene.predictions - mean(mash.gene.predictions)) * sd(meta[[predictor]]) / sd(mash.gene.predictions)

	data.frame(
		known = meta[[predictor]],
		predicted = mash.standardized.predictions[rownames(meta)],
		region = meta$Region,
		sex = meta[[sex.variable]],
		threshold = threshold,
		genes = length(prediction.genes)
	)
}

get.threshold.performance = function(results) {
	require(caret)
	do.call(rbind,lapply(results,function(x) {
		if (all(is.nan(x$predicted))) {
			data.frame(
				threshold=unique(x$threshold),
				genes=unique(x$genes),
				error = NA,
				rmse = NA,
				slope = NA,
				intercept = NA
			)
		} else {
			fit = as.numeric(coef(lm(predicted~known,data=x)))
			data.frame(
				threshold=unique(x$threshold),
				genes=unique(x$genes),
				error = mean(abs(with(x,predicted-known))),
				rmse = with(x,RMSE(predicted,known)),
				slope = fit[2],
				intercept = fit[1]
			)
		}
	}))
}

# Calculate predictions across threshold
mash.prediction.results.merged = lapply(coarse.thresholds,run.predictions)

# Repeat, but with tissue-specific predictors
mash.prediction.results.single = lapply(coarse.thresholds,run.predictions,merge.regions=FALSE)

# Get perfomrannce indicators
mash.threshold.performance.merged = get.threshold.performance(mash.prediction.results.merged)
mash.threshold.performance.single = get.threshold.performance(mash.prediction.results.single)

mash.threshold.performance.merged$predictor = 'region-agnostic'
mash.threshold.performance.single$predictor = 'region-specific'

mash.threshold.performance = rbind(mash.threshold.performance.merged,mash.threshold.performance.single)

best.threshold = with(mash.threshold.performance.merged,threshold[which.min(error)])

# # Temporarily set best threshold to default threshold
# best.threshold = fsr.cutoff

library(ggplot2)

# Run predictions with best threshold

mash.best.genes = get.prediction.genes(best.threshold)
mash.best.predictions.merged = run.predictions(best.threshold)
mash.best.predictions.single = run.predictions(best.threshold,merge.regions=FALSE)

mash.default.genes = get.prediction.genes(fsr.cutoff)
mash.default.predictions.merged = run.predictions(fsr.cutoff)
mash.default.predictions.single = run.predictions(fsr.cutoff,merge.regions=FALSE)

library(egg)

# Source: http://goo.gl/K4yh
lm_eqn = function(df){
    m = lm(y ~ x, df);
    eq = substitute(italic(y) == b * italic(x) + a*","~~italic(r)^2~"="~r2, 
            list(   a = as.character(format(coef(m)[1], digits = 2)), 
                    b = as.character(format(coef(m)[2], digits = 2)), 
                    r2 = as.character(format(summary(m)$r.squared, digits = 3))))
    as.character(as.expression(eq));                 
}

for (prediction.type in c('single','merged')) {

	tmp.df = get(paste0('mash.default.predictions.',prediction.type))
	p = ggplot(tmp.df,aes(known,predicted)) +
		geom_abline(slope=1,intercept=0,col='black',size=0.5) +
		geom_point() +
		geom_text(data = data.frame(label = lm_eqn(with(tmp.df,data.frame(x=known, y=predicted)))), aes(x=(ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / 4, y=(ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) * 0.95, label=label), size=6, parse=TRUE) +
		geom_smooth(method=lm,se=TRUE) +
		coord_fixed() +
		scale_x_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		scale_y_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		theme_classic(base_size=16) +
		xlab(paste0('Known ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_default_threshold_',prediction.type,'.pdf'),useDingbats=FALSE)

	# Blank plot for presentationss
	p = ggplot(tmp.df,aes(known,predicted)) +
		geom_abline(slope=1,intercept=0,col='black',size=0.5) +
#		geom_point() +
#		geom_smooth(method=lm,se=TRUE) +
		coord_fixed() +
		scale_x_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		scale_y_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		theme_classic(base_size=16) +
		xlab(paste0('Known ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_default_threshold_',prediction.type,'_blank.pdf'),useDingbats=FALSE)

	tmp.df = get(paste0('mash.threshold.performance.',prediction.type))
	p = ggplot(tmp.df,aes(-log10(threshold),error)) +
		geom_point(aes(size=genes)) + geom_text(aes(label=genes),nudge_y=0.05) +
		geom_point(aes(-log10(threshold),error),data=subset(tmp.df,threshold == best.threshold),color='red',shape=21,size=4) +
		scale_x_continuous(
			breaks = -log10(tmp.df$threshold),
			labels = gsub(' ','',format(sort(tmp.df$threshold,decreasing=TRUE),scientific=FALSE,drop0trailing=TRUE))
		) +
		scale_size_continuous(name='Num. genes',range=c(1,6),breaks=c(5,10,50,100,500,1000)) +
		xlab('LFSR threshold') + ylab(paste0('Mean error (',tolower(predictor.units),')')) +
		theme_classic(base_size=10)
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_threshold_performance_',prediction.type,'_error.pdf'),height=4,useDingbats=FALSE)

	p = ggplot(tmp.df,aes(-log10(threshold),slope)) +
		geom_point(aes(size=genes)) + geom_text(aes(label=genes),nudge_y=0.1) +
		geom_point(aes(-log10(threshold),slope),data=subset(tmp.df,threshold == best.threshold),color='red',shape=21,size=4) +
		scale_x_continuous(
			breaks = -log10(tmp.df$threshold),
			labels = gsub(' ','',format(sort(tmp.df$threshold,decreasing=TRUE),scientific=FALSE,drop0trailing=TRUE))
		) +
		scale_y_continuous(limits=c(min(c(0,min(tmp.df$slope))),max(c(max(tmp.df$slope),1)))) +
		scale_size_continuous(name='Num. genes',range=c(1,6),breaks=c(5,10,50,100,500,1000)) +
		xlab('LFSR threshold') + ylab(paste0('Slope')) +
		theme_classic(base_size=10)
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_threshold_performance_',prediction.type,'_slope.pdf'),height=4,useDingbats=FALSE)

	p = ggplot(tmp.df,aes(-log10(threshold),intercept)) +
		geom_point(aes(size=genes)) + geom_text(aes(label=genes),nudge_y=-0.5) +
		geom_point(aes(-log10(threshold),intercept),data=subset(tmp.df,threshold == best.threshold),color='red',shape=21,size=4) +
		scale_x_continuous(
			breaks = -log10(tmp.df$threshold),
			labels = gsub(' ','',format(sort(tmp.df$threshold,decreasing=TRUE),scientific=FALSE,drop0trailing=TRUE))
		) +
		scale_y_continuous(limits=c(min(c(0,min(tmp.df$intercept))),max(tmp.df$intercept))) +
		scale_size_continuous(name='Num. genes',range=c(1,6),breaks=c(5,10,50,100,500,1000)) +
		xlab('LFSR threshold') + ylab(paste0('Intercept')) +
		theme_classic(base_size=10)
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_threshold_performance_',prediction.type,'_intercept.pdf'),height=4,useDingbats=FALSE)

	p1 = ggplot(tmp.df,aes(-log10(threshold),error)) +
		geom_point(aes(size=genes)) + geom_text(aes(label=genes),nudge_y=0.2) +
		geom_point(aes(-log10(threshold),error),data=subset(tmp.df,threshold == best.threshold),color='red',shape=21,size=4) +
		scale_x_continuous(
			breaks = -log10(tmp.df$threshold),
			labels = gsub(' ','',format(sort(tmp.df$threshold,decreasing=TRUE),scientific=FALSE,drop0trailing=TRUE))
		) +
		scale_y_continuous(
			limits=c(with(tmp.df,min(error))-0.1,with(tmp.df,max(error))+0.25),
			labels=function(x) sprintf('%0.2f',x)
		) +
		scale_size_continuous(name='Num. genes',range=c(1,6),breaks=c(5,10,50,100,500,1000)) +
		xlab('LFSR threshold') + ylab(paste0('Mean error (',tolower(predictor.units),')')) +
		theme_classic(base_size=10) + theme(legend.position='none',axis.title.x=element_blank())

	p2 = ggplot(tmp.df,aes(-log10(threshold),slope)) +
		geom_point(aes(size=genes)) +
		geom_point(aes(-log10(threshold),slope),data=subset(tmp.df,threshold == best.threshold),color='red',shape=21,size=4) +
		scale_x_continuous(
			breaks = -log10(tmp.df$threshold),
			labels = gsub(' ','',format(sort(tmp.df$threshold,decreasing=TRUE),scientific=FALSE,drop0trailing=TRUE))
		) +
		scale_y_continuous(limits=c(min(c(0,min(tmp.df$slope))),max(c(max(tmp.df$slope),1))),labels=function(x) sprintf('%0.2f',x)) +
		scale_size_continuous(name='Num. genes',range=c(1,6),breaks=c(5,10,50,100,500,1000)) +
		xlab('LFSR threshold') + ylab(paste0('Slope')) +
		theme_classic(base_size=10) + theme(axis.title.x=element_blank())

	p3 = ggplot(tmp.df,aes(-log10(threshold),intercept)) +
		geom_point(aes(size=genes)) +
		geom_point(aes(-log10(threshold),intercept),data=subset(tmp.df,threshold == best.threshold),color='red',shape=21,size=4) +
		scale_x_continuous(
			breaks = -log10(tmp.df$threshold),
			labels = gsub(' ','',format(sort(tmp.df$threshold,decreasing=TRUE),scientific=FALSE,drop0trailing=TRUE))
		) +
		scale_y_continuous(
			limits=c(with(tmp.df,min(intercept))-0.5,with(tmp.df,max(intercept))+0.5),
			labels=function(x) sprintf('%0.2f',x)
		) +
		scale_size_continuous(name='Num. genes',range=c(1,6),breaks=c(5,10,50,100,500,1000)) +
		xlab('LFSR threshold') + ylab(paste0('Intercept')) +
		theme_classic(base_size=10) + theme(legend.position='none')

	pdf(file=paste0('figures/',social.prefix,'model_predictions_threshold_performance_',prediction.type,'_all.pdf'),useDingbats=FALSE,height=6,width=8)
		ggarrange(p1,p2,p3,ncol=1,nrow=3,heights=c(1,1,1),newpage=FALSE)
	dev.off()

	tmp.df = get(paste0('mash.best.predictions.',prediction.type))
	p = ggplot(tmp.df,aes(known,predicted)) +
		geom_abline(slope=1,intercept=0,col='black',size=0.5) +
		geom_point() +
		geom_smooth(method=lm,se=TRUE) +
		coord_fixed() +
		scale_x_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		scale_y_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		theme_classic(base_size=16) +
		xlab(paste0('Known ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_best_threshold_',prediction.type,'.pdf'),useDingbats=FALSE)

	p = ggplot(tmp.df,aes(known,predicted,color=region)) +
		geom_abline(slope=1,intercept=0,col='black',size=0.5) +
		geom_point() +
		geom_smooth(method=lm,se=FALSE) +
		coord_fixed() +
		scale_color_manual(name='Region',values=region.colors) +
		scale_x_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		scale_y_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		theme_classic(base_size=16) +
		guides(color = guide_legend(override.aes = list(size = 3,linetype=0))) +
		xlab(paste0('Known ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_best_threshold_',prediction.type,'_region.pdf'),useDingbats=FALSE)

	p = ggplot(tmp.df,aes(known,predicted,color=region)) +
		geom_abline(slope=1,intercept=0,col='black',size=0.1) +
		geom_point() +
		geom_smooth(method=lm,se=TRUE) +
		facet_wrap(~region,nrow=region.rows) +
		coord_fixed() +
		scale_color_manual(name='Region',values=region.colors) +
		scale_x_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		scale_y_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		theme_classic() +
		theme(legend.position='none') +
		xlab(paste0('Known ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_best_threshold_',prediction.type,'_region_facet.pdf'),useDingbats=FALSE)

	rm(tmp.df)
}

saveRDS(mash.threshold.performance,file=paste0('checkpoints/',social.prefix,'mashr_predictions_threshold_performance.rds'))

saveRDS(mash.best.genes,file=paste0('checkpoints/',social.prefix,'mashr_predictions_best_genes.rds'))
saveRDS(rowMeans(mash.beta[mash.best.genes,]),file=paste0('checkpoints/',social.prefix,'mashr_predictions_best_predictors.rds'))
