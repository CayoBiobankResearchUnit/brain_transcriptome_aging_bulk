#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

if (!length(arguments)) {
	kf = 0
	do.social = FALSE
	which.social = NULL
	message('No argument provided. Assuming leave-one-out cross-validation.')
} else if (length(arguments) > 3) {
	stop('Only 3 arguments supported, ',length(arguments),' provided')
} else if (length(arguments) == 3) {
	kf = as.integer(arguments[1])
	do.social = as.logical(as.integer(arguments[2]))
	which.social = as.integer(arguments[2])
	this.cutoff = as.numeric(arguments[3])
} else {
	kf = as.integer(arguments[1])
	do.social = FALSE
	which.social = NULL
	this.cutoff = NULL
}

source('scripts/_include_options.R')

social.prefix = if (do.social) paste0('social_',if (which.social > 1) 'universal_' else '',formatC(this.cutoff),'_') else ''

meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')
mash.best.genes = readRDS(paste0('checkpoints/',social.prefix,'mashr_predictions_best_genes.rds'))

animals = sort(unique(meta$Individual))

if (!kf) kf = length(animals)

# mash.checkpoints = paste('checkpoints',list.files('checkpoints',pattern=paste0(social.prefix,'predictions_cv_k_',formatC(kf,width=2,flag=0),'_i_[0-9]+.rds')),sep='/')
mash.checkpoints = file.path('checkpoints',paste0(social.prefix,'predictions_cv_k_',formatC(kf,width=2,flag=0),'_i_',formatC(1:kf,width=3,flag=0),'.rds'))

mash.gene.predictions = do.call(rbind,lapply(mash.checkpoints,readRDS))

mash.gene.predictions.merged = subset(mash.gene.predictions,type == 'region-agnostic')
mash.gene.predictions.single = subset(mash.gene.predictions,type == 'region-specific')

# Scale predictions

mash.standardized.predictions.merged = mean(meta[[predictor]]) + (mash.gene.predictions.merged$p - mean(mash.gene.predictions.merged$p)) * sd(meta[[predictor]]) / sd(mash.gene.predictions.merged$p)
mash.standardized.predictions.single = mean(meta[[predictor]]) + (mash.gene.predictions.single$p - mean(mash.gene.predictions.single$p)) * sd(meta[[predictor]]) / sd(mash.gene.predictions.single$p)

mash.prediction.results.merged = data.frame(
	x = mash.gene.predictions.merged$x,
	id = meta[mash.gene.predictions.merged$x,'Individual'],
	known = meta[mash.gene.predictions.merged$x,predictor],
	predicted = mash.standardized.predictions.merged,
	region = meta[mash.gene.predictions.merged$x,'Region'],
	sex = meta[mash.gene.predictions.merged$x,sex.variable],
	social.rank = meta[mash.gene.predictions.merged$x,rank.variable],
	social.rank.ordinal = meta[mash.gene.predictions.merged$x,rank.ordinal.variable],
	meta[mash.gene.predictions.merged$x,social.metrics],
	delta = mash.standardized.predictions.merged - meta[mash.gene.predictions.merged$x,predictor],
	residual = residuals(lm(predicted~known,data.frame(known=meta[mash.gene.predictions.merged$x,predictor],predicted=mash.standardized.predictions.merged))),
	threshold = NA,
	genes = length(mash.best.genes),
	i = mash.gene.predictions.merged$i,
	k = mash.gene.predictions.merged$k,
	type = 'region-agnostic',
	stringsAsFactors=FALSE
)

mash.prediction.results.single = data.frame(
	x = mash.gene.predictions.single$x,
	id = meta[mash.gene.predictions.single$x,'Individual'],
	known = meta[mash.gene.predictions.single$x,predictor],
	predicted = mash.standardized.predictions.single,
	region = meta[mash.gene.predictions.single$x,'Region'],
	sex = meta[mash.gene.predictions.single$x,sex.variable],
	social.rank = meta[mash.gene.predictions.single$x,rank.variable],
	social.rank.ordinal = meta[mash.gene.predictions.single$x,rank.ordinal.variable],
	meta[mash.gene.predictions.single$x,social.metrics],
	delta = mash.standardized.predictions.single - meta[mash.gene.predictions.single$x,predictor],
	residual = residuals(lm(predicted~known,data.frame(known=meta[mash.gene.predictions.single$x,predictor],predicted=mash.standardized.predictions.single))),
	threshold = NA,
	genes = length(mash.best.genes),
	i = mash.gene.predictions.single$i,
	k = mash.gene.predictions.single$k,
	type = 'region-specific',
	stringsAsFactors=FALSE
)

# Regress out the age:deltaAge relationship
#mash.prediction.results.merged = within(mash.prediction.results.merged,{
#	residual.adj = residuals(with(mash.prediction.results.merged,lm(residual~known:sex)))
#})
#
#mash.prediction.results.single = within(mash.prediction.results.single,{
#	residual.adj = residuals(with(mash.prediction.results.single,lm(residual~known:sex)))
#})

# Source: http://goo.gl/K4yh
lm_eqn = function(df){
    m = lm(y ~ x, df);
    eq = substitute(italic(y) == b * italic(x) + a*","~~italic(r)^2~"="~r2, 
            list(   a = as.character(format(coef(m)[1], digits = 2)), 
                    b = as.character(format(coef(m)[2], digits = 2)), 
                    r2 = as.character(format(summary(m)$r.squared, digits = 3))))
    as.character(as.expression(eq));                 
}

library(ggplot2)
library(RColorBrewer)
library(egg)

for (prediction.type in c('merged','single')) {

	tmp.df = get(paste0('mash.prediction.results.',prediction.type))
	levels(tmp.df$region) = gsub('ACCg','ACC',levels(tmp.df$region))
	p = ggplot(tmp.df,aes(known,predicted)) +
		geom_abline(slope=1,intercept=0,col='black',size=0.5) +
		geom_point() +
		geom_text(data = data.frame(label = lm_eqn(with(tmp.df,data.frame(x=known, y=predicted)))), aes(x=(ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / 4, y=(ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5), label=label), size=6, parse=TRUE) +
		geom_smooth(method=lm,se=TRUE) +
		coord_fixed() +
		scale_x_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		scale_y_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		theme_classic(base_size=24) +
		xlab(paste0('Chronological ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'.pdf'),useDingbats=FALSE)

	p = ggplot(tmp.df,aes(known,predicted)) +
		geom_abline(slope=1,intercept=0,col='black',size=0.5) +
		geom_point(alpha=0) +
		geom_smooth(method=lm,se=TRUE,fill=NA,color=NA,alpha=0) +
		coord_fixed() +
		scale_x_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		scale_y_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		theme_classic(base_size=24) +
		xlab(paste0('Chronological ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_blank.pdf'),useDingbats=FALSE)

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
		xlab(paste0('Chronological ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region.pdf'),useDingbats=FALSE)

	
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
		xlab(paste0('Chronological ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet.pdf'),useDingbats=FALSE,height=6,width=8)

	p = ggplot(within(tmp.df,{levels(region) = gsub('ACCg','ACC',levels(region))}),aes(known,predicted,color=region,fill=region)) +
		geom_abline(slope=1,intercept=0,col='black',size=0.1) +
		geom_point() +
		geom_smooth(method=lm,se=TRUE) +
		facet_wrap(~region,nrow=region.rows) +
		coord_fixed() +
		scale_color_manual(name='Region',values=region.colors) +
		scale_fill_manual(name='Region',values=region.colors) +
		scale_x_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		scale_y_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		theme_article(base_size=16) +
		theme(legend.position='none') +
		xlab(paste0('Chronological ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_article.pdf'),useDingbats=FALSE,height=6,width=8)

	p = ggplot(
			within(tmp.df,{
				region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(tmp.df,tmp.df$region),function(x) coef(summary(lm(predicted~known*sex,data=x)))[paste0('known:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
			}),
			aes_string('known','predicted',color='sex')
		) +
		geom_abline(slope=1,intercept=0,col='black',size=0.1) +
		geom_point() +
		geom_smooth(method=lm,se=TRUE) +
		facet_wrap(~region,nrow=region.rows) +
		coord_fixed() +
		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
		scale_x_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		scale_y_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		theme_classic() +
		theme() +
		guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
		xlab(paste0('Chronological ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_sex.pdf'),useDingbats=FALSE,height=6,width=8)

	p = ggplot(tmp.df,aes_string('known','predicted',color='social.rank')) +
		geom_abline(slope=1,intercept=0,col='black',size=0.1) +
		geom_point() +
		geom_smooth(method=lm,se=TRUE) +
		facet_wrap(~region,nrow=region.rows) +
		coord_fixed() +
		scale_x_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		scale_y_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		scale_color_distiller(name='Rank',palette='RdYlGn',breaks=seq(0,1,0.2),labels=c('','low','','','high','')) +
		theme_classic() +
		theme() +
		xlab(paste0('Chronological ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_rank.pdf'),useDingbats=FALSE,height=6,width=8)

	p = ggplot(
			within(tmp.df,{
				region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(tmp.df,tmp.df$region),function(x) coef(summary(lm(residual~known*sex,data=x)))[paste0('known:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
			}),
			aes_string('known','residual',color='sex')
		) +
		geom_abline(slope=0,intercept=0,col='black',size=0.2) +
		geom_point() +
		geom_smooth(method=lm,se=TRUE) +
		facet_wrap(~region,nrow=region.rows) +
		coord_fixed(ratio = 1) +
		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
		scale_x_continuous(limits=c(with(tmp.df,min(c(known,predicted,0))),ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5)) +
		scale_y_continuous(limits=c(with(tmp.df,-max(abs(residual))),with(tmp.df,max(abs(residual))))) +
		theme_classic() +
		theme() +
		guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
		xlab(paste0('Chronological ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_sex_',tolower(predictor.label),'_delta.pdf'),useDingbats=FALSE,height=6,width=8)

	p = ggplot(
			within(tmp.df,{
				region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(tmp.df,tmp.df$region),function(x) wilcox.test(residual~sex,data=x)$p.value < 0.05)),paste(region.levels,'*'),region.levels))
			}),
			aes_string('sex','residual',color='sex')
		) +
		geom_abline(slope=0,intercept=0,col='black',size=0.2) +
		geom_boxplot() +
		facet_wrap(~region,nrow=region.rows) +
		coord_fixed(xlim=c(1,2),ratio = 1 * 2 * (with(tmp.df,max(abs(residual))) * 2) / (ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / (with(tmp.df,max(abs(residual))) * 2)) +
		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
		scale_x_discrete(labels=sex.labels) +
		scale_y_continuous(limits=c(with(tmp.df,-max(abs(residual))),with(tmp.df,max(abs(residual))))) +
		theme_classic() +
		theme() +
		xlab('Sex') + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_sex_delta.pdf'),useDingbats=FALSE,height=6,width=8)

	p = ggplot(
			within(tmp.df,{
				region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(tmp.df,tmp.df$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
			}),
			aes_string('social.rank','residual',color='sex')
		) +
		geom_abline(slope=0,intercept=0,col='black',size=0.2) +
		geom_point() +
		geom_smooth(method=lm,se=FALSE) +
		facet_wrap(~region,nrow=region.rows) +
		coord_fixed(ratio = (with(tmp.df,max(abs(residual))) * 2) / (ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / (with(tmp.df,max(abs(residual))) * 2)) +
		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
		scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
		scale_y_continuous(limits=c(with(tmp.df,-max(abs(residual))),with(tmp.df,max(abs(residual))))) +
		theme_classic() +
		theme(axis.ticks.x = element_line(linetype=c(1,0,1,0,1))) +
		guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
		xlab(paste0('Scaled rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_rank_delta_lm.pdf'),useDingbats=FALSE,height=6,width=8)

	# Generate blank version for presentations
	p = ggplot(
			within(tmp.df,{
				region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(tmp.df,tmp.df$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
			}),
			aes_string('social.rank','residual',color='sex')
		) +
		geom_abline(slope=0,intercept=0,col='black',size=0.2) +
		geom_point(alpha=0) +
		geom_smooth(method=lm,se=FALSE,fill=NA,color=NA,alpha=0) +
		facet_wrap(~region,nrow=region.rows) +
		coord_fixed(ratio = (with(tmp.df,max(abs(residual))) * 2) / (ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / (with(tmp.df,max(abs(residual))) * 2)) +
		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
		scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
		scale_y_continuous(limits=c(with(tmp.df,-max(abs(residual))),with(tmp.df,max(abs(residual))))) +
		theme_classic() +
		theme(axis.ticks.x = element_line(linetype=c(1,0,1,0,1))) +
		guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3,alpha=1))) +
		xlab(paste0('Scaled rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_rank_delta_lm_blank.pdf'),useDingbats=FALSE,height=6,width=8)

	p = ggplot(
			within(tmp.df,{
				region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(tmp.df,tmp.df$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
			}),
			aes_string('social.rank','residual',color='sex')
		) +
		geom_abline(slope=0,intercept=0,col='black',size=0.2) +
		geom_point() +
		geom_smooth(method=loess,se=FALSE) +
		facet_wrap(~region,nrow=region.rows) +
		coord_fixed(ratio = (with(tmp.df,max(abs(residual))) * 2) / (ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / (with(tmp.df,max(abs(residual))) * 2)) +
		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
		scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
		scale_y_continuous(limits=c(with(tmp.df,-max(abs(residual))),with(tmp.df,max(abs(residual))))) +
		theme_classic() +
		theme(axis.ticks.x = element_line(linetype=c(1,0,1,0,1))) +
		guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
		xlab(paste0('Scaled rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_rank_delta_loess.pdf'),useDingbats=FALSE,height=6,width=8)

	p = ggplot(
			within(tmp.df,{
				region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(tmp.df,tmp.df$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
			}),
			aes_string('social.rank','residual',color='sex')
		) +
		geom_abline(slope=0,intercept=0,col='black',size=0.2) +
		geom_point() +
		geom_smooth(method=lm,se=FALSE) +
		coord_fixed(ratio = (with(tmp.df,max(abs(residual))) * 2) / (ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / (with(tmp.df,max(abs(residual))) * 2)) +
		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
		scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
		scale_y_continuous(limits=c(with(tmp.df,-max(abs(residual))),with(tmp.df,max(abs(residual))))) +
		theme_classic(base_size=16) +
		theme(axis.ticks.x = element_line(linetype=c(1,0,1,0,1))) +
		guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
		xlab(paste0('Scaled rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_rank_delta_lm.pdf'),useDingbats=FALSE,height=6,width=8)

	p = ggplot(
			within(tmp.df,{
				region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(tmp.df,tmp.df$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
			}),
			aes_string('social.rank.ordinal','residual',color='sex')
		) +
		geom_abline(slope=0,intercept=0,col='black',size=0.2) +
		geom_boxplot() +
		coord_fixed(xlim=c(1,3),ratio = 0.9695779 * 3 * (with(tmp.df,max(abs(residual))) * 2) / (ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / (with(tmp.df,max(abs(residual))) * 2)) +
		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
		scale_x_discrete(labels=rank.ordinal.labels) +
		scale_y_continuous(limits=c(with(tmp.df,-max(abs(residual))),with(tmp.df,max(abs(residual))))) +
		theme_classic(base_size=16) +
		theme() +
		guides() +
		xlab(paste0('Ordinal rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_rank_ordinal_delta_boxplot.pdf'),useDingbats=FALSE,height=6,width=8)

	p = ggplot(
			within(tmp.df,{
				region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(tmp.df,tmp.df$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
			}),
			aes_string('social.rank.ordinal','residual',color='sex')
		) +
		geom_abline(slope=0,intercept=0,col='black',size=0.2) +
		geom_violin(draw_quantiles=0.5) +
		coord_fixed(xlim=c(1,3),ratio = 0.9695779 * 3 * (with(tmp.df,max(abs(residual))) * 2) / (ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / (with(tmp.df,max(abs(residual))) * 2)) +
		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
		scale_x_discrete(labels=rank.ordinal.labels) +
		scale_y_continuous(limits=c(with(tmp.df,-max(abs(residual))),with(tmp.df,max(abs(residual))))) +
		theme_classic(base_size=16) +
		theme() +
		guides() +
		xlab(paste0('Ordinal rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_rank_ordinal_delta_violin.pdf'),useDingbats=FALSE,height=6,width=8)

	dodge = position_dodge(width = 0.9)
	p = ggplot(
			within(tmp.df,{
				region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(tmp.df,tmp.df$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
			}),
			aes_string('social.rank.ordinal','residual',color='sex')
		) +
		geom_abline(slope=0,intercept=0,col='black',size=0.2) +
		geom_violin(position=dodge,show.legend=FALSE,fill='#fcfcfc') +
		geom_boxplot(width=0.1,outlier.colour=NA,position=dodge,fill='#eeeeee') +
		coord_fixed(xlim=c(1,3),ratio = 0.9695779 * 3 * (with(tmp.df,max(abs(residual))) * 2) / (ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / (with(tmp.df,max(abs(residual))) * 2)) +
		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
		scale_x_discrete(labels=rank.ordinal.labels) +
		scale_y_continuous(limits=c(with(tmp.df,-max(abs(residual))),with(tmp.df,max(abs(residual))))) +
		theme_classic(base_size=16) +
		theme() +
		guides(color=guide_legend(override.aes=list(fill=NA))) +
		xlab(paste0('Ordinal rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_rank_ordinal_delta_violin_with_boxplot.pdf'),useDingbats=FALSE,height=6,width=8)

	p = ggplot(
			within(tmp.df,{
				region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(tmp.df,tmp.df$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
			}),
			aes_string('social.rank','residual',color='sex')
		) +
		geom_abline(slope=0,intercept=0,col='black',size=0.2) +
		geom_point() +
		geom_smooth(method=loess,se=FALSE) +
		coord_fixed(ratio = (with(tmp.df,max(abs(residual))) * 2) / (ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / (with(tmp.df,max(abs(residual))) * 2)) +
		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
		scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
		scale_y_continuous(limits=c(with(tmp.df,-max(abs(residual))),with(tmp.df,max(abs(residual))))) +
		theme_classic(base_size=16) +
		theme(axis.ticks.x = element_line(linetype=c(1,0,1,0,1))) +
		guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
		xlab(paste0('Scaled rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_rank_delta_loess.pdf'),useDingbats=FALSE,height=6,width=8)

	p = ggplot(tmp.df,aes_string('social.rank.ordinal','residual',color='sex')) +
		geom_abline(slope=0,intercept=0,col='black',size=0.2) +
		geom_boxplot() +
		facet_wrap(~region,nrow=region.rows) +
		coord_fixed(xlim=c(1,3),ratio = 0.9695779 * 3 * (with(tmp.df,max(abs(residual))) * 2) / (ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / (with(tmp.df,max(abs(residual))) * 2)) +
		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
		scale_x_discrete(labels=rank.ordinal.labels) +
		scale_y_continuous(limits=c(with(tmp.df,-max(abs(residual))),with(tmp.df,max(abs(residual))))) +
		theme_classic() +
		theme() +
		xlab(paste0('Ordinal rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_rank_ordinal_delta_boxplot.pdf'),useDingbats=FALSE,height=6,width=8)

	p = ggplot(tmp.df,aes_string('social.rank.ordinal','residual',color='sex')) +
		geom_abline(slope=0,intercept=0,col='black',size=0.2) +
		geom_violin(draw_quantiles=0.5) +
		facet_wrap(~region,nrow=region.rows) +
		coord_fixed(xlim=c(1,3),ratio = 0.9695779 * 3 * (with(tmp.df,max(abs(residual))) * 2) / (ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / (with(tmp.df,max(abs(residual))) * 2)) +
		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
		scale_x_discrete(labels=rank.ordinal.labels) +
		scale_y_continuous(limits=c(with(tmp.df,-max(abs(residual))),with(tmp.df,max(abs(residual))))) +
		theme_classic() +
		theme() +
		xlab(paste0('Ordinal rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_rank_ordinal_delta_violin.pdf'),useDingbats=FALSE,height=6,width=8)

	for (social.metric in social.metrics) {
		social.label = social.labels[match(social.metric,social.metrics)]

		p = ggplot(
				within(tmp.df,{
					region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(tmp.df,tmp.df$region),function(x) coef(summary(lm(as.formula(paste0('residual~',social.metric,'*sex')),data=x)))[paste0(social.metric,':sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
				}),
				aes_string(social.metric,'residual',color='sex')
			) +
			geom_abline(slope=0,intercept=0,col='black',size=0.2) +
			geom_point() +
			geom_smooth(method=lm,se=FALSE) +
			coord_fixed(ratio = (with(tmp.df,max(abs(residual))) * 2) / (ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / (with(tmp.df,max(abs(residual))) * 2)) +
			scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
			scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
			scale_y_continuous(limits=c(with(tmp.df,-max(abs(residual))),with(tmp.df,max(abs(residual))))) +
			theme_classic(base_size=16) +
			theme(axis.ticks.x = element_line(linetype=c(1,0,1,0,1))) +
			guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
			xlab(paste0(social.label)) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
		ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_',gsub(' ','_',tolower(social.label)),'_delta_lm.pdf'),useDingbats=FALSE,height=6,width=8)

		p = ggplot(
				within(tmp.df,{
					region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(tmp.df,tmp.df$region),function(x) coef(summary(lm(as.formula(paste0('residual~',social.metric,'*sex')),data=x)))[paste0(social.metric,':sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
				}),
				aes_string(social.metric,'residual',color='sex')
			) +
			geom_abline(slope=0,intercept=0,col='black',size=0.2) +
			geom_point() +
			geom_smooth(method=lm,se=FALSE) +
			facet_wrap(~region,nrow=region.rows) +
			coord_fixed(ratio = (with(tmp.df,max(abs(residual))) * 2) / (ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / (with(tmp.df,max(abs(residual))) * 2)) +
			scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
			scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
			scale_y_continuous(limits=c(with(tmp.df,-max(abs(residual))),with(tmp.df,max(abs(residual))))) +
			theme_classic(base_size=16) +
			theme(axis.ticks.x = element_line(linetype=c(1,0,1,0,1))) +
			guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
			xlab(paste0(social.label)) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
		ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_',gsub(' ','_',tolower(social.label)),'_delta_lm.pdf'),useDingbats=FALSE,height=6,width=8)

	}

#	p = ggplot(
#			within(tmp.df,{
#				region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(tmp.df,tmp.df$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
#			}),
#			aes_string('social.rank','residual.adj',color='sex')
#		) +
#		geom_abline(slope=0,intercept=0,col='black',size=0.2) +
#		geom_point() +
#		geom_smooth(method=lm,se=TRUE) +
#		facet_wrap(~region,nrow=region.rows) +
#		coord_fixed(ratio  = (with(tmp.df,max(abs(residual))) * 2) / (ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / (with(tmp.df,max(abs(residual))) * 2)) +
#		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
#		scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
#		scale_y_continuous(limits=c(with(tmp.df,-max(abs(residual))),with(tmp.df,max(abs(residual))))) +
#		theme_classic() +
#		theme(axis.ticks.x = element_line(linetype=c(1,0,1,0,1))) +
#		guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
#		xlab(paste0('Scaled rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
#	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_rank_delta_adjusted.pdf'),useDingbats=FALSE,height=6,width=8)
#
#	p = ggplot(tmp.df,aes_string('social.rank.ordinal','residual.adj',color='sex')) +
#		geom_abline(slope=0,intercept=0,col='black',size=0.2) +
#		geom_boxplot() +
#		facet_wrap(~region,nrow=region.rows) +
#		coord_fixed(xlim=c(1,3),ratio = 0.9695779 * 3 * (with(tmp.df,max(abs(residual))) * 2) / (ceiling(with(tmp.df,max(c(known,predicted))) / 5) * 5) / (with(tmp.df,max(abs(residual))) * 2)) +
#		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
#		scale_x_discrete(labels=rank.ordinal.labels) +
#		scale_y_continuous(limits=c(with(tmp.df,-max(abs(residual))),with(tmp.df,max(abs(residual))))) +
#		theme_classic() +
#		theme() +
#		xlab(paste0('Ordinal rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
#	ggsave(p,file=paste0('figures/',social.prefix,'model_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_rank_ordinal_delta_adjusted.pdf'),useDingbats=FALSE,height=6,width=8)

	saveRDS(tmp.df,file=paste0('checkpoints/mashr_model_predictions_',social.prefix,'cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'.rds'))
	
	rm(tmp.df)
}
