#!/usr/bin/env Rscript

source('scripts/_include_options.R')

arguments = commandArgs(trailing=TRUE)

if (length(arguments)) {
	do.social = as.logical(as.integer(arguments[1]))
	which.social = as.integer(arguments[1])
	this.cutoff = if (do.social) as.numeric(arguments[2]) else NULL
	if (length(arguments) > 2 & as.logical(as.integer(arguments[3]))) {
		do.all = TRUE
	} else {
		do.all = FALSE
	}
} else {
	do.social = FALSE
	which.social = NULL
	this.cutoff = NULL
	do.all = FALSE
}

social.prefix = if (do.social) paste0('social_',if (which.social > 1) 'universal_' else '',formatC(this.cutoff),'_') else ''

meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')

if (do.all) {
	prediction.files = with(unique(subset(meta,select='Individual')),
		file.path('checkpoints',paste0(social.prefix,'glmnet_predictions_',Individual,'_','all','.rds'))
	)
	
	feature.files = with(unique(subset(meta,select='Individual')),
		file.path('checkpoints',paste0(social.prefix,'glmnet_betas_',Individual,'_','all','.rds'))
	)
} else {
	prediction.files = with(meta,
		file.path('checkpoints',paste0(social.prefix,'glmnet_predictions_',Library,'_',Region,'.rds'))
	)

	feature.files = with(meta,
		file.path('checkpoints',paste0(social.prefix,'glmnet_betas_',Library,'_',Region,'.rds'))
	)
}

library(parallel)

glmnet.predictions = do.call(rbind,mclapply(prediction.files,readRDS))

glmnet.prediction.results = data.frame(
	meta[rownames(glmnet.predictions),],glmnet.predictions[,c('prediction','lambda')]
)

# Rename some columns to reuse old code
glmnet.prediction.results$known = glmnet.prediction.results[[predictor]]
glmnet.prediction.results$predicted = glmnet.prediction.results$prediction
glmnet.prediction.results$region = glmnet.prediction.results$Region
glmnet.prediction.results$social.rank = glmnet.prediction.results$rank.scaled
glmnet.prediction.results$social.rank.ordinal = glmnet.prediction.results$ordinal.rank
glmnet.prediction.results$delta = with(glmnet.prediction.results,predicted - known)
glmnet.prediction.results$residual = residuals(lm(predicted~known,data=glmnet.prediction.results))

# glmnet.prediction.results$residual = glmnet.prediction.results$delta

animals = sort(unique(meta$Individual))

kf = 0
if (!kf) kf = length(animals)

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

if (do.all) {
	prediction.type = 'merged'
} else {
	prediction.type = 'single'
}

p = ggplot(glmnet.prediction.results,aes(known,predicted)) +
	geom_abline(slope=1,intercept=0,col='black',size=0.5) +
	geom_point() +
	geom_text(data = data.frame(label = lm_eqn(with(glmnet.prediction.results,data.frame(x=known, y=predicted)))), aes(x=(ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5) / 4, y=(ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5) * 0.95, label=label), size=6, parse=TRUE) +
	geom_smooth(method=lm,se=TRUE) +
	coord_fixed() +
	scale_x_continuous(limits=c(with(glmnet.prediction.results,min(c(known,predicted,0))),ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5)) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,min(c(known,predicted,0))),ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5)) +
	theme_classic(base_size=24) +
	xlab(paste0('Known ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'.pdf'),useDingbats=FALSE)

p = ggplot(glmnet.prediction.results,aes(known,predicted)) +
	geom_abline(slope=1,intercept=0,col='black',size=0.5) +
	geom_point(alpha=0) +
	geom_smooth(method=lm,se=TRUE,fill=NA,color=NA,alpha=0) +
	coord_fixed() +
	scale_x_continuous(limits=c(with(glmnet.prediction.results,min(c(known,predicted,0))),ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5)) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,min(c(known,predicted,0))),ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5)) +
	theme_classic(base_size=24) +
	xlab(paste0('Known ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_blank.pdf'),useDingbats=FALSE)

p = ggplot(glmnet.prediction.results,aes(known,predicted,color=region)) +
	geom_abline(slope=1,intercept=0,col='black',size=0.5) +
	geom_point() +
	geom_smooth(method=lm,se=FALSE) +
	coord_fixed() +
	scale_color_manual(name='Region',values=region.colors) +
	scale_x_continuous(limits=c(with(glmnet.prediction.results,min(c(known,predicted,0))),ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5)) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,min(c(known,predicted,0))),ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5)) +
	theme_classic(base_size=16) +
	guides(color = guide_legend(override.aes = list(size = 3,linetype=0))) +
	xlab(paste0('Known ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region.pdf'),useDingbats=FALSE)

p = ggplot(glmnet.prediction.results,aes(known,predicted,color=region)) +
	geom_abline(slope=1,intercept=0,col='black',size=0.1) +
	geom_point() +
	geom_smooth(method=lm,se=TRUE) +
	facet_wrap(~region,nrow=region.rows) +
	coord_fixed() +
	scale_color_manual(name='Region',values=region.colors) +
	scale_x_continuous(limits=c(with(glmnet.prediction.results,min(c(known,predicted,0))),ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5)) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,min(c(known,predicted,0))),ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5)) +
	theme_classic() +
	theme(legend.position='none') +
	xlab(paste0('Known ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(
		within(glmnet.prediction.results,{
			region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(glmnet.prediction.results,glmnet.prediction.results$region),function(x) coef(summary(lm(predicted~known*sex,data=x)))[paste0('known:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
		}),
		aes_string('known','predicted',color='sex')
	) +
	geom_abline(slope=1,intercept=0,col='black',size=0.1) +
	geom_point() +
	geom_smooth(method=lm,se=TRUE) +
	facet_wrap(~region,nrow=region.rows) +
	coord_fixed() +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_continuous(limits=c(with(glmnet.prediction.results,min(c(known,predicted,0))),ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5)) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,min(c(known,predicted,0))),ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5)) +
	theme_classic() +
	theme() +
	guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
	xlab(paste0('Known ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_sex.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(glmnet.prediction.results,aes_string('known','predicted',color='social.rank')) +
	geom_abline(slope=1,intercept=0,col='black',size=0.1) +
	geom_point() +
	geom_smooth(method=lm,se=TRUE) +
	facet_wrap(~region,nrow=region.rows) +
	coord_fixed() +
	scale_x_continuous(limits=c(with(glmnet.prediction.results,min(c(known,predicted,0))),ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5)) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,min(c(known,predicted,0))),ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5)) +
	scale_color_distiller(name='Rank',palette='RdYlGn',breaks=seq(0,1,0.2),labels=c('','low','','','high','')) +
	theme_classic() +
	theme() +
	xlab(paste0('Known ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_rank.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(
		within(glmnet.prediction.results,{
			region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(glmnet.prediction.results,glmnet.prediction.results$region),function(x) coef(summary(lm(residual~known*sex,data=x)))[paste0('known:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
		}),
		aes_string('known','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_point() +
	geom_smooth(method=lm,se=TRUE) +
	facet_wrap(~region,nrow=region.rows) +
	coord_fixed(ratio = 1) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_continuous(limits=c(with(glmnet.prediction.results,min(c(known,predicted,0))),ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5)) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,-max(abs(residual))),with(glmnet.prediction.results,max(abs(residual))))) +
	theme_classic() +
	theme() +
	guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
	xlab(paste0('Known ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_sex_',tolower(predictor.label),'_delta.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(
		within(glmnet.prediction.results,{
			region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(glmnet.prediction.results,glmnet.prediction.results$region),function(x) wilcox.test(residual~sex,data=x)$p.value < 0.05)),paste(region.levels,'*'),region.levels))
		}),
		aes_string('sex','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_boxplot() +
	facet_wrap(~region,nrow=region.rows) +
	coord_fixed(xlim=c(1,2),ratio = 1 * 2 * (with(glmnet.prediction.results,max(abs(residual))) * 2) / (ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5) / (with(glmnet.prediction.results,max(abs(residual))) * 2)) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_discrete(labels=sex.labels) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,-max(abs(residual))),with(glmnet.prediction.results,max(abs(residual))))) +
	theme_classic() +
	theme() +
	xlab('Sex') + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_sex_delta.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(
		within(glmnet.prediction.results,{
			region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(glmnet.prediction.results,glmnet.prediction.results$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
		}),
		aes_string('social.rank','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_point() +
	geom_smooth(method=lm,se=FALSE) +
	facet_wrap(~region,nrow=region.rows) +
	coord_fixed(ratio = (with(glmnet.prediction.results,max(abs(residual))) * 2) / (ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5) / (with(glmnet.prediction.results,max(abs(residual))) * 2)) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,-max(abs(residual))),with(glmnet.prediction.results,max(abs(residual))))) +
	theme_classic() +
	theme(axis.ticks.x = element_line(linetype=c(1,0,1,0,1))) +
	guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
	xlab(paste0('Scaled rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_rank_delta_lm.pdf'),useDingbats=FALSE,height=6,width=8)

# Generate blank version for presentations
p = ggplot(
		within(glmnet.prediction.results,{
			region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(glmnet.prediction.results,glmnet.prediction.results$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
		}),
		aes_string('social.rank','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_point(alpha=0) +
	geom_smooth(method=lm,se=FALSE,fill=NA,color=NA,alpha=0) +
	facet_wrap(~region,nrow=region.rows) +
	coord_fixed(ratio = (with(glmnet.prediction.results,max(abs(residual))) * 2) / (ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5) / (with(glmnet.prediction.results,max(abs(residual))) * 2)) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,-max(abs(residual))),with(glmnet.prediction.results,max(abs(residual))))) +
	theme_classic() +
	theme(axis.ticks.x = element_line(linetype=c(1,0,1,0,1))) +
	guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3,alpha=1))) +
	xlab(paste0('Scaled rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_rank_delta_lm_blank.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(
		within(glmnet.prediction.results,{
			region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(glmnet.prediction.results,glmnet.prediction.results$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
		}),
		aes_string('social.rank','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_point() +
	geom_smooth(method=loess,se=FALSE) +
	facet_wrap(~region,nrow=region.rows) +
	coord_fixed(ratio = (with(glmnet.prediction.results,max(abs(residual))) * 2) / (ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5) / (with(glmnet.prediction.results,max(abs(residual))) * 2)) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,-max(abs(residual))),with(glmnet.prediction.results,max(abs(residual))))) +
	theme_classic() +
	theme(axis.ticks.x = element_line(linetype=c(1,0,1,0,1))) +
	guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
	xlab(paste0('Scaled rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_rank_delta_loess.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(
		within(glmnet.prediction.results,{
			region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(glmnet.prediction.results,glmnet.prediction.results$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
		}),
		aes_string('social.rank','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_point() +
	geom_smooth(method=lm,se=FALSE) +
	coord_fixed(ratio = (with(glmnet.prediction.results,max(abs(residual))) * 2) / (ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5) / (with(glmnet.prediction.results,max(abs(residual))) * 2)) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,-max(abs(residual))),with(glmnet.prediction.results,max(abs(residual))))) +
	theme_classic(base_size=16) +
	theme(axis.ticks.x = element_line(linetype=c(1,0,1,0,1))) +
	guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
	xlab(paste0('Scaled rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_rank_delta_lm.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(
		within(glmnet.prediction.results,{
			region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(glmnet.prediction.results,glmnet.prediction.results$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
		}),
		aes_string('social.rank.ordinal','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_boxplot() +
	coord_fixed(xlim=c(1,3),ratio = 0.9695779 * 3 * (with(glmnet.prediction.results,max(abs(residual))) * 2) / (ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5) / (with(glmnet.prediction.results,max(abs(residual))) * 2)) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_discrete(labels=rank.ordinal.labels) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,-max(abs(residual))),with(glmnet.prediction.results,max(abs(residual))))) +
	theme_classic(base_size=16) +
	theme() +
	guides() +
	xlab(paste0('Ordinal rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_rank_ordinal_delta_boxplot.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(
		within(glmnet.prediction.results,{
			region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(glmnet.prediction.results,glmnet.prediction.results$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
		}),
		aes_string('social.rank.ordinal','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_violin(draw_quantiles=0.5) +
	coord_fixed(xlim=c(1,3),ratio = 0.9695779 * 3 * (with(glmnet.prediction.results,max(abs(residual))) * 2) / (ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5) / (with(glmnet.prediction.results,max(abs(residual))) * 2)) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_discrete(labels=rank.ordinal.labels) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,-max(abs(residual))),with(glmnet.prediction.results,max(abs(residual))))) +
	theme_classic(base_size=16) +
	theme() +
	guides() +
	xlab(paste0('Ordinal rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_rank_ordinal_delta_violin.pdf'),useDingbats=FALSE,height=6,width=8)

dodge = position_dodge(width = 0.9)
p = ggplot(
		within(glmnet.prediction.results,{
			region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(glmnet.prediction.results,glmnet.prediction.results$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
		}),
		aes_string('social.rank.ordinal','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_violin(position=dodge,show.legend=FALSE,fill='#fcfcfc') +
	geom_boxplot(width=0.1,outlier.colour=NA,position=dodge,fill='#eeeeee') +
	coord_fixed(xlim=c(1,3),ratio = 0.9695779 * 3 * (with(glmnet.prediction.results,max(abs(residual))) * 2) / (ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5) / (with(glmnet.prediction.results,max(abs(residual))) * 2)) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_discrete(labels=rank.ordinal.labels) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,-max(abs(residual))),with(glmnet.prediction.results,max(abs(residual))))) +
	theme_classic(base_size=16) +
	theme() +
	guides(color=guide_legend(override.aes=list(fill=NA))) +
	xlab(paste0('Ordinal rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_rank_ordinal_delta_violin_with_boxplot.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(
		within(glmnet.prediction.results,{
			region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(glmnet.prediction.results,glmnet.prediction.results$region),function(x) coef(summary(lm(residual~social.rank*sex,data=x)))[paste0('social.rank:sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
		}),
		aes_string('social.rank','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_point() +
	geom_smooth(method=loess,se=FALSE) +
	coord_fixed(ratio = (with(glmnet.prediction.results,max(abs(residual))) * 2) / (ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5) / (with(glmnet.prediction.results,max(abs(residual))) * 2)) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,-max(abs(residual))),with(glmnet.prediction.results,max(abs(residual))))) +
	theme_classic(base_size=16) +
	theme(axis.ticks.x = element_line(linetype=c(1,0,1,0,1))) +
	guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
	xlab(paste0('Scaled rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_rank_delta_loess.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(glmnet.prediction.results,aes_string('social.rank.ordinal','residual',color='sex')) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_boxplot() +
	facet_wrap(~region,nrow=region.rows) +
	coord_fixed(xlim=c(1,3),ratio = 0.9695779 * 3 * (with(glmnet.prediction.results,max(abs(residual))) * 2) / (ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5) / (with(glmnet.prediction.results,max(abs(residual))) * 2)) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_discrete(labels=rank.ordinal.labels) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,-max(abs(residual))),with(glmnet.prediction.results,max(abs(residual))))) +
	theme_classic() +
	theme() +
	xlab(paste0('Ordinal rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_rank_ordinal_delta_boxplot.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(glmnet.prediction.results,aes_string('social.rank.ordinal','residual',color='sex')) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_violin(draw_quantiles=0.5) +
	facet_wrap(~region,nrow=region.rows) +
	coord_fixed(xlim=c(1,3),ratio = 0.9695779 * 3 * (with(glmnet.prediction.results,max(abs(residual))) * 2) / (ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5) / (with(glmnet.prediction.results,max(abs(residual))) * 2)) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_discrete(labels=rank.ordinal.labels) +
	scale_y_continuous(limits=c(with(glmnet.prediction.results,-max(abs(residual))),with(glmnet.prediction.results,max(abs(residual))))) +
	theme_classic() +
	theme() +
	xlab(paste0('Ordinal rank')) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_rank_ordinal_delta_violin.pdf'),useDingbats=FALSE,height=6,width=8)

for (social.metric in social.metrics) {
	social.label = social.labels[match(social.metric,social.metrics)]

	p = ggplot(
			within(glmnet.prediction.results,{
				region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(glmnet.prediction.results,glmnet.prediction.results$region),function(x) coef(summary(lm(as.formula(paste0('residual~',social.metric,'*sex')),data=x)))[paste0(social.metric,':sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
			}),
			aes_string(social.metric,'residual',color='sex')
		) +
		geom_abline(slope=0,intercept=0,col='black',size=0.2) +
		geom_point() +
		geom_smooth(method=lm,se=FALSE) +
		coord_fixed(ratio = (with(glmnet.prediction.results,max(abs(residual))) * 2) / (ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5) / (with(glmnet.prediction.results,max(abs(residual))) * 2)) +
		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
		scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
		scale_y_continuous(limits=c(with(glmnet.prediction.results,-max(abs(residual))),with(glmnet.prediction.results,max(abs(residual))))) +
		theme_classic(base_size=16) +
		theme(axis.ticks.x = element_line(linetype=c(1,0,1,0,1))) +
		guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
		xlab(paste0(social.label)) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
	ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_',gsub(' ','_',tolower(social.label)),'_delta_lm.pdf'),useDingbats=FALSE,height=6,width=8)

	p = ggplot(
			within(glmnet.prediction.results,{
				region = factor(region,levels=region.levels,labels=ifelse(unlist(lapply(split(glmnet.prediction.results,glmnet.prediction.results$region),function(x) coef(summary(lm(as.formula(paste0('residual~',social.metric,'*sex')),data=x)))[paste0(social.metric,':sex',with(x,levels(sex))[2]),'Pr(>|t|)'] < 0.05)),paste(region.levels,'*'),region.levels))
			}),
			aes_string(social.metric,'residual',color='sex')
		) +
		geom_abline(slope=0,intercept=0,col='black',size=0.2) +
		geom_point() +
		geom_smooth(method=lm,se=FALSE) +
		facet_wrap(~region,nrow=region.rows) +
		coord_fixed(ratio = (with(glmnet.prediction.results,max(abs(residual))) * 2) / (ceiling(with(glmnet.prediction.results,max(c(known,predicted))) / 5) * 5) / (with(glmnet.prediction.results,max(abs(residual))) * 2)) +
		scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
		scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
		scale_y_continuous(limits=c(with(glmnet.prediction.results,-max(abs(residual))),with(glmnet.prediction.results,max(abs(residual))))) +
		theme_classic(base_size=16) +
		theme(axis.ticks.x = element_line(linetype=c(1,0,1,0,1))) +
		guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3))) +
		xlab(paste0(social.label)) + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
	ggsave(p,file=paste0('figures/',social.prefix,'glmnet_predictions_cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'_region_facet_',gsub(' ','_',tolower(social.label)),'_delta_lm.pdf'),useDingbats=FALSE,height=6,width=8)
}

saveRDS(glmnet.prediction.results,file=paste0('checkpoints/glmnet_model_predictions_',social.prefix,'cv_k_',formatC(kf,width=2,flag=0),'_',prediction.type,'.rds'))
