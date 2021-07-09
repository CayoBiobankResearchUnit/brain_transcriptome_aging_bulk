#!/usr/bin/env Rscript

source('scripts/_include_options.R')

mash.main = readRDS('checkpoints/mashr_model_predictions_cv_k_36_merged.rds')
mash.0.0100 = readRDS('checkpoints/mashr_model_predictions_social_0.01_cv_k_36_merged.rds')
mash.0.0010 = readRDS('checkpoints/mashr_model_predictions_social_0.001_cv_k_36_merged.rds')
mash.0.0001 = readRDS('checkpoints/mashr_model_predictions_social_0.0001_cv_k_36_merged.rds')

glmnet.main = readRDS('checkpoints/glmnet_model_predictions_cv_k_36_merged.rds')
glmnet.0.2000 = readRDS('checkpoints/glmnet_model_predictions_social_0.2_cv_k_36_merged.rds')
glmnet.0.1000 = readRDS('checkpoints/glmnet_model_predictions_social_0.1_cv_k_36_merged.rds')
glmnet.0.0500 = readRDS('checkpoints/glmnet_model_predictions_social_0.05_cv_k_36_merged.rds')

mash.main$model = mash.0.0100$model = mash.0.0010$model = mash.0.0001$model = 'mash'

glmnet.main$model = glmnet.0.2000$model = glmnet.0.1000$model = glmnet.0.0500$model = 'glmnet'

mash.main$social.threshold = 1
mash.0.0100$social.threshold = 0.01
mash.0.0010$social.threshold = 0.001
mash.0.0001$social.threshold = 0.0001

glmnet.main$social.threshold = 1
glmnet.0.2000$social.threshold = 0.2
glmnet.0.1000$social.threshold = 0.1
glmnet.0.0500$social.threshold = 0.05

names(glmnet.main)[match(c('Individual','Library','Region'),names(glmnet.main))] = c('id','x','region')
names(glmnet.0.2000)[match(c('Individual','Library','Region'),names(glmnet.0.2000))] = c('id','x','region')
names(glmnet.0.1000)[match(c('Individual','Library','Region'),names(glmnet.0.1000))] = c('id','x','region')
names(glmnet.0.0500)[match(c('Individual','Library','Region'),names(glmnet.0.0500))] = c('id','x','region')

common.columns = Reduce(intersect,lapply(list(
	mash.main,
	mash.0.0100,
	mash.0.0010,
	mash.0.0001,
	glmnet.main,
	glmnet.0.2000,
	glmnet.0.1000,
	glmnet.0.0500
),function(x) colnames(x)))

mash.main     =     mash.main[common.columns]
mash.0.0100   =   mash.0.0100[common.columns]
mash.0.0010   =   mash.0.0010[common.columns]
mash.0.0001   =   mash.0.0001[common.columns]
glmnet.main   =   glmnet.main[common.columns]
glmnet.0.2000 = glmnet.0.2000[common.columns]
glmnet.0.1000 = glmnet.0.1000[common.columns]
glmnet.0.0500 = glmnet.0.0500[common.columns]

all.models = rbind(
	mash.main,
	mash.0.0100,
	mash.0.0010,
	mash.0.0001,
	glmnet.main,
	glmnet.0.2000,
	glmnet.0.1000,
	glmnet.0.0500
)

mash.models = subset(all.models,model == 'mash')
glmnet.models = subset(all.models,model == 'glmnet')

levels(mash.models$region) = gsub('ACCg','ACC',levels(mash.models$region))
levels(glmnet.models$region) = gsub('ACCg','ACC',levels(glmnet.models$region))

library(lme4)
library(lmerTest)
library(tidyverse)

# Report continuous
report.lm = function(threshold=NA,subset.sex=NA,social.var='social.rank',model='mash') {
	require(lme4)

	this = get(paste(model,'models',sep='.'))
	if (!is.na(threshold)) this = subset(this,social.threshold == threshold)
	if (!is.na(subset.sex)) this = subset(this,sex == subset.sex)
	this = droplevels(this)

 	if (!is.na(subset.sex)) {
 		f = residual~social.rank+factor(region)+(1|id)
 	} else {
 		f = residual~social.rank+sex+factor(region)+(1|id)
 	}
 	res = lmer(f,this)
 	out = data.frame(
 		model=model,
 		variable=social.var,
 		set=ifelse(is.na(subset.sex),'all',subset.sex),
 		threshold=ifelse(is.na(threshold),'none',threshold),
 		estimate=coef(summary(res))[social.var,'Estimate'],
 		df=coef(summary(res))[social.var,'df'],
 		t=coef(summary(res))[social.var,'t value'],
 		p=pt(coef(summary(res))[social.var,'t value'],coef(summary(res))[social.var,'df'],lower=TRUE)
 	)
 	out
}

# Report categorical
report.pw = function(threshold=NA,subset.sex=NA,social.var='social.rank.ordinal',model='mash') {
	require(lme4)

	this = get(paste(model,'models',sep='.'))
	if (!is.na(threshold)) this = subset(this,social.threshold == threshold)
	if (!is.na(subset.sex)) this = subset(this,sex == subset.sex)
	this = droplevels(this)
	
	this$social.rank.categorical = factor(as.character(this$social.rank.ordinal),levels=c('H','M','L'))

 	if (!is.na(subset.sex)) {
 		f = residual~social.rank.categorical+factor(region)+(1|id)
 	} else {
 		f = residual~social.rank.categorical+sex+factor(region)+(1|id)
 	}
 	res = lmer(f,mutate(this,social.rank.categorical = factor(as.character(social.rank.ordinal),levels=c('H','M','L'))))
 	res2 = lmer(f,mutate(this,social.rank.categorical = factor(as.character(social.rank.ordinal),levels=c('M','L','H'))))
 	out = data.frame(
 		model=model,
 		variable='social.rank.categorical',
 		set=ifelse(is.na(subset.sex),'all',subset.sex),
 		threshold=ifelse(is.na(threshold),'none',threshold),
 		estimate.hm=-coef(summary(res))['social.rank.categoricalM','Estimate'],
 		df.hm=coef(summary(res))['social.rank.categoricalM','df'],
 		t.hm=-coef(summary(res))['social.rank.categoricalM','t value'],
 		p.hm=pt(coef(summary(res))['social.rank.categoricalM','t value'],coef(summary(res))['social.rank.categoricalM','df'],lower=FALSE),
 		estimate.hl=-coef(summary(res))['social.rank.categoricalL','Estimate'],
 		df.hl=coef(summary(res))['social.rank.categoricalL','df'],
 		t.hl=-coef(summary(res))['social.rank.categoricalL','t value'],
 		p.hl=pt(coef(summary(res))['social.rank.categoricalL','t value'],coef(summary(res))['social.rank.categoricalL','df'],lower=FALSE),
 		estimate.ml=-coef(summary(res2))['social.rank.categoricalL','Estimate'],
 		df.ml=coef(summary(res2))['social.rank.categoricalL','df'],
 		t.ml=-coef(summary(res2))['social.rank.categoricalL','t value'],
 		p.ml=pt(coef(summary(res2))['social.rank.categoricalL','t value'],coef(summary(res))['social.rank.categoricalL','df'],lower=FALSE)
 	)
 	out
}

report.lm.wrapper = function(thresholds,subset.sex=NA,social.var='social.rank',model='mash') {
	out = do.call(rbind,lapply(thresholds,function(i) report.lm(i,subset.sex=subset.sex,social.var=social.var,model=model)))
	out$p.bonf = p.adjust(out$p,'bonferroni')
	subset(out,select=c('model','variable','set','threshold','estimate','df','t','p','p.bonf'))
}
report.pw.wrapper = function(thresholds,subset.sex=NA,social.var='social.rank.ordinal',model='mash') {
	out = do.call(rbind,lapply(thresholds,function(i) report.pw(i,subset.sex=subset.sex,social.var=social.var,model=model)))

 	out$p.hm.bonf = p.adjust(out$p.hm,'bonferroni')
 	out$p.hl.bonf = p.adjust(out$p.hl,'bonferroni')
 	out$p.ml.bonf = p.adjust(out$p.ml,'bonferroni')
 	subset(out,select=c('model','variable','set','threshold','estimate.hm','df.hm','t.hm','p.hm','p.hm.bonf','estimate.hl','df.hl','t.hl','p.hl','p.hl.bonf','estimate.ml','df.ml','t.ml','p.ml','p.ml.bonf'))
}

lm.results = rbind(
	report.lm.wrapper(c(1,0.01,0.001,0.0001),subset.sex=NA ,model='mash'),
	report.lm.wrapper(c(1,0.01,0.001,0.0001),subset.sex='f',model='mash'),
	report.lm.wrapper(c(1,0.01,0.001,0.0001),subset.sex='m',model='mash'),
	report.lm.wrapper(c(1,0.2,0.1,0.05),     subset.sex=NA ,model='glmnet'),
	report.lm.wrapper(c(1,0.2,0.1,0.05),     subset.sex='f',model='glmnet'),
	report.lm.wrapper(c(1,0.2,0.1,0.05),     subset.sex='m',model='glmnet')
)

pw.results = rbind(
	report.pw.wrapper(c(1,0.01,0.001,0.0001),subset.sex=NA ,model='mash'),
	report.pw.wrapper(c(1,0.01,0.001,0.0001),subset.sex='f',model='mash'),
	report.pw.wrapper(c(1,0.01,0.001,0.0001),subset.sex='m',model='mash'),
	report.pw.wrapper(c(1,0.2,0.1,0.05),     subset.sex=NA ,model='glmnet'),
	report.pw.wrapper(c(1,0.2,0.1,0.05),     subset.sex='f',model='glmnet'),
	report.pw.wrapper(c(1,0.2,0.1,0.05),     subset.sex='m',model='glmnet')
)

saveRDS(lm.results,file='checkpoints/comparison_model_predictions_lm_tests.rds')
saveRDS(pw.results,file='checkpoints/comparison_model_predictions_pw_tests.rds')

# unlist(lapply(c(1,0.01,0.001,0.0001),function(i) wilcox.test(subset(mash.models,social.threshold == i & social.rank.ordinal == 'H' & sex == 'f')$residual,subset(mash.models,social.threshold == i & social.rank.ordinal == 'L' & sex == 'f')$residual,alternative='less')$p.value))
# unlist(lapply(c(1,0.01,0.001,0.0001),function(i) wilcox.test(subset(mash.models,social.threshold == i & social.rank.ordinal == 'H' & sex == 'f')$residual,subset(mash.models,social.threshold == i & social.rank.ordinal == 'M' & sex == 'f')$residual,alternative='less')$p.value))
# unlist(lapply(c(1,0.01,0.001,0.0001),function(i) coef(summary(lmer(residual~social.rank+factor(region)+(1|id),data=subset(mash.models,social.threshold == i & sex == 'f'))))['social.rank','Pr(>|t|)']))
# unlist(lapply(c(1,0.01,0.001,0.0001),function(i) pt(coef(summary(lmer(residual~social.rank+factor(region)+(1|id),data=subset(mash.models,social.threshold == i & sex == 'f'))))['social.rank','t value'],coef(summary(lmer(residual~social.rank+factor(region)+(1|id),data=subset(mash.models,social.threshold == i & sex == 'f'))))['social.rank','df'],lower=TRUE)))

unlist(lapply(c(1,0.01,0.001,0.0001),function(i) pt(coef(summary(lmer(residual~social.rank+sex+factor(region)+(1|id),data=subset(mash.models,social.threshold == i))))['social.rank','t value'],coef(summary(lmer(residual~social.rank+sex+factor(region)+(1|id),data=subset(mash.models,social.threshold == i))))['social.rank','df'],lower=TRUE)))

unlist(lapply(c(1,0.01,0.001,0.0001),function(i) pt(coef(summary(lmer(residual~social.rank+factor(region)+(1|id),data=subset(mash.models,social.threshold == i & sex == 'f'))))['social.rank','t value'],coef(summary(lmer(residual~social.rank+factor(region)+(1|id),data=subset(mash.models,social.threshold == i & sex == 'f'))))['social.rank','df'],lower=TRUE)))

unlist(lapply(c(1,0.01,0.001,0.0001),function(i) pt(coef(summary(lmer(residual~social.rank+factor(region)+(1|id),data=subset(mash.models,social.threshold == i & sex == 'm'))))['social.rank','t value'],coef(summary(lmer(residual~social.rank+sex+factor(region)+(1|id),data=subset(mash.models,social.threshold == i))))['social.rank','df'],lower=TRUE)))
# unlist(lapply(c(1,0.2,0.1,0.05),function(i) wilcox.test(subset(glmnet.models,social.threshold == i & social.rank.ordinal == 'H' & sex == 'f')$residual,subset(glmnet.models,social.threshold == i & social.rank.ordinal == 'L' & sex == 'f')$residual,alternative='less')$p.value))
# unlist(lapply(c(1,0.2,0.1,0.05),function(i) wilcox.test(subset(glmnet.models,social.threshold == i & social.rank.ordinal == 'H' & sex == 'f')$residual,subset(glmnet.models,social.threshold == i & social.rank.ordinal == 'M' & sex == 'f')$residual,alternative='less')$p.value))
# unlist(lapply(c(1,0.2,0.1,0.05),function(i) anova(lmer(residual~social.rank+factor(region)+(1|id),data=subset(glmnet.models,social.threshold == i & sex == 'f')))['social.rank','Pr(>F)']))
# unlist(lapply(c(1,0.2,0.1,0.05),function(i) pt(coef(summary(lmer(residual~social.rank+factor(region)+(1|id),data=subset(glmnet.models,social.threshold == i & sex == 'f'))))['social.rank','t value'],coef(summary(lmer(residual~social.rank+factor(region)+(1|id),data=subset(glmnet.models,social.threshold == i & sex == 'f'))))['social.rank','df'],lower=TRUE)))

unlist(lapply(c(1,0.2,0.1,0.05),function(i) pt(coef(summary(lmer(residual~social.rank+sex+factor(region)+(1|id),data=subset(glmnet.models,social.threshold == i))))['social.rank','t value'],coef(summary(lmer(residual~social.rank+sex+factor(region)+(1|id),data=subset(glmnet.models,social.threshold == i))))['social.rank','df'],lower=TRUE)))

## MASH

# P value of high vs mid rank
# unlist(lapply(c(1,0.01,0.001,0.0001),function(i) pt(coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(mash.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('H','M','L'))),social.threshold == i & sex == 'f'))))['social.rank.categoricalM','t value'],coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(mash.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('H','M','L'))),social.threshold == i & sex == 'f'))))['social.rank.categoricalM','df'],lower.tail=FALSE)))
# unlist(lapply(c(1,0.01,0.001,0.0001),function(i) pt(coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(mash.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('H','M','L'))),social.threshold == i & sex == 'm'))))['social.rank.categoricalM','t value'],coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(mash.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('H','M','L'))),social.threshold == i & sex == 'm'))))['social.rank.categoricalM','df'],lower.tail=FALSE)))


# P value of high vs low rank
# unlist(lapply(c(1,0.01,0.001,0.0001),function(i) pt(coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(mash.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('H','M','L'))),social.threshold == i & sex == 'f'))))['social.rank.categoricalL','t value'],coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(mash.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('H','M','L'))),social.threshold == i & sex == 'f'))))['social.rank.categoricalL','df'],lower.tail=FALSE)))
# unlist(lapply(c(1,0.01,0.001,0.0001),function(i) pt(coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(mash.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('H','M','L'))),social.threshold == i & sex == 'm'))))['social.rank.categoricalL','t value'],coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(mash.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('H','M','L'))),social.threshold == i & sex == 'm'))))['social.rank.categoricalL','df'],lower.tail=FALSE)))

# P value of mid vs low rank
# unlist(lapply(c(1,0.01,0.001,0.0001),function(i) pt(coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(mash.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('M','L','H'))),social.threshold == i & sex == 'f'))))['social.rank.categoricalL','t value'],coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(mash.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('M','L','H'))),social.threshold == i & sex == 'f'))))['social.rank.categoricalL','df'],lower.tail=FALSE)))
# unlist(lapply(c(1,0.01,0.001,0.0001),function(i) pt(coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(mash.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('M','L','H'))),social.threshold == i & sex == 'm'))))['social.rank.categoricalL','t value'],coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(mash.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('M','L','H'))),social.threshold == i & sex == 'm'))))['social.rank.categoricalL','df'],lower.tail=FALSE)))


## GLMnet

# P value of high vs mid rank
# unlist(lapply(c(1,0.2,0.1,0.05),function(i) pt(coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(glmnet.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('H','M','L'))),social.threshold == i & sex == 'f'))))['social.rank.categoricalM','t value'],coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(glmnet.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('H','M','L'))),social.threshold == i & sex == 'f'))))['social.rank.categoricalM','df'],lower.tail=FALSE)))
# unlist(lapply(c(1,0.2,0.1,0.05),function(i) pt(coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(glmnet.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('H','M','L'))),social.threshold == i & sex == 'm'))))['social.rank.categoricalM','t value'],coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(glmnet.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('H','M','L'))),social.threshold == i & sex == 'm'))))['social.rank.categoricalM','df'],lower.tail=FALSE)))


# P value of high vs low rank
# unlist(lapply(c(1,0.2,0.1,0.05),function(i) pt(coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(glmnet.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('H','M','L'))),social.threshold == i & sex == 'f'))))['social.rank.categoricalL','t value'],coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(glmnet.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('H','M','L'))),social.threshold == i & sex == 'f'))))['social.rank.categoricalL','df'],lower.tail=FALSE)))
# unlist(lapply(c(1,0.2,0.1,0.05),function(i) pt(coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(glmnet.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('H','M','L'))),social.threshold == i & sex == 'm'))))['social.rank.categoricalL','t value'],coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(glmnet.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('H','M','L'))),social.threshold == i & sex == 'm'))))['social.rank.categoricalL','df'],lower.tail=FALSE)))

# P value of mid vs low rank
# unlist(lapply(c(1,0.2,0.1,0.05),function(i) pt(coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(glmnet.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('M','L','H'))),social.threshold == i & sex == 'f'))))['social.rank.categoricalL','t value'],coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(glmnet.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('M','L','H'))),social.threshold == i & sex == 'f'))))['social.rank.categoricalL','df'],lower.tail=FALSE)))
# unlist(lapply(c(1,0.2,0.1,0.05),function(i) pt(coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(glmnet.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('M','L','H'))),social.threshold == i & sex == 'm'))))['social.rank.categoricalL','t value'],coef(summary(lmer(residual~social.rank.categorical+factor(region)+(1|id),data=subset(mutate(glmnet.models,social.rank.categorical=factor(as.character(social.rank.ordinal),levels=c('M','L','H'))),social.threshold == i & sex == 'm'))))['social.rank.categoricalL','df'],lower.tail=FALSE)))


library(ggplot2)

dodge = position_dodge(width = 0.9)
p = ggplot(
		within(mash.models,{
			social.threshold.factor = factor(formatC(mash.models$social.threshold),levels=c('1','0.01','0.001','0.0001'))
			levels(social.threshold.factor) = c('atop("no","filter")','atop(LFSR[RANK]," < 0.01")','atop(LFSR[RANK]," < 0.001")','atop(LFSR[RANK]," < 0.0001")')
			box.category = paste(sex,social.rank.ordinal,formatC(social.threshold),sep='.')
			box.factor = factor(character(length(box.category)),levels=c('TRUE','FALSE'))
			box.factor[box.category %in% c('f.M.1','f.H.1','f.M.0.01','f.H.0.01','f.M.0.001','f.H.0.001','f.M.0.0001','f.H.0.0001')] = 'TRUE'
			box.factor[is.na(box.factor)] = 'FALSE'
		}),
		aes_string('social.rank.ordinal','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_violin(position=dodge,show.legend=FALSE,fill='#fcfcfc') +
	geom_boxplot(width=0.1,outlier.colour=NA,position=dodge,fill='#eeeeee') +
	facet_wrap(~social.threshold.factor,ncol=1,strip.position='right',labeller=label_parsed) +
	coord_fixed(xlim=c(1,3),ratio = 1/5 * 0.9693487 * 3 * (with(mash.models,max(abs(residual))) * 2) / (ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5) / (with(mash.models,max(abs(residual))) * 2)) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_discrete(labels=rank.ordinal.labels) +
	scale_y_continuous(
		limits=c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual)))),
		breaks=seq(-8,8,2),
		labels=c(-8,'','','',0,'','','',8)
	) +
	theme_classic(base_size=16) +
	theme(
		strip.text.y.right = element_text(angle = 0),
		strip.background=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(fill=NA))) +
	xlab(paste0('Ordinal rank')) + ylab(expression(italic(e)[scriptscriptstyle('AGE')])) # ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_mash_rank_ordinal_delta_violin_with_boxplot.pdf'),useDingbats=FALSE,height=6,width=8)
#ggsave(p,file=paste0('figures/comparison_model_predictions_mash_rank_ordinal_delta_violin_with_boxplot_sig.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(
		within(mash.models,{
			social.threshold.factor = factor(formatC(mash.models$social.threshold),levels=c('1','0.01','0.001','0.0001'))
			levels(social.threshold.factor) = c('atop("no","filter")','atop(LFSR[RANK]," < 0.01")','atop(LFSR[RANK]," < 0.001")','atop(LFSR[RANK]," < 0.0001")')
			box.category = paste(sex,social.rank.ordinal,formatC(social.threshold),sep='.')
			box.factor = factor(character(length(box.category)),levels=c('TRUE','FALSE'))
			box.factor[box.category %in% c('f.M.1','f.H.1','f.M.0.01','f.H.0.01','f.M.0.001','f.H.0.001','f.M.0.0001','f.H.0.0001')] = 'TRUE'
			box.factor[is.na(box.factor)] = 'FALSE'
		}),
		aes_string('social.rank.ordinal','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_violin(position=dodge,show.legend=FALSE,fill='#fcfcfc') +
	geom_boxplot(width=0.1,outlier.colour=NA,position=dodge,fill='#eeeeee') +
	facet_wrap(~social.threshold.factor,ncol=1,strip.position='right',labeller=label_parsed) +
	coord_fixed(xlim=c(1,3),ratio = 1/5 * 0.9693487 * 3 * (with(mash.models,max(abs(residual))) * 2) / (ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5) / (with(mash.models,max(abs(residual))) * 2)) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_discrete(labels=rank.ordinal.labels) +
	scale_y_continuous(
		limits=c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual)))),
		breaks=seq(-10,10,10)
	) +
	theme_classic(base_size=20) +
	theme(
		strip.text.y.right = element_text(angle = 0),
		strip.background=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(fill=NA))) +
	xlab(paste0('Ordinal rank')) + ylab(expression(italic(e)[scriptscriptstyle('AGE')])) # ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_mash_rank_ordinal_delta_violin_with_boxplot_presentation.pdf'),useDingbats=FALSE,height=5,width=8)



p = ggplot(
		subset(within(mash.models,{
			social.threshold.factor = factor(formatC(mash.models$social.threshold),levels=c('1','0.01','0.001','0.0001'))
			levels(social.threshold.factor) = c('atop("no","filter")','atop(LFSR[RANK]," < 0.01")','atop(LFSR[RANK]," < 0.001")','atop(LFSR[RANK]," < 0.0001")')
			box.category = paste(sex,social.rank.ordinal,formatC(social.threshold),sep='.')
			box.factor = factor(character(length(box.category)),levels=c('TRUE','FALSE'))
			box.factor[box.category %in% c('f.M.1','f.H.1','f.M.0.01','f.H.0.01','f.M.0.001','f.H.0.001','f.M.0.0001','f.H.0.0001')] = 'TRUE'
			box.factor[is.na(box.factor)] = 'FALSE'
		}),sex=='f'),
		aes_string('social.rank.ordinal','residual')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_violin(position=dodge,show.legend=FALSE,fill='#fcfcfc',color=sex.colors[1]) +
	geom_boxplot(width=0.1,outlier.colour=NA,position=dodge,fill='#eeeeee',color=sex.colors[1]) +
	facet_wrap(~social.threshold.factor,ncol=1,strip.position='right',labeller=label_parsed) +
	coord_fixed(xlim=c(1,3),ratio = 1/5 * 0.9693487 * 3 * (with(mash.models,max(abs(residual))) * 2) / (ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5) / (with(mash.models,max(abs(residual))) * 2)) +
#	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_discrete(labels=rank.ordinal.labels) +
	scale_y_continuous(
		limits=c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual)))),
		breaks=seq(-8,8,2),
		labels=c(-8,'','','',0,'','','',8)
	) +
	theme_classic(base_size=16) +
	theme(
		strip.text.y.right = element_text(angle = 0),
		strip.background=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(fill=NA))) +
	xlab(paste0('Ordinal rank')) + ylab(expression(italic(e)[scriptscriptstyle('AGE')])) # ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_mash_rank_ordinal_delta_violin_with_boxplot_f.pdf'),useDingbats=FALSE,height=6,width=8)
#ggsave(p,file=paste0('figures/comparison_model_predictions_mash_rank_ordinal_delta_violin_with_boxplot_sig.pdf'),useDingbats=FALSE,height=6,width=8)


p = ggplot(
		subset(within(mash.models,{
			social.threshold.factor = factor(formatC(mash.models$social.threshold),levels=c('1','0.01','0.001','0.0001'))
			levels(social.threshold.factor) = c('atop("no","filter")','atop(LFSR[RANK]," < 0.01")','atop(LFSR[RANK]," < 0.001")','atop(LFSR[RANK]," < 0.0001")')
			box.category = paste(sex,social.rank.ordinal,formatC(social.threshold),sep='.')
			box.factor = factor(character(length(box.category)),levels=c('TRUE','FALSE'))
			box.factor[box.category %in% c('f.M.1','f.H.1','f.M.0.01','f.H.0.01','f.M.0.001','f.H.0.001','f.M.0.0001','f.H.0.0001')] = 'TRUE'
			box.factor[is.na(box.factor)] = 'FALSE'
		}),sex=='f' & social.threshold == 0.001),
		aes_string('social.rank.ordinal','residual')
	) +
	geom_abline(slope=0,intercept=0,color='black',size=0.2) +
	geom_violin(position=dodge,show.legend=FALSE,fill='#fcfcfc',color=sex.colors[1]) +
	geom_boxplot(width=0.1,outlier.colour=NA,position=dodge,fill='#eeeeee',color=sex.colors[1]) +
#	facet_wrap(~social.threshold.factor,ncol=1,strip.position='right',labeller=label_parsed) +
#	coord_fixed(xlim=c(1,3),ratio = 1/5 * 0.9693487 * 3 * (with(mash.models,max(abs(residual))) * 2) / (ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5) / (with(mash.models,max(abs(residual))) * 2)) +
#	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_discrete(labels=paste(rank.ordinal.labels,'rank')) +
	scale_y_continuous(
		limits=with(subset(mash.models,sex=='f'&social.threshold==0.001),c(-max(abs(residual)),max(abs(residual)))),
		breaks=seq(-10,10,5),
#		labels=c(-8,'','','',0,'','','',8)
	) +
	theme_classic(base_size=24) +
	theme(
#		strip.text.y.right = element_text(angle = 0),
#		strip.background=element_blank(),
		axis.title.x = element_blank()
	) +
	guides(color=guide_legend(override.aes=list(fill=NA))) +
	ylab(expression(italic(e)[scriptscriptstyle('AGE')])) # ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_mash_rank_ordinal_delta_violin_with_boxplot_f_0.001.pdf'),useDingbats=FALSE,height=6,width=8)
#ggsave(p,file=paste0('figures/comparison_model_predictions_mash_rank_ordinal_delta_violin_with_boxplot_sig.pdf'),useDingbats=FALSE,height=6,width=8)

library(ggbeeswarm)

p = ggplot(
		subset(within(mash.models,{
			social.threshold.factor = factor(formatC(mash.models$social.threshold),levels=c('1','0.01','0.001','0.0001'))
			levels(social.threshold.factor) = c('atop("no","filter")','atop(LFSR[RANK]," < 0.01")','atop(LFSR[RANK]," < 0.001")','atop(LFSR[RANK]," < 0.0001")')
			box.category = paste(sex,social.rank.ordinal,formatC(social.threshold),sep='.')
			box.factor = factor(character(length(box.category)),levels=c('TRUE','FALSE'))
			box.factor[box.category %in% c('f.M.1','f.H.1','f.M.0.01','f.H.0.01','f.M.0.001','f.H.0.001','f.M.0.0001','f.H.0.0001')] = 'TRUE'
			box.factor[is.na(box.factor)] = 'FALSE'
		}),sex=='f' & social.threshold == 0.001),
		aes_string('social.rank.ordinal','residual')
	) +
	geom_abline(slope=0,intercept=0,color='black',size=0.2) +
	geom_boxplot(width=0.1,size=0.5,outlier.colour=NA,position=dodge,fill='#cccccc',color='#000000',alpha=0.5) +
	geom_beeswarm(size=2,cex=2,color=sex.colors[1]) +
#	geom_violin(position=dodge,show.legend=FALSE,fill='#fcfcfc',color=sex.colors[1]) +
#	geom_boxplot(width=0.1,outlier.colour=NA,position=dodge,fill='#eeeeee',color=sex.colors[1]) +
#	geom_jitter(alpha=0.5,size=1,col=sex.colors[1]) +
#	facet_wrap(~social.threshold.factor,ncol=1,strip.position='right',labeller=label_parsed) +
#	coord_fixed(xlim=c(1,3),ratio = 1/5 * 0.9693487 * 3 * (with(mash.models,max(abs(residual))) * 2) / (ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5) / (with(mash.models,max(abs(residual))) * 2)) +
#	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_discrete(labels=paste(rank.ordinal.labels,'rank')) +
	scale_y_continuous(
		limits=with(subset(mash.models,sex=='f'&social.threshold==0.001),c(-max(abs(residual)),max(abs(residual)))),
		breaks=seq(-10,10,5),
#		labels=c(-8,'','','',0,'','','',8)
	) +
	theme_classic(base_size=24) +
	theme(
#		strip.text.y.right = element_text(angle = 0),
#		strip.background=element_blank(),
		axis.title.x = element_blank()
	) +
	guides(color=guide_legend(override.aes=list(fill=NA))) +
	ylab(expression(italic(e)[scriptscriptstyle('AGE')])) # ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_mash_rank_ordinal_delta_beeswarm_f_0.001.pdf'),useDingbats=FALSE,height=5,width=7)
#ggsave(p,file=paste0('figures/comparison_model_predictions_mash_rank_ordinal_delta_violin_with_boxplot_sig.pdf'),useDingbats=FALSE,height=6,width=8)


p = ggplot(
		subset(within(mash.models,{
			social.threshold.factor = factor(formatC(mash.models$social.threshold),levels=c('1','0.01','0.001','0.0001'))
			levels(social.threshold.factor) = c('atop("no","filter")','atop(LFSR[RANK]," < 0.01")','atop(LFSR[RANK]," < 0.001")','atop(LFSR[RANK]," < 0.0001")')
			box.category = paste(sex,social.rank.ordinal,formatC(social.threshold),sep='.')
			box.factor = factor(character(length(box.category)),levels=c('TRUE','FALSE'))
			box.factor[box.category %in% c('f.M.1','f.H.1','f.M.0.01','f.H.0.01','f.M.0.001','f.H.0.001','f.M.0.0001','f.H.0.0001')] = 'TRUE'
			box.factor[is.na(box.factor)] = 'FALSE'
		}),sex=='m'),
		aes_string('social.rank.ordinal','residual')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_violin(position=dodge,show.legend=FALSE,fill='#fcfcfc',color=sex.colors[2]) +
	geom_boxplot(width=0.1,outlier.colour=NA,position=dodge,fill='#eeeeee',color=sex.colors[2]) +
	facet_wrap(~social.threshold.factor,ncol=1,strip.position='right',labeller=label_parsed) +
	coord_fixed(xlim=c(1,3),ratio = 1/5 * 0.9693487 * 3 * (with(mash.models,max(abs(residual))) * 2) / (ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5) / (with(mash.models,max(abs(residual))) * 2)) +
#	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_discrete(labels=rank.ordinal.labels) +
	scale_y_continuous(
		limits=c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual)))),
		breaks=seq(-8,8,2),
		labels=c(-8,'','','',0,'','','',8)
	) +
	theme_classic(base_size=16) +
	theme(
		strip.text.y.right = element_text(angle = 0),
		strip.background=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(fill=NA))) +
	xlab(paste0('Ordinal rank')) + ylab(expression(italic(e)[scriptscriptstyle('AGE')])) # ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_mash_rank_ordinal_delta_violin_with_boxplot_m.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(
		within(mash.models,{
			social.threshold.factor = factor(formatC(mash.models$social.threshold),levels=c('1','0.01','0.001','0.0001'))
			levels(social.threshold.factor) = c('atop("no","filter")','atop(LFSR[RANK]," < 0.01")','atop(LFSR[RANK]," < 0.001")','atop(LFSR[RANK]," < 0.0001")')
		}),
		aes_string('social.rank','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_point(size=0.5,alpha=0.25,shape=21) +
	geom_smooth(method=lm,se=FALSE,size=0.5) +
	facet_wrap(~social.threshold.factor,ncol=1,strip.position='right',labeller=label_parsed) +
	coord_fixed(
		ratio = 1/5 * (with(mash.models,max(abs(residual))) * 2) / (ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5) / (with(mash.models,max(abs(residual))) * 2),
		ylim=1 * c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual))))
	) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
	scale_y_continuous(
		limits=c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual)))),
		breaks=seq(-8,8,2),
		labels=c(-8,'','','',0,'','','',8)
	) +
	theme_classic(base_size=16) +
	theme(
		strip.text.y.right = element_text(angle = 0),
		strip.background=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(linetype=NA,fill=NA,size=3,alpha=1))) +
	xlab(paste0('Scaled rank')) + ylab(expression(italic(e)[scriptscriptstyle('AGE')])) # + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_mash_rank_scaled_delta_lm.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(
		within(mash.models,{
			social.threshold.factor = factor(formatC(mash.models$social.threshold),levels=c('1','0.01','0.001','0.0001'))
			levels(social.threshold.factor) = c('atop("no","filter")','atop(LFSR[RANK]," < 0.01")','atop(LFSR[RANK]," < 0.001")','atop(LFSR[RANK]," < 0.0001")')
		}),
		aes_string('social.rank','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
#	geom_point(size=0.5,alpha=0.25,shape=21) +
	geom_smooth(aes_string(fill='sex'),method=lm,se=TRUE,size=0,alpha=0.1) +
	geom_smooth(method=lm,se=FALSE,size=1) +
	facet_wrap(~social.threshold.factor,ncol=1,strip.position='right',labeller=label_parsed) +
	coord_fixed(
		ratio = 1/5 * 4 * (with(mash.models,max(abs(residual))) * 2) / (ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5) / (with(mash.models,max(abs(residual))) * 2),
		ylim=0.25 * c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual))))
	) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_fill_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
	scale_y_continuous(
		limits=c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual)))),
		breaks=0.25 * seq(-8,8,8)
	) +
	theme_classic(base_size=16) +
	theme(
		strip.text.y.right = element_text(angle = 0),
		strip.background=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(linetype=1,fill=NA,size=1,alpha=1))) +
	xlab(paste0('Scaled rank')) + ylab(expression(italic(e)[scriptscriptstyle('AGE')])) # # ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_mash_rank_scaled_delta_lm_clipped.pdf'),useDingbats=FALSE,height=6,width=8)




p = ggplot(
		within(mash.models,{
			social.threshold.factor = factor(formatC(mash.models$social.threshold),levels=c('1','0.01','0.001','0.0001'))
			levels(social.threshold.factor) = c('"no filter"','LFSR[RANK]~"< 0.01"','LFSR[RANK]~"< 0.001"','LFSR[RANK]~"< 0.0001"')
		}),
		aes_string('social.rank','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_point(size=0.5,alpha=0.25,shape=21,show.legend=FALSE) +
	geom_smooth(aes_string(fill='sex'),method=lm,se=TRUE,size=0,alpha=0.25) +
	geom_smooth(method=lm,se=FALSE,size=1) +
	facet_wrap(~social.threshold.factor,nrow=1,strip.position='top',labeller=label_parsed) +
	coord_fixed(
		ratio = 15 * (1/5 * (with(mash.models,max(abs(residual))) * 2) / (ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5) / (with(mash.models,max(abs(residual))) * 2)),
		ylim=1 * c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual))))
	) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_fill_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
	scale_y_continuous(
		limits=c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual)))),
		breaks=seq(-10,10,2),
		labels=c(-10,'','','','',0,'','','','',10)
	) +
	theme_classic(base_size=16) +
	theme(
		strip.text = element_text(size=10),
		strip.background=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(linetype=1,fill=NA,size=1,alpha=1))) +
	xlab(paste0('Scaled rank')) + ylab(expression(italic(e)[scriptscriptstyle('AGE')])) # + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_mash_rank_scaled_delta_lm_side.pdf'),useDingbats=FALSE,height=6,width=8)






# levels(mash.models$region) = gsub('ACCg','ACC',levels(mash.models$region))
# levels(glmnet.models$region) = gsub('ACCg','ACC',levels(glmnet.models$region))

p = ggplot(
		within(mash.models,{
			social.threshold.factor = factor(formatC(mash.models$social.threshold),levels=c('1','0.01','0.001','0.0001'))
			levels(social.threshold.factor) = c('atop("no","filter")','atop(LFSR[RANK]," < 0.01")','atop(LFSR[RANK]," < 0.001")','atop(LFSR[RANK]," < 0.0001")')
		}),
		aes(known,predicted,color=region)
	) +
	geom_abline(slope=1,intercept=0,col='black',size=0.1) +
	geom_point(size=0.25) +
	geom_smooth(method=lm,se=TRUE) +
	facet_grid(rows=vars(social.threshold.factor),cols=vars(region),labeller=label_parsed) +
	coord_fixed() +
	scale_color_manual(name='Region',values=region.colors) +
	scale_x_continuous(limits=c(with(mash.models,min(c(known,predicted,0))),ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5)) +
	scale_y_continuous(limits=c(with(mash.models,min(c(known,predicted,0))),ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5)) +
	theme_classic() +
	theme(
		legend.position='none',
		strip.text = element_text(size=8),
		strip.text.y.right = element_text(angle = 0),,
		strip.background=element_blank(),
		axis.text = element_blank()
	) +
	xlab(paste0('Chronological ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
ggsave(p,file=paste0('figures/comparison_model_predictions_mash_predictions.pdf'),useDingbats=FALSE,height=6,width=8)


p = ggplot(
		within(mash.models,{
			social.threshold.factor = factor(formatC(mash.models$social.threshold),levels=c('1','0.01','0.001','0.0001'))
			levels(social.threshold.factor) = c('"no filter"','LFSR[RANK]~"< 0.01"','LFSR[RANK]~"< 0.001"','LFSR[RANK]~"< 0.0001"')
		}),
		aes_string('social.rank','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_point(size=0.5,alpha=0.25,shape=21,show.legend=FALSE) +
	geom_smooth(aes_string(fill='sex'),method=lm,se=TRUE,size=0,alpha=0.25) +
	geom_smooth(method=lm,se=FALSE,size=1) +
	facet_wrap(~social.threshold.factor,nrow=1,strip.position='top',labeller=label_parsed) +
#	coord_fixed(
#		ratio = 15 * (1/5 * (with(mash.models,max(abs(residual))) * 2) / (ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5) / (with(mash.models,max(abs(residual))) * 2)),
#		ylim=1 * c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual))))
#	) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_fill_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
	scale_y_continuous(
		limits=c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual)))),
		breaks=seq(-10,10,5)
	) +
	theme_classic(base_size=20) +
	theme(
		strip.text = element_text(size=11),
		strip.background=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(linetype=1,fill=NA,size=1,alpha=1))) +
	xlab(paste0('Scaled rank')) + ylab(expression(italic(e)[scriptscriptstyle('AGE')])) # + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_mash_rank_scaled_delta_lm_side_presentation.pdf'),useDingbats=FALSE,height=5,width=8)


dodge = position_dodge(width = 0.9)
p = ggplot(
		within(glmnet.models,{
			social.threshold.factor = factor(formatC(glmnet.models$social.threshold),levels=c('1','0.2','0.1','0.05'))
			levels(social.threshold.factor) = c('atop("no","filter")','atop(LFSR[RANK]," < 0.2")','atop(LFSR[RANK]," < 0.1")','atop(LFSR[RANK]," < 0.05")')
		}),
		aes_string('social.rank.ordinal','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_violin(position=dodge,show.legend=FALSE,fill='#fcfcfc') +
	geom_boxplot(width=0.1,outlier.colour=NA,position=dodge,fill='#eeeeee') +
	facet_wrap(~social.threshold.factor,ncol=1,strip.position='right',labeller=label_parsed) +
	coord_fixed(
		xlim=c(1,3),
		ratio = 1/5 * 0.9693487 * 3 * 2 * (with(mash.models,max(abs(residual))) * 2) / (ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5) / (with(mash.models,max(abs(residual))) * 2)
	) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_discrete(labels=rank.ordinal.labels) +
	scale_y_continuous(
		limits=0.5*c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual)))),
		breaks=seq(-4,4,1),
		labels=c(-4,'','','',0,'','','',4)
	) +
	theme_classic(base_size=16) +
	theme(
		strip.text.y.right = element_text(angle = 0),
		strip.background=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(fill=NA))) +
	xlab(paste0('Ordinal rank')) + ylab(expression(italic(e)[scriptscriptstyle('AGE')])) + # ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_glmnet_rank_ordinal_delta_violin_with_boxplot.pdf'),useDingbats=FALSE,height=6,width=8)
ggsave(p,file=paste0('figures/comparison_model_predictions_glmnet_rank_ordinal_delta_violin_with_boxplot.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(
		within(glmnet.models,{
			social.threshold.factor = factor(formatC(glmnet.models$social.threshold),levels=c('1','0.2','0.1','0.05'))
			levels(social.threshold.factor) = c('atop("no","filter")','atop(LFSR[RANK]," < 0.2")','atop(LFSR[RANK]," < 0.1")','atop(LFSR[RANK]," < 0.05")')
		}),
		aes_string('social.rank.ordinal','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_violin(position=dodge,show.legend=FALSE,fill='#fcfcfc') +
	geom_boxplot(width=0.1,outlier.colour=NA,position=dodge,fill='#eeeeee') +
	facet_wrap(~social.threshold.factor,ncol=1,strip.position='right',labeller=label_parsed) +
	coord_fixed(
		xlim=c(1,3),
		ratio = 1/5 * 0.9693487 * 3 * 2 * (with(mash.models,max(abs(residual))) * 2) / (ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5) / (with(mash.models,max(abs(residual))) * 2)
	) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_discrete(labels=rank.ordinal.labels) +
	scale_y_continuous(
		limits=0.5*c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual)))),
		breaks=seq(-4,4,4)
	) +
	theme_classic(base_size=20) +
	theme(
		strip.text.y.right = element_text(angle = 0),
		strip.background=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(fill=NA))) +
	xlab(paste0('Ordinal rank')) + ylab(expression(italic(e)[scriptscriptstyle('AGE')])) + # ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_glmnet_rank_ordinal_delta_violin_with_boxplot_presentation.pdf'),useDingbats=FALSE,height=5,width=8)


p = ggplot(
		subset(within(glmnet.models,{
			social.threshold.factor = factor(formatC(glmnet.models$social.threshold),levels=c('1','0.2','0.1','0.05'))
			levels(social.threshold.factor) = c('atop("no","filter")','atop(LFSR[RANK]," < 0.2")','atop(LFSR[RANK]," < 0.1")','atop(LFSR[RANK]," < 0.05")')
		}),sex=='f'),
		aes_string('social.rank.ordinal','residual')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_violin(position=dodge,show.legend=FALSE,fill='#fcfcfc',color=sex.colors[1]) +
	geom_boxplot(width=0.1,outlier.colour=NA,position=dodge,fill='#eeeeee',color=sex.colors[1]) +
	facet_wrap(~social.threshold.factor,ncol=1,strip.position='right',labeller=label_parsed) +
	coord_fixed(
		xlim=c(1,3),
		ratio = 1/5 * 0.9693487 * 3 * 2 * (with(mash.models,max(abs(residual))) * 2) / (ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5) / (with(mash.models,max(abs(residual))) * 2)
	) +
#	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_discrete(labels=rank.ordinal.labels) +
	scale_y_continuous(
		limits=0.5*c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual)))),
		breaks=seq(-4,4,1),
		labels=c(-4,'','','',0,'','','',4)
	) +
	theme_classic(base_size=16) +
	theme(
		strip.text.y.right = element_text(angle = 0),
		strip.background=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(fill=NA))) +
	xlab(paste0('Ordinal rank')) + ylab(expression(italic(e)[scriptscriptstyle('AGE')])) + # ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_glmnet_rank_ordinal_delta_violin_with_boxplot_f.pdf'),useDingbats=FALSE,height=6,width=8)
ggsave(p,file=paste0('figures/comparison_model_predictions_glmnet_rank_ordinal_delta_violin_with_boxplot_f.pdf'),useDingbats=FALSE,height=6,width=8)

p = ggplot(
		subset(within(glmnet.models,{
			social.threshold.factor = factor(formatC(glmnet.models$social.threshold),levels=c('1','0.2','0.1','0.05'))
			levels(social.threshold.factor) = c('atop("no","filter")','atop(LFSR[RANK]," < 0.2")','atop(LFSR[RANK]," < 0.1")','atop(LFSR[RANK]," < 0.05")')
		}),sex=='m'),
		aes_string('social.rank.ordinal','residual')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_violin(position=dodge,show.legend=FALSE,fill='#fcfcfc',color=sex.colors[2]) +
	geom_boxplot(width=0.1,outlier.colour=NA,position=dodge,fill='#eeeeee',color=sex.colors[2]) +
	facet_wrap(~social.threshold.factor,ncol=1,strip.position='right',labeller=label_parsed) +
	coord_fixed(
		xlim=c(1,3),
		ratio = 1/5 * 0.9693487 * 3 * 2 * (with(mash.models,max(abs(residual))) * 2) / (ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5) / (with(mash.models,max(abs(residual))) * 2)
	) +
#	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_discrete(labels=rank.ordinal.labels) +
	scale_y_continuous(
		limits=0.5*c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual)))),
		breaks=seq(-4,4,1),
		labels=c(-4,'','','',0,'','','',4)
	) +
	theme_classic(base_size=16) +
	theme(
		strip.text.y.right = element_text(angle = 0),
		strip.background=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(fill=NA))) +
	xlab(paste0('Ordinal rank')) + ylab(expression(italic(e)[scriptscriptstyle('AGE')])) + # ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_glmnet_rank_ordinal_delta_violin_with_boxplot_m.pdf'),useDingbats=FALSE,height=6,width=8)
ggsave(p,file=paste0('figures/comparison_model_predictions_glmnet_rank_ordinal_delta_violin_with_boxplot_m.pdf'),useDingbats=FALSE,height=6,width=8)


p = ggplot(
		within(glmnet.models,{
			social.threshold.factor = factor(formatC(glmnet.models$social.threshold),levels=c('1','0.2','0.1','0.05'))
			levels(social.threshold.factor) = c('atop("no","filter")','atop(LFSR[RANK]," < 0.2")','atop(LFSR[RANK]," < 0.1")','atop(LFSR[RANK]," < 0.05")')
		}),
		aes_string('social.rank','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
#	geom_point(size=0.5,alpha=0.25,shape=21) +
	geom_smooth(aes_string(fill='sex'),method=lm,se=TRUE,size=0,alpha=0.1) +
	geom_smooth(method=lm,se=FALSE,size=1) +
	facet_wrap(~social.threshold.factor,ncol=1,strip.position='right',labeller=label_parsed) +
	coord_fixed(
		ratio = 1/5 * 4 * 2 * (with(mash.models,max(abs(residual))) * 2) / (ceiling(with(mash.models,max(c(known,predicted))) / 5) * 5) / (with(mash.models,max(abs(residual))) * 2),
		ylim=0.25 * 0.5 * c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual))))
	) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_fill_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
	scale_y_continuous(
		limits=0.5*c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual)))),
		breaks=0.25 * seq(-4,4,4)
	) +
	theme_classic(base_size=16) +
	theme(
		strip.text.y.right = element_text(angle = 0),
		strip.background=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(linetype=1,fill=NA,size=1,alpha=1))) +
	xlab(paste0('Scaled rank')) + ylab(expression(italic(e)[scriptscriptstyle('AGE')])) + # ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_glmnet_rank_scaled_delta_lm_clipped.pdf'),useDingbats=FALSE,height=6,width=8)
ggsave(p,file=paste0('figures/comparison_model_predictions_glmnet_rank_scaled_delta_lm_clipped.pdf'),useDingbats=FALSE,height=6,width=8)





p = ggplot(
		within(glmnet.models,{
			social.threshold.factor = factor(formatC(glmnet.models$social.threshold),levels=c('1','0.2','0.1','0.05'))
			levels(social.threshold.factor) = c('"no filter"','LFSR[RANK]~"< 0.2"','LFSR[RANK]~"< 0.1"','LFSR[RANK]~"< 0.05"')
		}),
		aes_string('social.rank','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_point(size=0.5,alpha=0.25,shape=21,show.legend=FALSE) +
	geom_smooth(aes_string(fill='sex'),method=lm,se=TRUE,size=0,alpha=0.25) +
	geom_smooth(method=lm,se=FALSE,size=1) +
	facet_wrap(~social.threshold.factor,nrow=1,strip.position='top',labeller=label_parsed) +
	coord_fixed(
		ratio = 3.2 * 15 * (1/5 * (with(glmnet.models,max(abs(residual))) * 2) / (ceiling(with(glmnet.models,max(c(known,predicted))) / 5) * 5) / (with(glmnet.models,max(abs(residual))) * 2)),
		ylim=1/4 * c(with(mash.models,-max(abs(residual))),with(mash.models,max(abs(residual))))
	) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_fill_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
	scale_y_continuous(
#		limits=1/4 * c(with(glmnet.models,-max(abs(residual))),with(glmnet.models,max(abs(residual)))),
		breaks=seq(-3,3,1),
		labels=c(-3,'','',0,'','',3)
	) +
	theme_classic(base_size=16) +
	theme(
		strip.text = element_text(size=10),
		strip.background=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(linetype=1,fill=NA,size=1,alpha=1))) +
	xlab(paste0('Scaled rank')) + ylab(expression(italic(e)[scriptscriptstyle('AGE')])) + # ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_glmnet_rank_scaled_delta_lm_side.pdf'),useDingbats=FALSE,height=6,width=8)
ggsave(p,file=paste0('figures/comparison_model_predictions_glmnet_rank_scaled_delta_lm_side.pdf'),useDingbats=FALSE,height=6,width=8)



p = ggplot(
		within(glmnet.models,{
			social.threshold.factor = factor(formatC(glmnet.models$social.threshold),levels=c('1','0.2','0.1','0.05'))
			levels(social.threshold.factor) = c('"no filter"','LFSR[RANK]~"< 0.2"','LFSR[RANK]~"< 0.1"','LFSR[RANK]~"< 0.05"')
		}),
		aes_string('social.rank','residual',color='sex')
	) +
	geom_abline(slope=0,intercept=0,col='black',size=0.2) +
	geom_point(size=0.5,alpha=0.25,shape=21,show.legend=FALSE) +
	geom_smooth(aes_string(fill='sex'),method=lm,se=TRUE,size=0,alpha=0.25) +
	geom_smooth(method=lm,se=FALSE,size=1) +
	facet_wrap(~social.threshold.factor,nrow=1,strip.position='top',labeller=label_parsed) +
#	coord_fixed(
#		ratio = 15 * (1/5 * (with(glmnet.models,max(abs(residual))) * 2) / (ceiling(with(glmnet.models,max(c(known,predicted))) / 5) * 5) / (with(glmnet.models,max(abs(residual))) * 2)),
#		ylim=1 * c(with(glmnet.models,-max(abs(residual))),with(glmnet.models,max(abs(residual))))
#	) +
	scale_color_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_fill_manual(name='Sex',values=sex.colors,labels=sex.labels) +
	scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=c('','low','','high','')) +
	scale_y_continuous(
	#	limits=c(with(glmnet.models,-max(abs(residual))),with(glmnet.models,max(abs(residual)))),
		breaks=seq(-4,4,2)
	) +
	theme_classic(base_size=24) +
	theme(
		strip.text = element_text(size=16),
		strip.background=element_blank()
	) +
	guides(color=guide_legend(override.aes=list(linetype=1,fill=NA,size=1,alpha=1))) +
	xlab(paste0('Scaled rank')) + ylab(expression(italic(e)[scriptscriptstyle('AGE')])) # + ylab(eval(substitute(expression(Delta*label),list('label'=paste(predictor.label,'(residuals)')))))
ggsave(p,file=paste0('figures/comparison_model_predictions_glmnet_rank_scaled_delta_lm_side_presentation.pdf'),useDingbats=FALSE,height=5,width=8)
ggsave(p,file=paste0('figures/comparison_model_predictions_glmnet_rank_scaled_delta_lm_side_presentation.pdf'),useDingbats=FALSE,height=5,width=15)





p = ggplot(
		within(glmnet.models,{
			social.threshold.factor = factor(formatC(glmnet.models$social.threshold),levels=c('1','0.2','0.1','0.05'))
			levels(social.threshold.factor) = c('atop("no","filter")','atop(LFSR[RANK]," < 0.2")','atop(LFSR[RANK]," < 0.1")','atop(LFSR[RANK]," < 0.05")')
		}),
		aes(known,predicted,color=region)
	) +
	geom_abline(slope=1,intercept=0,col='black',size=0.1) +
	geom_point(size=0.25) +
	geom_smooth(method=lm,se=TRUE) +
	facet_grid(rows=vars(social.threshold.factor),cols=vars(region),labeller=label_parsed) +
	coord_fixed() +
	scale_color_manual(name='Region',values=region.colors) +
	scale_x_continuous(limits=c(with(glmnet.models,min(c(known,predicted,0))),ceiling(with(glmnet.models,max(c(known,predicted))) / 5) * 5)) +
	scale_y_continuous(limits=c(with(glmnet.models,min(c(known,predicted,0))),ceiling(with(glmnet.models,max(c(known,predicted))) / 5) * 5)) +
	theme_classic() +
	theme(
		legend.position='none',
		strip.text = element_text(size=8),
		strip.text.y.right = element_text(angle = 0),
		strip.background=element_blank(),
		axis.text = element_blank()
	) +
	xlab(paste0('Chronological ',tolower(predictor.label),' (',tolower(predictor.units),')')) + ylab(paste0('Predicted ',tolower(predictor.label),' (',tolower(predictor.units),')'))
ggsave(p,file=paste0('figures/comparison_model_predictions_glmnet_predictions.pdf'),useDingbats=FALSE,height=6,width=8)
ggsave(p,file=paste0('figures/comparison_model_predictions_glmnet_predictions.pdf'),useDingbats=FALSE,height=6,width=8)


# FKBP5 alone

e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')

out = reshape2::melt(e.keep['ENSMMUG00000011520',])

out = data.frame(e=out$value,meta[rownames(out),])

p = ggplot() +
	geom_point(data=out,aes_string('exact_age_years','e',color='ordinal.rank')) +
	geom_smooth(data=out,aes_string('exact_age_years','e'),method=lm,se=FALSE,size=2) +
	geom_smooth(data=out,aes_string('exact_age_years','e',color='ordinal.rank'),method=lm,se=FALSE) +
	facet_wrap(~Region,nrow=3) +
	theme_classic()

m.social = readRDS('checkpoints/mashr_results_social_ordinal.rank.L.rds')