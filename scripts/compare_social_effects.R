#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

# arguments = c('exact_age_years','ordinal.rank.L')

source('scripts/_include_options.R')

aging.var = arguments[1]
social.var = arguments[2]

aging.m = readRDS(paste0('checkpoints/mashr_results.rds'))
social.m = readRDS(paste0('checkpoints/mashr_results_social_',social.var,'.rds'))

library(mashr)
library(tidyr)
library(ggplot2)

ignore.checkpoints = FALSE

if (aging.var == 'exact_age_years') {
	aging.pos = 'old age'
	aging.neg = 'young age'
	aging.label = 'age'
	aging.label.full = 'Aging'
	aging.pm = ashr::get_pm(aging.m)
}

if (social.var == 'ordinal.rank.L') {
	social.pos = 'low rank'
	social.neg = 'high rank'
	social.label = 'rank'
	social.label.full = 'Dominance rank'
	social.pm = ashr::get_pm(social.m) * -1 # Reverse sign so that positive signs correspond to increased expression with LOW rank
} else if (social.var == 'CSI') {
	social.pos = 'low integration'
	social.neg = 'high integration'
	social.label = 'integration'
	social.label.full = 'Social integration'
	social.pm = ashr::get_pm(social.m) * -1 # Reverse sign so that positive signs correspond to increased expression with LOW integration
} else if (social.var == 'CSIgroom') {
	social.pos = 'low integration'
	social.neg = 'high integration'
	social.label = 'integration'
	social.label.full = 'Social integration'
	social.pm = ashr::get_pm(social.m) * -1 # Reverse sign so that positive signs correspond to increased expression with LOW integration
} else if (social.var == 'numPartnersGroom') {
	social.pos = 'low integration'
	social.neg = 'high integration'
	social.label = 'integration'
	social.label.full = 'Social integration'
	social.pm = ashr::get_pm(social.m) * -1 # Reverse sign so that positive signs correspond to increased expression with LOW integration
}

aging.lfsr = ashr::get_lfsr(aging.m)
social.lfsr = ashr::get_lfsr(social.m)

colnames(aging.pm) = paste('beta',aging.var,colnames(aging.pm),sep='.')
colnames(social.pm) = paste('beta',social.var,colnames(social.pm),sep='.')
colnames(aging.lfsr) = paste('lfsr',aging.var,colnames(aging.lfsr),sep='.')
colnames(social.lfsr) = paste('lfsr',social.var,colnames(social.lfsr),sep='.')

main.comparisons = data.frame(ensembl_gene_id=rownames(aging.pm),as.data.frame(cbind(aging.pm,social.pm,aging.lfsr,social.lfsr)))

# Takes the data, summarizes macaque effects, filters to significant genes, and outputs data frame
get.matching.data = function(comparison.data,region,social.var=social.var,aging.var=aging.var,sig.cutoff=0.2,keep.criterion='both') {
	social.effect = paste('beta',social.var,region,sep='.')
	aging.effect = paste('beta',aging.var,region,sep='.')
	social.sig = paste('lfsr',social.var,region,sep='.')
	aging.sig = paste('lfsr',aging.var,region,sep='.')
	social.keep = comparison.data[,social.sig] <= sig.cutoff
	aging.keep = comparison.data[,aging.sig] <= sig.cutoff
	if (keep.criterion == 'both') {
		keep = aging.keep & social.keep
	} else if (keep.criterion == 'social') {
		keep = social.keep
	} else if (keep.criterion == 'aging') {
		keep = aging.keep
	} else if (keep.criterion == 'either') {
		keep = aging.keep | social.keep
	} else {
		stop('Error')
	}
	if (sum(keep)) {
		out.data = comparison.data[keep,]
		out.data$gene = rownames(out.data)
		out.data$social.effect = out.data[[social.effect]]
		out.data$aging.effect = out.data[[aging.effect]]
		out.data$social.sig = out.data[[social.sig]]
		out.data$aging.sig = out.data[[aging.sig]]
		out.data$region = factor(region,levels=region.levels)
		out.data[c('gene','region','aging.effect','social.effect','aging.sig','social.sig')]
	} else {
		data.frame(gene=character(0),region=factor(character(0),levels=region.levels),aging.effect=numeric(0),social.effect=numeric(0),aging.sig=numeric(0),social.sig=numeric(0))
	}
}

# matching.data = get.matching.data(main.comparisons,region='dmPFC',social.var=social.var,aging.var=aging.var,sig.cutoff=0.2,keep.criterion='both')

report.stat = function(matching.data,output='spearman') {
	if (output == 'spearman') {
		out = try(with(matching.data,cor(aging.effect,social.effect,use='complete.obs',method='spearman')),silent=TRUE)
		if (class(out) == 'try-error') out = NA
	} else if (output == 'codirectional') {
		out = sum(with(matching.data,(aging.effect > 0) == (social.effect > 0)),na.rm=TRUE) / sum(with(matching.data,!is.na(aging.effect) & !is.na(social.effect)))
		if (is.nan(out)) out = NA
	} else if (output == 'genes') {
		out = sum(with(matching.data,!is.na(aging.effect) & !is.na(social.effect)))
	}
	out
}

sig.data = do.call(rbind,lapply(region.levels,function(i) {
	get.matching.data(main.comparisons,region=i,social.var=social.var,aging.var=aging.var,sig.cutoff=0.2,keep.criterion='both')
}))

library(ggrastr)

levels(sig.data$region) = gsub('ACCg','ACC',levels(sig.data$region))

p = ggplot(
		within(sig.data,{shared = social.effect * aging.effect > 0}),
		aes(aging.effect,social.effect)
	) +
	geom_point_rast(aes(alpha=shared,color=shared),size=0.1) +
	geom_smooth(method=lm) +
	geom_hline(yintercept=0,linetype=2) +
	geom_vline(xintercept=0,linetype=2) +
	facet_wrap(~region,nrow=region.rows,scales='free') +
	scale_alpha_manual(values=c(0.1,0.1)) +
	scale_color_manual(values=c('#000000','#ff0000')) +
	ylab(social.label.full) + xlab(aging.label.full) +
	theme_classic(base_size=24) +
	theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.position='none')
ggsave(p,file=paste0('figures/social_correlation_',social.label,'_',aging.label,'.pdf'),useDingbats=FALSE)

p = ggplot(
		within(sig.data,{shared = social.effect * aging.effect > 0}),
		aes(aging.effect,social.effect)
	) +
	geom_point_rast(aes(alpha=shared,color=shared),size=0.1) +
	geom_smooth(method=lm) +
	geom_hline(yintercept=0,linetype=2) +
	geom_vline(xintercept=0,linetype=2) +
	facet_wrap(~region,nrow=region.rows,scales='free') +
	scale_alpha_manual(values=c(0.1,0.1)) +
	scale_color_manual(values=c('#000000','#ff0000')) +
	ylab(social.label.full) + xlab(aging.label.full) +
	theme_article(base_size=24) +
	theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.position='none')
ggsave(p,file=paste0('figures/social_correlation_',social.label,'_',aging.label,'_presentation.pdf'),useDingbats=FALSE)

library(parallel)

bootstrap.data = do.call(rbind,lapply(region.levels,function(j) {
	do.call(rbind,mclapply(seq(0,1,0.01),function(i) {
		this.data = get.matching.data(main.comparisons,region=j,social.var=social.var,aging.var=aging.var,sig.cutoff=i,keep.criterion='both')
		out = apply(replicate(1000,{
			this.bootstrap = this.data[sample(1:nrow(this.data),replace=TRUE),]
			c(report.stat(this.bootstrap,output='spearman'),report.stat(this.bootstrap,output='codirectional'))
		}),1,quantile,probs=c(0.025,0.975),na.rm=TRUE)
		colnames(out) = c('rho','shared')
		data.frame(
			region = factor(j,levels=region.levels),
			sig.cutoff = i,
			rho=report.stat(this.data,output='spearman'),
			rho.lb = out[1,'rho'],
			rho.ub = out[2,'rho'],
			shared = report.stat(this.data,output='codirectional'),
			shared.lb = out[1,'shared'],
			shared.ub = out[2,'shared'],
			n = report.stat(this.data,output='genes')		
		)
	},mc.cores=n.cores))
}))


library(egg)

levels(bootstrap.data$region) = gsub('ACCg','ACC',levels(bootstrap.data$region))



p = ggplot(bootstrap.data,aes(sig.cutoff,shared)) +
	geom_point(aes(size=n)) +
	geom_errorbar(aes(ymin=shared.lb,ymax=shared.ub),size=0.1,color='blue') +
	geom_hline(yintercept=0.5,linetype=2) +
	facet_wrap(~region,nrow=region.rows) +
	scale_x_continuous(breaks=seq(0,1,0.5)) +
	scale_y_continuous(limits=c(0,1)) +
	scale_size_continuous(name='Genes',range=c(0,1)) +
	theme_classic() +
	xlab('Threshold') +
	ylab('Fraction concordant')
ggsave(p,file=paste0('figures/social_correlation_threshold_shared_fraction_',social.label,'_',aging.label,'.pdf'),height=4,useDingbats=FALSE)

p = ggplot(bootstrap.data,aes(sig.cutoff,rho)) +
	geom_point(aes(size=n)) +
	geom_errorbar(aes(ymin=rho.lb,ymax=rho.ub),size=0.1,color='blue') +
	geom_hline(yintercept=0,linetype=2) +
	facet_wrap(~region,nrow=region.rows) +
	scale_x_continuous(breaks=seq(0,1,0.5)) +
	scale_y_continuous(limits=c(-1,1)) +
	scale_size_continuous(name='Genes',range=c(0,1)) +
	theme_classic() +
	xlab('Threshold') +
	ylab(expression('Spearman\'s'~rho))
ggsave(p,file=paste0('figures/social_correlation_threshold_spearman_rho_',social.label,'_',aging.label,'.pdf'),height=4,useDingbats=FALSE)

p = ggplot(bootstrap.data,aes(sig.cutoff,shared)) +
	geom_point(aes(size=n)) +
	geom_errorbar(aes(ymin=shared.lb,ymax=shared.ub),size=0.1,color='blue') +
	geom_hline(yintercept=0.5,linetype=2) +
	facet_wrap(~region,nrow=region.rows) +
	scale_x_continuous(breaks=seq(0,1,0.5)) +
	scale_y_continuous(limits=c(0,1)) +
	scale_size_continuous(name='Genes',range=c(0,1)) +
	theme_article() +
	xlab('Threshold') +
	ylab('Fraction concordant')
ggsave(p,file=paste0('figures/social_correlation_threshold_shared_fraction_',social.label,'_',aging.label,'_article.pdf'),height=4,useDingbats=FALSE)

p = ggplot(bootstrap.data,aes(sig.cutoff,rho)) +
	geom_point(aes(size=n)) +
	geom_errorbar(aes(ymin=rho.lb,ymax=rho.ub),size=0.1,color='blue') +
	geom_hline(yintercept=0,linetype=2) +
	facet_wrap(~region,nrow=region.rows) +
	scale_x_continuous(breaks=seq(0,1,0.5)) +
	scale_y_continuous(limits=c(-1,1)) +
	scale_size_continuous(name='Genes',range=c(0,1)) +
	theme_article() +
	xlab('Threshold') +
	ylab(expression('Spearman\'s'~rho))
ggsave(p,file=paste0('figures/social_correlation_threshold_spearman_rho_',social.label,'_',aging.label,'_article.pdf'),height=4,useDingbats=FALSE)



p = ggplot(bootstrap.data,aes(sig.cutoff,rho,color=region)) +
	geom_line() +
	geom_hline(yintercept=0,linetype=2) +
	scale_x_continuous(breaks=seq(0,1,0.5)) +
	scale_y_continuous(limits=c(-1,1)) +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=24) +
	theme(legend.position='none') +
	xlab('Threshold') +
	ylab(expression('Spearman\'s'~rho))
ggsave(p,file=paste0('figures/social_correlation_threshold_spearman_rho_',social.label,'_',aging.label,'_presentation.pdf'),height=4,useDingbats=FALSE)





test.data = do.call(rbind,lapply(region.levels,function(j) {
	do.call(rbind,mclapply(seq(0,1,0.01),function(i) {
		this.data = get.matching.data(main.comparisons,region=j,social.var=social.var,aging.var=aging.var,sig.cutoff=i,keep.criterion='both')

		rho = try(with(this.data,cor.test(aging.effect,social.effect,use='complete.obs',method='spearman')),silent=FALSE)
		rho.p = if (class(rho) == 'try-error') NA else rho$p.value

		shared = try(with(this.data,binom.test(sum((social.effect * aging.effect) > 0),nrow(this.data),p=0.5,alternative='greater')),silent=FALSE)
		shared.p = if (class(shared) == 'try-error') NA else shared$p.value

		data.frame(
			region = factor(j,levels=region.levels),
			sig.cutoff = i,
			rho=report.stat(this.data,output='spearman'),
			rho.p,
			shared = report.stat(this.data,output='codirectional'),
			shared.p,
			n = report.stat(this.data,output='genes')		
		)
	},mc.cores=n.cores))
}))


if (!file.exists(paste0('checkpoints/social_rnaseq_genes_entrez.rds'))) {
	library(biomaRt)

	mmul = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL',
		dataset='mmulatta_gene_ensembl')

	mmul.genes = getBM(
		attributes=c('ensembl_gene_id','entrezgene_id'),
		filters = 'ensembl_gene_id',
		values = rownames(main.comparisons),
		mart = mmul)
	saveRDS(mmul.genes,file=paste0('checkpoints/social_rnaseq_genes_entrez.rds'))
} else {
	message('Checkpoint found!\nLoading Entrez annotations from file.')

	mmul.genes = readRDS(paste0('checkpoints/social_rnaseq_genes_entrez.rds'))
}

if (!file.exists(paste0('checkpoints/kegg_pathways_info.rds'))) {

	# Do KEGG correlations
	library(KEGGREST)

	# kegg.genes = keggList('mcc')

	kegg.genes = keggConv('mcc','ncbi-geneid')

	entrez.gene.ids = unique(mmul.genes$entrezgene_id)
	entrez.gene.ids = entrez.gene.ids[!is.na(entrez.gene.ids)]

	entrez.kegg = intersect(paste0('ncbi-geneid:',entrez.gene.ids),names(kegg.genes))

	kegg.lookup = kegg.genes[entrez.kegg]

	kegg.pathways = keggLink('pathway','mcc')
	kegg.ko = keggLink('ko','mcc')

	pathways.all = data.frame(kegg_gene_id=names(kegg.pathways),kegg_pathway_id=as.character(kegg.pathways),stringsAsFactors=FALSE)

	genes.join = mmul.genes
	genes.join$kegg_gene_id = kegg.lookup[paste0('ncbi-geneid:',mmul.genes$entrezgene_id)]

	pathways.all = merge(pathways.all,genes.join,by='kegg_gene_id')

	pathways.info = data.frame(kegg_pathway_id=sort(unique(pathways.all$kegg_pathway_id)),kegg_pathway_name=NA)
	rownames(pathways.info) = pathways.info$kegg_pathway_id
	while (any(is.na(pathways.info$kegg_pathway_name))) {
		these.pathways = rownames(head(subset(pathways.info,is.na(kegg_pathway_name)),10))
		res = keggGet(these.pathways)
		keys = paste0('path:',unlist(lapply(res,function(x) x$ENTRY)))
		vals = unlist(lapply(res,function(x) x$NAME))
		pathways.info[keys,'kegg_pathway_name'] = vals
		pathways.info[setdiff(these.pathways,keys),'kegg_pathway_name'] = 'UNKNOWN'
	}
	
	saveRDS(pathways.info,file='checkpoints/kegg_pathways_info.rds')
	saveRDS(pathways.all,file='checkpoints/kegg_pathways_all.rds')

} else {
	message('Checkpoint found!\nLoading KEGG annotations from file.')

	pathways.info = readRDS('checkpoints/kegg_pathways_info.rds')
	pathways.all = readRDS('checkpoints/kegg_pathways_all.rds')
}

pathways.pass = subset(pathways.all,kegg_pathway_id %in% names(which(table(kegg_pathway_id) >= 10)))

pathways.split = split(pathways.pass,pathways.pass$kegg_pathway_id)

n.permutations = 10000

# Add the all regions category

main.comparisons[[paste0('beta.',aging.var,'.all')]] = rowMeans(do.call(cbind,lapply(region.levels,function(i) main.comparisons[[paste0('beta.',aging.var,'.',i)]])))
main.comparisons[[paste0('beta.',social.var,'.all')]] = rowMeans(do.call(cbind,lapply(region.levels,function(i) main.comparisons[[paste0('beta.',social.var,'.',i)]])))
main.comparisons[[paste0('lfsr.',aging.var,'.all')]] = rowMeans(do.call(cbind,lapply(region.levels,function(i) main.comparisons[[paste0('lfsr.',aging.var,'.',i)]])))
main.comparisons[[paste0('lfsr.',social.var,'.all')]] = rowMeans(do.call(cbind,lapply(region.levels,function(i) main.comparisons[[paste0('lfsr.',social.var,'.',i)]])))

kegg.test.results = do.call(rbind,lapply(c('all',region.levels),function(i) {

	all.rhesus = get.matching.data(main.comparisons,region=i,social.var=social.var,aging.var=aging.var,sig.cutoff=1,keep.criterion='either')
	successes.age = with(all.rhesus,
		as.integer(table((social.effect * aging.effect) > 0)[c('TRUE','FALSE')])
	)

	kegg.tests = do.call(rbind,mclapply(pathways.split,function(x) {
		x = merge(x,main.comparisons,by='ensembl_gene_id',all.x=FALSE,all.y=FALSE)
		sig.x = get.matching.data(x,region=i,social.var=social.var,aging.var=aging.var,sig.cutoff=1,keep.criterion='either')
		successes.this = with(sig.x,c(sum((social.effect * aging.effect) > 0),nrow(sig.x)-sum((social.effect * aging.effect) > 0)))
		# successes.this = with(sig.x,as.integer(table((social.effect * aging.effect) > 0)[c('TRUE','FALSE')]))
		contingency.this = rbind(successes.this,successes.age-successes.this)
		data.frame(
			region = factor(i,levels=region.levels),
			kegg_pathway_id = unique(x$kegg_pathway_id),
			kegg_pathway_name = pathways.info[unique(x$kegg_pathway_id),'kegg_pathway_name'],
			n = nrow(sig.x),
			cor.p = with(sig.x,cor.test(aging.effect,social.effect,method='spearman',alternative='greater'))$p.value,
			cor.est = with(sig.x,cor.test(aging.effect,social.effect,method='spearman',alternative='greater'))$estimate,
			shared.p = with(sig.x,fisher.test(contingency.this,alternative='greater'))$p.value,
			shared.est = with(sig.x,fisher.test(contingency.this,alternative='greater'))$estimate,
			shared.binom.p = with(sig.x,binom.test(sum((social.effect * aging.effect) > 0),nrow(sig.x),p=0.5,alternative='greater'))$p.value,
			shared.binom.est = with(sig.x,binom.test(sum((social.effect * aging.effect) > 0),nrow(sig.x),p=0.5,alternative='greater'))$estimate
		)
	},mc.cores=n.cores))

	permutation.null = mclapply(unique(sort(kegg.tests$n)),function(n) {
		out = replicate(n.permutations,{
			this = all.rhesus[sample(1:nrow(all.rhesus),n,replace=FALSE),]
			as.numeric(with(this,c(
				cor.test(aging.effect,social.effect,method='spearman',alternative='greater')$estimate,
				sum(social.effect * aging.effect > 0)/n
			)))
		})
		rownames(out) = c('cor','shared')
		out
	},mc.cores=n.cores)
	names(permutation.null) = unique(sort(kegg.tests$n))

	kegg.tests$cor.permute.p = sapply(1:nrow(kegg.tests),function(i) {
		x = kegg.tests[i,]
		n = x$n
		sum(permutation.null[[as.character(n)]]['cor',] > x$cor.est) / n.permutations
	})

	kegg.tests$shared.permute.p = sapply(1:nrow(kegg.tests),function(i) {
		x = kegg.tests[i,]
		n = x$n
		sum(permutation.null[[as.character(n)]]['shared',] > x$shared.binom.est) / n.permutations
	})


	kegg.tests = within(kegg.tests,{
		cor.q = p.adjust(cor.p,'fdr')	
		shared.q = p.adjust(shared.p,'fdr')	
		shared.binom.q = p.adjust(shared.binom.p,'fdr')
		cor.permute.q = p.adjust(cor.permute.p,'fdr')
		shared.permute.q = p.adjust(shared.permute.p,'fdr')
	})
	
	kegg.tests

}))

kegg.test.results$region = factor(kegg.test.results$region,levels=c('all',region.levels))
kegg.test.results$region[is.na(kegg.test.results$region)] = 'all'

kegg.test.results = kegg.test.results[order(kegg.test.results$region,kegg.test.results$cor.permute.p),]

# illustrative.pathways = kegg.tests[order(kegg.tests$cor.p),]$kegg_pathway_id[c(1:4,(nrow(kegg.tests)-3):nrow(kegg.tests))]
# illustrative.data = do.call(rbind,lapply(illustrative.pathways,function(i) {
# 	these.genes = subset(pathways.pass,kegg_pathway_id == i)$ensembl_gene_id
# 	out = get.matching.data(subset(main.comparisons,mmul_id %in% these.genes),social.day=social.day,social.var=social.var,aging.var=aging.var,sig.cutoff=1,summarize.proportion=0,keep.criterion='either')
# 	data.frame(kegg_pathway_id = i,out)
# }))
# illustrative.data$kegg_pathway_id = factor(illustrative.data$kegg_pathway_id,levels=illustrative.pathways)
# 
# p = ggplot(
# 		within(illustrative.data,{shared = social.effect * aging.effect > 0}),
# 		aes(aging.effect,social.effect)
# 	) +
# 	geom_point(aes(alpha=shared)) +
# 	geom_smooth(method=lm) +
# 	geom_hline(yintercept=0,linetype=2) +
# 	geom_vline(xintercept=0,linetype=2) +
# 	facet_wrap(~kegg_pathway_id,nrow=2,scales='free') +
# 	scale_alpha_manual(values=c(0.25,1)) +
# 	ylab(paste0('Increase with ',social.pos,' (ants)')) + xlab(paste0('Increase with ',aging.pos,' (rhesus)')) +
# 	theme_classic(base_size=24) +
# 	theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.position='none',strip.text=element_text(size=12))
# ggsave(p,file=paste0('figures/',file.prefix,'illustrative_correlation_d',social.day,'.',social.var,'_',aging.var,'.pdf'),useDingbats=FALSE,height=5)
# 
# message('KEGG tests with correlation FDR < 20%')
# out.print = subset(kegg.tests,cor.q < 0.2,select=c('kegg_pathway_name','n','cor.est','cor.p','cor.q'))
# print(out.print[order(out.print$cor.p),])
# message('KEGG tests with codirectional FDR < 20%')
# out.print = subset(kegg.tests,shared.q < 0.2,select=c('kegg_pathway_name','n','shared.est','shared.p','shared.q'))
# print(out.print[order(out.print$shared.p),])

kegg.all = subset(kegg.test.results,region == 'all')

illustrative.pathways = kegg.all[order(kegg.all$shared.permute.p + kegg.all$cor.permute.p),]$kegg_pathway_id[c(1:4,(nrow(kegg.all)-3):nrow(kegg.all))]
illustrative.data = do.call(rbind,lapply(illustrative.pathways,function(i) {
	these.genes = subset(pathways.pass,kegg_pathway_id == i)$ensembl_gene_id
	out = get.matching.data(subset(main.comparisons,ensembl_gene_id %in% these.genes),region='all',social.var=social.var,aging.var=aging.var,sig.cutoff=1,keep.criterion='either')
	data.frame(kegg_pathway_id = i,out)
}))
illustrative.data$kegg_pathway_id = factor(illustrative.data$kegg_pathway_id,levels=illustrative.pathways)

p = ggplot(
		within(illustrative.data,{shared = aging.effect * social.effect > 0}),
		aes(aging.effect,social.effect)
	) +
	geom_point(aes(alpha=shared)) +
	geom_smooth(method=lm) +
	geom_hline(yintercept=0,linetype=2) +
	geom_vline(xintercept=0,linetype=2) +
	facet_wrap(~kegg_pathway_id,nrow=2,scales='free') +
	scale_alpha_manual(values=c(0.25,1)) +
	ylab(paste0('Increase with ',social.pos)) + xlab(paste0('Increase with ',aging.pos)) +
	theme_classic(base_size=24) +
	theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.position='none',strip.text=element_text(size=12))
ggsave(p,file=paste0('figures/social_correlation_example_',social.label,'.',aging.label,'.pdf'),useDingbats=FALSE,height=5)



saveRDS(kegg.test.results,file=paste0('checkpoints/kegg_tests_',social.label,'_',aging.label,'.rds'))

save(list=ls(),file=paste0('checkpoints/image_',social.label,'_',aging.label,'.RData'))
