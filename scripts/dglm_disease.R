#!/usr/bin/env Rscript

source('scripts/_include_options.R')

mash.results = readRDS('checkpoints/dglm_mashr_results.rds')
keep.genes = readRDS('checkpoints/keep_genes.rds')

library(mashr)
mash.beta = get_pm(mash.results)
mash.lfsr = get_lfsr(mash.results)
mash.sbet = mash.beta / get_psd(mash.results)

ensembl.gene.names = dimnames(mash.beta)[[1]]
mashr.genes = rownames(mash.beta)
names(mashr.genes) = rownames(mash.beta)

do.data = read.table('data/human_disease_associations.tsv',
	sep='\t',
	quote='',
	col.names=c('protein_id','protein_name','do_id','do_name','z_score','confidence'),
	stringsAsFactors=FALSE)

do.def = unique(subset(do.data,select=c('do_id','do_name')))

do.mmul = readRDS('checkpoints/rnaseq_genes_do.rds')

all.region.fet = all.region.kst = numeric(length=length(ensembl.gene.names))
names(all.region.fet) = names(all.region.kst) = ensembl.gene.names

fraction.shared.cutoff=1/3

all.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
	(sum(mash.lfsr[x,] < fsr.cutoff) >= fraction.shared.cutoff * length(keep.genes)) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff] > 0) >= fraction.shared.cutoff * length(keep.genes)
}))))] = 1

# all.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
# 	(sum(mash.lfsr[x,] < fsr.cutoff) >= fraction.shared.cutoff * length(keep.genes)) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff] < 0) >= fraction.shared.cutoff * length(keep.genes)
# }))))] = -1

all.region.kst[rownames(mash.sbet)] = rowMeans(mash.sbet)

all.region.join = data.frame(ensembl_gene_id = ensembl.gene.names, direction = as.integer(all.region.fet), effect = as.numeric(all.region.kst))

all.region.do = merge(all.region.join, do.mmul, by='ensembl_gene_id')

# Toggle number
all.region.do.pass = subset(all.region.do,do_id %in% names(which(table(subset(all.region.do,confidence >= 0)$do_id) >= 10)))

all.region.do.gene.pass = unique(all.region.do.pass[c('ensembl_gene_id','direction')])

do.inc.total = as.integer(table(factor(all.region.do.gene.pass$direction > 0,levels=c('TRUE','FALSE'))))

all.region.do.split = split(all.region.do.pass,all.region.do.pass$do_id)

library(parallel)

all.region.do.test = do.call(rbind,mclapply(names(all.region.do.split),function(i) {
	x = all.region.do.split[[i]]

	this.inc.total = as.integer(table(factor(x$direction > 0,levels=c('TRUE','FALSE'))))
	contingency.matrix.inc = matrix(rbind(this.inc.total,do.inc.total - this.inc.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

	inc.fet.test = fisher.test(contingency.matrix.inc,alternative='greater')
	inc.kst.test = ks.test(x$effect,subset(all.region.do.pass,do_id != i)$effect,alternative='less')

	data.frame(
		do_id = unique(x$do_id),
		do.size = nrow(x),
		inc.n = this.inc.total[1],
		inc.fet.score = inc.fet.test$estimate,
		inc.kst.score = inc.kst.test$statistic,
		inc.fet.pval = inc.fet.test$p.value,
		inc.kst.pval = inc.kst.test$p.value
	)
},mc.cores=n.cores))

all.region.do.test = within(all.region.do.test,{
	inc.kst.qval = p.adjust(inc.kst.pval,'fdr')
	inc.fet.qval = p.adjust(inc.fet.pval,'fdr')
})

all.region.do.results = merge(all.region.do.test,do.def,by='do_id')
all.region.do.results$dataset = 'DISEASES'
all.region.do.results$set = 'union'
all.region.do.results$region = 'all'

saveRDS(all.region.do.results,file='checkpoints/dglm_disease_results.rds')

