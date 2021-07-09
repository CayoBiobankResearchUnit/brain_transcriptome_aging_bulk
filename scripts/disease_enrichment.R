#!/usr/bin/env Rscript

source('scripts/_include_options.R')

mash.results = readRDS('checkpoints/mashr_results.rds')
emma.results = readRDS('checkpoints/emma_results.rds')
keep.genes = readRDS('checkpoints/keep_genes.rds')

emma.beta = emma.results[,paste('beta',predictor,sep='.'),]
emma.pval = emma.results[,paste('pval',predictor,sep='.'),]
emma.qval = apply(emma.pval,2,function(x) p.adjust(x,'fdr'))

library(mashr)
mash.beta = get_pm(mash.results)
mash.lfsr = get_lfsr(mash.results)
mash.sbet = mash.beta / get_psd(mash.results)

ensembl.gene.names = dimnames(emma.results)[[1]]
mashr.genes = rownames(mash.beta)
names(mashr.genes) = rownames(mash.beta)

# Import disease associations from DISEASES dataset (i.e., "Disease Ontology")

do.data = read.table('data/human_disease_associations.tsv',
	sep='\t',
	quote='',
	col.names=c('protein_id','protein_name','do_id','do_name','z_score','confidence'),
	stringsAsFactors=FALSE)

do.def = unique(subset(do.data,select=c('do_id','do_name')))

if (ignore.checkpoints || !file.exists('checkpoints/disease_human_orthologs.rds')) {
	library(biomaRt)

	hsap = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
		dataset='hsapiens_gene_ensembl')

	hsap.info = getBM(
		attributes=c('ensembl_gene_id','ensembl_peptide_id','external_gene_name','mmulatta_homolog_ensembl_gene','mmulatta_homolog_orthology_type'),
		mart = hsap)

	saveRDS(hsap.info,file='checkpoints/disease_human_orthologs.rds')
} else {
	message('Checkpoint found!\nLoading human ortholog annotations from file.')

	hsap.info = readRDS('checkpoints/disease_human_orthologs.rds')
}

# For matching ensembl peptides, linking genes is straightforward
do.ensembl = subset(do.data,grepl('^ENSP[0-9]{11}',protein_id))
do.ensembl = merge(do.ensembl,hsap.info,by.x='protein_id','ensembl_peptide_id',all.x=FALSE,all.y=FALSE)
do.ensembl$ensembl_peptide_id = do.ensembl$protein_id

# For everything else
do.proteinname = subset(do.data,!grepl('^ENSP[0-9]{11}',protein_id))
do.proteinname = merge(do.proteinname,hsap.info,by.x='protein_name',by.y='external_gene_name',all.x=FALSE,all.y=FALSE)
do.proteinname$external_gene_name = do.proteinname$protein_name


do.all = rbind(do.ensembl[intersect(names(do.ensembl),names(do.proteinname))],do.proteinname[intersect(names(do.ensembl),names(do.proteinname))])

do.mmul = subset(do.all,
	nchar(mmulatta_homolog_ensembl_gene) > 0 & 
	mmulatta_homolog_orthology_type == 'ortholog_one2one',
	select=c('mmulatta_homolog_ensembl_gene','external_gene_name','ensembl_gene_id','ensembl_peptide_id','do_id','do_name','z_score','confidence'))

names(do.mmul) = c('ensembl_gene_id','external_gene_name','hsapiens_homolog_ensembl_gene','hsapiens_homolog_ensembl_peptide_id','do_id','do_name','z_score','confidence')

do.mmul.out = subset(do.mmul,ensembl_gene_id %in% ensembl.gene.names,select=c('ensembl_gene_id','do_id','do_name','z_score','confidence'))
rownames(do.mmul.out) = NULL
saveRDS(do.mmul.out,file='checkpoints/rnaseq_genes_do.rds')

all.region.fet = all.region.kst = numeric(length=length(ensembl.gene.names))
names(all.region.fet) = names(all.region.kst) = ensembl.gene.names

# Code upregulated genes as 1 and downregulated genes as -1 if they are significant and co-directional in at least a fraction [fraction.shared.cutoff] of regions

all.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
	(sum(mash.lfsr[x,] < fsr.cutoff) >= fraction.shared.cutoff * length(keep.genes)) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff] > 0) >= fraction.shared.cutoff * length(keep.genes)
}))))] = 1
all.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
	(sum(mash.lfsr[x,] < fsr.cutoff) >= fraction.shared.cutoff * length(keep.genes)) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff] < 0) >= fraction.shared.cutoff * length(keep.genes)
}))))] = -1

all.region.kst[rownames(mash.sbet)] = rowMeans(mash.sbet)

all.region.join = data.frame(ensembl_gene_id = ensembl.gene.names, direction = as.integer(all.region.fet), effect = as.numeric(all.region.kst))

all.region.do = merge(all.region.join, do.mmul, by='ensembl_gene_id')

# Toggle number
all.region.do.pass = subset(all.region.do,do_id %in% names(which(table(subset(all.region.do,confidence >= 0)$do_id) >= 10)))

all.region.do.gene.pass = unique(all.region.do.pass[c('ensembl_gene_id','direction')])

do.deg.total = as.integer(table(factor(all.region.do.gene.pass$direction != 0,levels=c('TRUE','FALSE'))))
do.inc.total = as.integer(table(factor(all.region.do.gene.pass$direction > 0,levels=c('TRUE','FALSE'))))
do.dec.total = as.integer(table(factor(all.region.do.gene.pass$direction < 0,levels=c('TRUE','FALSE'))))

all.region.do.split = split(all.region.do.pass,all.region.do.pass$do_id)

library(parallel)

all.region.do.test = do.call(rbind,mclapply(names(all.region.do.split),function(i) {
	x = all.region.do.split[[i]]

	this.deg.total = as.integer(table(factor(x$direction != 0,levels=c('TRUE','FALSE'))))
	contingency.matrix.deg = matrix(rbind(this.deg.total,do.deg.total - this.deg.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

	this.inc.total = as.integer(table(factor(x$direction > 0,levels=c('TRUE','FALSE'))))
	contingency.matrix.inc = matrix(rbind(this.inc.total,do.inc.total - this.inc.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

	this.dec.total = as.integer(table(factor(x$direction < 0,levels=c('TRUE','FALSE'))))
	contingency.matrix.dec = matrix(rbind(this.dec.total,do.dec.total - this.dec.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

	deg.fet.test = fisher.test(contingency.matrix.deg,alternative='greater')
	inc.fet.test = fisher.test(contingency.matrix.inc,alternative='greater')
	dec.fet.test = fisher.test(contingency.matrix.dec,alternative='greater')
	inc.kst.test = ks.test(x$effect,subset(all.region.do.pass,do_id != i)$effect,alternative='less')
	dec.kst.test = ks.test(x$effect,subset(all.region.do.pass,do_id != i)$effect,alternative='greater')

	data.frame(
		do_id = unique(x$do_id),
		do.size = sum(this.deg.total),
		deg.n = this.deg.total[1],
		inc.n = this.inc.total[1],
		dec.n = this.dec.total[1],
		deg.fet.score = deg.fet.test$estimate,
		inc.fet.score = inc.fet.test$estimate,
		dec.fet.score = dec.fet.test$estimate,
		inc.kst.score = inc.kst.test$statistic,
		dec.kst.score = dec.kst.test$statistic,
		deg.fet.pval = deg.fet.test$p.value,
		inc.fet.pval = inc.fet.test$p.value,
		dec.fet.pval = dec.fet.test$p.value,
		inc.kst.pval = inc.kst.test$p.value,
		dec.kst.pval = dec.kst.test$p.value
	)
},mc.cores=n.cores))

all.region.do.test = within(all.region.do.test,{
	dec.kst.qval = p.adjust(dec.kst.pval,'fdr')
	inc.kst.qval = p.adjust(inc.kst.pval,'fdr')
	dec.fet.qval = p.adjust(dec.fet.pval,'fdr')
	inc.fet.qval = p.adjust(inc.fet.pval,'fdr')
	deg.fet.qval = p.adjust(deg.fet.pval,'fdr')
})

all.region.do.results = merge(all.region.do.test,do.def,by='do_id')
all.region.do.results$dataset = 'DISEASES'
all.region.do.results$set = 'union'
all.region.do.results$region = 'all'


# Now repeat for DisGeNET

# library(disgenet2r)

if (ignore.checkpoints || !file.exists('checkpoints/disgenet_human_entrez_ids.rds')) {
	library(biomaRt)

	hsap = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
		dataset='hsapiens_gene_ensembl')

	hsap.entrez = getBM(
		attributes=c('ensembl_gene_id','entrezgene_id'),
		mart = hsap)

	saveRDS(hsap.entrez,file='checkpoints/disgenet_human_entrez_ids.rds')
} else {
	message('Checkpoint found!\nLoading human ortholog annotations from file.')

	hsap.entrez = readRDS('checkpoints/disgenet_human_entrez_ids.rds')
}

hsap.entrez = merge(hsap.entrez,unique(subset(hsap.info,mmulatta_homolog_orthology_type == 'ortholog_one2one' & mmulatta_homolog_ensembl_gene %in% ensembl.gene.names,select=c('ensembl_gene_id','external_gene_name','mmulatta_homolog_ensembl_gene','mmulatta_homolog_orthology_type'))),by='ensembl_gene_id')

dg.data = read.table(
	'data/curated_gene_disease_associations.tsv',
	sep='\t',
	quote='',
	header = TRUE,
	stringsAsFactors=FALSE)

dg.def = unique(subset(dg.data,select=c('diseaseId','diseaseName','diseaseType','diseaseClass','diseaseSemanticType')))
names(dg.def) = c('dg_id','dg_name','dg_type','dg_class','dg_semantic_type')

dg.mmul = merge(dg.data,hsap.entrez,by.x='geneId',by.y='entrezgene_id')

dg.mmul = subset(dg.mmul,
	select=c('mmulatta_homolog_ensembl_gene','external_gene_name','ensembl_gene_id','diseaseId','DSI','DPI','diseaseName','diseaseType','diseaseClass','diseaseSemanticType','score','EI'))

names(dg.mmul) = c('ensembl_gene_id','external_gene_name','hsapiens_homolog_ensembl_gene','dg_id','dsi','dpi','dg_name','dg_type','dg_class','dg_semantic_type','dg_score','ei')

dg.mmul.out = subset(dg.mmul,ensembl_gene_id %in% ensembl.gene.names,select=c('ensembl_gene_id','dg_id','dg_name','dg_class','dg_score'))
rownames(dg.mmul.out) = NULL
saveRDS(dg.mmul.out,file='checkpoints/rnaseq_genes_dg.rds')


all.region.dg = merge(all.region.join, dg.mmul, by='ensembl_gene_id')

# Toggle number
# Restrict to F branch
all.region.dg.pass = subset(all.region.dg,dg_id %in% names(which(table(subset(all.region.dg,dg_score >= 0 & grepl('F03',dg_class))$dg_id) >= 10)))

all.region.dg.gene.pass = unique(all.region.dg.pass[c('ensembl_gene_id','direction')])

dg.deg.total = as.integer(table(factor(all.region.dg.gene.pass$direction != 0,levels=c('TRUE','FALSE'))))
dg.inc.total = as.integer(table(factor(all.region.dg.gene.pass$direction > 0,levels=c('TRUE','FALSE'))))
dg.dec.total = as.integer(table(factor(all.region.dg.gene.pass$direction < 0,levels=c('TRUE','FALSE'))))

all.region.dg.split = split(all.region.dg.pass,all.region.dg.pass$dg_id)

library(parallel)

all.region.dg.test = do.call(rbind,mclapply(names(all.region.dg.split),function(i) {
	x = all.region.dg.split[[i]]

	this.deg.total = as.integer(table(factor(x$direction != 0,levels=c('TRUE','FALSE'))))
	contingency.matrix.deg = matrix(rbind(this.deg.total,dg.deg.total - this.deg.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

	this.inc.total = as.integer(table(factor(x$direction > 0,levels=c('TRUE','FALSE'))))
	contingency.matrix.inc = matrix(rbind(this.inc.total,dg.inc.total - this.inc.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

	this.dec.total = as.integer(table(factor(x$direction < 0,levels=c('TRUE','FALSE'))))
	contingency.matrix.dec = matrix(rbind(this.dec.total,dg.dec.total - this.dec.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

	deg.fet.test = fisher.test(contingency.matrix.deg,alternative='greater')
	inc.fet.test = fisher.test(contingency.matrix.inc,alternative='greater')
	dec.fet.test = fisher.test(contingency.matrix.dec,alternative='greater')
	inc.kst.test = ks.test(x$effect,subset(all.region.dg.pass,dg_id != i)$effect,alternative='less')
	dec.kst.test = ks.test(x$effect,subset(all.region.dg.pass,dg_id != i)$effect,alternative='greater')

	data.frame(
		dg_id = unique(x$dg_id),
		dg.size = sum(this.deg.total),
		deg.n = this.deg.total[1],
		inc.n = this.inc.total[1],
		dec.n = this.dec.total[1],
		deg.fet.score = deg.fet.test$estimate,
		inc.fet.score = inc.fet.test$estimate,
		dec.fet.score = dec.fet.test$estimate,
		inc.kst.score = inc.kst.test$statistic,
		dec.kst.score = dec.kst.test$statistic,
		deg.fet.pval = deg.fet.test$p.value,
		inc.fet.pval = inc.fet.test$p.value,
		dec.fet.pval = dec.fet.test$p.value,
		inc.kst.pval = inc.kst.test$p.value,
		dec.kst.pval = dec.kst.test$p.value
	)

},mc.cores=n.cores))

all.region.dg.test = within(all.region.dg.test,{
	dec.kst.qval = p.adjust(dec.kst.pval,'fdr')
	inc.kst.qval = p.adjust(inc.kst.pval,'fdr')
	dec.fet.qval = p.adjust(dec.fet.pval,'fdr')
	inc.fet.qval = p.adjust(inc.fet.pval,'fdr')
	deg.fet.qval = p.adjust(deg.fet.pval,'fdr')
})

all.region.dg.results = merge(all.region.dg.test,dg.def,by='dg_id')
all.region.dg.results$dataset = 'DisGeNET'
all.region.dg.results$set = 'union'
all.region.dg.results$region = 'all'

########################################################################
# One region
########################################################################

one.region.do.results = do.call(rbind,mclapply(names(keep.genes),function(i) {
	this.region.fet = this.region.kst = numeric(length=length(ensembl.gene.names))
	names(this.region.fet) = names(this.region.kst) = ensembl.gene.names

	# Code upregulated genes as 1 and downregulated genes as -1 if they are significant and co-directional in at most a fraction [fraction.unique.cutoff] of regions and up- or down-regulated in the current region
	this.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
		mash.lfsr[x,i] < fsr.cutoff & mash.beta[x,i] > 0
	}))))] = 1
	this.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
		mash.lfsr[x,i] < fsr.cutoff & mash.beta[x,i] < 0
	}))))] = -1

	this.region.kst[rownames(mash.sbet)] = mash.sbet[,i]

	this.region.join = data.frame(ensembl_gene_id = ensembl.gene.names, direction = as.integer(this.region.fet), effect = as.numeric(this.region.kst))

	this.region.do = merge(this.region.join, do.mmul, by='ensembl_gene_id')

	# Toggle number
	this.region.do.pass = subset(this.region.do,do_id %in% names(which(table(subset(this.region.do,confidence >= 0)$do_id) >= 10)))

	this.region.do.gene.pass = unique(this.region.do.pass[c('ensembl_gene_id','direction')])

	do.deg.total = as.integer(table(factor(this.region.do.gene.pass$direction != 0,levels=c('TRUE','FALSE'))))
	do.inc.total = as.integer(table(factor(this.region.do.gene.pass$direction > 0,levels=c('TRUE','FALSE'))))
	do.dec.total = as.integer(table(factor(this.region.do.gene.pass$direction < 0,levels=c('TRUE','FALSE'))))

	this.region.do.split = split(this.region.do.pass,this.region.do.pass$do_id)

	this.region.do.test = do.call(rbind,lapply(names(this.region.do.split),function(j) {
		x = this.region.do.split[[j]]

		this.deg.total = as.integer(table(factor(x$direction != 0,levels=c('TRUE','FALSE'))))
		contingency.matrix.deg = matrix(rbind(this.deg.total,do.deg.total - this.deg.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

		this.inc.total = as.integer(table(factor(x$direction > 0,levels=c('TRUE','FALSE'))))
		contingency.matrix.inc = matrix(rbind(this.inc.total,do.inc.total - this.inc.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

		this.dec.total = as.integer(table(factor(x$direction < 0,levels=c('TRUE','FALSE'))))
		contingency.matrix.dec = matrix(rbind(this.dec.total,do.dec.total - this.dec.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

		deg.fet.test = fisher.test(contingency.matrix.deg,alternative='greater')
		inc.fet.test = fisher.test(contingency.matrix.inc,alternative='greater')
		dec.fet.test = fisher.test(contingency.matrix.dec,alternative='greater')
		inc.kst.test = ks.test(x$effect,subset(this.region.do.pass,do_id != i)$effect,alternative='less')
		dec.kst.test = ks.test(x$effect,subset(this.region.do.pass,do_id != i)$effect,alternative='greater')

		data.frame(
			do_id = unique(x$do_id),
			do.size = sum(this.deg.total),
			deg.n = this.deg.total[1],
			inc.n = this.inc.total[1],
			dec.n = this.dec.total[1],
			deg.fet.score = deg.fet.test$estimate,
			inc.fet.score = inc.fet.test$estimate,
			dec.fet.score = dec.fet.test$estimate,
			inc.kst.score = inc.kst.test$statistic,
			dec.kst.score = dec.kst.test$statistic,
			deg.fet.pval = deg.fet.test$p.value,
			inc.fet.pval = inc.fet.test$p.value,
			dec.fet.pval = dec.fet.test$p.value,
			inc.kst.pval = inc.kst.test$p.value,
			dec.kst.pval = dec.kst.test$p.value
		)
	}))

	this.region.do.test = within(this.region.do.test,{
		dec.kst.qval = p.adjust(dec.kst.pval,'fdr')
		inc.kst.qval = p.adjust(inc.kst.pval,'fdr')
		dec.fet.qval = p.adjust(dec.fet.pval,'fdr')
		inc.fet.qval = p.adjust(inc.fet.pval,'fdr')
		deg.fet.qval = p.adjust(deg.fet.pval,'fdr')
	})

	this.region.do.test = merge(this.region.do.test,do.def,by='do_id')
	this.region.do.test$dataset = 'DISEASES'
	this.region.do.test$set = 'one'
	this.region.do.test$region = i
	this.region.do.test

},mc.cores=n.cores))

one.region.dg.results = do.call(rbind,mclapply(names(keep.genes),function(i) {
	this.region.fet = this.region.kst = numeric(length=length(ensembl.gene.names))
	names(this.region.fet) = names(this.region.kst) = ensembl.gene.names

	# Code upregulated genes as 1 and downregulated genes as -1 if they are significant and co-directional in at most a fraction [fraction.unique.cutoff] of regions and up- or down-regulated in the current region
	this.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
		mash.lfsr[x,i] < fsr.cutoff & mash.beta[x,i] > 0
	}))))] = 1
	this.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
		mash.lfsr[x,i] < fsr.cutoff & mash.beta[x,i] < 0
	}))))] = -1

	this.region.kst[rownames(mash.sbet)] = mash.sbet[,i]

	this.region.join = data.frame(ensembl_gene_id = ensembl.gene.names, direction = as.integer(this.region.fet), effect = as.numeric(this.region.kst))

	this.region.dg = merge(this.region.join, dg.mmul, by='ensembl_gene_id')

	# Toggle number
	this.region.dg.pass = subset(this.region.dg,dg_id %in% names(which(table(subset(this.region.dg,dg_score >= 0 & grepl('F',dg_class))$dg_id) >= 10)))

	this.region.dg.gene.pass = unique(this.region.dg.pass[c('ensembl_gene_id','direction')])

	dg.deg.total = as.integer(table(factor(this.region.dg.gene.pass$direction != 0,levels=c('TRUE','FALSE'))))
	dg.inc.total = as.integer(table(factor(this.region.dg.gene.pass$direction > 0,levels=c('TRUE','FALSE'))))
	dg.dec.total = as.integer(table(factor(this.region.dg.gene.pass$direction < 0,levels=c('TRUE','FALSE'))))

	this.region.dg.split = split(this.region.dg.pass,this.region.dg.pass$dg_id)

	this.region.dg.test = do.call(rbind,lapply(names(this.region.dg.split),function(j) {
		x = this.region.dg.split[[j]]

		this.deg.total = as.integer(table(factor(x$direction != 0,levels=c('TRUE','FALSE'))))
		contingency.matrix.deg = matrix(rbind(this.deg.total,dg.deg.total - this.deg.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

		this.inc.total = as.integer(table(factor(x$direction > 0,levels=c('TRUE','FALSE'))))
		contingency.matrix.inc = matrix(rbind(this.inc.total,dg.inc.total - this.inc.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

		this.dec.total = as.integer(table(factor(x$direction < 0,levels=c('TRUE','FALSE'))))
		contingency.matrix.dec = matrix(rbind(this.dec.total,dg.dec.total - this.dec.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

		deg.fet.test = fisher.test(contingency.matrix.deg,alternative='greater')
		inc.fet.test = fisher.test(contingency.matrix.inc,alternative='greater')
		dec.fet.test = fisher.test(contingency.matrix.dec,alternative='greater')
		inc.kst.test = ks.test(x$effect,subset(this.region.dg.pass,dg_id != i)$effect,alternative='less')
		dec.kst.test = ks.test(x$effect,subset(this.region.dg.pass,dg_id != i)$effect,alternative='greater')

		data.frame(
			dg_id = unique(x$dg_id),
			dg.size = sum(this.deg.total),
			deg.n = this.deg.total[1],
			inc.n = this.inc.total[1],
			dec.n = this.dec.total[1],
			deg.fet.score = deg.fet.test$estimate,
			inc.fet.score = inc.fet.test$estimate,
			dec.fet.score = dec.fet.test$estimate,
			inc.kst.score = inc.kst.test$statistic,
			dec.kst.score = dec.kst.test$statistic,
			deg.fet.pval = deg.fet.test$p.value,
			inc.fet.pval = inc.fet.test$p.value,
			dec.fet.pval = dec.fet.test$p.value,
			inc.kst.pval = inc.kst.test$p.value,
			dec.kst.pval = dec.kst.test$p.value
		)
	}))

	this.region.dg.test = within(this.region.dg.test,{
		dec.kst.qval = p.adjust(dec.kst.pval,'fdr')
		inc.kst.qval = p.adjust(inc.kst.pval,'fdr')
		dec.fet.qval = p.adjust(dec.fet.pval,'fdr')
		inc.fet.qval = p.adjust(inc.fet.pval,'fdr')
		deg.fet.qval = p.adjust(deg.fet.pval,'fdr')
	})

	this.region.dg.test = merge(this.region.dg.test,dg.def,by='dg_id')
	this.region.dg.test$dataset = 'DisGeNET'
	this.region.dg.test$set = 'one'
	this.region.dg.test$region = i
	this.region.dg.test

},mc.cores=n.cores))

########################################################################
# Single region
########################################################################

single.region.do.results = do.call(rbind,mclapply(names(keep.genes),function(i) {
	this.region.fet = this.region.kst = numeric(length=length(ensembl.gene.names))
	names(this.region.fet) = names(this.region.kst) = ensembl.gene.names

	# Code upregulated genes as 1 and downregulated genes as -1 if they are significant and co-directional in at most a fraction [fraction.unique.cutoff] of regions and up- or down-regulated in the current region
	this.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
		(sum((mash.lfsr[x,] < fsr.cutoff) & (mash.beta[x,] * mash.beta[x,i] > 0)) <= fraction.unique.cutoff * length(keep.genes)) && mash.lfsr[x,i] < fsr.cutoff && mash.beta[x,i] > 0
	}))))] = 1
	this.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
		(sum((mash.lfsr[x,] < fsr.cutoff) & (mash.beta[x,] * mash.beta[x,i] > 0)) <= fraction.unique.cutoff * length(keep.genes)) && mash.lfsr[x,i] < fsr.cutoff && mash.beta[x,i] < 0
	}))))] = -1

	this.region.kst[rownames(mash.sbet)] = mash.sbet[,i]

	this.region.join = data.frame(ensembl_gene_id = ensembl.gene.names, direction = as.integer(this.region.fet), effect = as.numeric(this.region.kst))

	this.region.do = merge(this.region.join, do.mmul, by='ensembl_gene_id')

	# Toggle number
	this.region.do.pass = subset(this.region.do,do_id %in% names(which(table(subset(this.region.do,confidence >= 0)$do_id) >= 10)))

	this.region.do.gene.pass = unique(this.region.do.pass[c('ensembl_gene_id','direction')])

	do.deg.total = as.integer(table(factor(this.region.do.gene.pass$direction != 0,levels=c('TRUE','FALSE'))))
	do.inc.total = as.integer(table(factor(this.region.do.gene.pass$direction > 0,levels=c('TRUE','FALSE'))))
	do.dec.total = as.integer(table(factor(this.region.do.gene.pass$direction < 0,levels=c('TRUE','FALSE'))))

	this.region.do.split = split(this.region.do.pass,this.region.do.pass$do_id)

	this.region.do.test = do.call(rbind,lapply(names(this.region.do.split),function(j) {
		x = this.region.do.split[[j]]

		this.deg.total = as.integer(table(factor(x$direction != 0,levels=c('TRUE','FALSE'))))
		contingency.matrix.deg = matrix(rbind(this.deg.total,do.deg.total - this.deg.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

		this.inc.total = as.integer(table(factor(x$direction > 0,levels=c('TRUE','FALSE'))))
		contingency.matrix.inc = matrix(rbind(this.inc.total,do.inc.total - this.inc.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

		this.dec.total = as.integer(table(factor(x$direction < 0,levels=c('TRUE','FALSE'))))
		contingency.matrix.dec = matrix(rbind(this.dec.total,do.dec.total - this.dec.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

		deg.fet.test = fisher.test(contingency.matrix.deg,alternative='greater')
		inc.fet.test = fisher.test(contingency.matrix.inc,alternative='greater')
		dec.fet.test = fisher.test(contingency.matrix.dec,alternative='greater')
		inc.kst.test = ks.test(x$effect,subset(this.region.do.pass,do_id != i)$effect,alternative='less')
		dec.kst.test = ks.test(x$effect,subset(this.region.do.pass,do_id != i)$effect,alternative='greater')

		data.frame(
			do_id = unique(x$do_id),
			do.size = sum(this.deg.total),
			deg.n = this.deg.total[1],
			inc.n = this.inc.total[1],
			dec.n = this.dec.total[1],
			deg.fet.score = deg.fet.test$estimate,
			inc.fet.score = inc.fet.test$estimate,
			dec.fet.score = dec.fet.test$estimate,
			inc.kst.score = inc.kst.test$statistic,
			dec.kst.score = dec.kst.test$statistic,
			deg.fet.pval = deg.fet.test$p.value,
			inc.fet.pval = inc.fet.test$p.value,
			dec.fet.pval = dec.fet.test$p.value,
			inc.kst.pval = inc.kst.test$p.value,
			dec.kst.pval = dec.kst.test$p.value
		)
	}))

	this.region.do.test = within(this.region.do.test,{
		dec.kst.qval = p.adjust(dec.kst.pval,'fdr')
		inc.kst.qval = p.adjust(inc.kst.pval,'fdr')
		dec.fet.qval = p.adjust(dec.fet.pval,'fdr')
		inc.fet.qval = p.adjust(inc.fet.pval,'fdr')
		deg.fet.qval = p.adjust(deg.fet.pval,'fdr')
	})

	this.region.do.test = merge(this.region.do.test,do.def,by='do_id')
	this.region.do.test$dataset = 'DISEASES'
	this.region.do.test$set = 'single'
	this.region.do.test$region = i
	this.region.do.test

},mc.cores=n.cores))

single.region.dg.results = do.call(rbind,mclapply(names(keep.genes),function(i) {
	this.region.fet = this.region.kst = numeric(length=length(ensembl.gene.names))
	names(this.region.fet) = names(this.region.kst) = ensembl.gene.names

	# Code upregulated genes as 1 and downregulated genes as -1 if they are significant and co-directional in at most a fraction [fraction.unique.cutoff] of regions and up- or down-regulated in the current region
	this.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
		(sum((mash.lfsr[x,] < fsr.cutoff) & (mash.beta[x,] * mash.beta[x,i] > 0)) <= fraction.unique.cutoff * length(keep.genes)) && mash.lfsr[x,i] < fsr.cutoff && mash.beta[x,i] > 0
	}))))] = 1
	this.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
		(sum((mash.lfsr[x,] < fsr.cutoff) & (mash.beta[x,] * mash.beta[x,i] > 0)) <= fraction.unique.cutoff * length(keep.genes)) && mash.lfsr[x,i] < fsr.cutoff && mash.beta[x,i] < 0
	}))))] = -1

	this.region.kst[rownames(mash.sbet)] = mash.sbet[,i]

	this.region.join = data.frame(ensembl_gene_id = ensembl.gene.names, direction = as.integer(this.region.fet), effect = as.numeric(this.region.kst))

	this.region.dg = merge(this.region.join, dg.mmul, by='ensembl_gene_id')

	# Toggle number
	this.region.dg.pass = subset(this.region.dg,dg_id %in% names(which(table(subset(this.region.dg,dg_score >= 0 & grepl('F',dg_class))$dg_id) >= 10)))

	this.region.dg.gene.pass = unique(this.region.dg.pass[c('ensembl_gene_id','direction')])

	dg.deg.total = as.integer(table(factor(this.region.dg.gene.pass$direction != 0,levels=c('TRUE','FALSE'))))
	dg.inc.total = as.integer(table(factor(this.region.dg.gene.pass$direction > 0,levels=c('TRUE','FALSE'))))
	dg.dec.total = as.integer(table(factor(this.region.dg.gene.pass$direction < 0,levels=c('TRUE','FALSE'))))

	this.region.dg.split = split(this.region.dg.pass,this.region.dg.pass$dg_id)

	this.region.dg.test = do.call(rbind,lapply(names(this.region.dg.split),function(j) {
		x = this.region.dg.split[[j]]

		this.deg.total = as.integer(table(factor(x$direction != 0,levels=c('TRUE','FALSE'))))
		contingency.matrix.deg = matrix(rbind(this.deg.total,dg.deg.total - this.deg.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

		this.inc.total = as.integer(table(factor(x$direction > 0,levels=c('TRUE','FALSE'))))
		contingency.matrix.inc = matrix(rbind(this.inc.total,dg.inc.total - this.inc.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

		this.dec.total = as.integer(table(factor(x$direction < 0,levels=c('TRUE','FALSE'))))
		contingency.matrix.dec = matrix(rbind(this.dec.total,dg.dec.total - this.dec.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

		deg.fet.test = fisher.test(contingency.matrix.deg,alternative='greater')
		inc.fet.test = fisher.test(contingency.matrix.inc,alternative='greater')
		dec.fet.test = fisher.test(contingency.matrix.dec,alternative='greater')
		inc.kst.test = ks.test(x$effect,subset(this.region.dg.pass,dg_id != i)$effect,alternative='less')
		dec.kst.test = ks.test(x$effect,subset(this.region.dg.pass,dg_id != i)$effect,alternative='greater')

		data.frame(
			dg_id = unique(x$dg_id),
			dg.size = sum(this.deg.total),
			deg.n = this.deg.total[1],
			inc.n = this.inc.total[1],
			dec.n = this.dec.total[1],
			deg.fet.score = deg.fet.test$estimate,
			inc.fet.score = inc.fet.test$estimate,
			dec.fet.score = dec.fet.test$estimate,
			inc.kst.score = inc.kst.test$statistic,
			dec.kst.score = dec.kst.test$statistic,
			deg.fet.pval = deg.fet.test$p.value,
			inc.fet.pval = inc.fet.test$p.value,
			dec.fet.pval = dec.fet.test$p.value,
			inc.kst.pval = inc.kst.test$p.value,
			dec.kst.pval = dec.kst.test$p.value
		)
	}))

	this.region.dg.test = within(this.region.dg.test,{
		dec.kst.qval = p.adjust(dec.kst.pval,'fdr')
		inc.kst.qval = p.adjust(inc.kst.pval,'fdr')
		dec.fet.qval = p.adjust(dec.fet.pval,'fdr')
		inc.fet.qval = p.adjust(inc.fet.pval,'fdr')
		deg.fet.qval = p.adjust(deg.fet.pval,'fdr')
	})

	this.region.dg.test = merge(this.region.dg.test,dg.def,by='dg_id')
	this.region.dg.test$dataset = 'DisGeNET'
	this.region.dg.test$set = 'single'
	this.region.dg.test$region = i
	this.region.dg.test

},mc.cores=n.cores))

all.do.results = rbind(all.region.do.results,one.region.do.results,single.region.do.results)
all.dg.results = rbind(all.region.dg.results,one.region.dg.results,single.region.dg.results)

# Melt and rename columns to make the structure similar to the GO results
d.x1 = data.frame(subset(all.do.results,select=c('do_id','dataset','do_name','inc.fet.score','inc.fet.pval','inc.fet.qval','set','region')),direction='increase',test='FET')
d.x2 = data.frame(subset(all.do.results,select=c('do_id','dataset','do_name','dec.fet.score','dec.fet.pval','dec.fet.qval','set','region')),direction='decrease',test='FET')
d.x3 = data.frame(subset(all.do.results,select=c('do_id','dataset','do_name','inc.kst.score','inc.kst.pval','inc.kst.qval','set','region')),direction='increase',test='KS')
d.x4 = data.frame(subset(all.do.results,select=c('do_id','dataset','do_name','dec.kst.score','dec.kst.pval','dec.kst.qval','set','region')),direction='decrease',test='KS')
d.y1 = data.frame(subset(all.dg.results,select=c('dg_id','dataset','dg_name','inc.fet.score','inc.fet.pval','inc.fet.qval','set','region')),direction='increase',test='FET')
d.y2 = data.frame(subset(all.dg.results,select=c('dg_id','dataset','dg_name','dec.fet.score','dec.fet.pval','dec.fet.qval','set','region')),direction='decrease',test='FET')
d.y3 = data.frame(subset(all.dg.results,select=c('dg_id','dataset','dg_name','inc.kst.score','inc.kst.pval','inc.kst.qval','set','region')),direction='increase',test='KS')
d.y4 = data.frame(subset(all.dg.results,select=c('dg_id','dataset','dg_name','dec.kst.score','dec.kst.pval','dec.kst.qval','set','region')),direction='decrease',test='KS')

names(d.x1) = names(d.x2) = names(d.x3) = names(d.x4) = names(d.y1) = names(d.y2) = names(d.y3) = names(d.y4) = c('do_id','dataset','do_name','score','pval','qval','set','region','direction','test')

disease.enrichment.results = subset(rbind(
	d.x1,
	d.x2,
	d.x3,
	d.x4,
	d.y1,
	d.y2,
	d.y3,
	d.y4
),select=c('do_id','dataset','do_name','score','pval','qval','direction','set','region','test'))

saveRDS(disease.enrichment.results,file='checkpoints/disease_enrichment_results.rds')

get.region.specificity = function(enrichment.results,annotation.dataset,tissue.set='one',effect.direction,statistical.test) {
	this.df = subset(enrichment.results,dataset==annotation.dataset & direction == effect.direction & set == tissue.set & test == statistical.test)
	this.df = subset(this.df,complete.cases(this.df))
	this.wide = tidyr::pivot_wider(droplevels(this.df),id_cols='do_id',names_from='region',values_from='score')
	this.mat = matrix(as.matrix(this.wide[,2:ncol(this.wide)]),nrow=nrow(this.wide),ncol=ncol(this.wide)-1,dimnames=list(this.wide$do_id,colnames(this.wide)[2:ncol(this.wide)]))
# 	apply(this.mat,1,function(x) max(x) - quantile(x,0.75))
	specificity = apply(this.mat,1,function(x) sort(x)[length(x)] - sort(x)[length(x)-1])
	specificity.null = do.call(cbind,parallel::mclapply(1:n.permutations,function(i) {
		apply(matrix(apply(this.mat,2,sample),nrow=nrow(this.mat),dimnames=dimnames(this.mat)),1,function(x) sort(x)[length(x)] - sort(x)[length(x)-1])
	},mc.cores=n.cores))
	data.frame(
		do_id = rownames(this.mat),
		dataset = annotation.dataset,
		do_name = as.character(unique(subset(this.df,select=c('do_id','do_name')))$do_name),
		direction = effect.direction,
		test = statistical.test,
		top.region = colnames(this.mat)[apply(this.mat,1,which.max)],
		specificity,
		specificity.pval = rowMeans(specificity < specificity.null),
		specificity.qval = p.adjust(rowMeans(specificity < specificity.null),'fdr')
	)
}

do.fet.inc.spec = get.region.specificity(disease.enrichment.results,'DISEASES','one','increase','FET')
do.fet.dec.spec = get.region.specificity(disease.enrichment.results,'DISEASES','one','decrease','FET')
do.kst.inc.spec = get.region.specificity(disease.enrichment.results,'DISEASES','one','increase','KS')
do.kst.dec.spec = get.region.specificity(disease.enrichment.results,'DISEASES','one','decrease','KS')
dg.fet.inc.spec = get.region.specificity(disease.enrichment.results,'DisGeNET','one','increase','FET')
dg.fet.dec.spec = get.region.specificity(disease.enrichment.results,'DisGeNET','one','decrease','FET')
dg.kst.inc.spec = get.region.specificity(disease.enrichment.results,'DisGeNET','one','increase','KS')
dg.kst.dec.spec = get.region.specificity(disease.enrichment.results,'DisGeNET','one','decrease','KS')

disease.specificity.results = rbind(
	do.fet.inc.spec,
	do.fet.dec.spec,
	do.kst.inc.spec,
	do.kst.dec.spec,
	dg.fet.inc.spec,
	dg.fet.dec.spec,
	dg.kst.inc.spec,
	dg.kst.dec.spec
)

saveRDS(disease.specificity.results,file='checkpoints/disease_specificity_results.rds')

save(list=c('all.do.results','all.dg.results'),file='checkpoints/disease_enrichment_results_separate.RData')