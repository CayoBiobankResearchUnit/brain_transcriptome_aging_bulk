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

if (ignore.checkpoints || !file.exists('checkpoints/rnaseq_genes_go.rds')) {
	library(biomaRt)

	mmul = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL',
		dataset='mmulatta_gene_ensembl')

	mmul.go = getBM(
		attributes=c('ensembl_gene_id','go_id'),
		filters = 'ensembl_gene_id',
		values = ensembl.gene.names,
		mart = mmul)
	saveRDS(mmul.go,file='checkpoints/rnaseq_genes_go.rds')
} else {
	message('Checkpoint found!\nLoading GO annotations from file.')

	mmul.go = readRDS('checkpoints/rnaseq_genes_go.rds')
}
gene2go = lapply(
	unique(mmul.go$ensembl_gene_id),
	function(x) {
		out = sort(mmul.go[mmul.go$ensembl_gene_id == x,'go_id'])
		out[out != '']
	}
)

names(gene2go) = unique(mmul.go$ensembl_gene_id)

library(topGO)

# Fetch names for GO terms identified as relations by topGO
if (ignore.checkpoints || !file.exists('checkpoints/rnaseq_go_names.rds')) {
	# Create dummy vector for topGO
	all.genes = numeric(length(gene2go))
	names(all.genes) = names(gene2go)

	# Initialize topGO object
	go.topgo = new('topGOdata',description='Simple session',ontology='BP',allGenes=all.genes,geneSelectionFun=function(x) x>0,nodeSize = 10,annotationFun = annFUN.gene2GO,gene2GO = gene2go)

	# Put together all relevant GO terms
	go.ids = sort(union(mmul.go$go_id,go.topgo@graph@nodes))

	# Exclude blank GO terms
	go.ids = go.ids[nchar(go.ids) == 10]

	# Grab names and namespaces from the GO database
	library(GO.db)

	go.terms = data.frame(
		go_id = go.ids,
		go_namespace = Ontology(go.ids),
		go_name = Term(go.ids),
		stringsAsFactors=FALSE
	)

	# Sometimes, there are names that are not in the GO database. Fetch these names from AmiGO

	library(XML)

	to.do = go.terms$go_id[!complete.cases(go.terms)]

	if (length(to.do)) {
		# Fetch missing GO metadata from AmiGO
		for (i in to.do) {
			message(i)
			search.url = paste0('http://amigo.geneontology.org/amigo/term/',i)

			amigo.data = xpathApply(htmlParse(search.url), '//dd', xmlValue)
			go.terms$go_namespace[match(i,go.terms$go_id)] = unlist(lapply(strsplit(amigo.data[[3]],'_'),function(x) paste(toupper(substr(x,1,1)),collapse='')))
			go.terms$go_name[match(i,go.terms$go_id)] = amigo.data[[2]]
		}
	}

	to.do = go.terms$go_id[!nchar(go.terms$go_name)]
	
	library(jsonlite)

	if (length(to.do)) {
		# Fetch still missing GO metadata from QuickGO
		for (i in to.do) {
			message(i)
			search.url = paste0('https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/',i)

			quickgo.data = read_json(search.url)
			go.terms$go_namespace[match(i,go.terms$go_id)] = c('molecular_function'='MF','biological_process'='BP','cellular_component'='CC')[quickgo.data$results[[1]]$aspect]
			go.terms$go_name[match(i,go.terms$go_id)] = quickgo.data$results[[1]]$name
		}
	}

	rownames(go.terms) = go.terms$go_id
	saveRDS(go.terms,file='checkpoints/rnaseq_go_names.rds')
} else {
	# If checkpoint is found, use it
	message('Checkpoint found!\nLoading GO metadata from file.')

	go.terms = readRDS('checkpoints/rnaseq_go_names.rds')
}

# Initialize results vector
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

# Initialize topGO data

all.region.fet.inc.topgo = new('topGOdata',
	description='Simple session',
	ontology='BP',
	allGenes=all.region.fet,
	geneSelectionFun=function(x) x > 0,
	nodeSize = 10,
	annotationFun = annFUN.gene2GO,
	gene2GO = gene2go
)

all.region.fet.dec.topgo = new('topGOdata',
	description='Simple session',
	ontology='BP',
	allGenes=all.region.fet,
	geneSelectionFun=function(x) x < 0,
	nodeSize = 10,
	annotationFun = annFUN.gene2GO,
	gene2GO = gene2go
)

all.region.kst.inc.topgo = new('topGOdata',
	description='Simple session',
	ontology='BP',
	allGenes=-all.region.kst,
	geneSelectionFun=function(x) x > 0,
	nodeSize = 10,
	annotationFun = annFUN.gene2GO,
	gene2GO = gene2go
)

all.region.kst.dec.topgo = new('topGOdata',
	description='Simple session',
	ontology='BP',
	allGenes=all.region.kst,
	geneSelectionFun=function(x) x < 0,
	nodeSize = 10,
	annotationFun = annFUN.gene2GO,
	gene2GO = gene2go
)

# Rewrite the following S4 methods to output the test statistics

# getMethod('GOKSTest',signature=c(object='classicScore'))
# getMethod('GOFisherTest',signature=c(object='classicCount'))

# But first ensure that the original methods are run first
setMethod('GOKSTest',signature(c(object='classicScore')),get.ks.pval)
setMethod('GOFisherTest',signature(c(object='classicCount')),get.fisher.pval)

# Run Fisher's exact test

all.region.fet.inc.test = runTest(all.region.fet.inc.topgo,algorithm='parentchild',statistic='fisher')
all.region.fet.dec.test = runTest(all.region.fet.dec.topgo,algorithm='parentchild',statistic='fisher')

# Run Kolmogorov–Smirnov test

all.region.kst.inc.test = runTest(all.region.kst.inc.topgo,algorithm='weight01',statistic='ks')
all.region.kst.dec.test = runTest(all.region.kst.dec.topgo,algorithm='weight01',statistic='ks')

library(reshape2)
all.region.fet.inc.results = melt(all.region.fet.inc.test@score)
all.region.fet.dec.results = melt(all.region.fet.dec.test@score)
all.region.kst.inc.results = melt(all.region.kst.inc.test@score)
all.region.kst.dec.results = melt(all.region.kst.dec.test@score)

# Write over the S4 methods for retrieving p-values from KS and Fisher tests.
# 	The topGO package outputs the p-values into the scores slot of the S4 results object.
# 	The methods for retrieving these scores are rewritten below to retrieve the model estimate instead
# 	Functions are in the `scripts/_include_options.R` script.
setMethod('GOKSTest',signature(c(object='classicScore')),get.ks.score)
setMethod('GOFisherTest',signature(c(object='classicCount')),get.fisher.score)

all.region.fet.inc.score = runTest(all.region.fet.inc.topgo,algorithm='classic',statistic='fisher')@score
all.region.fet.dec.score = runTest(all.region.fet.dec.topgo,algorithm='classic',statistic='fisher')@score
all.region.kst.inc.score = runTest(all.region.kst.inc.topgo,algorithm='weight01',statistic='ks')@score
all.region.kst.dec.score = runTest(all.region.kst.dec.topgo,algorithm='weight01',statistic='ks')@score

all.region.fet.inc.results = data.frame(
	go.terms[rownames(all.region.fet.inc.results),],
	score = all.region.fet.inc.score,
	pval = all.region.fet.inc.results$value,
	qval = p.adjust(all.region.fet.inc.results$value,'fdr'),
	direction = 'increase',
	set = 'union',
	region = 'all',
	test = 'FET',
	stringsAsFactors=FALSE
)

all.region.fet.dec.results = data.frame(
	go.terms[rownames(all.region.fet.dec.results),],
	score = all.region.fet.dec.score,
	pval = all.region.fet.dec.results$value,
	qval = p.adjust(all.region.fet.dec.results$value,'fdr'),
	direction = 'decrease',
	set = 'union',
	region = 'all',
	test = 'FET',
	stringsAsFactors=FALSE
)

all.region.kst.inc.results = data.frame(
	go.terms[rownames(all.region.kst.inc.results),],
	score = all.region.kst.inc.score,
	pval = all.region.kst.inc.results$value,
	qval = p.adjust(all.region.kst.inc.results$value,'fdr'),
	direction = 'increase',
	set = 'union',
	region = 'all',
	test = 'KS',
	stringsAsFactors=FALSE
)

all.region.kst.dec.results = data.frame(
	go.terms[rownames(all.region.kst.dec.results),],
	score = all.region.kst.dec.score,
	pval = all.region.kst.dec.results$value,
	qval = p.adjust(all.region.kst.dec.results$value,'fdr'),
	direction = 'decrease',
	set = 'union',
	region = 'all',
	test = 'KS',
	stringsAsFactors=FALSE
)

# Set topGO methods back to the originals (saving p-values into the scores S4 slot)
setMethod('GOKSTest',signature(c(object='classicScore')),get.ks.pval)
setMethod('GOFisherTest',signature(c(object='classicCount')),get.fisher.pval)

all.region.results = rbind(all.region.fet.inc.results,all.region.fet.dec.results,all.region.kst.inc.results,all.region.kst.dec.results)
rownames(all.region.results) = NULL

library(parallel)

# The function below was originally parallelized with mclapply. Turned this off to avoid conflicts in the writing/rewriting of topGO methods.
# one.region.results = do.call(rbind,mclapply(names(keep.genes),function(i) {
one.region.results = do.call(rbind,lapply(names(keep.genes),function(i) {
	this.region.fet = this.region.kst = numeric(length=length(ensembl.gene.names))
	names(this.region.fet) = names(this.region.kst) = ensembl.gene.names

	# Code upregulated genes as 1 and downregulated genes as -1 if they are significant and in the noted direction
	this.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
		mash.lfsr[x,i] < fsr.cutoff & mash.beta[x,i] > 0
	}))))] = 1
	this.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
		mash.lfsr[x,i] < fsr.cutoff & mash.beta[x,i] < 0
	}))))] = -1

	this.region.kst[rownames(mash.sbet)] = mash.sbet[,i]

	# Set topGO methods back to the originals (saving p-values into the scores S4 slot)
	setMethod('GOKSTest',signature(c(object='classicScore')),get.ks.pval)
	setMethod('GOFisherTest',signature(c(object='classicCount')),get.fisher.pval)

	this.region.fet.inc.topgo = new('topGOdata',
		description='Simple session',
		ontology='BP',
		allGenes=this.region.fet,
		geneSelectionFun=function(x) x > 0,
		nodeSize = 10,
		annotationFun = annFUN.gene2GO,
		gene2GO = gene2go
	)

	this.region.fet.dec.topgo = new('topGOdata',
		description='Simple session',
		ontology='BP',
		allGenes=this.region.fet,
		geneSelectionFun=function(x) x < 0,
		nodeSize = 10,
		annotationFun = annFUN.gene2GO,
		gene2GO = gene2go
	)

	this.region.kst.inc.topgo = new('topGOdata',
		description='Simple session',
		ontology='BP',
		allGenes=-this.region.kst,
		geneSelectionFun=function(x) x > 0,
		nodeSize = 10,
		annotationFun = annFUN.gene2GO,
		gene2GO = gene2go
	)

	this.region.kst.dec.topgo = new('topGOdata',
		description='Simple session',
		ontology='BP',
		allGenes=this.region.kst,
		geneSelectionFun=function(x) x < 0,
		nodeSize = 10,
		annotationFun = annFUN.gene2GO,
		gene2GO = gene2go
	)

	this.region.fet.inc.test = runTest(this.region.fet.inc.topgo,algorithm='parentchild',statistic='fisher')
	this.region.fet.dec.test = runTest(this.region.fet.dec.topgo,algorithm='parentchild',statistic='fisher')

	this.region.kst.inc.test = runTest(this.region.kst.inc.topgo,algorithm='weight01',statistic='ks')
	this.region.kst.dec.test = runTest(this.region.kst.dec.topgo,algorithm='weight01',statistic='ks')

	this.region.fet.inc.results = reshape2::melt(this.region.fet.inc.test@score)
	this.region.fet.dec.results = reshape2::melt(this.region.fet.dec.test@score)

	this.region.kst.inc.results = reshape2::melt(this.region.kst.inc.test@score)
	this.region.kst.dec.results = reshape2::melt(this.region.kst.dec.test@score)

	# Write over the S4 methods for retrieving p-values from KS and Fisher tests.
	setMethod('GOKSTest',signature(c(object='classicScore')),get.ks.score)
	setMethod('GOFisherTest',signature(c(object='classicCount')),get.fisher.score)

	this.region.fet.inc.score = runTest(this.region.fet.inc.topgo,algorithm='classic',statistic='fisher')@score
	this.region.fet.dec.score = runTest(this.region.fet.dec.topgo,algorithm='classic',statistic='fisher')@score
	this.region.kst.inc.score = runTest(this.region.kst.inc.topgo,algorithm='weight01',statistic='ks')@score
	this.region.kst.dec.score = runTest(this.region.kst.dec.topgo,algorithm='weight01',statistic='ks')@score

	this.region.fet.inc.results = data.frame(
		go.terms[rownames(this.region.fet.inc.results),],
		score = this.region.fet.inc.score,
		pval = this.region.fet.inc.results$value,
		qval = p.adjust(this.region.fet.inc.results$value,'fdr'),
		direction = 'increase',
		set = 'one',
		region = i,
		test = 'FET',
		stringsAsFactors=FALSE
	)

	this.region.fet.dec.results = data.frame(
		go.terms[rownames(this.region.fet.dec.results),],
		score = this.region.fet.dec.score,
		pval = this.region.fet.dec.results$value,
		qval = p.adjust(this.region.fet.dec.results$value,'fdr'),
		direction = 'decrease',
		set = 'one',
		region = i,
		test = 'FET',
		stringsAsFactors=FALSE
	)

	this.region.kst.inc.results = data.frame(
		go.terms[rownames(this.region.kst.inc.results),],
		score = this.region.kst.inc.score,
		pval = this.region.kst.inc.results$value,
		qval = p.adjust(this.region.kst.inc.results$value,'fdr'),
		direction = 'increase',
		set = 'one',
		region = i,
		test = 'KS',
		stringsAsFactors=FALSE
	)

	this.region.kst.dec.results = data.frame(
		go.terms[rownames(this.region.kst.dec.results),],
		score = this.region.kst.dec.score,
		pval = this.region.kst.dec.results$value,
		qval = p.adjust(this.region.kst.dec.results$value,'fdr'),
		direction = 'decrease',
		set = 'one',
		region = i,
		test = 'KS',
		stringsAsFactors=FALSE
	)
	rbind(this.region.fet.inc.results,this.region.fet.dec.results,this.region.kst.inc.results,this.region.kst.dec.results)
}))
# },mc.cores=n.cores))

rownames(one.region.results) = NULL

# # The function below was originally parallelized with mclapply. Turned this off to avoid conflicts in the writing/rewriting of topGO methods.
# # single.region.results = do.call(rbind,mclapply(names(keep.genes),function(i) {
# single.region.results = do.call(rbind,lapply(names(keep.genes),function(i) {
# 	this.region.fet = numeric(length=length(ensembl.gene.names))
# 	names(this.region.fet) = ensembl.gene.names
# 
# 	# Code upregulated genes as 1 and downregulated genes as -1 if they are significant and co-directional in at most a fraction [fraction.unique.cutoff] of regions and up- or down-regulated in the current region
# 	this.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
# 		(sum((mash.lfsr[x,] < fsr.cutoff) & (mash.beta[x,] * mash.beta[x,i] > 0)) <= fraction.unique.cutoff * length(keep.genes)) && mash.lfsr[x,i] < fsr.cutoff && mash.beta[x,i] > 0
# 	}))))] = 1
# 	this.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
# 		(sum((mash.lfsr[x,] < fsr.cutoff) & (mash.beta[x,] * mash.beta[x,i] > 0)) <= fraction.unique.cutoff * length(keep.genes)) && mash.lfsr[x,i] < fsr.cutoff && mash.beta[x,i] < 0
# 	}))))] = -1
# 
# 	# Set topGO methods back to the originals (saving p-values into the scores S4 slot)
# 	setMethod('GOFisherTest',signature(c(object='classicCount')),get.fisher.pval)
# 
# 	this.region.fet.inc.topgo = new('topGOdata',
# 		description='Simple session',
# 		ontology='BP',
# 		allGenes=this.region.fet,
# 		geneSelectionFun=function(x) x > 0,
# 		nodeSize = 10,
# 		annotationFun = annFUN.gene2GO,
# 		gene2GO = gene2go
# 	)
# 
# 	this.region.fet.dec.topgo = new('topGOdata',
# 		description='Simple session',
# 		ontology='BP',
# 		allGenes=this.region.fet,
# 		geneSelectionFun=function(x) x < 0,
# 		nodeSize = 10,
# 		annotationFun = annFUN.gene2GO,
# 		gene2GO = gene2go
# 	)
# 
# 	this.region.fet.inc.test = runTest(this.region.fet.inc.topgo,algorithm='parentchild',statistic='fisher')
# 	this.region.fet.dec.test = runTest(this.region.fet.dec.topgo,algorithm='parentchild',statistic='fisher')
# 
# 	this.region.fet.inc.results = reshape2::melt(this.region.fet.inc.test@score)
# 	this.region.fet.dec.results = reshape2::melt(this.region.fet.dec.test@score)
# 
# 	# Write over the S4 methods for retrieving p-values from KS and Fisher tests.
# 	setMethod('GOFisherTest',signature(c(object='classicCount')),get.fisher.score)
# 
# 	this.region.fet.inc.score = runTest(this.region.fet.inc.topgo,algorithm='classic',statistic='fisher')@score
# 	this.region.fet.dec.score = runTest(this.region.fet.dec.topgo,algorithm='classic',statistic='fisher')@score
# 
# 	this.region.fet.inc.results = data.frame(
# 		go.terms[rownames(this.region.fet.inc.results),],
# 		score = this.region.fet.inc.score,
# 		pval = this.region.fet.inc.results$value,
# 		qval = p.adjust(this.region.fet.inc.results$value,'fdr'),
# 		direction = 'increase',
# 		set = 'single',
# 		region = i,
# 		test = 'FET',
# 		stringsAsFactors=FALSE
# 	)
# 
# 	this.region.fet.dec.results = data.frame(
# 		go.terms[rownames(this.region.fet.dec.results),],
# 		score = this.region.fet.dec.score,
# 		pval = this.region.fet.dec.results$value,
# 		qval = p.adjust(this.region.fet.dec.results$value,'fdr'),
# 		direction = 'decrease',
# 		set = 'single',
# 		region = i,
# 		test = 'FET',
# 		stringsAsFactors=FALSE
# 	)
# 	rbind(this.region.fet.inc.results,this.region.fet.dec.results)
# }))
# # },mc.cores=n.cores))
# 
# rownames(single.region.results) = NULL

go.enrichment.results = rbind(all.region.results,one.region.results)
rownames(go.enrichment.results) = NULL

go.enrichment.results$dataset = 'GO'

saveRDS(go.enrichment.results,file='checkpoints/topgo_results.rds')

# n.permutations = 1000
# 
# get.region.specificity = function(enrichment.results,annotation.dataset=NULL,tissue.set='one',effect.direction,statistical.test) {
# 	this.df = subset(enrichment.results,dataset==annotation.dataset & direction == effect.direction & set == tissue.set & test == statistical.test)
# 	this.df = subset(this.df,complete.cases(this.df))
# 	this.wide = tidyr::pivot_wider(droplevels(this.df),id_cols='go_id',names_from='region',values_from='score')
# 	this.mat = matrix(as.matrix(this.wide[,2:ncol(this.wide)]),nrow=nrow(this.wide),ncol=ncol(this.wide)-1,dimnames=list(this.wide$go_id,colnames(this.wide)[2:ncol(this.wide)]))
# # 	apply(this.mat,1,function(x) max(x) - quantile(x,0.75))
# 	specificity = apply(this.mat,1,function(x) sort(x)[length(x)] - sort(x)[length(x)-1])
# 	specificity.null = do.call(cbind,parallel::mclapply(1:n.permutations,function(i) {
# 		apply(matrix(apply(this.mat,2,sample),nrow=nrow(this.mat),dimnames=dimnames(this.mat)),1,function(x) sort(x)[length(x)] - sort(x)[length(x)-1])
# 	},mc.cores=n.cores))
# 	data.frame(
# 		go_id = rownames(this.mat),
# 		dataset = annotation.dataset,
# 		go_name = as.character(unique(subset(this.df,select=c('go_id','go_name')))$go_name),
# 		direction = effect.direction,
# 		test = statistical.test,
# 		top.region = colnames(this.mat)[apply(this.mat,1,which.max)],
# 		specificity,
# 		specificity.pval = rowMeans(specificity < specificity.null),
# 		specificity.qval = p.adjust(rowMeans(specificity < specificity.null),'fdr')
# 	)
# }
# 
# go.fet.inc.spec = get.region.specificity(go.enrichment.results,'GO','one','increase','FET')
# go.fet.dec.spec = get.region.specificity(go.enrichment.results,'GO','one','decrease','FET')
# go.kst.inc.spec = get.region.specificity(go.enrichment.results,'GO','one','increase','KS')
# go.kst.dec.spec = get.region.specificity(go.enrichment.results,'GO','one','decrease','KS')
# 
# go.specificity.results = rbind(
# 	go.fet.inc.spec,
# 	go.fet.dec.spec,
# 	go.kst.inc.spec,
# 	go.kst.dec.spec
# )
# 
# saveRDS(go.specificity.results,file='checkpoints/go_specificity_results.rds')
# 