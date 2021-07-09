#!/usr/bin/env Rscript

source('scripts/_include_options.R')

mash.results = readRDS('checkpoints/mashr_results.rds')

keep.genes = readRDS('checkpoints/keep_genes.rds')

# Import mashr results
library(mashr)
mash.beta = get_pm(mash.results)
mash.lfsr = get_lfsr(mash.results)

# Calculate standardized beta for KS tests
mash.sbet = mash.beta / get_psd(mash.results)

mashr.genes = rownames(mash.beta)
names(mashr.genes) = rownames(mash.beta)

# Analyze enrichment of AD-related genes

# Data downloaded from AD Knowledge Portal
# AMP-AD meta-analysis of AD vs. control
# 7 brain regions:
## Human postmortem brain RNA-seq data were obtained from seven distinct regions: dorsolateral prefrontal cortex (DLPFC), temporal cortex (TCX), inferior frontal gyrus (IFG), superior temporal gyrus (STG), frontal pole (FP), parahippocampal gyrus (PHG), and cerebellum (CBE).
# File: Synapse:syn11914808
# https://doi.org/10.7303/syn11914606
# 
## MSBB
# The results published here are in whole or in part based on data obtained from the AD Knowledge Portal ( https://adknowledgeportal.org/ ). These data were generated from postmortem brain tissue collected through the Mount Sinai VA Medical Center Brain Bank and were provided by Dr. Eric Schadt from Mount Sinai School of Medicine.
## ROSMAP
# The results published here are in whole or in part based on data obtained from the AD Knowledge Portal ( https://adknowledgeportal.org ). Study data were provided by the Rush Alzheimer’s Disease Center, Rush University Medical Center, Chicago. Data collection was supported through funding by NIA grants P30AG10161 (ROS), R01AG15819 (ROSMAP; genomics and RNAseq), R01AG17917 (MAP), R01AG30146, R01AG36042 (5hC methylation, ATACseq), RC2AG036547 (H3K9Ac), R01AG36836 (RNAseq), R01AG48015 (monocyte RNAseq) RF1AG57473 (single nucleus RNAseq), U01AG32984 (genomic and whole exome sequencing), U01AG46152 (ROSMAP AMP-AD, targeted proteomics), U01AG46161(TMT proteomics), U01AG61356 (whole genome sequencing, targeted proteomics, ROSMAP AMP-AD), the Illinois Department of Public Health (ROSMAP), and the Translational Genomics Research Institute (genomic). Additional phenotypic data can be requested at www.radc.rush.edu.
# Cite:
# Genotype data: doi:10.1038/mp.2017.20.
# DNA methylation: doi:10.1038/nn.3786.
# RNAseq: doi:10.1038/s41593 -018-0154-9.
# ChIPseq-H3K9ac: doi:10.1101/273789.
# miRNA: doi:10.1186/s13024 -017-0191-y.
#
## MayoRNAseq
# The results published here are in whole or in part based on data obtained from the AD Knowledge Portal ( https://adknowledgeportal.org ). Study data were provided by the following sources: The Mayo Clinic Alzheimers Disease Genetic Studies, led by Dr. Nilufer Taner and Dr. Steven G. Younkin, Mayo Clinic, Jacksonville, FL using samples from the Mayo Clinic Study of Aging, the Mayo Clinic Alzheimers Disease Research Center, and the Mayo Clinic Brain Bank. Data collection was supported through funding by NIA grants P50 AG016574, R01 AG032990, U01 AG046139, R01 AG018023, U01 AG006576, U01 AG006786, R01 AG025711, R01 AG017216, R01 AG003949, NINDS grant R01 NS080820, CurePSP Foundation, and support from Mayo Foundation. Study data includes samples collected through the Sun Health Research Institute Brain and Body Donation Program of Sun City, Arizona. The Brain and Body Donation Program is supported by the National Institute of Neurological Disorders and Stroke (U24 NS072026 National Brain and Tissue Resource for Parkinsons Disease and Related Disorders), the National Institute on Aging (P30 AG19610 Arizona Alzheimers Disease Core Center), the Arizona Department of Health Services (contract 211002, Arizona Alzheimers Research Center), the Arizona Biomedical Research Commission (contracts 4001, 0011, 05-901 and 1001 to the Arizona Parkinson's Disease Consortium) and the Michael J. Fox Foundation for Parkinsons Research.

ad.genes = read.delim('data/meta.anlz.ad_cntrl.tsv')

# Pull from biomaRt to get ortholog info for AD genes
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

# Only keep genes that have mmul orthologs
hsap.info = unique(subset(hsap.info,nchar(mmulatta_homolog_ensembl_gene) > 0,select=c('ensembl_gene_id','external_gene_name','mmulatta_homolog_ensembl_gene','mmulatta_homolog_orthology_type')))

# Only keep genes that are 1:1 orthologs
hsap.one2one = subset(hsap.info,mmulatta_homolog_orthology_type == 'ortholog_one2one')

# Add mmul ortholog info to AD genes
ad.mmul = merge(ad.genes,subset(hsap.one2one,select=c('ensembl_gene_id','mmulatta_homolog_ensembl_gene')),by='ensembl_gene_id',all.x=TRUE,all.y=FALSE)

# Assign AD genes to categories (filter AD genes by FDR < 0.05) (com = combined, inc = up in AD, dec = down in AD)
ad.com = subset(ad.mmul,fdr.fixed <= 0.05 & abs(TE.fixed) > 0.2)
ad.inc = subset(ad.mmul,fdr.fixed <= 0.05 & abs(TE.fixed) > 0.2 & TE.fixed > 0)
ad.dec = subset(ad.mmul,fdr.fixed <= 0.05 & abs(TE.fixed) > 0.2 & TE.fixed < 0)

# Create list objects with each category of genes (only keeping mmul genes)
ad.deg = list(
	dis = subset(ad.com,mmulatta_homolog_ensembl_gene %in% mashr.genes)$mmulatta_homolog_ensembl_gene,
	inc = subset(ad.inc,mmulatta_homolog_ensembl_gene %in% mashr.genes)$mmulatta_homolog_ensembl_gene,
	dec = subset(ad.dec,mmulatta_homolog_ensembl_gene %in% mashr.genes)$mmulatta_homolog_ensembl_gene
)

# Function that takes a vector of logical values specifying if gene is a DEG in macaque or not (length must equal number of genes in dataset)
fet = function(x,direction,region) {
	# x is a named vector of logical values (i.e., DEG or not)

	# Subset macaque data to only genes that overlap AD analysis
	x = x[names(x) %in% ad.mmul$mmulatta_homolog_ensembl_gene]

	# Standardize the directions
	if (tolower(direction) %in% c('inc','increase','up')) {
		direction = 'inc'
	} else if (tolower(direction) %in% c('dec','decrease','down')) {
		direction = 'dec'
	} else if (tolower(direction) %in% c('var','dis','dispersion')) {
		direction = 'dis'
	} else {
		stop('direction must be "inc", "dec", or "dis".')
	}

	if (!sum(x)) {
		# In the unlikely scenario that there are no DEGs, return NAs
		data.frame(
			region = region,
			direction = direction,
			estimate = NA,
			test = 'FET',
			pval = NA
		)
	} else {

		do.deg.total = as.integer(table(factor(x,levels=c('TRUE','FALSE'))))

		this.deg.total = as.integer(table(factor(ad.deg[[direction]] %in% names(x[x]),levels=c('TRUE','FALSE'))))
		contingency.matrix.deg = matrix(rbind(this.deg.total,do.deg.total - this.deg.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))

		out = fisher.test(contingency.matrix.deg,alternative='greater')
	
		data.frame(
			region = region,
			direction = direction,
			estimate = out$estimate,
			test = 'FET',
			pval = out$p.value
		)
	}
}

kst = function(x,direction,region) {
	# x is a named vector of standardized beta values
	x = x[names(x) %in% ad.mmul$mmulatta_homolog_ensembl_gene]

	if (tolower(direction) %in% c('inc','increase','up')) {
		direction = 'inc'
	} else if (tolower(direction) %in% c('dec','decrease','down')) {
		direction = 'dec'
	} else if (tolower(direction) %in% c('var','dis','dispersion')) {
		direction = 'dis'
	} else {
		stop('direction must be "inc", "dec", or "dis".')
	}

	x1 = x[names(x) %in% ad.deg[[direction]]]
	x2 = x[!names(x) %in% ad.deg[[direction]]]

	out = if (direction == 'inc' || direction == 'dis') {
		ks.test(x1,x2,alternative='less')
	} else if (direction == 'dec') {
		ks.test(x1,x2,alternative='greater')
	}

	data.frame(
		region = region,
		direction = direction,
		estimate = out$statistic,
		test = 'KST',
		pval = out$p.value
	)
}

# Calculate "wbaDEGs" (whole-brain differentially expressed genes)
wbadeg.inc = mashr.genes %in% names(which(unlist(lapply(mashr.genes,function(x) {
	(sum(mash.lfsr[x,] < fsr.cutoff) >= fraction.shared.cutoff * length(keep.genes)) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff] > 0) >= fraction.shared.cutoff * length(keep.genes)
}))))

wbadeg.dec = mashr.genes %in% names(which(unlist(lapply(mashr.genes,function(x) {
	(sum(mash.lfsr[x,] < fsr.cutoff) >= fraction.shared.cutoff * length(keep.genes)) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff] < 0) >= fraction.shared.cutoff * length(keep.genes)
}))))
names(wbadeg.inc) = names(wbadeg.dec) = mashr.genes

# For KS, calculate average standardized beta across all genes
wbadeg.sbet = rowMeans(mash.sbet)

# Put together FET results for whole-brain and single-region analyses.
ad.enrich.fet = rbind(
	fet(wbadeg.inc,'inc','all'),
	fet(wbadeg.dec,'dec','all'),
	do.call(rbind,lapply(region.levels,function(i) {
		out = rbind(
			fet(mash.beta[,i] > 0 & mash.lfsr[,i] < fsr.cutoff,'inc',i),
			fet(mash.beta[,i] < 0 & mash.lfsr[,i] < fsr.cutoff,'dec',i)
		)
		rownames(out) = NULL
		out
	}))
)

# Put together KS results for whole-brain and single-region analyses.
ad.enrich.kst = rbind(
	kst(wbadeg.sbet,'inc','all'),
	kst(wbadeg.sbet,'dec','all'),
	do.call(rbind,lapply(region.levels,function(i) {
		out = rbind(
			kst(mash.sbet[,i],'inc',i),
			kst(mash.sbet[,i],'dec',i)
		)
		rownames(out) = NULL
		out
	}))
)

ad.enrich.fet$qval = p.adjust(ad.enrich.fet$pval,'fdr')
ad.enrich.kst$qval = p.adjust(ad.enrich.kst$pval,'fdr')

rownames(ad.enrich.fet) = rownames(ad.enrich.kst) = NULL

ad.enrich = rbind(ad.enrich.fet,ad.enrich.kst)

## Import DGLM variance results

dglm.results = readRDS('checkpoints/dglm_mashr_results.rds')

dglm.beta = get_pm(dglm.results)
dglm.lfsr = get_lfsr(dglm.results)
dglm.sbet = dglm.beta / get_psd(dglm.results)

dglm.genes = rownames(dglm.beta)
names(dglm.genes) = rownames(dglm.beta)

fraction.shared.cutoff = 1/3

wbadvg.inc = dglm.genes %in% names(which(unlist(lapply(dglm.genes,function(x) {
	(sum(dglm.lfsr[x,] < fsr.cutoff) >= fraction.shared.cutoff * length(keep.genes)) && sum(dglm.beta[x,][dglm.lfsr[x,] < fsr.cutoff] > 0) >= fraction.shared.cutoff * length(keep.genes)
}))))

wbadvg.sbet = rowMeans(dglm.sbet)

ad.dglm.enrich.fet = rbind(
	fet(wbadvg.inc,'dis','all'),
	do.call(rbind,lapply(region.levels,function(i) {
		out = fet(dglm.beta[,i] > 0 & dglm.lfsr[,i] < fsr.cutoff,'dis',i)
		rownames(out) = NULL
		out
	}))
)

ad.dglm.enrich.kst = rbind(
	kst(wbadvg.sbet,'dis','all'),
	do.call(rbind,lapply(region.levels,function(i) {
		out = kst(mash.sbet[,i],'dis',i)
		rownames(out) = NULL
		out
	}))
)

ad.dglm.enrich.fet$qval = p.adjust(ad.dglm.enrich.fet$pval,'fdr')
ad.dglm.enrich.kst$qval = p.adjust(ad.dglm.enrich.kst$pval,'fdr')

rownames(ad.dglm.enrich.fet) = rownames(ad.dglm.enrich.kst) = NULL

ad.dglm.enrich = rbind(ad.dglm.enrich.fet,ad.dglm.enrich.kst)

ad.dglm.enrich$sigFDR_0.05 = ad.dglm.enrich$qval < 0.05
ad.dglm.enrich$sigFDR_0.20 = ad.dglm.enrich$qval < 0.20


saveRDS(ad.enrich,file='checkpoints/alzherimers_results.rds')
saveRDS(ad.dglm.enrich,file='checkpoints/alzherimers_variance_results.rds')