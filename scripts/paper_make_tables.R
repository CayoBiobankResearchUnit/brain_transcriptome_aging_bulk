#!/usr/bin/env Rscript

# GO enrichment

go.traj = subset(readRDS('checkpoints_backup/gene_trajectories_hclust_go_results.rds'),k == 8)

library(GO.db)

go.traj$go_name = Term(go.traj$go_id)

go.traj = go.traj[order(go.traj$cluster,go.traj$pval),]


write.table(subset(do.call(rbind,lapply(split(go.traj,go.traj$cluster),head,n=20)),select=c('cluster','go_id','go_name','pval','qval')),file='tables/gene_trajectories_go.txt',sep='\t',row.names=FALSE,quote=FALSE)

go.enrichment.results = readRDS('checkpoints/topgo_results.rds')

go.enrichment.results$region = factor(go.enrichment.results$region,levels=region.levels)
levels(go.enrichment.results$region) = gsub('ACCg','ACC',levels(go.enrichment.results$region))

go.enrichment.results = go.enrichment.results[with(go.enrichment.results,order(
	region,
	test,
	direction,
	pval
)),]

# Increase in expression in all regions
subset(go.enrichment.results,direction == 'increase' & set =='union' & test == 'FET' & qval < 0.05,select=c('go_id','go_name','direction','pval','qval'))

# Decrease in expression in all regions
subset(go.enrichment.results,direction == 'decrease' & set =='union' & test == 'FET' & qval < 0.05)

write.table(subset(go.enrichment.results,set =='union' & test == 'FET' & qval < 0.05,select=c('go_id','go_name','direction','pval','qval')),file='tables/wbadegs_go.txt',sep='\t',row.names=FALSE,quote=FALSE)
write.table(subset(go.enrichment.results,set =='one' & test == 'FET' & qval < 0.05,select=c('region','go_id','go_name','direction','pval','qval')),file='tables/adegs_go.txt',sep='\t',row.names=FALSE,quote=FALSE)


go.specificity.results = readRDS('checkpoints/go_specificity_results.rds')
go.specificity.results = go.specificity.results[order(go.specificity.results$specificity.pval),]

disease.enrichment.results = readRDS('checkpoints/disease_enrichment_results.rds')

disease.enrichment.results$region = factor(disease.enrichment.results$region,levels=region.levels)
levels(disease.enrichment.results$region) = gsub('ACCg','ACC',levels(disease.enrichment.results$region))

disease.enrichment.results = disease.enrichment.results[with(disease.enrichment.results,order(
	region,
	dataset,
	test,
	direction,
	pval
)),]

write.table(subset(disease.enrichment.results,dataset == 'DISEASES' & set == 'union' & test == 'FET' & qval < 0.05,select=c('do_id','do_name','direction','pval','qval')),file='tables/wbadegs_diseases.txt',sep='\t',row.names=FALSE,quote=FALSE)
write.table(subset(disease.enrichment.results,dataset == 'DISEASES' & set == 'one' & test == 'FET' & qval < 0.05,select=c('region','do_id','do_name','direction','pval','qval')),file='tables/adegs_diseases.txt',sep='\t',row.names=FALSE,quote=FALSE)

disease.specificity.results = readRDS('checkpoints/disease_specificity_results.rds')

mash.results = readRDS('checkpoints/mashr_results.rds')

library(mashr)

mash.beta = get_pm(mash.results)
mash.lfsr = get_lfsr(mash.results)

sig.genes = names(c(
	which(unlist(lapply(rownames(mash.lfsr),function(x) { out = sum(mash.lfsr[x,] < 0.005 & mash.beta[x,] > 0) >= fraction.shared.cutoff * 15; names(out) = x; out }))),
	which(unlist(lapply(rownames(mash.lfsr),function(x) { out = sum(mash.lfsr[x,] < 0.005 & mash.beta[x,] < 0) >= fraction.shared.cutoff * 15; names(out) = x; out })))
))

wbadeg = names(c(
	which(unlist(lapply(rownames(mash.lfsr),function(x) { out = sum(mash.lfsr[x,] < fsr.cutoff & mash.beta[x,] > 0) >= fraction.shared.cutoff * 15; names(out) = x; out }))),
	which(unlist(lapply(rownames(mash.lfsr),function(x) { out = sum(mash.lfsr[x,] < fsr.cutoff & mash.beta[x,] < 0) >= fraction.shared.cutoff * 15; names(out) = x; out })))
))


library(biomaRt)

mmul = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL',
		dataset='mmulatta_gene_ensembl')

mmul.sig = getBM(attributes=c('ensembl_gene_id','external_gene_name','hsapiens_homolog_ensembl_gene'),filters='ensembl_gene_id',values=names(sig.genes),mart=mmul)
mmul.wbadeg = getBM(attributes=c('ensembl_gene_id','external_gene_name','hsapiens_homolog_ensembl_gene'),filters='ensembl_gene_id',values=wbadeg,mart=mmul)

wbadeg = unique(mmul.wbadeg[c('ensembl_gene_id','external_gene_name')])
rownames(wbadeg) = wbadeg$ensembl_gene_id

wbadeg.lfsr = mash.lfsr[rownames(wbadeg),]
# Order table by maximum LFSR across top 13 regions
sort.mash = rownames(wbadeg.lfsr)[order(apply(wbadeg.lfsr,1,function(x) {
	max(sort(x)[1:13])
}))]

wbadeg = wbadeg[sort.mash,]

lfsr.join = as.data.frame(mash.lfsr[rownames(wbadeg),])
beta.join = as.data.frame(mash.beta[rownames(wbadeg),])

names(lfsr.join) = paste0('LFSR ',names(lfsr.join))
names(beta.join) = paste0('beta ',names(beta.join))

wbadeg = data.frame(wbadeg,lfsr.join,beta.join)

write.table(wbadeg,file='tables/wbadegs.txt',sep='\t',row.names=FALSE,quote=FALSE)


genewalk.results = readRDS('checkpoints/genewalk_results.rds')
genewalk.sig = subset(genewalk.results,ensembl_id %in% mmul.sig$hsapiens_homolog_ensembl_gene)

genewalk.sig = genewalk.sig[order(genewalk.sig$hgnc_symbol,-genewalk.sig$global_padj),]

# R biomaRt down, but here are the results from web search
# Gene stable ID	Gene name
# ENSMMUG00000000857	ANXA4
# ENSMMUG00000005775	PLAAT3
# ENSMMUG00000008799	MAMU-F
# ENSMMUG00000011520	FKBP5
# ENSMMUG00000013467	AMACR
# ENSMMUG00000015959	ARPP21
# ENSMMUG00000017225	C20H16orf89
# ENSMMUG00000018687	CD99
# ENSMMUG00000020202	
# ENSMMUG00000020426	SERPINA1
# ENSMMUG00000032418	
# ENSMMUG00000052014	TMEM38A
# ENSMMUG00000063637	

# FKBP5: very interesting
# SERPINA1: interesting
# AMACR: maybe interesting
# ANXA4: maybe interesting
# CD99: maybe interesting
# ARPP21: interesting (neurodegenerative disease)
# TMEM38A: maybe interesting
do.data = read.table('data/human_disease_associations.tsv',
	sep='\t',
	quote='',
	col.names=c('protein_id','protein_name','do_id','do_name','z_score','confidence'),
	stringsAsFactors=FALSE)

dg.data = read.table(
	'data/curated_gene_disease_associations.tsv',
	sep='\t',
	quote='',
	header = TRUE,
	stringsAsFactors=FALSE)

subset(do.data,protein_name %in% foo)
subset(dg.data,geneSymbol %in% foo & grepl('F03',diseaseClass))

homer.1 = readRDS('checkpoints/homer_results_1.rds')
homer.2 = readRDS('checkpoints/homer_results_2.rds')
table(subset(homer.2,qval < 0.01)$tf.family) / sum(table(subset(homer.2,qval < 0.01)$tf.family))
(table(subset(homer.2,qval >= 0.01)$tf.family) / sum(table(subset(homer.2,qval >= 0.01)$tf.family)))[names(table(subset(homer.2,qval < 0.01)$tf.family) / sum(table(subset(homer.2,qval < 0.01)$tf.family)))]

# Split tf.family

homer.out = do.call(rbind,lapply(1:nrow(homer.2),function(i) {
	this = homer.2[i,]
	out = this[,names(homer.2)[names(homer.2) != 'tf.family']]
	out.families = toupper(with(this,unlist(strsplit(gsub('\\?','',tf.family),'[+,:]'))))
	if (length(out.families) == 1 & !nchar(out.families[1])) out.families = 'UNKNOWN'
	out.families = out.families[as.logical(nchar(out.families))]
	out = out[rep(1,length(out.families)),]
	out$tf.family = out.families
	out
}))

tf.family.enrichment = do.call(rbind,lapply(unique(homer.out$tf.family),function(i) {
	sig1.group1 = length(unique(subset(homer.out,tf.family == i & qval <  0.01)$motif))
	sig0.group1 = length(unique(subset(homer.out,tf.family == i & qval >= 0.01)$motif))
	sig1.group0 = length(unique(subset(homer.out,tf.family != i & qval <  0.01)$motif))
	sig0.group0 = length(unique(subset(homer.out,tf.family != i & qval >= 0.01)$motif))
	contingency.matrix = matrix(c(sig1.group1,sig1.group0,sig0.group1,sig0.group0),nrow=2,dimnames=list(c('in group','not in group'),c('sig','not sig')))
	test = fisher.test(contingency.matrix,alternative='greater')
	this.tf = gsub(paste0('.*?(',i,').*'),'\\1',names(sort(table(subset(homer.2,grepl(i,toupper(tf.family)))$tf.family),decreasing=TRUE))[1])
	data.frame(tf.family=this.tf,pval = test$p.value,sig1.group1,sig0.group1,sig1.group0,sig0.group0)
}))

tf.family.enrichment = tf.family.enrichment[order(tf.family.enrichment$pval),]
subset(tf.family.enrichment,sig1.group1 + sig0.group1 >= 10)

subset(homer.2,qval < 0.01)

write.table(subset(homer.2,qval < 0.01,select=c('tf.name','tf.family','motif','consensus','logP','qval')),file='tables/homer_results.txt',sep='\t',row.names=FALSE,quote=FALSE)



sample.cell.proportions = readRDS(paste0('../cayo_sc_brain/checkpoints/dlpfc_cell_proportions.rds'))

sample.cluster.proportions = readRDS(paste0('../cayo_sc_brain/checkpoints/dlpfc_cluster_proportions.rds'))

lapply(split(sample.cell.proportions,sample.cell.proportions$type),function(x) coef(summary(lm(prop~predictor+group,data=x)))['predictor',])

lapply(split(sample.cluster.proportions,sample.cluster.proportions$cluster_type),function(x) coef(summary(lm(prop~predictor+group,data=x)))['predictor',])

cluster.top.markers = readRDS(paste0('../cayo_sc_brain/checkpoints/dlpfc_clusters_top_markers.rds'))


mash.dglm = readRDS('checkpoints/dglm_mashr_results.rds')

m.dglm.beta = get_pm(mash.dglm)
m.dglm.lfsr = get_lfsr(mash.dglm)


wba.var = names(c(
	which(unlist(lapply(rownames(m.dglm.lfsr),function(x) { out = sum(m.dglm.lfsr[x,] < fsr.cutoff & m.dglm.beta[x,] > 0) >= 1/3 * 15; names(out) = x; out })))
))

mmul.var = getBM(attributes=c('ensembl_gene_id','external_gene_name'),filters='ensembl_gene_id',values=wba.var,mart=mmul)

wbadvg = unique(mmul.var[c('ensembl_gene_id','external_gene_name')])
rownames(wbadvg) = wbadvg$ensembl_gene_id

wbadvg.lfsr = m.dglm.lfsr[rownames(wbadvg),]
# Order table by maximum LFSR across top 13 regions
sort.mash = rownames(wbadvg.lfsr)[order(apply(wbadvg.lfsr,1,function(x) {
	max(sort(x)[1:5])
}))]

wbadvg = wbadvg[sort.mash,]

lfsr.join = as.data.frame(m.dglm.lfsr[rownames(wbadvg),])
beta.join = as.data.frame(m.dglm.beta[rownames(wbadvg),])

names(lfsr.join) = paste0('LFSR ',names(lfsr.join))
names(beta.join) = paste0('beta ',names(beta.join))

wbadvg = data.frame(wbadvg,lfsr.join,beta.join)

write.table(wbadvg,file='tables/wbadvgs.txt',sep='\t',row.names=FALSE,quote=FALSE)

dglm.topgo = readRDS('checkpoints/dglm_topgo_results.rds')

write.table(subset(dglm.topgo,method == 'MASH' & qval < 0.05 & set == 'union' & direction == 'increase',select=c('go_id','go_name','direction','pval','qval')),file='tables/wbadvgs_go.txt',sep='\t',row.names=FALSE,quote=FALSE)



sc.mash = readRDS('../cayo_sc_brain/checkpoints/dlpfc_model_emma_mashr.rds')
bk.mash = readRDS('../cayo_bulk_brain/checkpoints/mashr_results.rds')

sc.keep.genes = readRDS(paste0('../cayo_sc_brain/checkpoints/',dataset,'_keep_genes.rds'))


sc.lfsr = get_lfsr(sc.mash)
bk.lfsr = get_lfsr(bk.mash)

sc.beta = get_pm(sc.mash)
bk.beta = get_pm(bk.mash)

# Count of genes sig in any cell type
table(apply(sc.lfsr,1,function(x) any(x < 0.2)))

# Overlap between sig genes and bulk genes: 390 genes overlap
table(names(which(apply(sc.lfsr,1,function(x) any(x < 0.2)))) %in% rownames(bk.beta))

# 119 not identified by our bulk analysis (any tissue)
table(intersect(names(which(apply(sc.lfsr,1,function(x) any(x < 0.2)))),rownames(bk.beta)) %in% names(which(apply(bk.lfsr,1,function(x) any(x < 0.2)))))

# 239 not identified by our bulk analysis (dlPFC)
table(intersect(names(which(apply(sc.lfsr,1,function(x) any(x < 0.2)))),rownames(bk.beta)) %in% names(which(bk.lfsr[,'dlPFC'] < 0.2)))





library(mashr)

sc.mash.beta = get_pm(sc.mash)
sc.mash.lfsr = get_lfsr(sc.mash)

scdeg = names(c(
	which(unlist(lapply(rownames(sc.mash.lfsr),function(x) { out = sum(sc.mash.lfsr[x,] < fsr.cutoff & sc.mash.beta[x,] > 0) >= 1; names(out) = x; out }))),
	which(unlist(lapply(rownames(sc.mash.lfsr),function(x) { out = sum(sc.mash.lfsr[x,] < fsr.cutoff & sc.mash.beta[x,] < 0) >= 1; names(out) = x; out })))
))


library(biomaRt)

mmul = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL',
		dataset='mmulatta_gene_ensembl')

mmul.scdeg = getBM(attributes=c('ensembl_gene_id','external_gene_name','hsapiens_homolog_ensembl_gene'),filters='ensembl_gene_id',values=scdeg,mart=mmul)

scdeg = unique(mmul.scdeg[c('ensembl_gene_id','external_gene_name')])
rownames(scdeg) = scdeg$ensembl_gene_id

scdeg.lfsr = sc.mash.lfsr[rownames(scdeg),]
# Order table by maximum LFSR across top 13 regions
sort.mash = rownames(scdeg.lfsr)[order(apply(scdeg.lfsr,1,function(x) {
	max(sort(x)[1])
}))]

scdeg = scdeg[sort.mash,]

lfsr.join = as.data.frame(sc.mash.lfsr[rownames(scdeg),])
beta.join = as.data.frame(sc.mash.beta[rownames(scdeg),])

names(lfsr.join) = paste0('LFSR ',names(lfsr.join))
names(beta.join) = paste0('beta ',names(beta.join))

scdeg = data.frame(scdeg,lfsr.join,beta.join)

write.table(scdeg,file='tables/scdegs.txt',sep='\t',row.names=FALSE,quote=FALSE)













# Top markers
top.cell = readRDS('../cayo_sc_brain/checkpoints/dlpfc_celltypes_top_markers.rds')
top.clusters = readRDS('../cayo_sc_brain/checkpoints/dlpfc_clusters_top_markers.rds')


top.cell = top.cell %>% filter(fraction_expressing >= 0.1) %>% group_by(cell_group) %>% top_n(20, pseudo_R2)
top.clusters = top.clusters %>% filter(fraction_expressing >= 0.1) %>% group_by(cell_group) %>% top_n(10, pseudo_R2)


write.table(top.cell,file='tables/sc_celltype_markergenes.txt',sep='\t',row.names=FALSE,quote=FALSE)
write.table(top.clusters,file='tables/sc_cellclusters_markergenes.txt',sep='\t',row.names=FALSE,quote=FALSE)











sc.topgo = readRDS('../cayo_sc_brain/checkpoints/dlpfc_topgo_results_emma_mashr.rds')


subset(sc.topgo,pval < 0.001,select=c('go_id','go_name','pval','qval','direction','cell'))



dlpfc.meta = SummarizedExperiment::colData(readRDS('../cayo_sc_brain/checkpoints/dlpfc_classified_manual.rds'))
macmul.integrated.hits = readRDS('../cayo_sc_brain/checkpoints/dlpfc_singleR_integrated_results_matches.rds')
macmul.integrated.hits = macmul.integrated.hits[match(levels(dlpfc.meta$assigned_cell_cluster),macmul.integrated.hits$assigned_cell_cluster),]

write.table(subset(macmul.integrated.hits,select=c('assigned_cell_cluster','assigned_cell_type','cell_type_accession_label','cell_type_alias_label','score','class_label','subclass_label','class2','class.diff','subclass2','subclass.diff')),
	file='tables/singleR_assignments_mouse.txt',sep='\t',row.names=FALSE,quote=FALSE)

macmul.integrated.hits = readRDS('../cayo_sc_brain/checkpoints/dlpfc_homsap_singleR_integrated_results_matches.rds')

write.table(subset(macmul.integrated.hits,select=c('assigned_cell_cluster','assigned_cell_type','cell_type_accession_label','cell_type_alias_label','score','class_label','subclass_label','class2','class.diff','subclass2','subclass.diff')),
	file='tables/singleR_assignments_human.txt',sep='\t',row.names=FALSE,quote=FALSE)


# 

# pseuobulk.emma.sig = subset(pseudobulk.emma,qval.age < 0.05)
# 
# ensembl.gene.info = read.delim('data/mmul_gene_names.txt')
# names(ensembl.gene.info) = c('gene','gene_name')
# ensembl.gene.info = subset(ensembl.gene.info,nchar(gene_name) > 0)
# 
# pseuobulk.emma.sig = merge(pseuobulk.emma.sig,ensembl.gene.info,by='gene',all.x=TRUE)
# pseuobulk.emma.sig = pseuobulk.emma.sig[with(pseuobulk.emma.sig,order(cell,pval.age)),]
# 
# emma.topgo = readRDS('../cayo_sc_brain/checkpoints/dlpfc_topgo_results_emma.rds')
# emma.topgo = emma.topgo[with(emma.topgo,order(cell,direction,pval)),]
# 
# subset(emma.topgo,pval < 0.005,select=c('go_name','pval','qval','direction','cell'))


# Predictions performance
mash.predictions = readRDS(paste0('checkpoints/mashr_model_predictions_cv_k_36_merged.rds'))
glmnet.predictions = readRDS(paste0('checkpoints/glmnet_model_predictions_cv_k_36_merged.rds'))

checkpoint.prefixes = list('mashr_model_predictions_cv_k_36_merged',
'mashr_model_predictions_social_0.01_cv_k_36_merged',
'mashr_model_predictions_social_0.001_cv_k_36_merged',
'mashr_model_predictions_social_0.0001_cv_k_36_merged',
'glmnet_model_predictions_cv_k_36_merged',
'glmnet_model_predictions_social_0.2_cv_k_36_merged',
'glmnet_model_predictions_social_0.1_cv_k_36_merged',
'glmnet_model_predictions_social_0.05_cv_k_36_merged',
'glmnet_model_predictions_social_0.01_cv_k_36_merged',
'glmnet_model_predictions_social_0.001_cv_k_36_merged',
'glmnet_model_predictions_social_0.0001_cv_k_36_merged')

out = do.call(rbind,lapply(1:length(checkpoint.prefixes),function(i) {
	this = readRDS(paste0('checkpoints/',checkpoint.prefixes[[i]],'.rds'))
	with(this,
		data.frame(
			mean.error = mean(abs(delta)),
			mean.resid = mean(abs(residual)),
			rmse = sqrt(mean((predicted - known)^2)),
			slope = as.numeric(coef(lm(predicted~known))['known']),
			intercept = as.numeric(coef(lm(predicted~known))['(Intercept)']),
			r.squared = summary(lm(predicted~known))$adj.r.squared,
			cor.est = cor.test(known,predicted,method='pearson')$estimate,
			cor.p = cor.test(known,predicted,method='pearson')$p.value
		)
	)
}))


# Social comparisons

all.data = do.call(rbind,lapply(region.levels,function(i) {
	get.matching.data(main.comparisons,region=i,social.var=social.var,aging.var=aging.var,sig.cutoff=1,keep.criterion='both')
}))

out = data.frame(
	region=names(split(all.data,all.data$region)),
	p.cor = unlist(lapply(split(all.data,all.data$region),function(x) cor.test(x$aging.effect,x$social.effect,alternative='greater')$estimate)),
	p.p.value = unlist(lapply(split(all.data,all.data$region),function(x) cor.test(x$aging.effect,x$social.effect,alternative='greater')$p.value)),
	cor = unlist(lapply(split(all.data,all.data$region),function(x) cor.test(x$aging.effect,x$social.effect,alternative='greater',method='spearman')$estimate)),
	p.value = unlist(lapply(split(all.data,all.data$region),function(x) cor.test(x$aging.effect,x$social.effect,alternative='greater',method='spearman')$p.value))
)

# 


kegg.test.results = readRDS(paste0('checkpoints/kegg_tests_rank_age.rds'))

# Age predictions
mash.df = readRDS('checkpoints/mashr_model_predictions_cv_k_36_merged.rds')

glmnet.df = readRDS('checkpoints/glmnet_model_predictions_cv_k_36_merged.rds')



dlpfc.meta = readRDS('../cayo_sc_brain/checkpoints/dlpfc_metadata.rds')

write.table(
	subset(dlpfc.meta,select=c('sample','animal_id','sex','region','group','age','extraction_batch')),
	file='tables/sc_meta.txt',
	sep='\t',row.names=FALSE,quote=FALSE)
	
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')

# full.covariates = c('exact_age_years','ordinal.rank','sex','RIN','Sequencing.batch','reads_mapped')
meta$group = 'HH'
write.table(subset(meta,select=c('Library','Individual','sex','Region','group','exact_age_years','ordinal.rank','Sequencing.batch','RIN','reads_mapped')),
	file='tables/bk_meta.txt',
	sep='\t',row.names=FALSE,quote=FALSE)
	
	
	
	
# Results from social filtering predictions

lm.results = readRDS('checkpoints/comparison_model_predictions_lm_tests.rds')
pw.results = readRDS('checkpoints/comparison_model_predictions_pw_tests.rds')

write.table(lm.results,
	file='tables/social_lm_results.txt',
	sep='\t',row.names=FALSE,quote=FALSE)
write.table(pw.results,
	file='tables/social_pw_results.txt',
	sep='\t',row.names=FALSE,quote=FALSE)




