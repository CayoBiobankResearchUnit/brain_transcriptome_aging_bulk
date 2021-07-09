#!/usr/bin/env Rscript

# Number of cores to use for parallel computing (set to 0 to use all detected cores)
n.cores = 0

# Ignore checkpoints (if TRUE, will rerun some code even when a checkpoint file already exists)
ignore.checkpoints = FALSE

# Predictor of interest (set to column name in metadata file)
predictor = 'exact_age_years'

# Predictor of interest (set to preferred label)
predictor.label = 'Age'

# Predictor units of interest (set to preferred label)
predictor.units = 'Years'

# Main batch variable for plotting (set to column name in metadata file)
batch.variable = 'Sequencing.batch'

# Sex is a special case that is necessary for some plots. Set the name of the sex variable here
sex.variable = 'sex'

# Age is another special case that is necessary for some plots. Set the name of the age variable here
age.variable = 'exact_age_years'

# Rank is another special case
rank.variable = 'rank.scaled'
rank.ordinal.variable = 'ordinal.rank'

# Transcripts-per-million cutoff
tpm.cutoff = 10

# False sign rate cutoff (mashr)
fsr.cutoff = 0.2

# False sign rate cutoff (strict)
fsr.strict = 0.01

# Fraction of regions considered shared
fraction.shared.cutoff = 0.85 # 13 or more regions

# Fraction of regions considered unique
fraction.unique.cutoff = 0.5 # 7 or fewer regions

# Full covariates (used for calculating partial residuals)
full.covariates = c('exact_age_years','ordinal.rank','sex','RIN','Sequencing.batch','reads_mapped')

# Batch covariates
batch.covariates = c('RIN','Sequencing.batch','reads_mapped')

# Covariates to include in model(s)
model.covariates = c('exact_age_years','ordinal.rank','sex','RIN','Sequencing.batch')

# Set factor
region.levels = c('dmPFC', 'dlPFC', 'vmPFC', 'vlPFC', 'ACCg', 'M1', 'STS', 'V1', 'AMY', 'CA3', 'DG', 'CN', 'Pu', 'LGN', 'VMH')

region.set1 = c('dmPFC', 'dlPFC', 'vmPFC', 'vlPFC', 'ACCg', 'M1', 'STS', 'V1')
region.set2 = c('AMY', 'CA3', 'DG', 'CN', 'Pu', 'LGN', 'VMH')


# Set colors
region.colors = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666', '#ed1c24', '#00aeef', '#86328c', '#00aaad', '#ac7eaf', '#2e3192', '#6a3e14') # (extension of palette Dark2)

# Set shapes
region.shapes = c(21,24,22,16,17,15,21,24,22,16,17,15,21,24,22)

# Number of facet rows (3 works great for 15 regions)
region.rows = 3

# Set sex colors
sex.colors = c('#4daf4a','#984ea3') # (colors based on palette PRGn)

# Set sex labels (should be in alphabetical order)
sex.labels = c('Female','Male')

# Set rank labels (for ordinal rank only)
rank.ordinal.labels = c('low','mid','high')

# Social metrics
social.metrics = c(
	'percentrank',
	'GroomOUT',
	'GroomIN',
	'Groom.Recip',
	'CSIgroom',
	'isIsolated',
	'numPartnersGroom',
	'between.groom',
	'eig.cent.groom',
	'clusterCoeff.groom',
	'closeness.groom',
	'Kin',
	'propKin',
	'numKin',
	'CSIprox',
	'CSI',
	'AggOUT',
	'AggIN',
	'Agg.Recip',
	'CSIAgg',
	'tenor',
	'vig.ra',
	'sdb.ra',
	'loneliness',
	'social.interest',
	'social.attainment'
)
social.labels = c(
	'Social rank',
	'Active grooming index',
	'Passive grooming index',
	'Grooming reciprocity',
	'Composite grooming index',
	'Isolation',
	'Number of grooming partners',
	'Betweenness centrality',
	'Eigenvector centrality',
	'Clustering coefficient',
	'Grooming closeness',
	'Kin index',
	'Proportion of kin grooming',
	'Number of kin partners',
	'Composite proximity index',
	'Composite sociality index',
	'Active aggression index',
	'Passive aggression index',
	'Aggression reciprocity',
	'Composite aggression index',
	'Tenor',
	'Vigilance',
	'Self-directed behaviors index',
	'Loneliness',
	'Social interest',
	'Social attainment'
)

# Plot ensembl genes (genes to plot against predictor)
plot.ensembl.genes = c('ENSMMUG00000004791','ENSMMUG00000000857','ENSMMUG00000005775','ENSMMUG00000008799','ENSMMUG00000011520','ENSMMUG00000017225','ENSMMUG00000018687','ENSMMUG00000020426','ENSMMUG00000032418','ENSMMUG00000013467','ENSMMUG00000015959','ENSMMUG00000020202','ENSMMUG00000052014','ENSMMUG00000063637')
# Set random seed
seed = 42

# Clustering max predictor filter (set to Inf to include all animals)
max.predictor.range = 20

# Clustering min predictor filter (set to 0 to include all animals)
min.predictor.range = 0

# Number of permutations
n.permutations = 100

# Number of batches for parallelization
n.batches = 20

# Group 1 label
g1.label = 'young'

# Group 2 label
g2.label = 'old'

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
#                             End configurations
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #

# If not set, sets n.cores to detected cores
if (!exists('n.cores') || !n.cores) n.cores = ifelse('future' %in% rownames(installed.packages()),future::availableCores(),parallel::detectCores(logical=FALSE))

get.ks.pval = function (object) 
{
    N <- numAllMembers(object)
    na <- numMembers(object)
    if (na == 0 || na == N) 
        return(1)
    x.a <- rankMembers(object)
    return(ks.test(x.a, seq_len(N)[-x.a], alternative = "greater")$p.value)
}

get.ks.score = function (object) 
{
    N <- numAllMembers(object)
    na <- numMembers(object)
    if (na == 0 || na == N) 
        return(1)
    x.a <- rankMembers(object)
    return(ks.test(x.a, seq_len(N)[-x.a], alternative = "greater")$statistic)
}

get.fisher.pval = function (object) 
{
    contMat <- contTable(object)
    if (all(contMat == 0)) 
        p.value <- 1
    else p.value <- fisher.test(contMat, alternative = "greater")$p.value
    return(p.value)
}

get.fisher.score = function (object) 
{
    contMat <- contTable(object)
    return(as.numeric(fisher.test(contMat, alternative = "greater")$estimate))
}

# Modified CILP function (from https://github.com/AmandaJLea/differential_correlation)
CILP_withLM = function(input_matrix,predictor_DF,predictor_column,covariates_column_start,covariates_column_end,class1_name,class2_name,plot=T) {

	# get all possible pairwise combinations of genes
	genes=dim(input_matrix)[1]
	e_all = as.matrix(input_matrix) # Ensure matrix type
	predictor = predictor_DF[,predictor_column]
	pairs = as.data.frame(t(combn( 1:genes, 2)))

	print('scaling data within the classes to be compared')
	tmp = matrix(nrow=dim(e_all)[1],ncol=dim(e_all)[2])
	for (i in c(1:dim(tmp)[1])) {
		tmp[i,which(predictor==1)] = scale( matrix(e_all[i,which(predictor==1)],ncol=1) )
		tmp[i,which(predictor==0)] = scale( matrix(e_all[i,which(predictor==0)],ncol=1) )
	} 
	e_all_norm = tmp

	print('testing for differences in correlation')

	model = model.matrix(as.formula(paste('~',paste(c(predictor_column,names(predictor_DF)[cov.start:cov.end]),collapse='+'))),predictor_DF)

	results  =  as.data.frame(t(do.call(cbind,lapply(1:dim(pairs)[1], function(i) {
		tmp1 = (e_all_norm[pairs$V1[i],])
		tmp2 = (e_all_norm[pairs$V2[i],])
		return( c( cor.test( e_all_norm[pairs$V1[i],which(predictor==0)] , e_all_norm[pairs$V2[i],which(predictor==0)] ,method='spearman')$estimate, cor.test( e_all_norm[pairs$V1[i],which(predictor==1)] , e_all_norm[pairs$V2[i],which(predictor==1)] ,method='spearman')$estimate, summary(lm( scale(tmp1*tmp2) ~ model))$coefficients[2,4]))
	}))))
	names(results) = c(paste('cor',class1_name,sep='_'),paste('cor',class2_name,sep='_'),'p_value')
	print('saving results')
	write.table(results,'results.txt',row.names=F,quote=F,sep='\t')
	write.table(pairs,'pairs_tested.txt',row.names=F,quote=F,sep='\t')

	if(plot==T) {
		print('plotting results')
		par(mfrow=c(2,2))

		# qq-plot
		qqplot(-log10(runif(dim(pairs)[1])),-log10(results$p_value),pch=20,xlab='p-value, uniform distribution',ylab='p-value, differential correlation test',bty='n',xlim=c(0,max(-log10(results$p_value))),ylim=c(0,max(-log10(results$p_value))))
		x=c(0,1);y=c(0,1);abline(lm(y~x),lty=2)

	} 
	print('done')
}