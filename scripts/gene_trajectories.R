#!/usr/bin/env Rscript

source('scripts/_include_options.R')

e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')
keep.genes = readRDS('checkpoints/keep_genes.rds')
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')
e.regressed = readRDS('checkpoints/regressed_expression_matrix.rds')

# Read in annotations
mmul.go = readRDS('checkpoints/rnaseq_genes_go.rds')
mmul.do = readRDS('checkpoints/rnaseq_genes_do.rds')
mmul.dg = readRDS('checkpoints/rnaseq_genes_dg.rds')
go.terms = readRDS('checkpoints/rnaseq_go_names.rds')

gene2go = lapply(
	unique(mmul.go$ensembl_gene_id),
	function(x) {
		out = sort(mmul.go[mmul.go$ensembl_gene_id == x,'go_id'])
		out[out != '']
	}
)
names(gene2go) = unique(mmul.go$ensembl_gene_id)

ensembl.gene.names = rownames(e.regressed)

library(parallel)
library(doParallel)

# Filter metadata for predictions
meta = meta[meta[[predictor]] >= min.predictor.range & meta[[predictor]] <= max.predictor.range,]

# Gene trajectories
loess.predictor.span = 1

library(ggplot2)

animal.libraries = with(meta,split(Library,Individual))
region.libraries = with(meta,split(Library,Region))
animal.meta = unique(subset(meta,select=c('Individual',predictor)))
rownames(animal.meta) = animal.meta$Individual

# Z-score expression matrix (do this within each library)
# Within each region, extracts the expression matrix for libraries and calculates z-score
e.z = do.call(cbind,lapply(region.libraries,function(libraries) {
	t(apply(e.regressed[,libraries],1,function(x) {
		(x - mean(x)) / sd(x)
	}))
}))

# Global (median) expression matrix
# Calculate median expression across all libraries from an animal
e.m = do.call(cbind,mclapply(names(animal.libraries),function(x) {
	matrix(apply(e.z[,animal.libraries[[x]]],1,median),ncol=1,dimnames=list(rownames(e.z),x))
},mc.cores=n.cores))

# Set the predictor range (age range) for clustering
predictor.range.clustering = seq(ceiling(min(meta[[predictor]])*2)/2,floor(max(meta[[predictor]])*2)/2,0.5)
# predictor.range.clustering = seq(ceiling(quantile(unique(subset(meta,select=c('Individual',predictor)))[[predictor]],0.025)*2)/2,floor(quantile(unique(subset(meta,select=c('Individual',predictor)))[[predictor]],0.975)*2)/2,0.5)

# Now set a separate predictor range for plotting
predictor.range.plotting = seq(ceiling(min(meta[[predictor]])*2)/2,floor(max(meta[[predictor]])*2)/2,0.1)
# predictor.range.plotting = seq(ceiling(quantile(unique(subset(meta,select=c('Individual',predictor)))[[predictor]],0.025)*2)/2,floor(quantile(unique(subset(meta,select=c('Individual',predictor)))[[predictor]],0.975)*2)/2,0.1)

# Expression by time point (smoothed with loess)
e.t = do.call(rbind,mclapply(rownames(e.m),function(x) {
	this.df = data.frame(e = e.m[x,],a = animal.meta[colnames(e.m),][[predictor]])
	out = predict(loess(e~a,data=this.df),predictor.range.clustering,span=loess.predictor.span)
	names(out) = paste0('T',formatC(predictor.range.clustering,digits=1,format='f'))
	out
},mc.cores=n.cores))
rownames(e.t) = rownames(e.m)

# Initialize topGO object

all.clusters = numeric(length=length(ensembl.gene.names))
names(all.clusters) = ensembl.gene.names

# Initialize topGO data

library(topGO)

all.clusters.topgo = new('topGOdata',
	description='Simple session',
	ontology='BP',
	allGenes=all.clusters,
	geneSelectionFun=function(x) x == 0,
	nodeSize = 10,
	annotationFun = annFUN.gene2GO,
	gene2GO = gene2go
)

# k-means clustering
library(cluster)
library(clValid)

# Set number of iterations
n = 50

# Initialize results vectors
w.ss = numeric(n)
w.ss[1] = (nrow(e.t) - 1) * sum(apply(e.t,2,var))

avg.silhouette.values = numeric(n - 1)

gene.k.m.results = vector('list',n)
gene.k.m.go.results = gene.k.m.do.results = gene.k.m.dg.results = vector('list',n)

gene.h.c.results = vector('list',n)
gene.h.c.go.results = gene.h.c.do.results = gene.h.c.dg.results = vector('list',n)

# Set up hclust clustering

euclidean.distance = dist(e.t,method='euclidean')
h.clust = hclust(euclidean.distance,method='complete')

# Dunn index (https://www.datacamp.com/community/tutorials/hierarchical-clustering-R)
dunn.values = numeric(n-1)

gene.k.m.trajectories = gene.h.c.trajectories = vector('list',n)

for (i in 2:n) {
	message('Now running k = ',i)
	
	# Set number of rows for faceting
#                       01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30
	k.rows   = switch(i, 1,  1,  1,  1,  1,  2,  1,  2,  3,  2,  2,  2,  2,  2,  3,  2,  4,  3,  4,  4,  3,  4,  4,  4,  5,  4,  3,  4,  5,  5)
	k.cols   = switch(i, 1,  2,  3,  4,  5,  3,  7,  4,  3,  5,  6,  6,  7,  7,  5,  8,  5,  6,  5,  5,  7,  6,  6,  6,  5,  7,  9,  7,  6,  6)
	k.height = switch(i, 5,  4,  4,  4,  4,  4,  4,  5,  7,  4,  5,  5,  5,  5,  6,  5,  9,  7,  7,  7,  7,  7,  7,  7,  9,  7,  5,  7,  9,  9)
	k.width  = switch(i, 5,  8,  9, 11, 13,  7, 15,  7,  7,  9, 11, 11, 11, 11,  9, 11,  9, 11,  9,  9, 11, 11, 11, 11,  9, 11, 13, 11, 13, 13)

	if (i > 30) {
		k.rows = k.cols = NULL
		k.height = k.width = ceiling(i/3)
	}

	km.result = gene.k.m.results[[i]] = kmeans(e.t,centers=i,nstart=25)

	hc.result = gene.h.c.results[[i]] = cutree(h.clust,k=i)

	# Save elbow scores (within-cluster sum of squares)
	w.ss[i] = sum(km.result$withinss)

	# Save silhouette score
	avg.silhouette.values[i-1] = mean(silhouette(km.result$cluster, dist(e.t))[,3])

	km.clusters = data.frame(gene=rownames(e.t),cluster=hc.result,stringsAsFactors=FALSE)
	hc.clusters = data.frame(gene=rownames(e.t),cluster=cutree(h.clust,k=i),stringsAsFactors=FALSE)

	dunn.values[i-1] = dunn(euclidean.distance,hc.clusters$cluster)

	these.clusters = all.clusters
	these.clusters[km.clusters$gene] = km.clusters$cluster

	these.do = merge(km.clusters,mmul.do,by.x='gene',by.y='ensembl_gene_id')
	these.dg = merge(km.clusters,mmul.dg,by.x='gene',by.y='ensembl_gene_id')

	these.do = subset(these.do,do_id %in% names(which(table(subset(these.do,confidence >= 0)$do_id) >= 10)))
	these.dg = subset(these.dg,dg_id %in% names(which(table(subset(these.dg,dg_score >= 0 & grepl('F03',dg_class))$dg_id) >= 10)))

	gene.k.m.go.results[[i]] = do.call(rbind,mclapply(split(km.clusters, km.clusters$cluster),function(x) {
		this.cluster = unique(x$cluster)
		this.topgo = all.clusters.topgo
		this.topgo@allScores = as.integer(these.clusters[this.topgo@allGenes])
		this.topgo@geneSelectionFun=eval(parse(text=paste0('function(y) y==',this.cluster)))
		fet = runTest(this.topgo,algorithm='parentchild',statistic='fisher')
		data.frame(
			k=i,
			cluster=this.cluster,
			go_id=names(fet@score),
			pval=fet@score,
			qval=p.adjust(fet@score,'fdr'))
	},mc.cores=n.cores))

	gene.k.m.do.results[[i]] = do.call(rbind,mclapply(split(km.clusters, km.clusters$cluster),function(x) {
		this.cluster = unique(x$cluster)

		do.total = as.integer(table(factor(these.do$cluster == this.cluster,levels=c('TRUE','FALSE'))))

		out = do.call(rbind,lapply(split(these.do,these.do$do_id),function(y) {
			this.total = as.integer(table(factor(y$cluster == this.cluster,levels=c('TRUE','FALSE'))))
			contingency.matrix = matrix(rbind(this.total,do.total - this.total),nrow=2,dimnames=list(c('in group','not in group'),c('in cluster','not in cluster')))
			fet.test = fisher.test(contingency.matrix,alternative='greater')
			data.frame(
				k=i,
				cluster=this.cluster,
				do_id=unique(y$do_id),
				pval=fet.test$p.value
			)
		}))
		out$qval = p.adjust(out$pval,'fdr')
		rownames(out) = NULL
		out
	},mc.cores=n.cores))

	gene.k.m.dg.results[[i]] = do.call(rbind,mclapply(split(km.clusters, km.clusters$cluster),function(x) {
		this.cluster = unique(x$cluster)

		dg.total = as.integer(table(factor(these.dg$cluster == this.cluster,levels=c('TRUE','FALSE'))))

		out = do.call(rbind,lapply(split(these.dg,these.dg$dg_id),function(y) {
			this.total = as.integer(table(factor(y$cluster == this.cluster,levels=c('TRUE','FALSE'))))
			contingency.matrix = matrix(rbind(this.total,dg.total - this.total),nrow=2,dimnames=list(c('in group','not in group'),c('in cluster','not in cluster')))
			fet.test = fisher.test(contingency.matrix,alternative='greater')
			data.frame(
				k=i,
				cluster=this.cluster,
				dg_id=unique(y$dg_id),
				pval=fet.test$p.value
			)
		}))
		out$qval = p.adjust(out$pval,'fdr')
		rownames(out) = NULL
		out
	},mc.cores=n.cores))

	these.clusters = all.clusters
	these.clusters[hc.clusters$gene] = hc.clusters$cluster

	gene.h.c.go.results[[i]] = do.call(rbind,mclapply(split(hc.clusters, hc.clusters$cluster),function(x) {
		this.cluster = unique(x$cluster)
		this.topgo = all.clusters.topgo
		this.topgo@allScores = as.integer(these.clusters[this.topgo@allGenes])
		this.topgo@geneSelectionFun=eval(parse(text=paste0('function(y) y==',this.cluster)))
		fet = runTest(this.topgo,algorithm='parentchild',statistic='fisher')
		data.frame(
			k=i,
			cluster=this.cluster,
			go_id=names(fet@score),
			pval=fet@score,
			qval=p.adjust(fet@score,'fdr'))
	},mc.cores=n.cores))

	gene.h.c.do.results[[i]] = do.call(rbind,mclapply(split(hc.clusters, hc.clusters$cluster),function(x) {
		this.cluster = unique(x$cluster)

		do.total = as.integer(table(factor(these.do$cluster == this.cluster,levels=c('TRUE','FALSE'))))

		out = do.call(rbind,lapply(split(these.do,these.do$do_id),function(y) {
			this.total = as.integer(table(factor(y$cluster == this.cluster,levels=c('TRUE','FALSE'))))
			contingency.matrix = matrix(rbind(this.total,do.total - this.total),nrow=2,dimnames=list(c('in group','not in group'),c('in cluster','not in cluster')))
			fet.test = fisher.test(contingency.matrix,alternative='greater')
			data.frame(
				k=i,
				cluster=this.cluster,
				do_id=unique(y$do_id),
				pval=fet.test$p.value
			)
		}))
		out$qval = p.adjust(out$pval,'fdr')
		rownames(out) = NULL
		out
	},mc.cores=n.cores))

	gene.h.c.dg.results[[i]] = do.call(rbind,mclapply(split(hc.clusters, hc.clusters$cluster),function(x) {
		this.cluster = unique(x$cluster)

		dg.total = as.integer(table(factor(these.dg$cluster == this.cluster,levels=c('TRUE','FALSE'))))

		out = do.call(rbind,lapply(split(these.dg,these.dg$dg_id),function(y) {
			this.total = as.integer(table(factor(y$cluster == this.cluster,levels=c('TRUE','FALSE'))))
			contingency.matrix = matrix(rbind(this.total,dg.total - this.total),nrow=2,dimnames=list(c('in group','not in group'),c('in cluster','not in cluster')))
			fet.test = fisher.test(contingency.matrix,alternative='greater')
			data.frame(
				k=i,
				cluster=this.cluster,
				dg_id=unique(y$dg_id),
				pval=fet.test$p.value
			)
		}))
		out$qval = p.adjust(out$pval,'fdr')
		rownames(out) = NULL
		out
	},mc.cores=n.cores))

#	gene.clusters.long = do.call(rbind,lapply(split(km.clusters,km.clusters$cluster),function(x) {
#		out = reshape2::melt(e.t[x$gene,])
#		out$Var2 = as.numeric(gsub('^T','',out$Var2))
#		out$cluster = unique(x$cluster)
#		names(out) = c('gene',predictor.label,'e','cluster')
#		out
#	}))
#	gene.clusters.long = gene.clusters.long[with(gene.clusters.long,order(gene,Age)),]

	# If i is 2, plot k=1 as well
	if (i == 2) {
		region.trajectories = do.call(rbind,lapply(split(km.clusters,1),function(x) {
			this.cluster = 1
			e.by.region = do.call(rbind,mclapply(names(region.libraries),function(y) {
				this.df = data.frame(
					e = colMeans(e.z[rownames(subset(km.clusters,TRUE)),region.libraries[[y]]]),
					a = meta[region.libraries[[y]],][[predictor]],
					stringsAsFactors=FALSE
				)
				out = predict(loess(e~a,data=this.df),predictor.range.plotting)
				names(out) = paste0('T',formatC(predictor.range.plotting,digits=1,format='f'))
				out
			},mc.cores=n.cores))
			rownames(e.by.region) = names(region.libraries)
			output = reshape2::melt(e.by.region)
			names(output) = c('Region',predictor.label,'e')
			output$cluster = this.cluster
			output[[predictor.label]] = as.numeric(gsub('^T','',output[[predictor.label]]))
			output$Region = factor(output$Region,levels=region.levels)
			output
		}))
		region.trajectories$cluster.label = region.trajectories$cluster.label.size = factor(paste('Cluster',region.trajectories$cluster),levels=paste('Cluster',unique(sort(region.trajectories$cluster))))
		region.trajectories$cluster.size = as.integer(nrow(km.clusters))[region.trajectories$cluster]
		levels(region.trajectories$cluster.label.size) = with(region.trajectories,paste0('atop(textstyle("',levels(cluster.label),'\"),atop(scriptstyle("("*italic(n)~"="~',as.integer(nrow(km.clusters)),'*")")))'))
	#	levels(region.trajectories$cluster.label.size) = with(region.trajectories,paste0('textstyle("',levels(cluster.label),'\")~scriptstyle(\"("*italic(n)~"="~',as.integer(table(km.clusters$cluster)),'*")")'))

		region.trajectories = subset(region.trajectories,complete.cases(region.trajectories))

		gene.k.m.trajectories[[1]] = region.trajectories
		gene.h.c.trajectories[[1]] = region.trajectories

		p = ggplot() +
			geom_line(data=region.trajectories,aes_string(predictor.label,'e',color='Region'),size=0.2) +
			geom_smooth(data=region.trajectories,aes_string(predictor.label,'e'),color='#000000',size=1,method=loess,se=TRUE) +
			facet_wrap(~cluster.label.size,labeller=label_parsed,nrow=k.rows) +
			scale_color_manual(values=region.colors) +
			scale_x_continuous(
				limits=range(predictor.range.plotting),
				breaks=eval(parse(text=paste0('seq(',paste(round(range(predictor.range.plotting)/5) * 5,collapse=','),',5)')))
			) +
			scale_y_continuous(
				limits=c(-max(abs(range(region.trajectories$e))),max(abs(range(region.trajectories$e)))),
				breaks=seq(round(-max(abs(range(region.trajectories$e)))),round(max(abs(range(region.trajectories$e)))))
			) +
			coord_fixed(ratio = diff(range(predictor.range.plotting)) / (2*max(abs(range(region.trajectories$e))))) +
			theme_classic(base_size=16) +
			theme(
				strip.background=element_blank(),
				strip.text=element_text(lineheight=0,margin=margin(b=-10))
			) +
			guides(color = guide_legend(override.aes = list(size = 1))) +
			xlab(paste0(predictor.label,' (',tolower(predictor.units),')')) +
			ylab(expression(italic(Z)*'-score'))
		ggsave(p,file=paste0('figures/gene_trajectories_by_clusters_kmeans_k_',formatC(1,width=2,flag=0),'.pdf'),useDingbats=FALSE)
		ggsave(p,file=paste0('figures/gene_trajectories_by_clusters_hclust_k_',formatC(1,width=2,flag=0),'.pdf'),useDingbats=FALSE)
	}

	# Trajectories from k-means
	region.trajectories = do.call(rbind,lapply(split(km.clusters,km.clusters$cluster),function(x) {
		this.cluster = unique(x$cluster)
		e.by.region = do.call(rbind,mclapply(names(region.libraries),function(y) {
			this.df = data.frame(
				e = colMeans(matrix(e.z[rownames(subset(km.clusters,cluster == this.cluster)),region.libraries[[y]]],nrow=sum(km.clusters$cluster == this.cluster))),
				a = meta[region.libraries[[y]],][[predictor]],
				stringsAsFactors=FALSE
			)
			out = predict(loess(e~a,data=this.df),predictor.range.plotting)
			names(out) = paste0('T',formatC(predictor.range.plotting,digits=1,format='f'))
			out
		},mc.cores=n.cores))
		rownames(e.by.region) = names(region.libraries)
		output = reshape2::melt(e.by.region)
		names(output) = c('Region',predictor.label,'e')
		output$cluster = this.cluster
		output[[predictor.label]] = as.numeric(gsub('^T','',output[[predictor.label]]))
		output$Region = factor(output$Region,levels=region.levels)
		output
	}))
	region.trajectories$cluster.label = region.trajectories$cluster.label.size = factor(paste('Cluster',region.trajectories$cluster),levels=paste('Cluster',unique(sort(region.trajectories$cluster))))
	region.trajectories$cluster.size = as.integer(table(km.clusters$cluster))[region.trajectories$cluster]
	levels(region.trajectories$cluster.label.size) = with(region.trajectories,paste0('atop(textstyle("',levels(cluster.label),'\"),atop(scriptstyle("("*italic(n)~"="~',as.integer(table(km.clusters$cluster)),'*")")))'))
#	levels(region.trajectories$cluster.label.size) = with(region.trajectories,paste0('textstyle("',levels(cluster.label),'\")~scriptstyle(\"("*italic(n)~"="~',as.integer(table(km.clusters$cluster)),'*")")'))

	region.trajectories = subset(region.trajectories,complete.cases(region.trajectories))

	gene.k.m.trajectories[[i]] = region.trajectories

	p = ggplot() +
		geom_line(data=region.trajectories,aes_string(predictor.label,'e',color='Region'),size=0.2) +
		geom_smooth(data=region.trajectories,aes_string(predictor.label,'e'),color='#000000',size=1,method=loess,se=TRUE) +
		facet_wrap(~cluster.label.size,labeller=label_parsed,nrow=k.rows) +
		scale_color_manual(values=region.colors) +
		scale_x_continuous(
			limits=range(predictor.range.plotting),
			breaks=eval(parse(text=paste0('seq(',paste(round(range(predictor.range.plotting)/5) * 5,collapse=','),',5)')))
		) +
		scale_y_continuous(
			limits=c(-max(abs(range(region.trajectories$e))),max(abs(range(region.trajectories$e)))),
			breaks=seq(round(-max(abs(range(region.trajectories$e)))),round(max(abs(range(region.trajectories$e)))))
		) +
		coord_fixed(ratio = diff(range(predictor.range.plotting)) / (2*max(abs(range(region.trajectories$e))))) +
		theme_classic(base_size=16) +
		theme(
			strip.background=element_blank(),
			strip.text=element_text(lineheight=0,margin=margin(b=-10))
		) +
		guides(color = guide_legend(override.aes = list(size = 1))) +
		xlab(paste0(predictor.label,' (',tolower(predictor.units),')')) +
		ylab(expression(italic(Z)*'-score'))
	ggsave(p,file=paste0('figures/gene_trajectories_by_clusters_kmeans_k_',formatC(i,width=2,flag=0),'.pdf'),height=k.height,width=k.width,useDingbats=FALSE)


	# Trajectories from hclust
	region.trajectories = do.call(rbind,lapply(split(hc.clusters,hc.clusters$cluster),function(x) {
		this.cluster = unique(x$cluster)
		e.by.region = do.call(rbind,mclapply(names(region.libraries),function(y) {
			this.df = data.frame(
				e = colMeans(matrix(e.z[rownames(subset(hc.clusters,cluster == this.cluster)),region.libraries[[y]]],nrow=sum(hc.clusters$cluster == this.cluster))),
				a = meta[region.libraries[[y]],][[predictor]],
				stringsAsFactors=FALSE
			)
			out = predict(loess(e~a,data=this.df),predictor.range.plotting)
			names(out) = paste0('T',formatC(predictor.range.plotting,digits=1,format='f'))
			out
		},mc.cores=n.cores))
		rownames(e.by.region) = names(region.libraries)
		output = reshape2::melt(e.by.region)
		names(output) = c('Region',predictor.label,'e')
		output$cluster = this.cluster
		output[[predictor.label]] = as.numeric(gsub('^T','',output[[predictor.label]]))
		output$Region = factor(output$Region,levels=region.levels)
		output
	}))
	region.trajectories$cluster.label = region.trajectories$cluster.label.size = factor(paste('Cluster',region.trajectories$cluster),levels=paste('Cluster',unique(sort(region.trajectories$cluster))))
	region.trajectories$cluster.size = as.integer(table(hc.clusters$cluster))[region.trajectories$cluster]
	levels(region.trajectories$cluster.label.size) = with(region.trajectories,paste0('atop(textstyle("',levels(cluster.label),'\"),atop(scriptstyle("("*italic(n)~"="~',as.integer(table(hc.clusters$cluster)),'*")")))'))
#	levels(region.trajectories$cluster.label.size) = with(region.trajectories,paste0('textstyle("',levels(cluster.label),'\")~scriptstyle(\"("*italic(n)~"="~',as.integer(table(hc.clusters$cluster)),'*")")'))

	region.trajectories = subset(region.trajectories,complete.cases(region.trajectories))

	gene.h.c.trajectories[[i]] = region.trajectories

	p = ggplot() +
		geom_line(data=region.trajectories,aes_string(predictor.label,'e',color='Region'),size=0.2) +
		geom_smooth(data=region.trajectories,aes_string(predictor.label,'e'),color='#000000',size=1,method=loess,se=TRUE) +
		facet_wrap(~cluster.label.size,labeller=label_parsed,nrow=k.rows) +
		scale_color_manual(values=region.colors) +
		scale_x_continuous(
			limits=range(predictor.range.plotting),
			breaks=eval(parse(text=paste0('seq(',paste(round(range(predictor.range.plotting)/5) * 5,collapse=','),',5)')))
		) +
		scale_y_continuous(
			limits=c(-max(abs(range(region.trajectories$e))),max(abs(range(region.trajectories$e)))),
			breaks=seq(round(-max(abs(range(region.trajectories$e)))),round(max(abs(range(region.trajectories$e)))))
		) +
		coord_fixed(ratio = diff(range(predictor.range.plotting)) / (2*max(abs(range(region.trajectories$e))))) +
		theme_classic(base_size=16) +
		theme(
			strip.background=element_blank(),
			strip.text=element_text(lineheight=0,margin=margin(b=-10))
		) +
		guides(color = guide_legend(override.aes = list(size = 1))) +
		xlab(paste0(predictor.label,' (',tolower(predictor.units),')')) +
		ylab(expression(italic(Z)*'-score'))
	ggsave(p,file=paste0('figures/gene_trajectories_by_clusters_hclust_k_',formatC(i,width=2,flag=0),'.pdf'),height=k.height,width=k.width,useDingbats=FALSE)

}

best.silhouette = (2:n)[which.max(avg.silhouette.values)]
best.elbow = which.min(unlist(lapply(1:n,function(x) {
	d = diff(w.ss[(x-1):(x+1)])
	d[2]/d[1]
})))
best.dunn = (2:n)[which.max(dunn.values)]


pdf(file='figures/gene_trajectories_kmeans_genes_elbow.pdf',useDingbats=FALSE)
	plot(1:n, w.ss, type='b', pch = 19, frame = FALSE, xlab='Number of clusters', ylab='Within groups sum of squares',xlim=c(1,n),xaxt = 'n')
	points(w.ss[best.elbow] ~ best.elbow, col='red', cex=2)
	axis(side = 1, at=1:n)
dev.off()

pdf(file='figures/gene_trajectories_kmeans_genes_silhouette.pdf',useDingbats=FALSE)
	plot(2:n, avg.silhouette.values, type = 'b', pch = 19, frame = FALSE,  xlab = 'Number of clusters', ylab = 'Average silhoettes', xlim=c(1,n),xaxt = 'n')
	points(avg.silhouette.values[best.silhouette-1] ~ best.silhouette, col='red', cex=2)
	axis(side = 1, at=1:n)
dev.off()

pdf(file='figures/gene_trajectories_hclust_genes_dunn.pdf',useDingbats=FALSE)
	plot(2:n, dunn.values, type = 'b', pch = 19, frame = FALSE,  xlab = 'Number of clusters', ylab = 'Dunn index', xlim=c(1,n),xaxt = 'n')
	points(dunn.values[best.dunn-1] ~ best.dunn, col='red', cex=2)
	axis(side = 1, at=1:n)
dev.off()

gene.k.m.go.results = do.call(rbind,gene.k.m.go.results)
gene.k.m.do.results = do.call(rbind,gene.k.m.do.results)
gene.k.m.dg.results = do.call(rbind,gene.k.m.dg.results)
gene.h.c.go.results = do.call(rbind,gene.h.c.go.results)
gene.h.c.do.results = do.call(rbind,gene.h.c.do.results)
gene.h.c.dg.results = do.call(rbind,gene.h.c.dg.results)

saveRDS(gene.k.m.results,file='checkpoints/gene_trajectories_kmeans_results.rds')
saveRDS(gene.k.m.go.results,file='checkpoints/gene_trajectories_kmeans_go_results.rds')
saveRDS(gene.k.m.do.results,file='checkpoints/gene_trajectories_kmeans_do_results.rds')
saveRDS(gene.k.m.dg.results,file='checkpoints/gene_trajectories_kmeans_dg_results.rds')
saveRDS(gene.h.c.results,file='checkpoints/gene_trajectories_hclust_results.rds')
saveRDS(gene.h.c.go.results,file='checkpoints/gene_trajectories_hclust_go_results.rds')
saveRDS(gene.h.c.do.results,file='checkpoints/gene_trajectories_hclust_do_results.rds')
saveRDS(gene.h.c.dg.results,file='checkpoints/gene_trajectories_hclust_dg_results.rds')

saveRDS(gene.k.m.trajectories,file='checkpoints/gene_trajectories_kmeans_trajectories.rds')
saveRDS(gene.h.c.trajectories,file='checkpoints/gene_trajectories_hclust_trajectories.rds')
