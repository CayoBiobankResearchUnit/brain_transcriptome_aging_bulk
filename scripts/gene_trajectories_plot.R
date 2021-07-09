#!/usr/bin/env Rscript

source('scripts/_include_options.R')

# This is a temporary script simply for replotting trajectories without going through
# the expensive computation again

meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')

# Read in results
gene.k.m.trajectories = readRDS('checkpoints/gene_trajectories_kmeans_trajectories.rds')
gene.h.c.trajectories = readRDS('checkpoints/gene_trajectories_hclust_trajectories.rds')

# Filter metadata for predictions
meta = meta[meta[[predictor]] >= min.predictor.range & meta[[predictor]] <= max.predictor.range,]

# Gene trajectories
loess.predictor.span = 1

library(ggplot2)

# Now set a separate predictor range for plotting
predictor.range.plotting = seq(ceiling(min(meta[[predictor]])*2)/2,floor(max(meta[[predictor]])*2)/2,0.1)
# predictor.range.plotting = seq(ceiling(quantile(unique(subset(meta,select=c('Individual',predictor)))[[predictor]],0.025)*2)/2,floor(quantile(unique(subset(meta,select=c('Individual',predictor)))[[predictor]],0.975)*2)/2,0.1)

# Set skip to TRUE to skip repeating the clustering
skip = TRUE
if (!skip) {
	# Set number of iterations
	n = 30

	for (i in 1:n) {
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

		# K-means first
		region.trajectories = gene.k.m.trajectories[[i]]

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

		region.trajectories = gene.h.c.trajectories[[i]]

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
}

# Chosen hclust
chosen.k = 8
region.trajectories = gene.h.c.trajectories[[chosen.k]]

cluster.categories = c('signaling',
					   'metabolic process',
					   'metabolic process',
					   'immune response',
					   'morphogenesis',
					   'signaling',
					   'stress response',
					   'signaling')

#                          01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30
k.rows   = switch(chosen.k, 1,  1,  1,  1,  1,  2,  1,  2,  3,  2,  2,  2,  2,  2,  3,  2,  4,  3,  4,  4,  3,  4,  4,  4,  5,  4,  3,  4,  5,  5)
k.cols   = switch(chosen.k, 1,  2,  3,  4,  5,  3,  7,  4,  3,  5,  6,  6,  7,  7,  5,  8,  5,  6,  5,  5,  7,  6,  6,  6,  5,  7,  9,  7,  6,  6)
k.height = switch(chosen.k, 5,  4,  4,  4,  4,  4,  4,  5,  7,  4,  5,  5,  5,  5,  6,  5,  9,  7,  7,  7,  7,  7,  7,  7,  9,  7,  5,  7,  9,  9)
k.width  = switch(chosen.k, 5,  8,  9, 11, 13,  7, 15,  7,  7,  9, 11, 11, 11, 11,  9, 11,  9, 11,  9,  9, 11, 11, 11, 11,  9, 11, 13, 11, 13, 13)

region.trajectories$cluster.label.category = gsub('(.+?Cluster ([0-9]+).+?)\\(\\"\\*italic\\(n\\)~\\"=\\"~[0-9]+\\*\\"\\)(.+)','\\1<\\2>\\3',region.trajectories$cluster.label.size)
region.trajectories$cluster.label.blank = gsub('"<[0-9]+>"',paste0('phantom(',paste(unique(sort(unlist(strsplit(cluster.categories,'')))),collapse=''),')'),region.trajectories$cluster.label.category)

cluster.labels = unique(region.trajectories[c('cluster.label','cluster.label.blank','cluster.label.size','cluster.label.category')])
cluster.labels$cluster.label.size =  gsub('.+?,atop\\(scriptstyle\\(\\"\\(\\"\\*(.+?[0-9]+)\\*.+','\\1',cluster.labels$cluster.label.size)
cluster.labels$x = min(predictor.range.plotting)
cluster.labels$y = -max(abs(region.trajectories$e)) + max(abs(region.trajectories$e)) * 0.1

cluster.labels$cluster.category = cluster.categories

for (i in 1:chosen.k) {
	region.trajectories$cluster.label.category = gsub(paste0('"<',i,'>"'),paste0('bold("',cluster.categories[i],'")'),region.trajectories$cluster.label.category)
	cluster.labels$cluster.label.category = gsub(paste0('"<',i,'>"'),paste0('bold("',cluster.categories[i],'")'),cluster.labels$cluster.label.category)
}

p = ggplot() +
	geom_blank(data=region.trajectories,aes_string(predictor.label,'e',color='Region')) +
	geom_line(data=region.trajectories,aes_string(predictor.label,'e',color='Region'),size=0.2,alpha=0) +
#	geom_text(data=cluster.labels,aes(x,y,label=cluster.label.size),parse=TRUE, hjust = 0, size=4) +
	facet_wrap(~cluster.label.blank,labeller=label_parsed,nrow=k.rows) +
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
	guides(color = guide_legend(override.aes = list(alpha = 1, size=1))) +
	xlab(paste0(predictor.label,' (',tolower(predictor.units),')')) +
	ylab(expression(italic(Z)*'-score'))
ggsave(p,file=paste0('figures/gene_trajectories_by_clusters_hclust_best_k_',formatC(chosen.k,width=2,flag=0),'_blank.pdf'),height=k.height,width=k.width,useDingbats=FALSE)

p = ggplot() +
#	geom_blank(data=region.trajectories,aes_string(predictor.label,'e',color='Region')) +
	geom_smooth(data=region.trajectories,aes_string(predictor.label,'e'),color='#000000',size=1,method=loess,se=TRUE) +
	geom_line(data=region.trajectories,aes_string(predictor.label,'e',color='Region'),size=0.2) +
	geom_text(data=cluster.labels,aes(x,y,label=cluster.label.size),parse=TRUE, hjust = 0, size=4) +
	facet_wrap(~cluster.label.blank,labeller=label_parsed,nrow=k.rows) +
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
ggsave(p,file=paste0('figures/gene_trajectories_by_clusters_hclust_best_k_',formatC(chosen.k,width=2,flag=0),'_data.pdf'),height=k.height,width=k.width,useDingbats=FALSE)

p = ggplot() +
#	geom_blank(data=region.trajectories,aes_string(predictor.label,'e',color='Region')) +
	geom_smooth(data=region.trajectories,aes_string(predictor.label,'e'),color='#000000',size=1,method=loess,se=TRUE) +
	geom_line(data=region.trajectories,aes_string(predictor.label,'e',color='Region'),size=0.2) +
	geom_text(data=cluster.labels,aes(x,y,label=cluster.label.size),parse=TRUE, hjust = 0, size=4) +
	facet_wrap(~cluster.label.category,labeller=label_parsed,nrow=k.rows) +
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
ggsave(p,file=paste0('figures/gene_trajectories_by_clusters_hclust_best_k_',formatC(chosen.k,width=2,flag=0),'_categories.pdf'),height=k.height,width=k.width,useDingbats=FALSE)

p = ggplot() +
#	geom_blank(data=region.trajectories,aes_string(predictor.label,'e',color='Region')) +
	geom_smooth(data=region.trajectories,aes_string(predictor.label,'e'),color='#000000',size=1,method=loess,se=TRUE) +
	geom_line(data=region.trajectories,aes_string(predictor.label,'e',color='Region'),size=0.2) +
	geom_text(data=cluster.labels,aes(x,y,label=cluster.label.size),parse=TRUE, hjust = 0, size=6) +
	facet_wrap(~cluster.label.category,labeller=label_parsed,nrow=1) +
	scale_color_manual(values=region.colors) +
	scale_x_continuous(
		limits=range(predictor.range.plotting),
		breaks=eval(parse(text=paste0('seq(',paste(round(range(predictor.range.plotting)/5) * 5,collapse=','),',5)')))
	) +
	scale_y_continuous(
		limits=c(-max(abs(range(region.trajectories$e))),max(abs(range(region.trajectories$e)))),
		breaks=seq(round(-max(abs(range(region.trajectories$e)))),round(max(abs(range(region.trajectories$e)))))
	) +
	coord_fixed(ratio = diff(range(predictor.range.plotting)) / (2*max(abs(range(region.trajectories$e)))),clip='off') +
	theme_classic(base_size=24) +
	theme(
		strip.background=element_blank(),
		strip.text=element_text(lineheight=0,margin=margin(b=-10)),
		legend.position='none',
		panel.spacing = unit(2,'lines')
	) +
	guides(color = guide_legend(override.aes = list(size = 1))) +
	xlab(paste0(predictor.label,' (',tolower(predictor.units),')')) +
	ylab(expression(italic(Z)*'-score'))
ggsave(p,file=paste0('figures/gene_trajectories_by_clusters_hclust_best_k_',formatC(chosen.k,width=2,flag=0),'_categories_singlerow.pdf'),height=k.height,width=k.width*2,useDingbats=FALSE)

# gt = ggplot_gtable(ggplot_build(p))
# gt$layout$clip[TRUE] = 'off'
# 
# pdf(file=paste0('figures/gene_trajectories_by_clusters_hclust_best_k_',formatC(chosen.k,width=2,flag=0),'_categories_singlerow.pdf'),height=k.height,width=k.width*2,useDingbats=FALSE)
# 	grid::grid.draw(gt)
# dev.off()

# Pull cluster assignments
cluster.assignments = readRDS('checkpoints/gene_trajectories_hclust_results.rds')[[chosen.k]]

emma.results = readRDS('checkpoints/emma_results.rds')
mash.results = readRDS('checkpoints/mashr_results.rds')

keep.genes = readRDS('checkpoints/keep_genes.rds')

library(mashr)
mash.beta = get_pm(mash.results)
mash.lfsr = get_lfsr(mash.results)
mash.sbet = mash.beta / get_psd(mash.results)

ensembl.gene.names = dimnames(emma.results)[[1]]
mashr.genes = rownames(mash.beta)
names(mashr.genes) = rownames(mash.beta)

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


all.assign = table(cluster.assignments)
inc.beta.assign = table(cluster.assignments[names(all.region.fet[all.region.fet ==  1])])[names(all.assign)]
dec.beta.assign = table(cluster.assignments[names(all.region.fet[all.region.fet == -1])])[names(all.assign)]

names(inc.beta.assign) = names(dec.beta.assign) = names(all.assign)
inc.beta.assign[is.na(inc.beta.assign)] = 0
dec.beta.assign[is.na(dec.beta.assign)] = 0

cluster.info = data.frame(
	cluster.labels,
	total = as.integer(all.assign),
	inc = as.integer(inc.beta.assign),
	dec = as.integer(dec.beta.assign)
)

cluster.info = within(cluster.info,{
	prop.inc = inc / total
	prop.dec = dec / total
	per.inc = paste0(format(round(prop.inc * 100,1),nsmall=1),'%')
	per.dec = paste0(format(round(prop.dec * 100,1),nsmall=1),'%')
	cluster.category.break = gsub(' ','\n',cluster.category)
	cluster.category.break = gsub('morphogenesis','morpho-\ngenesis',cluster.category.break)
})

out = do.call(rbind,lapply(1:chosen.k,function(i) {
	this.prob = cluster.info$total[i]/sum(cluster.info$total)
	inc.outcomes = c(cluster.info$inc[i],sum(cluster.info$inc) - cluster.info$inc[i])
	dec.outcomes = c(cluster.info$dec[i],sum(cluster.info$dec) - cluster.info$dec[i])
	data.frame(
		cluster.label = cluster.info$cluster.label[i],
		inc.p = binom.test(inc.outcomes,p=this.prob,alternative='greater')$p.value,
		dec.p = binom.test(dec.outcomes,p=this.prob,alternative='greater')$p.value
	)
}))

out = within(out,{
	inc.q = p.adjust(inc.p,'fdr')
	dec.q = p.adjust(dec.p,'fdr')
})

cluster.info = merge(cluster.info,out,by='cluster.label')

p1 = ggplot() +
#	geom_blank(data=region.trajectories,aes_string(predictor.label,'e',color='Region')) +
	geom_smooth(data=region.trajectories,aes_string(predictor.label,'e'),color='#000000',size=1,method=loess,se=FALSE) +
	geom_line(data=region.trajectories,aes_string(predictor.label,'e',color='Region'),size=0.2) +
	facet_wrap(~cluster.label,nrow=1) +
	scale_color_manual(values=region.colors) +
	scale_x_continuous(
		limits=range(predictor.range.plotting)
	) +
	scale_y_continuous(
		limits=c(-max(abs(range(region.trajectories$e))),max(abs(range(region.trajectories$e))))
	) +
	coord_fixed(ratio = diff(range(predictor.range.plotting)) / (2*max(abs(range(region.trajectories$e))))) +
	theme_classic(base_size=16) +
	theme(
		strip.background=element_blank(),
		strip.text=element_text(size=14),
#		strip.text=element_blank(),
		axis.text=element_blank(),
		axis.title=element_blank(),
		axis.ticks=element_blank(),
		axis.line.y=element_blank(),
		legend.position='none'
	)


p2 = ggplot() +
#	geom_blank(data=region.trajectories,aes_string(predictor.label,'e',color='Region')) +
	geom_text(data=cluster.info,aes(1,1,label=cluster.category.break),parse=FALSE, hjust = 0.5, size=4, fontface='bold') +
	facet_wrap(~cluster.label,nrow=1) +
	theme_classic(base_size=16) +
	theme(
		strip.background=element_blank(),
		strip.text=element_blank(),
		axis.text=element_blank(),
		axis.title=element_blank(),
		axis.ticks=element_blank(),
		axis.line=element_blank(),
		legend.position='none'
	)

p3 = ggplot() +
#	geom_blank(data=region.trajectories,aes_string(predictor.label,'e',color='Region')) +
	geom_text(data=cluster.info,aes(1,1,label=cluster.label.size),parse=TRUE, hjust = 0.5, size=4) +
	facet_wrap(~cluster.label,nrow=1) +
	theme_classic(base_size=16) +
	theme(
		strip.background=element_blank(),
		strip.text=element_blank(),
		axis.text=element_blank(),
		axis.title=element_blank(),
		axis.ticks=element_blank(),
		axis.line=element_blank(),
		legend.position='none'
	)


library(ggtext)

cluster.long1 = reshape2::melt(cluster.info[,c('cluster.label','prop.inc','prop.dec')],id.vars='cluster.label')
cluster.long2 = reshape2::melt(cluster.info[,c('cluster.label','per.inc','per.dec')],id.vars='cluster.label')
cluster.long3 = reshape2::melt(cluster.info[,c('cluster.label','inc','dec')],id.vars='cluster.label')
cluster.long4 = reshape2::melt(cluster.info[,c('cluster.label','inc.q','dec.q')],id.vars='cluster.label')
names(cluster.long1) = c('cluster.label','variable','prop')
names(cluster.long2) = c('cluster.label','variable','per')
names(cluster.long3) = c('cluster.label','variable','dir')
names(cluster.long4) = c('cluster.label','variable','q')

cluster.long = data.frame(cluster.long1,cluster.long2[3],cluster.long3[3],cluster.long4[3])

cluster.long$dir.label = as.character(cluster.long$dir)
cluster.long$dir.label[cluster.long$dir.label == '0'] = ''

cluster.long$significant = ''
cluster.long$significant[cluster.long$q < fsr.cutoff] = '***'

cluster.long$variable = factor(cluster.long$variable,levels=c('prop.dec','prop.inc'))
# cluster.long$variable.emoji = cluster.long$variable
# levels(cluster.long$variable.emoji) = c('⬇️','⬆️')
# cluster.long$variable.emoji = as.character(cluster.long$variable.emoji)
# 
# 
# # https://www.hvitfeldt.me/blog/real-emojis-in-ggplot2/
# emoji_to_link <- function(x) {
#   paste0("https://emojipedia.org/emoji/",x) %>%
#     read_html() %>%
#     html_nodes("tr td a") %>%
#     .[1] %>%
#     html_attr("href") %>%
#     paste0("https://emojipedia.org/", .) %>%
#     read_html() %>%
#     html_node('div[class="vendor-image"] img') %>%
#     html_attr("src")
# }
# 
# link_to_img <- function(x, size = 25) {
#   paste0("<img src='", x, "' width='", size, "'/>")
# }
# 
# library(tidyverse)
# library(emo)
# library(xml2)
# library(rvest)
# 
# cluster.long = cluster.long %>%
# 	mutate(emoji=ji_extract_all(variable.emoji)) %>%
# 	unnest(cols = c(emoji)) %>%
# 	mutate(url = map_chr(emoji, slowly(~emoji_to_link(.x), rate_delay(1))), label = link_to_img(url))
# 
# cluster.long$label = factor(cluster.long$label,levels=unique(cluster.long$label)[order(!grepl('down',unique(cluster.long$label)))])
# levels(cluster.long$label) = gsub("width='25'","width='16'",levels(cluster.long$label))

cluster.long$label = cluster.long$variable

levels(cluster.long$label) = c(
	"<img src='https://emojipedia-us.s3.dualstack.us-west-1.amazonaws.com/thumbs/320/openmoji/252/down-arrow_2b07.png' width='16'/>",
	"<img src='http://emojipedia-us.s3.dualstack.us-west-1.amazonaws.com/thumbs/320/openmoji/252/up-arrow_2b06.png' width='16'/>"
#	"<img src='https://cdn.jsdelivr.net/gh/hfg-gmuend/openmoji/color/svg/2B07.svg' width='16'/>",
#	"<img src='https://cdn.jsdelivr.net/gh/hfg-gmuend/openmoji/color/svg/2B06.svg' width='16'/>"
)

p4 = ggplot(data=cluster.long,aes(label,prop,fill=variable)) +
	geom_bar(stat='identity') +
	geom_text(aes(label=significant),size=6, vjust=-0.1) +
	facet_wrap(~cluster.label,nrow=1) +
	scale_x_discrete() +
	scale_y_continuous(limits=c(0,1.2),breaks=seq(0,1,0.5),labels=c('0%','','100%')) +
	scale_fill_manual(values=c('#f1a340','#998ec3')) +
	theme_classic(base_size=16) +
	theme(
		strip.background=element_blank(),
		strip.text=element_blank(),
#		axis.text=element_blank(),
#		axis.text.x=element_text(angle=-60,hjust=0,vjust=1,size=12),
		axis.text.x=element_markdown(size=3),
		axis.text.y=element_text(size=12),
#		axis.text.y=element_blank(),
		axis.title=element_blank(),
#		axis.ticks.y=element_blank(),
		legend.position='none'
	)

library(egg)

pdf(file='figures/gene_trajectories_by_clusters_hclust_best_k_mashr_genes_overlap.pdf',width=11,height=5)
ggarrange(p1,p2,p3,p4,ncol=1,nrow=4,widths=15,heights=c(1,0.75,0.375,2),newpage=FALSE)
dev.off()
