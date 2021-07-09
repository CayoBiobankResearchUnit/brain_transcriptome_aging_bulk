#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

this = arguments[1]

homer.results = read.delim(paste0('results/homer_age_tf_',this,'/knownResults.txt'))

names(homer.results) = c('motif','consensus','pval','logP','qval','n.target','pct.target','n.background','pct.background')

homer.results$p.value = exp(homer.results$logP)
homer.results$q.value = p.adjust(homer.results$p.value,'fdr')

homer.results$tf.family = with(homer.results,gsub('^(.+?)\\((.+?)\\).+','\\2',motif))
homer.results$tf.name = with(homer.results,gsub('^(.+?)\\((.+?)\\).+','\\1',motif))

homer.results = homer.results[order(homer.results$logP),]

saveRDS(homer.results,file=paste0('checkpoints/homer_results_',this,'.rds'))
