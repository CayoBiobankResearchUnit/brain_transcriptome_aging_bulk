#!/usr/bin/env Rscript

# Read in fasta index (first 20 rows are autosomes)
fai = read.table('genomes/Mmul_10.dna.fa.fai',sep='\t')[1:20,]

regions = unlist(lapply(1:20,function(x) {
	# Calculate lower bounds of the regions (1e6 intervals)
	region.lower = seq(1,fai$V2[x],1e6)
	# Calculate upper bounds of the regions
	region.upper = c(region.lower[2:length(region.lower)] - 1,fai$V2[x])
	gsub(' ','',paste0(x,':',format(region.lower,scientific=FALSE),'-',format(region.upper,scientific=FALSE)))
}))

write(regions,file='genomes/angsd_regions.txt',sep='\t')
