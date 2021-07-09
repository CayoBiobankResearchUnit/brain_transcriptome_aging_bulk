#!/bin/bash

mkdir -p ngsrelate

# Calculate frequencies
zcat angsd/cayo_brain_region_*.mafs.gz | cut -f5 | sed 1d > ngsrelate/cayo_brain_freq.txt

# Make merged glf.gz file
zcat angsd/cayo_brain_region_*.glf.gz | gzip > angsd/cayo_brain_autosomes.glf.gz

tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | cut -f 1 > ngsrelate/cayo_brain_individuals.txt

ngsRelate -g angsd/cayo_brain_autosomes.glf.gz \
	-p ${SLURM_CPUS_ON_NODE} -z ngsrelate/cayo_brain_individuals.txt \
	-n $(wc -l ngsrelate/cayo_brain_individuals.txt | tr -s ' ' | cut -d ' ' -f 2) \
	-f ngsrelate/cayo_brain_freq.txt \
	-O ngsrelate/cayo_brain_relatedness.txt
