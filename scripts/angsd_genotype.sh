#!/bin/bash

region=`sed -n ${SLURM_ARRAY_TASK_ID}p genomes/angsd_regions.txt`

mkdir -p angsd

angsd -GL 2 -out angsd/cayo_brain_region_$(printf "%08d\n" $SLURM_ARRAY_TASK_ID) \
	-nThreads ${SLURM_CPUS_ON_NODE} \
	-doGlf 3 -doMajorMinor 1 -doMaf 1 -minmaf 0.05 -SNP_pval 1e-6 \
	-bam data/merged_genotypes.bamlist.txt \
	-r ${region}
