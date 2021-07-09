#!/bin/bash

mkdir -p lcmlkin

lcmlkin -i vcf/cayo_brain_all_thinned.vcf \
	-o lcmlkin/cayo_brain_all_thinned.txt \
	-l phred -g all -t ${SLURM_CPUS_ON_NODE}

# Copy results to data/ folder

cp lcmlkin/cayo_brain_all_thinned.txt data/cayo_brain_lcmlkin_results.txt
