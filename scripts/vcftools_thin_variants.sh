#!/bin/bash

vcftools --gzvcf vcf/cayo_brain_all_raw.vcf.gz \
	--remove-filtered-all --thin 100000 --remove-indels --maf 0.3 \
	--recode --recode-INFO-all --stdout --max-missing 0.9 > vcf/cayo_brain_all_thinned.vcf
