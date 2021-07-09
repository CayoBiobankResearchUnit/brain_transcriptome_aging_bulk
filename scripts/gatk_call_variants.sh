#!/bin/bash

chr=$(echo {1..20} | cut -d ' ' -f $SLURM_ARRAY_TASK_ID)

genome=genomes/Mmul_10.dna.fa

mkdir -p vcf

gatk --java-options '-Xmx64g' HaplotypeCaller \
	--reference ${genome} \
	$(ls bam/*.recalibrated.chr${chr}.bam | sed 's/^/--input /g') \
	--standard-min-confidence-threshold-for-calling 20.0 \
	--output vcf/cayo_brain_chr${chr}_raw.vcf \
	--intervals ${chr}
