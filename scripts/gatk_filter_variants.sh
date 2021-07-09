#!/bin/bash

chr=$(echo {1..20} | cut -d ' ' -f $SLURM_ARRAY_TASK_ID)

genome=genomes/Mmul_10.dna.fa

gatk --java-options '-Xmx16g' VariantFiltration \
	--reference ${genome} \
	--output vcf/cayo_brain_chr${chr}_filtered.vcf \
	--variant vcf/cayo_brain_chr${chr}_raw.vcf \
	--window 35 \
	--cluster 3 \
	--filter-name 'FS' \
	--filter 'FS > 30.0' \
	--filter-name 'QD' \
	--filter 'QD < 2.0'

# gatk VariantFiltration \
# 	--reference ${genome} \
# 	--output vcf/cayo_brain_chr${chr}_filtered.vcf \
# 	--variant vcf/cayo_brain_chr${chr}_raw.vcf \
# 	--filter-name "QDfilter" \
# 	--filter "QD < 2.0" \
# 	--filter-name "MQfilter" \
# 	--filter "MQ < 40.0" \
# 	--filter-name "FSfilter" \
# 	--filter "FS > 60.0" \
# 	--filter-name "HAPSCfilter" \
# 	--filter "HaplotypeScore > 13.0" \
# 	--filter-name "MQRSfilter" \
# 	--filter "MQRankSum < -12.5" \
# 	--filter-name "RPRSfilter" \
# 	--filter "ReadPosRankSum < -8.0" \
# 	--missing-values-evaluate-as-failing	
# 