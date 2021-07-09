#!/bin/bash

i=$(tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | cut -f 1 | sed -n $((($SLURM_ARRAY_TASK_ID-1)/20 + 1))p)
chr=$(( ($SLURM_ARRAY_TASK_ID-1) % 20 + 1))

genome=genomes/Mmul_10.dna.fa

mkdir -p reports

gatk --java-options '-Xmx64g' MarkDuplicates \
	--INPUT bam/${i}.chr${chr}.bam \
	--OUTPUT bam/${i}.dedupped.chr${chr}.bam  \
	--CREATE_INDEX true \
	--VALIDATION_STRINGENCY SILENT \
	--METRICS_FILE reports/${i}.dedupped.chr${chr}.metrics.txt
	
gatk --java-options '-Xmx64g' SplitNCigarReads \
	--reference ${genome} \
	--input bam/${i}.dedupped.chr${chr}.bam \
	--output bam/${i}.split.chr${chr}.bam

# Reference file data/cayo.pass.vcf required, indexed with the following:
# gatk IndexFeatureFile --feature-file data/cayo.pass.vcf

gatk --java-options '-Xms16g' BaseRecalibrator \
	--reference ${genome} \
	--input bam/${i}.split.chr${chr}.bam \
	--use-original-qualities \
	--output reports/${i}.chr${chr}.recalibration.report.txt \
	--known-sites data/cayo.pass.vcf

gatk --java-options "-Xms16g" ApplyBQSR \
   --add-output-sam-program-record \
   --reference ${genome} \
   --input bam/${i}.split.chr${chr}.bam \
   --use-original-qualities \
   --output bam/${i}.recalibrated.chr${chr}.bam \
   --bqsr-recal-file reports/${i}.chr${chr}.recalibration.report.txt
