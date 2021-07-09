#!/bin/bash

i=$(tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | cut -f 1 | sed -n ${SLURM_ARRAY_TASK_ID}p)

mkdir -p bam_lists

# Make list of bams of this genotype
tail -n+2  data/cayo_brain_bulk_metadata_technical.tsv | \
	awk '$2 == "'$i'"' | cut -f 1 | \
	sed 's/\(.*\)/star\/\1.bam/g' > \
	bam_lists/${i}.bamlist.txt

mkdir -p bam

# If the merged, sorted, and reheadered bam already exists, skip these steps
if [ $(samtools quickcheck bam/${i}.bam && echo 0 || echo 1) -eq 1 ]; then

	# Merge bams from the same genotype into one file
	samtools merge -@ ${SLURM_CPUS_ON_NODE} -rfb bam_lists/${i}.bamlist.txt bam/${i}.merged.bam;

	# Sort bam
	samtools sort -@ ${SLURM_CPUS_ON_NODE} bam/${i}.merged.bam -o bam/${i}.sorted.bam;

	# Add @RG info to other header metadata
	mkdir -p headers
	cat <(samtools view -H bam/${i}.sorted.bam | grep '@HD') \
		<(samtools view -H bam/${i}.sorted.bam | grep '@SQ') \
		<(cat bam_lists/${i}.bamlist.txt | sed 's:bam/\([A-Z0-9_]*\)\.bam:@RG\tID\:\1\tPL\:Illumina\tSM\:'${i}'\tLB\:\1:g') \
		<(samtools view -H bam/${i}.sorted.bam | grep '@PG') \
		<(samtools view -H bam/${i}.sorted.bam | grep '@CO') \
		> headers/${i}.rg.txt

	# Replace SAM header with updated @RG information
	samtools reheader headers/${i}.rg.txt bam/${i}.sorted.bam > bam/${i}.reheadered.bam

fi

# If the bam index already exists, skip this step
if [ ! -f bam/${i}.bam.bai ]; then
	samtools index -@ ${SLURM_CPUS_ON_NODE} -b bam/${i}.reheadered.bam
fi
