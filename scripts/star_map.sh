#!/bin/bash

lib=$(tail -n+2 data/cayo_brain_bulk_metadata_technical.tsv | cut -f 1 | sed -n ${SLURM_ARRAY_TASK_ID}p)

mkdir -p star

STAR --runThreadN $SLURM_CPUS_ON_NODE \
    --genomeDir genomes \
    --readFilesIn fastq/${lib}.R1.fastq.gz,fastq/${lib}.R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 31000000000 \
    --quantMode GeneCounts \
    --twopassMode Basic \
    --outFileNamePrefix star/${lib}.

# Rename bam file (for cleaner @RG info)
mv star/${lib}.Aligned.sortedByCoord.out.bam star/${lib}.bam

# Index bam file
samtools index -@ $SLURM_CPUS_ON_NODE star/${lib}.bam