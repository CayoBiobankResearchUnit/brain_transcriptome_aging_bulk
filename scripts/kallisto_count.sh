#!/bin/bash

index=genomes/Mmul_10.cdna.idx

# Assign a library ID
## samples.txt is a text file with each library ID listed (must be identical to prefixes used for fastq files)
lib=$(tail -n+2 data/cayo_brain_bulk_metadata_technical.tsv | cut -f 1 | sed -n ${SLURM_ARRAY_TASK_ID}p)

mkdir -p kallisto

# Map to transcriptome

## All reads are in the fastq/ directory with the naming convention ${lib}.R<read: 1 or 2>.fastq.gz

kallisto quant -i $index -t $SLURM_CPUS_ON_NODE -o kallisto/${lib} \
	fastq/${lib}.R1.fastq.gz fastq/${lib}.R2.fastq.gz