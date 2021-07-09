#!/bin/bash

mkdir -p genomes

# Fetch fasta from Ensembl
wget -O genomes/Mmul_10.dna.fa.gz \
	ftp://ftp.ensembl.org/pub/release-101/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz

# Unzip fasta
gunzip genomes/Mmul_10.dna.fa.gz

# Index fasta
samtools faidx genomes/Mmul_10.dna.fa

# Fetch gtf from Ensembl
wget -O genomes/Mmul_10.gtf.gz \
	ftp://ftp.ensembl.org/pub/release-101/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.101.gtf.gz

# Unzip gtf
gunzip genomes/Mmul_10.gtf.gz

STAR --runThreadN $SLURM_CPUS_ON_NODE \
     --runMode genomeGenerate \
     --genomeDir genomes \
     --genomeFastaFiles genomes/Mmul_10.dna.fa \
     --sjdbGTFfile genomes/Mmul_10.gtf \
     --sjdbOverhang 49

# Create sequence dictionary

gatk CreateSequenceDictionary --REFERENCE genomes/Mmul_10.dna.fa
