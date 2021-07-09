#!/bin/bash

mkdir -p genomes

# Fetch fasta from Ensembl
wget -O genomes/Mmul_10.cdna.fa.gz \
	ftp://ftp.ensembl.org/pub/release-99/fasta/macaca_mulatta/cdna/Macaca_mulatta.Mmul_10.cdna.all.fa.gz

# Unzip fasta
gunzip genomes/Mmul_10.cdna.fa.gz

# Index fasta
kallisto index -i genomes/Mmul_10.cdna.idx genomes/Mmul_10.cdna.fa