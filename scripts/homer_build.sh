#!/bin/bash

mkdir -p software

wget -O software/configureHomer.pl http://homer.ucsd.edu/homer/configureHomer.pl

cd software

chmod +x configureHomer.pl

perl software/configureHomer.pl -install

# Download genome
wget -O genomes/Mmul_10.fa.gz http://ftp.ensembl.org/pub/release-103/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna_sm.toplevel.fa.gz
wget -O genomes/Mmul_10.gtf.gz http://ftp.ensembl.org/pub/release-103/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.103.gtf.gz
# wget -O genomes/Mmul_10.gff3.gz http://ftp.ensembl.org/pub/release-103/gff3/macaca_mulatta/Macaca_mulatta.Mmul_10.103.gff3.gz

gunzip genomes/Mmul_10.fa.gz
gunzip genomes/Mmul_10.gtf.gz

export PATH=$(pwd)/software/bin:$PATH

loadGenome.pl -force -name mmul10 -org null \
	-fasta genomes/Mmul_10.fa \
	-gtf genomes/Mmul_10.gtf \
	-version 10 -gid -promoters mmul10
