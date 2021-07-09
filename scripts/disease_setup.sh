#!/bin/bash

wget -O data/human_disease_textmining_filtered.tsv http://download.jensenlab.org/human_disease_textmining_filtered.tsv

cut -f 1-6 data/human_disease_textmining_filtered.tsv > data/human_disease_associations.tsv

wget -O data/curated_gene_disease_associations.tsv.gz https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz

gunzip data/curated_gene_disease_associations.tsv.gz