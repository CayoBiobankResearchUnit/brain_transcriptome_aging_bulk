#!/bin/bash

ls -1 bam/*.bam | grep -v chr > data/merged_genotypes.bamlist.txt

scripts/angsd_calculate_regions.R