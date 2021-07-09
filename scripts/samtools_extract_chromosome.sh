#!/bin/bash

i=$1
chr=$2

samtools view -b bam/${i}.reheadered.bam $chr > bam/${i}.chr${chr}.bam
