#!/bin/bash

i=$(ls star/*.Aligned.sortedByCoord.out.bam | sed -n ${1}p | sed 's/star\/\([A-Z0-9_]*\)\.Aligned.sortedByCoord.out.bam/\1/g')

mkdir -p stats

samtools flagstat star/${i}.Aligned.sortedByCoord.out.bam > stats/${i}.flagstat.txt

fastq-stats fastq/${i}.R1.fastq.gz > stats/${i}.fastq-stats.txt
