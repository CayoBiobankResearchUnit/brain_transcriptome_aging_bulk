#!/bin/bash

i=$(tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | cut -f 1 | sed -n ${SLURM_ARRAY_TASK_ID}p)

srun="srun --exclusive -N1 -n1"

mkdir -p logs

parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog logs/parallel_${SLURM_JOB_ID}.log --resume"

$parallel "$srun scripts/samtools_extract_chromosome.sh {1} {2} > logs/parallel.${SLURM_JOB_ID}_{1}.log" ::: $i ::: {1..20}

exit
