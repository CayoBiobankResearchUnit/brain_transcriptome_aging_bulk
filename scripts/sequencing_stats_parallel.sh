#!/bin/bash

srun="srun --exclusive -N1 -n1"

mkdir -p logs

parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog logs/parallel_${SLURM_JOB_ID}.log --resume"

$parallel "$srun scripts/sequencing_stats.sh {1} > logs/parallel.${SLURM_JOB_ID}_{1}.log" ::: $(eval echo -e {1..$(ls star/*.Aligned.sortedByCoord.out.bam | wc -l)})

exit
