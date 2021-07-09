#!/bin/sh

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

do_single=$3

if [ $do_single -eq "1" ]; then
parallel -j $slots scripts/glmnet_predict.R {1} {2} {3} {4} ::: 0 ::: {1..36} ::: $1 ::: $2
else
parallel -j $slots scripts/glmnet_predict.R {1} {2} {3} {4} ::: {1..15} ::: {1..36} ::: $1 ::: $2
fi

exit