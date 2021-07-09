#!/bin/bash

scripts/mashr_predict_cv.R ${SLURM_ARRAY_TASK_ID} 0 $1 $2
