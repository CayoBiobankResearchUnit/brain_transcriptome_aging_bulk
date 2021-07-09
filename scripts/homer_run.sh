#!/bin/bash

export PATH=$(pwd)/software/bin:$PATH

this=$1

findMotifs.pl results/mmul10_best_genes_${this}.txt \
	mmul10 \
	results/homer_age_tf_${this} \
	-start -2000 -end 2000 \
	-bg results/mmul10_background_genes_${this}.txt \
	-mset vertebrates -nomotif -nogo
