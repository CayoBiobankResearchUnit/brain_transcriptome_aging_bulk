#!/bin/bash

rm -rf data/library_sequencing_stats.txt;
echo library$'\t'reads$'\t'read_length$'\t'duplicates$'\t'perc_duplicate$'\t'qual_min$'\t'qual_max$'\t'qual_mean$'\t'perc_A$'\t'perc_C$'\t'perc_G$'\t'perc_T$'\t'bam_total$'\t'bam_secondary$'\t'reads_mapped$'\t'perc_mapped$'\t'splices_total$'\t'reads_multimapped >> data/library_sequencing_stats.txt;
for i in $(ls star/*.Aligned.sortedByCoord.out.bam | sed 's/star\/\([A-Z0-9_]*\)\.Aligned.sortedByCoord.out.bam/\1/g'); do
paste <(echo $i) \
	<(grep reads stats/${i}.fastq-stats.txt | cut -f 2) \
	<(grep len$'\t' stats/${i}.fastq-stats.txt | cut -f 2) \
	<(grep dups stats/${i}.fastq-stats.txt | cut -f 2) \
	<(grep '%dup' stats/${i}.fastq-stats.txt | cut -f 2) \
	<(grep 'qual min' stats/${i}.fastq-stats.txt | cut -f 2) \
	<(grep 'qual max' stats/${i}.fastq-stats.txt | cut -f 2) \
	<(grep 'qual mean' stats/${i}.fastq-stats.txt | cut -f 2) \
	<(grep '%A' stats/${i}.fastq-stats.txt | cut -f 2) \
	<(grep '%C' stats/${i}.fastq-stats.txt | cut -f 2) \
	<(grep '%G' stats/${i}.fastq-stats.txt | cut -f 2) \
	<(grep '%T' stats/${i}.fastq-stats.txt | cut -f 2) \
	<(grep total stats/${i}.flagstat.txt | cut -d ' ' -f 1) \
	<(grep secondary stats/${i}.flagstat.txt | cut -d ' ' -f 1) \
	<(grep 'Uniquely mapped reads number' star_out/${i}Log.final.out | cut -f 2) \
	<(grep 'Uniquely mapped reads %' star_out/${i}Log.final.out | cut -f 2 | sed 's/%//g') \
	<(grep 'Number of splices: Total' star_out/${i}Log.final.out | cut -f 2) \
	<(grep 'Number of reads mapped to multiple loci' star_out/${i}Log.final.out | cut -f 2) \
	>> data/library_sequencing_stats.txt;
done;