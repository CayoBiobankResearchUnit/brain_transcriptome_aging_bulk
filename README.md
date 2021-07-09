# Brain bulk-tissue RNA-seq analysis pipeline

**Repository for the aging brain transcriptome project from Cayo Santiago macaques (bulk-tissue RNA-seq)**

This repository contains scripts used in the analysis of age effects in bulk-tissue RNA-seq data for the Cayo Santiago rhesus macaque population. Corresponding code use in the analysis of age effects in single-nucleus RNA-seq data can be found in the respository [CayoBiobankResearchUnit/brain\_transriptome\_aging\_sc](https://github.com/CayoBiobankResearchUnit/brain_transcriptome_aging_sc).

Note that we ran most steps on the University of Washington ([Mox](https://wiki.cac.washington.edu/display/hyakusers/Hyak+mox+Overview)) and Arizona State University ([Agave](https://cores.research.asu.edu/research-computing/user-guide)) high-performance computing clusters. We have aimed to generalize the code here by removing system-specific references to installed software and modules. Instead, we document required software and version numbers below (excluding standard Unix programs and R). For HPC systems, the required scripts and binaries must be in the PATH. The easiest way to do this is to use an existing module or to install your own. In these cases, the modules should be loaded prior to running the appropriate code below.

As Mox and Agave use the [slurm](https://slurm.schedmd.com/documentation.html) scheduler, most code below should run on slurm systems with little or no modification. For non-slurm HPC systems, slurm scripts and environmental variables will need to be adjusted, though hopefully without too much hassle.

We ran most analysis steps using [R](https://cran.r-project.org) (v4.1). We recommend the following utility or visualization packages to extend base R's functionality.

| Package                                                              | Description                                        |
| -----------                                                          | -----------                                        |
| [tidyverse](https://www.tidyverse.org/)                              | utilities for data manipulation and visualization  |
| [reshape2](https://cran.r-project.org/web/packages/reshape2)         | data manipulation                                  |
| [abind](https://cran.r-project.org/web/packages/abind)               | combining multidimensional arrays                  |
| [XML](https://cran.r-project.org/web/packages/XML)                   | parsing XML files                                  |
| [jsonlite](https://cran.r-project.org/web/packages/jsonlite)         | parsing JSON files                                 |
| [ggrastr](https://cran.r-project.org/web/packages/ggrastr)           | rasterizing big data visualizations                |
| [ggtext](https://cran.r-project.org/web/packages/ggtext)             | rendering text                                     |
| [ggrepel](https://cran.r-project.org/web/packages/ggrepel)           | avoiding overplotting of labels                    |
| [ggbeeswarm](https://cran.r-project.org/web/packages/ggbeeswarm)     | beeswarm plot support                              |
| [egg](https://cran.r-project.org/web/packages/egg)                   | advanced layouts for visualization                 |
| [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer) | color advice for plots                             |
| [viridis](https://cran.r-project.org/web/packages/viridis)           | perceptually uniform color gradients               |
| [doParallel](https://cran.r-project.org/web/packages/doParallel)     | support for parallel computing                     |
| [future](https://cran.r-project.org/web/packages/future)             | support for parallel computing                     |

More specialized R packages are listed with their specific scripts below.

# Inputs

The following files are expected:

* Demultiplexed pair-end fastq files should be compressed with gzip and placed in the `fastq/` folder with the naming convention `<library ID>.R1.fastq.gz` (read 1) and `<library ID>.R2.fastq.gz` (read 2).

* An animal metadata file should be placed in `data/cayo_brain_bulk_metadata_animals.tsv`

* A library metadata file should be placed in `data/cayo_brain_bulk_metadata_technical.tsv`

* An animal social metrics file should be placed in `data/social_metrics.csv`

# Pipeline

## Map reads with splice-aware aligner

* ***Required software***: [STAR](https://github.com/alexdobin/STAR) (v2.5), [SAMtools](http://www.htslib.org/) (v1.9), [GATK](https://gatk.broadinstitute.org/) (v4.1.2.0)

```
# Download and index genome
sbatch scripts/star_index.sh

# Map each library using STAR
sbatch --array=1-$(tail -n+2 data/cayo_brain_bulk_metadata_technical.tsv | wc -l | xargs) scripts/star_map.sh
```

## Merge alignments per genotype

* ***Required software***: [SAMtools](http://www.htslib.org/) (v1.9)

```
# Merge bam files for each genotype
sbatch --array=1-$(tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | wc -l | xargs) scripts/samtools_merge.sh
```

## Call and filter genotypes

* ***Required software***: [SAMtools](http://www.htslib.org/) (v1.9), [GATK](https://gatk.broadinstitute.org/) (v4.1.2.0), [VCFtools](https://vcftools.github.io/) (v0.1.16)

```
# Split variants into chromosomes (per genotype)
sbatch --array=1-$(tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | wc -l | xargs) scripts/samtools_split.sh

# Clean variants
sbatch --array=1-$(($(tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | wc -l | xargs)*20)) scripts/gatk_clean_reads.sh

# Call variants
sbatch --array=1-20 scripts/gatk_call_variants.sh

# Filter variants
sbatch --array=1-20 scripts/gatk_filter_variants.sh

# Concatenate variants across chromosomes
scripts/vcftools_concat.sh
```

## Compute relatedness with lcMLkin

* ***Required software***: [VCFtools](https://vcftools.github.io/) (v0.1.16), [lcMLkin](https://github.com/COMBINE-lab/maximum-likelihood-relatedness-estimation) (v20190218)

```
# Thin variants
scripts/vcftools_thin_variants.sh

# Call kinship
scripts/lcmlkin_call_kinship.sh
```

## Quantify transcripts with kallisto

* ***Required software***: [kallisto](https://pachterlab.github.io/kallisto) (v0.43.1)

```
# Index transcriptome for kallisto
scripts/kallisto_index.sh

# Count transcripts
sbatch --array=1-$(tail -n+2 data/cayo_brain_bulk_metadata_technical.tsv | wc -l | xargs) scripts/kallisto_count.sh
```

## Calculate sequencing stats

* ***Required software***: [SAMtools](http://www.htslib.org/) (v1.9), [ea-utils](https://expressionanalysis.github.io/ea-utils/) (v1.04.807), [GNU parallel](https://www.gnu.org/software/parallel/) (v20171122)

```
# Calculate sequencing stats across libraries
sbatch scripts/sequencing_stats_parallel.sh

# Summarize sequencing stats
scripts/sequencing_stats_summarize.sh
```

## Import and clean sample and library metadata

```
# Read, format, and clean animal and library metadata
scripts/clean_metadata.R
```

## Import expression data

* ***Key libraries***: [tximport](https://doi.org/doi:10.18129/B9.bioc.tximport), [rdf5](https://doi.org/doi:10.18129/B9.bioc.rhdf5), [biomaRt](https://doi.org/doi:10.18129/B9.bioc.biomaRt)

```
# Import kallisto results into R
scripts/kallisto_import.R
```

## Filter expression matrix

* ***Key libraries***: [biomaRt](https://doi.org/doi:10.18129/B9.bioc.biomaRt), [limma](https://doi.org/doi:10.18129/B9.bioc.limma)


```
# Apply filters to gene expression dataset
scripts/filter_expression.R
```

## Visualize data (pre-modeling)

* ***Key libraries***: [variancePartition](https://doi.org/doi:10.18129/B9.bioc.variancePartition), [umap](https://cran.r-project.org/web/packages/umap), [dendextend](https://cran.r-project.org/web/packages/dendextend), [ape](https://cran.r-project.org/web/packages/ape), [phangorn](https://cran.r-project.org/web/packages/phangorn)

```
# Visualize metadata
scripts/visualize_metadata.R

# Visualize data
scripts/visualize_expression.R
```

## Fit linear mixed model(s)

* ***Key libraries***: [EMMREML](https://cran.r-project.org/web/packages/EMMREML)

```
# Fit linear mixed effect model
scripts/emma_model.R
```

## Apply adaptive shrinkage

* ***Key libraries***: [mashr](https://cran.r-project.org/web/packages/mashr)

```
# Refine effects and significance values with adaptive shrinkage
scripts/mashr_model.R
```

## Enrichment of age effects

* ***Key libraries***: [biomaRt](https://doi.org/doi:10.18129/B9.bioc.biomaRt), [topGO](https://doi.org/doi:10.18129/B9.bioc.topGO), [GO.db](https://doi.org/doi:10.18129/B9.bioc.GO.db)

```
# Perform GO enrichment analysis
scripts/topgo_enrichment.R

# Setup disease enrichment analysis
scripts/disease_setup.R

# Perform disease enrichment analysis
scripts/disease_enrichment.R
```

## Transcription factor binding motif enrichment

* ***Required software***: [HOMER](http://homer.ucsd.edu/homer) (v4.11)

```
scripts/homer_prep.R
scripts/homer_build.sh
for i in {1..2}; do scripts/homer_run.sh $i; scripts/homer_summarize.R $i; done
```

## Visualize results (post-modeling)

```
# Visualize model results
scripts/visualize_model_results.R
```

## Rerun models following cell-type deconvolution

* ***Required input***: `checkpoints/cayo_bulkbrain_cell_proportions.rds` (file with libraries and cell-type proportions estimated using BRETTIGEA).

* ***Key libraries***: [EMMREML](https://cran.r-project.org/web/packages/EMMREML), [mashr](https://cran.r-project.org/web/packages/mashr)

```
# Run EMMA models controlling for cell-type proportions
scripts/cellprop_emma.R

# Refine effects and significance values with adaptive shrinkage
scripts/cellprop_mashr.R

# Summarize and visualize model results
scripts/cellprop_summarize.R
```

## Gene trajectory analysis

* ***Key libraries***: [cluster](https://cran.r-project.org/web/packages/cluster), [clValid](https://cran.r-project.org/web/packages/clValid), [biomaRt](https://doi.org/doi:10.18129/B9.bioc.biomaRt), [topGO](https://doi.org/doi:10.18129/B9.bioc.topGO)

```
# Identify and cluster gene trajectories
scripts/gene_trajectories.R

# Visualize gene trajectories
scripts/gene_trajectories_plot.R
```

## Age effects on variance of gene expression

* ***Key libraries***: [dglm](https://cran.r-project.org/web/packages/dglm), [mashr](https://cran.r-project.org/web/packages/mashr), [biomaRt](https://doi.org/doi:10.18129/B9.bioc.biomaRt), [topGO](https://doi.org/doi:10.18129/B9.bioc.topGO), [GO.db](https://doi.org/doi:10.18129/B9.bioc.GO.db)

```
# Fit double generalized linear model(s)
scripts/dglm_model.R

# Refine effects and significance values with adaptive shrinkage
scripts/dglm_mashr.R

# Gene Ontology enrichment analysis
scripts/dglm_topgo.R

# DISEASES enrichment analysis
scripts/dglm_disease.R

# Visualize results
scripts/dglm_visualize.R
```

## Enrichment of Alzheimer's disease differentially expressed genes

* ***Required input***: `data/meta.anlz.ad_cntrl.tsv` (data table obtained from the [AD Knowledge portal](https://adknowledgeportal.synapse.org), Synapse:syn11914808).

```
scripts/ad_enrichment.R
```

## Compare effects of age and dominance rank

* ***Key libraries***: [mashr](https://cran.r-project.org/web/packages/mashr)

```
# Refine effects and significance values for dominance rank with adaptive shrinkage
scripts/mashr_model_social.R ordinal.rank.L

# Compare effects of age and dominance rank
scripts/compare_social_effects.R exact_age_years ordinal.rank.L
```

## Make predictions using fit models

* ***Key libraries***: [mashr](https://cran.r-project.org/web/packages/mashr), [glmnet](https://cran.r-project.org/web/packages/glmnet)

```
# Predict on same dataset (wbaDEGs model)
scripts/mashr_predict.R

# Perform cross-validation (k set to 0, which in this case codes for LOOCV)
sbatch --array=1-$(echo $(tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | wc -l | xargs)) scripts/mashr_predict_cv.slurm 0

# Summarize and visualize wbaDEGs cross-validation results
scripts/visualize_predictions_cv.R 0

# Predict with glmnet (with LOOCV)
sbatch scripts/glmnet_predict.slurm 0 0 1

# Summarize and visualize glmnet model results
scripts/glmnet_summarize.R 0 0 1
```

## Make prediction after filtering dominance rank effects (wbaDEGs)

* ***Key libraries***: [mashr](https://cran.r-project.org/web/packages/mashr)

```
# Run predictions with three different thresholds
scripts/mashr_predict.R 1 0.01
scripts/mashr_predict.R 1 0.001
scripts/mashr_predict.R 1 0.0001

# For each animal (left-out of training but used for evaluation), predict transcriptional age using filtered genes
sbatch --array=1-$(echo $(tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | wc -l | xargs)) scripts/mashr_predict_cv.slurm 1 0.01
sbatch --array=1-$(echo $(tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | wc -l | xargs)) scripts/mashr_predict_cv.slurm 1 0.001
sbatch --array=1-$(echo $(tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | wc -l | xargs)) scripts/mashr_predict_cv.slurm 1 0.0001

# Visualize predictions at different dominance rank thresholds
scripts/visualize_predictions_cv.R 0 1 0.01
scripts/visualize_predictions_cv.R 0 1 0.001
scripts/visualize_predictions_cv.R 0 1 0.0001
```

## Make predictions after filtering dominance rank effects (glmnet)

* ***Key libraries***: [glmnet](https://cran.r-project.org/web/packages/glmnet)

```
# Run predictions with three different thresholds
sbatch scripts/glmnet_predict.slurm 1 0.2 1
sbatch scripts/glmnet_predict.slurm 1 0.1 1
sbatch scripts/glmnet_predict.slurm 1 0.05 1

# Summarize and visualize predictions at different dominance rank thresholds
scripts/glmnet_summarize.R 1 0.2
scripts/glmnet_summarize.R 1 0.1
scripts/glmnet_summarize.R 1 0.05
```

## Compare all dominance-rank-filtered predictions

* ***Key libraries***: [lme4](https://cran.r-project.org/web/packages/lme4), [lmerTest](https://cran.r-project.org/web/packages/lmerTest)

```
scripts/visualize_predictions_comparison.R
```
