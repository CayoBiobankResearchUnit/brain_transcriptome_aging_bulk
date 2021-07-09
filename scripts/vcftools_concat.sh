#!/bin/bash

vcf-concat $(echo vcf/cayo_brain_chr{1..20}_raw.vcf) | gzip -c > vcf/cayo_brain_all_raw.vcf.gz

vcf-concat $(echo vcf/cayo_brain_chr{1..20}_filtered.vcf) | gzip -c > vcf/cayo_brain_all_filtered.vcf.gz
