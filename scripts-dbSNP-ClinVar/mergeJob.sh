#!/bin/bash

#SBATCH --job-name=merge
#SBATCH --output=merge.out
#SBATCH --error=merge.out
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH --account=rrg-bourqueg-ad

module load mugqic/R_Bioconductor
Rscript mergeClinVarDbSNPGencode.R clinvar-grch38-deletion_20180701.tsv dbsnp-del.tsv gencode.v28.annotation.gtf.gz dbsnp-clinvar-deletion.tsv
