#!/bin/bash

#SBATCH --time=40:00:00
#SBATCH --mem-per-cpu=30000M
#SBATCH --account=rrg-bourqueg-ad

echo "Start:`date`"

module load python/2.7.14 mugqic/blast/2.7.1+

MHcut -var dbsnp-clinvar-deletion-100KforBenchmark.tsv -ref hg38.fa -out dbsnp-clinvar-deletion-100KforBenchmark-1CMM-blast-bestFl -nofilt -minvarL 1 -chunkS 10

echo "Finish:`date`"
