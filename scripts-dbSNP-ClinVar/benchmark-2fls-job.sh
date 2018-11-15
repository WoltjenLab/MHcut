#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=3000M
#SBATCH --account=rrg-bourqueg-ad

echo "Start:`date`"

module load mugqic/jellyfish/2.2.6 python/2.7.14

MHcut -var dbsnp-clinvar-deletion-100KforBenchmark.tsv -ref hg38.fa -jf mer_counts23.jf -out dbsnp-clinvar-deletion-100KforBenchmark-1CMM-jellyfish-2fls -nofilt -minvarL 1 -2fls

echo "Finish:`date`"
