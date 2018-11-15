#!/bin/bash

#SBATCH --job-name=makeblastdb
#SBATCH --output=makeblastdb.out
#SBATCH --error=makeblastdb.out
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=8000M
#SBATCH --account=rrg-bourqueg-ad

module load mugqic/blast/2.7.1+
makeblastdb -in hg38.fa -dbtype nucl -title hg38

