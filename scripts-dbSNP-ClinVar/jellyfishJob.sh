#!/bin/bash

#SBATCH --job-name=jellyfish
#SBATCH --output=jellyfish.out
#SBATCH --error=jellyfish.out
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=128000M
#SBATCH --account=rrg-bourqueg-ad

module load mugqic/jellyfish/2.2.6
jellyfish count -m 23 --out-counter-len 1 -C -s 1000M -o mer_counts23.jf hg38.fa
