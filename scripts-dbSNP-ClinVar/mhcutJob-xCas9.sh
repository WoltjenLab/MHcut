#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=3000M
#SBATCH --account=rrg-bourqueg-ad

echo "Start:`date`"

if [ $SLURM_ARRAY_TASK_ID -ne "1" ]; then
    head -1 $1.tsv > temp-mhcut-$1-$SLURM_ARRAY_TASK_ID.tsv
    split -d -n l/$SLURM_ARRAY_TASK_ID/600 $1.tsv >> temp-mhcut-$1-$SLURM_ARRAY_TASK_ID.tsv
else
    split -d -n l/$SLURM_ARRAY_TASK_ID/600 $1.tsv > temp-mhcut-$1-$SLURM_ARRAY_TASK_ID.tsv
fi

module load mugqic/jellyfish/2.2.6

MHcut -var temp-mhcut-$1-$SLURM_ARRAY_TASK_ID.tsv -ref hg38.fa -jf mer_counts23.jf -out mhcut-$1-xCas9-$SLURM_ARRAY_TASK_ID -nofilt -minvarL 1 -PAM NGN,GAA,GAT -indelphi -restart

rm temp-mhcut-$1-$SLURM_ARRAY_TASK_ID.tsv

echo "Finish:`date`"
