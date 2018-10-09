#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=3000M
#SBATCH --account=rrg-bourqueg-ad

echo "Start:`date`"

if [ $2 -ne "1" ]; then
    head -1 $1.tsv > temp-mhcut-$1-$2.tsv
    split -d -n l/$2/30 $1.tsv >> temp-mhcut-$1-$2.tsv
else
    split -d -n l/$2/30 $1.tsv > temp-mhcut-$1-$2.tsv
fi

module load jmonlong/python mugqic/jellyfish/2.2.6
python MHcut.py -var temp-mhcut-$1-$2.tsv -ref hg38.fa -jf mer_counts23.jf -out mhcut-$1-2MM-$2 -nofilt -minvarL 1 -maxConsMM 2

echo "Finish:`date`"
