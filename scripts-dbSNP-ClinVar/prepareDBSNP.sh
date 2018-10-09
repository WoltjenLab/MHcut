#!/bin/bash

#SBATCH --job-name=prepareDBSNP
#SBATCH --output=prepareDBSNP.out
#SBATCH --error=prepareDBSNP.out
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH --account=rrg-bourqueg-ad

module load jmonlong/python
python prepareVcfIndel.py -vcf dbsnp-del.vcf.gz -o dbsnp-del.tsv -info "RS|0,CAF|1,TOPMED|1,GENEINFO|0,PM|0" -addchr -infomerge 'MC|NSF,NSM,NSN,SYN,U3,U5,ASS,DSS,INT,R3,R5'
