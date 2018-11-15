## Download dbSNP deletions and transform to TSV
curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz | gunzip -c | grep "VC=DIV" | awk '{if(length($4)>length($5)){print $0}}' | gzip > dbsnp-del.vcf.gz
sbatch prepareDBSNP.sh

## Download ClinVar and select relevant information
wget ftp://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20180701.vcf.gz
python prepareVcfIndel.py -vcf clinvar_20180701.vcf.gz -o clinvar_20180701.tsv -info "AF_EXAC|0,AF_TGP|0,ALLELEID|0,CLNDN|-,CLNSIG|-,CLNVC|0,DBVARID|-,GENEINFO|1,MC|-,RS|0" -addchr
wget -O var_citations_17July2018.txt ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt
python prepareClinVar.py -tsv clinvar_20180701.tsv -cit var_citations_17July2018.txt -o clinvar-grch38-deletion_20180701.tsv

## Merge ClinVar and dbSNP, and annotate with gencode
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz
sbatch mergeJob.sh

## Check duplicates
wc -l dbsnp-clinvar-deletion.tsv
cut -f 1-4 dbsnp-clinvar-deletion.tsv | sort -u | wc -l

## Remove duplicates
head -1 dbsnp-clinvar-deletion.tsv > temp.tsv
sed 1d dbsnp-clinvar-deletion.tsv | sort -k1,4 -V -u >> temp.tsv
mv temp.tsv dbsnp-clinvar-deletion.tsv

## Get reference genome
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

## Precount khmer with JellyFish
sbatch jellyfishJob.sh

## Split input and run MHcut in parallel
INFILE=dbsnp-clinvar-deletion
for CHUNK in `seq 1 30`
do
    sbatch -J mhcut-$INFILE-$CHUNK -o mhcut-$INFILE-$CHUNK.out -e mhcut-$INFILE-$CHUNK.out mhcutJob.sh $INFILE $CHUNK
done

## If necessary, rerun the ones that didn't finish
REDO=`grep -c "100%" mhcut-dbsnp-clinvar-deletion-*out | awk 'BEGIN{FS=":"}{if($2==0){match($1, "deletion-(.*).out", a); print a[1]}}'`
INFILE=dbsnp-clinvar-deletion
for CHUNK in $REDO
do
    sbatch -J mhcut-$INFILE-$CHUNK -o mhcut-$INFILE-$CHUNK.out -e mhcut-$INFILE-$CHUNK.out mhcutJob.sh $INFILE $CHUNK
done

## Merge MHcut results
for TYPE in variants guides cartoons
do
    head -1 mhcut-dbsnp-clinvar-deletion-1-$TYPE.tsv > mhcut-dbsnp-clinvar-deletion-$TYPE.tsv
    for CHUNK in `seq 1 30`
    do
	sed 1d mhcut-dbsnp-clinvar-deletion-$CHUNK-$TYPE.tsv >> mhcut-dbsnp-clinvar-deletion-$TYPE.tsv
    done
    gzip -f mhcut-dbsnp-clinvar-deletion-$TYPE.tsv
done



#
##
### xCas9
##
#
## Split input and run MHcut in parallel
INFILE=dbsnp-clinvar-deletion
for CHUNK in `seq 1 30`
do
    sbatch -J mhcutxCas9-$INFILE-$CHUNK -o mhcutxCas9-$INFILE-$CHUNK.out -e mhcutxCas9-$INFILE-$CHUNK.out mhcutJob-xCas9.sh $INFILE $CHUNK
done

## Merge MHcut results
for TYPE in variants guides cartoons
do
    head -1 mhcut-dbsnp-clinvar-deletion-xCas9-1-$TYPE.tsv > mhcut-dbsnp-clinvar-deletion-xCas9-$TYPE.tsv
    for CHUNK in `seq 1 30`
    do
	sed 1d mhcut-dbsnp-clinvar-deletion-xCas9-$CHUNK-$TYPE.tsv >> mhcut-dbsnp-clinvar-deletion-xCas9-$TYPE.tsv
    done
    gzip -f mhcut-dbsnp-clinvar-deletion-xCas9-$TYPE.tsv
done




#
##
### Benchmark
##
#
# 1 vs 2 consecutive mismatches
# BLAST vs JellyFish
# Best flank vs both flanks
#

## Prepare BLAST
sbatch makeblastdbJob.sh

## Randomly select 100K variants
head -1 dbsnp-clinvar-deletion.tsv > dbsnp-clinvar-deletion-100KforBenchmark.tsv
shuf -n 100000 dbsnp-clinvar-deletion.tsv >> dbsnp-clinvar-deletion-100KforBenchmark.tsv

## Default values: 1 cons MM, JellyFish, best flank
sbatch -J mhcut-bmrk-default -o mhcut-bmrk-default.out -e mhcut-bmrk-default.out benchmark-default-job.sh

## 2 cons MM
sbatch -J mhcut-bmrk-2cmm -o mhcut-bmrk-2cmm.out -e mhcut-bmrk-2cmm.out benchmark-2cmm-job.sh

## BLAST
sbatch -J mhcut-bmrk-blast -o mhcut-bmrk-blast.out -e mhcut-bmrk-blast.out benchmark-blast-job.sh

## Both flank
sbatch -J mhcut-bmrk-2fls -o mhcut-bmrk-2fls.out -e mhcut-bmrk-2fls.out benchmark-2fls-job.sh
