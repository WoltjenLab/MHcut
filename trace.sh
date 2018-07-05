#### Prepare the genome file

## Download, unzip and index reference genome
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
python indexFasta.py hg38.fa

## Build Blast database
makeblastdb -in hg38.fa -dbtype nucl -title hg38



#### Deletions from ClinVar: small and larger deletions in the entire genome

## Download ClinVar variants
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

## Select pathogenic deletion with GRCh38 coordinates.
## Some records have 'na' for chromosome, don't take those.
## Also add a 'chr' prefix to the chromosome name to match the genome reference fasta.
## Keep the gene name, RS/dbSNP, dbVar and extra columns.
zcat variant_summary.txt.gz | awk 'BEGIN{FS="\t"; OFS="\t"; print "Chromosome\tStart\tStop\tGeneSymbol\tdbSNP\tdbVar\treviewStatus\tnbSubmitters"}{if($17=="GRCh38" && $19!="na" && $2=="deletion" && $7=="Pathogenic"){print "chr"$19,$20,$21,$5,$10,$11,$25,$26}}' | sort -u > clinvar-grch38-pathogenic-deletion.tsv

## Optional: check if duplicated variants (at the same position)
wc -l clinvar-grch38-pathogenic-deletion.tsv
cut -f 1-3 clinvar-grch38-pathogenic-deletion.tsv | sort -u | wc -l
## -> A few are duplicated, REMEMBER TO CHECK THAT IN MHCUT OUTPUT

## Run MHcut
python MHcut.py -var clinvar-grch38-pathogenic-deletion.tsv -ref hg38.fa -out clinvar-grch38-pathogenic-deletion

### To run MHcut on Janin's PC
cd ~/Documents/02\ Research\ project/07\ MHcut/01\ Data
python2.7 ../MHcutGitHub/MHcut.py -var variant_summary-grch38-pathogenic-deletion.tsv -ref hg38.fa -out ../"02 Output"/clinvar-grch38-pathogenic-deletion

#### Do not copy commands from here into Terminal, “ “ spaces are not formatted correctly


## Comprehensive run including MH down to 1 bp.
python MHcut.py -var clinvar-grch38-pathogenic-deletion.tsv -ref hg38.fa -out clinvar-grch38-pathogenic-deletion-1bp -minMHL 1 -minm1L 1 -minvarL 1


## Select all deletions with GRCh38 coordinates.
## Some records have 'na' for chromosome, don't take those.
## Also add a 'chr' prefix to the chromosome name to match the genome reference fasta.
## Keep the gene name, RS/dbSNP, dbVar and extra columns.
zcat variant_summary.txt.gz | awk 'BEGIN{FS="\t"; OFS="\t"; print "Chromosome\tStart\tStop\tGeneSymbol\tdbSNP\tdbVar\tClinicalSignificance\treviewStatus\tnbSubmitters"}{if($17=="GRCh38" && $19!="na" && $2=="deletion"){print "chr"$19,$20,$21,$5,$10,$11,$7,$25,$26}}' | sort -u > clinvar-grch38-all-deletion.tsv

python MHcut.py -var clinvar-grch38-all-deletion.tsv -ref hg38.fa -out clinvar-grch38-all-deletion-1bp -minMHL 1 -minm1L 1 -minvarL 1



#### Frameshift indels in dbSNP
curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz | gunzip -c | grep "VC=DIV" | grep -e ";NSF;" -e ";NSM;" -e ";NSN;" -e ";SYN;" | awk '{if(length($4)>length($5)){print $0}}' > dbsnp-nsf-div-del.vcf

## ADD PM/PMC info
python prepareVcfIndel.py -vcf dbsnp-nsf-div-del.vcf -o dbsnp-nsf-div-del.tsv -info "RS|0,CAF|1,GENEINFO|0" -addchr

python MHcut.py -var dbsnp-nsf-div-del.tsv -ref hg38.fa -out dbsnp-nsf-div-del-1bp -minMHL 1 -minm1L 1 -minvarL 1


### Exhaustive run testing both flanks (instead of the one with stronger homology)
python MHcut-2flanks.py -var clinvar-grch38-pathogenic-deletion.tsv -ref hg38.fa -out clinvar-grch38-pathogenic-deletion-2flanks
