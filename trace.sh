#### Prepare the genome file
## Download, unzip and index reference genome
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
python indexFasta.py hg38.fa

## Build Blast database
makeblastdb -in hg38.fa -dbtype nucl -title hg38

###
###
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

###
###
## Select all deletions with GRCh38 coordinates.
## Some records have 'na' for chromosome, don't take those.
## Also add a 'chr' prefix to the chromosome name to match the genome reference fasta.
## Keep the gene name, RS/dbSNP, dbVar and extra columns.
zcat variant_summary.txt.gz | awk 'BEGIN{FS="\t"; OFS="\t"; print "Chromosome\tStart\tStop\tGeneSymbol\tdbSNP\tdbVar\tClinicalSignificance\treviewStatus\tnbSubmitters"}{if($17=="GRCh38" && $19!="na" && $2=="deletion"){print "chr"$19,$20,$21,$5,$10,$11,$7,$25,$26}}' | sort -u > clinvar-grch38-all-deletion.tsv

python MHcut.py -var clinvar-grch38-all-deletion.tsv -ref hg38.fa -out clinvar-grch38-all-deletion-1bp -minMHL 1 -minm1L 1 -minvarL 1

###
###
##
### All dbSNP
curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz | gunzip -c | grep "VC=DIV" | awk '{if(length($4)>length($5)){print $0}}' | gzip > dbsnp-del.vcf.gz
python prepareVcfIndel.py -vcf dbsnp-del.vcf.gz -o dbsnp-gene.tsv -info "RS|0,CAF|1,GENEINFO|0,PM|0" -infomerge 'MC|NSF,NSM,NSN,SYN,U3,U5,ASS,DSS,INT,R3,R5' -addchr
gzip dbsnp-del.tsv

### All dbSNP aroung genes and ClinVar
curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz | gunzip -c | grep "VC=DIV"  | grep -E ";NSF;|;NSM;|;NSN;|;SYN;|;U3;|;U5;|;R3;|;R5;|;INT;|;ASS;|;RSS;|;REF;|;NSF$|;NSM$|;NSN$|;SYN$|;U3$|;U5$|;R3$|;R5$|;INT$|;ASS$|;RSS$|;REF$" | awk '{if(length($4)>length($5)){print $0}}' | gzip > dbsnp-gene-del.vcf.gz
python prepareVcfIndel.py -vcf dbsnp-gene-del.vcf.gz -o dbsnp-gene-del.tsv -info "RS|0,CAF|1,GENEINFO|0,PM|0" -addchr
gzip dbsnp-gene-del.tsv

wget ftp://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20180701.vcf.gz
python prepareVcfIndel.py -vcf clinvar_20180701.vcf.gz -o clinvar_20180701.tsv -info "AF_EXAC|0,AF_TGP|0,ALLELEID|0,CLNDN|-,CLNSIG|-,CLNVC|0,DBVARID|-,GENEINFO|1,MC|-,RS|0" -addchr
wget -O var_citations_17July2018.txt ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt
python prepareClinVar.py -tsv clinvar_20180701.tsv -cit var_citations_17July2018.txt -o clinvar-grch38-deletion_20180701.tsv

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz

# Merge the ClinVar and dbSNP and annotate with Gencode
Rscript mergeClinVarDbSNPGencode.R clinvar-grch38-deletion.tsv.gz dbsnp-gene-del.tsv gencode.v28.annotation.gtf.gz dbsnp-clinvar-deletion.tsv

## Run the faster version that uses JellyFish khmer counting
jellyfish count -m 20 -s 100M -t 10 hg38.fa ## Prepare count table

python MHcut.py -var dbsnp-clinvar-deletion.tsv.gz -ref hg38.fa -jf mer_counts.jf -out dbsnp-clinvar-deletion

###
###
### Exhaustive run testing both flanks (instead of the one with stronger homology)
python MHcut-2flanks.py -var clinvar-grch38-pathogenic-deletion.tsv -ref hg38.fa -out clinvar-grch38-pathogenic-deletion-2flanks
