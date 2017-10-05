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
## Keep the gene name, RS/dbSNP and dbVar extra columns.
zcat variant_summary.txt.gz | awk 'BEGIN{FS="\t"; OFS="\t"; print "Chromosome\tStart\tStop\tGeneSymbol\tdbSNP\tdbVar\treviewStatus\tnbSubmitters"}{if($17=="GRCh38" && $19!="na" && $2=="deletion" && $7=="Pathogenic"){print "chr"$19,$20,$21,$5,$10,$11,$25,$26}}' | sort -u > clinvar-grch38-pathogenic-deletion.tsv

## Optional: check if duplicated variants (at the same position)
wc -l clinvar-grch38-pathogenic-deletion.tsv
cut -f 1-3 clinvar-grch38-pathogenic-deletion.tsv | sort -u | wc -l
## -> A few are duplicated, REMEMBER TO CHECK THAT IN MHCUT OUTPUT

## Run MHcut
python MHcut.py -var clinvar-grch38-pathogenic-deletion.tsv -ref hg38.fa -out clinvar-grch38-pathogenic-deletion

### To run MHcut on Janin's PC
cd ~/Documents/"02 Research project"/"07 MHcut”/“01 Data”
python2.7 ../MHcutGitHub/MHcut.py -var variant_summary-grch38-pathogenic-deletion.tsv -ref hg38.fa -out ../"02 Output"/clinvar-grch38-pathogenic-deletion

#### Do not copy commands from here into Terminal, “ “ spaces are not formatted correctly
