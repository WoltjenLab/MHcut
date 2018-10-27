Overlap gene annotation with MHcut's output
===========================================

There are several ways of comparing variants and gene annotations. Sometimes there are some gene information in the input data, like in our dbSNP/ClinVar analysis which had a *GENEINFO* (from dbSNP) and *GENEINFO.ClinVar*. These columns could be used directly to identify variants around a specific gene or count how many genes have variants around them. However these gene information might not be exactly the one we want. Another approach is to directly overlap variant location with a gene annotation (e.g. Gencode). This report shows how to do both.

Let say the goal is to **count how many genes have an exonic variant**.

First we load packages and MHcut's output.

``` r
## For data manipulation
library(data.table)  # good for big data
library(dplyr)  # To manipulate/transform data.frames
library(magrittr)  # pipes (e.g. %>%)
## For tables
library(knitr)

## Read the first row (headers) to remind us the order of each column
var = read.table("../data/mhcut-dbsnp-clinvar-deletion-variants.tsv.gz", nrows = 1)
## Import variants coordinates and colums mh1L, GENEINFO and geneloc (fread
## is data.table command)
var = fread("gunzip -c ../data/mhcut-dbsnp-clinvar-deletion-variants.tsv.gz", 
    select = c(1:3, 7, 19, 22))
head(var)
```

    ##     chr start   end          GENEINFO    geneloc mh1L
    ## 1: chr1 10109 10114 DDX11L1:100287102 intergenic    6
    ## 2: chr1 10110 10114 DDX11L1:100287102 intergenic    1
    ## 3: chr1 10120 10120 DDX11L1:100287102 intergenic    0
    ## 4: chr1 10129 10147 DDX11L1:100287102 intergenic   19
    ## 5: chr1 10132 10132 DDX11L1:100287102 intergenic    0
    ## 6: chr1 10133 10137 DDX11L1:100287102 intergenic    0

*geneloc* can take the following values:

``` r
unique(var$geneloc)
```

    ## [1] "intergenic" "intronic"   "exonic"

Let's say we are interested in exonic variants with at least 3 bp of mh1L.

``` r
var3ex = subset(var, mh1L > 2 & geneloc == "exonic")
```

Working on the GENEINFO column
------------------------------

We could split the GENEINFO to get gene names for our subset of variants and count the unique gene names.

``` r
var3ex.genes.l = strsplit(var3ex$GENEINFO, "\\|")
var3ex.genes.uniq = unique(unlist(var3ex.genes.l))
length(var3ex.genes.uniq)
```

    ## [1] 30700

### Problem 1: some exonic variants have empty GENEINFO

GENEINFO comes from dbSNP while the *geneloc* annotation from overlapping Gencode annotation so some inconsistencies exists. For example:

``` r
subset(var3ex, GENEINFO == "-") %>% head %>% kable
```

| chr  |   start|     end| GENEINFO | geneloc |  mh1L|
|:-----|-------:|-------:|:---------|:--------|-----:|
| chr1 |   63031|   63038| -        | exonic  |     7|
| chr1 |   63034|   63036| -        | exonic  |     3|
| chr1 |   63736|   63738| -        | exonic  |     3|
| chr1 |   63736|   63738| -        | exonic  |     3|
| chr1 |  134720|  134723| -        | exonic  |     3|
| chr1 |  297438|  297441| -        | exonic  |     3|

### Problem 2: multiple genes in GENEINFO

When there are multiple genes in GENEINFO, are they all the ones with coding regions overlapping the variant?

Overlap the variant with Gencode annotation
-------------------------------------------

GENEINFO might be useful to quickly look for genes in the large dataset, but to compute numbers or if we know exactly which annotation we want to use we might better compare the annotation directly with the variants. In the following we load Gencode annotation, and overlap the exons with the variants selected earlier.

First, download Gencode v28 (if not already there) and import it.

``` r
if (!file.exists("../data/gencode.v28.annotation.gtf.gz")) {
    download.file("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz", 
        "../data/gencode.v28.annotation.gtf.gz")
}
genc = read.table("../data/gencode.v28.annotation.gtf.gz", as.is = TRUE, sep = "\t")
colnames(genc) = c("chr", "source", "type", "start", "end", "score", "strand", 
    "phase", "attributes")
genc %>% head(1) %>% kable
```

| chr  | source | type |  start|    end| score | strand | phase | attributes                                                                                                                                   |
|:-----|:-------|:-----|------:|------:|:------|:-------|:------|:---------------------------------------------------------------------------------------------------------------------------------------------|
| chr1 | HAVANA | gene |  11869|  14409| .     | +      | .     | gene\_id ENSG00000223972.5; gene\_type transcribed\_unprocessed\_pseudogene; gene\_name DDX11L1; level 2; havana\_gene OTTHUMG00000000961.2; |

We might want to parse the gene type and gene name out of the *attributes* column.

``` r
genc$genetype = gsub(".*gene_type ([^;]*);.*", "\\1", genc$attributes)
genc$genename = gsub(".*gene_name ([^;]*);.*", "\\1", genc$attributes)
```

Now let's create a GRanges object with the exons of *protein\_coding* genes, keeping the columns (to keep the gene names that we'll need later).

``` r
library(GenomicRanges)
exons.gr = genc %>% filter(type == "exon", genetype == "protein_coding") %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
length(unique(exons.gr$genename))
```

    ## [1] 19883

We have exon coordinates for 19883 protein-coding genes. We can now filter these exons to keep only the ones overlapping the variants we selected previously.

``` r
var3ex.gr = makeGRangesFromDataFrame(var3ex)
exons3ex.gr = subsetByOverlaps(exons.gr, var3ex.gr)
length(unique(exons3ex.gr$genename))
```

    ## [1] 19348

19348 genes have a variant in exonic region.
