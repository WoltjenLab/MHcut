Overlap known SNPs with MHcut's output
======================================

``` r
## For data manipulation
library(data.table)  # good for big data
library(dplyr)  # To manipulate/transform data.frames
library(magrittr)  # pipes (e.g. %>%)
## For tables
library(knitr)

library(GenomicRanges)

## Read the first row (headers) to remind us the order of each column
var = read.table("../scripts-dbSNP-ClinVar/mhcut-dbsnp-clinvar-deletion-variants.tsv.gz", 
    nrows = 1)
var = fread("gunzip -c ../scripts-dbSNP-ClinVar/mhcut-dbsnp-clinvar-deletion-variants.tsv.gz", 
    select = c(1:4, 20, 21, 23, 24))
head(var)
```

    ##     chr start   end         RS varL flank mhL mh1L
    ## 1: chr1 10109 10114 1377973775    6     2   6    6
    ## 2: chr1 10110 10114 1462685959    5     1   5    1
    ## 3: chr1 10120 10120 1156821933    1     2   0    0
    ## 4: chr1 10129 10147 1457723673   19     2  19   19
    ## 5: chr1 10132 10132 1289482855    1     2   0    0
    ## 6: chr1 10133 10137 1390118706    5     2   0    0

Subset of variants
------------------

Let's work with a random subset of 100K variants. This should be enough to produce good estimates.

Let's say we are interested in variants with at least 3 bp of mh1L.

``` r
var = var[sample.int(nrow(var), 1e+05)]
var.s = as.data.frame(subset(var, mh1L > 2))
var.s %<>% mutate(id = paste(chr, start, end, RS, sep = "_"))
```

Micro-homology regions
----------------------

``` r
fl1 = with(var.s, GRanges(var.s$chr, IRanges(ifelse(flank == 1, start - mhL, 
    start), ifelse(flank == 1, start, start + mhL)), id = id, flank = 1))
fl2 = with(var.s, GRanges(var.s$chr, IRanges(ifelse(flank == 1, end - mhL, end), 
    ifelse(flank == 1, end, end + mhL)), id = id, flank = 2))
fl.gr = c(fl1, fl2)
```

Common SNPs in the human population
-----------------------------------

``` r
snps = fread("gunzip -c ../scripts-dbSNP-ClinVar/snp151Common-formatted.tsv.gz")
snps = as.data.frame(snps)
colnames(snps) = c("chr", "start", "end", "RS", "vtype", "freq")
dim(snps)
```

    ## [1] 15175044        6

``` r
snps = subset(snps, !(RS %in% paste0("rs", var.s$RS)))
dim(snps)
```

    ## [1] 15174576        6

``` r
snps = makeGRangesFromDataFrame(snps, keep.extra.columns = TRUE)
```

Overlap
-------

``` r
sum.df = lapply(unique(snps$vtype), function(vt) {
    fl.gr$snp.ol = overlapsAny(fl.gr, subset(snps, vtype == vt))
    fl.id = fl.gr %>% as.data.frame %>% group_by(id) %>% summarize(snp.ol = any(snp.ol))
    fl.id %>% ungroup %>% summarize(prop.ol = mean(snp.ol)) %>% mutate(vtype = vt)
})
sum.df = do.call(rbind, sum.df)

sum.df %>% arrange(desc(prop.ol)) %>% kable
```

|    prop.ol| vtype          |
|----------:|:---------------|
|  0.0880993| deletion       |
|  0.0824078| single         |
|  0.0151773| insertion      |
|  0.0034386| in-del         |
|  0.0021738| microsatellite |
|  0.0000000| mnp            |
|  0.0000000| named          |
