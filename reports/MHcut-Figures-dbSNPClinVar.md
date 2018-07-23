Figure sandbox
==============

First load the packages:

``` r
## For data manipulation
library(data.table)  # good for big data
library(dplyr)  # To manipulate/transform data.frames
library(magrittr)  # pipes (e.g. %>%)

## For graphs and tables
library(ggplot2)
library(knitr)

## Functions
winsor <- function(x, u = 10) {
    ## Winsorize: if larger than u, force value to u
    if (any(x > u)) 
        x[x > u] = u
    x
}
```

Then we’ll read the input TSV file (with 43M variants) using the `fread`
function. This will create a *data.table* object which uses disk space
rather than memory (if I understand correctly) so it’s good for big
data. To reduce the amount of data loaded in R, let’s only load some
specific columns (using `select=`).

``` r
## Read the first row (headers) to remind us the order of each column
read.table("../data/mhcut-dbsnp-clinvar-deletion-variants.tsv.gz", nrows = 1)
```

    ##    V1    V2  V3 V4  V5     V6       V7 V8 V9     V10    V11      V12   V13
    ## 1 chr start end RS CAF TOPMED GENEINFO PM MC AF_EXAC AF_TGP ALLELEID CLNDN
    ##      V14     V15              V16        V17      V18     V19  V20 V21
    ## 1 CLNSIG DBVARID GENEINFO.ClinVar MC.ClinVar citation geneloc varL mhL
    ##    V22 V23  V24    V25    V26    V27    V28     V29        V30         V31
    ## 1 mh1L hom nbMM mhDist MHseq1 MHseq2 pamMot pamUniq guidesNoOT guidesMinOT

``` r
## Import variants and colums varL, mhL, mh1L, mhDist
var = fread("gunzip -c ../data/mhcut-dbsnp-clinvar-deletion-variants.tsv.gz", 
    select = c(20:22, 25))
```

This takes about 1 min. Now the `var` object is a *data.table* with 43M
values. Before feeding this to ggplot2, the most efficient is to compute
the summary statistics using the `data.table` functions, then convert to
*data.frame* and call ggplot.

Let’s say we want the number of variants for each pair of {variant size,
MH length}, the *data.table* command is:

``` r
varmh = var[, .N, by = .(varL, mhL)]
head(varmh)
```

    ##    varL mhL       N
    ## 1:    6   6  577270
    ## 2:    5   4  316808
    ## 3:    1   0 6275487
    ## 4:   19  19   12577
    ## 5:    5   0  349498
    ## 6:   13  13   66297

``` r
nrow(varmh)
```

    ## [1] 2972

This is now just ~3,000 values so we can convert to *data.frame* and use
dplyr/ggplot to make the graph.

``` r
varmh.df = as.data.frame(varmh)

## Create a new column with the mhL class
varmh.df$mhL.class = cut(varmh.df$mhL, breaks = c(-1, 0, 1, 2, 3, 4, 10, Inf), 
    labels = c(0:4, "5-10", ">10"))

## After winsorizing at varL 30, we compute the total number of variant in
## each class for each variant size using dplyr.
bar.df = varmh.df %>% mutate(varL = winsor(varL, 20)) %>% group_by(varL, mhL.class) %>% 
    summarize(N = sum(N))
head(bar.df)
```

    ## # A tibble: 6 x 3
    ## # Groups:   varL [3]
    ##    varL mhL.class        N
    ##   <dbl> <fct>        <int>
    ## 1  1.00 0          6275487
    ## 2  1.00 1         10790647
    ## 3  2.00 0          1781047
    ## 4  2.00 1          1761704
    ## 5  2.00 2          3629809
    ## 6  3.00 0           769328

``` r
## Bar plot using ggplot2
ggplot(bar.df, aes(x = varL, y = N, fill = mhL.class)) + geom_bar(stat = "identity") + 
    theme_bw() + scale_fill_brewer(name = "mH length", palette = "Set1") + ylab("variants") + 
    xlab("variant size (bp)") + scale_x_continuous(breaks = 1:20, labels = c(1:19, 
    "20+")) + theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 
    1))
```

![](MHcut-Figures-dbSNPClinVar_files/figure-markdown_github/unnamed-chunk-4-1.png)

Wow, there are not that many variants with no micro-homology. I expected
much more.

Just a safety check, are we really looking at ~43M variants?

``` r
sum(varmh.df$N)
```

    ## [1] 43567375

Yep.

Of note, there is a one-liner way of doing this using the power of
pipes. It’s a long line but can be read from left to right:

``` r
varmh %>% as.data.frame %>% mutate(mhL.class = cut(mhL, breaks = c(-1, 0, 1, 
    2, 3, 4, 10, Inf), labels = c(0:4, "5-10", ">10")), varL = winsor(varL, 
    20)) %>% group_by(varL, mhL.class) %>% summarize(N = sum(N)) %>% ggplot(aes(x = varL, 
    y = N, fill = mhL.class)) + geom_bar(stat = "identity") + theme_bw() + scale_fill_brewer(name = "mH length", 
    palette = "Set1") + ylab("variants") + xlab("variant size (bp)") + scale_x_continuous(breaks = 1:20, 
    labels = c(1:19, "20+")) + theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 
    1))
```

![](MHcut-Figures-dbSNPClinVar_files/figure-markdown_github/unnamed-chunk-6-1.png)
