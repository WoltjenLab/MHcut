dbSNP/ClinVar exploration
=========================

`GENEINFO` columns the same?
----------------------------

There is a *GENEINFO* column that comes from the dbSNP VCF and another one from the ClinVar. Are they the same? Or should we really keep both?

    ## [1] 0.8838559

|     | GENEINFO                    | GENEINFO.ClinVar      |
|-----|:----------------------------|:----------------------|
| 3   | PLCH2:9651|PEX10:5192       | PEX10:5192|PLCH2:9651 |
| 10  | MIR6808:102466740|DVL1:1855 | DVL1:1855             |
| 13  | B3GALT6:126792|SDF4:51150   | B3GALT6:126792        |
| 15  | B3GALT6:126792|SDF4:51150   | B3GALT6:126792        |
| 16  | B3GALT6:126792|SDF4:51150   | B3GALT6:126792        |
| 17  | B3GALT6:126792|SDF4:51150   | B3GALT6:126792        |

The columns are the same only 88% of the time. Often a gene is "missing" in the ClinVar column. Or is it because the ClinVar column only shows the disease-relevant gene?

Consistent frequencies?
-----------------------

Are the frequency from the 1000GP, ExAC and TOPMED consistent?

![](dbSNPClinVarExploration_files/figure-markdown_github/unnamed-chunk-3-1.png)![](dbSNPClinVarExploration_files/figure-markdown_github/unnamed-chunk-3-2.png)![](dbSNPClinVarExploration_files/figure-markdown_github/unnamed-chunk-3-3.png)

Quite consistent.

No frequency means rare?
------------------------

When the frequency information is missing (e.g. 1000GP frequency), is the variant rare?

We would expect so: there is no information because these variants were not seen in the study because they are rare. Let's check that the frequency in the larger studies (TOPMED, ExAC) is on the rare side.

![](dbSNPClinVarExploration_files/figure-markdown_github/unnamed-chunk-4-1.png)![](dbSNPClinVarExploration_files/figure-markdown_github/unnamed-chunk-4-2.png)

Yep.

Citations
---------

The *citation* column is specific to ClinVar variants.

Do most variants in ClinVar have associated citations? How many non-pathogenic variants have a citation associated with them?

| CLNSIG                                          |  nb.var|  prop.wth.citation|
|:------------------------------------------------|-------:|------------------:|
| Pathogenic                                      |   10550|          0.5813270|
| Uncertain\_significance                         |    4812|          0.2711970|
| Likely\_pathogenic                              |    2430|          0.4164609|
| Likely\_benign                                  |    1932|          0.3048654|
| Benign                                          |    1340|          0.4000000|
| not\_provided                                   |    1225|          0.2318367|
| Conflicting\_interpretations\_of\_pathogenicity |     455|          0.8131868|
| Pathogenic/Likely\_pathogenic                   |     436|          0.8853211|
| Benign/Likely\_benign                           |     266|          0.7857143|
| -                                               |      50|          0.8200000|

Around 60% of pathogenic variants (the largest class of variants in ClinVar) have at least one associated publication.

Duplicates
----------

Any variants annotated several times? (Defining a variant in term of its location/coordinates)

![](dbSNPClinVarExploration_files/figure-markdown_github/dup-1.png)

Most of them are fine but a few are annotated twice. Some are even annotated 10 times !?

Variant size
------------

![](dbSNPClinVarExploration_files/figure-markdown_github/unnamed-chunk-6-1.png)![](dbSNPClinVarExploration_files/figure-markdown_github/unnamed-chunk-6-2.png)![](dbSNPClinVarExploration_files/figure-markdown_github/unnamed-chunk-6-3.png)![](dbSNPClinVarExploration_files/figure-markdown_github/unnamed-chunk-6-4.png)
