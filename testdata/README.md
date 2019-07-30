This folder contains a small dataset to test MHcut's installation.
It's a subset of the first 1 bp from chr 20 and 30 randomly selected deletions fomr dbSNP.

Note that the off-target estimation will not be correct in this test because the reference doesn't the rest of the genome.
It's just here for testing purpose.

```sh
jellyfish count --out-counter-len 1 -C -m 23 -s 100M chr20-1Mbp.fa
MHcut -var test-chr20-1Mbp.tsv -ref chr20-1Mbp.fa -jf mer_counts.jf -out test
```
