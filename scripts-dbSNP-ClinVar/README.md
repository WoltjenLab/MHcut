Scripts used to prepare the dbSNP and ClinVar databases, and run MHcut.
See commands in `trace.sh`.

Briefly:

1. The input TSV file is prepared by merging dbSNP and ClinVar deletions.
1. The reference genome and JellyFish index are prepared.
1. MHcut is run in parallel on (600) chunks of the input TSV using Slurm's Job Arrays.
1. Eventually, jobs that died or hit the walltime are restarted.
1. Finally the chunked outputs are merged.

The only difference between the main results and the *xCas9* results are the PAMs used (`-PAM NGN,GAA,GAT` in [mhcutJob-xCas9.sh](mhcutJob-xCas9.sh)).

For the benchmarking experiment, 100K variants were randomly selected and MHcut was run with different parameters.
