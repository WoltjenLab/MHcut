# MHcut


## Microhomology search

    Flank 1:             |||x|||        |||x|||
	AGTGCCGTTAATCAGAGGTC-GGGCTGTGATGGTC-GGGGTGTTGTCGTTGACGTC
	Flank 2:        ||||           ||||

For each flank the microhomology is extended until the end of the variant as long as:

- first base is a match.
- no 2 consecutive mismatches.

The MH can be chosen from 2 flanks. MHcut uses a score to choose the "best" flank. 
The score is currently the number of matches + number of consecutive first matches.
In the example above, *Flank 1* is chosen (score: 9 vs 8).

## PAM cut search

                         |||x|||        |||x|||
	AGTGCCGTTAATCAGAGGTC-GGGCTGTGATGGTC-GGGGTGTTGTCGTTGACGTC
	                        <---------->

PAM cuts are searched between the end of the first exact match stretch of the MH and the variant boundary. 
In the example above, the first valid cut is between the G and C, and the last valid cut between the C and G.

PAM cuts are enumerated in both strands. For each valid cut, the protospacer sequence is retrieved.

Protospacers are blasted to the genome and the top 5 positions are parsed. *mm0* represents the number of position with full alignment and no mismatch. If *mm0* is equal to one, i.e. a unique match in the genome, the PAM cut is kept.

For each protospacer/cut, we also list other MHs that flank the cut and could be used preferentially instead of the one desired.
Only exact MHs of at least 3 bp are considered and if at least as close from each other as the target MH.
Among other, the output contains information about the best off-target MH (shortened to *bot*) defined as the off-target MH with the highest pattern score ([Bae et al 2014](http://www.nature.com.proxy3.library.mcgill.ca/articles/nmeth.3015)).

## Install

What you need:

- Blast. Executables available at [ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST).
- Python 2.7 or higher (but not Python 3)
- [pyfaidx](https://pypi.python.org/pypi/pyfaidx) module. Install with `pip install pyfaidx` or find alternatives [here](https://pypi.python.org/pypi/pyfaidx).
- Optional: JellyFish. Installation instructions on [the official webpage](http://www.genome.umd.edu/jellyfish.html).

## Preparing the reference genome

First download and unzip the reference genome, for example:

```shell
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

Eventually, you can index the genome using the following command:

```shell
python indexFasta.py hg38.fa
```

Otherwise this indexing will be done automatically the first time that MHcut is run (might take a few extra minutes).

By default MHcut uses BLAST to test if the protospacer sequence is unique in the genome. 
As a faster alternative, it can also use JellyFish to search for exact matches in the genome.

### Using BLAST

To build the BLAST database:

```shell
makeblastdb -in hg38.fa -dbtype nucl -title hg38
```

The files produced will be used automatically when providing the reference genome using `-ref` (see Usage below).

### Using JellyFish

JellyFish is a much faster alternative but can only find exact matches.
To count 23-mers in the reference genome:

```shell
jellyfish count -m 23 -s 100M hg38.fa
```

Note: Use `-t` to use multiple cores, e.g. `-t 10` to use 10 cores.

The output file `mer_counts.jf` will later be given to MHcut using `-jf` (see Usage below).

## Usage

	python MHcut.py -var NCBI_Variation_Viewer_data_uniq.tsv -ref hg38.fa -out MHcut-NCBI-chrX

The required parameters are:

- *-var* a TAB delimited file starting with chr/start/end columns and then whatever else (e.g. rsid, gene). The first row contains the column names.
- *-ref* a fasta file with the reference genome. It must be indexed and prepared for Blast (see previous section).
- *-out* the prefix for the output files (TSV files). 

Other optional parameters:

- *-minvarL* the minimum length for a variant to be considered. Default is `3`.
- *-minMHL* the minimum length of the MH. Default is `3`.
- *-maxConsMM*. The maximum number of consecutive mismatches allowed when extending the MH. Default is `1`.
- *-maxTail* the maximum distance betweem the MHs and a PAM cut to be considered valid. Currently default at `50`. Relevant for large variants.
- *-minhom* the minimum ratio of homology in the whole microhomology. Default is `0.8`.
- *-minm1L* the minimum length of the first stretch if the microhomology. Default is `3`.
- *-PAM* the PAM sequence. Default is `NGG`.
- *-PAMcut* the cut position relative to the PAM motif. Default is `-3`
- *-minMHLot* the minimum length of off-target MH to be considered in the off-target check. Default is `3`.
- *-nofilt* Don't filter variants without MH. all input variants will be present in the output *-variants* file. If used, the following parameters will NOT be taken into account: -minMHL, -minhom, -minm1L.
- *-jf* The 23-mers count file created by JellyFish. If provided, JellyFish will be used to test protospacer uniqueness instead of BLAST.

## Output

### The "variant" file

Named `PREFIX-variants.tsv`, the "variant" file  has one line per input variant with information about the MH found and if a valid PAM cut is available or not.

Currently the columns of the output are:

- The columns of the input `.tsv` file. For example: *chr*, *start*, *end*, *rsid*, *gene*.
- *mhL*: MH length.
- *mh1L*: number of first consecutive matches.
- *hom*: proportion of matches.
- *nbMM*: number of mismatches.
- *mhDist*: distance between the end of MH and the variant boundary.
- *MHseq1*/*MHseq2*: sequences of the MH.
- *pamMot*: the number of PAM motives in a valid location, no matter how unique the protospacer is.
- *pamUniq*: the number of PAM motives in a valid location and with unique the protospacer.
- *guidesNoOT*: the number of guides that have no off-target MH.
- *guidesMinOT*: the number of off-target MH for the guide which have the least amount of off-target MH.

### The "guide" file

Named `PREFIX-guides.tsv`, the "guide" file has one line per protospacer. 
This means that the same variant can be present several times if several valid PAM cuts are available. 
Only valid protospacers are returned, i.e. between micro-homologies and unique in the genome.

Currently the columns of the output are the same as for the "variant" file with the following additional columns:

- *protospacer* the sequence of the protospacer.
- *mm0* the number of position in the genome where the sequence align with no mismatches.
- *mm1* the number of position in the genome where the sequence align with 1 mismatches.
- *mm2* the number of position in the genome where the sequence align with 2 mismatches.
- *m1Dist1* and *m1Dist2*: the distance between the PAM cut the left or right stretch of perfect match, respectively.
- *mhDist1* and *mhDist2*: the distance between the PAM cut the left or right micro-homology, respectively.
- *nbOffTgt* the number of off-target MH.
- *largestOffTgt* the size of the largest off-target MH.
- *botScore* the MMEJ score of the **b**est **o**ff-**t**arget MH (best defined as the highest MMEJ score).
- *botSize* the MH length of the **b**est **o**ff-**t**arget MH (best defined as the highest MMEJ score).
- *botVarL* the length of the variant created by the **b**est **o**ff-**t**arget MH (best defined as the highest MMEJ score).
- *botGC* the GC content of the **b**est **o**ff-**t**arget MH (best defined as the highest MMEJ score).
- *botSeq* the sequence of the **b**est **o**ff-**t**arget MH (best defined as the highest MMEJ score).

### The "cartoon" file

Named `PREFIX-cartoons.tsv`, the "cartoon" file has one paragraph per variant with:

1. The corresponding line from the "variant" file (e.g. MH metrics, PAM found or not).
1. The location of the micro-homology. `|` means a match, `x` a mismatch.
1. The sequence of the flanks and variant. `-` to stress the variant's limits. Eventually `...` is the middle part of a large variant that is not shown for clarity purpose.
1. The location of valid PAM cuts. `\` and `/` depending on the strand of the PAM motif. `X` means there are cuts on each side of the base in opposite strand.
1. A list of valid protospacers if any.

For example, a 3 bases perfect MH with 3 valid PAM cuts:

```
chr8	41725834	41725853	ANK1	True	3	3	1.0	0	17	GCG	GCG
                     |||                  |||
GCGTGTCGTCGTTGCGGGCC-GCGATGTGCAGGGCCGGGAG-GCGCACCTTCCCCTTGGTGC
____________________ _______\___\\_______ ____________________
Protospacers:
GTTGCGGGCCGCGATGTGCA
CGGGCCGCGATGTGCAGGGC
GGGCCGCGATGTGCAGGGCC
```

## Next

- Add a *"two-cuts"* mode for large variants.
- Docker container to avoid installation nightmare.
- Tutorials: how to use a different PAM, how to use the Docker container, how was the dbSNP/ClinVar analysis for the paper ran.
