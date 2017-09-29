# MHcut


## Microhomology search

    Flank 1:             |||x|||        |||x|||
	AGTGCCGTTAATCAGAGGTC-GGGCTGTGATGGTC-GGGGTGTTGTCGTTGACGTC
	Flank 2:        ||||           ||||

For each flank the microhomology is extended until the end of the variant as long as:

- first base is a match.
- no 2 consecutive mismatches.

The MH can be chosen from 2 flanks. MHcut uses a score to choose the "best" flank. The score is the number of matches + number of consecutive first matches.

## PAM cut search

                         |||x|||        |||x|||
	AGTGCCGTTAATCAGAGGTC-GGGCTGTGATGGTC-GGGGTGTTGTCGTTGACGTC
	                            <------>

PAM cuts are searched between the MH and the variant boundary. In the example above, the first valid cut is between the T and G, and the last valid cut between the C and G.

PAM cuts are enumerated in both strands. For each valid cut, the protospacer sequence is retrieved.

Protospacers are blasted to the genome and the top 5 positions are parsed. *mm0* represents the number of position with full alignment and no mismatch. If *mm0* is equal to one, i.e. a unique match in the genome, the PAM cut is kept.

## Install

What you need:

- Blast. Executables available at [ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST).
- Python 2.7 or higher (but not Python 3)
- [pyfaidx](https://pypi.python.org/pypi/pyfaidx) module. Install with `pip install pyfaidx` or find alternatives [here](https://pypi.python.org/pypi/pyfaidx).

## Preparing the reference genome

The following commands download the reference genome, index it and build a Blast database.

```shell
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
python indexFasta.py hg38.fa
makeblastdb -in hg38.fa -dbtype nucl -title hg38
```

## Usage

	python MHcut.py -var NCBI_Variation_Viewer_data_uniq.tsv -ref hg38.fa -out MHcut-NCBI-chrX

The required parameters are:

- *-var* a TAB delimited file starting with chr/start/end columns and then whatever else (e.g. rsid, gene). The first row contains the column names.
- *-ref* a fasta file with the reference genome. It must be indexed and prepared for Blast (see previous section).
- *-out* the prefix for the output files (TSV files). 

Other optional parameters:

- *-minvarL* the minimum length for a variant to be considered. Default is 3.
- *-minMHL* the minimum length of the MH. Default is 3.
- *-maxTail* the maximum distance betweem the MHs and a PAM cut to be considered valid. Currently default at 50. Relevant for large variants.
- *-minhom* the minimum ratio of homology in the whole microhomology. Default is 0.8.
- *-minm1L* the minimum length of the first stretch if the microhomology. Default is 2.

## Output

### The "variant" file

Named `PREFIX-variants.tsv`, the "variant" file  has one line per input variant with information about the MH found and if a valid PAM cut is available or not.

Currently the columns of the output are:

- The columns of the input `.tsv` file. For example: *chr*, *start*, *end*, *rsid*, *gene*.
- *pam*: True is there is a valid PAM cut with unique protospacer. False otherwise.
- *mhL*: MH length.
- *mh1L*: number of first consecutive matches.
- *hom*: proportion of matches.
- *nbMM*: number of mismatches.
- *mhDist*: distance between the end of MH and the variant boundary.
- *MHseq1*/*MHseq2*: sequences of the MH.
- *pamMot*: the number of PAM motives in a valid location, no matter how unique the protospacer is.

### The "guide" file

Named `PREFIX-guides.tsv`, the "guide" file has one line per protospacer. It means the same variant can be present several times if several valid PAM cuts are available. Only valid protospacers are returned, i.e. between micro-homologies and unique in the genome.

Currently the columns of the output are the same as for the "variant" file with the following additional columns:

- *protospacer* the sequence of the protospacer.
- *mm0* the number of position in the genome where the sequence align with no mismatches.
- *mm1* the number of position in the genome where the sequence align with 1 mismatches.
- *mm0* the number of position in the genome where the sequence align with 2 mismatches.

### The "cartoon" file

Named `PREFIX-cartoons.tsv`, the "cartoon" file has one paragraph per variant with:

1. The corresponding line from the "variant" file (e.g. MH metrics, PAM found or not).
1. The location of the micro-homology. `|` means a match, `x` a mismatch.
1. The sequence of the flanks and variant. `-` to stress the variant's limits. Eventually `...` is the middle part of a large variant that is not shown for clarity purpose.
1. The location of valid PAM cuts. `\` and `/` depending on the strand of the PAM motif.
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
