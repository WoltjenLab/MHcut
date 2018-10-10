import argparse

# Define arguments
parser = argparse.ArgumentParser(description='Prepare ClinVar information. ' +
                                 'Example: python prepareClinVar.py -vcf clinvar_20180701.tsv' +
                                 '-cit var_citations.txt -o clinvar-grch38-deletion.tsv')
parser.add_argument('-tsv', dest='tsvfile', required=True,
                    help='the tab-delimited prod>uced by prepareVcfIndel.py on ClinVar VCF.')
parser.add_argument('-cit', dest='citfile', default='',
                    help='the tab-delimited variant citation file (*.txt)')
parser.add_argument('-o', dest='outfile', required=True,
                    help='the output file')
args = parser.parse_args()

# Save citation information (AlleleID -> unique citations as a string)
cit = {}
for line in open(args.citfile):
    line = line.rstrip().split('\t')
    citsrc = 'PM'
    if(line[4] != 'PubMed'):
        citsrc = 'PMC'
        if(line[4] == 'NCBIBookShelf'):
            citsrc = ''
    citid = citsrc + line[5]
    if(line[0] in cit):
        if(citid not in cit[line[0]]):
            cit[line[0]].append(citid)
    else:
        cit[line[0]] = [citid]
for al in cit:
    cit[al] = ';'.join(cit[al])

# Open output file, write headers
# and figure out index of relevant columns
outf = open(args.outfile, 'w')
tsvf = open(args.tsvfile)
header = tsvf.next()
header = header.rstrip().split('\t')
typecol = header.index('CLNVC')
idcol = header.index('ALLELEID')
header.pop(typecol)
header.append('citation')
outf.write('\t'.join(header) + '\n')

# Read TSV file line by line, keep deletions,
# add citation column and write to output file
for line in tsvf:
    line = line.rstrip().split('\t')
    vtype = line.pop(typecol)
    if(vtype == 'Deletion'):
        if(line[idcol] in cit):
            line.append(cit[line[idcol]])
        else:
            line.append('-')
        outf.write('\t'.join(line) + '\n')
outf.close()
