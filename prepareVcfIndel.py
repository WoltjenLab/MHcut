import vcf
import argparse


# Define arguments
parser = argparse.ArgumentParser(description='Convert VCF into TSV. ' +
                                 'Example: python prepareVcfIndel.py -vcf var.vcf' +
                                 ' -o var.tsv -info "RS|0,CAF|1,GENEINFO|0"')
parser.add_argument('-vcf', dest='vcffile', required=True,
                    help='the VCF file with the variants location (BED-TSV format with header)')
parser.add_argument('-o', dest='outfile', required=True,
                    help='the output BED-TSV format with header')
parser.add_argument('-addchr', dest='addchr', action='store_true',
                    help='Add "chr" prefix to chrom.')
parser.add_argument('-rmchr', dest='rmchr', action='store_true',
                    help='Remove "chr" prefix from chrom.')
parser.add_argument('-info', dest='infos', default='',
                    help='the INFO fields to retrieve in the form ' +
                    '"FIELD|ID,FIELD|ID...". FIELD being the field ' +
                    'name and ID the element position (if list). ' +
                    'If ID is "-", all the elements are taken')
parser.add_argument('-infomerge', dest='infosmerge', default='',
                    help='the boolean INFO fields to retrieve and merge into a columns' +
                    '"COL|FIELD,FIELD...". FIELD being the field ' +
                    'name and COL the name of the new column.')
args = parser.parse_args()

# Parse 'info' argument
infos = []
for info in args.infos.split(','):
    info = info.split('|')
    if(info[1] == '-'):
        infos.append([info[0], -1])
    else:
        infos.append([info[0], int(info[1])])
if(args.infosmerge != ''):
    fields_merge = args.infosmerge.split('|')[1].split(',')
    merge_col = args.infosmerge.split('|')[0]

# Open input and output files
vcf_reader = vcf.Reader(open(args.vcffile, 'r'))
outf = open(args.outfile, 'w')
# Write header of the output file
outline = ['chr\tstart\tend']
for info in infos:
    outline.append(info[0])
if(args.infosmerge != ''):
    outline.append(merge_col)
outf.write('\t'.join(outline) + '\n')
# Read input VCF one line at a time and write the corresponding output line.
for record in vcf_reader:
    # Filter out duplications (where the alternative allele is several bp long)
    if(len(record.ALT) == 1):
        # Add/remove 'chr' if asked by '-addchr' or '-rmchr'
        if(args.addchr):
            record.CHROM = 'chr' + record.CHROM
        if(args.rmchr):
            record.CHROM = record.CHROM.lstrip('chr')
        outline = [record.CHROM]
        # Calculate the deletion coordinates
        outline.append(str(record.POS + 1))
        outline.append(str(record.POS + len(record.REF) - 1))
        # Look for each of the desired info fields
        for info in infos:
            if(info[0] not in record.INFO):
                # If this information is missing
                outline.append('-')
            else:
                rec = record.INFO[info[0]]
                if(type(rec) is bool):
                    # If a boolean output 'T', otherwise '-' (see above)
                    rec = 'T'
                else:
                    if(type(rec) == list):
                        # If list either get the relevant element...
                        if(info[1] >= 0):
                            rec = rec[info[1]]
                        else:
                            # ...or join all of them together
                            rec = ';'.join(rec)
                if(rec is None):
                    # Edge case: sometimes info is there but with no value
                    outline.append('-')
                else:
                    # Otherwise add the value to the output line
                    outline.append(str(rec))
        # Look for the (boolean) info fields that should be merge into one column
        if(args.infosmerge != ''):
            rec = []
            for field in fields_merge:
                if(field in record.INFO):
                    rec.append(field)
            if(len(rec) == 0):
                outline.append('-')
            else:
                outline.append(','.join(rec))
        # Write the output line to the output file
        outf.write('\t'.join(outline) + '\n')
outf.close()
