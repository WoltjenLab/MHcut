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
                    'name and ID the element position (if list).')
args = parser.parse_args()

# Parse 'info' argument
infos = []
for info in args.infos.split(','):
    info = info.split('|')
    infos.append([info[0], int(info[1])])

outf = open(args.outfile, 'w')
vcf_reader = vcf.Reader(open(args.vcffile, 'r'))
outline = ['chr\tstart\tend']
for info in infos:
    outline.append(info[0])
outf.write('\t'.join(outline) + '\n')
for record in vcf_reader:
    if(len(record.ALT) == 1):
        if(args.addchr):
            record.CHROM = 'chr' + record.CHROM
        if(args.rmchr):
            record.CHROM = record.CHROM.lstrip('chr')
        outline = [record.CHROM]
        outline.append(str(record.POS + 1))
        outline.append(str(record.POS + len(record.REF) - 1))
        for info in infos:
            if(info[0] not in record.INFO):
                outline.append('-')
            else:
                rec = record.INFO[info[0]]
                if(type(rec) is bool):
                   outline.append('T')
                else:
                   outline.append(rec[info[1]])
        outf.write('\t'.join(outline) + '\n')
outf.close()