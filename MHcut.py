import argparse
from pyfaidx import Fasta
import subprocess
import os
# import re

Spy_PAM = "GG"
Spy_cut = -4
Spy_PAM_rev = "CC"


def mhTest(seq1, seq2, debug=False):
    '''Test for presence of microhomology between two sequences.'''
    res = {'score': 0, 'm1L': 0, 'mhL': 0, 'hom': 0, 'al': '', 'seq1': '', 'seq2': ''}
    # ALignment
    al = []
    for pos in range(len(seq1)):
        al.append(seq1[pos] == seq2[pos])
    # First base must match
    if(not al[0]):
        return(res)
    # Cut alignment if two consecutive mismatches
    alcut = al
    for pos in range(len(al)-1):
        if(not al[pos] and not al[pos+1]):
            alcut = al[:pos]
            break
    # Cut potential last mismatch
    if(not alcut[-1]):
        alcut = alcut[:-1]
    # Count consecutive matches in the beginning
    m1L = 0
    for pos in range(len(alcut)):
        if(alcut[pos]):
            m1L += 1
        else:
            break
    # Other alignment stats
    nb_match = sum(alcut)
    hom = float(nb_match) / len(alcut)
    # Compute a score
    score = m1L + nb_match
    # Print debug message
    alverb = ''.join(['|' if all else 'x' for all in alcut])
    res = {'score': score, 'm1L': m1L, 'mhL': len(alcut), 'hom': hom, 'al': alverb,
           'nbMM': len(alcut) - nb_match,
           'mhdist': len(al) - len(alcut),
           'seq1': seq1[0:len(alcut)], 'seq2': seq2[0:len(alcut)]}
    return(res)


def alignPamsBlast(pams, reffile):
    '''Align protospacers and return an updated 'pams' object.'''
    fasta_file = 'tempMHcut.fasta'
    ff = open(fasta_file, 'w')
    pams_hash = {}
    for pam in pams:
        pamid = str(pam['cut']) + '_' + pam['strand']
        pam['mm0'] = 0
        pam['mm1'] = 0
        pam['mm2'] = 0
        pams_hash[pamid] = pam
        ff.write('>' + pamid + '\n' + pam['proto'] + '\n')
    ff.close()
    dump = open('/dev/null')
    blast_cmd = ['blastn', '-db', reffile,  '-query', fasta_file, '-outfmt', '6',
                 '-word_size', '10', '-max_target_seqs', '20']
    blast_out = subprocess.check_output(blast_cmd, stderr=dump)
    dump.close()
    blast_out = blast_out.split('\n')
    for line in blast_out:
        line = line.split('\t')
        if(len(line) > 1 and int(line[3]) == len(pams_hash[line[0]]['proto'])):
            if(int(line[4]) == 0):
                pams_hash[line[0]]['mm0'] += 1
            if(int(line[4]) == 1):
                pams_hash[line[0]]['mm1'] += 1
            if(int(line[4]) == 2):
                pams_hash[line[0]]['mm2'] += 1
    os.remove(fasta_file)
    return pams


# def alignPamsBWA(pams, reffile):
#     '''Align protospacers and return an updated 'pams' object.'''
#     fastq_file = 'tempMHcut.fastq'
#     ff = open(fastq_file, 'w')
#     for pam in pams:
#         ff.write('@' + str(pam['cut']) + '_' + pam['strand'] + '\n' +
#                  pam['proto'] + '\n+\n' + '~'*len(pam['proto']) + '\n')
#     ff.close()
#     sai_file = 'tempMHcut.sai'
#     dump = open('/dev/null')
#     aln_cmd = ['bwa', 'aln', reffile, fastq_file]
#     aln_out = subprocess.check_output(aln_cmd, stderr=dump)
#     ff = open(sai_file, 'w')
#     ff.write(aln_out)
#     ff.close()
#     samse_cmd = ['bwa', 'samse', reffile, sai_file, fastq_file]
#     samse_out = subprocess.check_output(samse_cmd, stderr=dump)
#     dump.close()
#     samse_out = samse_out.split('\n')
#     map = {}
#     for line in samse_out:
#         if(line[0] != '@'):
#             mism = re.search('NM:i:(\d*)', line)
#             best = re.search('X0:i:(\d*)', line)
#             line = line.split('\t')
#             map[line[0]] = {'mm': mism.group(1), 'best': best.group(1)}
#     for pam in pams:
#         pamap = map[str(pam['cut']) + '_' + pam['strand']]
#         pam['mm'] = pamap['mm']
#         pam['best'] = pamap['best']
#     # os.remove(fastq_file)
#     # os.remove(sai_file)
#     return pams


def findPAM(seq, flank_size):
    '''Find PAM motives in a sequence.'''
    res = []
    for pos in xrange(len(seq)-1):
        cut = proto_seq = strand = real_cut = False
        if(seq[pos:pos+2] == Spy_PAM):
            cut = pos + Spy_cut
            if(cut > flank_size - 1 and cut < len(seq) - flank_size + 1):
                real_cut = cut - 1
                proto_seq = seq[(cut - 20 - Spy_cut - 1):(cut - Spy_cut - 1)]  # TODO
                strand = '+'
        if(seq[pos:pos+2] == Spy_PAM_rev):
            cut = pos - Spy_cut + len(Spy_PAM_rev) - 1
            if(cut > flank_size - 2 and cut < len(seq) - flank_size):
                real_cut = cut
                proto_seq = seq[(cut + Spy_cut + 2):(cut + 20 + Spy_cut + 2)]  # TODO
                strand = '-'
        if(strand):
            res.append({'cut': cut, 'realcut': real_cut, 'strand': strand,
                        'proto': proto_seq})
    return(res)

# Define arguments
parser = argparse.ArgumentParser(description='Find regions with microhomology and a cut position.')
parser.add_argument('-var', dest='varfile', required=True,
                    help='the file with the variants location (BED-TSV format with header)')
parser.add_argument('-ref', dest='reffile', required=True, help='the reference genome fasta file')
parser.add_argument('-minvarL', dest='minvarL', default=4, help='the minimum variant length')
parser.add_argument('-maxvarL', dest='maxvarL', default=50, help='the maximum variant length')
parser.add_argument('-minMHL', dest='minMHL', default=3, help='the minimum microhomology length')
parser.add_argument('-maxTail', dest='maxTail', default=20, help='the maximum hanging tail allowed')
parser.add_argument('-debug', dest='debug', action='store_true', help='debug mode')
parser.add_argument('-vout', dest='voutfile', required=True, help='the variant output file')
parser.add_argument('-gout', dest='goutfile', required=True, help='the guide output file')
args = parser.parse_args()

# Open connection to reference genome
reffa = Fasta(args.reffile)

# Read each line of the variant file and analyze it and write in output file
var_lines = open(args.varfile, 'r')
inhead = var_lines.next().rstrip('\n')
outf = open(args.voutfile, 'w')
outhead = inhead + '\tpam\tmhL\tmh1L\thom\tnbMM\tmhDist\tseq1\tseq2'
outf.write(outhead + '\n')
goutf = open(args.goutfile, 'w')
gouthead = outhead + '\tseq\tmm0\tmm1\tmm2\n'
goutf.write(gouthead)
if(args.debug):
    print outhead
for var_line in var_lines:
    var_line_raw = var_line.rstrip('\n')
    var_line = var_line_raw.split('\t')
    vstart = int(var_line[1])
    vend = int(var_line[2])
    vsize = vend - vstart + 1
    if(vsize < args.minvarL or vsize > args.maxvarL):
        # If variant is too small or too big, skip and jump to next iteration
        continue
    if(var_line[0] not in reffa.keys()):
        # If chromosome name not in reference, skip variant
        continue
    # print var_line_raw
    # Retrieve sequences. Careful with position and shift
    flsize = max(vsize, 20)
    varseq = str(reffa[var_line[0]][(vstart-1):vend])
    fl1seq = str(reffa[var_line[0]][(vstart-flsize-1):(vstart-1)])
    fl2seq = str(reffa[var_line[0]][vend:(vend+flsize)])
    # Test MH in each flank (reverse for flank 1) and save best MH
    mhfl1 = mhTest(varseq[::-1], fl1seq[::-1])
    mhfl2 = mhTest(varseq, fl2seq)
    fl = 1
    if(mhfl1['score'] > mhfl2['score']):
        mhfl = mhfl1
        mhfl['al'] = mhfl['al'][::-1]
    else:
        mhfl = mhfl2
        fl = 2
    # If not MH, skip
    if(mhfl['score'] == 0 or mhfl['mhL'] < args.minMHL):
        continue
    # Find PAM motives
    allseq = fl1seq + varseq + fl2seq
    pams = findPAM(allseq, len(fl1seq))
    # Find PAM in correct location and filter others
    goodMin = flsize - 1
    goodMax = flsize + vsize - mhfl['mhL']
    if(fl == 2):
        goodMin = flsize + mhfl['mhL'] - 1
        goodMax = flsize + vsize
    for pam in pams:
        if(pam['realcut'] >= goodMin and pam['realcut'] < goodMin + args.maxTail
           and pam['realcut'] < goodMax and pam['realcut'] >= goodMax - args.maxTail):
            pam['goodpam'] = True
        else:
            pam['goodpam'] = False
    pams_filter = []
    for pam in pams:
        if pam['goodpam']:
            pams_filter.append(pam)
    pams = pams_filter
    # Align protospacers and keep unique ones
    if(len(pams) > 0):
        pams = alignPamsBlast(pams, args.reffile)
        pams_filter = []
        for pam in pams:
            if pam['mm0'] == 1:
                pams_filter.append(pam)
        pams = pams_filter
    # Write in output files
    outline = var_line_raw + '\t' + str(len(pams) > 0) + '\t' + str(mhfl['mhL']) + '\t'
    outline += str(mhfl['m1L']) + '\t' + str(round(mhfl['hom'], 2)) + '\t' + str(mhfl['nbMM'])
    outline += '\t' + str(mhfl['mhdist']) + '\t' + mhfl['seq1'] + '\t' + mhfl['seq2']
    outf.write(outline + '\n')
    for pam in pams:
        goutline = outline + '\t' + pam['proto'] + '\t' + str(pam['mm0'])
        goutline += '\t' + str(pam['mm1']) + '\t' + str(pam['mm2']) + '\n'
        goutf.write(goutline)
    # Debug 'cartoon'
    if(args.debug):
        print outline
        spaces = flsize + 1
        space1 = flsize - mhfl['mhL']
        if(fl == 1):
            spaces = space1
        spaces = ' ' * spaces
        print spaces + mhfl['al'] + ' ' * (vsize - mhfl['mhL'] + 1) + mhfl['al']
        print fl1seq + '-' + varseq + '-' + fl2seq
        pamdebug = ['_' for i in range(2*flsize + vsize)]
        for pam in pams:
            if(pam['strand'] == '+'):
                pamdebug[pam['cut']] = '\\'
            else:
                pamdebug[pam['cut']] = '/'
        pamdebug = ''.join(pamdebug)
        pamdebug = pamdebug[:flsize] + ' ' + pamdebug[flsize:(flsize + vsize)] + ' ' + pamdebug[(flsize + vsize):]
        print pamdebug
        if(len(pams) > 0):
            print 'Protospacers:'
            for pam in pams:
                print pam['proto']
        print '\n'

var_lines.close()
outf.close()
goutf.close()
