import argparse
from pyfaidx import Fasta
import subprocess
import os
import sys
# import re

Spy_PAM = "GG"
Spy_cut = -4
Spy_PAM_rev = "CC"


def mhTest(var_seq, fl_seq):
    '''Test for presence of microhomology between two sequences.'''
    res = {'score': 0, 'm1L': 0, 'mhL': 0, 'hom': 0, 'cartoon': '', 'seq1': '', 'seq2': ''}
    # ALignm the two sequences
    al_full = []
    for pos in range(min(len(var_seq), len(fl_seq))):
        al_full.append(var_seq[pos] == fl_seq[pos] and var_seq[pos] != 'N')
    # First base must match, otherwise return the 'res' as is
    if(not al_full[0]):
        return(res)
    # Trim the end of the alignment if two consecutive mismatches
    al_trimmed = al_full
    for pos in range(len(al_full)-1):
        if(not al_full[pos] and not al_full[pos+1]):
            al_trimmed = al_full[:pos]
            break
    # Cut potential last mismatch
    if(not al_trimmed[-1]):
        al_trimmed = al_trimmed[:-1]
    # Count consecutive matches in the beginning
    for pos in range(len(al_trimmed)):
        if(al_trimmed[pos]):
            res['m1L'] += 1
        else:
            break
    # Other alignment stats: length, homology, nb of mismatches
    nb_match = sum(al_trimmed)
    res['mhL'] = len(al_trimmed)
    res['hom'] = float(nb_match) / res['mhL']
    res['nbMM'] = res['mhL'] - nb_match
    res['mhdist'] = len(al_full) - res['mhL']  # this is the two hologous sequences
    # Compute a score, later used to choose which flank has the best MH
    res['score'] = res['m1L'] + nb_match
    # Cartoon of the MH (e.g. ||x|)
    res['cartoon'] = ''.join(['|' if al else 'x' for al in al_trimmed])
    # The two homologous sequences
    res['seq1'] = var_seq[0:res['mhL']]
    res['seq2'] = fl_seq[0:res['mhL']]
    return(res)


def findPAM(varseq, fl1seq, fl2seq, mhfl, maxTail):
    '''Look for PAM cuts between the MH regions.'''
    seq = fl1seq + varseq + fl2seq
    search_range = [len(fl1seq) - 1, len(fl1seq) + len(varseq) - mhfl['mhL']]
    if(mhfl['flank'] == 2):
        search_range = [len(fl1seq) + mhfl['mhL'] - 1, len(fl1seq) + len(varseq)]
    # Test each position: if it matched the motif and in the search range, add to list
    pams = []
    for pos in xrange(len(seq)-1):
        proto_seq = strand = cut_pos = False
        if(seq[pos:pos+2] == Spy_PAM):
            cut_pos = pos + Spy_cut - 1
            strand = '+'
        if(seq[pos:pos+2] == Spy_PAM_rev):
            cut_pos = pos - Spy_cut + len(Spy_PAM_rev) - 1
            strand = '-'
        if(strand and cut_pos >= search_range[0] and cut_pos < search_range[0] + maxTail
           and cut_pos < search_range[1] and cut_pos >= search_range[1] - maxTail):
            if(strand == '+'):
                proto_seq = seq[(cut_pos - 20 - Spy_cut):(cut_pos - Spy_cut)]
            else:
                proto_seq = seq[(cut_pos + Spy_cut + 2):(cut_pos + 20 + Spy_cut + 2)]
            pams.append({'cutPosition': cut_pos, 'strand': strand, 'proto': proto_seq})
    return(pams)


def alignPamsBlast(pams, reffile):
    '''Align protospacers and return an updated version of the input "pams" list.'''
    fasta_file = 'tempMHcut.fasta'
    ff = open(fasta_file, 'w')
    pams_hash = {}
    for pam in pams:
        pamid = str(pam['cutPosition']) + '_' + pam['strand']
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

# Define arguments
parser = argparse.ArgumentParser(description='Find regions with microhomology and a cut position.')
parser.add_argument('-var', dest='varfile', required=True,
                    help='the file with the variants location (BED-TSV format with header)')
parser.add_argument('-ref', dest='reffile', required=True, help='the reference genome fasta file')
parser.add_argument('-minvarL', dest='minvarL', default=3, help='the minimum variant length')
parser.add_argument('-minMHL', dest='minMHL', default=3, help='the minimum microhomology length')
parser.add_argument('-maxTail', dest='maxTail', default=50, help='the maximum hanging tail allowed')
parser.add_argument('-out', dest='outprefix', required=True, help='the prefix for the output files')
parser.add_argument('-minhom', dest='minhom', default=0.8, help='the minimum homology ratio')
parser.add_argument('-minm1L', dest='minm1L', default=2, help='the minimum length of first microhomology stretch')
args = parser.parse_args()

# Open connection to reference genome
reffa = Fasta(args.reffile)

# Open connection to output files
variant_output_file = open(args.outprefix + '-variants.tsv', 'w')
guide_output_file = open(args.outprefix + '-guides.tsv', 'w')
cartoon_output_file = open(args.outprefix + '-cartoons.tsv', 'w')

# Open connection to input file
variant_input_file = open(args.varfile, 'r')

# Read input file once to get number of line for progress bar
input_nb_lines = 0
for line in variant_input_file:
    input_nb_lines += 1
# Reopen
variant_input_file.close()
variant_input_file = open(args.varfile, 'r')

# Read firt line of input file (header) and write header in output file
# Change colunm names here.
# Add/remove columns here but also in the "Write in output files" section
inhead = variant_input_file.next().rstrip('\n')
outhead = inhead + '\tvarL\tpam\tmhL\tmh1L\thom\tnbMM\tmhDist\tMHseq1\tMHseq2\tpamMot'
variant_output_file.write(outhead + '\n')
gouthead = outhead + '\tprotospacer\tmm0\tmm1\tmm2\n'
guide_output_file.write(gouthead)
cartoon_output_file.write(outhead + '\n\n')

# Start progress bar
sys.stdout.write('[' + ' ' * 100 + ']')

# Read each line of the input file
line_cpt = 0
input_nb_lines = input_nb_lines / 100
for input_line in variant_input_file:
    line_cpt += 1
    if(line_cpt % input_nb_lines == 0):
        percent = line_cpt / input_nb_lines
        sys.stdout.write('\r[' + '*' * percent + ' ' * (100 - percent) + ']')
    input_line_raw = input_line.rstrip('\n')
    input_line = input_line_raw.split('\t')
    vstart = int(input_line[1])
    vend = int(input_line[2])
    vsize = vend - vstart + 1
    if(vsize < args.minvarL):
        # If variant is too small or too big, skip and jump to next iteration
        continue
    if(input_line[0] not in reffa.keys()):
        # If chromosome name not in reference, skip variant
        continue
    # Retrieve sequences. Careful with position and shift
    flsize = max(vsize, 20)
    varseq = str(reffa[input_line[0]][(vstart-1):vend])
    fl1seq = str(reffa[input_line[0]][(vstart-flsize-1):(vstart-1)])
    fl2seq = str(reffa[input_line[0]][vend:(vend+flsize)])
    # Test MH in each flank (reverse for flank 1) and save best MH
    mhfl1 = mhTest(varseq[::-1], fl1seq[::-1])
    mhfl2 = mhTest(varseq, fl2seq)
    if(mhfl1['score'] > mhfl2['score']):  # Using the alignment score, the best flank is chosen
        mhfl = mhfl1
        mhfl['flank'] = 1  # This is to remember which flank is chosen
        mhfl['cartoon'] = mhfl['cartoon'][::-1]
    else:
        mhfl = mhfl2
        mhfl['flank'] = 2
    # If no MH or too small, or too low MH ratio or too short first microhomology stretch, jump to the next input line
    if(mhfl['score'] == 0 or mhfl['mhL'] < args.minMHL or mhfl['hom'] < args.minhom or mhfl['m1L'] < args.minm1L):
        continue
    # Find PAM motives
    pams = findPAM(varseq, fl1seq, fl2seq, mhfl, args.maxTail)
    # Map protospacers to the genome and keep unique ones
    nb_pam_motives = len(pams)
    if(nb_pam_motives > 0):
        pams = alignPamsBlast(pams, args.reffile)
        pams_filter = []
        for pam in pams:
            # This where to define how unique the protospacer must be
            # Here there must be only one position in the genome aligning perfectly
            if pam['mm0'] == 1:
                pams_filter.append(pam)
        pams = pams_filter
    # Write in output files
    # Add/remove columns here (without forgetting the header)
    voutline = input_line_raw + '\t' + str(vsize) + '\t' + str(len(pams) > 0)
    voutline += '\t' + str(mhfl['mhL']) + '\t'
    voutline += str(mhfl['m1L']) + '\t' + str(round(mhfl['hom'], 2)) + '\t' + str(mhfl['nbMM'])
    voutline += '\t' + str(mhfl['mhdist']) + '\t' + mhfl['seq1'] + '\t' + mhfl['seq2'] + '\t' + str(nb_pam_motives)
    variant_output_file.write(voutline + '\n')
    for pam in pams:
        guide_output_file.write(voutline + '\t' + pam['proto'] + '\t' + str(pam['mm0']))
        guide_output_file.write('\t' + str(pam['mm1']) + '\t' + str(pam['mm2']) + '\n')
    # Write the cartoon
    cartoon_output_file.write(voutline + '\n')
    cartoon_output_lines = ['', '', '']
    # Cartoon: alignment line
    white_spaces_before = len(fl1seq) - mhfl['mhL']
    if(mhfl['flank'] == 2):
        white_spaces_before += mhfl['mhL'] + 1
    white_spaces_before = ' ' * white_spaces_before
    white_spaces_between = ' ' * (vsize - mhfl['mhL'] + 1)
    cartoon_output_lines[0] = white_spaces_before + mhfl['cartoon']
    cartoon_output_lines[0] += white_spaces_between + mhfl['cartoon']
    # Cartoon: sequence line
    cartoon_output_lines[1] = fl1seq + '-' + varseq + '-' + fl2seq
    # Cartoon: PAM position line
    pam_cartoon = ['_' for i in range(len(fl1seq) + len(fl2seq) + vsize)]
    for pam in pams:
        if(pam['strand'] == '+'):
            pam_cartoon[pam['cutPosition'] + 1] = '\\'
        else:
            pam_cartoon[pam['cutPosition']] = '/'
    pam_cartoon = ''.join(pam_cartoon)
    cartoon_output_lines[2] = pam_cartoon[:len(fl1seq)] + ' '
    cartoon_output_lines[2] += pam_cartoon[len(fl1seq):(len(fl1seq) + vsize)]
    cartoon_output_lines[2] += ' ' + pam_cartoon[(len(fl1seq) + vsize):]
    # If the line is too long (large variants), trim the ends and middle
    flank_buffer = 10  # How many extra bases to show on the flanks (outside of the MH region)
    if(len(cartoon_output_lines[1]) > 2 * (flank_buffer + mhfl['mhL'] + args.maxTail)):
        cartoon_part1 = [len(fl1seq) - mhfl['mhL'] - flank_buffer, len(fl1seq) + args.maxTail]
        cartoon_part2 = [len(fl1seq) + vsize - mhfl['mhL'] - args.maxTail, len(fl1seq) + vsize + flank_buffer]
        # If second flank was used, shift the positions
        if(mhfl['flank'] == 2):
            cartoon_part1 = [p + mhfl['mhL'] + 1 for p in cartoon_part1]
            cartoon_part2 = [p + mhfl['mhL'] + 1 for p in cartoon_part2]
        # Sanity checks that the position are in the correct range
        cartoon_part1[0] = max(cartoon_part1[0], 0)
        # Update lines (if the two part overlaps merge them into one)
        cartoon_output_lines_trimmed = []
        if(cartoon_part1[1] + 1 > cartoon_part2[0]):
            cartoon_part1[1] = cartoon_part2[1]
            for line in cartoon_output_lines:
                line = line[cartoon_part1[0]:min(cartoon_part1[1], len(line))]
                cartoon_output_lines_trimmed.append(line)
        else:
            for line in cartoon_output_lines:
                full_line = line
                line = full_line[cartoon_part1[0]:cartoon_part1[1]]
                line += '...'
                line += full_line[cartoon_part2[0]:min(cartoon_part2[1], len(full_line))]
                cartoon_output_lines_trimmed.append(line)
        cartoon_output_lines = cartoon_output_lines_trimmed
    # Write cartoon lines
    for line in cartoon_output_lines:
        cartoon_output_file.write(line + '\n')
    # Cartoon: protospacers sequence
    if(len(pams) > 0):
        cartoon_output_file.write('Protospacers:\n')
        for pam in pams:
            cartoon_output_file.write(pam['proto'] + '\n')
    cartoon_output_file.write('\n\n')

print '\nDone.\n'

variant_input_file.close()
variant_output_file.close()
guide_output_file.close()
cartoon_output_file.close()
