from pyfaidx import Fasta
import subprocess
import os
import sys
from math import exp
# import re

dnacomp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'W': 'W', 'S': 'S',
           'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R', 'V': 'B', 'B': 'V',
           'H': 'D', 'D': 'H', 'N': 'N'}


def GC(seq):
    seq = list(seq)
    gc_count = seq.count('C') + seq.count('G')
    gc_prop = float(gc_count) / len(seq)
    return round(gc_prop, 3)


def mhTest(var_seq, fl_seq, maxConsMM=1):
    '''Test for presence of microhomology between two sequences.'''
    res = {'score': 0, 'm1L': 0, 'mhL': 0, 'hom': 0, 'nbMM': 0, 'cartoon': '',
           'seq1': '', 'seq2': ''}
    # ALignm the two sequences
    al_full = []
    for pos in range(min(len(var_seq), len(fl_seq))):
        al_full.append(var_seq[pos] == fl_seq[pos] and var_seq[pos] != 'N')
    # First base must match, otherwise return the 'res' as is
    if(not al_full[0]):
        return(res)
    # Trim the end of the alignment if X consecutive mismatches
    al_trimmed = al_full
    consMM = 0
    for pos in range(len(al_full)-1):
        if(not al_full[pos]):
            consMM += 1
        else:
            consMM = 0
        if(consMM > maxConsMM):
            al_trimmed = al_full[:(pos-consMM+1)]
            break
    # Cut potential last mismatch
    while(not al_trimmed[-1]):
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
    res['mhdist'] = len(al_full) - res['mhL']  # dist btw homologous sequences
    # Compute a score, later used to choose which flank has the best MH
    res['score'] = res['m1L'] + nb_match
    # Cartoon of the MH (e.g. ||x|)
    res['cartoon'] = ''.join(['|' if al else 'x' for al in al_trimmed])
    # The two homologous sequences
    res['seq1'] = var_seq[0:res['mhL']]
    res['seq2'] = fl_seq[0:res['mhL']]
    # Maximum GC content of the homologous sequences
    gc1 = GC(res['seq1'])
    gc2 = GC(res['seq2'])
    res['gc'] = max(gc1, gc2)
    return(res)


def aligned(seq, pamseq):
    '''Is there a match between seq and pamseq?
    Takes NWSMKRYVBHD into account.'''
    align = True
    pos = 0
    while(align and pos < len(pamseq)):
        if(pamseq[pos] != 'N'):
            align = pamseq[pos] == seq[pos]
        if(pamseq[pos] == 'W'):
            align = seq[pos] == 'A' or seq[pos] == 'T'
        elif(pamseq[pos] == 'S'):
            align = seq[pos] == 'G' or seq[pos] == 'C'
        elif(pamseq[pos] == 'M'):
            align = seq[pos] == 'A' or seq[pos] == 'C'
        elif(pamseq[pos] == 'K'):
            align = seq[pos] == 'G' or seq[pos] == 'T'
        elif(pamseq[pos] == 'R'):
            align = seq[pos] == 'A' or seq[pos] == 'G'
        elif(pamseq[pos] == 'Y'):
            align = seq[pos] == 'T' or seq[pos] == 'C'
        elif(pamseq[pos] == 'V'):
            align = seq[pos] == 'G' or seq[pos] == 'C' or seq[pos] == 'A'
        elif(pamseq[pos] == 'B'):
            align = seq[pos] == 'T' or seq[pos] == 'C' or seq[pos] == 'G'
        elif(pamseq[pos] == 'H'):
            align = seq[pos] == 'T' or seq[pos] == 'C' or seq[pos] == 'A'
        elif(pamseq[pos] == 'D'):
            align = seq[pos] == 'T' or seq[pos] == 'G' or seq[pos] == 'A'
        pos += 1
    return align


def revComp(seq):
    '''Return the reverse complememt of a sequence'''
    res = ''
    for pos in xrange(len(seq)):
        res = dnacomp[seq[pos]] + res
    return res


def findPAM(varseq, fl1seq, fl2seq, mhfl, maxTail, pamseq, pamcut):
    '''Look for PAM cuts between the MH regions.'''
    pamseq_rev = revComp(pamseq)
    seq = fl1seq + varseq + fl2seq
    # Minimum (exact) MH length that must remains after a cut.
    min_mh_cut = min(mhfl['m1L'], 3)
    search_range = [len(fl1seq) - 1, len(fl1seq) + len(varseq) - min_mh_cut]
    reduced_search_range = [len(fl1seq) - 1,
                            len(fl1seq) + len(varseq) - mhfl['mhL']]
    if(mhfl['flank'] == 2):
        search_range = [len(fl1seq) + min_mh_cut - 1,
                        len(fl1seq) + len(varseq)]
        reduced_search_range = [len(fl1seq) + mhfl['mhL'] - 1,
                                len(fl1seq) + len(varseq)]
    # Test each position: if it matched the motif and
    # in the search range, add to list
    pams = []
    for pos in xrange(len(seq)-len(pamseq)+1):
        proto_seq = strand = cut_pos = False
        if(aligned(seq[pos:pos+len(pamseq)], pamseq)):
            cut_pos = pos + pamcut - 1
            strand = '+'
        if(aligned(seq[pos:pos+len(pamseq)], pamseq_rev)):
            cut_pos = pos - pamcut + len(pamseq_rev) - 1
            strand = '-'
        if(strand and cut_pos >= search_range[0] and cut_pos < search_range[1]
           and (cut_pos < search_range[0] + maxTail or
                cut_pos >= search_range[1] - maxTail)):
            if(strand == '+'):
                proto_seq = seq[(cut_pos-19-pamcut):(cut_pos-pamcut+1)]
            else:
                proto_seq = seq[(cut_pos+pamcut+1):(cut_pos+20+pamcut+1)]
            pam_info = {'cutPosition': cut_pos, 'strand': strand,
                        'proto': proto_seq}
            # Distance to MH on each side using first stretch of
            # perfect match or the extended MH
            pam_info['m1Dist1'] = cut_pos - search_range[0]
            pam_info['m1Dist2'] = search_range[1] - cut_pos - 1
            pam_info['mhDist1'] = cut_pos - reduced_search_range[0]
            pam_info['mhDist2'] = reduced_search_range[1] - cut_pos - 1
            pam_info['pamseq'] = pamseq
            pams.append(pam_info)
    return(pams)


def enumN(seq):
    '''Enumerate sequences by replacing Ns by ATCG.'''
    if 'N' not in seq:
        return [seq]
    # Finds the position of the first N
    n_pos = seq.find('N')
    seq = list(seq)
    res = []
    for nuc in ['A', 'T', 'C', 'G']:
        seq2 = seq
        # Replace N
        seq2[n_pos] = nuc
        seq2 = ''.join(seq2)
        # Recursive call in case other Ns are present
        res.extend(enumN(seq2))
    return res


def alignPamsBlast(pams, reffile, include_pam=True, prefix='blast'):
    '''Align protospacers and return an updated version of
    the input "pams" list.'''
    fasta_file = prefix + '_tempMHcut.fasta'
    ff = open(fasta_file, 'w')
    pams_hash = {}
    for pam in pams:
        pam['mm0'] = 0
        pam['mm1'] = 0
        pam['mm2'] = 0
        if('N' in pam['proto']):
            pam['mm0'] = 'NA'
            pam['mm1'] = 'NA'
            pam['mm2'] = 'NA'
            continue
        if(include_pam):
            if(pam['strand'] == '+'):
                protoguide = pam['proto'] + pam['pamseq']
            else:
                protoguide = revComp(pam['pamseq']) + pam['proto']
        else:
            protoguide = pam['proto']
        protoguides = enumN(protoguide)
        for ii in xrange(len(protoguides)):
            pamid = '{}_{}_{}'.format(pam['cutPosition'], pam['strand'], ii)
            pams_hash[pamid] = pam
            ff.write('>' + pamid + '\n' + protoguides[ii] + '\n')
    ff.close()
    dump = open('/dev/null')
    blast_cmd = ['blastn', '-db', reffile,  '-query', fasta_file, '-outfmt',
                 '6', '-word_size', '17', '-max_target_seqs', '20']
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


def alignPamsJellyfish(pams, jffile, include_pam=True, prefix='jf'):
    '''Align protospacers and return an updated version of
    the input "pams" list.'''
    pams_hash = {}
    fasta_file = prefix + '_tempMHcut.fasta'
    ff = open(fasta_file, 'w')
    cpt = 0
    jellyfish_cmd = ['jellyfish', 'query', '-L', jffile, '-s', fasta_file]
    for pam in pams:
        pam['mm0'] = 0
        pam['mm1'] = 'NA'
        pam['mm2'] = 'NA'
        if('N' in pam['proto']):
            pam['mm0'] = 'NA'
            continue
        if(include_pam):
            if(pam['strand'] == '+'):
                protoguide = pam['proto'] + pam['pamseq']
            else:
                protoguide = revComp(pam['pamseq']) + pam['proto']
        else:
            protoguide = pam['proto']
        protoguides = enumN(protoguide)
        for pg in protoguides:
            ff.write('>' + str(cpt) + '\n' + pg + '\n')
            cpt += 1
            if(pg in pams_hash):
                pams_hash[pg].append(pam)
            else:
                pams_hash[pg] = [pam]
            pg_rc = revComp(pg)
            if(pg_rc in pams_hash):
                pams_hash[pg_rc].append(pam)
            else:
                pams_hash[pg_rc] = [pam]
    ff.close()
    if(cpt > 0):
        dump = open('/dev/null')
        jellyfish_out = subprocess.check_output(jellyfish_cmd, stderr=dump)
        dump.close()
        jellyfish_out = jellyfish_out.rstrip('\n').split('\n')
        for line in jellyfish_out:
            line = line.split(' ')
            for pam in pams_hash[line[0]]:
                pam['mm0'] += int(line[1])
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


class RegionExactMH:
    '''A region where exact micro-homology are searched.

    This class is used to list other nested (exact) MHs that might be used
    during the NHEJ. It takes an input sequence and build an alignment array
    with all MH in the sequence. Then for a particular "cut position" it
    quickly lists MH that are appropriate.
    '''
    def __init__(self, seq):
        self.seq = seq
        self.seql = len(seq)
        # Build the alignment array (lower-triangle only using a hash).
        self.ktab = {}
        for ii in range(self.seql):
            self.ktab['-1_' + str(ii-1)] = -1
        for ii in range(self.seql - 1):
            for jj in range(ii + 1, self.seql):
                ii_jj = '{}_{}'.format(ii, jj)
                if(seq[ii] == seq[jj]):
                    iip_jjp = '{}_{}'.format(ii-1, jj-1)
                    if(self.ktab[iip_jjp] != -1):
                        self.ktab[ii_jj] = self.ktab[iip_jjp]
                    else:
                        self.ktab[ii_jj] = ii
                else:
                    self.ktab[ii_jj] = -1

    def listmh(self, cutpos, minmh_size=3):
        res = {}
        # Search for matches in the sub-array defined by the cut position
        for ii in range(cutpos+1)[::-1]:
            for jj in range(cutpos+1, self.seql):
                ii_jj = '{}_{}'.format(ii, jj)
                if self.ktab[ii_jj] != -1:
                    start2 = self.ktab[ii_jj]
                    size = ii - start2 + 1
                    start1 = jj - size + 1
                    # Trim if it crosses the cut
                    if start1 <= cutpos:
                        start_shift = cutpos - start1 + 2
                        start1 += start_shift
                        start2 += start_shift
                        size = ii - start2 + 1
                    start12 = '{}_{}'.format(start1, start2)
                    if(start12 not in res and size >= minmh_size):
                        mh = {'size': size,
                              'seq': self.seq[start1:(start1+size)]}
                        mh['dist'] = start1 - start2 - size
                        mh['vsize'] = start1 - start2
                        mh['startD'] = start1
                        mh['startU'] = start2
                        res[start12] = mh
        # Compute the MMEJ score from Bae et al
        for mh in res:
            length_weight = 20.0
            length_factor = round(1/exp((res[mh]['size'])/(length_weight)), 3)
            num_GC = 0
            for ii in range(res[mh]['size']):
                if(res[mh]['seq'][ii] == 'G' or res[mh]['seq'][ii] == 'C'):
                    num_GC += 1
            score = 100 * length_factor * ((res[mh]['size'] - num_GC) +
                                           (num_GC * 2))
            res[mh]['score'] = score
            res[mh]['gc'] = float(num_GC) / res[mh]['size']
        return res

    def printKtab(self):
        toprint = '    '
        for ii in range(self.seql):
            if(ii > 9):
                toprint += ' ' + str(ii)
            else:
                toprint += '  ' + str(ii)
        toprint += '\n'
        toprint += '0     ' + self.seq[0] + '\n'
        for jj in range(1, self.seql):
            if(jj > 9):
                toprint += str(jj) + ' ' + self.seq[jj]
            else:
                toprint += str(jj) + '  ' + self.seq[jj]
            for ii in range(jj):
                elt = self.ktab[str(ii) + '_' + str(jj)]
                if(elt == -1):
                    elt = ' .'
                else:
                    if(elt > 9):
                        elt = str(elt)
                    else:
                        elt = ' ' + str(elt)
                toprint += ' ' + elt
            toprint += '  ' + self.seq[jj] + '\n'
        print toprint


def mhcut(args):
    # Open connection to reference genome
    print "Check if reference is indexed and index it if not"
    print "...might take a minute..."
    reffa = Fasta(args.reffile)
    print "...Done."

    if(args.varfile == '' and args.outprefix == ''):
        print "Use -var and -out to run MHcut. "
        print "Run 'MHcut -h' for more info."
        sys.exit(0)

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
    outhead = inhead + '\tvarL\tmhL\tmh1L\thom\tnbMM\tmhDist\tMHseq1\tMHseq2\tGC'
    outhead += '\tpamMot\tpamUniq\tguidesNoNMH\tguidesMinNMH\tmax2cutsDist'
    variant_output_file.write(outhead + '\n')
    gouthead = outhead + '\tprotospacer\tpamSeq\tmm0\tmm1\tmm2\tm1Dist1\tm1Dist2'
    gouthead += '\tmhDist1\tmhDist2\tnbNMH\tlargestNMH\tnmhScore\tnmhSize\tnmhVarL'
    gouthead += '\tnmhGC\tnmhSeq\n'
    guide_output_file.write(gouthead)
    cartoon_output_file.write(outhead + '\n\n')

    # Start progress bar
    sys.stdout.write('Completed: 0%')

    # Read each line of the input file
    line_cpt = 0
    if(input_nb_lines > 100):
        input_nb_lines = input_nb_lines / 100
    for input_line in variant_input_file:
        line_cpt += 1
        if(line_cpt % input_nb_lines == 0):
            percent = line_cpt / input_nb_lines
            sys.stdout.write('\rCompleted: ' + str(percent) + '%')
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
        varseq = str(reffa[input_line[0]][(vstart-1):vend]).upper()
        fl1seq = str(reffa[input_line[0]][(vstart-flsize-1):(vstart-1)]).upper()
        fl2seq = str(reffa[input_line[0]][vend:(vend+flsize)]).upper()
        # Test MH in each flank (reverse for flank 1) and save best MH
        mhfl1 = mhTest(varseq[::-1], fl1seq[::-1], args.maxConsMM)
        # If no MH or too small, or too low MH ratio or
        # too short first microhomology stretch
        if(not args.nofilter
           and (mhfl1['mhL'] < args.minMHL or mhfl1['hom'] < args.minhom or
                mhfl1['m1L'] < args.minm1L)):
            mhfl1['score'] = 0
        # Same for other flank
        mhfl2 = mhTest(varseq, fl2seq, args.maxConsMM)
        if(not args.nofilter
           and (mhfl2['mhL'] < args.minMHL or mhfl2['hom'] < args.minhom or
                mhfl2['m1L'] < args.minm1L)):
            mhfl2['score'] = 0
        # Using the alignment score, the best flank is chosen
        if(mhfl1['score'] > mhfl2['score']):
            mhfl = mhfl1
            mhfl['flank'] = 1  # This is to remember which flank is chosen
            mhfl['cartoon'] = mhfl['cartoon'][::-1]
        else:
            mhfl = mhfl2
            mhfl['flank'] = 2
        # If a score of 0, either no MH or didn't satisfy criteria above,
        # jump to the next input line
        if(mhfl['score'] == 0):
            if(args.nofilter):
                # Write line
                voutline = '\t'.join([input_line_raw, str(vsize), str(mhfl['mhL']),
                                      str(mhfl['m1L']), str(round(mhfl['hom'], 2)),
                                      str(mhfl['nbMM'])])
                voutline += '\tNA' * 9
                variant_output_file.write(voutline + '\n')
            continue
        # Find PAM motives
        pamseqs = args.pamseq.split(',')
        pams = []
        for pamseq in pamseqs:
            pams.extend(findPAM(varseq, fl1seq, fl2seq, mhfl,
                                args.maxTail, pamseq, args.pamcut))
        # Map protospacers to the genome and keep unique ones
        nb_pam_motives = len(pams)
        max2cutsDist = 'NA'
        mincut = maxcut = ''
        if(nb_pam_motives > 0):
            if(args.jffile == ''):
                pams = alignPamsBlast(pams, args.reffile, args.outprefix)
            else:
                pams = alignPamsJellyfish(pams, args.jffile, args.outprefix)
            pams_filter = []
            for pam in pams:
                # This is where to define how unique the protospacer must be
                # With mm0=1, there must be only one position
                # in the genome aligning perfectly
                if pam['mm0'] == 1:
                    pams_filter.append(pam)
                    if(mincut == ''):
                        mincut = maxcut = pam['cutPosition']
                    else:
                        mincut = min(mincut, pam['cutPosition'])
                        maxcut = max(maxcut, pam['cutPosition'])
            pams = pams_filter
            if(mincut != '' and len(pams) > 1):
                max2cutsDist = maxcut - mincut
        # Search for other MH that could be used by the MMEJ
        if len(pams) > 0 and vsize < args.maxTail*2:
            other_mh = RegionExactMH(fl1seq + varseq + fl2seq)
        else:
            other_mh = False
        no_offtargets = 'NA'  # Number of guides with no off target MH
        min_offtargets = 'NA'
        for pam in pams:
            # Number of off target and maximum size (no matter the score)
            pam['nmh_nb'] = 'NA'
            pam['nmh_maxL'] = 'NA'
            # Best nested MH stats
            pam['bnmh_score'] = 'NA'
            pam['bnmh_size'] = 'NA'
            pam['bnmh_vsize'] = 'NA'
            pam['bnmh_gc'] = 'NA'
            pam['bnmh_seq'] = 'NA'
            if other_mh:
                pam['nmh_nb'] = 0
                pam['nmh_maxL'] = 0
                pam['bnmh_score'] = 0
                mhhet = other_mh.listmh(pam['cutPosition'], args.minLnmh)
                for mho in mhhet:
                    data = mhhet[mho]
                    # Only consider other MH that at least as close
                    # from each other as our target MH.
                    if(data['vsize'] <= vsize and
                       (data['startU'] != flsize or
                        data['startD'] != flsize + vsize)
                       and (data['startU'] + data['size'] != flsize or
                            data['startD'] + data['size'] != flsize + vsize)):
                        pam['nmh_nb'] += 1
                        pam['nmh_maxL'] = max(pam['nmh_maxL'], data['size'])
                        if data['score'] > pam['bnmh_score']:
                            pam['bnmh_score'] = data['score']
                            pam['bnmh_size'] = data['size']
                            pam['bnmh_vsize'] = data['vsize']
                            pam['bnmh_gc'] = round(data['gc'], 3)
                            pam['bnmh_seq'] = data['seq']
                if no_offtargets == 'NA':
                    no_offtargets = 0
                if pam['nmh_nb'] == 0:
                    no_offtargets += 1
                if min_offtargets == 'NA':
                    min_offtargets = pam['nmh_nb']
                else:
                    min_offtargets = min(min_offtargets, pam['nmh_nb'])
        # Write in output files
        # Add/remove columns here (without forgetting the header)
        voutline = '\t'.join([input_line_raw, str(vsize), str(mhfl['mhL']),
                              str(mhfl['m1L']), str(round(mhfl['hom'], 2)),
                              str(mhfl['nbMM']), str(mhfl['mhdist']),
                              mhfl['seq1'], mhfl['seq2'], str(mhfl['gc']),
                              str(nb_pam_motives), str(len(pams)),
                              str(no_offtargets), str(min_offtargets),
                              str(max2cutsDist)])
        variant_output_file.write(voutline + '\n')
        for pam in pams:
            goutline = '\t'.join([voutline, pam['proto'], pam['pamseq'],
                                  str(pam['mm0']), str(pam['mm1']),
                                  str(pam['mm2']), str(pam['m1Dist1']),
                                  str(pam['m1Dist2']), str(pam['mhDist1']),
                                  str(pam['mhDist2']), str(pam['nmh_nb']),
                                  str(pam['nmh_maxL']), str(pam['bnmh_score']),
                                  str(pam['bnmh_size']), str(pam['bnmh_vsize']),
                                  str(pam['bnmh_gc']), pam['bnmh_seq']])
            guide_output_file.write(goutline + '\n')
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
                if(pam_cartoon[pam['cutPosition'] + 1] != '_'):
                    pam_cartoon[pam['cutPosition'] + 1] = 'X'
                else:
                    pam_cartoon[pam['cutPosition'] + 1] = '\\'
            else:
                if(pam_cartoon[pam['cutPosition']] != '_'):
                    pam_cartoon[pam['cutPosition']] = 'X'
                else:
                    pam_cartoon[pam['cutPosition']] = '/'
        pam_cartoon = ''.join(pam_cartoon)
        cartoon_output_lines[2] = pam_cartoon[:len(fl1seq)] + ' '
        cartoon_output_lines[2] += pam_cartoon[len(fl1seq):(len(fl1seq) + vsize)]
        cartoon_output_lines[2] += ' ' + pam_cartoon[(len(fl1seq) + vsize):]
        # If the line is too long (large variants), trim the ends and middle
        # How many extra bases to show on the flanks (outside of the MH region)
        flank_buffer = 10
        if(len(cartoon_output_lines[1]) > 2 *
           (flank_buffer + mhfl['mhL'] + args.maxTail)):
            cartoon_part1 = [len(fl1seq) - mhfl['mhL'] - flank_buffer,
                             len(fl1seq) + args.maxTail]
            cartoon_part2 = [len(fl1seq) + vsize - mhfl['mhL'] - args.maxTail,
                             len(fl1seq) + vsize + flank_buffer]
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
                    line += full_line[cartoon_part2[0]:min(cartoon_part2[1],
                                                           len(full_line))]
                    cartoon_output_lines_trimmed.append(line)
            cartoon_output_lines = cartoon_output_lines_trimmed
        # Write cartoon lines
        for line in cartoon_output_lines:
            cartoon_output_file.write(line + '\n')
        # Cartoon: protospacers sequence
        if(len(pams) > 0):
            cartoon_output_file.write('Protospacers, m1Dist1, m1Dist2, mhDist1, '
                                      'mhDist2, bestOffTarget:\n')
            for pam in pams:
                coutline = '\t'.join([pam['proto'], str(pam['m1Dist1']),
                                      str(pam['m1Dist2']), str(pam['mhDist1']),
                                      str(pam['mhDist2']), str(pam['bnmh_seq'])])
                cartoon_output_file.write(coutline + '\n')
        cartoon_output_file.write('\n\n')

    print '\nDone.\n'

    variant_input_file.close()
    variant_output_file.close()
    guide_output_file.close()
    cartoon_output_file.close()
