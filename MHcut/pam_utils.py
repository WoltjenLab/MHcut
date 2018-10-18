import seq_utils
import subprocess
import os
from math import exp


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


class PAM():
    '''Contains all the information about one PAM.'''
    def __init__(self, pamseq, cutPosition, strand, proto, m1Dist1,
                 m1Dist2, mhDist1, mhDist2):
        self.pamseq = pamseq
        self.cutPosition = cutPosition
        self.strand = strand
        self.proto = proto
        self.m1Dist1 = m1Dist1
        self.m1Dist2 = m1Dist2
        self.mhDist1 = mhDist1
        self.mhDist2 = mhDist2
        self.uniq = False
        # Protospacer alignment
        self.mm0 = 0
        self.mm1 = 0
        self.mm2 = 0
        # Nested MH
        # Number of nested MH and maximum size (no matter the score)
        self.nmh_nb = 'NA'
        self.nmh_maxL = 'NA'
        # Best nested MH stats
        self.bnmh_score = 'NA'
        self.bnmh_size = 'NA'
        self.bnmh_vsize = 'NA'
        self.bnmh_gc = 'NA'
        self.bnmh_seq = 'NA'


class PAMs():
    '''List of PAMs. Useful for bulk processing (protospacer alignment).'''
    def __init__(self, var_fl, pam_seqs, cut_offset, max_tail=50):
        self.pams = []
        for pamseq in pam_seqs:
            self.findPAMs(var_fl, pamseq, cut_offset, max_tail)
        # Global stats
        self.nb_pam_motives = 0
        # Number of guides with no off target MH
        self.no_offtargets = 'NA'
        self.min_offtargets = 'NA'

    def findPAMs(self, var_fl, pamseq, cut_offset, max_tail):
        pamseq_rev = seq_utils.revComp(pamseq)
        seq = var_fl.var.fl1seq + var_fl.var.varseq + var_fl.var.fl2seq
        # Minimum (exact) MH length that must remains after a cut.
        min_mh_cut = min(var_fl.m1L, 3)
        fl_L = len(var_fl.var.fl1seq)
        var_L = len(var_fl.var.varseq)
        search_range = [fl_L - 1, fl_L + var_L - min_mh_cut]
        reduced_search_range = [fl_L - 1,
                                fl_L + var_L - var_fl.mhL]
        if(var_fl.flank == 2):
            search_range = [fl_L + min_mh_cut - 1,
                            fl_L + var_L]
            reduced_search_range = [fl_L + var_fl.mhL - 1,
                                    fl_L + var_L]
        # Test each position: if it matched the motif and
        # in the search range, add to list
        for pos in xrange(len(seq)-len(pamseq)+1):
            proto_seq = strand = cut_pos = False
            if(seq_utils.aligned(seq[pos:pos+len(pamseq)], pamseq)):
                cut_pos = pos + cut_offset - 1
                strand = '+'
            if(seq_utils.aligned(seq[pos:pos+len(pamseq)], pamseq_rev)):
                cut_pos = pos - cut_offset + len(pamseq_rev) - 1
                strand = '-'
            if(strand and cut_pos >= search_range[0]
               and cut_pos < search_range[1]
               and (cut_pos < search_range[0] + max_tail or
                    cut_pos >= search_range[1] - max_tail)):
                if(strand == '+'):
                    proto_start = cut_pos - 19 - cut_offset
                    proto_end = cut_pos - cut_offset + 1
                else:
                    proto_start = cut_pos + cut_offset + 1
                    proto_end = cut_pos + 20 + cut_offset + 1
                proto_seq = seq[proto_start:proto_end]
                # Info about the PAM found, including
                # Distance to MH on each side using first stretch of
                # perfect match or the extended MH
                pam = PAM(pamseq=pamseq, cutPosition=cut_pos,
                          strand=strand, proto=proto_seq,
                          m1Dist1=cut_pos - search_range[0],
                          m1Dist2=search_range[1] - cut_pos - 1,
                          mhDist1=cut_pos - reduced_search_range[0],
                          mhDist2=reduced_search_range[1] - cut_pos - 1)
                self.pams.append(pam)

    def nbPAMs(self):
        return(len(self.pams))

    def getMax2cutsDist(self, uniq_pam_only=True):
        mincut = maxcut = ''
        nb_pams = 0
        for pam in self.pams:
            if(uniq_pam_only and not pam.uniq):
                continue
            else:
                nb_pams += 1
            if(mincut == ''):
                mincut = maxcut = pam.cutPosition
            else:
                mincut = min(mincut, pam.cutPosition)
                maxcut = max(maxcut, pam.cutPosition)
        if(mincut != '' and nb_pams > 1):
            return(maxcut - mincut)
        else:
            return('NA')

    def annotateUnique(self, max_mm0=1, max_mm1=float('inf'),
                       max_mm2=float('inf')):
        # This is where to define how unique the protospacer must be
        # By default max_mm0=1, i.e. there must be only one position
        # in the genome aligning perfectly
        for pam in self.pams:
            # if(pam.mm0 > 0 and pam.mm0 <= max_mm0 and
            #    pam.mm1 <= max_mm1 and pam.mm2 <= max_mm2):
            if(pam.mm0 == 1):
                pam.uniq = True

    def alignPamsBlast(self, reffile, include_pam=True, prefix='blast'):
        '''Align protospacers and return an updated version of
        the input "pams" list.'''
        fasta_file = prefix + '_tempMHcut.fasta'
        ff = open(fasta_file, 'w')
        pams_hash = {}
        for pam in self.pams:
            if('N' in pam.proto):
                pam.mm0 = 'NA'
                pam.mm1 = 'NA'
                pam.mm2 = 'NA'
                continue
            if(include_pam):
                if(pam.strand == '+'):
                    protoguide = pam.proto + pam.pamseq
                else:
                    protoguide = seq_utils.revComp(pam.pamseq) + pam.proto
            else:
                protoguide = pam.proto
            protoguides = seq_utils.enumN(protoguide)
            for ii in xrange(len(protoguides)):
                pamid = '{}_{}_{}'.format(pam.cutPosition, pam.strand, ii)
                pams_hash[pamid] = pam
                ff.write('>' + pamid + '\n' + protoguides[ii] + '\n')
        ff.close()
        dump = open('/dev/null')
        blast_cmd = ['blastn', '-db', reffile,  '-query', fasta_file,
                     '-outfmt', '6', '-word_size', '10',
                     '-max_target_seqs', '20']
        blast_out = subprocess.check_output(blast_cmd, stderr=dump)
        dump.close()
        blast_out = blast_out.split('\n')
        for line in blast_out:
            line = line.split('\t')
            if(len(line) > 1):
                pam_cand = pams_hash[line[0]]
                if(int(line[3]) == len(pam_cand.proto)):
                    if(int(line[4]) == 0):
                        pam_cand.mm0 += 1
                    if(int(line[4]) == 1):
                        pam_cand.mm1 += 1
                    if(int(line[4]) == 2):
                        pam_cand.mm2 += 1
        os.remove(fasta_file)

    def alignPamsJellyfish(self, jffile, include_pam=True, prefix='jf'):
        '''Align protospacers and return an updated version of
        the input "pams" list.'''
        pams_hash = {}
        fasta_file = prefix + '_tempMHcut.fasta'
        ff = open(fasta_file, 'w')
        cpt = 0
        jellyfish_cmd = ['jellyfish', 'query', '-L', jffile, '-s', fasta_file]
        for pam in self.pams:
            pam.mm0 = 0
            pam.mm1 = 'NA'
            pam.mm2 = 'NA'
            if('N' in pam.proto):
                pam.mm0 = 'NA'
                continue
            if(include_pam):
                if(pam.strand == '+'):
                    protoguide = pam.proto + pam.pamseq
                else:
                    protoguide = seq_utils.revComp(pam.pamseq) + pam.proto
            else:
                protoguide = pam.proto
            protoguides = seq_utils.enumN(protoguide)
            for pg in protoguides:
                ff.write('>' + str(cpt) + '\n' + pg + '\n')
                cpt += 1
                if(pg in pams_hash):
                    pams_hash[pg].append(pam)
                else:
                    pams_hash[pg] = [pam]
                pg_rc = seq_utils.revComp(pg)
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
                    pam.mm0 += int(line[1])
        os.remove(fasta_file)

    def findNestedMH(self, var, max_tail=50, min_l_nmh=3,
                     uniq_pam_only=False):
        # Search for other MH that could be used by the MMEJ
        pams = []
        if(uniq_pam_only):
            for pam in self.pams:
                if(pam.uniq):
                    pams.append(pam)
        else:
            pams = self.pams
        if len(pams) > 0 and var.vsize < max_tail*2:
            other_mh = RegionExactMH(var.fl1seq + var.varseq + var.fl2seq)
        else:
            other_mh = False
        for pam in pams:
            if other_mh:
                pam.nmh_nb = 0
                pam.nmh_maxL = 0
                pam.bnmh_score = 0
                mhhet = other_mh.listmh(pam.cutPosition, min_l_nmh)
                for mho in mhhet:
                    data = mhhet[mho]
                    # Only consider other MH that at least as close
                    # from each other as our target MH.
                    v_fl_size = var.flsize + var.vsize
                    if(data['vsize'] <= var.vsize and
                       (data['startU'] != var.flsize or
                        data['startD'] != v_fl_size)
                       and (data['startU'] + data['size'] != var.flsize or
                            data['startD'] + data['size'] != v_fl_size)):
                        pam.nmh_nb += 1
                        pam.nmh_maxL = max(pam.nmh_maxL, data['size'])
                        if data['score'] > pam.bnmh_score:
                            pam.bnmh_score = data['score']
                            pam.bnmh_size = data['size']
                            pam.bnmh_vsize = data['vsize']
                            pam.bnmh_gc = round(data['gc'], 3)
                            pam.bnmh_seq = data['seq']
                if self.no_offtargets == 'NA':
                    self.no_offtargets = 0
                if pam.nmh_nb == 0:
                    self.no_offtargets += 1
                if self.min_offtargets == 'NA':
                    self.min_offtargets = pam.nmh_nb
                else:
                    self.min_offtargets = min(self.min_offtargets, pam.nmh_nb)

    def toStringVariants(self):
        # Add/remove columns here (without forgetting the header)
        pams_uniq = 0
        for pam in self.pams:
            if(pam.uniq):
                pams_uniq += 1
        tostr = [self.nbPAMs(), pams_uniq, self.no_offtargets,
                 self.min_offtargets, self.getMax2cutsDist()]
        tostr = '\t'.join([str(ii) for ii in tostr])
        return tostr

    def toStringGuides(self, voutline, uniq_pam_only=True):
        tostr = ''
        for pam in self.pams:
            if(uniq_pam_only and not pam.uniq):
                continue
            pam_str = [pam.proto, pam.pamseq, pam.mm0, pam.mm1, pam.mm2,
                       pam.m1Dist1, pam.m1Dist2, pam.mhDist1,  pam.mhDist2,
                       pam.nmh_nb, pam.nmh_maxL, pam.bnmh_score, pam.bnmh_size,
                       pam.bnmh_vsize, pam.bnmh_gc, pam.bnmh_seq]
            pam_str = '\t'.join([str(ii) for ii in pam_str])
            tostr += voutline + '\t' + pam_str + '\n'
        return tostr


def headersVariants(NAs=False):
    headers = ['pamMot', 'pamUniq', 'guidesNoNMH', 'guidesMinNMH',
               'max2cutsDist']
    if(NAs):
        return(['NA' for ii in range(len(headers))])
    else:
        return(headers)


def headersGuides():
    return ['protospacer', 'pamSeq', 'mm0', 'mm1', 'mm2', 'm1Dist1', 'm1Dist2',
            'mhDist1', 'mhDist2', 'nbNMH', 'largestNMH', 'nmhScore', 'nmhSize',
            'nmhVarL', 'nmhGC', 'nmhSeq']
