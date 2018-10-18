def GC(seq):
    seq = list(seq)
    gc_count = seq.count('C') + seq.count('G')
    gc_prop = float(gc_count) / len(seq)
    return round(gc_prop, 3)


class VarFlank():
    '''A pair of variant and flanking sequences.'''
    def __init__(self, var, flank=1):
        self.var = var
        self.flank = flank
        # Init MH stats
        self.score = 0
        self.m1L = 0
        self.mhL = 0
        self.hom = 0
        self.nbMM = 0
        self.mh_cartoon = '',
        self.inner_mh_seq = 'NA'
        self.outer_mh_seq = 'NA'
        self.mhdist = 'NA'
        self.gc = 'NA'

    def findMH(self, max_cons_mhh=1):
        '''Test for presence of microhomology between the two sequences.'''
        if(self.flank == 2):
            mht_varseq = self.var.varseq
            mht_flseq = self.var.fl2seq
        else:
            mht_varseq = self.var.varseq[::-1]
            mht_flseq = self.var.fl1seq[::-1]
        # ALign the two sequences
        al_full = []
        for pos in range(min(len(mht_varseq), len(mht_flseq))):
            al_full.append(mht_varseq[pos] == mht_flseq[pos]
                           and mht_varseq[pos] != 'N')
        # First base must match
        if(al_full[0]):
            # Trim the end of the alignment if X consecutive mismatches
            al_trimmed = al_full
            consMM = 0
            for pos in range(len(al_full)-1):
                if(not al_full[pos]):
                    consMM += 1
                else:
                    consMM = 0
                if(consMM > max_cons_mhh):
                    al_trimmed = al_full[:(pos-consMM+1)]
                    break
            # Cut potential last mismatch
            while(not al_trimmed[-1]):
                al_trimmed = al_trimmed[:-1]
            # Count consecutive matches in the beginning
            for pos in range(len(al_trimmed)):
                if(al_trimmed[pos]):
                    self.m1L += 1
                else:
                    break
            # Other alignment stats: length, homology, nb of mismatches
            nb_match = sum(al_trimmed)
            self.mhL = len(al_trimmed)
            self.hom = float(nb_match) / self.mhL
            self.nbMM = self.mhL - nb_match
            # dist btw homologous sequences
            self.mhdist = len(al_full) - self.mhL
            # Compute a score, later used to choose which flank has the best MH
            self.score = self.m1L + nb_match
            # Cartoon of the MH (e.g. ||x|)
            cartoon = ['|' if al else 'x' for al in al_trimmed]
            if(self.flank == 1):
                self.mh_cartoon = ''.join(cartoon[::-1])
            else:
                self.mh_cartoon = ''.join(cartoon)
            # The two homologous sequences
            self.inner_mh_seq = mht_varseq[0:self.mhL]
            self.outer_mh_seq = mht_flseq[0:self.mhL]
            # Maximum GC content of the homologous sequences
            gc1 = GC(self.inner_mh_seq)
            gc2 = GC(self.outer_mh_seq)
            self.gc = max(gc1, gc2)

    def toString(self, flank_info=False):
        tostr = [self.mhL, self.m1L, round(self.hom, 2), self.nbMM,
                 self.mhdist, self.inner_mh_seq, self.outer_mh_seq, self.gc]
        if flank_info:
            tostr = [self.flank, self.score] + tostr
        tostr = '\t'.join([str(ii) for ii in tostr])
        return tostr


def headers(flank_info=False):
    headers = ['mhL', 'mh1L', 'hom', 'nbMM', 'mhDist', 'MHseq1', 'MHseq2',
               'GC']
    if flank_info:
        headers = ['flank', 'mhScore'] + headers
    return headers
