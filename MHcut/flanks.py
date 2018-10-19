'''
Class to represent a MH configuration for a variant.

Two flanking configurations exists: outer-inner (1) and inner-outer (2).
If flank1-variant-flank2, configuration 1 is MH between flank1 and variant,
configuration 2 is MH between flank2 and variant. In other words:
Outer-inner: MH between 3' variant sequence and 5' flanking sequence.
Inner-outer: MH between 5' variant sequence and 3' flanking sequence.
'''

import seq_utils


class VarFlank():
    '''
    A pair of variant and flanking sequences.
    Finds MH and contains MH metrics.
    '''
    def __init__(self, var, flank=1):
        self.var = var
        # The flank configuration, see details above
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
        if self.flank == 2:
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
        if al_full[0]:
            # Trim the end of the alignment if X consecutive mismatches
            al_trimmed = al_full
            consMM = 0
            for pos in range(len(al_full)-1):
                if not al_full[pos]:
                    consMM += 1
                else:
                    consMM = 0
                if consMM > max_cons_mhh:
                    al_trimmed = al_full[:(pos-consMM+1)]
                    break
            # Cut potential last mismatch
            while(not al_trimmed[-1]):
                al_trimmed = al_trimmed[:-1]
            # Count consecutive matches in the beginning
            for pos in range(len(al_trimmed)):
                if al_trimmed[pos]:
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
            if self.flank == 1:
                self.mh_cartoon = ''.join(cartoon[::-1])
            else:
                self.mh_cartoon = ''.join(cartoon)
            # The two homologous sequences
            self.inner_mh_seq = mht_varseq[0:self.mhL]
            self.outer_mh_seq = mht_flseq[0:self.mhL]
            # Maximum GC content of the homologous sequences
            gc1 = seq_utils.GC(self.inner_mh_seq)
            gc2 = seq_utils.GC(self.outer_mh_seq)
            self.gc = max(gc1, gc2)

    def toString(self, flank_info=False):
        '''The relevant string to write in the "variants" output.'''
        # IF YOU CHANGE SOMETHING HERE, CHANGE THE HEADERS TOO (below)
        tostr = [self.mhL, self.m1L, round(self.hom, 2), self.nbMM,
                 self.mhdist, self.inner_mh_seq, self.outer_mh_seq, self.gc]
        if flank_info:
            tostr = [self.flank, self.score] + tostr
        tostr = '\t'.join([str(ii) for ii in tostr])
        return tostr


def headers(flank_info=False):
    '''
    The headers corresponding to the "variants" output from the
    "toString" output.
    '''
    # IF YOU CHANGE SOMETHING HERE, CHANGE THE toString TOO (above)
    headers = ['mhL', 'mh1L', 'hom', 'nbMM', 'mhDist', 'MHseq1', 'MHseq2',
               'GC']
    if flank_info:
        headers = ['flank', 'mhScore'] + headers
    return headers
