'''
Functions to manipulate sequences.
E.g. reverse complement, alignment, GC content.
'''

from math import exp

# DNA complement
dnacomp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'W': 'W', 'S': 'S',
           'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R', 'V': 'B', 'B': 'V',
           'H': 'D', 'D': 'H', 'N': 'N'}


def GC(seq):
    '''Computes the GC content of a sequence.'''
    seq = list(seq)
    gc_count = seq.count('C') + seq.count('G')
    gc_prop = float(gc_count) / len(seq)
    return round(gc_prop, 3)


def revComp(seq):
    '''Return the reverse complement of a sequence'''
    res = ''
    for pos in xrange(len(seq)):
        res = dnacomp[seq[pos]] + res
    return res


def aligned(seq, pamseq):
    '''Are seq and pamseq matching?
    Takes NWSMKRYVBHD into account.'''
    align = True
    pos = 0
    while(align and pos < len(pamseq)):
        if pamseq[pos] != 'N':
            align = pamseq[pos] == seq[pos]
        if pamseq[pos] == 'W':
            align = seq[pos] == 'A' or seq[pos] == 'T'
        elif pamseq[pos] == 'S':
            align = seq[pos] == 'G' or seq[pos] == 'C'
        elif pamseq[pos] == 'M':
            align = seq[pos] == 'A' or seq[pos] == 'C'
        elif pamseq[pos] == 'K':
            align = seq[pos] == 'G' or seq[pos] == 'T'
        elif pamseq[pos] == 'R':
            align = seq[pos] == 'A' or seq[pos] == 'G'
        elif pamseq[pos] == 'Y':
            align = seq[pos] == 'T' or seq[pos] == 'C'
        elif pamseq[pos] == 'V':
            align = seq[pos] == 'G' or seq[pos] == 'C' or seq[pos] == 'A'
        elif pamseq[pos] == 'B':
            align = seq[pos] == 'T' or seq[pos] == 'C' or seq[pos] == 'G'
        elif pamseq[pos] == 'H':
            align = seq[pos] == 'T' or seq[pos] == 'C' or seq[pos] == 'A'
        elif pamseq[pos] == 'D':
            align = seq[pos] == 'T' or seq[pos] == 'G' or seq[pos] == 'A'
        pos += 1
    return align


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
                if seq[ii] == seq[jj]:
                    iip_jjp = '{}_{}'.format(ii-1, jj-1)
                    if self.ktab[iip_jjp] != -1:
                        self.ktab[ii_jj] = self.ktab[iip_jjp]
                    else:
                        self.ktab[ii_jj] = ii
                else:
                    self.ktab[ii_jj] = -1

    def listmh(self, cutpos, minmh_size=3):
        '''List MH on each side of a cut position.'''
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
                    if start12 not in res and size >= minmh_size:
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
                if res[mh]['seq'][ii] == 'G' or res[mh]['seq'][ii] == 'C':
                    num_GC += 1
            score = 100 * length_factor * ((res[mh]['size'] - num_GC) +
                                           (num_GC * 2))
            res[mh]['score'] = score
            res[mh]['gc'] = float(num_GC) / res[mh]['size']
        return res

    def printKtab(self):
        '''Prints a table (used during implementation and debugging).'''
        toprint = '    '
        for ii in range(self.seql):
            if ii > 9:
                toprint += ' ' + str(ii)
            else:
                toprint += '  ' + str(ii)
        toprint += '\n'
        toprint += '0     ' + self.seq[0] + '\n'
        for jj in range(1, self.seql):
            if jj > 9:
                toprint += str(jj) + ' ' + self.seq[jj]
            else:
                toprint += str(jj) + '  ' + self.seq[jj]
            for ii in range(jj):
                elt = self.ktab[str(ii) + '_' + str(jj)]
                if elt == -1:
                    elt = ' .'
                else:
                    if elt > 9:
                        elt = str(elt)
                    else:
                        elt = ' ' + str(elt)
                toprint += ' ' + elt
            toprint += '  ' + self.seq[jj] + '\n'
        print toprint
