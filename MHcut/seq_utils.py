dnacomp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'W': 'W', 'S': 'S',
           'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R', 'V': 'B', 'B': 'V',
           'H': 'D', 'D': 'H', 'N': 'N'}


def revComp(seq):
    '''Return the reverse complememt of a sequence'''
    res = ''
    for pos in xrange(len(seq)):
        res = dnacomp[seq[pos]] + res
    return res


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
