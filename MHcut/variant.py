class Variant():
    '''A variant with flanking sequence.'''
    def __init__(self, v_chr, vstart, vend, reffa):
        self.vsize = vend - vstart + 1
        # Retrieve sequences. Careful with position and shift
        self.flsize = max(self.vsize, 20)
        self.varseq = reffa[v_chr][(vstart-1):vend]
        self.fl1seq = reffa[v_chr][(vstart-self.flsize-1):(vstart-1)]
        self.fl2seq = reffa[v_chr][vend:(vend+self.flsize)]
        # Force to uppercase strings
        self.varseq = str(self.varseq).upper()
        self.fl1seq = str(self.fl1seq).upper()
        self.fl2seq = str(self.fl2seq).upper()
