class Variant():
    '''A variant with flanking sequence.'''
    # A bit of a useless class but nicer than passing multiple objects
    # or a dictionnary around.
    def __init__(self, v_chr, vstart, vend, reffa):
        self.vchr = v_chr
        self.vstart = vstart
        self.vend = vend
        self.reffa = reffa
        self.vsize = vend - vstart + 1
        # Retrieve sequences. Careful with position and shift
        self.flsize = max(self.vsize, 60)
        self.varseq = reffa[v_chr][(vstart-1):vend]
        self.fl1seq = reffa[v_chr][(vstart-self.flsize-1):(vstart-1)]
        self.fl2seq = reffa[v_chr][vend:(vend+self.flsize)]
        # Force to uppercase strings
        self.varseq = str(self.varseq).upper()
        self.fl1seq = str(self.fl1seq).upper()
        self.fl2seq = str(self.fl2seq).upper()

    def right_shift_possible(self):
        return(self.varseq[0] == self.fl2seq[0])

    def left_shift_possible(self):
        return(self.varseq[-1] == self.fl1seq[-1])


def list_shifted_variants(var):
    '''Lists all possible variants that don't change the deletion'''
    vars = [var]
    # Shift right
    cur_var = var
    while cur_var.right_shift_possible():
        cur_var = Variant(cur_var.vchr, cur_var.vstart+1, cur_var.vend+1,
                          cur_var.reffa)
        vars.append(cur_var)
    # Shift left
    cur_var = var
    while cur_var.left_shift_possible():
        cur_var = Variant(cur_var.vchr, cur_var.vstart-1, cur_var.vend-1,
                          cur_var.reffa)
        vars.append(cur_var)
    return vars
