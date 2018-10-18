def drawCartoon(var, var_fl, pams, max_tail=50):
    out_lines = ['', '', '']
    # Cartoon: alignment line
    white_spaces_before = len(var.fl1seq) - var_fl.mhL
    if(var_fl.flank == 2):
        white_spaces_before += var_fl.mhL + 1
    white_spaces_before = ' ' * white_spaces_before
    white_spaces_between = ' ' * (var.vsize - var_fl.mhL + 1)
    out_lines[0] = white_spaces_before + var_fl.mh_cartoon
    out_lines[0] += white_spaces_between + var_fl.mh_cartoon
    # Cartoon: sequence line
    out_lines[1] = var.fl1seq + '-' + var.varseq + '-' + var.fl2seq
    # Cartoon: PAM position line
    size_cartoon = len(var.fl1seq) + len(var.fl2seq) + var.vsize
    # Init with '_' everywhere
    pam_cartoon = ['_' for i in range(size_cartoon)]
    for pam in pams.pams:
        if(pam.uniq):
            if(pam.strand == '+'):
                if(pam_cartoon[pam.cutPosition + 1] != '_'):
                    pam_cartoon[pam.cutPosition + 1] = 'X'
                else:
                    pam_cartoon[pam.cutPosition + 1] = '\\'
            else:
                if(pam_cartoon[pam.cutPosition] != '_'):
                    pam_cartoon[pam.cutPosition] = 'X'
                else:
                    pam_cartoon[pam.cutPosition] = '/'
    pam_cartoon = ''.join(pam_cartoon)
    out_lines[2] = pam_cartoon[:len(var.fl1seq)] + ' '
    out_lines[2] += pam_cartoon[len(var.fl1seq):(len(var.fl1seq) + var.vsize)]
    out_lines[2] += ' ' + pam_cartoon[(len(var.fl1seq) + var.vsize):]
    # If the line is too long (large variants), trim the ends and middle
    # How many extra bases to show on the flanks (outside of the MH region)
    flank_buffer = 10
    if(len(out_lines[1]) > 2 *
       (flank_buffer + var_fl.mhL + max_tail)):
        cartoon_part1 = [len(var.fl1seq) - var_fl.mhL - flank_buffer,
                         len(var.fl1seq) + max_tail]
        cartoon_part2 = [len(var.fl1seq) + var.vsize - var_fl.mhL - max_tail,
                         len(var.fl1seq) + var.vsize + flank_buffer]
        # If second flank was used, shift the positions
        if(var_fl.flank == 2):
            cartoon_part1 = [p + var_fl.mhL + 1 for p in cartoon_part1]
            cartoon_part2 = [p + var_fl.mhL + 1 for p in cartoon_part2]
        # Sanity checks that the position are in the correct range
        cartoon_part1[0] = max(cartoon_part1[0], 0)
        # Update lines (if the two part overlaps merge them into one)
        out_lines_trimmed = []
        if(cartoon_part1[1] + 1 > cartoon_part2[0]):
            cartoon_part1[1] = cartoon_part2[1]
            for line in out_lines:
                line = line[cartoon_part1[0]:min(cartoon_part1[1], len(line))]
                out_lines_trimmed.append(line)
        else:
            for line in out_lines:
                full_line = line
                line = full_line[cartoon_part1[0]:cartoon_part1[1]]
                line += '...'
                line += full_line[cartoon_part2[0]:min(cartoon_part2[1],
                                                       len(full_line))]
                out_lines_trimmed.append(line)
        out_lines = out_lines_trimmed
    # Cartoon: protospacers sequence
    if(len(pams.pams) > 0):
        proto_head = 'Protospacers, m1Dist1, m1Dist2, mhDist1, '
        proto_head += 'mhDist2, bestOffTarget:'
        out_lines.append(proto_head)
        for pam in pams.pams:
            coutline = '\t'.join([pam.proto, str(pam.m1Dist1),
                                  str(pam.m1Dist2), str(pam.mhDist1),
                                  str(pam.mhDist2), str(pam.bnmh_seq)])
            out_lines.append(coutline)
    return('\n'.join(out_lines))
