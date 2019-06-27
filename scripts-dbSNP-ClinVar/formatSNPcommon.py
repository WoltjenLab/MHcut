import fileinput


for line in fileinput.input():
    line = line.rstrip().split('\t')
    rsid = line[4]
    ref = line[7]
    vtype = line[11]
    als = line[22].split(',')
    freqs = line[24].split(',')
    altfreq = 0
    for ii in range(len(freqs)):
        if als[ii] != ref and als[ii] != '':
            altfreq += float(freqs[ii])
    out = line[1:4] + [rsid, vtype, str(altfreq)]
    print '\t'.join(out)
