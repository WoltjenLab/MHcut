'''
Align guides and count number hits with 0, 1, 2 mismatches.
'''

import h5py
import seq_utils
from Bio import SeqIO


class SeededGuides:
    def __init__(self, file_name='seededguides.hdf5', k=23, ns=3, sl=7,
                 PAM='NGG'):
        self.file_name = file_name
        self.ns = ns
        self.sl = sl
        self.k = k
        self.PAM = PAM
        # Init HDF5 file
        hf = h5py.File(self.file_name, 'w')
        hf.create_group('seeds')
        hf.create_dataset('seqs', (1, 0), dtype='S23', maxshape=(1, None))
        hf.close()

    def scanRef(self, ref_file,
                chrs=['chr' + str(ii) for ii in range(1, 23) + ['X', 'Y']]):
        for record in SeqIO.parse(ref_file, "fasta"):
            if record.id in chrs:
                self.scanChunk(record.seq, record.id)

    def scanChunk(self, seq, chr_name, idx_start=0, idx_end=-1):
        if idx_end == -1:
            idx_end = len(seq)
        seqs = []
        seeds = {}
        PAMrc = seq_utils.revComp(self.PAM)
        hf = h5py.File(self.file_name, 'a')
        hfseqs = hf['seqs']
        cpt = hfseqs.shape[1]
        for ii in xrange(len(seq)-self.k):
            seq_ii = seq[ii:(ii+self.k)]
            seq_ii = str(seq_ii).upper()
            if 'N' not in seq_ii:
                al = False
                if seq_utils.aligned(seq_ii[(self.k-3):], self.PAM):
                    al = True
                elif seq_utils.aligned(seq_ii[:3], PAMrc):
                    seq_ii = seq_utils.revComp(seq_ii)
                    al = True
                if al:
                    seqs.append(seq_ii)
                    for seedi in range(self.ns):
                        sstart = (seedi*self.k/self.ns)
                        seed = seq_ii[sstart:(sstart+self.sl)]
                        if seed not in seeds:
                            seeds[seed] = [cpt]
                        else:
                            seeds[seed].append(cpt)
                    cpt += 1
        # Update sequences in HDF5 file
        cpt = hfseqs.shape[1]
        hfseqs.resize((1, hfseqs.shape[1] + len(seqs)))
        hfseqs[0, cpt:] = seqs
        # Update seeds in HDF5 file
        for seed in seeds:
            seed_pos = seeds[seed]
            if 'seeds/' + seed in hf:
                hfs = hf['seeds/' + seed]
                cpt = hfs.shape[1]
                hfs.resize((1, hfs.shape[1] + len(seed_pos)))
                hfs[0][cpt:] = seed_pos
            else:
                hf.create_dataset('seeds/' + seed, dtype=int,
                                  shape=(1, len(seed_pos)), data=seed_pos,
                                  maxshape=(1, None))
        hf.close()

    def querySeq(self, inseq):
        hf = h5py.File(self.file_name, 'r')
        # Unique positions of candidate sequences
        seq_cands = {}
        for seedi in range(self.ns):
            sstart = (seedi*self.k/self.ns)
            seed = inseq[sstart:(sstart+self.sl)]
            if 'seeds/' + seed in hf:
                for ss in hf['seeds/' + seed][0]:
                    seq_cands[ss] = True
        # Align and compute number of mismatches
        nm_sum = [0, 0, 0]
        hfseqs = hf['seqs'][0]
        for seqid in seq_cands:
            seq = str(hfseqs[seqid])
            nm = 0
            for ii in range(len(seq)):
                if seq[ii] != inseq[ii]:
                    nm += 1
            if nm == 0:
                nm_sum[0] += 1
            if nm == 1:
                nm_sum[1] += 1
            if nm == 2:
                nm_sum[2] += 1
        hf.close()
        return(nm_sum)


# Optimization
# - Access position of arrays using h5py
# - Auto-resize array when adding sequences


sg = SeededGuides()
sg.scanRef('/home/jmonlong/Documents/genomecenter/mhcut/ref.fa')

sg.querySeq('GATATGCCTCCCCCTACCTTGGG')
sg.querySeq('GATATGCATCCTCCTACCTTGGG')

hf = h5py.File(sg.file_name, 'a')
hf.keys()
hf['seqs'].shape
hf['seqs'][0][0]
len(hf['seeds'].keys())
hf.close()
