'''
Align guides and count number hits with 0, 1, 2 mismatches.
'''

import h5py
import seq_utils
from Bio import SeqIO
import tqdm


class SeededGuides:
    def __init__(self, file_name='seededguides.hdf5', k=23, ns=3, sl=7,
                 PAM='NGG'):
        self.file_name = file_name
        self.ns = ns
        self.sl = sl
        self.k = k
        self.PAM = PAM

    def initFile(self):
        # Init HDF5 file
        hf = h5py.File(self.file_name, 'w')
        hf.create_group('seeds')
        hf.create_dataset('seqs', (0, 1), dtype='S23', maxshape=(None, 1))
        hf.close()

    def scanRef(self, ref_file, verbose=True,
                chrs=['chr' + str(ii) for ii in range(1, 23) + ['X', 'Y']]):
        '''Scan a reference FASTA file, by chr and chunks.'''
        for record in SeqIO.parse(ref_file, "fasta"):
            if record.id in chrs:
                if verbose:
                    print 'Scanning ' + record.id + '...'
                # split into chunks if larger than 100 Mbp
                chunk_size = 100000000
                if len(record.seq) > chunk_size:
                    nb_chunks = (len(record.seq) / chunk_size) + 1
                    chunk_size = int(len(record.seq)/nb_chunks)
                    for chunk in range(nb_chunks):
                        chunk_start = chunk*chunk_size
                        if chunk == nb_chunks - 1:
                            chunk_end = -1
                        else:
                            chunk_end = (chunk+1)*chunk_size
                        if verbose and nb_chunks > 1:
                            print 'Chunk ' + str(chunk) + '...'
                        self.scanChunk(record.seq, verbose=verbose,
                                       idx_start=chunk_start,
                                       idx_end=chunk_end)

    def scanChunk(self, seq, idx_start=0, idx_end=-1, verbose=True):
        '''Scans a sequence chunk between start and end indexes.'''
        if verbose:
            # Start progress bar
            pbar = tqdm.tqdm(total=len(seq))
            # When to update progress
            if len(seq) < 100:
                update_pb = 1
            else:
                update_pb = int(len(seq)/100)
        # if not specified scan to the end of the seq
        if idx_end == -1:
            idx_end = len(seq)
        # adjust the last position to consider
        if idx_end > len(seq)-self.k:
            idx_end = len(seq)-self.k
        # prepare objects and open hdf5 file
        seqs = []
        seeds = {}
        PAMrc = seq_utils.revComp(self.PAM)
        hf = h5py.File(self.file_name, 'a')
        hfseqs = hf['seqs']
        cpt = hfseqs.shape[0]
        # Loop over kmers
        line_cpt = 0
        for ii in xrange(idx_start, idx_end):
            seq_ii = seq[ii:(ii+self.k)]
            seq_ii = str(seq_ii).upper()
            # Discard kmer with a N in the reference
            if 'N' not in seq_ii:
                # align the PAM in both strands
                al = False
                if seq_utils.aligned(seq_ii[(self.k-3):], self.PAM):
                    al = True
                elif seq_utils.aligned(seq_ii[:3], PAMrc):
                    seq_ii = seq_utils.revComp(seq_ii)
                    al = True
                # if aligned extract seeds for this sequence
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
            # Update progress bar every couple of variants
            line_cpt += 1
            if verbose and line_cpt % update_pb == 0:
                pbar.update(update_pb)
        if verbose:
            # Update and close progress bar
            if pbar.total > pbar.last_print_n:
                pbar.update(pbar.total - pbar.last_print_n)
            pbar.close()
        # Update sequences in HDF5 file
        cpt = hfseqs.shape[0]
        hfseqs.resize((hfseqs.shape[0] + len(seqs), 1))
        hfseqs[cpt:, 0] = seqs
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
                    if ss in seq_cands:
                        seq_cands[ss] += 1
                    else:
                        seq_cands[ss] = 1
        # Keep only positions with 2 seeds
        seqs2 = []
        for seqid in seq_cands:
            if seq_cands[seqid] > 1:
                seqs2.append(seqid)
        # Align and compute number of mismatches
        nm_sum = [0, 0, 0]
        hfseqs = hf['seqs']
        for seqid in seqs2:
            seq = str(hfseqs[seqid][0])
            nm = 0
            for ii in range(len(seq)):
                if seq[ii] != inseq[ii]:
                    nm += 1
                    if nm > 2:
                        break
            if nm == 0:
                nm_sum[0] += 1
            if nm == 1:
                nm_sum[1] += 1
            if nm == 2:
                nm_sum[2] += 1
        hf.close()
        return(nm_sum)


# # FOR TESTING
# sg = SeededGuides()
# sg.scanRef('../data/smallref.fa')

# sg.querySeq('ACCGTAAGATACATTTACAGCGG')
# sg.querySeq('ATCGTAAGATACATTTACAGGGG')

# hf = h5py.File(sg.file_name, 'a')
# hf.keys()
# hf['seqs'].shape
# hf['seqs'][0][0]
# len(hf['seeds'].keys())
# hf.close()
