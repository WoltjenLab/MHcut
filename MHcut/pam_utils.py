'''
PAM classes with functions to work on PAMs.
'''

import seq_utils
import subprocess
import os
import inDelphi.inDelphi


class PAM():
    '''Contains all the information about one PAM.'''
    def __init__(self, pamseq, cutPosition, strand, proto, m1Dist1,
                 m1Dist2, mhDist1, mhDist2):
        self.pamseq = pamseq
        self.cutPosition = cutPosition
        self.strand = strand
        self.proto = proto
        self.m1Dist1 = m1Dist1
        self.m1Dist2 = m1Dist2
        self.mhDist1 = mhDist1
        self.mhDist2 = mhDist2
        self.uniq = False
        # Protospacer alignment
        self.mm0 = 0
        self.mm1 = 0
        self.mm2 = 0
        # Nested MH
        # Number of nested MH and maximum size (no matter the score)
        self.nmh_nb = 'NA'
        self.nmh_maxL = 'NA'
        # Best nested MH stats
        self.bnmh_score = 'NA'
        self.bnmh_size = 'NA'
        self.bnmh_vsize = 'NA'
        self.bnmh_gc = 'NA'
        self.bnmh_seq = 'NA'
        # inDelphi predictions
        self.indelphi_freq = 'NA'
        self.indelphi_freq_size = 'NA'


class PAMs():
    '''List of PAMs. Useful for bulk processing (protospacer alignment).'''
    def __init__(self, var_fl, pam_seqs, cut_offset, max_tail=50):
        # The list of PAM objects
        self.pams = []
        # Look for PAMs for each PAM sequence
        for pamseq in pam_seqs:
            self.findPAMs(var_fl, pamseq, cut_offset, max_tail)
        # Global stats
        self.nb_pam_motives = 0
        # Number of guides with no off target MH (now called nested MH)
        self.no_offtargets = 'NA'
        self.min_offtargets = 'NA'
        # inDelphi stats at the variant level
        self.max_indelphi_freq = 'NA'
        self.max_indelphi_freq_size = 'NA'

    def findPAMs(self, var_fl, pamseq, cut_offset, max_tail):
        '''Looks for PAM for a PAM sequence and cut offset.'''
        # Also look for the reverse complement
        pamseq_rev = seq_utils.revComp(pamseq)
        # Full sequence
        seq = var_fl.var.fl1seq + var_fl.var.varseq + var_fl.var.fl2seq
        # Minimum (exact) MH length that must remains after a cut.
        min_mh_cut = min(var_fl.m1L, 3)
        fl_L = len(var_fl.var.fl1seq)
        var_L = len(var_fl.var.varseq)
        # Where we should look for cuts
        search_range = [fl_L - 1, fl_L + var_L - min_mh_cut]
        # Interval between the full MHs
        reduced_search_range = [fl_L - 1,
                                fl_L + var_L - var_fl.mhL]
        if var_fl.flank == 2:
            search_range = [fl_L + min_mh_cut - 1,
                            fl_L + var_L]
            reduced_search_range = [fl_L + var_fl.mhL - 1,
                                    fl_L + var_L]
        # Test each position: if it matched the motif and
        # in the search range, add to list
        for pos in xrange(len(seq)-len(pamseq)+1):
            proto_seq = strand = cut_pos = False
            if seq_utils.aligned(seq[pos:pos+len(pamseq)], pamseq):
                cut_pos = pos + cut_offset - 1
                strand = '+'
            if seq_utils.aligned(seq[pos:pos+len(pamseq)], pamseq_rev):
                cut_pos = pos - cut_offset + len(pamseq_rev) - 1
                strand = '-'
            # If there is a match in the correct interval
            if(strand and cut_pos >= search_range[0]
               and cut_pos < search_range[1]
               and (cut_pos < search_range[0] + max_tail or
                    cut_pos >= search_range[1] - max_tail)):
                # Get protospacer sequence
                if strand == '+':
                    proto_start = cut_pos - 19 - cut_offset
                    proto_end = cut_pos - cut_offset + 1
                else:
                    proto_start = cut_pos + cut_offset + 1
                    proto_end = cut_pos + 20 + cut_offset + 1
                proto_seq = seq[proto_start:proto_end]
                # Info about the PAM found, including
                # distance to MH on each side using first stretch of
                # perfect match or the extended MH
                pam = PAM(pamseq=pamseq, cutPosition=cut_pos,
                          strand=strand, proto=proto_seq,
                          m1Dist1=cut_pos - search_range[0],
                          m1Dist2=search_range[1] - cut_pos - 1,
                          mhDist1=cut_pos - reduced_search_range[0],
                          mhDist2=reduced_search_range[1] - cut_pos - 1)
                self.pams.append(pam)

    def nbPAMs(self):
        '''Number of PAMs in the list.'''
        return(len(self.pams))

    def getMax2cutsDist(self, uniq_pam_only=True):
        '''Computes the maximum distance between 2 valid cuts'''
        mincut = maxcut = ''
        nb_pams = 0
        for pam in self.pams:
            if uniq_pam_only and not pam.uniq:
                continue
            else:
                nb_pams += 1
            if mincut == '':
                mincut = maxcut = pam.cutPosition
            else:
                mincut = min(mincut, pam.cutPosition)
                maxcut = max(maxcut, pam.cutPosition)
        # If at least two cuts considered
        if nb_pams > 1:
            return(maxcut - mincut)
        else:
            # Otherwise return NA
            return('NA')

    def annotateUnique(self, max_mm0=1, max_mm1=float('inf'),
                       max_mm2=float('inf')):
        '''
        Annotate each PAM in the list according to their protospacer
        alignment.

        This is where to define how unique the protospacer must be
        By default max_mm0=1, i.e. there must be only one position
        in the genome aligning perfectly.
        '''
        for pam in self.pams:
            # Handle NAs (e.g. when using JellyFish)
            if pam.mm1 == 'NA':
                mm1 = 0
            else:
                mm1 = pam.mm1
            if pam.mm2 == 'NA':
                mm2 = 0
            else:
                mm2 = pam.mm2
            # Test if unique
            if(pam.mm0 > 0 and pam.mm0 <= max_mm0 and
               mm1 <= max_mm1 and mm2 <= max_mm2):
                pam.uniq = True

    def alignPamsBlast(self, reffile, include_pam=True, prefix='blast',
                       chunk_size=30):
        '''Align protospacers and update the PAMs.'''
        # Chunk PAMs to avoid memory explosion
        pams_chunks = []
        for ii in range(0, len(self.pams), chunk_size):
            pams_chunks.append(self.pams[ii:(ii + chunk_size)])
        fasta_file = prefix + '_tempMHcut.fasta'
        for pams in pams_chunks:
            ff = open(fasta_file, 'w')
            pams_hash = {}
            for pam in pams:
                if 'N' in pam.proto:
                    pam.mm0 = 'NA'
                    pam.mm1 = 'NA'
                    pam.mm2 = 'NA'
                    continue
                if include_pam:
                    if pam.strand == '+':
                        protoguide = pam.proto + pam.pamseq
                    else:
                        protoguide = seq_utils.revComp(pam.pamseq) + pam.proto
                else:
                    protoguide = pam.proto
                protoguides = seq_utils.enumN(protoguide)
                for ii in xrange(len(protoguides)):
                    pamid = '{}_{}_{}'.format(pam.cutPosition, pam.strand, ii)
                    pams_hash[pamid] = pam
                    ff.write('>' + pamid + '\n' + protoguides[ii] + '\n')
            ff.close()
            dump = open('/dev/null')
            blast_cmd = ['blastn', '-db', reffile,  '-query', fasta_file,
                         '-outfmt', '6', '-word_size', '10',
                         '-max_target_seqs', '20']
            blast_out = subprocess.check_output(blast_cmd, stderr=dump)
            dump.close()
            blast_out = blast_out.split('\n')
            for line in blast_out:
                line = line.split('\t')
                if len(line) > 1:
                    pam_cand = pams_hash[line[0]]
                    # If that the blast hit is as long as the full sequence
                    # and has no gap.
                    protoguide_length = len(pam_cand.proto)
                    if include_pam:
                        protoguide_length += len(pam_cand.pamseq)
                    if(int(line[3]) == protoguide_length and
                       int(line[5]) == 0):
                        # Update appropriate mismatch count
                        if int(line[4]) == 0:
                            pam_cand.mm0 += 1
                        if int(line[4]) == 1:
                            pam_cand.mm1 += 1
                        if int(line[4]) == 2:
                            pam_cand.mm2 += 1
            os.remove(fasta_file)

    def alignPamsJellyfish(self, jffile, include_pam=True, prefix='jf'):
        '''Align protospacers and update the PAMs.'''
        pams_hash = {}
        # Temporary file with the protospacer sequence to align
        fasta_file = prefix + '_tempMHcut.fasta'
        ff = open(fasta_file, 'w')
        cpt = 0
        for pam in self.pams:
            # Init PAM values
            pam.mm0 = 0
            pam.mm1 = 'NA'
            pam.mm2 = 'NA'
            # If Ns in the protospacer, skip
            if 'N' in pam.proto:
                pam.mm0 = 'NA'
                continue
            # Should we look for protospacer + PAM
            if include_pam:
                if pam.strand == '+':
                    protoguide = pam.proto + pam.pamseq
                else:
                    protoguide = seq_utils.revComp(pam.pamseq) + pam.proto
            else:
                protoguide = pam.proto
            # List sequence to remove the Ns (from the PAMs)
            protoguides = seq_utils.enumN(protoguide)
            for pg in protoguides:
                # Write the sequence
                ff.write('>' + str(cpt) + '\n' + pg + '\n')
                cpt += 1
                # Save the corresponding PAM object
                if pg in pams_hash:
                    pams_hash[pg].append(pam)
                else:
                    pams_hash[pg] = [pam]
                # Same for the reverse complement
                pg_rc = seq_utils.revComp(pg)
                if pg_rc in pams_hash:
                    pams_hash[pg_rc].append(pam)
                else:
                    pams_hash[pg_rc] = [pam]
        ff.close()
        # If we wrote some sequence to align, run JellyFish
        if cpt > 0:
            # JellyFish command
            jellyfish_cmd = ['jellyfish', 'query', '-L', jffile,
                             '-s', fasta_file]
            dump = open('/dev/null')
            jellyfish_out = subprocess.check_output(jellyfish_cmd, stderr=dump)
            dump.close()
            # Parse jellyfish output
            jellyfish_out = jellyfish_out.rstrip('\n').split('\n')
            for line in jellyfish_out:
                line = line.split(' ')
                # For each protoguide, update the count in corresponding PAM
                for pam in pams_hash[line[0]]:
                    pam.mm0 += int(line[1])
        # Remove temporary file
        os.remove(fasta_file)

    def inDelphi(self, var, uniq_pam_only=False):
        '''Run inDelphi to predict repair outcome.'''
        pams = []
        # If we want only unique PAMs, retrieve them
        if uniq_pam_only:
            for pam in self.pams:
                if pam.uniq:
                    pams.append(pam)
        else:
            # Otherwise analyze at all the PAMs
            pams = self.pams
        # Prepare input and target sequences
        full_seq = var.fl1seq + var.varseq + var.fl2seq
        target_seq = var.fl1seq + var.fl2seq
        full_seq_rc = seq_utils.revComp(full_seq)
        target_seq_rc = seq_utils.revComp(target_seq)
        # If there is a N in the sequence, inDelphi raises an error.
        # Let's skip these rare variants
        if 'N' in full_seq:
            return()
        # Otherwise, we init the info and go over the cuts
        if len(pams) > 0:
            # Init frequencies to 0
            self.max_indelphi_freq = 0
            self.max_indelphi_freq_size = 0
        # For each PAM, run inDelphi and update PAM/variant stats
        for pam in pams:
            # Init frequencies to 0
            pam.indelphi_freq = 0
            pam.indelphi_freq_size = 0
            # Prepare cut position and sequences
            if pam.strand == '-':
                cut_pos = len(full_seq) - pam.cutPosition - 1
                delphi_input = full_seq_rc
                delphi_comp = target_seq_rc
            else:
                cut_pos = pam.cutPosition + 1
                delphi_input = full_seq
                delphi_comp = target_seq
            # Run inDelphi
            pred_df, stats = inDelphi.inDelphi.predict(delphi_input, cut_pos)
            pred_df = inDelphi.inDelphi.add_genotype_column(pred_df, stats)
            # Loop over prediction looking for target sequence/size
            for row in pred_df.iterrows():
                if row[1]['Genotype'] == delphi_comp:
                    pam.indelphi_freq = round(row[1]['Predicted frequency'], 3)
                if(row[1]['Length'] == var.vsize and
                   row[1]['Genotype position'] != 'e'):
                    pam.indelphi_freq_size += row[1]['Predicted frequency']
            pam.indelphi_freq_size = round(pam.indelphi_freq_size, 3)
            # Update variant-level stats
            self.max_indelphi_freq_size = max(self.max_indelphi_freq_size,
                                              pam.indelphi_freq_size)
            self.max_indelphi_freq = max(self.max_indelphi_freq,
                                         pam.indelphi_freq)

    def findNestedMH(self, var, max_tail=50, min_l_nmh=3,
                     uniq_pam_only=False):
        '''Look for nested MH for each valid cut.'''
        # Search for other MH that could be used by the MMEJ
        pams = []
        # If we want only unique PAMs, retrieve them
        if uniq_pam_only:
            for pam in self.pams:
                if pam.uniq:
                    pams.append(pam)
        else:
            # Otherwise analyze at all the PAMs
            pams = self.pams
        # Look for nested MH on variants that are small enough
        if len(pams) > 0 and var.vsize < max_tail*2:
            full_seq = var.fl1seq + var.varseq + var.fl2seq
            other_mh = seq_utils.RegionExactMH(full_seq)
            # Go over PAMs and look for nested MH
            for pam in pams:
                # Init PAM values
                pam.nmh_nb = 0
                pam.nmh_maxL = 0
                pam.bnmh_score = 0
                # List nested MH for this cut position
                mhhet = other_mh.listmh(pam.cutPosition, min_l_nmh)
                for mho in mhhet:
                    data = mhhet[mho]
                    # Only consider other MH that are at least as close
                    # from each other as our target MH.
                    v_fl_size = var.flsize + var.vsize
                    if(data['vsize'] <= var.vsize and
                       (data['startU'] != var.flsize or
                        data['startD'] != v_fl_size)
                       and (data['startU'] + data['size'] != var.flsize or
                            data['startD'] + data['size'] != v_fl_size)):
                        pam.nmh_nb += 1
                        pam.nmh_maxL = max(pam.nmh_maxL, data['size'])
                        # Check and saves info about the strongest nested MH
                        if data['score'] > pam.bnmh_score:
                            pam.bnmh_score = data['score']
                            pam.bnmh_size = data['size']
                            pam.bnmh_vsize = data['vsize']
                            pam.bnmh_gc = round(data['gc'], 3)
                            pam.bnmh_seq = data['seq']
                # Update global info
                # Number of PAMs with no nested MH
                if self.no_offtargets == 'NA':
                    self.no_offtargets = 0
                if pam.nmh_nb == 0:
                    self.no_offtargets += 1
                # Number of nested MH in the PAM with the least
                if self.min_offtargets == 'NA':
                    self.min_offtargets = pam.nmh_nb
                else:
                    self.min_offtargets = min(self.min_offtargets, pam.nmh_nb)

    def toStringVariants(self):
        '''The relevant string to write in the "variants" output.'''
        # IF YOU CHANGE SOMETHING HERE, CHANGE THE headersVariants TOO (below)
        # Count number of unique PAMs
        pams_uniq = 0
        for pam in self.pams:
            if pam.uniq:
                pams_uniq += 1
        # Create string to output
        tostr = [self.nbPAMs(), pams_uniq, self.no_offtargets,
                 self.min_offtargets, self.getMax2cutsDist(),
                 self.max_indelphi_freq, self.max_indelphi_freq_size]
        tostr = '\t'.join([str(ii) for ii in tostr])
        return tostr

    def toStringGuides(self, voutline, uniq_pam_only=True):
        '''The relevant string to write in the "guides" output.'''
        # IF YOU CHANGE SOMETHING HERE, CHANGE THE headersGuides TOO (below)
        tostr = ''
        for pam in self.pams:
            # Skip non-unique PAMs
            if uniq_pam_only and not pam.uniq:
                continue
            # Create string to output for this guide
            pam_str = [pam.proto, pam.pamseq, pam.mm0, pam.mm1, pam.mm2,
                       pam.m1Dist1, pam.m1Dist2, pam.mhDist1,  pam.mhDist2,
                       pam.nmh_nb, pam.nmh_maxL, pam.bnmh_score, pam.bnmh_size,
                       pam.bnmh_vsize, pam.bnmh_gc, pam.bnmh_seq,
                       pam.indelphi_freq, pam.indelphi_freq_size]
            pam_str = '\t'.join([str(ii) for ii in pam_str])
            tostr += voutline + '\t' + pam_str + '\n'
        return tostr


def headersVariants(NAs=False):
    '''
    The headers corresponding to the "variants" output from the
    "toStringVariants" output.
    '''
    # IF YOU CHANGE SOMETHING HERE, CHANGE THE toStringVariants TOO (above)
    headers = ['pamMot', 'pamUniq', 'guidesNoNMH', 'guidesMinNMH',
               'max2cutsDist', 'maxInDelphiFreqDel', 'maxInDelphiFreqSize']
    # Should it returns NAs instead of the headers
    # (useful when skipping variants)
    if NAs:
        return(['NA' for ii in range(len(headers))])
    else:
        return(headers)


def headersGuides():
    '''
    The headers corresponding to the "variants" output from the
    "toStringGuides" output.
    '''
    # IF YOU CHANGE SOMETHING HERE, CHANGE THE toStringGuides TOO (above)
    return ['protospacer', 'pamSeq', 'mm0', 'mm1', 'mm2', 'm1Dist1', 'm1Dist2',
            'mhDist1', 'mhDist2', 'nbNMH', 'largestNMH', 'nmhScore', 'nmhSize',
            'nmhVarL', 'nmhGC', 'nmhSeq', 'inDelphiFreqDel',
            'inDelphiFreqSize']
