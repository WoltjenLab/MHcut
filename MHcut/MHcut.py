'''
Main workflow:
- Opening the input files.
- Reading each variant.
- Looking for MH.
- Aligning protospacers.
- Looking for nested MH.
- Writing in the output files.
'''

from pyfaidx import Fasta
import sys
import os
import shutil
import flanks
import pam_utils
import variant
import cartoon
import tqdm
import inDelphi.inDelphi


def mhcut(args):
    '''Main workflow for input arguments "args".'''

    # Open connection to reference genome
    fa_indexed = os.path.isfile(args.reffile + '.fai')
    if not fa_indexed:
        print "The reference file doesn't seem to be indexed."
        print "Indexing now (might take a minute)..."
    reffa = Fasta(args.reffile)
    if not fa_indexed:
        print "Indexing completed."

    if(args.varfile == ''):
        print "Use -ref and -var to run MHcut. "
        print "Run 'MHcut -h' for more info."
        sys.exit(0)

    # Open connection to input file
    variant_input_file = open(args.varfile, 'r')

    # Read input file once to get number of line for progress bar
    input_nb_lines = 0
    for line in variant_input_file:
        input_nb_lines += 1
    # Reopen
    variant_input_file.close()
    variant_input_file = open(args.varfile, 'r')

    # Read first line of input file (header) and write header in output file
    # Change colunm names here.
    # Add/remove columns here but also in the "Write in output files" section
    inhead = variant_input_file.next().rstrip('\n')
    outhead = inhead + '\tvarL'
    outhead += '\t' + '\t'.join(flanks.headers(flank_info=True))
    outhead += '\t' + '\t'.join(pam_utils.headersVariants())
    gouthead = outhead + '\t' + '\t'.join(pam_utils.headersGuides())

    # If restart mode: copy and open unfinished files
    # Also read the header line
    in_nbcols = len(inhead.split('\t'))
    out_nbcols = len(outhead.split('\t'))
    gout_nbcols = len(gouthead.split('\t'))
    restart_mode = False
    test_restart = True
    r_cur_variant_line = ''
    r_cur_guide_line = ''
    if args.restart and os.path.isfile(args.outprefix + '-variants.tsv'):
        print 'Restart mode.'
        shutil.copyfile(args.outprefix + '-variants.tsv',
                        args.outprefix + '-restart-variants.tsv')
        shutil.copyfile(args.outprefix + '-guides.tsv',
                        args.outprefix + '-restart-guides.tsv')
        shutil.copyfile(args.outprefix + '-cartoons.tsv',
                        args.outprefix + '-restart-cartoons.tsv')
        r_variant_outfile = open(args.outprefix + '-restart-variants.tsv')
        r_variant_outfile.next()
        r_cur_variant_line = r_variant_outfile.next()
        r_guide_outfile = open(args.outprefix + '-restart-guides.tsv')
        r_guide_outfile.next()
        r_cur_guide_line = r_guide_outfile.next()
        r_cartoon_outfile = open(args.outprefix + '-restart-cartoons.tsv')
        r_cartoon_outfile.next()
        restart_mode = True

    # Open connection to output files
    variant_outfile = open(args.outprefix + '-variants.tsv', 'w')
    guide_outfile = open(args.outprefix + '-guides.tsv', 'w')
    cartoon_outfile = open(args.outprefix + '-cartoons.tsv', 'w')
    variant_outfile.write(outhead + '\n')
    guide_outfile.write(gouthead + '\n')
    cartoon_outfile.write(outhead + '\n\n')

    # Init inDelphi model if necessary
    idmodels = {}
    if args.indelphi:
        for ct in ['mESC', 'U2OS', 'HEK293', 'HCT116', 'K562']:
            idmodels[ct] = inDelphi.inDelphi.init_model(celltype=ct)

    # Start progress bar
    pbar = tqdm.tqdm(total=input_nb_lines)
    # When to update progress
    if input_nb_lines < 100:
        update_pb = 1
    else:
        update_pb = int(input_nb_lines/100)

    # Read each line of the input file
    line_cpt = 0
    for input_line in variant_input_file:
        # Update progress bar every couple of variants
        line_cpt += 1
        if line_cpt % update_pb == 0:
            pbar.update(update_pb)
        # Parse line from the variant input
        input_line_raw = input_line.rstrip('\n')
        input_line = input_line_raw.split('\t')
        vstart = int(input_line[1])
        vend = int(input_line[2])
        if vend - vstart + 1 < args.minvarL:
            # If variant is too small or too big,
            # skip and jump to next iteration
            continue
        if input_line[0] not in reffa.keys():
            # If chromosome name not in reference, skip variant
            continue
        # Restart mode
        if restart_mode and test_restart:
            r_vline = r_cur_variant_line.rstrip('\n').split('\t')
            r_vline_raw = '\t'.join(r_vline[:in_nbcols])
            try:
                r_cur_variant_line = r_variant_outfile.next()
            except StopIteration:
                test_restart = False
            if input_line_raw == r_vline_raw and len(r_vline) == out_nbcols:
                # Write guides: as long as it matches the variant
                # and not end of file
                r_gline = ''
                r_cgline = r_cur_guide_line.rstrip('\n').split('\t')
                while '\t'.join(r_cgline[:in_nbcols]) == input_line_raw and \
                      len(r_cgline) == gout_nbcols:
                    r_gline += r_cur_guide_line
                    try:
                        r_cur_guide_line = r_guide_outfile.next()
                    except StopIteration:
                        r_cur_guide_line = ''
                        test_restart = False
                    r_cgline = r_cur_guide_line.rstrip('\n').split('\t')
                # Write cartoon: 5 lines
                r_cline = ''
                for ii in range(5):
                    try:
                        r_cline += r_cartoon_outfile.next()
                    except StopIteration:
                        test_restart = False
                # Skip this variant
                if test_restart:
                    variant_outfile.write('\t'.join(r_vline) + '\n')
                    guide_outfile.write(r_gline)
                    cartoon_outfile.write(r_cline)
                    continue
            else:
                # Don't skip and don't test restart files anymore
                test_restart = False
        # Create a variant object with the reference sequence
        var = variant.Variant(input_line[0], vstart, vend, reffa)
        if args.no_opt_shift:
            vars = [var]
        else:
            # List variant when shifting represent the same deletion
            vars = variant.list_shifted_variants(var)
        # To record which shifted variant has the best score
        var = vars[0]
        var_opt_score = -1
        var_fl_1 = var_fl_2 = False
        for var_c in vars:
            # Test MH in flank configuration 1
            var_fl_c_1 = flanks.VarFlank(var_c, flank=1)
            var_fl_c_1.findMH(args.maxConsMM)
            # If no MH or too small, or too low MH ratio or
            # too short first microhomology stretch
            if(not args.nofilter
               and (var_fl_c_1.mhL < args.minMHL or
                    var_fl_c_1.hom < args.minhom or
                    var_fl_c_1.m1L < args.minm1L)):
                var_fl_c_1.score = 0
            # Test MH in flank configuration 1
            var_fl_c_2 = flanks.VarFlank(var_c, flank=2)
            var_fl_c_2.findMH(args.maxConsMM)
            if(not args.nofilter
               and (var_fl_c_2.mhL < args.minMHL or
                    var_fl_c_2.hom < args.minhom or
                    var_fl_c_2.m1L < args.minm1L)):
                var_fl_c_2.score = 0
            # Test if better than previous variants, save if so
            if var_fl_c_2.score > var_opt_score or \
               var_fl_c_1.score > var_opt_score:
                var = var_c
                var_fl_1 = var_fl_c_1
                var_fl_2 = var_fl_c_2
                var_opt_score = max(var_fl_c_1.score, var_fl_c_2.score)

        if args.twofls:
            # Two-flanks mode, i.e. both flank configurations are tested
            # (inner-outer and outer-inner)
            flank_cfgs = [var_fl_1, var_fl_2]
        else:
            # Using the alignment score, the best flank configuration is chosen
            if var_fl_1.score > var_fl_2.score:
                flank_cfgs = [var_fl_1]
            else:
                flank_cfgs = [var_fl_2]

        # Continue analysis on each flank configuration of interest
        for var_fl in flank_cfgs:
            # If a score of 0, either no MH or didn't satisfy criteria above,
            # so skip the rest
            if var_fl.score == 0:
                # If nofilter mode, print the variant output line
                if args.nofilter:
                    voutline = '\t'.join([input_line_raw, str(var.vsize),
                                          var_fl.toString(True)] +
                                         pam_utils.headersVariants(NAs=True))
                    variant_outfile.write(voutline + '\n')
                continue

            # Find PAM motives
            pamseqs = args.pamseq.split(',')
            pams = pam_utils.PAMs(var_fl, pamseqs, cut_offset=args.pamcut,
                                  max_tail=args.maxTail)

            # Align protospacers
            if pams.nbPAMs() > 0:
                if args.jffile == '':
                    if args.bwa:
                        pams.alignPamsBwa(args.reffile, prefix=args.outprefix)
                    else:
                        pams.alignPamsBlast(args.reffile,
                                            prefix=args.outprefix,
                                            chunk_size=args.chunkS)
                else:
                    pams.alignPamsJellyfish(args.jffile, prefix=args.outprefix)

            # We define a unique pam if the protospacer maps exactly
            # at maximum one position in the genome.
            pams.annotateUnique(max_mm0=1)

            # Looking for nested MH only in unique PAMs
            pams.findNestedMH(var, max_tail=args.maxTail,
                              min_l_nmh=args.minLnmh,
                              uniq_pam_only=False)

            if args.indelphi:
                # Running inDelphi
                pams.inDelphi(idmodels, var, uniq_pam_only=False)

            # Prepare variant output line
            voutline = '\t'.join([input_line_raw, str(var.vsize),
                                  var_fl.toString(flank_info=True),
                                  pams.toStringVariants()])

            # Prepare guides output lines
            goutline = pams.toStringGuides(voutline=voutline,
                                           uniq_pam_only=False)

            # Write the cartoon and other files
            cartoon_outfile.write(voutline + '\n')
            cartoon_outfile.write(cartoon.drawCartoon(var, var_fl, pams,
                                                      max_tail=args.maxTail))
            cartoon_outfile.write('\n\n')
            guide_outfile.write(goutline)
            variant_outfile.write(voutline + '\n')
            # variant output last to make sure restart mode works.

    # Update and close progress bar
    if pbar.total > pbar.last_print_n:
        pbar.update(pbar.total - pbar.last_print_n)
    pbar.close()
    print '\nDone.\n'

    # Close connections to input and output files
    variant_input_file.close()
    variant_outfile.close()
    guide_outfile.close()
    cartoon_outfile.close()

    if restart_mode:
        r_variant_outfile.close()
        r_guide_outfile.close()
        r_cartoon_outfile.close()
        os.remove(args.outprefix + '-restart-variants.tsv')
        os.remove(args.outprefix + '-restart-guides.tsv')
        os.remove(args.outprefix + '-restart-cartoons.tsv')
