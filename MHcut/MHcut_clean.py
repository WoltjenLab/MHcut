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
import flanks
import pam_utils
import variant
import cartoon
import tqdm


def mhcut(args):
    '''Main workflow for input arguments "args".'''
    
    # Open connection to reference genome
    print "Check if reference is indexed and index it if not"
    print "...might take a minute..."
    reffa = Fasta(args.reffile)
    print "...Done."

    if(args.varfile == ''):
        print "Use -ref and -var to run MHcut. "
        print "Run 'MHcut -h' for more info."
        sys.exit(0)

    # Open connection to output files
    variant_outfile = open(args.outprefix + '-variants.tsv', 'w')
    guide_outfile = open(args.outprefix + '-guides.tsv', 'w')
    cartoon_outfile = open(args.outprefix + '-cartoons.tsv', 'w')

    # Open connection to input file
    variant_input_file = open(args.varfile, 'r')

    # Read input file once to get number of line for progress bar
    input_nb_lines = 0
    for line in variant_input_file:
        input_nb_lines += 1
    # Reopen
    variant_input_file.close()
    variant_input_file = open(args.varfile, 'r')

    # Read firt line of input file (header) and write header in output file
    # Change colunm names here.
    # Add/remove columns here but also in the "Write in output files" section
    inhead = variant_input_file.next().rstrip('\n')
    outhead = inhead + '\tvarL'
    outhead += '\t' + '\t'.join(flanks.headers(flank_info=args.twofls))
    outhead += '\t' + '\t'.join(pam_utils.headersVariants())
    variant_outfile.write(outhead + '\n')
    gouthead = outhead + '\t' + '\t'.join(pam_utils.headersGuides())
    guide_outfile.write(gouthead + '\n')
    cartoon_outfile.write(outhead + '\n\n')

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
        # Update progress bar every 100 variants
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
        # Create a variant object with the reference sequence
        var = variant.Variant(input_line[0], vstart, vend, reffa)

        # Test MH in flank configuration 1
        var_fl_1 = flanks.VarFlank(var, flank=1)
        var_fl_1.findMH(args.maxConsMM)
        # If no MH or too small, or too low MH ratio or
        # too short first microhomology stretch
        if(not args.nofilter
           and (var_fl_1.mhL < args.minMHL or var_fl_1.hom < args.minhom or
                var_fl_1.m1L < args.minm1L)):
            var_fl_1.score = 0
        # Test MH in flank configuration 1
        var_fl_2 = flanks.VarFlank(var, flank=2)
        var_fl_2.findMH(args.maxConsMM)
        if(not args.nofilter
           and (var_fl_2.mhL < args.minMHL or var_fl_2.hom < args.minhom or
                var_fl_2.m1L < args.minm1L)):
            var_fl_2.score = 0

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
                                          var_fl.toString(args.twofls)] +
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
                    pams.alignPamsBlast(args.reffile, prefix=args.outprefix,
                                        chunk_size=args.chunkS)
                else:
                    pams.alignPamsJellyfish(args.jffile, prefix=args.outprefix)

            # We define a unique pam if the protospacer maps exactly
            # at maximum one position in the genome.
            pams.annotateUnique(max_mm0=1)

            # Looking for nested MH only in unique PAMs
            pams.findNestedMH(var, max_tail=args.maxTail,
                              min_l_nmh=args.minLnmh,
                              uniq_pam_only=True)

            # Write variant output line
            voutline = '\t'.join([input_line_raw, str(var.vsize),
                                 var_fl.toString(args.twofls),
                                  pams.toStringVariants()])
            variant_outfile.write(voutline + '\n')

            # Write guides output lines
            goutline = pams.toStringGuides(voutline=voutline,
                                           uniq_pam_only=True)
            guide_outfile.write(goutline)

            # Write the cartoon
            cartoon_outfile.write(voutline + '\n')
            cartoon_outfile.write(cartoon.drawCartoon(var, var_fl, pams,
                                                      max_tail=args.maxTail))
            cartoon_outfile.write('\n\n')

    # Update and close progress bar
    pbar.update(pbar.maxinterval-pbar.last_print_n)
    pbar.close()
    print '\nDone.\n'

    # Close connections to input and output files
    variant_input_file.close()
    variant_outfile.close()
    guide_outfile.close()
    cartoon_outfile.close()
