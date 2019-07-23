import argparse
import sys
import MHcut
import MHcut_old


def main():
    # Define arguments
    parser = argparse.ArgumentParser(
        description='Find regions with microhomology and a cut position.')
    parser.add_argument('-ref', dest='reffile', required=True,
                        help='the reference genome fasta file. Required.')
    parser.add_argument('-var', dest='varfile', default='', required=True,
                        help='the file with the variants location '
                        '(BED-TSV format with header).')
    parser.add_argument('-jf', dest='jffile', default='', required=True,
                        help='the jellyfish file of the reference genome.')
    parser.add_argument('-out', dest='outprefix', default='mhcutres',
                        help='the prefix for the output files. '
                        'Default: "mhcutres".')
    parser.add_argument('-2fls', dest='twofls', action='store_true',
                        help='report results for both flank configurations '
                        'instead of the one with strongest MH.')
    parser.add_argument('-noShift', dest='no_opt_shift',
                        action='store_true',
                        help='use input coordinates without trying to '
                        'shift the variant to find the best MH')
    parser.add_argument('-v0', dest='v0', action='store_true',
                        help='use older version.')
    parser.add_argument('-restart', dest='restart', action='store_true',
                        help='restart mode. Will skip variants already in'
                        ' output and append new output. Works only with '
                        'nofilt mode. It does not work if -2fls is used.')

    # Parameters related to protospacer alignment
    align_pars = parser.add_argument_group('Protospacer alignment')
    align_pars.add_argument('-chunkS', type=int, default=30,
                            help='if using BLAST, the number of PAMs per '
                            'chunks. Default: 30. Pick lower value if BLAST '
                            'uses too much memory.')
    align_pars.add_argument('-bwa', dest='bwa',
                            action='store_true',
                            help='use BWA for protospacer alignment')

    # Parameters about filtering
    filt_pars = parser.add_argument_group('Filters')
    filt_pars.add_argument('-minvarL', dest='minvarL', default=3, type=int,
                           help='the minimum variant length. Default: 3.')
    filt_pars.add_argument('-minMHL', dest='minMHL', default=3, type=int,
                           help='the minimum microhomology length. '
                           'Default: 3.')
    filt_pars.add_argument('-minm1L', dest='minm1L', default=3, type=int,
                           help='the minimum length of first.'
                           'microhomology stretch. Default: 3.')
    filt_pars.add_argument('-minhom', dest='minhom', default=0, type=float,
                           help='the minimum homology ratio. Default: 0.')
    filt_pars.add_argument('-nofilt', dest='nofilter', action='store_true',
                           help="Don't filter variants without MH. Ignore "
                           "-minMHL, -minhom and -minm1L.")

    # Parameters about the PAMs
    pam_pars = parser.add_argument_group('PAM specification')
    pam_pars.add_argument('-PAM', dest='pamseq', default='NGG',
                          help='the PAM motif. Default is NGG. '
                          'Possibly several separated by ",".')
    pam_pars.add_argument('-PAMcut', dest='pamcut', default=-3, type=int,
                          help='the cut position relative to the PAM motif. '
                          'Default is -3.')
    pam_pars.add_argument('-indelphi', dest='indelphi', action='store_true',
                          help="Run inDelphi to predict repair for each PAM.")

    # Parameters for MH search
    mh_pars = parser.add_argument_group('Micro-homology specification')
    mh_pars.add_argument('-maxConsMM', dest='maxConsMM', default=1, type=int,
                         help='the maximum number of consecutive mismatches '
                         'allowed when extending the MH. Default: 1.')
    mh_pars.add_argument('-maxTail', dest='maxTail', default=50, type=int,
                         help='the maximum hanging tail allowed. Default: 50.')
    mh_pars.add_argument('-minLnhm', dest='minLnmh', default=3, type=int,
                         help='the minimum length of the nested '
                         'microhomology. Default: 3.')

    # Parse arguments
    args = parser.parse_args(sys.argv[1:])

    # Warnings that filtering arguments are not used in nofilter mode
    if args.nofilter:
        print('no filter mode (-nofilt): all variants will be kept and the ' +
              'following parameters will NOT be taken into account: ' +
              '-minMHL, -minhom, -minm1L')

    if args.v0:
        MHcut_old.mhcut(args)
    else:
        MHcut.mhcut(args)


if __name__ == "__main__":
    main()
    # try:
    #     main()
    # except Exception as e:
    #     print e.message
    #     sys.exit(1)
