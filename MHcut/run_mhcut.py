import argparse
import sys
import MHcut
import MHcut_clean


def main():
        # Define arguments
    parser = argparse.ArgumentParser(
        description='Find regions with microhomology and a cut position.')
    parser.add_argument('-ref', dest='reffile', required=True,
                        help='the reference genome fasta file')
    parser.add_argument('-var', dest='varfile', default='',
                        help='the file with the variants location '
                        '(BED-TSV format with header)')
    parser.add_argument('-out', dest='outprefix', default='',
                        help='the prefix for the output files')
    parser.add_argument('-jf', dest='jffile', default='',
                        help='the jellyfish file of the reference genome')
    parser.add_argument('-minvarL', dest='minvarL', default=3, type=int,
                        help='the minimum variant length')
    parser.add_argument('-minMHL', dest='minMHL', default=3, type=int,
                        help='the minimum microhomology length')
    parser.add_argument('-maxConsMM', dest='maxConsMM', default=1, type=int,
                        help='the maximum number of consecutive mismatches' +
                             'allowed when extending the MH.')
    parser.add_argument('-maxTail', dest='maxTail', default=50, type=int,
                        help='the maximum hanging tail allowed')
    parser.add_argument('-minhom', dest='minhom', default=0, type=float,
                        help='the minimum homology ratio')
    parser.add_argument('-minm1L', dest='minm1L', default=3, type=int,
                        help='the minimum length of first'
                        'microhomology stretch')
    parser.add_argument('-PAM', dest='pamseq', default='NGG',
                        help='the PAM motif. Possibly several'
                        'separated by ","')
    parser.add_argument('-PAMcut', dest='pamcut', default=-3, type=int,
                        help='the cut position relative to the PAM motif')
    parser.add_argument('-minLnhm', dest='minLnmh', default=3, type=int,
                        help='the minimum length of the nested microhomology')
    parser.add_argument('-nofilt', dest='nofilter', action='store_true',
                        help="Don't filter variants without MH.")
    parser.add_argument('-2fls', dest='twofls', action='store_true',
                        help='Report results for both flank configuration'
                        'instead of the one with highest MH.')
    parser.add_argument('-v2', dest='v2', action='store_true',
                        help='clean version')
    args = parser.parse_args(sys.argv[1:])

    if(args.nofilter):
        print('no filter mode (-nofilt): all variants will be kept and the ' +
              'following parameters will NOT be taken into account: ' +
              '-minMHL, -minhom, -minm1L')

    if(args.v2):
        MHcut_clean.mhcut(args)
    else:
        MHcut.mhcut(args)


if __name__ == "__main__":
    main()
    # try:
    #     main()
    # except Exception as e:
    #     print e.message
    #     sys.exit(1)
