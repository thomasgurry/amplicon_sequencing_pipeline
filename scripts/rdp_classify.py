#!/usr/bin/env python
#
# author: scott olesen <swo@mit.edu>
# modified by: thomas gurry <tg296@mit.edu>

import argparse, subprocess

if __name__ == '__main__':
    p = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('fasta', help='input fasta')
    p.add_argument('output', help='output tab-delimited file')
    p.add_argument('--gene', '-g', choices=['16srrna', 'fungalits_unite'], default='16srrna', help='Marker gene')
    p.add_argument('--format', '-f', choices=['fixrank', 'allrank'], default='fixrank', help='output format')
    p.add_argument('--jar', '-j', default='/home/ubuntu/tools/RDPTools/classifier.jar', metavar='/path/to/classifier.jar', help='path to RDP classifier jar')
    p.add_argument('--java', '-x', default='/usr/bin/java', metavar='/path/to/java', help='java binary')
    p.add_argument('--verbose', '-v', action='store_true', help='show commands being executed?')
    args = p.parse_args()

    cmd = [args.java, '-Xmx1g', '-jar', args.jar, 'classify', '-o', args.output, args.fasta, '-g', args.gene]

    if args.verbose:
        print(" ".join(cmd))

    subprocess.check_call(cmd)
