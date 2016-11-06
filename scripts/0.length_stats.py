# Print read length percentiles (Min, 5%, 10%, 15%, 20%, 25%, Max)

import argparse, util
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-f', help = 'Input FASTA file', default = '')
parser.add_argument('-q', help = 'Input FASTQ file', default = '')
parser.add_argument('-n', help = 'Number of sequences to evaluate', default = 100000, type = int)
args = parser.parse_args()

if args.f != '':
    fn = args.f
    iter_seq = util.iter_fst

if args.q != '':
    fn = args.q
    iter_seq = util.iter_fsq

x = []
for record in iter_seq(fn):
    sid, seq = record[:2]
    if args.n > 0:
        args.n = args.n - 1
    else:
        break
    x.append(len(seq))

x = np.array(x)
print '\nMin = %d' %(min(x))
print 'Max = %d\n' %(max(x))
print '95%% = %d' %(np.percentile(x, 95))
print '90%% = %d' %(np.percentile(x, 90))
print '75%% = %d' %(np.percentile(x, 75))
print '50%% = %d' %(np.percentile(x, 50))
print '5%% = %d\n' %(np.percentile(x, 5))
