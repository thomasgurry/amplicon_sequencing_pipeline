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
print '5%% = %d' %(np.percentile(x, 5))
print '10%% = %d' %(np.percentile(x, 10))
print '15%% = %d' %(np.percentile(x, 15))
print '20%% = %d' %(np.percentile(x, 20))
print '25%% = %d\n' %(np.percentile(x, 25))
