# Convert a derep mapping file to an OTU table

import pandas as pd
import argparse
import numpy as np
import util

parser = argparse.ArgumentParser()
parser.add_argument('--fst', help='Input fasta sequences (optional)', default='')
parser.add_argument('--map', help='Input mapping file', required=True)
parser.add_argument('--min_count', help='Minimum read count', type=int)
parser.add_argument('--min_samples', help='Minimum number of samples', type=int)
parser.add_argument('--out', help='Output counts matrix')

# Parse command line arguments
args = parser.parse_args()

# Load valid fst seqs
keep = {}
if args.fst:
    for [otu, seq] in util.iter_fst(args.fst):
        otu = otu[1:]
        keep[otu] = 1

# Keep track of samples and otus
samples = {}
otus = {}

# For every line in the mapping file
for line in open(args.map):
    # Load otu name and table of sample counts
    otu, table = line.rstrip().split('\t')
    if len(keep) > 0 and otu not in keep:
        continue
    entries = table.split(' ')
    count = sum([int(entry.split(':')[1]) >= args.min_count for entry in entries])
    if count < args.min_samples:
        continue
    if otu not in otus:
        otus[otu] = len(otus)
    for entry in entries:
        sample, count = entry.split(':')
        if sample not in samples:
            samples[sample] = len(samples)

x = np.zeros([len(samples), len(otus)])

for line in open(args.map):
    otu, table = line.rstrip().split('\t')
    if len(keep) > 0 and otu not in keep:
        continue
    if otu in otus:
        for entry in table.split(' '):
            sample, count = entry.split(':')
            i = samples[sample]
            j = otus[otu]
            x[i,j] += int(count)

x = pd.DataFrame(x)

sort_otus = sorted(otus.keys(), key=lambda a: otus[a])
sort_samp = sorted(samples.keys(), key=lambda a: samples[a])
keep = ((x > 0).sum(axis=0) > args.min_samples)
x = x.ix[:, keep]

out = open(args.out, 'w')
out.write( 'sample\t' + '\t'.join(sort_otus) + '\n')
for sample in sort_samp:
    i = samples[sample]
    out.write( '%s\t%s' %(sample, '\t'.join(['%d' %(xi) for xi in x.ix[i,:]])) + '\n')
