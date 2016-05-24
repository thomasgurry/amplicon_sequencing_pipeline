"""
Script to take a fasta_filemap.txt file, relabel the fastas with their sampleID and sequence number, and concatenate into one big fasta file.
"""

import util
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fmap', help='tab delimited fasta<->sample map')
parser.add_argument('-d', '--dataset', default='myDataset', help='Dataset ID')

args = parser.parse_args()
filemap = args.fmap

with open(filemap, 'r') as f:
    lines = f.readlines()

lines = [line.strip().split('\t') for line in lines]
sids = [line[1] for line in lines]
fastas = [line[0] for line in lines]

## Relabel sequences in individual fastas. Save each one as *.sb
## Also concatenate all of these fastas into one large fasta, named myDataset.raw_concat.fasta

concat_fasta = args.dataset + '.raw_concat.fasta'
with open(concat_fasta, 'w') as concat:
    for fasta, newsid in zip(fastas, sids):
        new_fasta = fasta + '.sb'
        counter = 0
        with open(new_fasta, 'w') as f:
            print(new_fasta)
            for oldsid, seq in util.iter_fst(fasta):
                counter += 1
                f.write('>' + newsid + '_' + str(counter) + '\n')
                f.write(seq + '\n')
                
                # @Thomas: I think concatenating at the same time as relabeling will be fastest, but I'm not sure!
                concat.write('>' + newsid + '_' + str(counter) + '\n')
                concat.write(seq + '\n')

