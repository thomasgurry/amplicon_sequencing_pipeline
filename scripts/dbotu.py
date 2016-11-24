#!/usr/bin/env python3
#
# author: scott olesen <swo@mit.edu>

from __future__ import print_function
import argparse, sys
import pandas as pd, numpy as np
import Levenshtein
from Bio import SeqIO
import scipy.stats

class OTU:
    '''
    Object for keeping track of an OTU's distribution and computing genetic distances
    '''
    def __init__(self, name, sequence, counts):
        '''
        name: str
          OTU ID
        sequence: str
          OTU's nucleotide sequence
        counts: numpy.Array
          length of sequence should be the number of samples
        '''
        # make this assertion so that lists of counts don't get concatenated
        self.name = name
        self.sequence = sequence
        self.counts = np.array(counts)

        self.abundance = sum(self.counts)

    def __eq__(self, other):
        return self.name == other.name and self.sequence == other.sequence and all(self.counts == other.counts)

    def __repr__(self):
        return "OTU(name={}, sequence={}, counts={})".format(repr(self.name), repr(self.sequence), repr(self.counts))

    def absorb(self, other):
        '''
        Add another OTU's counts to this one

        other: OTU

        returns: nothing
        '''
        self.counts += other.counts
        self.abundance += other.abundance

    def distance_to(self, other):
        '''
        Length-adjusted Levenshtein "distance" to other OTU

        other: OTU
          distance to this OTU

        returns: float
        '''
        return Levenshtein.distance(self.sequence, other.sequence) / (0.5 * (len(self.sequence) + len(other.sequence)))

    @staticmethod
    def _D_helper(x):
        '''A helper function for the _D method'''
        x = np.array(x)
        x = x[x > 0]
        return np.sum(x * np.log(x)) - (np.sum(x) * np.log(np.sum(x)))

    @classmethod
    def _D(cls, x, y):
        '''
        Statistic for the likelihood ratio test. See docs for mathematical derivation.
        '''
        x = np.array(x)
        y = np.array(y)
        return -2.0 * (cls._D_helper(x + y) - cls._D_helper(x) - cls._D_helper(y))

    @classmethod
    def _distribution_test_pval(cls, x, y):
        '''
        P-value from the likelihood ratio test comparing the distribution of the abundances
        of two taxa (x and y). See docs for explanation of the test.
        '''
        assert len(x) == len(y)
        df = len(x) - 1
        return scipy.stats.chi2.sf(cls._D(x, y), df=df)

    def distribution_pval(self, other):
        '''
        P-value from the likelihood ratio test comparing the distribution of the abundances
        of two OTU objects. See docs for explanation of the test.

        other: OTU

        returns: float
        '''
        return self._distribution_test_pval(self.counts, other.counts)


class DBCaller:
    '''
    Object for processing the sequence table and distance matrix into an OTU table.
    '''
    def __init__(self, seq_table, records, max_dist, min_fold, threshold_pval, log=None):
        '''
        seq_table: pandas.DataFrame
          Samples on the columns; sequences on the rows
        records: index of Bio.Seq
          Indexed, unaligned input sequences. This could come from BioPython's
          SeqIO.to_dict or SeqIO.index.
        max_dist: float
          genetic distance cutoff above which a sequence will not be merged into an OTU
        min_fold: float
          Multiply the sequence's abundance by this fold to get the minimum abundance
          of an OTU for merging
        threshold_pval: float
          P-value below which a sequence will not be merged into an OTU
        log: filehandle
          Log file reporting the abundance, genetic, and distribution checks.
        '''
        self.seq_table = seq_table
        self.records = records
        self.max_dist = max_dist
        self.min_fold = min_fold
        self.threshold_pval = threshold_pval
        self.log = log

        # get a list of the names of the sequences in order of their (decreasing) abundance
        self.seq_abunds = self.seq_table.sum(axis=1).sort_values(ascending=False)

        # check that all sequence IDs in the table are in the fasta
        missing_ids = [seq_id for seq_id in self.seq_abunds.index if seq_id not in self.records]
        if len(missing_ids) > 0:
            raise RuntimeError("{} sequence IDs found in the sequence table but not in the fasta: {}".format(len(missing_ids), missing_ids))

        # initialize OTU information
        self.membership = {}
        self.otus = []

    def ga_matches(self, candidate):
        '''
        OTUs that meet the genetic and abundance criteria

        candidate: OTU
          sequence to evaluate
        '''

        # find abundance matches
        min_abundance = self.min_fold * candidate.abundance
        abundance_matches = [otu for otu in self.otus if otu.abundance > min_abundance]

        if self.log is not None:
            print(candidate.name, 'abundance_check', *[otu.name for otu in abundance_matches], sep='\t', file=self.log)

        if len(abundance_matches) == 0:
            return []
        else:
            # find genetic matches (in order of decreasing genetic distance)
            matches_distances = [(otu.distance_to(candidate), otu) for otu in abundance_matches]
            matches_distances.sort(key=lambda x: (x[0], -x[1].abundance, x[1].name))
            matches = [otu for dist, otu in matches_distances if dist < self.max_dist]

            if self.log is not None:
                print(candidate.name, 'genetic_check', *[otu.name for otu in matches], sep='\t', file=self.log)

            return matches

    def _process_record(self, record_id):
        '''
        Process the next sequence: run the genetic, abundance, and distribution checks, either
        merging the sequence into an existing OTU or creating a new OTU.
        '''
        assert record_id in self.seq_table.index
        record = self.records[record_id]

        candidate = OTU(record.id, str(record.seq), self.seq_table.loc[record.id])

        if self.log is not None:
            print('seq', candidate.name, sep='\t', file=self.log)

        merged = False
        for otu in self.ga_matches(candidate):
            test_pval = candidate.distribution_pval(otu)

            if self.log is not None:
                print(candidate.name, 'distribution_check', otu.name, test_pval, sep='\t', file=self.log)

            if test_pval > self.threshold_pval:
                otu.absorb(candidate)
                self.membership[otu.name].append(candidate.name)
                merged = True
                break

        if not merged:
            # form own otu
            self.otus.append(candidate)
            self.membership[candidate.name] = [candidate.name]

    def generate_otu_table(self):
        '''
        Process all the input sequences to make an OTU table.

        returns: pandas.DataFrame
          OTU table (which can also be found at instance.otu_table)
        '''
        for record_id in self.seq_abunds.index:
            self._process_record(record_id)

        self.otus.sort(key=lambda otu: otu.abundance, reverse=True)
        self.otu_table = pd.DataFrame([otu.counts for otu in self.otus], index=[otu.name for otu in self.otus])
        self.otu_table.columns = self.seq_table.columns

        return self.otu_table

    def write_otu_table(self, output):
        '''
        Write the QIIME-style OTU table to a file.

        output: filehandle
        '''
        self.otu_table.to_csv(output, sep='\t')

    def write_membership(self, output):
        '''
        Write the QIIME-style OTU mapping information to a file.

        output: filehandle
        '''
        for otu in self.otus:
            print(otu.name, *self.membership[otu.name], sep='\t', file=output)


def read_sequence_table(fn):
    '''
    Read in a table of sequences. Expect a header and the sequence IDs in the
    first column. Samples are on the columns.

    fn: filename (or handle)

    returns: pandas.DataFrame
    '''
    df = pd.read_table(fn, index_col=0, header=0).astype(int)
    df.index = df.index.astype(str)
    return df

def call_otus(seq_table_fh, fasta_fh, output_fh, dist_crit, abund_crit, pval_crit, log=None, membership=None):
    '''
    Read in input files, call OTUs, and return output.

    seq_table_fh: filehandle
      sequence count table
    fasta_fh: filehandle or filename
      sequences fasta
    output_fh: filehandle
      place to write main output OTU table
    dist_crit, abund_crit, pval_crit: float
      threshold values for distance, abundance, and pvalue
    log, membership: filehandles
      places to write supplementary output
    '''

    # read in the sequences table
    seq_table = read_sequence_table(seq_table_fh)

    # set up the input fasta records
    records = SeqIO.index(fasta_fh, 'fasta')

    # generate the caller object
    caller = DBCaller(seq_table, records, dist_crit, abund_crit, pval_crit, log)
    caller.generate_otu_table()
    caller.write_otu_table(output_fh)

    if membership is not None:
        caller.write_membership(membership)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description='dbOTU3: Distribution-based OTU calling')
    p.add_argument('table', type=argparse.FileType('r'), help='sequence count table')
    p.add_argument('fasta', help='sequences (trimmed if unpaired, merged if paired-end; always unaligned)')

    g = p.add_argument_group(title='criteria')
    g.add_argument('--dist', '-d', type=float, default=0.1, metavar='D', help='maximum sequence dissimilarity when two sequences (default: 0.1)')
    g.add_argument('--abund', '-a', type=float, default=10.0, metavar='A', help='minimum fold difference for comparing two OTUs (0=no abundance criterion; default 10.0)')
    g.add_argument('--pval', '-p', type=float, default=0.0005, metavar='P', help='minimum p-value for merging OTUs (default: 0.0005)')

    g = p.add_argument_group(title='output options')
    g.add_argument('--output', '-o', default=sys.stdout, metavar='FILE', help='OTU table output (default: stdout)')
    g.add_argument('--membership', '-m', default=None, type=argparse.FileType('w'), metavar='FILE', help='QIIME-style OTU mapping file output')
    g.add_argument('--log', '-l', default=None, type=argparse.FileType('w'), metavar='FILE', help='log output')
    args = p.parse_args()

    assert args.dist >= 0
    assert args.abund >= 0.0
    assert args.pval >= 0.0 and args.pval <= 1.0

    call_otus(args.table, args.fasta, args.output, args.dist, args.abund, args.pval, args.log, args.membership)
