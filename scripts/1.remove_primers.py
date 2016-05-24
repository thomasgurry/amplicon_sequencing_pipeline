import argparse, primer, sys, util

# Remove primer from the beginning of every sequence


def parse_args():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', default='', help='Input FASTA file')
    parser.add_argument('-q', default='', help='Input FASTQ file')
    parser.add_argument('-p', default='', help='Primer sequence')
    parser.add_argument('-l', default='', help='Primer list')
    parser.add_argument('-d', default=1, type=int, help='Max primer differences')
    parser.add_argument('-w', default=35, type=int, help='Size of search window (bp)')
    parser.add_argument('-o', default='', help='Output file (FASTA or FASTQ)')
    args = parser.parse_args()
    
    # Check for consistency
    if args.f and args.q:
        quit('Error: cannot specify both FASTA and FASTQ file')
    if args.p and args.l:
        quit('Error: cannot specify both primer and primer list')
    
    return args


def mismatches(seq, p, w):
    # Calculate the number of mismatches between a sequence and a primer
    # Searches in a sliding window from 0:w
    best_i = 0 # index of best match
    best_d = len(seq) # edit distance of best match
    # for every start position
    for i in range(w):
        # calculate edit distance to the given primer
        d = primer.MatchPrefix(seq[i:], p)
        # keep track of the best match
        if d < best_d:
            best_i = i
            best_d = d
    # Return the index and edit distance of the best match
    return [best_i, best_d]


def find_best_match(seq, primers, w, max_dist):
    # For a given sequence, find the best matching primer
    # If edit distance > max_dist, return empty match
    best_i = '' # index of best match
    best_p = '' # best matching primer
    best_d = len(seq) # edit distance of best match
    # Calculate edit distance to every primer
    for p in primers:
        [i,d] = mismatches(seq, p, w) # get index and edit distance
        if d < best_d:
            best_i = i
            best_p = p
            best_d = d
    # Return [index, edit distance, primer] of best match
    if best_d <= max_dist:
        return [best_i, best_d, best_p]
    else:
        return ['', '', '']


def run():
    # Remove primers from FASTA/FASTQ file
    
    # Get command line arguments
    args = parse_args()
    
    # Get primer sequences
    if args.p:
        primers = [args.p]
    elif args.l:
        primers = [line.rstrip() for line in open(args.l)]
    else:
        quit('Error: must specify primer or primer list')
    
    # Get FASTA/FASTQ iterators
    if args.f:
        fn = args.f
        iter_fst = util.iter_fst
    elif args.q:
        fn = args.q
        iter_fst = util.iter_fsq
    else:
        quit('Error: must specify FASTA or FASTQ file')
    
    # Iterate through FASTA/FASTQ file
    n_seqs = 0
    n_keep = 0
    out = open(args.o, 'w')
    for record in iter_fst(fn):
        n_seqs += 1
        seq = record[1]
        qual = ''
        if len(record) > 2:
            qual = record[3]
        [i,d,p] = find_best_match(seq, primers, args.w, args.d)
        if i != '':
            n_keep += 1
            seq = seq[i+len(p):]
            record[1] = seq
            if qual:
                qual = qual[i+len(p):]
                record[3] = qual
            out.write('\n'.join(record) + '\n')
    out.close()
    
    # Print statistics
    print 'Successfully removed primers from %d of %d total sequences %.2f' %(n_keep, n_seqs, 100.*n_keep/n_seqs)


run()
