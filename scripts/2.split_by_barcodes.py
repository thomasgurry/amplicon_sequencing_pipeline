import argparse, primer, sys, util
from string import maketrans

# Demultiplex FASTA/FASTQ file

def parse_args():
    
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', default='', help='Input FASTA file')
    parser.add_argument('-q', default='', help='Input FASTQ file')
    parser.add_argument('-b', default='', help='Barcodes file (samples -> barcodes)', required=True)
    parser.add_argument('-B', default='', help='Barcodes file format', choices=['fasta', 'tab'], required=True)
    parser.add_argument('-i', default='', help='Index file (seqs -> barcodes)')
    parser.add_argument('-I', default='', help='Index file format', choices=['fasta', 'fastq', 'tab'])
    parser.add_argument('-d', default=0, help='Max barcode differences', type=int)
    parser.add_argument('-w', default=5, help='Search positions 1-w for barcode', type=int)
    parser.add_argument('--mode', default=1, type=int, help='Barcodes in [1] seqids, [2] seqs, [3] index file', choices=[1,2,3], required=True)
    parser.add_argument('--rc', default=False, action='store_true', help='Reverse complement barcodes?')
    parser.add_argument('-o', help='Output file', required=True)
    args = parser.parse_args()
    
    # Check for consistency
    if not args.f and not args.q:
        quit('Error: must specify FASTA or FASTQ file')
    if args.i and not args.I:
        quit('Error: must specify index file format')
    return args


rctab = maketrans('ACGTacgt','TGCAtgca')
def reverse_complement(x):
    # Reverse complement a sequence
    return x[::-1].translate(rctab)


def parse_barcodes_file(map_fn, format='fasta', rc=False):
    # Map barcodes to samples
    b2s = {} # maps barcodes to samples
    # Case 1: barcodes file is FASTA format
    if format == 'fasta':
        for [s,b] in util.iter_fst(map_fn):
            if rc == True:
                seq = reverse_complement(s)
            b2s[b] = s
    # Case 2: barcodes file is tab-delimited
    elif format == 'tab':
        with open(map_fn) as fid:
            all_lines = fid.readlines()
        for line in open(map_fn):
            [s,b] = line.rstrip().split()
            if rc == True:
                b = reverse_complement(b)
            b2s[b] = s
    # Return map of barcodes to samples
    return b2s


def parse_index_file(index_fn, format='fasta'):
    # Map FASTQ sequences to their barcodes
    s2b = {} # maps sequences to barcodes
    # Case 1: index file is FASTA format
    if format=='fasta':
        for [s,b] in util.iter_fst(index_fn):
            # note: I'm pretty sure this won't work for downstream, because you need
            # to remove the first character from sequence ID
            s2b[s] = b  
    # Case 2: index file is tab-delimited
    elif format=='tab':
        for line in open(index_fn):
            [s,b] = line.rstrip().split()
            s2b[s] = b
    # Case 3: index file is FASTQ format
    elif format=='fastq':
        for [s,b,_,_] in util.iter_fsq(index_fn):
            # If sequence ID has :Y:0: thing at the end (standard Illumina format), remove it
            # For this kind of fastq line: @SL-MAJ:AY3TB170104:AY3TB:1:1101:10000:7854 :N:0:
            s = s.rsplit(' ', 1)[0]
            s2b[s[1:]] = b
    return s2b


def extract_barcode_from_id(line):
    # for this type of fasta line:
    # @MISEQ:1:1101:14187:1716#ATAGGTGG/1
    bcode = line.split('#')[-1].split('/')[0]
    return bcode


def mismatches(seq, subseq, w):
    # Calculate the number of mismatches between a sequence and a given subsequence
    # Searches in a sliding window that starts at position 1 and ends at position w
    best_i = 0 # index (start position)
    best_d = len(seq) # edit distance
    # for every start position
    for i in range(w):
        # calculate edit distance to the given subsequence
        if len(seq[i:]) < len(subseq):
            continue
        d = primer.MatchPrefix(seq[i:], subseq)
        # keep track of the best index and edit distance
        if d < best_d:
            best_i = i
            best_d = d
    return [best_i, best_d]


def find_best_match(seq, b2s, w, max_diff):
    # Find the sample with the best matching barcode
    best_i = ''
    best_b = '' # barcode
    best_d = len(seq) # edit distance
    # Calculate edit distance to every barcode
    for b in b2s:
        [i,d] = mismatches(seq, b, w) # index, edit distance
        if d < best_d:
            best_i = i
            best_b = b
            best_d = d
    # Return [index, edit distance, barcode, sample id] of best match
    if best_d <= max_diff:
        return [best_i, best_d, best_b, b2s[best_b]]
    else:
        return ['', '', '', '']


def run():
    # Maps FASTQ sequences to samples by finding the best matching barcodes
    # Create new sequence ids of the form: sample_count, for compatibility with QIIME etc.
    
    # Initialize variables
    args = parse_args()
    b2s = parse_barcodes_file(args.b, format=args.B, rc=args.rc) # barcodes to samples
    s2b = parse_index_file(args.i, format=args.I) # samples to barcodes
    s2c = {} # count number for a given sample
    out = open(args.o, 'w')
    
    # Get FASTA, FASTQ filenames and iterators
    if args.f:
        fn = args.f
        iter_fst = util.iter_fst
    if args.q:
        fn = args.q
        iter_fst = util.iter_fsq
    
    # For every record in FASTA/FASTQ file...
    for record in iter_fst(fn):
        sid = record[0][1:] # id
        seq = record[1] # sequence
        qual = ''
        if len(record) > 2:
            qual = record[3]
        
        # Case 1: barcodes are in the sample IDs
        if args.mode == 1:
            # Extract barcode from sequence id
            barcode = extract_barcode_from_id(record[0])
            # Find best matching sample
            [i,d,b,s] = find_best_match(barcode, b2s, args.w, args.d)
        
        # Case 2: barcodes are in the sequences
        elif args.mode == 2:
            # Search sequence for best barcode
            [i,d,b,s] = find_best_match(seq, b2s, args.w, args.d)
            # Trim barcode from sequence
            if i != '':
                record[1] = seq[i+len(b):]
                if qual != '':
                    record[3] = qual[i+len(b):]
        
        # Case 3: barcodes are in index file
        elif args.mode == 3:
            # Just in case this is Illumina data with the :Y:0: suffix,
            # remove it from the sequence ID. (When we made s2b from the
            # index file, we removed that suffix too, since merged reads
            # automatically have it removed by usearch.)
            sid = sid.rsplit(' ', 1)[0]
            # Get barcode from index file
            b = s2b[sid]
            # Find best matching sample
            [i,d,b,s] = find_best_match(b, b2s, args.w, args.d)
        
        # If sample found, replace seqid with new seqid
        if s:
            s2c[s] = s2c.get(s, 0) + 1 # Increment sample count
            new_sid = '@%s_%d' %(s, s2c[s])
            record[0] = new_sid
            out.write('\n'.join(record) + '\n')
    out.close()

if __name__ == "__main__":
    run()
