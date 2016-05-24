# Dereplicate sequences in fasta file

import util, argparse, re, os, seqdb

def parse_args():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='Input FASTA file', default='')
    parser.add_argument('-q', help='Input FASTQ file', default='')
    parser.add_argument('-s', help='Sample ID separator', required=True)
    parser.add_argument('-o', help='Output mapping file', required=True)
    parser.add_argument('-d', help='OTU database', default='', required=True)
    parser.add_argument('-P', help='Processing summary', default='', required=True)
    parser.add_argument('-M', help='Min count', default=10, type=int)
    parser.add_argument('-S', help='Min samples', default=1, type=int)
    parser.add_argument('-l', help='Trim length', type=int, default=0)
    args = parser.parse_args()
    return args


def dereplicate(fst='', fsq='', sep='', trim_len=''):
    # Dereplicate sequences
    # NOTE:
    #      Separator for barcodes must be specified in summary file. 
    #      e.g. 'SRR230982_142' the separator is '_'
    #      
    x = {}
    if fst:
        fn = fst
        iter_fst = util.iter_fst
    if fsq:
        fn = fsq
        iter_fst = util.iter_fsq

    for record in iter_fst(fn):
        [sid, seq] = record[:2]
        sid = sid[1:]
        #sa = re.search('(.*?)%s' %(sep), sid).group(1)
#        sa = sid.split(sep)[0]
        sa = sid.split(sep)
        sa = sep.join(sa[:len(sa)-1])
        if trim_len:
            if len(seq) >= trim_len:
                seq = seq[:trim_len]
            else:
                continue
        if seq not in x:
            x[seq] = {}
#        if sid not in x[seq]:
#            x[seq][sid] = 0
#        x[seq][sid] += 1
        if sa not in x[seq]:
            x[seq][sa] = 0
        x[seq][sa] += 1

    return x

def write_output(x, map_fn, db_fn, min_size=1, min_samples=1, processing_summary_file=None):
    # Write output (database + mapping file)
    # Load SeqDB
    min_samples = 1
    db = seqdb.SeqDB(fn=db_fn)
    out = open(map_fn, 'w')
    nthrown_out = 0
    for seq in x:
        size = sum(x[seq].values())
        if size < min_size:
            nthrown_out += size
            continue
        if len(x[seq]) < min_samples:
            continue
        db.add_seq(seq, size=size)
        out.write('%s\t%s\n' %(db.db[:seq], ' '.join(['%s:%d' %(sa, x[seq][sa]) for sa in x[seq]])))
    if processing_summary_file != None:
        proc_sum = open(processing_summary_file, 'w')
        proc_sum.write('Total number of reads thrown out during dereplication (falling under minimum total read count across samples of ' + str(min_size) + '):\t' + str(nthrown_out) + '\n\n')
        proc_sum.close()
    out.close()
    db.write(db_fn)

args = parse_args()
if args.l == 0:
    args.l = ''
x = dereplicate(fst=args.f, fsq=args.q, sep=args.s, trim_len=args.l)
write_output(x, map_fn=args.o, db_fn=args.d, min_size=args.M, min_samples=args.S, processing_summary_file=args.P)
