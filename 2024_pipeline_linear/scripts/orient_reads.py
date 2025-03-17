# This script takes a pair of fastq files (read 1 and read 2) which contain reads
# in a mix of directions (due to ligation-based library prep) and orients all
# reads to be in the same direction given known start and end sequences.

import gzip
from itertools import islice
from rapidfuzz.distance import Levenshtein

COMP_TABLE = str.maketrans("ACTGN", "TGACN")
        
def revcomp(seq):
    '''Return the reverse complement of a DNA sequence'''
    return seq.translate(COMP_TABLE)[::-1]

def read_fastq(fastqfile):
    '''Parse a fastq-formatted file, yielding a (header, sequence, quality) tuple'''
    fastqiter = (l.strip('\n') for l in fastqfile)  # strip trailing newlines
    fastqiter = filter(lambda l: l, fastqiter)  # skip blank lines
    while True:
        fqlines = list(islice(fastqiter, 4))
        if len(fqlines) == 4:
            header1, seq, header2, qual = fqlines
        elif len(fqlines) == 0:
            return
        else:
            raise EOFError("Failed to parse four lines from fastq file!")

        if header1.startswith('@') and header2.startswith('+'):
            yield header1, seq, qual
        else:
            raise ValueError("Invalid header lines: %s and %s for seq %s" % (header1, header2, seq))
        
def get_orientation(r1_seq, r2_seq, up_start, down_start, threshold=3):
    '''
    Returns "fwd" if read orientation matches fwd_start sequence argument,
    "rc" if the read is in reverse complement orientation, or None if the start
    sequence doesn't match either orientation within the threshold 
    (i.e. edit distance is greater than threshold).
    '''

    up_dist_1 = Levenshtein.distance(up_start, r1_seq[:len(up_start)])
    up_dist_2 = Levenshtein.distance(up_start, r2_seq[:len(up_start)])

    if up_dist_1 <= threshold:
        return "fwd_up"
    elif up_dist_2 <= threshold:
        return "rc_up"
    else:
        down_dist_1 = Levenshtein.distance(down_start, r1_seq[:len(down_start)])
        down_dist_2 = Levenshtein.distance(down_start, r2_seq[:len(down_start)])
        if down_dist_1 <= threshold:
            return "rc_down"
        elif down_dist_2 <= threshold:
            return "fwd_down"
        else:
            return None
        
up_start = snakemake.params.fwd_start
down_start = snakemake.params.rev_start
r1_in_fp = snakemake.input[0]
r2_in_fp = snakemake.input[1]
r1_out_up_fp = snakemake.output[0]
r2_out_up_fp = snakemake.output[1]
r1_out_down_fp = snakemake.output[2]
r2_out_down_fp = snakemake.output[3]
log_fp = snakemake.log[0]


with (
    gzip.open(r1_in_fp, 'rt') as in_f1,
    gzip.open(r2_in_fp, 'rt') as in_f2,
    open(r1_out_up_fp, 'wt') as out_f1,
    open(r2_out_up_fp, 'wt') as out_f2,
    open(r1_out_down_fp, 'wt') as out_f3,
    open(r2_out_down_fp, 'wt') as out_f4,
    open(log_fp, 'wt') as log_f
):
    r1_reader = read_fastq(in_f1)
    r2_reader = read_fastq(in_f2)

    for (r1_head, r1_seq, r1_qual), (r2_head, r2_seq, r2_qual) in zip(r1_reader, r2_reader):
        if r1_head.split('/')[0] != r1_head.split('/')[0]:
            raise ValueError(f"Read 1 and read 2 headers don't match: {r1_head} and {r2_head}")
        
        orientation = get_orientation(r1_seq, r2_seq, up_start, down_start, threshold=3)

        if not orientation:
            log_f.write(f"No match: {r1_head.split('/')[0]}\n")
        
        elif orientation == "fwd_up":
            out_f1.write(f"{r1_head}\n{r1_seq}\n+\n{r1_qual}\n")
            out_f2.write(f"{r2_head}\n{r2_seq}\n+\n{r2_qual}\n")

        elif orientation == "rc_up":
            out_f1.write(f"{r1_head}\n{revcomp(r1_seq)}\n+\n{r1_qual[::-1]}\n")
            out_f2.write(f"{r2_head}\n{revcomp(r2_seq)}\n+\n{r2_qual[::-1]}\n")

        elif orientation == "fwd_down":
            out_f3.write(f"{r1_head}\n{r1_seq}\n+\n{r1_qual}\n")
            out_f4.write(f"{r2_head}\n{r2_seq}\n+\n{r2_qual}\n")

        else:
            out_f3.write(f"{r1_head}\n{revcomp(r1_seq)}\n+\n{r1_qual[::-1]}\n")
            out_f4.write(f"{r2_head}\n{revcomp(r2_seq)}\n+\n{r2_qual[::-1]}\n")
