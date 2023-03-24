# This script generates two test fastq.gz files (read 1 and read 2) to run through the PLAB pipeline

import gzip
import random

with open('test_seqs.txt', 'r') as f:
    seqs = [line.strip() for line in f]

def rev_comp(seq):
    comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return ''.join([comp[base] for base in seq[::-1]])

# shuffle seqs
random.shuffle(seqs)

seqs_revcomp = [rev_comp(seq) for seq in seqs]

# add r1 upstream and downstream to seqs and r2 upstream and downstream to seqs_revcomp
r1_upstream = "TGCATTTTTTTCACATC"
r1_downstream = "GGTTACGGCTGTT"
r2_upstream = "AACAGCCGTAACC"
r2_downstream = "GATGTGAAAAAAATGCA"

seqs = [r1_upstream + seq + r1_downstream for seq in seqs]
seqs_revcomp = [r2_upstream + seq + r2_downstream for seq in seqs_revcomp]

def generate_fastq(seqs, filename, read):
    with gzip.open(filename, 'wt') as f:
        for i, seq in enumerate(seqs):
            f.write(f"@NOVASEQ1:272:HNJKLDSX3:2:1101:{i} {read}:N:0:NNNNNNNN+NNNNNNNN\n")
            f.write(f"{seq}\n")
            f.write(f"+\n")
            f.write(f"{'F' * len(seq)}\n")

generate_fastq(seqs, 'test_r1.fastq.gz', '1')
generate_fastq(seqs_revcomp, 'test_r2.fastq.gz', '2')