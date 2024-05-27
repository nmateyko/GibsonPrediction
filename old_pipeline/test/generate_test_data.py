# This script generates four test fastq.gz files (read 1 and read 2 for input and assembled)
# to run through the GibsPred pipeline. The read structure matches the PCR-adaptered samples,
# not the ligation-adaptered samples.

# Input should have 5 A30, 8 T30, 6 G30, and 4 C30 sequences.
# Assembled should have 5 A30, 5 T30, 12 G30, and 4 C30 sequences.

import gzip
import random

with open('test_seqs_input.txt', 'r') as f:
    seqs_input = [line.strip() for line in f]

with open('test_seqs_assembled.txt', 'r') as f:
    seqs_assembled = [line.strip() for line in f]

def rev_comp(seq):
    comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return ''.join([comp[base] for base in seq[::-1]])

# shuffle seqs
random.shuffle(seqs_input)
random.shuffle(seqs_assembled)

revcomp_input = [rev_comp(seq) for seq in seqs_input]
revcomp_assembled = [rev_comp(seq) for seq in seqs_assembled]

# add r1 upstream and downstream to seqs and r2 upstream and downstream to revcomp
r1_upstream = "CGCCAGCTCTTC"
r1_downstream = "GAAGAGCCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
r2_upstream = "ATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGGCTCTTC"
r2_downstream = "GAAGAGCTGGCGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

seqs_input = [r1_upstream + seq + r1_downstream for seq in seqs_input]
revcomp_input = [r2_upstream + seq + r2_downstream for seq in revcomp_input]
seqs_assembled = [r1_upstream + seq + r1_downstream for seq in seqs_assembled]
revcomp_assembled = [r2_upstream + seq + r2_downstream for seq in revcomp_assembled]

def generate_fastq(seqs, filename, read):
    with gzip.open(filename, 'wt') as f:
        for i, seq in enumerate(seqs):
            f.write(f"@NOVASEQ1:272:HNJKLDSX3:2:1101:{i} {read}:N:0:NNNNNNNN+NNNNNNNN\n")
            f.write(f"{seq}\n")
            f.write(f"+\n")
            f.write(f"{'F' * len(seq)}\n")

generate_fastq(seqs_input, 'test_input_r1.fastq.gz', '1')
generate_fastq(revcomp_input, 'test_input_r2.fastq.gz', '2')
generate_fastq(seqs_assembled, 'test_assembled_r1.fastq.gz', '1')
generate_fastq(revcomp_assembled, 'test_assembled_r2.fastq.gz', '2')