# Takes two fastq files (read 1 and read 2) and combines reads into one sequence with strandedness of read 1.
# Assumes the fastqs have been trimmed of adapter sequences and that read 2 is aligned with read 1
# (i.e. they fully overlap), but there may be some substitution errors.
# Mostly taken from https://github.com/Carldeboer/CisRegModels/blob/master/alignFastqsIntoSeqs.py

from itertools import islice

# stolen from https://gist.github.com/jakebiesinger/759018/1b7d6bd6967780a8bbae743760c37885bdf86467
def readFastq(fastqfile):
    "parse a fastq-formatted file, yielding a (header, sequence, quality) tuple"
    fastqiter = (l.strip('\n') for l in fastqfile)  # strip trailing newlines 
    fastqiter = filter(lambda l: l, fastqiter)  # skip blank lines
    while True:
        fqlines = list(islice(fastqiter, 4))
        if len(fqlines) == 4:
            header1,seq,header2,qual = fqlines
        elif len(fqlines) == 0:
            return
        else:
            raise EOFError("Failed to parse four lines from fastq file!")

        if header1.startswith('@') and header2.startswith('+'):
            yield header1[1:], seq, qual
        else:
            raise ValueError("Invalid header lines: %s and %s for seq %s" % (header1, header2, seq))


def revcomp(seq):
    if not set(seq).issubset({'A', 'C', 'G', 'T', 'N'}):
        raise ValueError(f"Sequence {seq} must only contain ACTGN")
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join([comp[i] for i in seq[::-1]])


def testAlignment(seq1, seq2, offset):
    score = 0.0
    for i in range(offset, len(seq1)):
        if seq1[i] == seq2[i-offset]:
            score += 1
        #else:
        #	print("mismatch at %i"%(i));
    return score/(len(seq1)-offset) 


def align(seq1, seq2, offset, overlapRange):
    for i in range(0,overlapRange+1):
        if testAlignment(seq1,seq2, offset + i) > 0.75:
            return (offset + i)
        elif testAlignment(seq1,seq2, offset - i) > 0.75:
            return (offset - i)


def getConsensus(read1, read2, offset):
    consensus  = []
    #print("from %i to %i"%(offset,len(read1[1])))
    for i in range(offset, len(read1[1])):
        if read1[2][i] > read2[2][i-offset]: #converts to ASCII, corresponds to quality; if equal, doesn't matter which is better since base is the same
            consensus.append(read1[1][i])
        else:
            consensus.append(read2[1][i-offset])
    #print("read 2 from %i"%((len(read1[1])-offset)))
    return read1[1][0:(offset)] + "".join(consensus) + read2[1][(len(read1[1])-offset):]


r1_fastq_fp = snakemake.input[0]
r2_fastq_fp = snakemake.input[1]
out_fp = snakemake.output[0]
log_fp = snakemake.log[0]

with open(r1_fastq_fp, 'rt') as r1, \
     open(r2_fastq_fp, 'rt') as r2, \
     open(out_fp, 'wt') as out_file, \
     open(log_fp, 'wt') as log_file:

    r1_reader = readFastq(r1)
    r2_reader = readFastq(r2)
    
    for r1_read, r2_read in zip(r1_reader, r2_reader):
        r1_seq = r1_read[1]
        r2_seq = r2_read[1]

        r2_read = (r2_read[0], revcomp(r2_read[1]), r2_read[2][::-1])
        if testAlignment(r1_read[1], r2_read[1], 0) > 0.8:
            consensus = getConsensus(r1_read, r2_read, 0)
            out_file.write(f"{consensus}\n")
        else:
            log_file.write(f"Skipping {r1_read[0]}/{r2_read[0]} because they don't align within given parameters\nr1: {r1_seq}\nr2: {r2_seq}\n")