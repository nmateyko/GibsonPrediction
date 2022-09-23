#!/usr/bin/env python

# Takes two fastq files (read 1 and read 2) and combines reads into one sequence with strandedness of read 1.
# Assumes the fastqs have been trimmed of adapter sequences and that read 2 is aligned with read 1
# (i.e. they fully overlap), but there may be some substitution errors.
# Mostly taken from https://github.com/Carldeboer/CisRegModels/blob/master/alignFastqsIntoSeqs.py

import gzip
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i1', dest='inFP1',	metavar='<inFile>', help='Input file of fastq.gz for read1', required=True)
parser.add_argument('-i2', dest='inFP2',	metavar='<inFile>', help='Input file of fastq.gz for read2', required=True)
parser.add_argument('-o', dest='outFP', metavar='<outFile>', help='Where to output results', required=True)
parser.add_argument('-l', dest='logFP', metavar='<logFile>', help='Where to output errors/warnings', required=True)
parser.add_argument('-v', dest='verbose', action='count', help='Verbose output?', required=False, default=0)
args = parser.parse_args()

def revcomp(seq):
	if not set(seq).issubset({'A', 'C', 'G', 'T', 'N'}):
		raise ValueError(f"Sequence {seq} must only contain ACTGN")
	comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
	return "".join([comp[i] for i in seq[::-1]])

def getNextRead(inFile):
	name = inFile.readline().rstrip()
	seq = inFile.readline().rstrip()
	plus = inFile.readline().rstrip()
	quality = inFile.readline().rstrip()
	return (name,seq,quality)

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
		if read1[2][i]> read2[2][i-offset]: #converts to ASCII, corresponds to quality; if equal, doesn't matter which is better since base is the same
			consensus.append(read1[1][i])
		else:
			consensus.append(read2[1][i-offset])
	#print("read 2 from %i"%((len(read1[1])-offset)))
	return read1[1][0:(offset)] + "".join(consensus) + read2[1][(len(read1[1])-offset):]


with gzip.open(args.inFP1,'rt') as inFile1, \
	   gzip.open(args.inFP2,'rt') as inFile2, \
		 gzip.open(args.outFP,'wt') as outFile, \
		 gzip.open(args.logFP,'wt') as logFile:

	notDone = True
	i = 0
	while notDone:
		i += 1
		d1 = getNextRead(inFile1)
		d2 = getNextRead(inFile2)

		if d1[0] == "" or d2[0] == "":
			if d1[0] == "" and d2[0] == "":
				notDone = False
			else:
				raise Exception(f"One of the files ended prematurely at fastq entry {i} (line {i * 4})")
		else:
			d2 = (d2[0], revcomp(d2[1]), d2[2][::-1])
			if testAlignment(d1[1], d2[1], 0) > 0.8:				
				consensus = getConsensus(d1, d2, 0)
				outFile.write(f"{consensus}\n")
			else: 
				logFile.write(f"Skipping {d1[1]}/{d2[1]} because they don't align within given parameters\n"); 