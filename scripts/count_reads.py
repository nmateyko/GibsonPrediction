#!/usr/bin/env python

# Takes a list of sequences and counts the occurrences of each unique sequence. Similar sequences are grouped together to account
# for sequencing errors or mutations during plasmid amplification/library prep.

import gzip
import argparse
from collections import Counter
import numpy as np
import random
from umi_tools import UMIClusterer

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='inFP',	metavar='<inFile>', help='Input file', required=True)
parser.add_argument('-o', dest='outFP', metavar='<outFile>', help='Where to output results', required=True)
parser.add_argument('-d', dest='downsample', type=int, metavar='<downsample n>', help='Number of reads to downsample to')
parser.add_argument('-l', dest='logFP', metavar='<logFile>', help='Where to output errors/warnings', required=True)
parser.add_argument('-t', dest='threshold', type=int, default=2, help='Threshold for umi clustering')
parser.add_argument('-m', dest='method', type=str, default='directional', choices=['directional', 'cluster', 'adjacency'], help='Algorithm for umi clustering')
args = parser.parse_args()

clusterer = UMIClusterer(cluster_method=args.method)

with gzip.open(args.inFP, 'rt') as inFile:
  seqs = inFile.read().splitlines()
seqs = [seq.encode() for seq in seqs] # UMI_tools requires binary encoded sequences
if args.downsample: 
  seqs = random.sample(seqs, args.downsample)
counts = dict(Counter(seqs))
clustered = clusterer(counts, threshold=args.threshold)

# For each cluster of similar sequences, sum the number of occurrences
# of each sequence in the cluster to get the total count for that cluster.
# The first sequence in the cluster should be the true sequence.
cluster_counts = []
for cluster in clustered:
  sum = 0
  for seq in cluster:
    sum += counts[seq]
  cluster_counts.append((cluster[0], sum))

with open(f'{args.outFP}_counts.npy', 'wb') as f:
  np.save(f, cluster_counts)

with open(f'{args.outFP}_clusters.npy', 'wb') as f:
  np.save(f, np.array(clustered, dtype=object))
