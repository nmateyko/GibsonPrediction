import csv
import pickle
from collections import Counter

def starcode_clusters_to_counter(clustered_path, unclustered_path):
    with open(unclustered_path, 'r') as f:
        reads = [line.strip() for line in f]
    with open(clustered_path, 'r') as f:
        clustered_lines = [tuple(line) for line in csv.reader(f, delimiter='\t')]

    cluster_counts = {}

    for cluster in clustered_lines:
        indices = [int(i) for i in cluster[3].split(',')]
        cluster_seqs = [reads[i - 1] for i in indices] # Sequence IDs start at 1
        cluster_counts[cluster[0]] = Counter(cluster_seqs)

    return cluster_counts

fastq_path = snakemake.input[0]
cluster_path = snakemake.input[1]
cluster_counter = starcode_clusters_to_counter(cluster_path, fastq_path)

with open(snakemake.output[0], 'wb') as f:
    pickle.dump(cluster_counter, f, protocol=pickle.HIGHEST_PROTOCOL)