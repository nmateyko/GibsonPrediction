import csv
import pandas as pd
import pickle
import subprocess
from collections import Counter


def read_starcode_results(clustered_path, unclustered_path, input_count):
  '''
  Reads starcode output file and generates a dictionary of Counter objects for both
  the input and assembled libraries. Input and assembled reads are separated by their seq-id
  in the starcode output, which is the line number from the starcode input file.

  Parameters:
    clustered_path: path to the starcode output file
    unclustered_path: path to the starcode input file; this is the input and assembled reads
      concatenated into a single file, with one sequence per line
    input_count: the number of sequences in the input sample, so that input and assembled reads
      can be separated

  Returns:
    cluster_counts_input and cluster_counts_assembled: both are dictionaries with cluster centers as keys
    and Counter objects as values. Each Counter object contains the counts of each sequence that belongs to the
    cluster center (the dictionary key) in either the input sample or assembled sample.

  Assumptions:
    Both files should be uncompressed. The starcode output must be generated with both the
    --print-clusters and --seq-id flags. 
  '''
  with open(unclustered_path, 'r') as f:
    unclustered_seqs = [line.rstrip() for line in f]
  with open(clustered_path, 'r') as f:
    clustered_data = [tuple(line) for line in csv.reader(f, delimiter='\t')]

  cluster_counts_input = {}
  cluster_counts_assembled = {}

  for cluster in clustered_data:
    indices = [int(i) for i in cluster[3].split(',')]
    input_indices = [i for i in indices if i <= input_count]
    assembled_indices = [i for i in indices if i > input_count]
    input_seqs = [unclustered_seqs[i - 1] for i in input_indices] # Sequence IDs start at 1
    assembled_seqs = [unclustered_seqs[i - 1] for i in assembled_indices] # Sequence IDs start at 1
    cluster_counts_input[cluster[0]] = Counter(input_seqs)
    cluster_counts_assembled[cluster[0]] = Counter(assembled_seqs)

  return cluster_counts_input, cluster_counts_assembled


def cluster_counter_to_count_df(input_clusters, assembled_clusters):
  '''
  Converts two dictionaries of Counter objects (input clusters and assembled clusters) to
  a dataframe with sequence, input count, and assembled count columns.

  Parameters:
    input_clusters: a dictionary containing cluster centers as keys and Counter objects as values.
                    Each Counter object contains the counts of each sequence that belongs to the
                    cluster center (the dictionary key) in the input sample.
    input_clusters: a dictionary containing cluster centers as keys and Counter objects as values.
                    Each Counter object contains the counts of each sequence that belongs to the
                    cluster center (the dictionary key) in the assembled sample.

  Returns:
    A pandas dataframe with three columns: sequence, input_count, and assembled_count.
  '''
  all_centers = input_clusters.keys() | assembled_clusters.keys() # extract all keys
  all_counts = []

  for seq in all_centers:
    if seq in input_clusters:
      input_count = input_clusters[seq].total()
    else:
      input_count = 0

    if seq in assembled_clusters:
      assembled_count = assembled_clusters[seq].total()
    else:
      assembled_count = 0

    all_counts.append((seq, input_count, assembled_count))

  return pd.DataFrame(all_counts, columns=['sequence', 'input_count', 'assembled_count'])


unclustered_seqs_path = snakemake.input[0]
starcode_output_path = snakemake.input[1]
input_seqs_path = snakemake.input[2]

input_count = int(subprocess.check_output(['wc', '-l', input_seqs_path]).split()[0])
cluster_counts_input, cluster_counts_assembled = read_starcode_results(starcode_output_path, unclustered_seqs_path, input_count)
counts_df = cluster_counter_to_count_df(cluster_counts_input, cluster_counts_assembled)

with open(snakemake.output[0], 'wb') as f:
    pickle.dump(cluster_counts_input, f, protocol=pickle.HIGHEST_PROTOCOL)
with open(snakemake.output[1], 'wb') as f:
    pickle.dump(cluster_counts_assembled, f, protocol=pickle.HIGHEST_PROTOCOL)
with open(snakemake.output[2], 'wb') as f:
    pickle.dump(counts_df, f, protocol=pickle.HIGHEST_PROTOCOL)