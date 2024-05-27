import csv
import pandas as pd
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


def create_all_seqs_dict(cluster_counts):
  '''
  Parameters:
    cluster_counts: a dictionary containing cluster centers as keys and Counter objects as values.
                    Each Counter object contains the counts of each sequence that belongs to the
                    cluster center (the dictionary key)

  Returns:
    A dictionary where every sequence present in all Counter objects in cluster_counts is a key, and
    the value is a tuple of the count of that sequence and the cluster center sequence.

  Assumptions:
    Assumes the same sequence isn't found in two different Counter objects in the cluster_counts dict.
    This should not occur with starcode; i.e. each sequence is found in only one cluster.
  '''
  all_seqs_dict = {}
  for center in cluster_counts:
    for seq, count in cluster_counts[center].items():
      all_seqs_dict[seq] = (count, center)
  return all_seqs_dict


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


def main():
  print("Running read_starcode_results tests")

  input, assembled = read_starcode_results("test/test_clustered.txt", "test/test_unclustered.txt", 8)
  print("Input counts:\n", input)
  print("Assembled counts:\n", assembled)

  assert input['AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'].total() == 4
  assert input['GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'].total() == 2
  assert input['TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'].total() == 2
  assert assembled['AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'].total() == 4
  assert assembled['GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'].total() == 4
  assert assembled['TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'].total() == 2
  print("read_starcode_results tests passed\n")

  
  all_seqs_test = {'AAAA': Counter({'AAAA': 4, 'AAAT': 2, 'ATAA': 1}), 'TTTT': Counter({'TTTT': 5, 'TTAT': 3, 'TTAA': 2})}
  all_seqs_expected = {'AAAA': (4, 'AAAA'), 'AAAT': (2, 'AAAA'), 'ATAA': (1, 'AAAA'), 'TTTT': (5, 'TTTT'), 'TTAT': (3, 'TTTT'), 'TTAA': (2, 'TTTT')}
  all_seqs_wrong = {'AAAA': (4, 'AAAA'), 'AAAT': (2, 'AAAA'), 'ATAA': (1, 'TTTT'), 'TTTT': (5, 'TTTT'), 'TTAT': (3, 'TTTT'), 'TTAA': (2, 'TTTT')}

  print('Running create_all_seqs_dict tests')
  all_seqs_test_out = create_all_seqs_dict(all_seqs_test)
  assert all_seqs_test_out == all_seqs_expected, 'Output dict does not match expected'
  assert all_seqs_test_out != all_seqs_wrong, 'Output dict matches incorrect dict'
  print('create_all_seqs_dict tests passed\n')


  to_df_input = {'AAAA': Counter({'AAAA': 4, 'AAAT': 2, 'ATAA': 1}), 'TTTT': Counter({'TTTT': 5, 'TTAT': 3, 'TTAA': 2})}
  to_df_assembled_1 = {'AAAA': Counter({'AAAA': 3, 'AAAT': 2, 'ATAA': 1}), 'TTTT': Counter({'TTTT': 5, 'TTAT': 4, 'TTAA': 2})}
  to_df_assembled_2 = {'AAAA': Counter({'AAAA': 4, 'AAAT': 2, 'ATAA': 1})}

  to_df_expected_1 = pd.DataFrame([('AAAA', 7, 6), ('TTTT', 10, 11)], columns=['sequence', 'input_count', 'assembled_count'])
  to_df_expected_2 = pd.DataFrame([('AAAA', 7, 7), ('TTTT', 10, 0)], columns=['sequence', 'input_count', 'assembled_count'])

  print('Running cluster_counter_to_count_df tests')
  # row order may change so can't directly compare dfs; just look at them instead
  df1 = cluster_counter_to_count_df(to_df_input, to_df_assembled_1)
  df2 = cluster_counter_to_count_df(to_df_input, to_df_assembled_2)
  print("First df should have these rows (order may be different):")
  print(f'{to_df_expected_1}\n')
  print("First df:")
  print(f'{df1}\n')
  print("Second df should have these rows (order may be different):")
  print(f'{to_df_expected_2}\n')
  print("Second df:")
  print(f'{df2}\n')
  print('cluster_counter_to_count_df tests finished\n')


if __name__ == '__main__':
  main()