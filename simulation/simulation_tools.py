import random
import numpy as np
from collections import Counter

def mutated_sequence(seq, mut_rate):
  '''
  Returns a mutated version of the sequence with a mutation rate per base of mut_rate.
  '''
  mut_seq = []
  for base in seq:
    if random.random() <= mut_rate:
      mut_seq.append(random.choice(['A', 'C', 'G', 'T']))
    else:
      mut_seq.append(base)
  return "".join(mut_seq)

def amplify_library(unique_seqs, length, sample_size, cycles, bias_mean, bias_sigma, dup_mean, dup_sigma, mut_rate, error_rate):
  '''
  Simulates amplifying and then sequencing a library of random DNA.
  '''
  starting_seqs = ["".join(random.choices(['A', 'C', 'G', 'T'], k=length)) for i in range(unique_seqs)]
  biases = np.random.normal(bias_mean, bias_sigma, unique_seqs)
  amplified_seqs = [[seq] for seq in starting_seqs]

  # Simulate amplification by PCR or growth in bacteria.
  for i in range(cycles):
    prev_seqs = amplified_seqs
    amplified_seqs = []
    for seq_group, bias in zip(prev_seqs, biases):
      amplified_group = []
      for seq in seq_group:
        amplified_group.append(seq)
        p = np.random.normal(dup_mean, dup_sigma) + bias
        if random.random() <= p:
          amplified_group.append(mutated_sequence(seq, mut_rate))
      amplified_seqs.append(amplified_group)
  
  # Simulate sequencing of library with a sampling step and a error-introducing step.
  if sum([len(seq_group) for seq_group in amplified_seqs]) > sample_size:
    expanded = [(seq, group_index) for group_index, group in enumerate(amplified_seqs) for seq in group]
    expanded_sample = random.sample(expanded, sample_size)
    sampled_seqs = [[] for i in range(len(amplified_seqs))]
    for seq, group_index in expanded_sample:
      sampled_seqs[group_index].append(seq)
  else:
    sampled_seqs = amplified_seqs
    print("Fewer sequences than sample_size")
  sequenced = []
  for seq_group in sampled_seqs:
    sequenced_group = [mutated_sequence(seq, error_rate) for seq in seq_group]
    sequenced.append(sequenced_group)

  sequenced_counts = {seq: Counter(seq_group) for seq, seq_group in zip(starting_seqs, sequenced)}
  
  return sequenced_counts


def get_consensus(seqs):
  '''
  Gets consensus of sequences. Ties are broken by the base that was present in the
  first sequence in seqs.
  '''
  consensus = []
  for bases in zip(*seqs):
    base_counter = Counter(bases)
    consensus.append(base_counter.most_common(1)[0][0])
  return "".join(consensus)

def hamming_distance(s1, s2):
  return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def compare_counter_dicts(ref, clustered):
  center_not_found = []
  count_diff = []
  for true_seq in ref:
    if true_seq in clustered:
      cluster = clustered[true_seq]
      count = cluster.total()
      true_count = ref[true_seq].total()
      count_diff.append(count / true_count)
    else:
      center_not_found.append(true_seq)
  return count_diff, center_not_found

def find_closest_match(query, seqs, dist_fun):
  distances = [dist_fun(query, seq) for seq in seqs]
  min_index = min(range(len(distances)), key=lambda i: distances[i])
  return seqs[min_index], distances[min_index]

def main():
  '''Tests for functions in this file.'''
  q = 'AAAA'
  s = ['ATAT', 'AAAT', 'CCCC']
  closest, dist = find_closest_match(q, s, hamming_distance)
  assert closest == 'AAAT'
  assert dist == 1
  print('All tests passed.')

if __name__ == "__main__":
  main()