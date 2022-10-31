# Calculate the distribution of minimum Hamming distance between each sequence
# and all other sequences in that set.

import random
import numpy as np
from sklearn.metrics.pairwise import pairwise_distances_chunked
import matplotlib.pyplot as plt

LENGTH = 30

def second_smallest(a):
    return np.partition(a, 1)[1]

assert second_smallest(np.array([1, 4, 3, 1])) == 1
assert second_smallest(np.array([1, 4, 3, 2])) == 2

def reduce_func(D_chunk, start):
  lowest_dists = [second_smallest(d) for d in D_chunk]
  return lowest_dists

# test_array1 = np.array([[0,0,0,0,0], [0,0,0,0,1], [0,1,0,1,1], [1,0,1,1,1]])
# D_chunk = next(pairwise_distances_chunked(test_array1, metric='hamming'))
# print(D_chunk * 5)
# reduced = pairwise_distances_chunked(test_array1, metric='hamming', reduce_func=reduce_func)
# print([j * 5 for i in reduced for j in i])

random_seqs = np.array([np.array(random.choices([0, 1, 2, 3], k=LENGTH)) for i in range(10000)])
# print(random_seqs[0])
reduced = pairwise_distances_chunked(random_seqs, metric='hamming', reduce_func=reduce_func)
min_dists = [j * LENGTH for i in reduced for j in i]
plt.hist(min_dists, bins=range(int(max(min_dists) + 1)))
plt.show()