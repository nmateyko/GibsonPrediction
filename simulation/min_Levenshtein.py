import random
from tqdm import tqdm
import matplotlib.pyplot as plt
from Levenshtein import distance

LENGTH = 30
SIZE = 100000

random_seqs = ["".join(random.choices(['A', 'C', 'G', 'T'], k=LENGTH)) for i in range(SIZE)]

min_dists = []
for i, seq1 in tqdm(enumerate(random_seqs)):
  min = LENGTH
  for j, seq2 in enumerate(random_seqs):
    if i != j:
      d = distance(seq1, seq2)
      if d < min:
        min = d
  min_dists.append(min)

plt.hist(min_dists, bins=(range(max(min_dists) + 1)))
plt.show()
