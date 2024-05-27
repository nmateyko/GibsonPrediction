import numpy as np
import pandas as pd
import random
from farmhash import fingerprint64


def hash_split(df, split_col, percentages, hash_func=fingerprint64, keep_bin=False):
    '''
    Splits a pandas dataframe by hashing the value in split_col and assigning the row
    to a split based on the hash value. Number and size of bins are defined by the percentages tuple passed.
    Hashing guarantees that a given row will end up in the same split as long as the value in the split_col
    remains the same, regardless of row order or other values in the row. Split sizes will not be exactly
    equal to the percentages passed.

    Parameters:
        df: a pandas dataframe to be split
        split_col: the df column used as an ID to assign each row to a split
        percentages: a tuple of the percent of rows assigned to each split
        hash_func: function used for hashing the values in split_col
        keep_bin: True to keep the bin column in the returned dfs, False to remove the column

    Returns:
        A list of splits (pandas dataframes), with percentages in the same order as the passed percentages tuple.
    '''
    if sum(percentages) != 100:
        raise ValueError("Sum of percentages must equal 100")
    
    # Calculate cumulative percentages to define bin range for each split
    cum_percentages = np.cumsum(percentages)
    cum_percentages = np.insert(cum_percentages, 0, 0)

    splits = []

    df = df.copy()

    # The lambda function will evenly distribute rows into bins from 0-99. Rows
    # are then assigned to splits using the bin value.
    df["bin"] = df[split_col].apply(lambda x: hash_func(x) % 100)
    for i in range(len(cum_percentages) - 1):
        splits.append(df[df['bin'].between(cum_percentages[i], cum_percentages[i + 1], inclusive="left")])
    
    if not keep_bin:
        for split in splits:
            split.drop('bin', axis=1, inplace=True)

    return splits


def main():

    print("Running hash_split tests")

    # Find a number that ends up in each bin from 0-99
    test_values = []
    for i in range(100):
        value = str(random.random())
        while fingerprint64(value) % 100 != i:
            value = str(random.random())
        test_values.append(value)

    # Place these values in a dataframe and split it. Rows 0-79 should
    # be in split 1, 80-89 in split 2, and 90-99 in split 3.
    test_df = pd.DataFrame(test_values, columns=['value'])

    split_1, split_2, split_3 = hash_split(test_df, 'value', (80, 10, 10), keep_bin=True)

    assert list(split_1['bin']) == list(range(80))
    assert list(split_2['bin']) == list(range(80, 90))
    assert list(split_3['bin']) == list(range(90, 100))
    print("hash_split tests passed")


if __name__ == '__main__':
    main()