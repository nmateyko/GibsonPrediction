import pandas as pd
from itertools import islice


def get_fold_change(count_df, input_pseudocount, assembled_pseudocount):
    '''
    Calculates fold change for each sequence (cluster center) between input and assembled (assembled/input).
    Normalizes for number of reads in each sample before computing fold change.
    Parameters:
        count_df: a pandas dataframe containing three columns: sequence, input_count, and assembled_count
        input_pseudocount: added to the input count before calculating frequency (normalized count) and fold change;
        must be > 0 or division will fail
        assembled_pseudocount: added to the assembled count before calculating frequency (normalized count) and fold change

    Returns:
        A pandas Series containing the fold change values
    '''
    input_total = count_df['input_count'].sum()
    assembled_total = count_df['assembled_count'].sum()
    input_frequency = (count_df['input_count'] + input_pseudocount) / input_total
    assembled_frequency = (count_df['assembled_count'] + assembled_pseudocount) / assembled_total
    fold_change = assembled_frequency / input_frequency
    return fold_change


def reverse_complement(seq):
    '''
    Returns the reverse complement of a DNA sequence.
    
        Parameters:
        seq (str): A DNA sequence.
        
        Returns:
        result (str): The reverse complement of seq.
    '''
    complement_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    result = []
    for i in range(len(seq)-1, -1, -1):
        result.append(complement_dict[seq[i]])
    return "".join(result)


def contains_RE_site(seq, site):
    '''
    Returns whether a DNA sequence contains a restriction site in either strand.
    
        Parameters:
        seq (str): A DNA sequence.
        
        Returns:
        (bool): True if seq contains the restriction site, False otherwise.
    '''
    if seq.find(site) == -1 and reverse_complement(seq).find(site) == -1:
        return False
    else:
        return True


# stolen from https://gist.github.com/jakebiesinger/759018/1b7d6bd6967780a8bbae743760c37885bdf86467
def read_fastq(fastqfile):
    '''parse a fastq-formatted file, yielding a (header, sequence, quality) tuple'''
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
  

def main():
    print("Testing contains_RE_site")
    site = "AATT"
    pos = "CGAATTCGA"
    pos_rc = "TCGAATTCG"
    neg = "AATAGCTA"

    assert contains_RE_site(pos, site)
    assert contains_RE_site(pos_rc, site)
    assert not contains_RE_site(neg, site)
    print("contains_RE_site tests passed\n")

    
    print('Running get_fold_change tests')
    fc_test1_df = pd.DataFrame([('AAAA', 10, 7), ('TTTT', 8, 30)], columns=['sequence', 'input_count', 'assembled_count'])
    fc_test1_expected = pd.Series([(8/37) / (11/18), (31/37) / (9/18)])
    assert get_fold_change(fc_test1_df, 1, 1).equals(fc_test1_expected), "Fold change does not match expected"
    print('get_fold_change tests passed')


if __name__ == '__main__':
    main()
