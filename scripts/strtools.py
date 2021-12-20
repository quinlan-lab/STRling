from Bio.Seq import Seq
import numpy as np

def circular_permuted(x):
    """
    Args:
        x (iterator)
    Returns:
        list: All cicular permutations of x
    """
    return([x[i:] + x[:i] for i in range(len(x))])

def normalise_str(in_dna):
    """Find all possible eqivalent STR sequences.
    And return the first alphabetically.
    For example, TA = AT. But would return AT.
    """
    if in_dna == None or type(in_dna) != str or len(in_dna) == 0:
        return ''
    all_possible = []
    # Circularly permute original sequence and reverse complement
    for permuted_seq in circular_permuted(in_dna): # Switch to faster permutation (6)
        all_possible.append(permuted_seq)
    # Sort and take the first
    all_possible.sort()
    return(all_possible[0])

