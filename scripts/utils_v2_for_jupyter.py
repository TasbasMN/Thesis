# this module contains helper functions for the v2.0 pipeline

import itertools

from scripts.nucleotide_toolkit import *

import pandas as pd
import numpy as np

#######################################################################################################################
### small utilities

def sliding_window(sequence, win_size):
    """
    A generator that yields subsequences of length `win_size` from the input
    `sequence`.

    Args:
        sequence (list): A list of elements.
        win_size (int): The length of the subsequences to yield.

    Yields:
        list: A subsequence of length `win_size` from the input `sequence`.
    """
    for i in range(len(sequence) - win_size + 1):
        yield sequence[i: i + win_size]


def parse_fasta_file_into_rna(file_path):
    """
    Parses a FASTA file and returns the DNA sequence as a string with thymine replaced by uracil.

    :param file_path: The path to the FASTA file.
    :return: The DNA sequence as a string with thymine replaced by uracil.
    """
    sequence = ""
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                continue
            sequence += line
    return thymine_to_uracil(sequence)


#######################################################################################################################
### find matches v2.0
def align_sequences(seq1, seq2, allow_wobbles=False):
    """
    Aligns two sequences of the same size and creates a string of 1's and 0's representing nucleotide matches.
    If allow_wobbles is True, it will also check for G:U wobbles and appends 2.

    Args:
        seq1 (str): sequence 1
        seq2 (str): sequence 2
        allow_wobbles (bool, optional): if allowed, check for G:U wobbles. Defaults to False.

    Returns:
        Tuple: alignment_string containing the alignment string, pair_count containing the number of nucleotide matches, 
        and wobble_count containing the number of G:U wobbles (if True)
    """

    alignment_string = []
    pair_count = 0
    wobble_count = 0

    for nucleotide_1, nucleotide_2 in zip(seq1, seq2):
        if nucleotide_1 == nucleotide_2:
            alignment_string.append("1")
            pair_count += 1
        elif allow_wobbles and ((nucleotide_1 == "G" and nucleotide_2 == "U") or (nucleotide_1 == "U" and nucleotide_2 == "G")):
            alignment_string.append("2")
            wobble_count += 1
        else:
            alignment_string.append("0")

    return "".join(alignment_string), pair_count, wobble_count

def find_k_consecutive_bps(df, k=8):
    """Finds kmers located anywhere in the mRNA - miRNA complex.

    Args:
        df (pd.DataFrame): DataFrame containing results.
        k (int, optional): k value. Defaults to 8.

    Returns:
        pd.DataFrame: DataFrame with additional column.
    """

    kmer_colname = f"{k}_consecutive_bps"
    kmer_regex = "1{" + str(k) + "}"

    kmer_mask = df["alignment_string"].str.contains(kmer_regex)

    df[kmer_colname] = kmer_mask.astype(int)

    return df

###############################################################################################################################################

def slide_and_compare_sequences(mirna_seq, mrna_seq, start_pos, min_matches):
    """
    Slide a window of length len(mirna_seq) along mrna_seq, and check for matches
    between the two sequences. Return the starting and ending positions of windows
    with at least min_matches matches, along with the matching window
    sequences and alignment strings.
    """

    # Reverse complement miRNA sequence
    mirna_seq = reverse_complement(mirna_seq)

    # Convert sequences to arrays
    mirna_arr = np.array(list(mirna_seq))
    mrna_arr = np.array(list(mrna_seq))

    # Create sliding window view of mRNA sequence
    sliding_window = np.lib.stride_tricks.sliding_window_view(mrna_arr, len(mirna_arr))

    # Get matching window sequences
    mrna_windows = np.apply_along_axis("".join, 1, sliding_window)

    # Get number of matches for each window
    no_of_matches = np.sum(sliding_window == mirna_arr, axis=1)

    # Get starting and ending positions of matching windows
    starts = start_pos + np.arange(len(mrna_arr) - len(mirna_arr) + 1) 
    ends = starts + len(mirna_arr) - 1

    # Generate alignment strings for each window
    alignments = ["".join(["1" if i == 1 else "0" for i in row]) for row in (sliding_window == mirna_arr)]

    # Filter out windows with too few no_of_matches
    mask = no_of_matches < min_matches
    
    starts = np.delete(starts, np.where(mask))
    ends = np.delete(ends, np.where(mask))
    alignments = np.delete(alignments, np.where(mask))
    no_of_matches = np.delete(no_of_matches, np.where(mask))
    mrna_windows = np.delete(mrna_windows, np.where(mask))
    
    # calculating matches in seed regions
    match_counts_in_seed_regions = [i[-7:-1].count("1") for i in alignments]

    return starts, ends, alignments, no_of_matches, mrna_windows, match_counts_in_seed_regions


def find_matches(df):
    """
    Find all matches between miRNA and mRNA sequences in a given DataFrame.
    """
    
    ids = []
    
    pred_starts = []
    true_starts = []
    
    pred_ends = []
    true_ends = []
    
    pred_no_of_bp = []
    true_no_of_bp = []
    
    pred_no_of_bp_in_seed = []
    true_no_of_bp_in_seed = []
    
    alignment_strings = []
    mrna_sequences = []
    mirna_sequences = []
    
    true_seed_types = []
    folding_classes = []
    
    mirna_names = []
    ensgs = []
    
    accessions = []

    

    for _, row in df.iterrows():
        start, end, alignment_string, no_of_bp, mrna_sequence, pred_matches_in_seed = slide_and_compare_sequences(row.mirna_sequence, row.mrna_sequence, row.true_start, 7)
        
        # appending results of slide_and_compare
        pred_starts.extend(start.tolist())
        true_starts.extend([row.true_start] * len(start))
        
        pred_ends.extend(end.tolist())
        true_ends.extend([row.true_end] * len(start))
        
        pred_no_of_bp.extend(no_of_bp.tolist())
        true_no_of_bp.extend([row.num_basepairs] * len(start))
        
        pred_no_of_bp_in_seed.extend(pred_matches_in_seed)
        true_no_of_bp_in_seed.extend([row.seed_basepairs] * len(start))    
        
        alignment_strings.extend(alignment_string)
        mrna_sequences.extend(mrna_sequence)
        
        # datas that are exploded
        ids.extend([row.id] * len(start))        
        mirna_sequences.extend([row.mirna_sequence] * len(start))
        true_seed_types.extend([row.true_seed_type] * len(start))
        folding_classes.extend([row.true_folding_class] * len(start))
        
        mirna_names.extend([row.mirna_name] * len(start))
        ensgs.extend([row.ensg] * len(start))
        
        accessions.extend([row.mirna_accession] * len(start))
        
        

    return pd.DataFrame(
        {
            "id": ids,
            "pred_start": pred_starts,
            "true_start": true_starts,
            
            "pred_end": pred_ends,
            "true_end": true_ends,
            
            "pred_no_of_bp": pred_no_of_bp,
            "true_no_of_bp": true_no_of_bp,
            
            "pred_no_of_bp_in_seed": pred_no_of_bp_in_seed,
            "true_no_of_bp_in_seed": true_no_of_bp_in_seed,
            
            "alignment_string": alignment_strings,
            
            "mirna_sequence": mirna_sequences,
            "mrna_sequence": mrna_sequences,
            
            "true_seed_type": true_seed_types,
            "true_folding_class": folding_classes,
            
            "mirna_name": mirna_names,
            "ensg": ensgs,
            "accession": accessions
        }
    )


def ohe_true_match_types(df):
    
        
    roman_dict = {
        "I": 1,
        "II": 2,
        "III": 3,
        "IV": 4,
        "V": 5,}
    
    # replacing romans into numerical
    df["true_folding_class"] = df["true_folding_class"].replace(roman_dict)

    # creating ohe labels
    labels_df = pd.get_dummies(df["true_folding_class"])
    
    # concatenating labels + original df
    df = pd.concat([df, labels_df], axis=1)
    
    # dropping original column
    df.drop(columns="true_folding_class", inplace=True)
    
    rename_dict = {1: "true_type_1",
                2: "true_type_2",
                3: "true_type_3",
                4: "true_type_4",
                5: "true_type_5"}

    df.rename(columns=rename_dict, inplace=True)
    
    return df


def smith_waterman(sequence1, sequence2, match_score, mismatch_score, gap_open_penalty, gap_extension_penalty):
    # sourcery skip: low-code-quality
    # Initialize the scoring matrix and traceback matrix
    rows = len(sequence1) + 1
    cols = len(sequence2) + 1
    score_matrix = [[0] * cols for _ in range(rows)]
    traceback_matrix = [[0] * cols for _ in range(rows)]

    # Fill the scoring matrix
    for i, j in itertools.product(range(1, rows), range(1, cols)):
        match = score_matrix[i - 1][j - 1] + (match_score if sequence1[i - 1] == sequence2[j - 1] else mismatch_score)
        delete = score_matrix[i - 1][j] - (gap_open_penalty if traceback_matrix[i - 1][j] != 'u' else gap_extension_penalty)
        insert = score_matrix[i][j - 1] - (gap_open_penalty if traceback_matrix[i][j - 1] != 'l' else gap_extension_penalty)
        score_matrix[i][j] = max(0, match, delete, insert)

        # Update the traceback matrix
        if score_matrix[i][j] == match:
            traceback_matrix[i][j] = 'd'  # Diagonal
        elif score_matrix[i][j] == delete:
            traceback_matrix[i][j] = 'u'  # Up
        elif score_matrix[i][j] == insert:
            traceback_matrix[i][j] = 'l'  # Left

    # Find the cell with the highest score
    max_score = 0
    max_i = 0
    max_j = 0
    for i, j in itertools.product(range(rows), range(cols)):
        if score_matrix[i][j] > max_score:
            max_score = score_matrix[i][j]
            max_i = i
            max_j = j

    # Traceback from the highest-scoring cell to build the alignment
    alignment1 = ""
    alignment2 = ""
    alignment_string = ""
    i, j = max_i, max_j
    while score_matrix[i][j] != 0:
        if traceback_matrix[i][j] == 'd':  # Diagonal
            alignment1 = sequence1[i - 1] + alignment1
            alignment2 = sequence2[j - 1] + alignment2
            alignment_string = ("1" if sequence1[i - 1] == sequence2[j - 1] else "0") + alignment_string
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 'u':  # Up
            alignment1 = sequence1[i - 1] + alignment1
            alignment2 = '-' + alignment2
            alignment_string = "0" + alignment_string
            i -= 1
        elif traceback_matrix[i][j] == 'l':  # Left
            alignment1 = '-' + alignment1
            alignment2 = sequence2[j - 1] + alignment2
            alignment_string = "0" + alignment_string
            j -= 1

    # Calculate starting and ending indices of the alignment
    start_index1 = i + 1
    end_index1 = i + len(alignment1)
    start_index2 = j + 1
    end_index2 = j + len(alignment2)

    # Print the alignment with start and end indices
    print(f"{alignment1}")
    print(alignment_string)
    print(f"{alignment2}")

    # Return the start and end indices
    
    return alignment_string
    # return start_index1, end_index1, start_index2, end_index2
    
    
    
    
def clash():
    return pd.read_csv("data/processed/clash/clash_parsed.csv")