# this module contains helper functions for the v2.0 pipeline



from scripts.nucleotide_toolkit import *

import pandas as pd

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
        

# def import_fasta(fasta):
#     """
#     Parses a FASTA file and returns a string of transcripts with thymine
#     bases converted to uracil.

#     :param fasta: The path to the FASTA file to parse.
#     :type fasta: str
#     :return: A string of transcripts with thymine bases converted to uracil.
#     :rtype: str
#     """
#     with open(fasta) as f:
#         records = SeqIO.parse(f, "fasta")
#         transcripts = [thymine_to_uracil(str(rec.seq)) for rec in records]
#         return "".join(transcripts)

def parse_fasta_file(file_path):
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

def find_matches_slice(sequence_slice, targetscan_df, start_pos=0, allow_wobbles=False, minimum_matches=7):
    
    # Preparing stuff
    names = targetscan_df["name"].tolist()
    mirna_sequences = targetscan_df["sequence"].tolist()
    name_results, starts, alignment_strings, pair_counts, wobble_counts = [], [], [], [], []

    # For each miRNA
    for i, mirna_sequence in enumerate(mirna_sequences):

        # Create a generator for sliding windows
        generator = sliding_window(thymine_to_uracil(sequence_slice), len(mirna_sequence))

        # For each window
        for c, window in enumerate(generator, start=1):
            # Align the window with the miRNA sequence
            alignment_string, pair_count, wobble_count = align_sequences(window, mirna_sequence, allow_wobbles=allow_wobbles)
            # Check if the result is poor
            if pair_count < minimum_matches:
                continue
            # Add the results to the lists
            name_results.append(names[i])
            starts.append(start_pos + c - 1)
            alignment_strings.append(alignment_string)
            pair_counts.append(pair_count)
            wobble_counts.append(wobble_count)

    return pd.DataFrame(
        {
            "names": name_results,
            "start": starts,
            "alignment_string": alignment_strings,
            "no_of_base_pairs": pair_counts,
            "no_of_wobbles": wobble_counts,
            "no_of_bp+wobbles": [
                sum(x) for x in zip(pair_counts, wobble_counts)
            ],
        }
    )
    
    
    
def find_matches(sequence, mirna_df, allow_wobbles=False, minimum_matches=7):

    # Preparing stuff
    df_names = mirna_df["name"].tolist()
    
    df_sequences = mirna_df["sequence"].tolist()
    
    names, starts, mirna_sequences, mrna_sequences, alignment_strings, pair_counts, wobble_counts = [], [], [], [], [], [], []

    # For each miRNA
    for i, mirna_sequence in enumerate(df_sequences):
        # Create a generator for sliding windows
        generator = sliding_window(
            thymine_to_uracil(sequence), len(mirna_sequence))

        # For each window
        for c, window in enumerate(generator, start=1):
            # Align the window with the miRNA sequence
            alignment_string, pair_count, wobble_count = align_sequences(
                window, mirna_sequence, allow_wobbles=allow_wobbles)
            # Check if the result is poor
            if pair_count < minimum_matches:
                continue
            # Add the results to the lists
            names.append(df_names[i])
            starts.append(c)
            mirna_sequences.append(mirna_sequence)
            mrna_sequences.append(window)
            alignment_strings.append(alignment_string)
            pair_counts.append(pair_count)
            wobble_counts.append(wobble_count)

    df = pd.DataFrame(
        {
            "name": names,
            "start": starts,
            "mirna_sequence": mirna_sequences,
            "mrna_sequence": mrna_sequences,
            "alignment_string": alignment_strings,
            "no_of_base_pairs": pair_counts,
            "no_of_wobbles": wobble_counts,
            "no_of_bp+wobbles": [
                sum(x) for x in zip(pair_counts, wobble_counts)
            ],
        }
    )
    
    # adding seed column
    df["seed"] = df["mirna_sequence"].str[-8:-1]
    
    return df


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

