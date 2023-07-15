import subprocess, random, itertools, editdistance

from itertools import pairwise, combinations


import pandas as pd

## Finding matches with RNADuplex
def invoke_rnaduplex(long_sequence: str, short_sequence: str, energy_range: float = 5.0,
                     rnaduplex_location: str = "/usr/bin/RNAduplex") -> tuple:

    input_sequence = f"{long_sequence}\n{short_sequence}".encode()

    rnaduplex_subprocess = subprocess.Popen(
        [rnaduplex_location, "-e", f"{energy_range}", "-s"],
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    output, error = rnaduplex_subprocess.communicate(input=input_sequence)
    rnaduplex_subprocess.wait()

    first_line = output.decode().split("\n")[0].split()

    dot_bracket_long, dot_bracket_short = first_line[0].split("&")
    start_long, end_long = map(int, first_line[1].split(","))
    start_short, end_short = map(int, first_line[3].split(","))
    energy = float(first_line[-1].strip("()"))

    return start_long, end_long, dot_bracket_long, start_short, end_short, dot_bracket_short, energy



def find_matches_with_rnaduplex(df):

    mrna_starts = []
    mrna_ends = []
    # mrna_dot_brackets = []
    mirna_starts = []
    mirna_ends = []
    mirna_dot_brackets = []
    energies = []

    for _, row in df.iterrows():
        start_long, end_long, dot_bracket_long, start_short, end_short, dot_bracket_short, energy = invoke_rnaduplex(
            row.mrna_sequence, row.mirna_sequence)

        mrna_starts.append(start_long)
        mrna_ends.append(end_long)
        # mrna_dot_brackets.append(dot_bracket_long)
        mirna_starts.append(start_short)
        mirna_ends.append(end_short)
        mirna_dot_brackets.append(dot_bracket_short)
        energies.append(energy)

    df = pd.DataFrame({"mrna_start": mrna_starts,
                       "mrna_end": mrna_ends,
                       "pred_energy": energies,
                       "mirna_start": mirna_starts,
                       "mirna_end": mirna_ends,
                       "mirna_dot_bracket_5to3": mirna_dot_brackets,
                       })

    return df


def find_most_different_string(original_string, column, cache=None):
    if cache is None:
        cache = {}
    max_distance = float('-inf')
    most_different_string = ''
    for sequence in column:
        if sequence in cache:
            distance = cache[sequence]
        else:
            distance = editdistance.eval(original_string, sequence)
            cache[sequence] = distance
        if distance > max_distance:
            max_distance = distance
            most_different_string = sequence
    return most_different_string

####################################################################
# random small methods

def clash():
    return pd.read_csv("data/processed/clash/clash_parsed.csv")
####################################################################
# feature columns

def generate_alignment_string_from_dot_bracket(df):

    full_strings = []

    for _, row in df.iterrows():

        start_string = (row.mirna_start - 1) * "0"
        mid_string = row["mirna_dot_bracket_5to3"].replace(
            ".", "0").replace(")", "1")
        end_string = (len(row.mirna_sequence) - row.mirna_end) * "0"

        full_string = start_string + mid_string + end_string

        full_strings.append(full_string)

    df["alignment_string"] = full_strings

    return df


def generate_match_count_columns(df):

    def count_ones(str, seed=False):
        return str[1:7].count("1") if seed else str.count("1")

    df["pred_num_basepairs"] = df["alignment_string"].apply(count_ones)

    df["pred_seed_basepairs"] = df["alignment_string"].apply(
        count_ones, seed=True)

    return df




######

# deprecated functions


def smith_waterman(sequence1, sequence2, match_score, mismatch_score, gap_open_penalty, gap_extension_penalty):

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


#############
# random functions


def generate_random_seq(length, gc_content=0.5):
    """
    Generate a random DNA sequence of the given length and GC content.

    Args:
        length (int): The length of the sequence.
        gc_content (float): The desired GC content, between 0 and 1.

    Returns:
        str: The randomly generated sequence.

    Example:
        >>> generate_random_seq(10, 0.6)
        'GCGCTAGCTA'
    """
    gc_nuc_length = round(gc_content * length)
    gc_seq = random.choices("GC", k=gc_nuc_length)
    at_seq = random.choices("AT", k=length - gc_nuc_length)
    dna_seq = gc_seq + at_seq
    random.shuffle(dna_seq)
    return "".join(dna_seq)


################################
# features 

def create_midpoint_dict(df):
    """creates a dictionary with miRNA names as keys and their average MRE positions as values, only for miRNAs that have >= 2 MREs

    Args:
        df (pd.DataFrame): results from find_matches()

    Returns:
        dict: dictionary with miRNA names as keys and their average MRE positions as values
    """
    # calculate midpoint for each MRE
    df["midpoint"] = df["mrna_start"] + 4

    # group by miRNA name and collect midpoints for each group in a list
    grouped = df.groupby("mirna_accession")["midpoint"].apply(list)

    # filter out miRNAs with less than two MREs
    filtered = grouped[grouped.apply(len) >= 2]

    return filtered.to_dict()




def generate_close_proximity_column(df: pd.DataFrame, m: int = 13, n: int = 35) -> pd.DataFrame:
    """
    Generates a new column that checks for miRNA response elements (MREs) of the same miRNA
    in close proximity (i.e., with a distance between 13 and 35 nucleotides).

    Args:
        df (pd.DataFrame): The dataframe output of the 'find_matches' function.
        m (int): The minimum distance between MREs.
        n (int): The maximum distance between MREs.

    Returns:
        pd.DataFrame: The result dataframe with an additional 'close_proximity' column.
    """

    # Add a midpoint column to the dataframe.
    df["midpoint"] = df["mrna_start"] + df["alignment_string"].str.len() / 2

    # Create a dictionary of midpoints for each miRNA.
    midpoint_dict = create_midpoint_dict(df)

    # Create a matrix of distances between all pairs of midpoints.
    distances = [[abs(a - b) for (a, b) in combinations(midpoint_dict[mirna], 2)]
                 for mirna in midpoint_dict]

    # Create a dictionary of miRNA names that fit the proximity criteria.
    temp_dict = {list(midpoint_dict)[x]: distances[x][y]
                 for x, line in enumerate(distances)
                 for y, i in enumerate(line)
                 if m + 4 <= i <= n + 4}

    # Find MRE positions having the desired distance from each other.
    results = []
    for name, temp_distance in temp_dict.items():
        all_distances = midpoint_dict[name]
        results.extend(
            [name, [v1, v2]]
            for v1, v2 in pairwise(all_distances)
            if abs(v2 - v1) == temp_distance
        )

    # Create a boolean mask to identify the matching rows.
    df_filter = df.apply(lambda row: any((row["mirna_accession"] == i[0]) and (row["midpoint"] == j)
                                         for i in results for j in i[1]), axis=1)

    # initiating close proximity column
    df["close_proximity"] = 0
    
    # Set the corresponding values to 1.
    df.loc[df_filter, "close_proximity"] = 1

    return df