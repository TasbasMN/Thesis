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

# -1's here convert biological coordinates into 0-index coordinates
    return start_long-1, end_long-1, dot_bracket_long, start_short-1, end_short-1, dot_bracket_short, energy



def find_matches_with_rnaduplex(df):

    mrna_starts = []
    mrna_ends = []
    mirna_starts = []
    mirna_ends = []
    mirna_dot_brackets = []
    energies = []

    for _, row in df.iterrows():
        start_long, end_long, dot_bracket_long, start_short, end_short, dot_bracket_short, energy = invoke_rnaduplex(
            row.mrna_sequence, row.mirna_sequence)

        mrna_starts.append(start_long)
        mrna_ends.append(end_long)
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


def reverse_complement_rna_to_dna(rna_sequence):
    complement = {'A': 'T', 'U': 'A', 'C': 'G', 'G': 'C'}
    reverse_seq = rna_sequence[::-1]
    return ''.join(complement[base] for base in reverse_seq)

def reverse_complement_dna_to_rna(rna_sequence):
    complement = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_seq = rna_sequence[::-1]
    return ''.join(complement[base] for base in reverse_seq)


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



################
# 7_process_vcf.ipynb

def get_nucleotides_in_interval(chrom, start, end):
    file_path = f"fasta/grch37/Homo_sapiens.GRCh37.dna.chromosome.{chrom}.fa"
    with open(file_path, 'r') as file:
        file.readline()
        byte_position = file.tell()
        line_length = len(file.readline().strip())
        start_offset = start - 1
        end_offset = end - 1
        num_start_new_lines = start_offset // line_length
        num_end_new_lines = end_offset // line_length
        start_byte_position = byte_position + start_offset + num_start_new_lines
        end_byte_position = byte_position + end_offset + num_end_new_lines
        file.seek(start_byte_position)

        # Read the nucleotides in the interval
        nucleotides = file.read(end_byte_position - start_byte_position + 1)

    # Remove newlines from the nucleotides
    nucleotides = nucleotides.replace('\n', '')

    return nucleotides


def get_nucleotide_at_position(chrom, position):
    file_path = f"fasta/grch37/Homo_sapiens.GRCh37.dna.chromosome.{chrom}.fa"
    with open(file_path, 'r') as file:
        file.readline()
        byte_position = file.tell()
        line_length = len(file.readline().strip())
        offset = position - 1
        num_new_lines = offset // line_length
        byte_position = byte_position + offset + num_new_lines
        file.seek(byte_position)

        # Read the nucleotide at the position
        nucleotide = file.read(1)
    return nucleotide

def invoke_rnaduplex_for_vcfs(long_sequence: str, short_sequence: str, energy_range: float = 5.0,
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

# -1's here convert biological coordinates into 0-index coordinates
    return start_long-1, end_long, dot_bracket_long, start_short-1, end_short, dot_bracket_short, energy

def find_matches_for_vcfs(df, mutated=False):
    
    mirna_df = pd.read_csv("data/processed/mirbase/mirbase22.csv")
    mirna_df["sequence"] = mirna_df["sequence"].apply(reverse_complement_rna_to_dna)
    
    mrna_starts = []
    mrna_ends = []
    mirna_starts = []
    mirna_ends = []
    mrna_dot_brackets =[]
    mirna_dot_brackets = []
    energies = []
    identifiers = []
    
    mrna_sequences = []
    mirna_sequences = []
    
    mirna_accessions = []

    for _, row in df.iterrows():
        for _, mirna_row in mirna_df.iterrows():
            mirna_sequence = mirna_row['sequence']
            mrna_sequence = row['mutated_sequence'] if mutated else row['sequence']
            start_long, end_long, dot_bracket_long, start_short, end_short, dot_bracket_short, energy = invoke_rnaduplex_for_vcfs(
                mrna_sequence, mirna_sequence)

            mrna_starts.append(start_long)
            mrna_ends.append(end_long)
            mirna_starts.append(start_short)
            mirna_ends.append(end_short)
            mrna_dot_brackets.append(dot_bracket_long)
            mirna_dot_brackets.append(dot_bracket_short)
            energies.append(energy)
            
            mrna_sequences.append(mrna_sequence)
            mirna_sequences.append(mirna_sequence)
            
            mirna_accessions.append(mirna_row.accession)
            
            identifiers.append(f"{row.id}_{mirna_row.accession}") 

            print(f"Processing row {row.name} in df with mirna_row {mirna_row.accession} in mirna_df")

    df = pd.DataFrame({
        "id": identifiers,
        "mrna_start": mrna_starts,
        "mrna_end": mrna_ends,
        "mrna_sequence": mrna_sequences,
        "mirna_accession": mirna_accessions,
        "mirna_start": mirna_starts,
        "mirna_end": mirna_ends,
        "mirna_sequence": mirna_sequences,
        "mrna_dot_bracket_5to3": mrna_dot_brackets,
        "mirna_dot_bracket_5to3": mirna_dot_brackets,
        "pred_energy": energies
        
    })

    return df



def augment_vcf(vcf_df, sequence_offset=30):
    """
    Augments a VCF DataFrame by adding additional columns.

    Args:
        vcf_df (pandas.DataFrame): The input VCF DataFrame.
        sequence_offset (int, optional): The offset used to fetch nucleotides in the sequence. Defaults to 30.

    Returns:
        pandas.DataFrame: The augmented VCF DataFrame.
    """

    # Add an 'id' column by combining 'chr', 'pos', 'ref', and 'alt' columns
    vcf_df.loc[:, "id"] = vcf_df.apply(lambda row: f"{row['chr']}_{row['pos']}_{row['ref']}_{row['alt']}", axis=1)

    # Add a 'fetched_nucleotides' column by fetching nucleotides at the specified position
    vcf_df.loc[:, 'fetched_nucleotides'] = vcf_df.apply(lambda x: get_nucleotide_at_position(x['chr'], x['pos']), axis=1)

    # Add an 'is_nucleotides_same' column to check if fetched nucleotides are the same as 'ref'
    vcf_df.loc[:, 'is_nucleotides_same'] = vcf_df["fetched_nucleotides"] == vcf_df["ref"]

    # Add a 'sequence' column by fetching nucleotides in the interval [pos-sequence_offset, pos+sequence_offset]
    vcf_df.loc[:, 'sequence'] = vcf_df.apply(lambda x: get_nucleotides_in_interval(x['chr'], x['pos']-sequence_offset, x["pos"]+sequence_offset), axis=1)

    # Add a 'mutated_sequence' column by replacing the nucleotide at the sequence_offset position with 'alt'
    vcf_df.loc[:,'mutated_sequence'] = vcf_df.apply(lambda row: row['sequence'][:sequence_offset] + row['alt'] + row['sequence'][sequence_offset+1:], axis=1)

    return vcf_df


################
# 8_adding_feature_cols

def generate_alignment_string_from_dot_bracket(df):
    full_strings = []
    for _, row in df.iterrows():
        start_string = (row.mirna_start) * "0"
        mid_string = row["mirna_dot_bracket_5to3"].replace(".", "0").replace(")", "1")
        end_string = (len(row.mirna_sequence) - row.mirna_end -1) * "0"
        
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

def add_ta_sps_columns(df):
    # Generate temporary seed column
    df["seed"] = df["mirna_sequence"].str[1:8].str.replace("T", "U")
    # Read ta sps data
    ta_sps_df = pd.read_csv("data/processed/ta_sps/ta_sps.csv", usecols=["seed_8mer", "ta_log10", "sps_mean"])
    ta_sps_df = ta_sps_df.rename(columns={"seed_8mer": "seed"})
    # Merge dataframes on seed column
    df = df.merge(ta_sps_df, on="seed", how="left")
    # Drop temporary column
    df.drop(columns=["seed"], inplace=True)

    return df

def add_mirna_conservation_column(df):
    targetscan = pd.read_csv("data/processed/targetscan/targetscan.csv")
    targetscan = targetscan.rename(columns={"accession": "mirna_accession", "conservation": "mirna_conservation"})
    targetscan = targetscan[["mirna_accession", "mirna_conservation"]]
    df = df.merge(targetscan, on="mirna_accession", how="left")
    return df

def find_seed_type(df):
    
    def check_nth_character(row):
        sequence = row["mrna_sequence"]
        mirna_end = row["mrna_end"] - 1
        return "x" if mirna_end >= len(sequence) else int(sequence[mirna_end] == "A")

    df["anchor_a"] = df.apply(check_nth_character, axis=1)

    df["6mer_seed"] = (df["alignment_string"].str[1:7].str.count("0") == 0).astype(int)
    df["match_8"] = (df["alignment_string"].str[7] == "1").astype(int)
    df["6mer_seed_1_mismatch"] = (df["alignment_string"].str[1:7].str.count("0") == 1).astype(int)
    
    df["compensatory_site"] = (df["alignment_string"].str[12:17].str.count("0") == 0).astype(int)
    
    df["supplementary_site"] = (df["alignment_string"].str[12:16].str.count("0") == 0).astype(int)
    df["supplementary_site_2"] = (df["alignment_string"].str[16:21].str.count("0") == 0).astype(int)
    df["empty_seed"] = (df["alignment_string"].str[1:8].str.count("1") == 0).astype(int)
    
    
    df["9_consecutive_match_anywhere"] = (df["alignment_string"]
                                          .str
                                          .contains("1{" + str(9) + ",}")
                                          .astype(int))
    
    
    
    return df


def generate_seed_type_columns(df):
    df['seed_8mer'] = ((df['anchor_a'] == 1) & (df['6mer_seed'] == 1) & (df['match_8'] == 1)).astype(int)
    df['seed_7mer_a1'] = ((df['anchor_a'] == 1) & (df['6mer_seed'] == 1) & (df['match_8'] == 0)).astype(int)
    df['seed_7mer_m8'] = ((df['anchor_a'] == 0) & (df['6mer_seed'] == 1) & (df['match_8'] == 1) & (df['supplementary_site'] == 0) & (df['supplementary_site_2'] == 0)).astype(int)
    df['seed_compensatory'] = ((df['compensatory_site'] == 1) & (df['6mer_seed_1_mismatch'] == 1) & (df['match_8'] == 1)).astype(int)

    df['seed_clash_2'] = ((df['supplementary_site'] == 1) & (df['6mer_seed'] == 1) & (df['match_8'] == 1)).astype(int)
    df['seed_clash_3'] = ((df['supplementary_site_2'] == 1) & (df['6mer_seed'] == 1) & (df['match_8'] == 1)).astype(int)
    df['seed_clash_4'] = ((df['empty_seed'] == 1) & (df['9_consecutive_match_anywhere'] == 1)).astype(int)
    df['seed_clash_5'] = ((df['pred_num_basepairs'] > 10) & (df['6mer_seed'] == 0)).astype(int)

    return df