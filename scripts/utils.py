import subprocess
from itertools import pairwise, combinations

import pandas as pd
import numpy as np
import editdistance

from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score


#####################################################
### small functions

def reverse_complement_rna_to_dna(rna_sequence):
    complement = {'A': 'T', 'U': 'A', 'C': 'G', 'G': 'C'}
    reverse_seq = rna_sequence[::-1]
    return ''.join(complement[base] for base in reverse_seq)

def reverse_complement_dna_to_rna(rna_sequence):
    complement = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_seq = rna_sequence[::-1]
    return ''.join(complement[base] for base in reverse_seq)

def calculate_au_content(sequence):
    au_count = sequence.count(
        'A') + sequence.count('T') + sequence.count('U')
    return None if len(sequence) == 0 else au_count / len(sequence)


def slice_column(row):
    return row['full_sequence_of_transcript'][row['temp_start']:row['temp_end']]


def import_pyensembl(grch):
    if grch not in [37, 38]:
        raise ValueError("grch must be either 37 or 38")

    ens_release = 75 if grch == 37 else 109
    
    from pyensembl import EnsemblRelease
    import os
    os.environ['PYENSEMBL_CACHE_DIR'] = "../data"
    assembly = EnsemblRelease(ens_release)
    assembly.download()
    assembly.index()
    
    return assembly



####################################################
# step 1
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


#########################
# step 3

def clash():
    return pd.read_csv("data/processed/clash/clash_parsed.csv")

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



###################### 
# step 5

def generate_ta_sps_columns(df):
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

def generate_mirna_conservation_column(df):
    mirna_df = (pd.read_csv("data/processed/mirbase/mirbase22.csv", 
                            usecols=["accession", "conservation"])
                            .rename(columns={"accession": "mirna_accession", "conservation": "mirna_conservation"})
                            [["mirna_accession", "mirna_conservation"]])

    df = df.merge(mirna_df, on="mirna_accession", how="left")
    return df

def generate_important_sites(df):
    df["anchor_a"] = (df["mre_region"].str[-1] == "A").astype(int)
    df["6mer_seed"] = (
        df["alignment_string"].str[1:7].str.count("0") == 0).astype(int)
    df["match_8"] = (df["alignment_string"].str[7] == "1").astype(int)
    df["6mer_seed_1_mismatch"] = (
        df["alignment_string"].str[1:7].str.count("0") == 1).astype(int)

    df["compensatory_site"] = (
        df["alignment_string"].str[12:17].str.count("0") == 0).astype(int)

    df["supplementary_site"] = (
        df["alignment_string"].str[12:16].str.count("0") == 0).astype(int)
    df["supplementary_site_2"] = (
        df["alignment_string"].str[16:21].str.count("0") == 0).astype(int)
    df["empty_seed"] = (
        df["alignment_string"].str[1:8].str.count("1") == 0).astype(int)

    df["9_consecutive_match_anywhere"] = (df["alignment_string"]
                                          .str
                                          .contains("1{" + str(9) + ",}")
                                          .astype(int))

    return df

def generate_seed_type_columns(df):
    df['seed_8mer'] = ((df['anchor_a'] == 1) & (
        df['6mer_seed'] == 1) & (df['match_8'] == 1)).astype(int)
    df['seed_7mer_a1'] = ((df['anchor_a'] == 1) & (
        df['6mer_seed'] == 1) & (df['match_8'] == 0)).astype(int)
    df['seed_7mer_m8'] = ((df['anchor_a'] == 0) & (df['6mer_seed'] == 1) & (df['match_8'] == 1) & (
        df['supplementary_site'] == 0) & (df['supplementary_site_2'] == 0)).astype(int)
    df['seed_compensatory'] = ((df['compensatory_site'] == 1) & (
        df['6mer_seed_1_mismatch'] == 1) & (df['match_8'] == 1)).astype(int)

    df['seed_clash_2'] = ((df['supplementary_site'] == 1) & (
        df['6mer_seed'] == 1) & (df['match_8'] == 1)).astype(int)
    df['seed_clash_3'] = ((df['supplementary_site_2'] == 1) & (
        df['6mer_seed'] == 1) & (df['match_8'] == 1)).astype(int)
    df['seed_clash_4'] = ((df['empty_seed'] == 1) & (
        df['9_consecutive_match_anywhere'] == 1)).astype(int)
    df['seed_clash_5'] = ((df['pred_num_basepairs'] > 10)
                          & (df['6mer_seed'] == 0)).astype(int)

    return df

def generate_mre_au_content_column(df):
    df["mre_au_content"] = df['mre_region'].apply(calculate_au_content)

    return df

def generate_local_au_content_column(df):

    # clip sets negative start indices to 0
    df["temp_start"] = (df.mre_start - 30).clip(0)

    # this handles mre_ends extending off the mRNA transcript
    df["temp_extended_end"] = df.mre_end + 30
    df["temp_transcript_length"] = df.full_sequence_of_transcript.str.len()
    df["temp_end"] = df[["temp_transcript_length",
                         "temp_extended_end"]].min(axis=1)

    # temp col for calculating au content
    df["temp_col_for_calculating_au_content"] = df.apply(slice_column, axis=1)

    df["local_au_content"] = df['temp_col_for_calculating_au_content'].apply(
        calculate_au_content)

    df.drop(["temp_start", "temp_extended_end", "temp_transcript_length",
            "temp_end", "temp_col_for_calculating_au_content"], axis=1, inplace=True)

    return df

def find_mres_in_close_proximity(df):
    # Calculate the middle of MRE
    df["middle_of_mre"] = ((df["mre_start"] + df["mre_end"]) // 2).astype(int)

    # Find the indices where the length of each group is greater than or equal to 2
    group_indices = df.groupby(["enst", "mirna_accession"]).indices
    filtered_indices = [indices for indices in group_indices.values() if len(indices) >= 2]

    # Create a dictionary to store the results
    result = {}

    # Loop through the filtered indices and perform the computation
    for indices in filtered_indices:
        # Get the values corresponding to the indices
        values = df.iloc[indices]["middle_of_mre"].values

        # Create a mask for the condition
        mask = np.logical_and(
            np.abs(values[:, None] - values) > (13 + 4),
            np.abs(values[:, None] - values) < (35 + 4)
        )

        # Get the indices where the condition is satisfied
        i, j = np.where(mask)

        # Add the pairs to the result dictionary
        for k in range(len(i)):
            # Create a key using the 'enst' and 'mirna_accession' values from the DataFrame
            key = (df.iloc[indices[i[k]]]['enst'], df.iloc[indices[i[k]]]['mirna_accession'])

            # If the key is not in the result dictionary, add it with an empty list as the value
            if key not in result:
                result[key] = []

            # Get the pair values
            a, b = values[i[k]], values[j[k]]

            # Check if the pair (a, b) or (b, a) is already in the list associated with the key
            if (a, b) not in result[key] and (b, a) not in result[key]:
                # Append the pair of values to the list associated with the key
                result[key].append((a, b))

    # Create a set of tuples for faster membership checking
    data_set = {(enst, mirna, i) for (enst, mirna), values in result.items() for i in values[0]}

    # Use vectorized operations to check membership and assign the result
    df['another_mre_in_close_proximity'] = np.where(
        df[['enst', 'mirna_accession', 'middle_of_mre']].apply(tuple, axis=1).isin(data_set),
        1,
        0
    )

    return df



#############################
# step 6 xgb

def report_performance(model, X, y):
    # Make predictions on the input data
    y_pred = model.predict(X)

    # Calculate performance metrics
    accuracy = accuracy_score(y, y_pred)
    precision = precision_score(y, y_pred)
    recall = recall_score(y, y_pred)
    f1 = f1_score(y, y_pred)
    roc_auc = roc_auc_score(y, y_pred)

    return {
        'Accuracy': accuracy,
        'Precision': precision,
        'Recall': recall,
        'F1-Score': f1,
        'ROC AUC': roc_auc,
    }
    
def calculate_metrics(metrics, y_true, y_pred):
    scores = {}
    for metric_name, metric_func in metrics.items():
        score = metric_func(y_true, y_pred)
        scores[metric_name] = score
    return scores

def drop_column_and_score(X, y):
    # List of metrics to calculate
    metrics = {
        'Accuracy': accuracy_score,
        'Precision': precision_score,
        'Recall': recall_score,
        'F1 Score': f1_score
    }

    # Split the data into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y)

    # Train an initial XGBoost model
    model = XGBClassifier()
    model.fit(X_train, y_train)
    initial_scores = calculate_metrics(metrics, y_test, model.predict(X_test))
    scores = {'Initial Model': initial_scores}
    # Drop each column one by one and train a new model
    for column in X.columns:
        # Create a new X without the current column
        X_dropped = X.drop(column, axis=1)

        # Split the modified data into train and test sets
        X_train_dropped, X_test_dropped, y_train, y_test = train_test_split(
            X_dropped, y, test_size=0.2, random_state=42, stratify=y)

        # Train a new XGBoost model without the current column
        model_dropped = XGBClassifier()
        model_dropped.fit(X_train_dropped, y_train)

        # Calculate the scores of the new model
        column_scores = calculate_metrics(
            metrics, y_test, model_dropped.predict(X_test_dropped))
        scores[f'Dropped {column}'] = column_scores

    return pd.DataFrame(scores).T

def stylize_df(df, base_intensity=0.2, difference_intensity=0.8):
    # Define a function to apply color based on cell values
    def apply_color(row):
        color = []  # List to store color values for each cell in the row
        for i, cell in enumerate(row):
            # Check if the cell value is equal to the first cell value in the DataFrame
            if cell == df.iloc[0, i]:
                color.append('background-color: black')  # If equal, set background color to black
            else:
                # Calculate the difference between the cell value and the first cell value
                diff = cell - df.iloc[0, i]
                # Calculate the maximum difference in the DataFrame
                max_diff = df.values.max() - df.values.min()
                if diff > 0:
                    # If the difference is positive, calculate the intensity of green color based on the difference
                    intensity = min(1.0, base_intensity + difference_intensity * (diff / max_diff))
                    color.append(f'background-color: rgba(0, 255, 0, {intensity:.2f})')  # Green with intensity
                else:
                    # If the difference is negative, calculate the intensity of red color based on the difference
                    intensity = min(1.0, base_intensity - difference_intensity * (diff / max_diff))
                    color.append(f'background-color: rgba(255, 0, 0, {intensity:.2f})')  # Red with intensity
        return color

    # Apply the color function to each row of the DataFrame using the `.style.apply()` method
    return df.style.apply(apply_color, axis=1)

def find_most_correlated_features(data, threshold=0.5):
    # Compute the correlation matrix
    correlation_matrix = data.corr().abs()

    # Extract the upper triangle of the correlation matrix
    upper_triangle = correlation_matrix.where(np.triu(np.ones(correlation_matrix.shape), k=1).astype(bool))

    # Find feature pairs with correlation above the threshold
    correlated_pairs = upper_triangle.unstack().sort_values(ascending=False)
    correlated_pairs = correlated_pairs[correlated_pairs > threshold]

    return correlated_pairs


def get_average_correlations(df: pd.DataFrame) -> pd.Series:
    return df.corr().abs().mean()


###########################
#step 7

def get_nucleotides_in_interval(chrom, start, end):
    # sourcery skip: extract-method
    file_path = f"data/fasta/grch37/Homo_sapiens.GRCh37.dna.chromosome.{chrom}.fa"
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
    file_path = f"data/fasta/grch37/Homo_sapiens.GRCh37.dna.chromosome.{chrom}.fa"
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

def augment_cancer_vcf(vcf_df, sequence_offset=30):
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


##############################################
# step 8

def generate_positions_from_id(vcf_df):
    vcf_df['chr'] = vcf_df['id'].str.split('_').str[0]

    vcf_df['start_coordinate'] = vcf_df['id'].str.split('_').str[1].astype(int) - 30 + vcf_df["mrna_start"]
    vcf_df['end_coordinate'] = vcf_df['id'].str.split('_').str[1].astype(int) - 30 + vcf_df["mrna_end"]
    
    return vcf_df

def generate_mre_sequence_for_vcf(vcf_df):

    def slice_column(row):
        return row["mrna_sequence"][row["mre_start"]:row["mre_end"]]
    
    # getting mirna length
    vcf_df["mirna_length"] = vcf_df["mirna_sequence"].str.len()

    # using mirna length to figure out mre coordinates
    vcf_df["mre_end"] = vcf_df["mrna_end"] + vcf_df["mirna_start"]
    vcf_df["mre_start"] = vcf_df["mre_end"] - vcf_df["mirna_length"]

    # some start values might be lower than zero, so we need to adjust
    vcf_df["mre_start"] = vcf_df["mre_start"].apply(lambda x: max(x, 0))

    # creating mre sequence column
    vcf_df["mre_region"] = vcf_df.apply(slice_column, axis=1)

    # dropping temp column
    vcf_df.drop(columns=["mirna_length"], inplace=True)
    
    return vcf_df

def generate_au_content_column_for_vcf(vcf_df):
    
    def calculate_au_content(sequence):
        au_count = sequence.count('A') + sequence.count('T') + sequence.count('U')
        return None if len(sequence) == 0 else au_count / len(sequence)

    vcf_df["au_content_sequence"] = vcf_df.apply(lambda x: get_nucleotides_in_interval(x['chr'], x['start_coordinate']-30, x["end_coordinate"]+30), axis=1)

    vcf_df["local_au_content"] = vcf_df['au_content_sequence'].apply(calculate_au_content)
    
    return vcf_df
