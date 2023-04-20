# this file contains feature column generator functions
import re
from itertools import pairwise, combinations

import pandas as pd


######################################################################################################################
### CLASH type finder functions working with 2 regexes

def find_2_regex_matches(df, regex1, regex2, new_column_name):
    """
    Finds rows in a DataFrame that match both regex patterns and adds a new column 
    with the results.

    Args:
        df: pandas DataFrame
            DataFrame to search for matches.
        regex1: str
            Regular expression pattern to match.
        regex2: str
            Regular expression pattern to match.
        new_column_name: str
            Name of the new column to add with the match results.

    Returns:
        pandas DataFrame
            The input DataFrame with the new column added.

    Raises:
        None
    """
    regex1_compiled = re.compile(regex1)
    regex2_compiled = re.compile(regex2)

    matches = []
    flag_column = df["flag_column"].tolist()
    alignment_string = df["alignment_string"].tolist()

    for i, string in enumerate(alignment_string):
        if flag_column[i] == 1:
            # Skip this row and set the match to 0
            matches.append(0)
            continue

        # Check for regex matches
        match1 = regex1_compiled.search(string)
        match2 = regex2_compiled.search(string)
        matches.append(1 if match1 and match2 else 0)

        # If both regexes match, append 1 to the flag column for this row
        if matches[i] == 1:
            df.at[i, "flag_column"] = 1

    df[new_column_name] = matches

    return df


def find_CLASH_II_sites(df):
    """
    Finds CLASH II sites in a dataframe by looking for 7mer-m8 patterns using regex_pattern1 and
    matches in positions 13-16 using regex_pattern2. Adds a new column named 'CLASH_II' to the 
    dataframe with the resulting matches.

    Args:
    - df: pandas DataFrame with data to search for CLASH II sites

    Returns:
    - pandas DataFrame: original DataFrame with a new 'CLASH_II' column containing CLASH II site matches
    """
    # looks for 7mer-m8
    regex_pattern1 = r"1{7}0$"
    # looks for matches in positions 13-16
    regex_pattern2 = r"1{4}[01]{12}$"
    new_column_name = "CLASH_II"

    return find_2_regex_matches(df, regex_pattern1, regex_pattern2, new_column_name)


def find_CLASH_III_sites(df):
    """Finds CLASH-III sites in a DataFrame of sequences.

    Args:
        df (pandas.DataFrame): A DataFrame with a "sequence" column containing
            DNA sequences to search for CLASH-III sites.

    Returns:
        pandas.DataFrame: The original DataFrame with an additional "CLASH_III"
        column containing 1 where a CLASH-III site was found and 0 otherwise.

    Notes:
        A CLASH-III site is defined as a sequence that matches both of the following
        two regex patterns:
        - r"1{7}0$": a sequence that ends with 7 consecutive 1's and a 0
        - r"1{4}[01]{16}$": a sequence that has 4 consecutive 1's followed by 16
          positions that can be either 0 or 1, and ends with either 0 or 1.
        The function uses the helper function `find_2_regex_matches` to find
        matches for both patterns in the "sequence" column of the input DataFrame.
    """

    # looks for 7mer-m8
    regex_pattern1 = r"1{7}0$"
    # looks for matches in positions 17-21
    regex_pattern2 = r"1{4}[01]{16}$"
    new_column_name = "CLASH_III"

    return find_2_regex_matches(df, regex_pattern1, regex_pattern2, new_column_name)


def find_CLASH_IV_sites(df):
    """Finds CLASH-IV sites in a DataFrame of sequences.

    Args:
        df (pandas.DataFrame): A DataFrame with a "sequence" column containing
            DNA sequences to search for CLASH-IV sites.

    Returns:
        pandas.DataFrame: The original DataFrame with an additional "CLASH_IV"
        column containing 1 where a CLASH-IV site was found and 0 otherwise.
    """
    # looks for non matching seed
    regex_pattern1 = r"0{7}0$"
    # looks for 9 consecutive matches
    regex_pattern2 = r"1{9}"
    new_column_name = "CLASH_IV"

    return find_2_regex_matches(df, regex_pattern1, regex_pattern2, new_column_name)


def find_compensatory_sites(df):
    """Finds compensatory sites in a DataFrame of sequences.

    Args:
        df (pandas.DataFrame): A DataFrame with a "sequence" column containing
            DNA sequences to search for compensatory sites.

    Returns:
        pandas.DataFrame: The original DataFrame with an additional "CLASH_IV"
        column containing 1 where a compensatory site was found and 0 otherwise.
    """
    # looks for 7mer with 1 mismatch
    regex_pattern1 = r"11{5}0{1}0$|11{4}0{1}1{1}0$|11{3}0{1}1{2}0$|11{2}0{1}1{3}0$|11{1}0{1}1{4}0$|11{0}0{1}1{5}0$|10{1}1{5}0$"
    # looks for matches in positions 13-16
    regex_pattern2 = r"1{4}[01]{12}$"
    new_column_name = "compensatory"

    return find_2_regex_matches(df, regex_pattern1, regex_pattern2, new_column_name)

# CLASH type finder functions working with 1 regex
def find_1_regex_matches(df, regex_pattern, new_column_name):
    """
    Finds matches of a regular expression pattern in a DataFrame column and
    adds a new column with the results.

    Args:
        df (pandas.DataFrame): The DataFrame to operate on.
        regex_pattern (str): The regular expression pattern to search for.
        new_column_name (str): The name of the new column to add.

    Returns:
        pandas.DataFrame: The modified DataFrame.

    Example:
        >>> df = pd.DataFrame({'flag_column': [0, 1, 0],
        ...                    'alignment_string': ['abc', 'def', 'ghi']})
        >>> find_1_regex_matches(df, r'[a-z]', 'matches')
           flag_column alignment_string  matches
        0            0              abc        1
        1            1              def        0
        2            0              ghi        1
    """
    regex = re.compile(regex_pattern)
    matches = []
    flag_column = df["flag_column"].tolist()
    alignment_string = df["alignment_string"].tolist()

    for i, string in enumerate(alignment_string):

        # Check if the flag column for this row is 1
        if flag_column[i] == 1:
            # Skip this row and set the match to 0
            matches.append(0)

        else:
            # Check for regex match
            matches.append(1 if regex.search(string) else 0)

            # If regex matches, append 1 to the flag column for this row
            if matches[i] == 1:
                df.at[i, "flag_column"] = 1

    df[new_column_name] = matches

    return df


def find_8mer_sites(df):
    """Find 8-mer sites in a DataFrame column.

    Args:
        df (pandas.DataFrame): The DataFrame to search for 8-mer sites.

    Returns:
        pandas.DataFrame: A new DataFrame with an additional column named "8mer",
        which contains the 8-mer sites found in the input DataFrame.

    """
    # looks for 8mers
    regex_pattern = r"1{8}$"
    new_column_name = "8mer"

    return find_1_regex_matches(df, regex_pattern, new_column_name)


def find_7mer_a1_sites(df):
    """
    Finds 7mer-A1 sites in a given DataFrame.

    Args:
    - df: pandas DataFrame

    Returns:
    - pandas DataFrame with new column "7mer-a1" containing the matches
    """
    regex_pattern = r"01{7}$"
    new_column_name = "7mer-a1"

    return find_1_regex_matches(df, regex_pattern, new_column_name)



def find_7mer_m8_sites(df):
    """
    Finds 7mer-m8 sites in a given dataframe by applying a regex pattern.

    Parameters:
    df (pandas.DataFrame): The input dataframe to search for matches.

    Returns:
    pandas.DataFrame: A copy of the input dataframe with a new column added, called '7mer-m8',
    that contains the matches found by applying the regex pattern '1{7}0$'.
    """
    # looks for 7mer-m8s
    regex_pattern = r"1{7}0$"
    new_column_name = "7mer-m8"

    return find_1_regex_matches(df, regex_pattern, new_column_name)


def find_seed_with_1_mismatch(df):
    """Find seeds with exactly 1 mismatch in a given DataFrame.

    Args:
        df (pandas.DataFrame): The input DataFrame to search for seeds.

    Returns:
        pandas.DataFrame: A new DataFrame with a new column added
            indicating which rows contain a seed with exactly 1 mismatch.
    """
    # looks for seeds with 1 mismatch
    regex_pattern = r"11{5}0{1}0$|11{4}0{1}1{1}0$|11{3}0{1}1{2}0$|11{2}0{1}1{3}0$|11{1}0{1}1{4}0$|11{0}0{1}1{5}0$|10{1}1{5}0$"
    new_column_name = "seed_with_1_mismatch"

    return find_1_regex_matches(df, regex_pattern, new_column_name)


def find_centered_site(df):
    """
    Returns a new DataFrame with an additional column called 'centered_site'
    that contains sequences of centered sites. 

    Args:
        df (pandas.DataFrame): The DataFrame to search for centered sites.

    Returns:
        pandas.DataFrame: A new DataFrame with an additional column called
        'centered_site' that contains the centered site sequences that match
        the given regex pattern.
    """
    # looks for seeds with 1 mismatch
    regex_pattern = r"1{11}[01]{4}$|1{11}[01]{3}$"
    new_column_name = "centered_site"

    return find_1_regex_matches(df, regex_pattern, new_column_name)


def find_CLASH_V_sites(df, threshold=11):
    """
    Find CLASH V sites in a DataFrame.

    Args:
        df (pandas.DataFrame): DataFrame containing the columns "flag_column" and "no_of_base_pairs".
        threshold (int): Threshold value for the "no_of_base_pairs" column.

    Returns:
        pandas.DataFrame: The input DataFrame with an additional "CLASH_V" column that contains 1 if the
        corresponding row meets the threshold condition and the "flag_column" value is None, and 0 otherwise.
    """

    mask = (df["no_of_base_pairs"] >= threshold) & df["flag_column"].isnull()

    df["CLASH_V"] = mask.astype(int)

    return df


def create_midpoint_dict(df):
    """creates a dictionary with miRNA names as keys and their average MRE positions as values, only for miRNAs that have >= 2 MREs

    Args:
        df (pd.DataFrame): results from find_matches()

    Returns:
        dict: dictionary with miRNA names as keys and their average MRE positions as values
    """
    # calculate midpoint for each MRE
    df["midpoint"] = df["start"] + 4

    # group by miRNA name and collect midpoints for each group in a list
    grouped = df.groupby("name")["midpoint"].apply(list)

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
    df["midpoint"] = df["start"] + df["alignment_string"].str.len() / 2

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
    df_filter = df.apply(lambda row: any((row["name"] == i[0]) and (row["midpoint"] == j)
                                         for i in results for j in i[1]), axis=1)

    # Set the corresponding values to 1.
    df.loc[df_filter, "close_proximity"] = 1

    return df



def generate_TA_column(df):
    """
    Adds TA values to every seed in result dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing the 'seed' column.

    Returns
    -------
    pandas.DataFrame
        A new DataFrame with the 'ta' column merged from the TA dataset.
    """
    ta_df = pd.read_csv("../data/supplementary_files/ta_sps.tsv", sep="\t", usecols=["mrna_seed", "ta"])

    ta_df = ta_df.rename(columns={"mrna_seed": "seed"})
    
    return pd.merge(df, ta_df, on="seed", how="left")


def generate_SPS_column(df):
    """
    Given a Pandas DataFrame containing miRNA target predictions, adds a column to the DataFrame indicating the SPS score
    for each target site. The SPS score is a measure of the strength of the interaction between the miRNA and the target mRNA.
    This function reads the SPS scores from a TSV file and merges them with the input DataFrame based on the seed sequence
    of each target site.

    Args:
        df (pandas.DataFrame): A DataFrame containing miRNA target predictions with a "seed" column.

    Returns:
        pandas.DataFrame: A copy of the input DataFrame with an additional "7mer_sps" column indicating the SPS score
        for each target site, or None if the site does not have a valid SPS score.
    """

    # Read in the SPS scores from a TSV file
    sps_df = pd.read_csv("../data/supplementary_files/ta_sps.tsv", sep="\t", usecols=["mrna_seed", "7mer_sps"])

    # Rename the "mrna_seed" column to match the "seed" column in the input DataFrame
    sps_df = sps_df.rename(columns={"mrna_seed": "seed"})

    # Filter the input DataFrame to only include rows where "7mer-m8" or "7mer-a1" is 1
    mask = (df["7mer-m8"] == 1) | (df["7mer-a1"] == 1)

    # Merge the input DataFrame with the SPS scores based on the "seed" column
    merged_df = pd.merge(df, sps_df, on="seed", how="left")

    # Set the "7mer_sps" value to None for rows that do not meet the filter criteria
    merged_df.loc[~mask, "7mer_sps"] = None

    return merged_df


######################################################################################################################
### driver function that generates CLASH type columns

def find_clash_types(df, find_non_CLASH_types=False, drop_flag_column=True):
    """
    Finds CLASH sites in the input DataFrame `df`.

    Parameters
    ----------
    df : pandas.DataFrame
        The input DataFrame to search for CLASH sites.
    find_non_CLASH_types : bool, optional
        If True, also finds non-CLASH sites in the DataFrame, by calling
        `find_compensatory_sites`, `find_seed_with_1_mismatch`, and
        `find_centered_site`. Defaults to False.
    drop_flag_column : bool, optional
        If True, drops the 'flag_column' column from the DataFrame before
        returning it. Defaults to True.

    Returns
    -------
    pandas.DataFrame
        The input DataFrame with additional columns indicating the presence of
        different types of CLASH sites (see `find_8mer_sites`, `find_7mer_a1_sites`,
        `find_7mer_m8_sites`, `find_CLASH_II_sites`, `find_CLASH_III_sites`,
        `find_CLASH_IV_sites`, and `find_CLASH_V_sites` for details).
    """
    df["flag_column"] = None

    df = find_8mer_sites(df)
    df = find_7mer_a1_sites(df)
    df = find_7mer_m8_sites(df)
    df = find_CLASH_II_sites(df)
    df = find_CLASH_III_sites(df)
    df = find_CLASH_IV_sites(df)
    df = find_CLASH_V_sites(df)

    if find_non_CLASH_types:
        df = find_compensatory_sites(df)
        df = find_seed_with_1_mismatch(df)
        df = find_centered_site(df)

    if drop_flag_column:
        df = df.drop("flag_column", axis=1)

    return df

######################################################################################################################
