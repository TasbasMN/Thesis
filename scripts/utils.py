from collections import defaultdict
from itertools import pairwise, combinations

from Bio import SeqIO
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib_venn as venn


def find_matches(sequence, mirna_df, ignore_first_15_nucleotides=True, find_6mers=False, comparison_columns=False):
    # sourcery skip: low-code-quality
    """function that find matches between a sequence and a mirna dataframe

    Args:
        sequence (str): sequence string
        mirna_df (df): miRNA dataframe
        ignore_first_15_nucleotides (bool, optional): if True, removes matches found in the first 15 nucleotides. Defaults to True.
        find_6mers (bool, optional): if True, finds 6mer matches too. Defaults to False.
    Returns:
        df: df containing matches
    """

    # unpacking mirna_df input
    db_names = mirna_df["name"].values.tolist()
    db_sequences = mirna_df["sequence"].values.tolist()
    db_seeds = [seq[-7:-1] for seq in db_sequences]

    # step 1: seed dict
    seed_dict = defaultdict(list)
    for i in range(len(db_names)):
        seed_dict[db_seeds[i]].append((db_names[i], db_sequences[i]))

    # step 2: creating sequence generator
    generator = sliding_window(dna_to_mrna(sequence), 8)

    # step 3: initializing lists
    names = []
    start_coords = []
    end_coords = []
    seed_matches = []
    match_types = []
    mirna_sequences = []

    # step 4: for loop
    for c, chunk in enumerate(generator, start=1):

        seed = chunk[1:7]
        if seed in seed_dict:

            # pythonic
            m8 = chunk[0]
            anchor_a = chunk[7] == "A"

            c2 = 0
            for record in seed_dict[seed]:
                mirna_sequence = seed_dict[seed][c2][1]
                m8_match = m8 == mirna_sequence[-8]
                c2 += 1

                if anchor_a and m8_match:
                    match_types.append("8mer")
                    start_coords.append(c)
                    end_coords.append(c + 7)
                    seed_matches.append(chunk)

                elif anchor_a:
                    match_types.append("7mer-A1")
                    start_coords.append(c + 1)
                    end_coords.append(c + 7)
                    seed_matches.append(chunk[1:])

                elif m8_match:
                    match_types.append("7mer-m8")
                    start_coords.append(c)
                    end_coords.append(c + 6)
                    seed_matches.append(chunk[:7])

                elif find_6mers:
                    match_types.append("6mer")
                    start_coords.append(c + 1)
                    end_coords.append(c + 6)
                    seed_matches.append(chunk[1:7])
                else:
                    continue

                names.append(record[0])
                mirna_sequences.append(mirna_sequence)

    # step 5: creating result dataframe
    df = pd.DataFrame(
        list(
            zip(
                names,
                start_coords,
                end_coords,
                seed_matches,
                match_types,
                mirna_sequences
            )
        ),
        columns=[
            "name",
            "start",
            "end",
            "seed match",
            "match_type",
            "mirna_sequence"
        ],
    )

    # step 6: extra option for UTR sequences

    if ignore_first_15_nucleotides:

        df = df[df["start"] > 15]

    # step 7: generating new columns
    if comparison_columns:

        df["coordinates"] = [f"{start}-{end}" for start,
                             end in zip(df["start"], df["end"])]

        df["comparison"] = [f"{a}_{b}_{c}" for a,
                            b, c in zip(df.name, df.start, df.end)]

    return df


def import_fasta(fasta):
    """imports fasta into a string

    Args:
        fasta (file): fasta file

    Returns:
        str: sequence string
    """

    with open(fasta) as f:
        records = SeqIO.parse(f, "fasta")
        transcripts = [dna_to_mrna(str(rec.seq)) for rec in records]
        return "".join(transcripts)


def sliding_window(sequence, win_size):
    """function that creates a sliding window
    """
    for i in range(len(sequence) - win_size + 1):
        yield sequence[i: i + win_size]


########################################################################################################################
# string manipulation functions

def complement(string):
    """computes the complement of a nucleic acid sequence
    """
    result = ""
    for nuc in string:
        if nuc == "A":
            result += "U"
        elif nuc == "C":
            result += "G"
        elif nuc == "G":
            result += "C"
        elif nuc in ["U", "T"]:
            result += "A"
    return result


def uracil_to_thymine(string):
    """changes uracils into thymines
    """
    return "".join("T" if nuc == "U" else nuc for nuc in string)


def mirna_to_mrna(string):
    """converts a miRNA string into a mRNA string
    """
    return complement(string[::-1])


def dna_to_mrna(string):
    """converts a dna string into a mrna string
    """
    return "".join("U" if nuc == "T" else nuc for nuc in string)


########################################################################################################################
# db comparison tools


def compare_mirna_dataframes(mirbase, targetscan):
    """function that compares 2 dataframes
    """

    mirbase_copy = mirbase.copy()
    targetscan_copy = targetscan.copy()

    # generating comparison columns
    mirbase_copy["comparison"] = [
        f"{a}_{b}" for a, b in zip(mirbase_copy.name, mirbase_copy.sequence)
    ]
    targetscan_copy["comparison"] = [
        f"{a}_{b}" for a, b in zip(targetscan_copy.name, targetscan_copy.sequence)
    ]

    mirbase_set = set(mirbase_copy["comparison"].values.tolist())
    targetscan_set = set(targetscan_copy["comparison"].values.tolist())

    total0 = len(mirbase_set.union(targetscan_set))
    venn.venn2(
        [mirbase_set, targetscan_set],
        set_labels=("miRBase", "TargetScan"),
        subset_label_formatter=lambda x: str(
            x) + "\n(" + f"{(x/total0):1.0%}" + ")",
    )
    plt.title("Fig 0: Differences between miRNA databases")
    plt.show()


def create_venn_diagrams(mirbase_set, targetscan_set, canonical_set):
    """function that creates venn diagrams
    """

    total1 = len(mirbase_set.union(targetscan_set))
    venn.venn2(
        [mirbase_set, targetscan_set],
        set_labels=("miRBase miRNAs", "TargetScan miRNAs"),
        subset_label_formatter=lambda x: str(
            x) + "\n(" + f"{(x/total1):1.0%}" + ")",
    )
    plt.title(
        "Fig 1: Results using miRBase miRNAs vs. Results using TargetScan miRNAs")
    plt.show()

    total2 = len(mirbase_set.union(canonical_set))
    venn.venn2(
        [mirbase_set, canonical_set],
        set_labels=("miRBase miRNAs", "TargetScan Canonical Results"),
        subset_label_formatter=lambda x: str(
            x) + "\n(" + f"{(x/total2):1.0%}" + ")",
    )
    plt.title("Fig 2: Results using miRBase miRNAs vs. TargetScan Canonical Results")
    plt.show()

    total3 = len(targetscan_set.union(canonical_set))
    venn.venn2(
        [targetscan_set, canonical_set],
        set_labels=("TargetScan miRNAs", "TargetScan Canonical Results"),
        subset_label_formatter=lambda x: str(
            x) + "\n(" + f"{(x/total3):1.0%}" + ")",
    )
    plt.title(
        "Fig 3: Results using TargetScan miRNAs vs. TargetScan Canonical Results")
    plt.show()

    triple = [mirbase_set, targetscan_set, canonical_set]
    venn.venn3(
        triple,
        set_labels=(
            "miRBase miRNAs",
            "TargetScan miRNAs",
            "TargetScan Canonical Results",
        ),
    )
    plt.title("Fig 4: Comparison of all results")
    plt.show()


def plot_venn_diagram(s1, s2, s1_name="set1", s2_name="set2"):
    """function that plots a venn diagram
    """

    total_number = len(s1.union(s2))
    venn.venn2(
        [s1, s2],
        set_labels=(s1_name, s2_name),
        subset_label_formatter=lambda x: str(x)
        + "\n("
        + f"{(x/total_number):1.0%}"
        + ")",
    )
    # plt.title("Fig 3: Results using TargetScan miRNAs vs. TargetScan Canonical Results")
    plt.show()


def parse_targetscan_result_file(f):
    """function that parses a targetscan result file
    """

    df = pd.read_csv(f, sep="\t")

    df.columns = [
        "name",
        "coordinates",
        "match_type",
        "c++_score",
        "c++_score_percentile",
        "weighted_c++_score",
        "conserved_branch_length",
        "pct",
        "relative_kd",
    ]

    df.drop(df.index[df["name"] == "Conserved sites"], inplace=True)
    df.drop(df.index[df["name"] == "Poorly conserved sites"], inplace=True)

    return df


def analyze_differences_between_sets(s1, s2):
    """function that analyzes the differences between two sets
    """

    diff = s1.difference(s2)

    s1_names = []
    s1_start_coords = []
    s1_end_coords = []

    # unpacks previously concatenated mirna names and predicted positions
    for i in diff:
        temp = i.split("_")
        s1_names.append(temp[0])

        temp2 = temp[1].split("-")
        s1_start_coords.append(temp2[0])
        s1_end_coords.append(temp2[1])

    s2_start_coords = []
    s2_end_coords = []

    for s1_name in s1_names:

        # pythonic
        start = targetscan_results[targetscan_results["name"]
                                   == s1_name].start.values
        end = targetscan_results[targetscan_results["name"]
                                 == s1_name].end.values

        s2_start_coords.append(start)
        s2_end_coords.append(end)

    # result df
    column_names = ["name", "canonical_start",
                    "canonical_end", "our_start", "our_end"]
    df = pd.DataFrame(
        list(
            zip(
                s1_names, s1_start_coords, s1_end_coords, s2_start_coords, s2_end_coords
            )
        ),
        columns=column_names,
    )

    # exploding results
    df = df.explode("our_start")
    df = df.explode("our_end")

    # post explode
    df.fillna(0, inplace=True)

    # changing values into int
    df["canonical_start"] = df["canonical_start"].astype(int)
    df["canonical_end"] = df["canonical_end"].astype(int)
    df["our_start"] = df["our_start"].astype(int)
    df["our_end"] = df["our_end"].astype(int)

    # explode causes unwanted coordinate pairs, this part removes them
    # flags rows to be dropped with True
    df["flag"] = abs(df["our_end"] - df["our_start"]) > 10

    # drops True & flag column
    df = df[df["flag"] == False]
    df = df.drop("flag", axis=1)

    df["difference_between_starts"] = df["canonical_start"] - df["our_start"]
    df["difference_between_ends"] = df["canonical_end"] - df["our_end"]

    # additional flags
    df["missed"] = df["our_end"] == 0

    return df


def create_comparison_set(df):
    """function that creates a comparison set
    """

    names = df["name"].values.tolist()
    coords = df["coordinates"].values.tolist()

    # zipping names and coordinates with underscore
    zipped_strings = [f"{m}_{n}" for m, n in zip(names, coords)]

    return set(zipped_strings)


########################################################################################################################
# TargetScan style printing

def pretty_print(df, summary=True):
    """prints a dataframe in monospaced format
    """
    if summary:
        return df.head().style.set_properties(**{'text-align': 'left', 'white-space': 'pre-wrap'}).set_table_styles([dict(selector="", props=[("font-size", "12pt"), ("font-family", 'Courier')])])
    else:
        return df.style.set_properties(**{'text-align': 'left', 'white-space': 'pre-wrap'}).set_table_styles([dict(selector="", props=[("font-size", "12pt"), ("font-family", 'Courier')])])


def print_results(results_df, sequence, summary=True):
    """prints match results in TargetScan style

    Args:
        df (dataframe): df containing find_matches() results
        sequence (str): sequence string

    Returns:
        dataframe: df containing multiline strings in TargetScan style
    """

    names, match_types, mirna_sequences, starts, ends = unpack_results_df(
        results_df)

    result_strings = []

    # creates DNA and miRNA string
    for i, _ in enumerate(names):

        # pythonic
        mirna_sequence = mirna_sequences[i]
        dna_start = (starts[i] - 24) if starts[i] > 25 else 0
        dna_end = ends[i]
        sequence_slice = sequence[dna_start:dna_end]
        whitespace_length = len(sequence_slice) - len(mirna_sequence)

        # creating sequence strings
        dna_string = f"5' {uracil_to_thymine(sequence_slice)} 3' position: {dna_start+1}-{dna_end} of DNA"
        mirna_string = "3' " + whitespace_length * \
            (" ") + complement(mirna_sequences[i]) + " 5' " + names[i]

        # creating middle string
        match_string = ""

        c = -1
        for nucleotide in reversed(mirna_sequence):
            if match_types[i] in ["7mer-m8", "6mer"]:
                match_string += "|" if nucleotide == sequence_slice[c+1] else " "
            else:
                match_string += "|" if nucleotide == sequence_slice[c] else " "
            c -= 1

        # +3 in whitespace_length corresponds to the "5' " part of the sequence strings
        middle_string = (whitespace_length+3)*(" ") + match_string[::-1]
        final_string = f"{dna_string}\n{middle_string}\n{mirna_string}"
        result_strings.append(final_string)

        # df that contains final results
        stylized_df = pd.DataFrame(list(zip(result_strings, match_types)),
                                   columns=["results", "match_types"])

    return pretty_print(stylized_df, summary)


########################################################################################################################
# feature column generator functions

def unpack_results_df(results_df):
    """unpacks results dataframe in a tidy way

    Args:
        results_df (pd.DataFrame): results of find_matches()
    """
    names = results_df.name.tolist()
    mirna_sequences = results_df.mirna_sequence.tolist()
    match_types = results_df.match_type.tolist()
    starts = results_df.start.tolist()
    ends = results_df.end.tolist()

    return names, match_types, mirna_sequences, starts, ends


def generate_3utr_length_column(sequence, results_df):
    """
        generates 3' UTR length column

    Args:
        sequence (str): sequence
        results_df (pd.DataFrame): results dataframe
    """
    results_df["3utr_length"] = len(sequence)
    return results_df


def generate_3_supplementary_pairing_column(sequence, results_df, extended=False):
    """generates 3' supplementary pairing column

    Args:
        sequence (string): DNA sequence under analysis
        results_df (pd.DataFrame): results of find_matches() function
        extended (bool, optional): if extended, 3' supplementary pairing is looked at positions 12-17. If not, positions 13-16 are checked. Defaults to False.

    Returns:
        pd.DataFrame: result dataframe with an extra column
    """

    names, match_types, mirna_sequences, starts, ends = unpack_results_df(
        results_df)

    results = []

    # main loop
    for i in range(len(names)):
        if extended:

            supplementary_start = (
                starts[i] - 11) if match_types[i] in ["8mer", "7mer-m8"] else (starts[i] - 12)
            supplementary_end = starts[i] - 5 if match_types[i] in [
                "8mer", "7mer-m8"] else (starts[i] - 6)

            supplementary_mirna_sequence = mirna_sequences[i][-17:-11]

        else:
            supplementary_start = (
                starts[i] - 10) if match_types[i] in ["8mer", "7mer-m8"] else (starts[i] - 11)
            supplementary_end = starts[i] - 6 if match_types[i] in [
                "8mer", "7mer-m8"] else (starts[i] - 7)

            supplementary_mirna_sequence = mirna_sequences[i][-16:-12]

        if (sequence[supplementary_start:supplementary_end] == supplementary_mirna_sequence):
            results.append(1)
        else:
            results.append(0)

        # supplementary_dna_sequence = (sequence[supplementary_start:supplementary_end])

    if extended:
        results_df["extended_supplementary_pairing"] = results
    else:
        results_df["3_supplementary_pairing"] = results

    return results_df


def create_abundance_dict(results_df):
    """creates a dictionary with miRNA names as keys and their average MRE positions as values, only for miRNAs that have >= 2 MREs

    Args:
        results_df (pd.DataFrame): results from find_matches()

    Returns:
        dict: dictionary with miRNA names as keys and their average MRE positions as values
    """

    names = results_df.name.tolist()
    starts = results_df.start.tolist()
    ends = results_df.end.tolist()

    abundance_dict = {}

    # populating the dict
    for c, mirna in enumerate(names):

        avg_position = int((starts[c] + ends[c]) / 2)

        if mirna in abundance_dict:
            abundance_dict[mirna].append(avg_position)

        else:
            abundance_dict[mirna] = [avg_position]

    # dropping miRNAs having only one match
    # list() is needed because the dictionary is getting changed during the loop
    # if not, "dictionary changed size during iteration" error is raised
    for i in list(abundance_dict):
        if len(abundance_dict[i]) < 2:
            abundance_dict.pop(i)

    return abundance_dict


def generate_avg_position_column(results_df):
    """generates column that contains the average value of the MREs' start & end positions

    Args:
        results_df (pd.DataFrame): output of find_matches()

    Returns:
        pd.DataFrame: result dataframe with an additional column
    """
    results_df["avg_position"] = (
        ((results_df.start + results_df.end) / 2)).astype(int)

    return results_df


def generate_close_proximity_column(results_df):
    """generates column that checks for MREs of the same miRNA in close proximity (13<dist<35)

    Args:
        results_df (pd.DataFrame): output of find_matches()

    Returns:
        pd.DataFrame: result dataframe with an additional column
    """

    # creating abundance dict
    abundance_dict = create_abundance_dict(results_df)

    # creating results matrix
    distances = []
    for mirna in abundance_dict:

        # pythonic
        tmp = abundance_dict[mirna]

        dist_combinations = [abs(a - b)
                             for (a, b) in combinations(tmp, 2)]
        distances.append(dist_combinations)

    # populating true_dict that contains miRNA names that fit the criteria as keys and their distances as values
    true_dict = {}
    for line in distances:
        for i in line:
            if 21 <= i <= 43:
                # 3' UTR abundance checks for MREs that are close to each other (when distance is between 13 and 35 nt)
                # given that we averaged MRE start & end positions and average length of MREs are 4, +8 is added to both
                # 13 and 35.
                x = distances.index(line)
                y = line.index(i)
                mirna_name = list(abundance_dict)[x]
                true_dict[mirna_name] = distances[x][y]

    # finding MRE positions having that specific distance from each other
    results = []
    for name in true_dict:

        # pythonic
        true_distance = true_dict.get(name)
        all_distances = abundance_dict.get(name)

        results.extend(
            [name, [v1, v2]]
            for v1, v2 in pairwise(all_distances)
            if abs(v2 - v1) == true_distance
        )

    # generating zeros column
    results_df["close_proximity"] = 0

    # finding which columns to write "1"
    # ones = []
    for i in results:
        for j in range(len(results[1])):

            # pythonic
            df_filter = (results_df["name"] == i[0]) & (
                results_df["avg_position"] == i[1][j])

            matching_row_index = results_df.loc[df_filter].index.values.astype(int)[
                0]
            results_df.at[matching_row_index, "close_proximity"] = 1
            # ones.append(matching_row_index)

    # # writing "1" to the corresponding columns
    # for i in ones:
    #     results_df.at[i, "3utr_abundance"] = 1
    return results_df


def generate_total_no_of_pairs_column(sequence, results_df):

    names, match_types, mirna_sequences, starts, ends = unpack_results_df(
        results_df)

    results = []

    for i, _ in enumerate(names):

        dna_start = (starts[i] - 24) if starts[i] > 25 else 0
        dna_end = ends[i]
        sequence_slice = sequence[dna_start:dna_end]
        mirna_sequence = mirna_sequences[i]

        no_of_matches = 0

        c = -1
        for nucleotide in reversed(mirna_sequence):
            if (
                match_types[i] in ["7mer-m8", "6mer"]
                and nucleotide == sequence_slice[c + 1]
                or match_types[i] not in ["7mer-m8", "6mer"]
                and nucleotide == sequence_slice[c]
            ):
                no_of_matches += 1
            c -= 1

        results.append(no_of_matches)

    results_df["total_no_of_pairs"] = results

    return results_df


def generate_position_in_utr_column(sequence, results_df):

    utr_length = len(sequence)

    names = results_df["name"].tolist()
    avg_positions = results_df["avg_position"].tolist()

    results = []

    for i, _ in enumerate(names):
        # round to 3 decimal places
        position_in_utr = round(avg_positions[i] / utr_length, 3)
        results.append(position_in_utr)

    results_df["position_in_utr"] = results

    return results_df
