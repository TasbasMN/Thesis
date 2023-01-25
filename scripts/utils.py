from collections import defaultdict

from Bio import SeqIO
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib_venn as venn


def find_matches(sequence, mirna_df, ignore_first_15_nucleotides=True):
    """function that find matches between a sequence and a mirna dataframe

    Args:
        sequence (str): sequence string
        mirna_df (df): miRNA dataframe
        ignore_first_15_nucleotides (bool, optional): if True, removes matches found in the first 15 nucleotides. Defaults to True.

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


def complement(string):
    """function that computes the complement of a string
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
    """function that converts a mirna string into a mrna string
    """

    return complement(string[::-1])


def dna_to_mrna(string):
    """function that converts a dna string into a mrna string
    """
    return "".join("U" if nuc == "T" else nuc for nuc in string)


def sliding_window(sequence, win_size):
    """function that creates a sliding window
    """
    for i in range(len(sequence) - win_size + 1):
        yield sequence[i: i + win_size]


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
    # df = df.drop([df.index[0], df.index[1], df.index[2]])

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


def pretty_print(df):
    """prints a dataframe in monospaced format
    """
    return df.style.set_properties(**{'text-align': 'left', 'white-space': 'pre-wrap'}).set_table_styles([dict(selector="", props=[("font-size", "12pt"), ("font-family", 'Courier')])])


def print_results(df, sequence):
    """prints match results in TargetScan style

    Args:
        df (dataframe): df containing find_matches() results
        sequence (str): sequence string

    Returns:
        dataframe: df containing multiline strings in TargetScan style
    """
    # pythonic
    names = df["name"].tolist()
    match_types = df["match_type"].tolist()
    starts = df["start"].tolist()
    ends = df["end"].tolist()
    mirna_sequences = df["mirna_sequence"].tolist()

    result_strings = []

    # creates DNA and miRNA string

    for i, _ in enumerate(names):

        # for i in range(len(names)):

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
        middle_string = ""

        c = -1
        for nucleotide in reversed(mirna_sequence):
            if match_types[i] == "7mer-m8":
                middle_string += "|" if nucleotide == sequence_slice[c+1] else " "
            else:
                middle_string += "|" if nucleotide == sequence_slice[c] else " "
            c -= 1

        # +3 in whitespace_length corresponds to the "5' " part of the sequence strings
        middle_string = (whitespace_length+3)*(" ") + middle_string[::-1]
        final_string = f"{dna_string}\n{middle_string}\n{mirna_string}"
        result_strings.append(final_string)

        # df that contains final results
        df = pd.DataFrame(list(zip(result_strings, match_types)),
                          columns=["results", "match_types"])

    return pretty_print(df)
