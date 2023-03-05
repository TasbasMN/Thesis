from collections import defaultdict
from itertools import pairwise, combinations

from Bio import SeqIO
import pandas as pd

import subprocess
import uuid
import os
import contextlib
import math


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
    seeds_without_anchor_a = []
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
                seeds_without_anchor_a.append(chunk[:7])

    # step 5: creating result dataframe
    df = pd.DataFrame(
        list(
            zip(
                names,
                start_coords,
                end_coords,
                seeds_without_anchor_a,
                seed_matches,
                match_types,
                mirna_sequences
            )
        ),
        columns=[
            "name",
            "start",
            "end",
            "seed_without_anchor_a",
            "seed_match",
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


def get_dna_seq_from_coordinates(sequence, start, end):
    """
    Get DNA sequence from coordinates, with indices starting from 1 (biological) instead of 0

    Args:
        sequence (str): DNA sequence
        start (int): start position
        end (int): end position

    Returns:
        str: DNA sequence
    """
    return sequence[start-1:end-1]


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

    names = results_df["name"].values.tolist()
    avg_positions = results_df["avg_position"].values.tolist()

    results = []

    for i, _ in enumerate(names):
        # round to 3 decimal places
        position_in_utr = round(avg_positions[i] / utr_length, 3)
        results.append(position_in_utr)

    results_df["position_in_utr"] = results

    return results_df


def generate_flanking_dinucleotides_columns(sequence, results_df):

    names = results_df["name"].values.tolist()
    starts = results_df["start"].values.tolist()
    ends = results_df["end"].values.tolist()

    results_a5, results_u5, results_g5, results_c5, results_a3, results_u3, results_g3, results_c3 = ([
    ] for _ in range(8))

    for i, _ in enumerate(names):

        a5 = u5 = g5 = c5 = a3 = u3 = g3 = c3 = 0

        nucleotides_5end = get_dna_seq_from_coordinates(
            sequence, starts[i]-3, starts[i]-1)

        nucleotides_3end = get_dna_seq_from_coordinates(
            sequence, ends[i]+1, ends[i]+3)

        for j in nucleotides_5end:
            if j == "A":
                a5 += 1
            elif j == "U":
                u5 += 1
            elif j == "G":
                g5 += 1
            elif j == "C":
                c5 += 1

        for k in nucleotides_3end:
            if k == "A":
                a3 += 1
            elif k == "U":
                u3 += 1
            elif k == "G":
                g3 += 1
            elif k == "C":
                c3 += 1

        results_a5.append(a5)
        results_u5.append(u5)
        results_g5.append(g5)
        results_c5.append(c5)

        results_a3.append(a3)
        results_u3.append(u3)
        results_g3.append(g3)
        results_c3.append(c3)

    results_df["adenine_5end"] = results_a5
    results_df["cytosine_5end"] = results_c5
    results_df["guanine_5end"] = results_g5
    results_df["uracil_5end"] = results_u5

    results_df["adenine_3end"] = results_a3
    results_df["cytosine_3end"] = results_c3
    results_df["guanine_3end"] = results_g3
    results_df["uracil_3end"] = results_u3

    return results_df


def generate_local_au_content_column(sequence, results_df):

    names = results_df["name"].values.tolist()
    starts = results_df["start"].values.tolist()
    ends = results_df["end"].values.tolist()

    # generating weights that correspond to each position in the 30nt window
    weight_dict_5end = {i: (1/(32-i)) for i in range(1, 31)}
    weight_dict_3end = {i: (1/(i+1)) for i in range(1, 31)}

    final_scores_5end = []
    final_scores_3end = []

    # for each miRNA;
    for i, _ in enumerate(names):

        start_5end = (starts[i]-30) if starts[i] > 30 else 1
        end_5end = starts[i]-1

        start_3end = ends[i]+1
        end_3end = ends[i]+31

        # getting flanking 30nt sequences
        window_5end = get_dna_seq_from_coordinates(
            sequence, start_5end, end_5end)
        window_3end = get_dna_seq_from_coordinates(
            sequence, start_3end, end_3end)

        score_5end = sum(
            weight_dict_5end[j + 1]
            for j, nucleotide in enumerate(window_5end)
            if nucleotide in ["A", "U"]
        )
        score_3end = sum(
            weight_dict_3end[k + 1]
            for k, nucleotide in enumerate(window_3end)
            if nucleotide in ["A", "U"]
        )
        final_scores_5end.append(score_5end)
        final_scores_3end.append(score_3end)

    results_df["5end_au_content_score"] = final_scores_5end
    results_df["3end_au_content_score"] = final_scores_3end

    return results_df


def generate_ta_column(results_df, ta_sps_df):

    bartel_seeds = ta_sps_df["mrna_seed"].values.tolist()
    bartel_ta = ta_sps_df["ta"].values.tolist()

    seed_to_ta_dict = {bartel_seeds[i]: bartel_ta[i]
                       for i in range(len(bartel_seeds))}

    results_df["ta"] = results_df["seed_without_anchor_a"].map(seed_to_ta_dict)

    return results_df


def generate_sps_column(results_df, ta_sps_df):

    # unpacking results df
    names = results_df["name"].values.tolist()
    seed_without_anchor_as = results_df["seed_without_anchor_a"].values.tolist(
    )
    match_types = results_df["match_type"].values.tolist()

    # unpacking bartel df
    bartel_seeds = ta_sps_df["mrna_seed"].values.tolist()
    bartel_6mer_sps = ta_sps_df["6mer_sps"].values.tolist()
    bartel_7mer_sps = ta_sps_df["7mer_sps"].values.tolist()

    seeds_to_6mer_sps_dict = {
        bartel_seeds[i]: bartel_6mer_sps[i] for i in range(len(bartel_seeds))}
    seeds_to_7mer_sps_dict = {
        bartel_seeds[i]: bartel_7mer_sps[i] for i in range(len(bartel_seeds))}

    # initiating empty new column
    # results_df["ta_6mer"] = [0 for _ in range(len(results_df))]

    results = [0 for _ in range(len(results_df))]

    for i, _ in enumerate(names):

        if match_types[i] in ["6mer", "7mer-A1"] and seed_without_anchor_as[i] in seeds_to_6mer_sps_dict:
            results[i] = bartel_6mer_sps[i]

        elif match_types[i] in ["7mer-m8", "8mer"] and seed_without_anchor_as[i] in seeds_to_7mer_sps_dict:
            results[i] = bartel_7mer_sps[i]

    results_df["sps"] = results

    return results_df


def one_hot_encode_match_types(results_df):

    match_types = results_df["match_type"].values.tolist()

    m8 = []
    seed_match = [1 for _ in range(len(match_types))]
    a1 = []

    for i, _ in enumerate(match_types):
        if match_types[i] == "8mer":
            m8.append(1)
            a1.append(1)
        elif match_types[i] == "7mer-A1":
            m8.append(0)
            a1.append(1)
        elif match_types[i] == "7mer-m8":
            m8.append(1)
            a1.append(0)
        else:
            m8.append(0)
            a1.append(0)

    results_df["m8"] = m8
    results_df["seedmatch"] = seed_match
    results_df["a1"] = a1

    return results_df


def get_pair_probabilities_of_sequence(sequence, delete_temp_files=True, window_size=80, bp_range=40, length_of_region=16, log10_output=True):

    # define RNAplfold location
    RNAplfold = "/usr/local/bin/RNAplfold"

    # getting unique session ID and adding it to sequence to uniquely rename each output file
    unique_id = f'fn_{str(uuid.uuid1())}'

    seq_with_unique_id = f'>{(unique_id)}' + '\n' + sequence
    seq_with_unique_id = seq_with_unique_id.encode()

    # running RNAplfold
    plfold_subprocess = subprocess.Popen([RNAplfold, "-W", f"{window_size}", "-L", f"{bp_range}", "-u", f"{length_of_region}", "-o"],
                                         stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)

    output, error = plfold_subprocess.communicate(input=seq_with_unique_id)
    plfold_subprocess.wait()

    # parsing output
    with open(f"{unique_id}_lunp", 'r') as f:
        lines = f.readlines()

    # removing whitespaces
    lines = [item.strip() for item in lines]

    # dropping empty lines
    lines = [line for line in lines if len(line) != 0]

    # dropping first 2 lines
    lines = lines[2:]

    # splitting lines by tabs
    lines = [item.split('\t') for item in lines]

    # getting last element of each line
    unpaired_probabilities = [line[-1] for line in lines]

    # deleting temp files
    if delete_temp_files:
        with contextlib.suppress(Exception):
            os.remove(f"{unique_id}_basepairs")
            os.remove(f"{unique_id}_lunp")

    if not log10_output:
        return unpaired_probabilities

    log10_results = []
    for i in unpaired_probabilities:
        if i == "NA":
            log10_results.append(0)
        else:
            log10_results.append(round(math.log10(float(i)), 3))

    return log10_results


def get_accessibility_column(sequence, results_df, log10_results=True):

    ends = results_df["end"].values.tolist()

    log10_results = get_pair_probabilities_of_sequence(
        sequence, length_of_region=16, log10_output=log10_results)

    results = []
    # for each miRNA;
    for i, _ in enumerate(ends):

        end = ends[i]

        results.append(log10_results[end-1])

    results_df["accessibility"] = results

    return results_df
