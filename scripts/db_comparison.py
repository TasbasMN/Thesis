import matplotlib.pyplot as plt
import matplotlib_venn as venn

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