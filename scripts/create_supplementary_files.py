import csv
import pandas as pd
from Bio import SeqIO


def mirna_to_mrna(string):
    return rna_complement(string[::-1])


def rna_complement(string):

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


def create_targetscan_db():

    df = pd.read_csv("data/raw_supplementary_files/targetscan_mirnas.txt",
                     sep="\t", header=0)

    # dropping non-human miRNAs
    df = df[df["Species ID"] == 9606]
    # dropping unnecessary columns
    df = df.drop(["miR family", "Species ID"], axis=1)

    # reordering columns
    df = df.iloc[:, [1, 2, 0, 3, 4]]

    # renaming columns
    renaming_dict = {
        "Seed+m8": "seed",
        "MiRBase ID": "name",
        "Mature sequence": "sequence",
        "Family Conservation?": "conservation",
        "MiRBase Accession": "accession"
    }

    df.rename(columns=renaming_dict, inplace=True)

    df["sequence"] = df["sequence"].apply(mirna_to_mrna)
    df["seed"] = df["seed"].apply(mirna_to_mrna)

    # trimming 7mer seed into 6mer
    df["seed"] = df["seed"].str[1:]

    return df


def create_mirbase_db():
    with open("data/raw_supplementary_files/mature_miRNAs.fa") as f:

        # initializing lists to zip into a df
        names = []
        sequences = []

        for record in SeqIO.parse(f, "fasta"):

            name = str(record.id)
            if name.startswith("hsa"):  # drops all non-human entries
                sequence = str(record.seq)

                # adding name
                names.append(name)

                # adding mRNA sequence
                sequences.append(mirna_to_mrna(sequence))

    return pd.DataFrame(list(zip(names, sequences)), columns=["name", "sequence"])


def read_ta_sps_data(file_path):
    df = pd.read_excel(file_path)

    df = df.drop(["TA percentile rank", "Mean SPS percentile rank",
                 "Mean SPS (kcal/mol)", "miRNA"], axis=1)

    renaming_dict = {"Seed + nt 8": "seed",
                     "7mer-m8 site": "mrna_seed",
                     "TA (log10)": "ta",
                     "6mer SPS (kcal/mol)": "6mer_sps",
                     "7mer-m8 SPS (kcal/mol)": "7mer_sps"
                     }

    df.rename(columns=renaming_dict, inplace=True)

    return df


def parse_clash_data(filename="../data/raw_supplementary_files/mmc1.txt"):
    """
    Parse the CLASH data into 2 csv's
    """
    with open(filename) as f:
        lines = f.readlines()

    columns = []
    data = []
    for line in lines:
        if line.startswith("#"):
            line = line[1:]
            row = line.strip().split("\t")

            columns.append(row)

        else:
            row = line.strip().split("\t")
            data.append(row)

    # removes header text
    columns.pop(0)

    with open("../data/supplementary_files/clash_column_details.tsv", "w") as f:
        csv.writer(f, delimiter="\t").writerows(columns)

    with open("../data/supplementary_files/clash_data_parsed.tsv", "w") as f:
        csv.writer(f, delimiter="\t").writerows(data)


targetscan_df = create_targetscan_db()
targetscan_df.to_csv("data/supplementary_files/targetscan.tsv",
                     sep="\t", index=False)

mirbase_df = create_mirbase_db()
mirbase_df.to_csv("data/supplementary_files/mirbase.tsv",
                  sep="\t", index=False)

ta_sps_df = read_ta_sps_data(
    "data/raw_supplementary_files/bartel_2011_ta_sps_data.xlsx")
ta_sps_df.to_csv("data/supplementary_files/ta_sps.tsv", sep="\t", index=False)

parse_clash_data()
