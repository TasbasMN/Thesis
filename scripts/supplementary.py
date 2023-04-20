 
import csv
import os

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


def parse_clash_data(filename="data/raw_supplementary_files/mmc1.txt"):
    """
    Parse Clash data from the given file and write the results to output files.

    Args:
        filename (str): The path to the input file. Default is "data/raw_supplementary_files/mmc1.txt".

    Raises:
        FileNotFoundError: If the input file does not exist.

    Returns:
        None.
    """

    with open(filename) as f:
        
        # iter(f) is the same as (line for line in f)
        lines = iter(f)

        columns = []
        data = []
        for line in lines:
            if line.startswith("#"):
                line = line[1:]
                row = next(csv.reader([line], delimiter='\t'))
                columns.append(row)
            else:
                row = next(csv.reader([line], delimiter='\t'))
                data.append(row)

        # removes header text
        columns.pop(0)

        # create the output directory if it doesn't exist
        output_dir = "data/supplementary_files"
        os.makedirs(output_dir, exist_ok=True)

        # write the output CSV files
        with open(os.path.join(output_dir, "clash_column_details.tsv"), "w") as f1, \
             open(os.path.join(output_dir, "clash.tsv"), "w") as f2:
            writer1 = csv.writer(f1, delimiter='\t')
            writer2 = csv.writer(f2, delimiter='\t')
            writer1.writerows(columns)
            writer2.writerows(data)
        

parse_clash_data()
        

parse_clash_data()
