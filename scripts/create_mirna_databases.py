import pandas as pd
from Bio import SeqIO


def mirna_to_mrna(string):
    return complement(string[::-1])


def complement(string):

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

    df = pd.read_csv("sequences/mirna/targetscan_mirnas.txt",
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
    with open("sequences/mirna/mature_mirnas.fa") as f:

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


targetscan_df = create_targetscan_db()
targetscan_df.to_csv("data/mirna_databases/targetscan.tsv",
                     sep="\t", index=False)

mirbase_df = create_mirbase_db()
mirbase_df.to_csv("data/mirna_databases/mirbase.tsv", sep="\t", index=False)
