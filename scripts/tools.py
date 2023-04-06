# this script contains useful standalone functions for general use

import random


### sequence actions

def complement(seq):
    """computes the complement of a nucleic acid seq"""
    complements = str.maketrans("ACTG", "TGAC")
    return seq.translate(complements)

def rna_complement(seq):
    """ complement but with uracil instead of thymine
    """
    complements = str.maketrans("ACTG", "UGAC")
    return seq.translate(complements)

def reverse_complement(seq):
    """computes the reverse complement of a nucleic acid seq"""
    complements = str.maketrans("ACTG", "TGAC")
    return seq.translate(complements)[::-1]

########################################################################################################################
### sequence manipulation in vivo

def mirna_to_mrna(seq):
    """ converts a miRNA string into a mRNA string
    """
    complements = str.maketrans("AUCG", "UAGC")
    return seq.translate(complements)[::-1]

def uracil_to_thymine(seq):
    """changes uracils into thymines
    """
    return seq.replace("U", "T")

def thymine_to_uracil(seq):
    """changes thymines into uracils
    """
    return seq.replace("T", "U")

########################################################################################################################

def seed_match(mirna, gene):
    """quick seed match with built in reverse complement
    """
    alphabet = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    seed = mirna[1:7]
    seed = ''.join([alphabet[s] for s in seed][::-1])
    return seed in gene



def generate_random_seq(length, gc_content=0.5):
    
    gc_nuc_length = round(gc_content*length)
    
    gc_seq = random.choices("GC", k=gc_nuc_length)
    at_seq = random.choices("AT", k=length-gc_nuc_length)
    dna_seq = gc_seq + at_seq
    random.shuffle(dna_seq)
    return ''.join(dna_seq)