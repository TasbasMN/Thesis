# this module contains standalone nucleotide manipulation functions

import random

########################################################################################################################
# complements


def complement(seq):
    """
    Returns the complement of a DNA sequence.

    Args:
        seq (str): A DNA sequence.

    Returns:
        str: The complement DNA sequence.

    Examples:
        >>> complement("ACTG")
        'TGAC'
    """
    complements = str.maketrans("ACTG", "TGAC")
    return seq.translate(complements)


def rna_complement(seq):
    """
    Translates a DNA sequence into its RNA complement.

    Args:
        seq (str): A DNA sequence to be translated.

    Returns:
        str: The RNA complement of the input sequence.

    Examples:
        >>> rna_complement("TACGGT")
        'AUGCCA'
    """
    complements = str.maketrans("ACTG", "UGAC")
    return seq.translate(complements)


def reverse_complement(seq):
    """
    Returns the reverse complement of a DNA sequence.

    Arguments:
    - seq -- a string representing a DNA sequence

    Returns:
    - a string representing the reverse complement of the input sequence

    Example:
    >>> reverse_complement("ATCG")
    'CGAT'
    """
    complements = str.maketrans("ACTG", "TGAC")
    return seq.translate(complements)[::-1]


########################################################################################################################
# thematic sequence manipulation


def mirna_to_mrna(string):
    """
    Convert a miRNA sequence to its complementary mRNA sequence.

    Args:
    - string (str): The miRNA sequence to convert.

    Returns:
    - str: The complementary mRNA sequence.

    Example usage:
    >>> mirna_to_mrna('UGAGGUAGUAGGUUGUAUAGUU')
    'UAUACAACCACUACUCCAUCA'
    """
    complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
    return ''.join(complement.get(base, base) for base in string)[::-1]


def uracil_to_thymine(seq):
    """
    Convert all uracil nucleotides (U) in a DNA or RNA sequence to thymine (T).

    Args:
        seq (str): The DNA or RNA sequence to convert.

    Returns:
        str: The converted sequence with all U's replaced by T's.
    """
    return seq.replace("U", "T")


def thymine_to_uracil(seq):
    """
    Convert a DNA sequence to an RNA sequence by replacing all occurrences of
    thymine (T) with uracil (U).

    Args:
        seq (str): The DNA sequence to convert.

    Returns:
        str: The RNA sequence with all T's replaced with U's.
    """
    return seq.replace("T", "U")

#######################################################################################################################
# misc

def get_dna_seq_from_coordinates(sequence, start, end):
    """
    Get DNA sequence from coordinates, with indices starting from 1 (biological) instead of 0

    Args:
    - sequence (str): The input DNA sequence.
    - start (int): The starting coordinate (1-indexed) for the DNA sequence to be extracted.
    - end (int): The ending coordinate (1-indexed, inclusive) for the DNA sequence to be extracted.

    Returns:
    - str: The DNA sequence from `start` to `end` coordinates within the input `sequence`.
    """
    return sequence[start-1:end-1]

def seed_match(mirna, sequence):
    """
    Check if a given miRNA seed sequence matches a target sequence.

    Args:
        mirna (str): The miRNA seed sequence to match.
        sequence (str): The target sequence to match against.

    Returns:
        bool: True if the miRNA seed sequence matches the target sequence, False otherwise.
    """
    alphabet = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    seed = mirna[1:7]
    seed = ''.join([alphabet[s] for s in seed][::-1])
    return seed in sequence


def generate_random_seq(length, gc_content=0.5):
    """
    Generate a random DNA sequence of the given length and GC content.

    Args:
        length (int): The length of the sequence.
        gc_content (float): The desired GC content, between 0 and 1.

    Returns:
        str: The randomly generated sequence.

    Example:
        >>> generate_random_seq(10, 0.6)
        'GCGCTAGCTA'
    """
    gc_nuc_length = round(gc_content * length)
    gc_seq = random.choices("GC", k=gc_nuc_length)
    at_seq = random.choices("AT", k=length - gc_nuc_length)
    dna_seq = gc_seq + at_seq
    random.shuffle(dna_seq)
    return "".join(dna_seq)
