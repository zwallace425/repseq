# Imports
from typing import List
from Bio import SeqIO


def parse_fasta(fasta_file: str) -> List[str]:
    """
    Loads sequences from a .fasta file to a list of strings.

    :param fasta_file: <str> Path to the .fasta file to load the data from.

    :return: <List[str]> .fasta file sequences stored in a list of strings.
    """   
    seq_list = []
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):
        seq = str(seq_record.seq)
        seq_list.append(seq)
        
    return(seq_list)


def write_fasta(
    sequences: List[str],
    output_file: str,
    original_file: str,
):
    """
    :param sequences: <List[str]> List of sequences to be stored.
    :param output_file: <str> File path to write the output to.
    :param original_file: <str> Path to the file where the sequences were extracted from.
                                It is used to recover sequence metadata.
    """
    raise NotImplementedError

