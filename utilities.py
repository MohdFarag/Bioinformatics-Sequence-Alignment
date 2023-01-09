#!/usr/bin/env python

# Import libraries
import os
import textwrap
# BioPython
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
# Math
import numpy as np
import pandas as pd
# Plotting
from matplotlib.colors import ListedColormap

#####################################################################
# Global variables
#####################################################################

# Reference of sequences letters:
# http://web.mit.edu/meme_v4.11.4/share/doc/alphabets.html
LETTERS_OF_DNA = ['A', 'G', 'T', 'C', 'N', 'X']
LETTERS_OF_RNA = ['A', 'G', 'U', 'C', 'N', 'X']
LETTERS_OF_PROTEIN = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X', 'B', 'Z', 'J']

OUTPUT_LOC = r"./src/output.fasta"
INPUT_LOC = r"./src/input.fasta"

#####################################################################
# .fasta file
#####################################################################

# Read .fasta file
def read_fasta(path: str, concat: bool = False):
    """Reads a fasta file and returns a dictionary with the sequences

    Args:
        path (str): path to the fasta file
        concat (bool, optional): If you want to concatenate list. Defaults to False.

    Returns:
        list: list of sequences
    """
    try:
        with open(path, "r", encoding='utf-8-sig') as file:
            # Read the text from a file
            sequences = file.read()
            sequences = sequences.split(">")
    except Exception as e:
        if "No such file or directory" in e.__str__():
            raise Exception("**ERROR** No such file or directory")

    sequences_list = {}
    for all_sequence in sequences:
        lines = all_sequence.splitlines()
        if len(lines) != 0:
            if concat == True:
                sequence = "".join(lines[1:-1])
            else:
                sequence = lines[1:-1]

            sequences_list[lines[0]] = sequence

    return sequences_list

# Write .fasta file
def write_fasta(sequences: dict):
    """
    Takes a dictionary and writes it to a fasta file
    Must specify the filename when calling the function
    """

    with open(INPUT_LOC, "w") as outfile:
        for key, value in sequences.items():
            outfile.write(f">{key}" + "\n")
            outfile.write("\n".join(textwrap.wrap(value, 60)))
            outfile.write("\n")

#####################################################################
# Check sequence validity
#####################################################################
def check_sequence(sequence: str, type_of_sequence: str):
    """
    Check validity of sequences
    """
    sequence = sequence.upper()
    type_of_sequence = type_of_sequence.upper()

    if type_of_sequence == "DNA":
        return is_dna(sequence)
    elif type_of_sequence == "RNA":
        return is_rna(sequence)
    elif type_of_sequence == "PROTEIN":
        return is_protein(sequence)

    return False

def is_dna(sequence):
    """
    Check if sequence is DNA
    """

    for letter in sequence:
        if letter not in LETTERS_OF_DNA:
            return False
    return True

def is_rna(sequence):
    """
    Check if sequence is RNA
    """

    for letter in sequence:
        if letter not in LETTERS_OF_RNA:
            return False
    return True

def is_protein(sequence):
    """
    Check if sequence is protein
    """
   
    for letter in sequence:
        if letter not in LETTERS_OF_PROTEIN:
            return False
    return True

#####################################################################
# Alignments
#####################################################################

# Pairwise global alignment (Needleman-Wunsch algorithm)
def pairwise_global_alignment(sequence_a:str, sequence_b:str, match:int=1, mismatch: int=0, gap:int=-1):
    """
    The Needleman-Wunsch Algorithm
    ==============================
    This is a dynamic programming algorithm for finding the optimal global alignment of
    two sequences.

    Arguments
    -------
    sequence_a: First sequence
    sequence_b: Second sequence
    match: Matching score
    mismatch: Mismatching score
    gap: Gapping score

    Returns
    -------
    results: dict
    results = {
        "matrix":matching matrix,
        "color: color matrix,
        "score": score,
        "alignments": alignments
    }

    Example
    -------
        >>> x = "GATTACA"
        >>> y = "GCATGCU"
        >>> results = pairwise_global_alignment(x, y)
        >>> results["score"]
        3.0
        >>> results["matrix"]
        [[ 0., -1., -2., -3., -4., -5., -6., -7.],
         [-1.,  1.,  0., -1., -2., -3., -4., -5.],
         [-2.,  0.,  1.,  1.,  0., -1., -2., -3.],
         [-3., -1.,  0.,  1.,  2.,  1.,  0., -1.],
         [-4., -2., -1.,  0.,  2.,  2.,  1.,  0.],
         [-5., -3., -2.,  0.,  1.,  2.,  2.,  1.],
         [-6., -4., -2., -1.,  0.,  1.,  3.,  2.],
         [-7., -5., -3., -1., -1.,  0.,  2.,  3.]]
        >>> results["alignments"]
        [[['G', 'A', 'T', 'T', 'A', 'C', 'A'], 
          ['G', 'C', 'A', 'T', 'G', 'C', 'U']]]
    """

    if isinstance(sequence_a, list):
        sequence_a = "".join(sequence_a)

    if isinstance(sequence_b, list):
        sequence_b = "".join(sequence_b)

    if isinstance(match, str):
        match = int(match)

    if isinstance(mismatch, str):
        mismatch = int(mismatch)

    if isinstance(gap, str):
        gap = int(gap)

    # STEP 1: Initialization of matrix
    match_matrix = np.zeros((len(sequence_a)+1, len(sequence_b)+1))
    color_matrix = np.zeros((len(sequence_a)+1, len(sequence_b)+1))

    for i in range(1, len(sequence_a)+1):
        match_matrix[i, 0] = match_matrix[i-1, 0] + gap

    for j in range(1, len(sequence_b)+1):
        match_matrix[0, j] = match_matrix[0, j-1] + gap

    # STEP 2: Filling to the matrix
    for i in range(1, len(sequence_a)+1):
        for j in range(1, len(sequence_b)+1):
            arr = []
            if sequence_a[i-1] == sequence_b[j-1]:
                arr.append(match_matrix[i-1, j-1] + match)
            elif sequence_a[i-1] != sequence_b[j-1]:
                arr.append(match_matrix[i-1, j-1] + mismatch)

            arr.append(match_matrix[i-1, j] + gap)
            arr.append(match_matrix[i, j-1] + gap)
            match_matrix[i, j] = max(arr)

    score = match_matrix[-1, -1]

    # STEP 3: Backtracing
    alignments = []
    i, j = match_matrix.shape[0] - 1, match_matrix.shape[1] - 1

    alignment_a = []
    alignment_b = []
    while not (i == 0 and j == 0):
        letter_a = sequence_a[i-1]
        letter_b = sequence_b[j-1]
        curr = match_matrix[i, j]
        color_matrix[i, j] = 1

        corner = match_matrix[i-1, j-1]
        top = match_matrix[i-1, j]
        left = match_matrix[i, j-1]
        if curr - match == corner and sequence_a[i-1] == sequence_b[j-1] and i != 0 and j != 0:
            alignment_a.append(letter_a)
            alignment_b.append(letter_b)
            i -= 1
            j -= 1
        elif curr - mismatch == corner and sequence_a[i-1] != sequence_b[j-1] and i != 0 and j != 0:
            alignment_b.append(letter_b)
            i -= 1
            j -= 1
        elif curr - gap == top and i != 0:
            alignment_a.append(letter_a)
            alignment_b.append("-")
            i -= 1

        elif curr - gap == left and j != 0:
            alignment_a.append("-")
            alignment_b.append(letter_b)
            j -= 1

    color_matrix[i, j] = 1
    alignments.append([alignment_a[::-1], alignment_b[::-1]])

    results = {
        "matrix": match_matrix,
        "color": color_matrix,
        "score": score,
        "alignments": alignments
    }

    return results

# Pairwise Local Alignment (smith-waterman algorithm)
def pairwise_local_alignment(sequence_a:str, sequence_b:str, match:int=1, mismatch: int=0, gap:int=-1):
    """The Smith-Waterman Algorithm
    ==============================
    This is a dynamic programming algorithm for finding the optimal local alignment of
    two sequences.

    Arguments
    -------
    sequence_a: First sequence
    sequence_b: Second sequence
    match: Matching score
    mismatch: Mismatching score
    gap: Gapping score

    Returns
    -------
    results: dict
    results = {
        "matrix":matching matrix,
        "color: color matrix,
        "score": score,
        "alignments": alignments
    }

    Example
    -------
        >>> x = "GATTACA"
        >>> y = "GCATGCU"
        >>> results = pairwise_global_alignment(x, y)
        >>> results["score"]
        3.0
        >>> results["matrix"]
        [[0., 0., 0., 0., 0., 0., 0., 0.],
         [0., 1., 0., 0., 0., 1., 0., 0.],
         [0., 0., 1., 1., 0., 0., 1., 0.],
         [0., 0., 0., 1., 2., 1., 0., 1.],
         [0., 0., 0., 0., 2., 2., 1., 0.],
         [0., 0., 0., 1., 1., 2., 2., 1.],
         [0., 0., 1., 0., 1., 1., 3., 2.],
         [0., 0., 0., 2., 1., 1., 2., 3.]]
        >>> results["alignments"]
        [[['G', 'A', 'T', 'T', 'A', 'C', 'A'],
          ['G', 'C', 'A', 'T', 'G', 'C', 'U']]]
    """

    if isinstance(sequence_a, list):
        sequence_a = "".join(sequence_a)

    if isinstance(sequence_b, list):
        sequence_b = "".join(sequence_b)

    if isinstance(match, str):
        match = int(match)

    if isinstance(mismatch, str):
        mismatch = int(mismatch)

    if isinstance(gap, str):
        gap = int(gap)

    # STEP 1: Initialization of matrix
    match_matrix = np.zeros((len(sequence_a)+1, len(sequence_b)+1))
    color_matrix = np.zeros((len(sequence_a)+1, len(sequence_b)+1))

    # STEP 2: Filling to the matrix
    for i in range(1, len(sequence_a)+1):
        for j in range(1, len(sequence_b)+1):
            arr = [0]
            if sequence_a[i-1] == sequence_b[j-1]:
                arr.append(match_matrix[i-1, j-1] + match)
            elif sequence_a[i-1] != sequence_b[j-1]:
                arr.append(match_matrix[i-1, j-1] + mismatch)

            arr.append(match_matrix[i-1, j] + gap)
            arr.append(match_matrix[i, j-1] + gap)
            match_matrix[i, j] = max(arr)
    score = match_matrix.max()

    # STEP 3: Backtracing
    alignments = []

    index_max = np.unravel_index(
        np.argmax(match_matrix, axis=None), match_matrix.shape)
    i, j = index_max[0], index_max[1]
    alignment_a = []
    alignment_b = []

    while match_matrix[i, j] != 0:
        curr = match_matrix[i, j]
        color_matrix[i, j] = 1

        corner = match_matrix[i-1, j-1]
        top = match_matrix[i-1, j]
        left = match_matrix[i, j-1]

        letter_a = sequence_a[i-1]
        letter_b = sequence_b[j-1]

        if curr - match == corner and sequence_a[i-1] == sequence_b[j-1]:
            alignment_a.append(letter_a)
            alignment_b.append(letter_b)
            i -= 1
            j -= 1
        elif curr - mismatch == corner and sequence_a[i-1] != sequence_b[j-1]:
            alignment_a.append(letter_a)
            alignment_b.append(letter_b)
            i -= 1
            j -= 1
        elif curr - gap == top:
            alignment_a.append(letter_a)
            alignment_b.append("-")
            i -= 1
        elif curr - gap == left:
            alignment_a.append("-")
            alignment_b.append(letter_b)
            j -= 1
        else:
            break
    color_matrix[i, j] = 1
    alignments.append([alignment_a[::-1], alignment_b[::-1]])

    results = {
        "matrix": match_matrix,
        "color": color_matrix,
        "score": score,
        "alignments": alignments
    }

    return results

# Multiple sequence alignment
def multiple_sequence_alignment(path:str=INPUT_LOC):
    """Generate multiple sequence alignment using Clustal Omega

    Args:
        path (str, optional): path of desired .fasta file. Defaults to INPUT_LOC.

    Returns:
        dict: dictionary of sequences
    """
    try:
        os.remove(OUTPUT_LOC)  # Delete file before sequencing
    except Exception as _:
        print("No file to delete")

    clustal_omega_cline = ClustalOmegaCommandline(
        infile=path, outfile=OUTPUT_LOC, verbose=True, auto=True)
    clustal_omega_cline()

    sequences = dict()
    sequences_alignment = SeqIO.to_dict(SeqIO.parse(OUTPUT_LOC, "fasta"))

    for key, value in sequences_alignment.items():
        sequences[key] = value.seq._data
    return sequences

#####################################################################
# Drawing plots
#####################################################################

# Draw the matrix of alignment
def draw_match_matrix(fig, axis, sequence_a: str, sequence_b: str, match_matrix: np.ndarray, color_matrix: np.ndarray):
    """Draw the matrix of alignment
    """
    # Draw the map
    axis.clear()

    rows, cols = len(sequence_a) + 1, len(sequence_b) + 1
    table = pd.DataFrame(match_matrix, columns=list(" " + sequence_b))
    colors = []
    for j in range(0, rows):
        color_row = []
        for i in range(0, cols):
            if color_matrix[j, i] == 0:
                color_row.append("#171717")
            else:
                color_row.append("#00a1c9")

        colors.append(color_row)

    colors = np.array(colors)

    row_label = " " + sequence_a
    # Color the map
    axis.table(cellText=table.values,
               cellColours=colors,
               cellLoc='center',
               colWidths=[0.05 for x in table.columns],
               rowLabels=row_label,
               rowColours=["#a1a1a1" for _ in row_label],
               rowLoc='center',
               colLabels=table.columns,
               colColours=["#a1a1a1" for _ in table.columns],
               colLoc='center',
               loc='center')

    fig.patch.set_visible(False)
    axis.axis('off')
    axis.axis('tight')

# Draw the multiple sequence alignment
def draw_multiple_sequence_alignment(fig, axis, sequences_dict: dict, max_size:int=100):
    """Draw the multiple sequence alignment

    Args:
        fig (matplotlib figure): figure of the plot
        axis (matplotlib axis): axis of the plot
        sequences_dict (dict): dictionary of sequences
    """
    
    # Draw the map
    axis.clear()
      
    sequences = []
    for _, sequence in sequences_dict.items():
        sequences.append(list(sequence[:min(len(sequence),max_size)].decode()))
   
    replace_dict = {'A': 0, 'G': 1, 'T': 2, 'C': 3, '-': 4}
    sequences = replace_repeated_elements(sequences, replace_dict)
    print(sequences)
    
    # Define colormap fpr every number (label2rgb)
    label2rgb_cmap = ListedColormap(
        ['#FFFFFF', '#000083', '#80FF80', '#830000', '#0080FF'],
        N=5
    )

    # Draw matrix
    axis.imshow(sequences, cmap=label2rgb_cmap)
    fig.savefig("multiple_sequence_alignment.png")

#####################################################################
# Metrics for alignment
#####################################################################

# Calculating Percent identity
def percent_identity(sequences: list):
    identical_pairs = 0

    for i in range(len(sequences[0])):
        if sequences[0][i] == sequences[1][i]:
            identical_pairs += 1

    # Calculate the total number of pairs in the multiple sequence alignment
    total_pairs = len(sequences[0])

    # Calculate the percent identity
    percent_identity = identical_pairs / total_pairs * 100

    return percent_identity

# Calculating Mutual information
def mutual_information(sequences: list, normalized=False):
    sequences += sequences
    normalized += normalized
    # residue_freq = {}

    # # Iterate over the MSA and count the number of times each residue appears
    # for sequence in sequences:
    #     for residue in sequence:
    #         if residue in residue_freq:
    #             residue_freq[residue] += 1
    #         else:
    #             residue_freq[residue] = 1
    # total_residues = sum(residue_freq.values())

    # # Initialize a variable to store the mutual information
    # mutual_information = 0

    # # Iterate over the MSA and calculate the mutual information for each pair of aligned residues
    # for i in range(len(sequences[0])):
    #     residue_1 = sequences[0][i]
    #     residue_2 = sequences[1][i]
    #     p_residue_1_residue_2 = residue_freq[residue_1+residue_2] / total_residues
    #     p_residue_1 = residue_freq[residue_1] / total_residues
    #     p_residue_2 = residue_freq[residue_2] / total_residues
    #     mutual_information += log(p_residue_1_residue_2 / (p_residue_1 * p_residue_2))

    # # Normalize the mutual information by the number of residues in the MSA
    # if normalized:
    #     mutual_information /= log(total_residues)

    # return mutual_information
    return 0

# Calculating Sum of pairs
def sum_of_pairs(sequences: list):
    """
    Calculate the sum of pairs score for a multiple sequence alignment.
    """
    
    sequences += sequences
    # # Initialize a variable to store the sum of pairs score
    # sum_of_pairs = 0

    # # Define the pair scores
    # pair_scores = {
    #     "AA":  5, "AC": -1, "AG": -2, "AT": -1,
    #     "CA": -1, "CC": 5, "CG": -3, "CT": -2,
    #     "GA": -2, "GC": -3, "GG": 5, "GT": -2,
    #     "TA": -1, "TC": -2, "TG": -2, "TT": 5
    # }

    # # Iterate over the MSA and calculate the sum of pairs score for each pair of aligned residues
    # for i in range(len(sequences[0])):
    #     residue_1 = sequences[0][i]
    #     residue_2 = sequences[1][i]
    #     sum_of_pairs += pair_scores[residue_1+residue_2]

    # return sum_of_pairs
    return 0

#####################################################################
# Other functions
#####################################################################

# Set colors of letters in multiple sequence
def color_of_letter(letter):
    if letter == "A":
        return "red"
    elif letter == "C":
        return "blue"
    elif letter == "G":
        return "green"
    elif letter == "T":
        return "yellow"
    else:
        return "black"

# Fill word with spaces
def fill_word(word, length):
    while len(word) < length:
        word += " "
    return word

# Transfer list of lists to list of strings
def lists_to_strings(list_of_lists: list):
    list_of_strings = []
    for single_list in list_of_lists:
        new_string = "".join(single_list)
        list_of_strings.append(new_string)
    return list_of_strings

# Replace repeated element in 2d list
def replace_repeated_elements(list_of_lists: list, replace_dict: dict):
    """Replace repeated elements by dictionary

    Args:
        list_of_lists (list): 2d list
        replace_dict (dict): Map of replacing elements

    Returns:
        list: modified list
    """
    
    new_list = []
    for single_list in list_of_lists:
        new_list.append([replace_dict[element] for element in single_list])
    return new_list