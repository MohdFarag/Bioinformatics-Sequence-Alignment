#!/usr/bin/env python
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

# Reference of sequences letters: 
# http://web.mit.edu/meme_v4.11.4/share/doc/alphabets.html
LETTERS_OF_DNA = ['A','G','T','C','N','X']
LETTERS_OF_RNA = ['A','G','U','C','N','X']
LETTERS_OF_PROTEIN = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X', 'B', 'Z', 'J']

# Read fasta file
def read_fasta(path: str, concat=False):
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
            
            print(sequence)

            sequences_list[lines[0]] = sequence

    return sequences_list

def check_sequence(sequence:str, type:str):
    sequence = sequence.upper() 
    type = type.upper()
    if type == "DNA":
        return is_dna(sequence)
    elif type == "RNA":
        return is_rna(sequence)
    elif type == "PROTEIN":
        return is_protein(sequence)

    return False

def is_dna(sequence):
    # Valid letters in DNA sequence
    for letter in sequence:
        if letter not in LETTERS_OF_DNA:
            return False
    return True

def is_rna(sequence):
    # Valid letters in RNA sequence
    for letter in sequence:
        if letter not in LETTERS_OF_RNA:
            return False
    return True

def is_protein(sequence):
    # Valid letters in protein sequence
    for letter in sequence:
        if letter not in LETTERS_OF_PROTEIN:
            return False
    return True

def pairwise_global_alignment(sequence_a:str, sequence_b:str, match:int=1, mismatch:int=0, gap:int=-1):
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

    if type(sequence_a) != str:
        sequence_a = "".join(sequence_a)
    
    if type(sequence_b) != str:
        sequence_b = "".join(sequence_b)

    if type(match) == str:
        match = int(match)

    if type(mismatch) == str:
        mismatch = int(mismatch)

    if type(gap) == str:
        gap = int(gap)

    # STEP 1: Initialization of matrix
    match_matrix = np.zeros((len(sequence_a)+1, len(sequence_b)+1))
    color_matrix = np.zeros((len(sequence_a)+1, len(sequence_b)+1))
    
    for i in range (1,len(sequence_a)+1):
        match_matrix[i,0] = match_matrix[i-1,0] + gap

    for j in range (1,len(sequence_b)+1):
        match_matrix[0,j] = match_matrix[0,j-1] + gap

    # STEP 2: Filling to the matrix
    for i in range(1,len(sequence_a)+1):
        for j in range(1,len(sequence_b)+1):
            arr = []
            if sequence_a[i-1] == sequence_b[j-1]:
                arr.append(match_matrix[i-1,j-1] + match)
            elif sequence_a[i-1] != sequence_b[j-1]:
                arr.append(match_matrix[i-1,j-1] + mismatch)
            
            arr.append(match_matrix[i-1,j] + gap)
            arr.append(match_matrix[i,j-1] + gap)
            match_matrix[i,j] = max(arr)

    score = match_matrix[-1,-1]
    print(score)
    print(match_matrix)
    # STEP 3: Backtracing
    alignments = []
    i, j = match_matrix.shape[0] - 1, match_matrix.shape[1] - 1

    alignment_a = []
    alignment_b = []
    while not(i == 0 and j == 0):
        letter_a = sequence_a[i-1]
        letter_b = sequence_b[j-1]
        curr = match_matrix[i,j]
        color_matrix[i,j] = 1

        corner = match_matrix[i-1,j-1]
        top = match_matrix[i-1,j]
        left = match_matrix[i,j-1]
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


    alignments.append([alignment_a[::-1],alignment_b[::-1]])

    results = {
        "matrix":match_matrix,
        "color": color_matrix,
        "score": score,
        "alignments": alignments
    }

    return results

def pairwise_local_alignment(sequence_a:str, sequence_b:str, match:int=1, mismatch:int=0, gap:int=-1):
    """
    The Smith-Waterman Algorithm
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
    
    if type(sequence_a) != str:
        sequence_a = "".join(sequence_a)
    
    if type(sequence_b) != str:
        sequence_b = "".join(sequence_b)

    if type(match) == str:
        match = int(match)

    if type(mismatch) == str:
        mismatch = int(mismatch)

    if type(gap) == str:
        gap = int(gap)

    # STEP 1: Initialization of matrix
    match_matrix = np.zeros((len(sequence_a)+1, len(sequence_b)+1))
    color_matrix = np.zeros((len(sequence_a)+1, len(sequence_b)+1))

    # STEP 2: Filling to the matrix
    for i in range(1,len(sequence_a)+1):
        for j in range(1,len(sequence_b)+1):
            arr = [0]
            if sequence_a[i-1] == sequence_b[j-1]:
                arr.append(match_matrix[i-1,j-1] + match)
            elif sequence_a[i-1] != sequence_b[j-1]:
                arr.append(match_matrix[i-1,j-1] + mismatch)
            
            arr.append(match_matrix[i-1,j] + gap)
            arr.append(match_matrix[i,j-1] + gap)
            match_matrix[i,j] = max(arr)
    score = match_matrix.max()

    # STEP 3: Backtracing
    alignments = []

    index_max = np.unravel_index(np.argmax(match_matrix, axis=None), match_matrix.shape)
    i, j = index_max[0], index_max[1]
    alignment_a = []
    alignment_b = []
    
    while match_matrix[i,j] != 0:
        curr = match_matrix[i,j]
        color_matrix[i,j] = 1

        corner = match_matrix[i-1,j-1]
        top = match_matrix[i-1,j]
        left = match_matrix[i,j-1]

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

    alignments.append([alignment_a[::-1],alignment_b[::-1]])

    results = {
        "matrix":match_matrix,
        "color": color_matrix,
        "score": score,
        "alignments": alignments
    }

    return results

def draw_match_matrix(sequence_a:str, sequence_b:str, match_matrix:np.ndarray, color_matrix:np.ndarray):
    N = len(sequence_a) + 1
    M = len(sequence_b) + 1
    
    # Draw the map
    fig, ax = plt.subplots()
    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')

    df = pd.DataFrame(match_matrix, columns=list(" " + sequence_a))

    color_matrix[color_matrix==0] = "w"
    color_matrix[color_matrix==1] = "g"

    # Color the map
    ax.table(cellText=df.values, colLabels=df.columns,
            loc='center', rowLabels=" " + sequence_b,
            cellColours=color_matrix, colWidths=[0.05 for x in df.columns], )

    fig.tight_layout()
    plt.show()

def multiple_sequence_alignment(sequences:list, match:int=1, mismatch:int=0, gap:int=-1):
    pass