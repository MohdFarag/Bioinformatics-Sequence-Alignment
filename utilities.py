#!/usr/bin/env python
import numpy as np
import os

def read_fasta_sequences(path: str):
    try:
        with open(path, "r", encoding='utf-8-sig') as file:
            # Read the text from a file
            sequences = file.read()
            print(sequences)
    except Exception as e:
        if "No such file or directory" in e.__str__():
            raise Exception("**ERROR** No such file or directory")

    sequences_group = {}
    for line in sequences:
        if(line.startswith(">")):
            key = line

def global_alignment(seq1:str, seq2:str, match:int=1, mismatch:int=1, gap:int=1):
    """
    The Needleman-Wunsch Algorithm
    ==============================
    This is a dynamic programming algorithm for finding the optimal alignment of
    two strings.
    
    Arguments
    -------
    seq1: first sequence (Must be combination of A,C,G,T)
    seq2: seconde sequence (Must be combination of A,C,G,T)
    match: matching score
    mismatch: mismatching score
    gap: gapping score

    Returns
    -------
    rx: first sequence alignment
    ry: second sequence alignment

    Example
    -------
        >>> x = "GATTACA"
        >>> y = "GCATGCU"
        >>> global_alignment(x, y)
        G-ATTACA
        GCA-TGCU
    """

    if type(seq1) != str:
        seq1 = "".join(seq1)
    
    if type(seq2) != str:
        seq2 = "".join(seq2)

    if type(match) == str:
        match = int(match)

    if type(mismatch) == str:
        mismatch = int(mismatch)

    if type(gap) == str:
        gap = int(gap)

    nx = len(seq1)
    ny = len(seq2)

    # Optimal score at each possible pair of characters.
    alignment_sequence = np.zeros((nx + 1, ny + 1))
    alignment_sequence[:,0] = np.linspace(0, -nx * gap, nx + 1)
    alignment_sequence[0,:] = np.linspace(0, -ny * gap, ny + 1)
    # Pointers to trace through an optimal alignment.
    P = np.zeros((nx + 1, ny + 1))
    P[:,0] = 3
    P[0,:] = 4
    
    # Temporary scores.
    t = np.zeros(3)
    for i in range(nx):
        for j in range(ny):
            if seq1[i] == seq2[j]:
                t[0] = alignment_sequence[i,j] + match
            else:
                t[0] = alignment_sequence[i,j] - mismatch
            t[1] = alignment_sequence[i,j+1] - gap
            t[2] = alignment_sequence[i+1,j] - gap
            tmax = np.max(t)
            alignment_sequence[i+1,j+1] = tmax
            if t[0] == tmax:
                P[i+1,j+1] += 2
            if t[1] == tmax:
                P[i+1,j+1] += 3
            if t[2] == tmax:
                P[i+1,j+1] += 4

    # Trace through an optimal alignment.
    i = nx
    j = ny
    rx = []
    ry = []
    while i > 0 or j > 0:
        if P[i,j] in [2, 5, 6, 9]:
            rx.append(seq1[i-1])
            ry.append(seq2[j-1])
            i -= 1
            j -= 1
        elif P[i,j] in [3, 5, 7, 9]:
            rx.append(seq1[i-1])
            ry.append('-')
            i -= 1
        elif P[i,j] in [4, 6, 7, 9]:
            rx.append('-')
            ry.append(seq2[j-1])
            j -= 1
    
    # Reverse the strings.
    rx = rx[::-1]
    ry = ry[::-1]
    return rx, ry

def local_alignment(seq1:str, seq2:str, match:int=1, mismatch:int=1, gap:int=1):
    """
    The Smith-Waterman Algorithm
    ==============================
    This is a dynamic programming algorithm for finding the optimal alignment of
    two strings.
    
    Arguments
    -------
    seq1: first sequence (Must be combination of A,C,G,T)
    seq2: second sequence (Must be combination of A,C,G,T)
    match: matching score
    mismatch: mismatching score
    gap: gapping score

    Returns
    -------
    rx: first sequence alignment
    ry: second sequence alignment

    Example
    -------
        >>> x = "GATTACA"
        >>> y = "GCATGCU"
        >>> global_alignment(x, y)
        G-ATTACA
        GCA-TGCU
    """

    if type(seq1) != str:
        seq1 = "".join(seq1)
    
    if type(seq2) != str:
        seq2 = "".join(seq2)

    nx = len(seq1)
    ny = len(seq2)