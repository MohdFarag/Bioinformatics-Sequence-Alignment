from math import *

def mutual_information(sequences:list, normalized=False):
    residue_freq = {}

    # Iterate over the MSA and count the number of times each residue appears
    for seq in sequences:
        for residue in seq:
            if residue not in residue_freq:
                residue_freq[residue] = 1
            else:
                residue_freq[residue] += 1
    
    # Calculate the total number of residues
    total_residues = sum(residue_freq.values())

    # Calculate the probability of each residue
    residue_prob = {}
    for residue in residue_freq:
        residue_prob[residue] = residue_freq[residue] / total_residues

    return residue_prob

seq = ['AGCAGA-----', 'AGCAGAAAAAG']

print(mutual_information(seq, normalized=True))