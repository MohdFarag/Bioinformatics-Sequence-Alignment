import matplotlib.pyplot as plt
import numpy as np

def draw_multiple_sequence_alignment(fig, ax, sequences:dict):
    """
    Draw plot for the multiple sequence alignment
    """
    ax.set_title("Multiple Sequence Alignment")
    ax.set_xlabel("Sequence")
    ax.set_ylabel("Position")
    ax.set_xticks(np.arange(len(sequences)))
    ax.set_xticklabels([x for x in sequences.keys()])
    ax.set_yticks(np.arange(len(sequences[list(sequences.keys())[0]])))
    ax.set_yticklabels([x for x in sequences[list(sequences.keys())[0]]])
    ax.imshow(np.array([list(x) for x in sequences.values()]), cmap="Greys")

fig, ax = plt.subplots()
sequences = {
    1:"AGCAGAGACA",
    2:"TACTGTATCA"
}
draw_multiple_sequence_alignment(fig, ax, sequences)