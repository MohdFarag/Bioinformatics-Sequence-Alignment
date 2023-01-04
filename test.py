#!/usr/bin/env python
import numpy as np
import pandas as pd
import os
import matplotlib
import matplotlib.pyplot as plt
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
import textwrap
from Bio import SeqIO




path = "./Data/Group_4.fasta"
alignment = multiple_sequence_alignment(path)
print(alignment)