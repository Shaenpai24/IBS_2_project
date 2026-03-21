import numpy as np
from Bio import SeqIO

amino_acids = "ACDEFGHIKLMNPQRSTVWY"
aa_dict = {aa: i for i, aa in enumerate(amino_acids)}

max_len = 260  # pad to uniform length

X = []

for record in SeqIO.parse("prnp_proteins_clean.fasta", "fasta"):
    seq = str(record.seq)

    # Pad sequence
    seq = seq.ljust(max_len, "X")

    one_hot = np.zeros((max_len, 20))

    for i, aa in enumerate(seq):
        if aa in aa_dict:
            one_hot[i, aa_dict[aa]] = 1

    X.append(one_hot.flatten())

X = np.array(X)

print("Shape of dataset:", X.shape)
