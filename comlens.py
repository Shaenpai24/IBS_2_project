from Bio import SeqIO

lengths = []

for record in SeqIO.parse("prnp_proteins_clean.fasta", "fasta"):
    lengths.append(len(record.seq))

print(set(lengths))
