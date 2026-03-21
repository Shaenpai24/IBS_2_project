from Bio import SeqIO

lengths = []

for record in SeqIO.parse("sequence.txt", "fasta"):
    lengths.append(len(record.seq))

print("Min length:", min(lengths))
print("Max length:", max(lengths))
print("Average length:", sum(lengths)/len(lengths))
