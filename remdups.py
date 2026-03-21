from Bio import SeqIO

unique = {}
for record in SeqIO.parse("prnp_clean_cds.fasta", "fasta"):
    seq = str(record.seq)
    if seq not in unique:
        unique[seq] = record

with open("prnp_clean_nr.fasta", "w") as out:
    for record in unique.values():
        SeqIO.write(record, out, "fasta")

print("Unique sequences:", len(unique))
