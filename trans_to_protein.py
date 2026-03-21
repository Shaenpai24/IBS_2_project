from Bio import SeqIO

count = 0

with open("prnp_proteins_clean.fasta", "w") as out:
    for record in SeqIO.parse("prnp_clean_nr.fasta", "fasta"):
        seq = record.seq

        # Keep only sequences divisible by 3
        if len(seq) % 3 == 0:
            protein_seq = seq.translate(to_stop=True)

            # Filter realistic PRNP length
            if 240 <= len(protein_seq) <= 260:
                record.seq = protein_seq
                SeqIO.write(record, out, "fasta")
                count += 1

print("Clean protein sequences:", count)
