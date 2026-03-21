from Bio import SeqIO

input_file = "sequence.txt"
output_file = "prnp_clean_cds.fasta"

count = 0

with open(output_file, "w") as out:
    for record in SeqIO.parse(input_file, "fasta"):
        if 700 <= len(record.seq) <= 900:
            SeqIO.write(record, out, "fasta")
            count += 1

print("Filtered sequences:", count)
