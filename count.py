count = 0

with open("sequence.txt", "r") as f:
    for line in f:
        if line.startswith(">"):
            count += 1

print("Total sequences:", count)
