import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# -----------------------------
# Helper Feature Functions
# -----------------------------

def longest_qn_run(seq):
    max_run = 0
    current = 0
    for aa in seq:
        if aa in ["Q", "N"]:
            current += 1
            max_run = max(max_run, current)
        else:
            current = 0
    return max_run


def hydrophobic_stretch(seq, window=5):
    hydrophobic = set("AILMFWYV")
    count = 0
    for i in range(len(seq) - window + 1):
        segment = seq[i:i+window]
        if all(aa in hydrophobic for aa in segment):
            count += 1
    return count


def aliphatic_index(seq):
    A = seq.count("A") / len(seq)
    V = seq.count("V") / len(seq)
    I = seq.count("I") / len(seq)
    L = seq.count("L") / len(seq)
    return (A * 100) + (2.9 * V * 100) + (3.9 * (I + L) * 100)


def net_charge(seq):
    pos = seq.count("K") + seq.count("R") + seq.count("H")
    neg = seq.count("D") + seq.count("E")
    return pos - neg


def motif_count(seq, motif="PHGGGWGQ"):
    return seq.count(motif)


def clean_sequence(seq):
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    return "".join([aa for aa in seq if aa in valid_aa])


def compute_features(seq):
    seq = clean_sequence(seq)

    if len(seq) == 0:
        return None

    analysed = ProteinAnalysis(seq)
    aa_percent = analysed.amino_acids_percent  # new method (no deprecation)

    features = {}

    features["length"] = len(seq)
    features["molecular_weight"] = analysed.molecular_weight()
    features["aromaticity"] = analysed.aromaticity()
    features["instability_index"] = analysed.instability_index()
    features["gravy"] = analysed.gravy()
    features["pI"] = analysed.isoelectric_point()
    features["aliphatic_index"] = aliphatic_index(seq)
    features["net_charge"] = net_charge(seq)

    # Q/N relevant
    features["Q_percent"] = seq.count("Q") / len(seq)
    features["N_percent"] = seq.count("N") / len(seq)
    features["QN_longest_run"] = longest_qn_run(seq)

    # Aggregation proxies
    features["hydrophobic_stretch_count"] = hydrophobic_stretch(seq)
    features["aromatic_density"] = (seq.count("F")+seq.count("W")+seq.count("Y")) / len(seq)
    features["glycine_density"] = seq.count("G") / len(seq)
    features["cysteine_count"] = seq.count("C")
    features["proline_count"] = seq.count("P")
    features["octapeptide_repeat_count"] = motif_count(seq)

    # Composition categories
    polar = set("STNQ")
    nonpolar = set("AVLIMFWY")
    tiny = set("AGST")
    bulky = set("FWY")

    features["polar_percent"] = sum(seq.count(x) for x in polar) / len(seq)
    features["nonpolar_percent"] = sum(seq.count(x) for x in nonpolar) / len(seq)
    features["tiny_percent"] = sum(seq.count(x) for x in tiny) / len(seq)
    features["bulky_percent"] = sum(seq.count(x) for x in bulky) / len(seq)

    # Disorder proxies
    disorder_promoting = set("GPQSENK")
    order_promoting = set("WFYILV")

    features["disorder_promoting_fraction"] = sum(seq.count(x) for x in disorder_promoting) / len(seq)
    features["order_promoting_fraction"] = sum(seq.count(x) for x in order_promoting) / len(seq)

    # 20 AA composition
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        features[f"AA_{aa}"] = aa_percent.get(aa, 0)

    # Codon 129
    if len(seq) > 128:
        residue_129 = seq[128]
    else:
        residue_129 = None

    features["residue_129"] = residue_129
    features["is_M129"] = 1 if residue_129 == "M" else 0
    features["is_V129"] = 1 if residue_129 == "V" else 0

    return features


# -----------------------------
# Load FASTA & Build Dataset
# -----------------------------

data = []

for record in SeqIO.parse("prnp_proteins_clean.fasta", "fasta"):
    seq = str(record.seq)
    row = {"sequence_id": record.id, "sequence": seq}
    row.update(compute_features(seq))
    data.append(row)

df = pd.DataFrame(data)

df.to_csv("prnp_feature_dataset_extended.csv", index=False)

print("Final dataset shape:", df.shape)
print("Saved as prnp_feature_dataset_extended.csv")
