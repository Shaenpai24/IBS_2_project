import pandas as pd

features = [
    ["delta_volume", "float", "difference", "computed", "Change in amino acid volume"],
    ["delta_hydropathy", "float", "difference", "computed", "Change in hydrophobicity"],
    ["blosum62", "int", "score", "BLOSUM62", "Evolutionary substitution score"],
    ["delta_charge", "int", "difference", "computed", "Charge change"],
    ["position", "int", "index", "sequence", "Mutation position"],
    ["conservation_wt_freq", "float", "frequency", "MSA", "Wild-type conservation"],
]

df_reg = pd.DataFrame(features, columns=["feature", "type", "unit", "source", "description"])
df_reg.to_csv("feature_registry.csv", index=False)