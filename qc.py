import pandas as pd

# Load dataset
df = pd.read_csv("prnp_mutation_features_tier2.csv")

print("Initial shape:", df.shape)

# -----------------------------
# 1. Duplicate Check
# -----------------------------
dup_count = df.duplicated(subset=["mutation"]).sum()
print("Duplicate mutations:", dup_count)

df = df.drop_duplicates(subset=["mutation"])

# -----------------------------
# 2. Missing Values Check
# -----------------------------
missing = df.isnull().sum()
missing = missing[missing > 0]

print("\nColumns with missing values:")
print(missing if len(missing) > 0 else "No missing values ")

# Optional: drop rows with missing values
df = df.dropna()

# -----------------------------
# 3. Valid Amino Acid Check
# -----------------------------
valid_aa = set("ACDEFGHIKLMNPQRSTVWY")

invalid_rows = df[
    (~df["wt_aa"].isin(valid_aa)) |
    (~df["mut_aa"].isin(valid_aa))
]

print("\nInvalid amino acid rows:", len(invalid_rows))

# -----------------------------
# 4. Mutation Consistency Check
# -----------------------------
def check_mutation(row):
    mutation = row["mutation"]   # e.g. M1A
    wt = row["wt_aa"]
    mut = row["mut_aa"]
    pos = str(row["position"])

    try:
        return mutation == f"{wt}{pos}{mut}"
    except:
        return False

invalid_mutations = df[~df.apply(check_mutation, axis=1)]

print("Invalid mutation format rows:", len(invalid_mutations))

# -----------------------------
# 5. Position Consistency Check
# -----------------------------
invalid_positions = df[df["position"] <= 0]

print("Invalid positions:", len(invalid_positions))

# -----------------------------
# 6. Final Clean Dataset
# -----------------------------
df_clean = df.copy()

print("\nFinal cleaned shape:", df_clean.shape)

# Save cleaned dataset
df_clean.to_csv("prnp_mutation_features_tier2_clean.csv", index=False)

print("Saved cleaned dataset ")