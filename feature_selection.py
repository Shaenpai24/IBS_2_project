import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier

# -----------------------------
# Load cleaned dataset
# -----------------------------
df = pd.read_csv("prnp_mutation_features_tier2_clean.csv")

print("Initial shape:", df.shape)

# -----------------------------
# Keep only numeric columns
# -----------------------------
df_model = df.select_dtypes(include=[np.number]).copy()

print("Numeric feature shape:", df_model.shape)

# -----------------------------
# Create proxy labels
# -----------------------------
df_model["label"] = (df["grantham"] > 100).astype(int)

# Separate features and labels
X = df_model.drop(columns=["label", "grantham"])
y = df_model["label"]

# -----------------------------
# 1. Correlation Filtering
# -----------------------------
corr_matrix = X.corr().abs()

upper = corr_matrix.where(
    np.triu(np.ones(corr_matrix.shape), k=1).astype(bool)
)

to_drop = [col for col in upper.columns if any(upper[col] > 0.9)]

print("\nHighly correlated features to drop:", len(to_drop))

X_reduced = X.drop(columns=to_drop)

print("Shape after correlation filter:", X_reduced.shape)

# -----------------------------
# 2. Train model for importance
# -----------------------------
model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_reduced, y)

importances = pd.Series(model.feature_importances_, index=X_reduced.columns)
importances = importances.sort_values(ascending=False)

print("\nTop 20 important features:\n")
print(importances.head(20))

# -----------------------------
# 3. Select top features (CORE)
# -----------------------------
top_features = importances.head(20).index.tolist()

core_df = df[["mutation"] + top_features]

# -----------------------------
# 4. Extended dataset
# -----------------------------
extended_df = df.copy()

# -----------------------------
# Save outputs
# -----------------------------
core_df.to_csv("dataset_v1_core.csv", index=False)
extended_df.to_csv("dataset_v1_extended.csv", index=False)

print("\nSaved:")
print("dataset_v1_core.csv")
print("dataset_v1_extended.csv")