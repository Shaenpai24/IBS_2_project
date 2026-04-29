import pandas as pd
import shap
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier

# -----------------------------
# Load datasets
# -----------------------------
core_df = pd.read_csv("dataset_v1_core.csv")
full_df = pd.read_csv("prnp_mutation_features_tier2_clean.csv")

# -----------------------------
# Create labels
# -----------------------------
core_df["label"] = (full_df["grantham"] > 100).astype(int)

# -----------------------------
# Prepare data
# -----------------------------
X = core_df.drop(columns=["mutation", "label"])
y = core_df["label"]

# -----------------------------
# Train model
# -----------------------------
model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X, y)

print("Model trained successfully.")

# -----------------------------
# SHAP analysis 
# -----------------------------
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X)

# Handle different SHAP output formats
if isinstance(shap_values, list):
    shap_vals = shap_values[1]   # classification → class 1
else:
    shap_vals = shap_values

# -----------------------------
#  SHAP Summary Plot (PNG)
# -----------------------------
plt.figure()
shap.summary_plot(shap_vals, X, show=False)
plt.tight_layout()
plt.savefig("shap_summary.png", dpi=300)
plt.close()

print("Saved shap_summary.png ")

# -----------------------------
#  SHAP Bar Plot
# -----------------------------
plt.figure()
shap.summary_plot(shap_vals, X, plot_type="bar", show=False)
plt.tight_layout()
plt.savefig("shap_bar.png", dpi=300)
plt.close()

print("Saved shap_bar.png ")

