import pandas as pd

df = pd.read_csv("dataset_v1_core.csv")

confidence = []

for col in df.columns:
    if col in ["delta_volume", "blosum62", "delta_hydropathy"]:
        confidence.append("high")
    elif "delta" in col:
        confidence.append("medium")
    else:
        confidence.append("low")

meta = pd.DataFrame({
    "feature": df.columns,
    "confidence": confidence
})

meta.to_csv("feature_confidence.csv", index=False)