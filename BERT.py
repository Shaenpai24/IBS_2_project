import torch
import pandas as pd
import numpy as np
from transformers import AutoTokenizer, AutoModel

# Load your dataset
df = pd.read_csv("prnp_feature_dataset_extended.csv")

sequences = df["sequence"].tolist()

# Load ESM-2 model
model_name = "facebook/esm2_t12_35M_UR50D"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModel.from_pretrained(model_name)

model.eval()

embeddings = []

with torch.no_grad():
    for seq in sequences:
        inputs = tokenizer(seq, return_tensors="pt", truncation=True)
        outputs = model(**inputs)

        # Mean pooling (average over tokens)
        last_hidden = outputs.last_hidden_state
        mean_embedding = last_hidden.mean(dim=1).squeeze().numpy()

        embeddings.append(mean_embedding)

embeddings = np.array(embeddings)

print("Embedding shape:", embeddings.shape)

np.save("prnp_esm_embeddings.npy", embeddings)
bert_df = pd.DataFrame(
    embeddings,
    columns=[f"bert_{i}" for i in range(embeddings.shape[1])]
)

df_final = pd.concat([df, bert_df], axis=1)

df_final.to_csv("prnp_full_dataset_with_bert.csv", index=False)

print("Final dataset shape:", df_final.shape)
