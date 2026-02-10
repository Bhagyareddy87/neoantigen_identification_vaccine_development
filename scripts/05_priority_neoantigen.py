import pandas as pd
import numpy as np

# Load Step-7 output
df = pd.read_csv("results/kras_G12D_final_neoantigens.csv")

# ---- Step 8: Neoantigen Prioritization ----

# 1. Normalize IC50 (lower is better → invert)
df["norm_ic50"] = 1 - (
    np.log10(df["IC50_nM"] + 1) - np.log10(df["IC50_nM"].min() + 1)
) / (
    np.log10(df["IC50_nM"].max() + 1) - np.log10(df["IC50_nM"].min() + 1)
)

# 2. Normalize TPM
df["norm_TPM"] = (
    df["TPM"] - df["TPM"].min()
) / (
    df["TPM"].max() - df["TPM"].min()
)

# 3. Normalize VAF
df["norm_vaf"] = (
    df["vaf"] - df["vaf"].min()
) / (
    df["vaf"].max() - df["vaf"].min()
)

# 4. Composite neoantigen priority score
df["neoantigen_score"] = (
    0.5 * df["norm_ic50"] +
    0.3 * df["norm_vaf"] +
    0.2 * df["norm_TPM"]
)

# 5. Rank candidates
df = df.sort_values("neoantigen_score", ascending=False)

# Save Step-8 output
df.to_csv("results/kras_G12D_priority_neoantigens.csv", index=False)

print("Step 8 completed: Neoantigens ranked and saved!")
