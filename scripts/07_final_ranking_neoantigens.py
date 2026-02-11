import pandas as pd

# Load Step 9 output
df = pd.read_csv("results/kras_G12D_immunogenic_neoantigens.csv")

# -----------------------------
# 1. Ensure numeric columns
# -----------------------------
score_cols = ["norm_ic50", "norm_TPM", "norm_vaf", "immunogenicity_score"]

for col in score_cols:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# -----------------------------
# 2. Fix TPM normalization
#    (all TPM values identical)
# -----------------------------
if df["TPM"].nunique() == 1:
    df["norm_TPM"] = 1
else:
    df["norm_TPM"] = (
        (df["TPM"] - df["TPM"].min()) /
        (df["TPM"].max() - df["TPM"].min())
    )

# -----------------------------
# 3. Handle any remaining NaNs
# -----------------------------
df[score_cols] = df[score_cols].fillna(0)

# -----------------------------
# 4. Final composite score
# -----------------------------
df["final_score"] = (
    0.35 * df["norm_ic50"] +
    0.25 * df["norm_TPM"] +
    0.20 * df["norm_vaf"] +
    0.20 * df["immunogenicity_score"]
)

# -----------------------------
# 5. Rank neoantigens
# -----------------------------
df = df.sort_values("final_score", ascending=False)

# -----------------------------
# 6. Save final output
# -----------------------------
df.to_csv(
    "results/kras_G12D_final_ranked_neoantigens.csv",
    index=False
)

# -----------------------------
# 7. Quick verification
# -----------------------------
print(df[["peptide", "final_score"]])

