import pandas as pd

# Load Step 8 output
df = pd.read_csv("results/kras_G12D_priority_neoantigens.csv")

# Load Step 9 immunogenicity scores
imm = pd.read_csv("results/immunogenicity_scores.csv")

# Rename column if needed
imm = imm.rename(columns={"score": "immunogenicity_score"})

# Merge on peptide
df = df.merge(
    imm[["peptide", "immunogenicity_score"]],
    on="peptide",
    how="left"
)

# Save Step 9 output
df.to_csv(
    "results/kras_G12D_immunogenic_neoantigens.csv",
    index=False
)

print("Step 9 completed correctly")


