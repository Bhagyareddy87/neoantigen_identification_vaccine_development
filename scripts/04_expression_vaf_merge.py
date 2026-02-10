import pandas as pd

# 1. Load files
binders = pd.read_csv("results/kras_G12D_strong_binders.csv")
vaf = pd.read_csv("data/processed/vaf_simulated.csv")
rna = pd.read_csv("data/processed/rna_final.csv")

# 🔧 Clean column names
binders.columns = binders.columns.str.strip()
vaf.columns = vaf.columns.str.strip()
rna.columns = rna.columns.str.strip()

# 2. Merge binders + VAF on peptide + allele
merged = binders.merge(
    vaf,
    on=["peptide", "allele"],
    how="left"
)

# 3. Get KRAS RNA expression
kras_tpm = rna.loc[rna["gene"] == "KRAS", "TPM"].values[0]

merged["TPM"] = kras_tpm
merged["gene"] = "KRAS"

# 4. Apply biological filters
final_candidates = merged[
    (merged["IC50_nM"] < 500) &
    (merged["TPM"] >= 1) &
    (merged["vaf"] >= 0.05)
]

# 5. Save result
final_candidates.to_csv(
    "results/kras_G12D_final_neoantigens.csv",
    index=False
)

print("Final neoantigen candidates saved!")


print(binders.columns.tolist())


import pandas as pd

binders = pd.read_csv("results/kras_G12D_strong_binders.csv")
binders.columns = binders.columns.str.strip()

print(binders.columns.tolist())
