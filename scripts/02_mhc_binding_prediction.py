import pandas as pd
from mhcflurry import Class1AffinityPredictor

# Load peptides
peptides = []
with open("data/processed/kras_G12D_neoantigens.txt") as f:
    for line in f:
        peptides.append(line.strip())

# HLA alleles
alleles = ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02"]

# Expand peptide × allele
peptide_list = []
allele_list = []
for pep in peptides:
    for allele in alleles:
        peptide_list.append(pep)
        allele_list.append(allele)

# Load predictor
predictor = Class1AffinityPredictor.load()

# Run prediction
results = predictor.predict(
    peptides=peptide_list,
    alleles=allele_list
)

# Convert to DataFrame
df = pd.DataFrame(results)
df["peptide"] = peptide_list
df["allele"] = allele_list

# Save output
df.to_csv(
    "results/kras_G12D_mhcflurry_affinity.csv",
    index=False
)

print("MHCflurry affinity prediction saved successfully")








import pandas as pd

# Load MHCflurry results
df = pd.read_csv("results/kras_G12D_mhcflurry_affinity.csv")

# Rename first column (clean)
df.columns = ["ic50_nm", "peptide", "allele"]

# Classify binding strength
df["binding_class"] = "non-binder"
df.loc[df["ic50_nm"] < 500, "binding_class"] = "binder"
df.loc[df["ic50_nm"] < 50, "binding_class"] = "strong_binder"

# Keep only binders
binders = df[df["ic50_nm"] < 500].copy()

# Rank by affinity (best first)
binders = binders.sort_values("ic50_nm")

# Save
binders.to_csv(
    "results/kras_G12D_strong_binders.csv",
    index=False
)

print("Filtered and ranked binders saved")
