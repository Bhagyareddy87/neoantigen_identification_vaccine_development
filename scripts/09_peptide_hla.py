import pandas as pd

# Load your file
data = pd.read_csv("results/kras_G12D_peptide_hla.csv")

# Convert long format → wide format
heatmap_data = data.pivot(index="peptide",
                          columns="HLA_allele",
                          values="IC50")

# Save to check
heatmap_data.to_csv("results/kras_G12D_heatmap_matrix.csv")

print(heatmap_data.head())

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# ==================================================
# 0. Ensure results folder exists
# ==================================================

os.makedirs("results", exist_ok=True)

# ==================================================
# 1. Load original file
# ==================================================

data = pd.read_csv("results/kras_G12D_peptide_hla.csv")
data.columns = data.columns.str.strip()

# ==================================================
# 2. Convert long → wide format
# ==================================================

heatmap_data = data.pivot_table(
    index="peptide",
    columns="HLA_allele",
    values="IC50",
    aggfunc="min"
)

heatmap_data = heatmap_data.fillna(100000)

# ==================================================
# 3. Select top 25 strongest peptides
# ==================================================

top_peptides = heatmap_data.min(axis=1).sort_values().head(25).index
heatmap_top = heatmap_data.loc[top_peptides]

# ==================================================
# 4. Log transform
# ==================================================

heatmap_log = np.log10(heatmap_top)

# ==================================================
# 5. Heatmap (Main Figure)
# ==================================================

sns.set_style("white")

plt.figure(figsize=(10, 14))
sns.heatmap(
    heatmap_log,
    cmap="viridis",
    linewidths=0.4,
    linecolor='white',
    cbar_kws={'label': 'log10(IC50 nM)'}
)

plt.xticks(rotation=45, ha='right')
plt.title("KRAS G12D Peptide–HLA Binding Landscape", fontsize=14)

plt.tight_layout()
plt.savefig("results/KRAS_G12D_Heatmap.png", dpi=300, bbox_inches='tight')
plt.close()

# ==================================================
# 6. Strong Binder Count per HLA (Moderate cutoff)
# ==================================================

strong_threshold = 2.7  # log10(500 nM)

strong_counts = (heatmap_log < strong_threshold).sum()

plt.figure(figsize=(8,6))
strong_counts.sort_values(ascending=False).plot(kind='bar')

plt.ylabel("Number of Peptides (IC50 < 500 nM)")
plt.xlabel("HLA Allele")
plt.title("Moderate/Strong Binder Distribution per HLA")
plt.xticks(rotation=45)

plt.tight_layout()
plt.savefig("results/strong_binder_distribution.png", dpi=300)
plt.close()

# ==================================================
# 7. Global Binding Distribution
# ==================================================

plt.figure(figsize=(8,6))
sns.histplot(heatmap_log.values.flatten(), bins=30, kde=True)

plt.xlabel("log10(IC50 nM)")
plt.title("Global Distribution of Peptide-HLA Binding")

plt.tight_layout()
plt.savefig("results/global_binding_distribution.png", dpi=300)
plt.close()

# ==================================================
# 8. Top 10 Vaccine Candidate Peptides
# ==================================================

best_binding = heatmap_log.min(axis=1).sort_values().head(10)

plt.figure(figsize=(8,6))
best_binding.plot(kind='barh')

plt.xlabel("Best log10(IC50)")
plt.title("Top 10 KRAS G12D Vaccine Candidates")

plt.tight_layout()
plt.savefig("results/top10_candidates.png", dpi=300)
plt.close()

# ==================================================
# 9. Population Coverage Style Plot
# ==================================================

coverage = (heatmap_log < strong_threshold).sum(axis=1)

plt.figure(figsize=(8,6))
coverage.sort_values(ascending=False).head(10).plot(kind='bar')

plt.ylabel("Number of HLAs Covered (IC50 < 500 nM)")
plt.title("Top Peptides by HLA Coverage")

plt.tight_layout()
plt.savefig("results/population_coverage.png", dpi=300)
plt.close()

print("✅ All plots generated successfully inside 'results/' folder.")


