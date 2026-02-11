import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# ===============================
# 1. Load Data
# ===============================

data = pd.read_csv("results/kras_G12D_IC50_immunogenicity.csv")

# Clean column names (good practice)
data.columns = data.columns.str.strip()

# ===============================
# 2. Define Binder Type
# ===============================

data["Binder_Type"] = data["IC50"].apply(
    lambda x: "Valid (<500 nM)" if x < 500 else "Weak (≥500 nM)"
)

# ===============================
# 3. Print Summary Statistics
# ===============================

print("\nGroup Summary Statistics:")
print(data.groupby("Binder_Type")["Immunogenicity"].describe())

# ===============================
# 4. Mann–Whitney U Test
# ===============================

valid = data[data["IC50"] < 500]["Immunogenicity"]
weak = data[data["IC50"] >= 500]["Immunogenicity"]

stat, p_value = mannwhitneyu(valid, weak, alternative="two-sided")

print("\nMann–Whitney U Test")
print("U statistic:", stat)
print("p-value:", p_value)

# ===============================
# 5. Create Publication-Style Boxplot
# ===============================

plt.figure(figsize=(7,6))

sns.boxplot(
    x="Binder_Type",
    y="Immunogenicity",
    data=data,
    order=["Valid (<500 nM)", "Weak (≥500 nM)"],
    palette=["red", "blue"]
)

# Add p-value text on plot
plt.text(
    0.5,
    max(data["Immunogenicity"]),
    f"Mann–Whitney p = {p_value:.4f}",
    horizontalalignment="center",
    fontsize=10
)

plt.title("Immunogenicity Comparison: Valid vs Weak Binders", fontsize=14)
plt.xlabel("Binder Type", fontsize=12)
plt.ylabel("Immunogenicity Score", fontsize=12)

plt.tight_layout()
plt.savefig("results/Validation_Binder_Comparison_Final.png", dpi=300)
plt.show()
