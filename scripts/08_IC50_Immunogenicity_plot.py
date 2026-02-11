import pandas as pd
import matplotlib.pyplot as plt

# Load data
data = pd.read_csv("results/kras_G12D_IC50_immunogenicity.csv")

# Define valid peptides (IC50 < 500 nM)
valid = data[data["IC50"] < 500]
others = data[data["IC50"] >= 500]

plt.figure(figsize=(7,6))

# Plot other peptides (blue)
plt.scatter(others["IC50"], others["Immunogenicity"],
            color="blue", alpha=0.5, label="Weak Binders (IC50 ≥ 500 nM)")

# Plot valid peptides (red)
plt.scatter(valid["IC50"], valid["Immunogenicity"],
            color="red", edgecolor="black", s=90,
            label="Valid Peptides (IC50 < 500 nM)")

# Add threshold line
plt.axvline(x=500, color="black", linestyle="--", linewidth=1)

plt.xscale("log")

plt.xlabel("IC50 (nM)", fontsize=12)
plt.ylabel("Immunogenicity Score", fontsize=12)
plt.title("KRAS G12D Neoantigen Binding Landscape", fontsize=14)

plt.legend()
plt.grid(alpha=0.2)

for _, row in valid.iterrows():
    plt.text(row["IC50"], row["Immunogenicity"],
             row["peptide"], fontsize=8)


plt.axvspan(1, 500, color='red', alpha=0.05)



plt.tight_layout()
plt.savefig("results/KRAS_G12D_IC50_Immunogenicity.png", dpi=300)
plt.show()





