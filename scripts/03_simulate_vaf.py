import pandas as pd

# Strong binders you already identified
data = [
    {
        "peptide": "KLVVVGADGV",
        "allele": "HLA-A*02:01",
        "vaf": 0.30
    },
    {
        "peptide": "VGADGVGKSAL",
        "allele": "HLA-C*07:02",
        "vaf": 0.15
    }
]

total_depth = 100

for row in data:
    row["t_alt_count"] = int(row["vaf"] * total_depth)
    row["t_ref_count"] = total_depth - row["t_alt_count"]

df = pd.DataFrame(data)

df.to_csv("data/processed/vaf_simulated.csv", index=False)

print("Toy VAF data generated successfully")























































































