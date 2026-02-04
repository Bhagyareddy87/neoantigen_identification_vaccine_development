#Check where Python is currently looking
import os
print(os.getcwd())
print(os.listdir())

##“Compressed GENCODE annotation files 
# were decompressed using Python’s gzip module.”##

import gzip, shutil

input_file = "gencode.v46.pc_translations.fa.gz"
output_file = "gencode.v46.pc_translations.fa"

with gzip.open(input_file, "rb") as f_in:
    with open(output_file, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

print("GENCODE unzipped")

## define mutations 
import pandas as pd
# Filter only missense mutations
maf_final = pd.read_csv("maf_final.csv")
print(maf_final.shape)
maf_final = maf_final[
    maf_final["Variant_Classification"] == "Missense_Mutation"
].copy()






import re
import pandas as pd

def parse_protein_change(pchange):
    if not isinstance(pchange, str):
        return None, None, None
    m = re.match(r"p\.([A-Z])(\d+)([A-Z])", pchange)
    if not m:
        return None, None, None
    return m.group(1), int(m.group(2)), m.group(3)

maf_final[["WT_AA", "Position", "Mut_AA"]] = (
    maf_final["Protein_Change"]
    .apply(lambda x: pd.Series(parse_protein_change(x)))
)

maf_final = maf_final.dropna(subset=["WT_AA"])

print(maf_final)

## check the gencode.v46.pc_translations.fa file description 
from Bio import SeqIO

for record in SeqIO.parse("gencode.v46.pc_translations.fa", "fasta"):
    print(record.description)
    break



##Load GENCODE protein FASTA
from Bio import SeqIO

protein_sequences = {}

for record in SeqIO.parse("gencode.v46.pc_translations.fa", "fasta"):
    parts = record.description.split("|")
    if len(parts) >= 7:
        gene = parts[6].strip()
        seq = str(record.seq)

        # keep the longest isoform per gene
        if gene not in protein_sequences or len(seq) > len(protein_sequences[gene]):
            protein_sequences[gene] = seq

print("Total proteins loaded:", len(protein_sequences))
print("KRAS length:", len(protein_sequences["KRAS"]))










# --- get KRAS protein ---
seq = protein_sequences["KRAS"]

# sanity check
print("WT AA at pos 12:", seq[11])  # should be G


# --- apply mutation ---
def apply_mutation(seq, pos, mut):
    pos = pos - 1  # 1-based → 0-based
    return seq[:pos] + mut + seq[pos+1:]

mut_seq = apply_mutation(seq, 12, "D")
print("Mut AA at pos 12:", mut_seq[11])  # should be D


# --- generate 8–11mer peptides ---
def generate_peptides(seq, pos, lengths=range(8,12)):
    peptides = set()
    pos = pos - 1
    for l in lengths:
        for i in range(max(0, pos-l+1), min(len(seq)-l+1, pos+1)):
            peptides.add(seq[i:i+l])
    return peptides

wt_peptides  = generate_peptides(seq, 12)
mut_peptides = generate_peptides(mut_seq, 12)


# --- true neoantigens ---
neoantigen_peptides = mut_peptides - wt_peptides

print("Number of mutant-only peptides:", len(neoantigen_peptides))
for p in sorted(neoantigen_peptides):
    print(p)



with open("kras_G12D_neoantigens.txt", "w") as f:
    for p in sorted(neoantigen_peptides):
        f.write(p + "\n")

print("Saved neoantigen peptides to kras_G12D_neoantigens.txt")














