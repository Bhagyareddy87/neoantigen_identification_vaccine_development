if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  BiocManager::install("TCGAbiolinks")
}
library(TCGAbiolinks)
# Get mutation (MAF) open samples
query_mut <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)
mut_samples <- unique(substr(getResults(query_mut)$cases,1,12))

# Get RNA-seq open samples
query_rna <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  access = "open"
)
rna_samples <- unique(substr(getResults(query_rna)$cases,1,12))

# Intersection = VALID STUDY PATIENTS
valid_patients <- intersect(mut_samples, rna_samples)
valid_patients[1:20]


valid_patients_df <- data.frame(patient_id = valid_patients)

write.csv(valid_patients_df, "TCGA_COAD_valid_patients.csv", row.names = FALSE)
