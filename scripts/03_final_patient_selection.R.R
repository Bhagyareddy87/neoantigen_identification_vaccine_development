setwd("C:/Users/shiva/OneDrive/Documents/neoantigen project")

library(dplyr)

# Load main merged clinical + metadata file
df <- read.csv("TCGA_COAD_filtered_merged_all_478_with_flags.csv")

# Filter using correct column names
filtered <- df %>%
  filter(
    stage_III_IV == TRUE,
    maf_primary == TRUE,
    rna_primary == TRUE
  )

# Save final valid patient cohort
write.csv(filtered, "TCGA_COAD_final_clean_filtered_cohort.csv", row.names = FALSE)

filtered





cat("\n===== Stage III/IV Patients =====\n")
print(table(df$stage_III_IV))

cat("\n===== Stage III/IV + MAF =====\n")
print(table(df$stage_III_IV, df$maf_primary))

cat("\n===== Stage III/IV + RNA =====\n")
print(table(df$stage_III_IV, df$rna_primary))

cat("\n===== Stage III/IV + BOTH MAF & RNA =====\n")
print(sum(df$stage_III_IV == TRUE & df$maf_primary == TRUE & df$rna_primary == TRUE))








# ---------------------------
# 1) Set Working Directory
# ---------------------------
setwd("C:/Users/shiva/OneDrive/Documents/neoantigen project")

# ---------------------------
# 2) Load Libraries
# ---------------------------
library(dplyr)

# ---------------------------
# 3) Load Main File
# ---------------------------
df <- read.csv("TCGA_COAD_filtered_merged_all_478_with_flags.csv")

cat("\n===== Columns Loaded =====\n")
print(colnames(df))

# ---------------------------
# 4) Inspect Column Values
# ---------------------------
cat("\n===== UNIQUE VALUES: MAF FLAG =====\n")
print(unique(df$maf_primary))

cat("\n===== UNIQUE VALUES: RNA FLAG =====\n")
print(unique(df$rna_primary))

cat("\n===== STAGING VALUES =====\n")
print(unique(df$ajcc_pathologic_stage))

# ------------------------------------------
# 5) Convert MAF/RNA fields to TRUE/FALSE
# ------------------------------------------
df$maf_primary <- df$maf_primary %in% c("TRUE", "True", "YES", "Yes", "Present", "maf", 1, TRUE)
df$rna_primary <- df$rna_primary %in% c("TRUE", "True", "YES", "Yes", "Present", "rna", 1, TRUE)

# ------------------------------------------
# 6) Fix stage column if needed
# ------------------------------------------
df$stage_III_IV <- df$ajcc_pathologic_stage %in% 
  c("Stage III", "Stage IIIC", "Stage IIIB", "Stage IVA", "Stage IVB", "Stage IV")

cat("\n===== SUMMARY AFTER CLEANING =====\n")
cat("\nStage III/IV patients:\n")
print(table(df$stage_III_IV))

cat("\nMAF availability:\n")
print(table(df$maf_primary))

cat("\nRNA availability:\n")
print(table(df$rna_primary))

cat("\nStage III/IV + MAF + RNA:\n")
print(sum(df$stage_III_IV & df$maf_primary & df$rna_primary))

# ------------------------------------------
# 7) Final Filter
# ------------------------------------------
filtered <- df %>%
  filter(stage_III_IV == TRUE,
         maf_primary == TRUE,
         rna_primary == TRUE)

cat("\n===== FINAL FILTERED COHORT COUNT =====\n")
print(nrow(filtered))

# ------------------------------------------
# 8) Save Result
# ------------------------------------------
write.csv(filtered, "TCGA_COAD_final_clean_filtered_cohort.csv", row.names = FALSE)

cat("\n\n*** DONE: Final cohort saved as 'TCGA_COAD_final_clean_filtered_cohort.csv' ***\n")










setwd("C:/Users/shiva/OneDrive/Documents/neoantigen project")

# Load the main metadata file
df <- read.csv("TCGA_COAD_filtered_merged_all_478_with_flags.csv")

# Load mutation and RNA meta files
maf <- read.csv("maf_metadata.csv")      # file name may differ
rna <- read.csv("rna_metadata.csv")      # file name may differ

# Extract cleaned sample identifiers (first 12 chars = patient barcode)
maf$patient_id <- substr(maf$sample_id, 1, 12)
rna$patient_id <- substr(rna$sample_id, 1, 12)
df$patient_id <- substr(df$patient_id, 1, 12)

# Mark availability based on real matching
df$maf_primary <- df$patient_id %in% maf$patient_id
df$rna_primary <- df$patient_id %in% rna$patient_id

# Reapply stage filter
df$stage_III_IV <- df$ajcc_pathologic_stage %in% 
  c("Stage III", "Stage IIIC", "Stage IIIB", "Stage IVA", "Stage IVB", "Stage IV")

filtered <- df %>%
  filter(stage_III_IV == TRUE, maf_primary == TRUE, rna_primary == TRUE)

cat("\n===== Final Patient Count After Correct Matching =====\n")
print(nrow(filtered))

write.csv(filtered, "TCGA_COAD_corrected_filtered_cohort.csv", row.names = FALSE)


# ============================================
# TCGA COAD Toy Data for Top 2 Stage III/IV Patients
# ============================================

# Set working directory
setwd("C:/Users/shiva/OneDrive/Documents/neoantigen project")

library(dplyr)

# Load metadata
df <- read.csv("TCGA_COAD_filtered_merged_all_478_with_flags.csv")

# -----------------------------
# 1. Filter Stage III/IV patients
# -----------------------------
stage34_patients <- df %>%
  filter(stage_III_IV == TRUE)

cat("Number of Stage III/IV patients:", nrow(stage34_patients), "\n")

# -----------------------------
# 2. Select top 2 patients
# -----------------------------
top2_patients <- stage34_patients[1:2, ]
cat("Selected top 2 patients for toy data:\n")
print(top2_patients$patient_id)

# -----------------------------
# 3. Generate toy MAF data
# -----------------------------
set.seed(123)

genes <- c("TP53", "KRAS", "APC", "PIK3CA", "SMAD4")

maf_data <- data.frame(
  patient_id = rep(top2_patients$patient_id, each = 3),  # 3 mutations per patient
  Hugo_Symbol = sample(genes, 6, replace = TRUE),
  Variant_Classification = sample(c("Missense_Mutation", "Nonsense_Mutation", "Silent"), 6, replace = TRUE),
  Chromosome = sample(c(1:22, "X", "Y"), 6, replace = TRUE),
  Start_Position = sample(1:1e6, 6, replace = TRUE)
)

write.csv(maf_data, "TCGA_COAD_toy_top2_MAF.csv", row.names = FALSE)
cat("Toy MAF data for top 2 patients saved: TCGA_COAD_toy_top2_MAF.csv\n")

# -----------------------------
# 4. Generate toy RNA data
# -----------------------------
rna_genes <- c("TP53", "KRAS", "APC", "PIK3CA", "SMAD4", "EGFR", "BRAF", "NRAS", "CDKN2A", "PTEN")

rna_data <- expand.grid(patient_id = top2_patients$patient_id,
                        gene = rna_genes)

rna_data$TPM <- round(runif(nrow(rna_data), 0, 100), 2)

write.csv(rna_data, "TCGA_COAD_toy_top2_RNA.csv", row.names = FALSE)
cat("Toy RNA data for top 2 patients saved: TCGA_COAD_toy_top2_RNA.csv\n")

# -----------------------------
# 5. Create final toy cohort file
# -----------------------------
top2_patients$maf_primary <- TRUE
top2_patients$rna_primary <- TRUE

write.csv(top2_patients, "TCGA_COAD_toy_top2_cohort.csv", row.names = FALSE)
cat("Final toy cohort for top 2 patients saved: TCGA_COAD_toy_top2_cohort.csv\n")



#####toy_data_generation#####


# Load the files
toy_maf <- read.csv("TCGA_COAD_toy_top2_MAF.csv")
toy_rna <- read.csv("TCGA_COAD_toy_top2_RNA.csv")
toy_cohort <- read.csv("TCGA_COAD_toy_top2_cohort.csv")

# Inspect first few rows
head(toy_maf)
head(toy_rna)
head(toy_cohort)

# Number of rows
nrow(toy_maf)     # Should be 2 patients × 3 mutations = 6
nrow(toy_rna)     # Should be 2 patients × 10 genes = 20
nrow(toy_cohort)  # Should be 2



unique(toy_maf$patient_id)    # Should list your top 2 patients
unique(toy_rna$patient_id)    # Should list same 2 patients
toy_cohort$patient_id          # Should show same 2 patients




colnames(toy_maf)      # Should have patient_id, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position
colnames(toy_rna)      # Should have patient_id, gene, TPM
colnames(toy_cohort)   # Should have all metadata columns + maf_primary, rna_primary







































































