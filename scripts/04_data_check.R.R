setwd("C:/Users/shiva/OneDrive/Documents/neoantigen project")
# Load the toy files
toy_maf <- read.csv("TCGA_COAD_toy_top2_MAF.csv")
toy_rna <- read.csv("TCGA_COAD_toy_top2_RNA.csv")
toy_cohort <- read.csv("TCGA_COAD_toy_top2_cohort.csv")

# Get unique patients
patients <- unique(toy_cohort$patient_id)

# Loop: create separate files per patient
for (p in patients) {
  
  # Create folder for patient
  dir.create(paste0("Patient_", p), showWarnings = FALSE)
  
  # Filter data for the patient
  maf_p <- toy_maf[toy_maf$patient_id == p, ]
  rna_p <- toy_rna[toy_rna$patient_id == p, ]
  cohort_p <- toy_cohort[toy_cohort$patient_id == p, ]
  
  # Save files inside patient folder
  write.csv(maf_p, paste0("Patient_", p, "/", p, "_maf.csv"), row.names = FALSE)
  write.csv(rna_p, paste0("Patient_", p, "/", p, "_rna.csv"), row.names = FALSE)
  write.csv(cohort_p, paste0("Patient_", p, "/", p, "_clinical.csv"), row.names = FALSE)
}

cat("🎉 Files successfully separated per patient!\n")
