# Full production-grade patient filtering script for neoantigen pipeline
# Uses TCGAbiolinks for metadata only (no big downloads)
# Author: ChatGPT (adapted for your project)
# Date: 2025-12-01

# -----------------------------
# Libraries
# -----------------------------
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  BiocManager::install("TCGAbiolinks")
}
library(TCGAbiolinks)
library(dplyr)
library(readr)
library(stringr)
library(rlang)

# -----------------------------
# Parameters / filenames
# -----------------------------
patient_file <- "TCGA_COAD_valid_patients.csv"   # your 478 patient IDs (one column)
out_prefix <- "TCGA_COAD_filtered"               # prefix for output files

# -----------------------------
# Step 0: Read patient list
# -----------------------------
patients <- read_csv(patient_file, col_names = TRUE)
# ensure column name standardized
if(!"patient_id" %in% colnames(patients)) {
  colnames(patients)[1] <- "patient_id"
}
# normalize patient ids to first 12 chars (TCGA-XX-YYYY)
patients <- patients %>%
  mutate(patient_id = str_trim(patient_id),
         patient_id_12 = substr(patient_id, 1, 12))

cat("Loaded", nrow(patients), "patient IDs (first 12 chars used as patient key)\n")

# -----------------------------
# Step 1: Get clinical metadata (small)
# -----------------------------
cat("Downloading clinical metadata (small)...\n")
clinical <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")

cat("Clinical rows:", nrow(clinical), "\n")

# -----------------------------
# Helper: robust column detection
# -----------------------------
detect_first <- function(candidates, df) {
  found <- intersect(candidates, colnames(df))
  if(length(found) == 0) return(NA_character_)
  found[1]
}

# detect plausible patient id column in clinical
patient_col_candidates <- c("bcr_patient_barcode", "submitter_id", "patient_id", "case_submitter_id")
patient_clin_col <- detect_first(patient_col_candidates, clinical)
if(is.na(patient_clin_col)) {
  # try to guess: some clinical tables use 'cases' nested; we'll handle later
  cat("Warning: could not find an explicit patient-id column in clinical. We'll attempt to merge using alternate approaches.\n")
} else {
  cat("Detected patient column in clinical:", patient_clin_col, "\n")
}

# detect age, stage and msi columns
age_col_candidates <- c("age_at_diagnosis","age_at_initial_diagnosis",
                        "age_at_initial_pathologic_diagnosis",
                        "demographic.age_at_index", "diagnosis.age_at_diagnosis",
                        "age_at_index")
age_col <- detect_first(age_col_candidates, clinical)
stage_col_candidates <- c("clinical_stage","ajcc_pathologic_stage","tumor_stage","pathologic_stage")
stage_col <- detect_first(stage_col_candidates, clinical)
msi_col_candidates <- c("msi_status","microsatellite_instability","mismatch_repair_status","msi")
msi_col <- detect_first(msi_col_candidates, clinical)

cat("Detected columns -> age:", age_col, " stage:", stage_col, " msi:", msi_col, "\n")

# -----------------------------
# Step 2: Merge patients + clinical safely
# -----------------------------
# Try to merge by common patient id columns, with several fallbacks
if(!is.na(patient_clin_col)) {
  merged <- merge(patients, clinical, by.x = "patient_id_12", by.y = patient_clin_col, all.x = TRUE)
} else {
  # fallback: try merge by 12-char prefix if clinical contains any column with TCGA- patterns
  # Create a 12-char patient id column from clinical if possible
  created_key <- FALSE
  for (cn in colnames(clinical)) {
    if(any(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", clinical[[cn]], perl = TRUE), na.rm = TRUE)) {
      clinical <- clinical %>% mutate(clin_patient_12 = substr(as.character(!!sym(cn)),1,12))
      merged <- merge(patients, clinical, by.x = "patient_id_12", by.y = "clin_patient_12", all.x = TRUE)
      created_key <- TRUE
      break
    }
  }
  if(!created_key) {
    # Last resort: join by exact 12-char patient ids if clinical has a column named 'patient_id' but different shape
    clinical$clin_tmp_id <- NA
    merged <- merge(patients, clinical, by.x = "patient_id_12", by.y = "clin_tmp_id", all.x = TRUE)
    cat("Warning: merge may have failed to attach clinical fields. Inspect 'merged' manually.\n")
  }
}

cat("Merged rows (should equal number of patients):", nrow(merged), "\n")

# -----------------------------
# Step 3: Query GDC metadata to check MAF and STAR counts availability (metadata only)
# -----------------------------
cat("Querying GDC for MAF metadata (masked somatic mutations)...\n")
q_maf <- GDCquery(project = "TCGA-COAD",
                  data.category = "Simple Nucleotide Variation",
                  data.type = "Masked Somatic Mutation")
maf_meta <- tryCatch(getResults(q_maf), error = function(e) {
  cat("Error retrieving MAF metadata: ", e$message, "\n"); NULL
})

cat("Querying GDC for RNA-Seq (STAR - Counts / STAR) metadata...\n")
q_rna <- GDCquery(project = "TCGA-COAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")
rna_meta <- tryCatch(getResults(q_rna), error = function(e) {
  cat("Error retrieving RNA metadata: ", e$message, "\n"); NULL
})

# Helper to robustly extract patient submitter IDs from the metadata returned by getResults()
extract_patient_ids_from_meta <- function(meta_df) {
  if(is.null(meta_df)) return(character(0))
  ids <- rep(NA_character_, nrow(meta_df))
  
  # Approach 1: nested 'cases' list-column -> attempt to extract submitter_id
  if("cases" %in% colnames(meta_df)) {
    safe_extract <- function(x) {
      tryCatch({
        if(is.list(x) && length(x) >= 1) {
          # first case's submitter_id if present
          case <- x[[1]]
          if(!is.null(case$submitter_id)) return(as.character(case$submitter_id))
          # try nested sample
          if(!is.null(case$case_id)) return(as.character(case$case_id))
        }
        NA_character_
      }, error = function(e) NA_character_)
    }
    ids <- sapply(meta_df$cases, safe_extract)
  }
  
  # Approach 2: file_name may contain TCGA barcode (fallback)
  if(all(is.na(ids)) || sum(!is.na(ids)) < length(ids)/2) {
    if("file_name" %in% colnames(meta_df)) {
      # look for a TCGA-XX-XXXX pattern
      extracted <- str_extract(meta_df$file_name, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")
      ids[is.na(ids)] <- extracted[is.na(ids)]
    }
  }
  
  # Approach 3: submitter_id or cases.submitter_id column
  candidate_cols <- c("submitter_id","cases.submitter_id","cases.submitter_id.1")
  for(col in candidate_cols) {
    if(col %in% colnames(meta_df)) {
      ids[is.na(ids)] <- as.character(meta_df[[col]][is.na(ids)])
    }
  }
  
  # final normalize: first 12 chars for patient id
  ids <- ifelse(is.na(ids), NA_character_, substr(as.character(ids), 1, 12))
  ids[!nzchar(ids)] <- NA_character_
  unique(na.omit(ids))
}

maf_patients <- extract_patient_ids_from_meta(maf_meta)
rna_patients <- extract_patient_ids_from_meta(rna_meta)

cat("Number of unique patients with MAF metadata:", length(maf_patients), "\n")
cat("Number of unique patients with STAR counts metadata:", length(rna_patients), "\n")

# -----------------------------
# Step 4: Determine which patients have required files & primary tumor sample
# -----------------------------
# We will attempt to detect sample type 'Primary Tumor' from the metadata (nested structure) if present

extract_primary_sample_flag <- function(meta_df) {
  # Returns a named logical vector: names = patient_id_12
  if(is.null(meta_df)) return(logical(0))
  patient_flags <- list()
  for(i in seq_len(nrow(meta_df))) {
    rec <- meta_df[i, ]
    pid <- NA_character_
    
    # try cases -> sample_type
    if("cases" %in% colnames(meta_df)) {
      cs <- rec$cases[[1]]
      pid <- tryCatch(cs$submitter_id, error = function(e) NA_character_)
      # sample type may be nested under cs[[1]]$sample_type
      sample_type <- NA_character_
      if(!is.null(cs) && length(cs) >= 1) {
        st <- tryCatch(cs[[1]]$sample_type, error = function(e) NA_character_)
        sample_type <- st
      }
      if(!is.na(pid)) {
        pid12 <- substr(as.character(pid), 1, 12)
        if(is.null(patient_flags[[pid12]])) patient_flags[[pid12]] <- FALSE
        if(!is.na(sample_type) && grepl("Primary", sample_type, ignore.case = TRUE)) patient_flags[[pid12]] <- TRUE
      }
    }
    
    # fallback: if file_name contains sample barcode TCGA-XX-XXXX-01...
    if(is.na(pid) && "file_name" %in% colnames(meta_df)) {
      fn <- as.character(rec$file_name)
      bc <- str_extract(fn, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[A-Z0-9]{2}[A-Z0-9]")
      if(!is.na(bc)) {
        pid12 <- substr(bc, 1, 12)
        sample_code <- substr(bc, 14, 15)
        if(is.null(patient_flags[[pid12]])) patient_flags[[pid12]] <- FALSE
        # sample code "01" means Primary Tumor
        if(sample_code == "01") patient_flags[[pid12]] <- TRUE
      }
    }
  }
  # convert list to named logical vector
  if(length(patient_flags) == 0) return(logical(0))
  v <- unlist(patient_flags)
  names(v) <- names(patient_flags)
  v
}

maf_primary_flags <- extract_primary_sample_flag(maf_meta)
rna_primary_flags <- extract_primary_sample_flag(rna_meta)

# -----------------------------
# Step 5: Attach availability flags to 'merged' patient table
# -----------------------------
merged <- merged %>%
  mutate(patient_id_12 = substr(patient_id, 1, 12),
         has_maf_meta = patient_id_12 %in% maf_patients,
         has_rna_meta = patient_id_12 %in% rna_patients,
         maf_primary = ifelse(patient_id_12 %in% names(maf_primary_flags),
                              as.logical(maf_primary_flags[patient_id_12]), FALSE),
         rna_primary = ifelse(patient_id_12 %in% names(rna_primary_flags),
                              as.logical(rna_primary_flags[patient_id_12]), FALSE)
  )

cat("Patients with MAF metadata in your 478:", sum(merged$has_maf_meta, na.rm = TRUE), "\n")
cat("Patients with STAR counts metadata in your 478:", sum(merged$has_rna_meta, na.rm = TRUE), "\n")
cat("Patients with MAF primary sample flag:", sum(merged$maf_primary, na.rm = TRUE), "\n")
cat("Patients with RNA primary sample flag:", sum(merged$rna_primary, na.rm = TRUE), "\n")

# -----------------------------
# Step 6: Compute stage-based filtering (Stage III/IV)
# -----------------------------
# Use the detected stage column; if not available, try common alternative columns
if(is.na(stage_col)) {
  # Attempt to find any column containing "stage" in its name
  cand <- grep("stage", colnames(merged), ignore.case = TRUE, value = TRUE)
  if(length(cand) > 0) stage_col <- cand[1]
}

if(is.na(stage_col)) {
  cat("Warning: no stage column found. Stage filtering not applied. Inspect merged columns manually.\n")
  # create a placeholder column to avoid errors
  merged$stage_for_filter <- NA_character_
  stage_col <- "stage_for_filter"
} else {
  # ensure stage column exists in merged (may be named differently after merge)
  if(!stage_col %in% colnames(merged)) {
    # try with prefix from clinical
    possible <- intersect(colnames(clinical), colnames(merged))
    stage_col <- detect_first(stage_col_candidates, merged)
  }
}

# Normalize clinical stage values to common names (best-effort)
normalize_stage <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_
  x <- trimws(x)
  x <- ifelse(grepl("^I\\b", x, ignore.case = TRUE), "Stage I", x)
  x <- ifelse(grepl("^II\\b", x, ignore.case = TRUE), "Stage II", x)
  x <- ifelse(grepl("^III\\b", x, ignore.case = TRUE), "Stage III", x)
  x <- ifelse(grepl("^IV\\b", x, ignore.case = TRUE), "Stage IV", x)
  x
}

merged <- merged %>%
  mutate(stage_norm = normalize_stage(!!sym(stage_col)))

# Filter valid stage entries and Stage III/IV
valid_stage_values <- c("Stage I","Stage II","Stage III","Stage IV")
merged <- merged %>%
  mutate(stage_valid = stage_norm %in% valid_stage_values,
         stage_III_IV = stage_norm %in% c("Stage III","Stage IV"))

cat("Patients with non-missing stage:", sum(merged$stage_valid, na.rm = TRUE), "\n")
cat("Patients Stage III/IV:", sum(merged$stage_III_IV, na.rm = TRUE), "\n")

# -----------------------------
# Step 7: MSI flagging (if column exists)
# -----------------------------
if(!is.na(msi_col) && msi_col %in% colnames(merged)) {
  merged <- merged %>%
    mutate(msi_flag = case_when(
      grepl("MSI-H|MSI High|MSI high|microsatellite instability-high|dMMR", as.character(!!sym(msi_col)), ignore.case = TRUE) ~ "MSI-H",
      grepl("MSS|MSI-L|MSI-Low|microsatellite stable", as.character(!!sym(msi_col)), ignore.case = TRUE) ~ "MSS",
      TRUE ~ NA_character_
    ))
  cat("MSI column found and MSI flags created.\n")
} else {
  merged$msi_flag <- NA_character_
  cat("MSI column not found; msi_flag left NA. You can add MSI info later.\n")
}

# -----------------------------
# Step 8: Final selection rules (Option B logic)
# -----------------------------
# Keep Stage III/IV patients, prefer primary tumor samples, require MAF + RNA metadata available and primary sample if possible
final_candidates <- merged %>%
  filter(stage_III_IV == TRUE) %>%
  mutate(has_both_meta = has_maf_meta & has_rna_meta,
         primary_sample_both = maf_primary & rna_primary)

cat("Stage III/IV count (before MAF/RNA availability):", nrow(merged %>% filter(stage_III_IV == TRUE)), "\n")
cat("Stage III/IV with both MAF and RNA metadata available:", sum(final_candidates$has_both_meta, na.rm = TRUE), "\n")
cat("Stage III/IV with primary sample for both MAF and RNA:", sum(final_candidates$primary_sample_both, na.rm = TRUE), "\n")

# Priority cohort: Stage III/IV + MAF + RNA + primary sample both
priority_cohort <- final_candidates %>% filter(primary_sample_both == TRUE)

# Broader analysis cohort: Stage III/IV + MAF + RNA (even if sample-type unknown)
analysis_cohort <- final_candidates %>% filter(has_both_meta == TRUE)

# Stage III/IV but missing MAF or RNA (to inspect)
missing_data <- final_candidates %>% filter(!has_both_meta)

# -----------------------------
# Step 9: Save outputs
# -----------------------------
write_csv(merged, paste0(out_prefix, "_merged_all_478_with_flags.csv"))
write_csv(final_candidates, paste0(out_prefix, "_stageIII_IV_with_meta_flags.csv"))
write_csv(priority_cohort, paste0(out_prefix, "_priority_cohort_stageIII_IV_primary_MAF_RNA.csv"))
write_csv(analysis_cohort, paste0(out_prefix, "_analysis_cohort_stageIII_IV_MAF_RNA.csv"))
write_csv(missing_data, paste0(out_prefix, "_stageIII_IV_missing_MAF_or_RNA.csv"))

cat("Saved CSV outputs:\n",
    paste0(out_prefix, "_merged_all_478_with_flags.csv"), "\n",
    paste0(out_prefix, "_stageIII_IV_with_meta_flags.csv"), "\n",
    paste0(out_prefix, "_priority_cohort_stageIII_IV_primary_MAF_RNA.csv"), "\n",
    paste0(out_prefix, "_analysis_cohort_stageIII_IV_MAF_RNA.csv"), "\n",
    paste0(out_prefix, "_stageIII_IV_missing_MAF_or_RNA.csv"), "\n"
)

# -----------------------------
# Done: Short summary printed
# -----------------------------
cat("\nSummary:\n")
cat("Total input patients:", nrow(patients), "\n")
cat("Stage III/IV patients:", nrow(merged %>% filter(stage_III_IV == TRUE)), "\n")
cat("Stage III/IV with both MAF & RNA metadata:", nrow(analysis_cohort), "\n")
cat("Priority cohort (Stage III/IV + primary samples present + MAF + RNA):", nrow(priority_cohort), "\n")
cat("Review", paste0(out_prefix, "_merged_all_478_with_flags.csv"), "for detailed flags and missing fields.\n\n")

# End of script
