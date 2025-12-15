# neoantigen_identification_vaccine_development

Project Overview
This project implements a computational workflow for identifying candidate tumor neoantigens from cancer genomics data, with the long-term goal of supporting personalized cancer vaccine development.
The pipeline integrates:
* Somatic mutation data
* Tumor RNA-seq expression data
* Patient and clinical metadata
* MHC class I binding prediction
At the current stage, the project focuses on data acquisition, preprocessing, and validation of a prototype neoantigen discovery pipeline using toy datasets, with clear extension points for full patient-level analysis.

Scientific Rationale
Tumor neoantigens arise from somatic mutations that create novel peptides absent in normal tissues. These peptides can trigger an immune response when:
* A nonsynonymous mutation generates a novel peptide
* The peptide binds patient-specific MHC molecules
* The gene is expressed in the tumor
* The mutation is present at sufficient variant allele frequency
This project follows standard immunogenomics strategies used in contemporary neoantigen discovery studies while maintaining modularity and reproducibility.

Workflow Summary
Somatic mutations ? mutant peptide generation ? MHC-I binding prediction ? expression & VAF filtering ? candidate prioritization

Current Project Status
? Data Collection and Preprocessing (Completed)
* Cancer genomics data were downloaded and processed from TCGA using TCGAbiolinks (R)
* Data types collected:
o Somatic mutation data (MAF-derived)
o Tumor RNA-seq expression data
o Clinical and sample metadata
* Controlled-access data limitations were handled appropriately
* Toy datasets were generated where direct patient-level data could not be shared
Prepared output files:
* mut_final.csv — curated somatic mutation table
* tumor_rna_final.csv — tumor gene expression data
* final_patients.csv — patient/sample metadata
* clinical_final.csv — clinical annotations

? Workflow Design (Completed)
A complete conceptual and computational workflow has been designed following standard neoantigen discovery pipelines.
The workflow is modular and allows future integration of:
* Patient-specific HLA typing
* Immunogenicity prediction models
* Machine learning–based candidate re-ranking

? Prototype Neoantigen Pipeline (Implemented with Toy Data)
The following steps have been implemented and validated using simplified datasets:
1. Mutant Peptide Generation
o Parsed protein-level mutation annotations (e.g., p.R248Q)
o Generated overlapping 8–11mer mutant peptides
o Output: candidate_peptides.csv
2. MHC Class I Binding Prediction
o Implemented MHC-I binding prediction using mhcflurry
o Supports patient-specific HLA alleles when available
o Common HLA alleles used for demonstration where needed
o Output: binding_predictions.csv
3. Expression and Variant Allele Frequency Filtering
o Integrated RNA expression (TPM) and mutation VAF
o Applied biologically motivated thresholds
o Output: filtered_candidates.csv
4. Candidate Prioritization
o Implemented a composite scoring system combining:
* Binding strength
* Tumor expression
* Variant allele frequency
o Output: priority_candidates.csv

What This Project Demonstrates
* Understanding of neoantigen biology and cancer immunogenomics
* Integration of multi-omics cancer datasets
* Practical implementation of MHC binding prediction
* Modular and reproducible pipeline design
* Responsible handling of controlled-access data using toy datasets

Planned Extensions (Future Work)
* Patient-specific HLA typing (OptiType / seq2HLA)
* Immunogenicity prediction (PRIME, DeepHLApan)
* Integration of IEDB epitope validation data
* Machine learning–based neoantigen re-ranking
* Visualization and cohort-level summaries

Tools & Technologies Used
* R: TCGAbiolinks (TCGA data acquisition)
* Python: pandas, numpy, biopython
* MHC binding prediction: mhcflurry
* Data sources: TCGA, UniProt

Reproducibility and Data Ethics
* No raw controlled-access patient data is shared
* Toy datasets are used to demonstrate pipeline logic
* All steps are scripted and reproducible
* Designed in compliance with TCGA data usage guidelines

Author Notes
This project represents an ongoing effort to build a research-grade, portfolio-ready neoantigen discovery pipeline. The current implementation focuses on validated core logic and workflow design, with future iterations planned as data access and computational resources expand.

