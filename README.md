# Neoantigen\_identification\_vaccine\_development

Project Overview



This repository implements a computational workflow for identifying and prioritizing tumor neoantigens from cancer genomics data, with the long-term goal of supporting personalized cancer vaccine development.



The pipeline integrates multiple data modalities, including:



Somatic mutation data



Tumor RNA-seq gene expression data



Variant allele frequency (VAF) information



MHC class I binding predictions



Due to controlled-access restrictions on patient-level genomic data, the current implementation uses toy datasets generated with PI approval to validate pipeline logic. The workflow is designed to be modular, reproducible, and directly extensible to real patient data once access is available.



Scientific Rationale



Tumor neoantigens arise from somatic, non-synonymous mutations that generate peptides absent in normal tissues. These peptides can elicit an immune response when:



A protein-altering mutation generates a novel peptide



The peptide binds to patient-specific MHC molecules



The corresponding gene is expressed in the tumor



The mutation is present at sufficient clonal frequency (VAF)



This project follows standard immunogenomics strategies used in contemporary neoantigen discovery studies, while emphasizing clarity, modularity, and data ethics.



Workflow Summary



Somatic mutations → Mutant peptide generation → MHC class I binding prediction → Expression \& VAF filtering → Neoantigen prioritization



Current Project Status

Data Collection \& Preprocessing (Completed)



Cancer genomics data were accessed using TCGAbiolinks (R)



Data types considered:



Somatic mutation data (MAF-derived)



Tumor RNA-seq expression data



Clinical and sample metadata



Controlled-access constraints were handled appropriately



Toy datasets were generated where patient-level data could not be shared



Prepared input files include:



mut\_final.csv — curated somatic mutation table



tumor\_rna\_final.csv — tumor gene expression (TPM)



final\_patients.csv — patient/sample metadata



clinical\_final.csv — clinical annotations



Workflow Design (Completed)



A complete conceptual and computational workflow was designed following established neoantigen discovery pipelines. The design supports future integration of:



Patient-specific HLA typing



Immunogenicity prediction models



Machine learning–based candidate re-ranking



Prototype Neoantigen Pipeline (Implemented with Toy Data)



The following core steps have been fully implemented and validated:



Mutant Peptide Generation



Parsed protein-level mutation annotations (e.g., p.R248Q)



Generated overlapping 8–11mer mutant peptides using canonical protein sequences (gencode.v46.pc\_translations.fa)



Output: candidate\_peptides.csv



MHC Class I Binding Prediction



Implemented MHC-I binding prediction using mhcflurry



Supports patient-specific HLA alleles when available



Common HLA alleles used for demonstration where required



Output: binding\_predictions.csv



Expression \& Variant Allele Frequency Filtering



Integrated RNA expression (TPM) and mutation VAF data



Applied biologically motivated thresholds (TPM, VAF, IC50)



Output: filtered\_candidates.csv



Neoantigen Prioritization



Implemented a composite neoantigen scoring system combining:



MHC binding strength



Tumor gene expression



Variant allele frequency



Ranked candidates to identify highest-priority neoantigens



Output: priority\_candidates.csv



What This Project Demonstrates



Understanding of neoantigen biology and cancer immunogenomics



Integration of multi-omics cancer datasets



Practical implementation of MHC class I binding prediction



End-to-end neoantigen filtering and prioritization logic



Modular, reproducible pipeline design



Responsible handling of controlled-access data using toy datasets



Planned Extensions (Future Work)



Patient-specific HLA typing (OptiType / seq2HLA)



Immunogenicity prediction (PRIME, DeepHLApan)



Integration of IEDB epitope validation data



Machine learning–based neoantigen re-ranking



Visualization and cohort-level summaries



Tools \& Technologies Used



R: TCGAbiolinks (TCGA data acquisition)



Python: pandas, numpy, biopython



MHC binding prediction: mhcflurry



Data sources: TCGA, GENCODE, UniProt



Reproducibility \& Data Ethics



No raw controlled-access patient data is shared



Toy datasets are used to demonstrate pipeline logic



All analysis steps are scripted and reproducible



Designed in compliance with TCGA data usage guidelines



Author Notes



This project represents an ongoing effort to build a research-grade, portfolio-ready neoantigen discovery pipeline.

The current implementation focuses on validated core logic and prioritization, with future iterations planned as data access and computational resources expand.

