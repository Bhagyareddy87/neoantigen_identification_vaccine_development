\# Neoantigen Identification \& Vaccine Development Pipeline



\## Project Overview



This repository implements a computational workflow for identifying and prioritizing tumor neoantigens from cancer genomics data, with the long-term goal of supporting personalized cancer vaccine development.



The pipeline integrates:



\- Somatic mutation data

\- Tumor RNA-seq gene expression data

\- Variant allele frequency (VAF) information

\- MHC class I binding predictions



Due to controlled-access restrictions on patient-level genomic data, the current implementation uses toy datasets generated with PI approval to validate pipeline logic. The workflow is designed to be modular, reproducible, and directly extensible to real patient data once access is available.



---



\## Scientific Rationale



Tumor neoantigens arise from somatic, non-synonymous mutations that generate novel peptides absent in normal tissues. These peptides can elicit an immune response when:



\- A protein-altering mutation generates a novel peptide

\- The peptide binds to patient-specific MHC molecules

\- The corresponding gene is expressed in the tumor

\- The mutation is present at sufficient variant allele frequency (VAF)



This project follows standard immunogenomics strategies used in contemporary neoantigen discovery studies, while emphasizing clarity, reproducibility, and responsible data handling.



---



\## Workflow Summary



\*\*Somatic mutations → Mutant peptide generation → MHC class I binding prediction → Expression \& VAF filtering → Neoantigen prioritization\*\*



---



\## HLA Typing Assumptions



Due to controlled-access restrictions on patient-level FASTQ/BAM files, patient-specific HLA typing could not be performed in the current implementation. For pipeline validation and demonstration purposes, common HLA class I alleles (e.g., HLA-A\*02:01, HLA-C\*07:02) were used.



This assumption is explicitly documented and can be replaced with true patient-specific HLA calls (e.g., OptiType or seq2HLA) without structural changes to the pipeline.



---



\## Current Project Status



\### Data Collection \& Preprocessing (Completed)



\- Cancer genomics data were accessed using \*\*TCGAbiolinks (R)\*\*

\- Data types considered:

&nbsp; - Somatic mutation data (MAF-derived)

&nbsp; - Tumor RNA-seq expression data

&nbsp; - Clinical and sample metadata

\- Controlled-access constraints were handled appropriately

\- Toy datasets were generated where patient-level data could not be shared



Prepared input files include:

\- `mut\_final.csv` — curated somatic mutation table

\- `tumor\_rna\_final.csv` — tumor gene expression (TPM)

\- `final\_patients.csv` — patient and sample metadata

\- `clinical\_final.csv` — clinical annotations



---



\### Workflow Design (Completed)



A complete conceptual and computational workflow was designed following established neoantigen discovery pipelines.



---



\## Implemented Neoantigen Discovery Pipeline (Toy Data)



\### Mutant Peptide Generation

\- Parsed protein-level mutation annotations (e.g., `p.R248Q`)

\- Generated overlapping 8–11mer mutant peptides using canonical protein sequences

\- Output: `candidate\_peptides.csv`



\### MHC Class I Binding Prediction

\- Implemented using \*\*mhcflurry\*\*

\- Output: `binding\_predictions.csv`



\### Expression \& Variant Allele Frequency Filtering

\- Integrated RNA expression (TPM) and mutation VAF data

\- Output: `filtered\_candidates.csv`

&nbsp; - Example: `kras\_G12D\_final\_neoantigens.csv`



\### Neoantigen Prioritization

\- Composite scoring using binding, expression, and VAF

\- Output: `priority\_candidates.csv`

&nbsp; - Example: `kras\_G12D\_priority\_neoantigens.csv`



---



\## Planned Extensions (Future Work)



\- Patient-specific HLA typing (OptiType / seq2HLA)

\- Immunogenicity prediction (PRIME, DeepHLApan)

\- IEDB epitope integration

\- Machine learning–based re-ranking

\- Cohort-level visualization



---



\## Tools \& Technologies Used



\- \*\*R\*\*: TCGAbiolinks

\- \*\*Python\*\*: pandas, numpy, biopython

\- \*\*MHC binding\*\*: mhcflurry

\- \*\*Data sources\*\*: TCGA, GENCODE, UniProt



---



\## Reproducibility \& Data Ethics



\- No controlled-access data is shared

\- Toy datasets used for demonstration

\- All steps are scripted and reproducible



---



\## Author Notes



This project represents a research-grade, portfolio-ready neoantigen discovery pipeline.



