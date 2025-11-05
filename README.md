# RNA Editing and edQTL Analysis Pipeline (v5/v2) README (Revised)

This pipeline is designed to robustly identify, quantify, filter, and normalize A-to-I RNA editing sites across $\sim\mathbf{6,000}$ individuals and $\mathbf{31}$ cell types for expression quantitative trait locus (edQTL) mapping.

# I. Conceptual Pipeline Overview

The pipeline is organized into six distinct phases, ensuring data quality and statistical rigor at every stage.

Phase,Conceptual Goal,Input Data,Output Data
P1,Reference Setup,"Genome FASTA, GTF",Annotated Blacklists
P2,Raw Calling,"∼6,000 individuals × 31 Cell Type BAMs","∼372,000 Raw Call Files (RED-ML & REDItools)"
P3,Consensus & Annotation,P2 Raw Calls,Individual Consensus Matrix (Raw ERs + Annotation)
P4,Final Filtering & QC,"P3 Matrix, Raw BAMs, Germline VCF",Final QC'd Editing Matrix (Known & UTR3 Sites)
P5,Feature Selection,"P4 Final Matrices (6,000 files)",Population Feature Matrix (Most Active Site ERs)
P6,Normalization & Prep,"P5 Feature Matrix, Covariate Files",Final FastQTL Phenotype Matrix (INT)

# II. Quality Control (QC) Measures (Revised)

The final site list must satisfy all criteria, including those enforced at the initial REDItools calling stage.

Type,Measure,Phase,Script / Tool
Pre-Filtering (REDItools),"Read Quality ≥20, Base Quality ≥20, Read Depth ≥10, and ≥3 Edits supporting the variant.",P2,call_reditools_v4.sh
Artifact Filtering,Blacklisting of Repeats (Alu/Simple) using pre-generated BED files.,"P1, P2","01_setup_references_mask_only_v2.sh, RED-ML/REDItools"
False Positive Control,Dual-Caller Consensus (Sites must be called by both RED-ML and REDItools).,P3,run_phase3_individual_processing_v5.py
Allele Specificity,Filter for canonical A-to-I (A > G or T > C) events only.,P3,run_phase3_individual_processing_v5.py
Site Specificity,Filter for Known Sites (REDIPortal) AND sites within the 3' UTR region.,P4,04_per_cell_type_quantification_and_filter_v3.sh
Genomic Confounding,Exclusion of sites overlapping a germline SNP (using VCF).,P4,run_phase4_quantification_v3.py
Data Reliability (Final),Final minimum Read Coverage ≥10 check during re-quantification.,P4,run_phase4_quantification_v3.py
Population Bias,Inverse Normal Transformation (INT) to standardize phenotype distribution.,P6,normalize_and_covariate_phase6.py
Hidden Variation,"Correction using multiple Covariates (PCs, PEER, AEI, Cell Props).",P6,normalize_and_covariate_phase6.py

# III. Pipeline Scripts and Core Tasks

## A. RNA Editing Calling Pipeline (Master Script: 000_Wraper_v17.sh)

Conceptual Step,Core Task,SLURM Wrapper,Python Script / Tool
P1 Setup,Generate reference blacklists and annotation files.,01_setup_references_mask_only_v2.sh,(Custom Shell/Tool)
P2 Prep,Index and chunk 31 cell-type BAMs per individual.,02_index_pseudobamfiles_v2.sh,(Shell/Samtools)
P2 Calling,Parallel A-to-I Calling (31 tasks per individual) with Strict QC enforced by REDItools.,call_red_ml_v3.sh / call_reditools_v4.sh,"RED-ML, REDItools3"
P3 Consensus,"Aggregate raw calls, enforce Dual Consensus, annotate features.",03_phase3_master_site_discovery_v5.sh,run_phase3_individual_processing_v5.py
P4 Filtering,"Final Filtering (Known/UTR3), Germline Exclusion, Re-quantification.",04_per_cell_type_quantification_and_filter_v3.sh,run_phase4_quantification_v3.py

## B. edQTL Statistical Pipeline (Master Script: master_edQTL_pipeline_v2.sh)

Conceptual Step,Core Task,SLURM Wrapper,Python Script
P5 Selection,"Collate ∼6,000 matrices and select the Most Active Site (highest median ER) per (Gene, CellType).",run_phase5_collation.sh,collate_and_select_phase5.py
P6 Normalization,"Apply Inverse Normal Transformation (INT) row-wise and merge all covariates (PEER, PC, AEI, Cell Props).",run_phase6_processing.sh,normalize_and_covariate_phase6.py