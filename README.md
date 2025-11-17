# RNA Editing and edQTL Analysis Pipeline (v5/v2) README (Revised)

This pipeline is designed to robustly identify, quantify, filter, and normalize A-to-I RNA editing sites across 6,000 individuals and 31 cell types for expression quantitative trait locus (edQTL) mapping.

### Overall Project Rationale (edQTL Importance)

The decision to map edQTLs is strongly supported by the literature:


Disease Relevance: edQTLs are highly enriched in GWAS signals for complex traits, particularly autoimmune and immune-related diseases (e.g., IBD, lupus, rheumatoid arthritis).


Effect Size: edQTLs showed effects of larger magnitude than eQTLs or sQTLs in many contexts and were more enriched in disease heritability.


Regulatory Mechanism: Genetic variants (edQTL SNPs) can alter editing levels by affecting ADAR-binding strength and the formation of RNA secondary structures (like cis-NATs), which are often immunogenic


# I. Conceptual Pipeline Overview

The pipeline is organized into six distinct phases, ensuring data quality and statistical rigor at every stage.

## Phase,Conceptual Goal,Input Data,Output Data

### P1,Reference Setup,"Genome FASTA, GTF",Annotated Blacklists
### P2,Raw Calling,"∼6,000 individuals × 31 Cell Type BAMs","∼372,000 Raw Call Files (RED-ML & REDItools)"
#### RED-ML: https://github.com/BGIRED/RED-ML
#### REDITools: https://github.com/BioinfoUNIBA/REDItools3
### P3,Consensus & Annotation,P2 Raw Calls,Individual Consensus Matrix (Raw ERs + Annotation)
#### REDIPortal: 04/11/25 http://srv00.recas.ba.infn.it/atlas/ MISSING!!!!
### P4,Final Filtering, QC ^ cleanup,"P3 Matrix, Raw BAMs, Germline VCF",Final QC'd Editing Matrix (Known & UTR3 Sites)
### P5,Feature Selection,"P4 Final Matrices (6,000 files)",Population Feature Matrix (Most Active Site ERs)
### P6,Normalization & Prep,"P5 Feature Matrix, Covariate Files",Final FastQTL Phenotype Matrix (INT)
### Association Mapping	P6 Phenotype & Covariates, Genotype VCF	Raw FastQTL Results (Nominal & Permutation P-values)
### FDR Correction & Lead SNP ID	P7 Raw Results	Final Significant Lead edQTL List (Q-value < 0.05)


# II. Quality Control (QC) Measures (Revised)

The final site list must satisfy all criteria, including those enforced at the initial REDItools calling stage.

## Type,Measure,Phase,Script / Tool

### Pre-Filtering (REDItools),"Read Quality ≥20, Base Quality ≥20, Read Depth ≥10, and ≥3 Edits supporting the variant.",P2,call_reditools_v4.sh
### Artifact Filtering,Blacklisting of Repeats (Alu/Simple) using pre-generated BED files.,"P1, P2","01_setup_references_mask_only_v2.sh, RED-ML/REDItools"
### False Positive Control,Dual-Caller Consensus (Sites must be called by both RED-ML and REDItools).,P3,run_phase3_individual_processing_v5.py
### Allele Specificity,Filter for canonical A-to-I (A > G or T > C) events only.,P3,run_phase3_individual_processing_v5.py
### Site Specificity,Filter for Known Sites (REDIPortal) AND sites within the 3' UTR region.,P4,04_per_cell_type_quantification_and_filter_v3.sh
### Genomic Confounding,Exclusion of sites overlapping a germline SNP (using VCF).,P4,run_phase4_quantification_v3.py
### Data Reliability (Final),Final minimum Read Coverage ≥10 check during re-quantification.,P4,run_phase4_quantification_v3.py
### Population Bias,Inverse Normal Transformation (INT) to standardize phenotype distribution.,P6,normalize_and_covariate_phase6.py
### Hidden Variation,"Correction using multiple Covariates (PCs, PEER, AEI, Cell Props).",P6,normalize_and_covariate_phase6.py
### Statistical FilteringFalse Discovery Rate (FDR) control using Q-value (FDR $\le 0.05$).P8process_fastqtl_results_p8.py

# III. Pipeline Scripts and Core Tasks

## A. RNA Editing Calling Pipeline (Master Script: 000_Wraper_v17.sh)

### Conceptual Step,Core Task,SLURM Wrapper,Python Script / Tool

#### P1 Setup,Generate reference blacklists and annotation files.,01_setup_references_mask_only_v2.sh,(Custom Shell/Tool)
#### P2 Prep,Index and chunk 31 cell-type BAMs per individual.,02_index_pseudobamfiles_v2.sh,(Shell/Samtools)
#### P2 Calling,Parallel A-to-I Calling (31 tasks per individual) with Strict QC enforced by REDItools.,call_red_ml_v3.sh / call_reditools_v4.sh,"RED-ML, REDItools3"
#### P3 Consensus,"Aggregate raw calls, enforce Dual Consensus, annotate features.",03_phase3_master_site_discovery_v5.sh,run_phase3_individual_processing_v5.py
#### P4 Filtering,"Final Filtering (Known/UTR3), Germline Exclusion, Re-quantification.",04_per_cell_type_quantification_and_filter_v3.sh,run_phase4_quantification_v3.py
#### A final Cleanup Phase (run_cleanup_p1_p4.sh) was added as the last job in the 000_Wraper_v17.sh. This job is responsible for deleting large, unnecessary intermediate files (including RED-ML's variation.sites.feature.txt and mut.txt.gz) to preserve disk quota before the statistical pipeline begins.

## B. edQTL Statistical Pipeline (Master Script: master_edQTL_pipeline_v2.sh)

### Phase,Script,Description,Dependencies
#### 5,run_phase5_collation_v2.sh,Collates raw editing calls into the final phenotype matrix.,→ (None)
#### AEI-Calc,run_AEI_calculation_array.sh,Calculates the raw AEI (Alu Editing Index) covariate.,→ P5
#### 6,run_phase6_processing_v2.sh,"Normalization & Covariate Merge. Applies INT to edQTL phenotypes. Merges all covariates (PCs, PEER, AEI).",→ P5 AND AEI-Calc
#### 7,run_phase7_edqtl_mapping_v2.sh,edQTL Mapping (FastQTL). Maps variants to INT-normalized editing sites.,→ P6
#### 7b,run_phase7b_aeiqtl_mapping_v3.sh,AEI-QTL Mapping (FastQTL). Maps variants to the AEI covariate as a trans-effect (runs parallel to P7).,→ P6
#### 8,run_phase8_qvalue_filter_v2.sh,"FDR Correction. Performs combined Benjamini-Hochberg FDR correction on P7 and P7b results, and identifies lead SNPs.",→ P7 AND P7b. Phase 8 is configured with a dual dependency (afterok:P7_ID:P7B_ID) to ensure all results are available for a single, unified FDR correction, maintaining statistical rigor across both QTL analyses.

# IV. Justification of Pipeline Decisions

## Justification for Dual Callers and Initial QC (P2 & P3)

### Pipeline Component,Literature Justification,Source

#### P2: Dual Caller Approach (REDItools ∩ RED-ML),"Enforcing a consensus between multiple robust callers is a critical strategy to filter out sequencing artifacts and caller-specific errors, which is essential for discovering true ADAR editing sites with high confidence.","Wang, F., et al. (2023). BMC Biol. 21, 160."
#### "P2: Strict Pre-Filtering (BaseQ ≥20, ReadQ ≥20, Min Edits ≥3)","Strict filtering by base quality, read quality, and minimum edited read count is explicitly required to distinguish true editing events from background noise and low-confidence candidates, a method championed by recent high-stringency studies.","Wang, F., et al. (2023). BMC Biol. 21, 160."
#### P3: Allele Specificity (A > G or T > C),"Confining the analysis to the canonical A-to-I editing type is necessary, as non-canonical events are frequently indicative of sequencing artifacts or somatic mutations.",General RNA Editing Standards

## Justification for Filtering and Site Selection (P1, P4, & P5)

### Pipeline Component,Literature Justification,Source

#### P1/P2: Blacklisting Repeats (Alu/Simple),More than 99% of human editing sites are located in inverted repeat Alu elements. General blacklisting controls for alignment ambiguity outside of these structured regions.,"Li, Q., et al. (2022). Nature 608"
#### P4: Germline SNP Exclusion (VCF),Excluding candidate sites that overlap with known SNPs (MAF>1%) is standard practice to prevent misidentification of germline polymorphisms as novel RNA editing events.,"General Genomic QC Principles (Informed by SComatic, etc.)"
#### P4: Known Sites AND 3' UTR Focus,"This decision is based on two key findings: 1) The validation ratio for annotated (Known) sites is significantly higher than unannotated sites (e.g., 38.5% vs. 7.3%). 2) The 3' UTR is a region known to be enriched in editing annotations and is critical for miRNA regulation.","Wang, F., et al. (2023). BMC Biol. 21, 160.; Gabay, O., et al. (2022). Nat Commun 13, 1184."
#### P5: Most Active Site Selection,"Selecting the most actively edited site per gene/cell-type combination simplifies association testing while reliably capturing the locus's primary regulatory signal, as multiple editing sites are often co-regulated.","Li, Q., et al. (2022). Nature 608"

## Justification for Normalization and Covariates (P6)

### Pipeline Component,Literature Justification,Source

#### P6: Inverse Normal Transformation (INT),"The raw editing ratios must be transformed to a standard, near-normal distribution to meet the assumptions of linear models used in edQTL software (like FastQTL).",Standard eQTL/edQTL Methodology
#### P6: Covariates (PCs, AEI, PEER, Cell Props)","Including Genotype PCs and metrics of global ADAR activity (like the Alu Editing Index or AEI) is essential to regress out global editing variation, as ADAR1 expression is highly correlated with the top principal components of editing variance.","Li, Q., et al. (2022). Nature 608"
#### P7: $\text{Cis-Window} (1\text{ Mb})$The $1\text{ Mb}$ cis-window size is a standard in QTL studies and aligns with the methodology used to discover thousands of edQTLs in human tissues.Li, Q., et al. (2022). Nature 608
#### P7: Permutation TestingThe use of permutation testing (e.g., 1,000 permutations) is standard practice in QTL analysis to accurately derive empirical P-values, which correctly accounts for the number of genetic variants tested per feature in the cis-window.Standard eQTL/edQTL Methodology
#### P8: FDR Correction (Q-value)Applying the Benjamini–Hochberg procedure to the empirical P-values is necessary to control the False Discovery Rate (FDR) across the millions of feature-SNP tests performed across all cell types, preventing the inflation of false positives.Standard eQTL/edQTL Methodology (Informed by Li, Q., et al. (2022))

# V original ideas:

## 1. AEI-QTL Mapping (Phase 7b): 

Novel Phenotype AnalysisThe most novel aspect is the decision to run a separate Quantitative Trait Loci (QTL) study using the Alu Editing Index (AEI) as the phenotype.Rationale: The AEI is a genome-wide metric of overall ADAR enzymatic activity. By mapping genetic variants to this metric, you are looking for trans-acting regulators—genes or pathways other than ADAR1 or ADAR2—that influence the global scale of RNA editing.Originality: While Li et al. (2022) focus on finding edQTLs (variants affecting specific editing sites), they do not report this systematic attempt to find AEI-QTLs. This analysis provides a unique insight into the upstream genetic control of ADAR activity, which is highly relevant to common inflammatory diseases mentioned in the literature.

Check https://www.youtube.com/watch?v=I1Na06UdW-E

## 2. Rigorous Handling of the AEI (Methodological Refinement)

While not strictly an original idea, your pipeline's design rigorously addresses a methodological challenge explicitly raised in the Li et al. (2022) introduction:Context from Li et al. (2022)Your Pipeline Implementation"...the reduced editing of immunogenic dsRNAs leads to interferon responses, which may subsequently induce expression of ADAR1 and affect the overall editing levels..."Phase 6: AEI as a Covariate. By including the AEI as a covariate in the primary edQTL model, you statistically decouple the specific variant-site association from the global, systemic effects of ADAR induction (like those caused by an IFN response). This ensures your edQTL results are more precise and less confounded by changes in ADAR expression.

## New database repository

https://github.com/TrinityCTAT/ctat-genome-lib-builder/wiki
https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/MUTATION_LIB_SUPPLEMENT/rna_editing/

$ wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/MUTATION_LIB_SUPPLEMENT/rna_editing/GRCh38.RNAediting.vcf.gz
