# LSCC MYBL2 Reproducible Analysis

This repository contains the **R code** used for the reproducible analysis of **MYBL2 in lung squamous cell carcinoma (LSCC)**.

The workflow integrates:

- Single-cell RNA sequencing
- Bulk transcriptomic analysis
- Machine-learning analysis
- Immune-microenvironment profiling
- Functional-programme analysis
- Cell–cell communication analysis
- Trajectory inference
- Sensitivity and benchmark analyses

## Main Analysis Workflow

Run the primary analyses in the following order.

### 1. Single-Cell RNA-Sequencing Analysis

**Script:** `01_scRNAseq_high_CNV_marker_analysis.R`

This script performs:

- Seurat preprocessing
- Cell-type annotation
- inferCNV analysis
- Malignant-cell identification
- High-CNV marker detection

### 2. Bulk RNA-Sequencing and Machine-Learning Analysis

**Script:** `02_bulk_RNAseq_ML_pipeline.R`

This script performs:

- Bulk RNA-seq preprocessing
- Differential-expression analysis
- Weighted gene co-expression network analysis (WGCNA)
- Integration of single-cell, bulk, and WGCNA results
- LASSO feature selection
- Internal receiver operating characteristic (ROC) analysis
- External machine-learning validation

### 3. Immune-Microenvironment Analysis

**Script:** `03_immune_microenvironment_analysis.R`

This script characterizes the **MYBL2-associated immune microenvironment** using:

- Single-sample gene set enrichment analysis (ssGSEA)
- ESTIMATE analysis

### 4. Functional-Programme Analysis

**Script:** `04_functional_programme_analysis.R`

This script calculates functional-programme scores across high-CNV malignant-cell subclusters.

### 5. Continuous MYBL2 Hallmark GSEA

**Script:** `05_MYBL2_continuous_Hallmark_GSEA.R`

This script performs Hallmark gene set enrichment analysis using continuous **MYBL2 expression** in bulk tumor samples.

### 6. Cluster-Signature Validation

**Script:** `06_cluster_signature_ssGSEA_validation.R`

This script validates the single-cell signatures of malignant **Cluster 0** and **Cluster 8** in bulk transcriptomic cohorts using ssGSEA.

### 7. Cell–Cell Communication Analysis

**Script:** `07_CellChat_clusters_0_8.R`

This script uses CellChat to analyze outgoing cell–cell signaling from malignant **Clusters 0 and 8**.

### 8. Monocle 2 Input Preparation

**Script:** `08_prepare_Monocle2_input.R`

This script prepares balanced Monocle 2 input from the saved high-CNV Seurat objects.

### 9. Pseudotime Analysis

**Script:** `09_Monocle2_pseudotime_analysis.R`

This script performs:

- Monocle 2 DDRTree trajectory inference
- Pseudotime ordering
- Evaluation of MYBL2 expression across pseudotime

## Sensitivity and Benchmark Analyses

The following scripts assess the robustness, stability, and comparative performance of the main findings.

### 10. Marker-Threshold Sensitivity Analysis

**Script:** `10_marker_threshold_sensitivity.R`

This script evaluates permissive, primary, and stringent single-cell marker definitions, followed by repeated intersections with bulk differentially expressed genes and WGCNA modules.

### 11. WGCNA Soft-Threshold Sensitivity Analysis

**Script:** `11_WGCNA_beta_sensitivity.R`

This script assesses WGCNA network stability across alternative soft-thresholding powers.

### 12. Multi-Seed Stability Analysis

**Script:** `12_multiseed_stability_analysis.R`

This script evaluates the stability of feature selection and model performance across multiple random seeds.

### 13. Multigene-Panel Benchmark

**Script:** `13_multigene_panel_benchmark.R`

This script compares:

- MYBL2
- Other LASSO-selected genes
- A five-gene logistic-regression model

### 14. Random Forest versus SVM Comparison

**Script:** `14_RF_vs_SVM_DeLong_test.R`

This script performs a paired DeLong test comparing the external area under the ROC curve values of the random forest and support vector machine models.

### 15. Established-Biomarker Benchmark

**Script:** `15_established_biomarker_benchmark.R`

This script benchmarks **MYBL2** against the following established biomarkers:

- **CDKN2A**
- **EGFR**
- **MKI67**

### 16. Purity-Adjusted Immune Sensitivity Analysis

**Script:** `16_purity_adjusted_immune_sensitivity.R`

This script performs tumor-purity- and dataset-adjusted immune sensitivity analyses.

