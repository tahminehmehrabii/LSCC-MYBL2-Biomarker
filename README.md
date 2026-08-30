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

## Repository Structure

```text
LSCC_MYBL2_Reproducible_Code/
├── README.md
├── requirements.R
├── config/
│   └── paths.example.R
└── scripts/
    ├── 01_scRNAseq_high_CNV_marker_analysis.R
    ├── 02_bulk_RNAseq_ML_pipeline.R
    ├── 03_immune_microenvironment_analysis.R
    ├── 04_functional_programme_analysis.R
    ├── 05_MYBL2_continuous_Hallmark_GSEA.R
    ├── 06_cluster_signature_ssGSEA_validation.R
    ├── 07_CellChat_clusters_0_8.R
    ├── 08_prepare_Monocle2_input.R
    ├── 09_Monocle2_pseudotime_analysis.R
    ├── 10_marker_threshold_sensitivity.R
    ├── 11_WGCNA_beta_sensitivity.R
    ├── 12_multiseed_stability_analysis.R
    ├── 13_multigene_panel_benchmark.R
    ├── 14_RF_vs_SVM_DeLong_test.R
    ├── 15_established_biomarker_benchmark.R
    └── 16_purity_adjusted_immune_sensitivity.R
```

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

## Configuration

The original analyses were conducted on Windows using absolute paths under:

```text
D:/LSCC
E:/LSCC
```

Before running the scripts:

1. Copy `config/paths.example.R` to `config/paths.R`.
2. Edit `config/paths.R` to define the project paths for your local system.
3. Update the path block at the beginning of each analysis script to match the configured paths.
4. Confirm that all required input files and intermediate objects are available in the specified directories.

Absolute paths are retained in the scripts to preserve the provenance of the reported analyses.


## Software Requirements

The trajectory workflow was developed using:

- **R version 4.4.3**
- **Monocle version 2.34.0**

Other scripts should be run using the package versions recorded in their generated `sessionInfo.txt` files.

Run the following script to identify missing R packages:

```r
source("requirements.R")
```

The `requirements.R` script reports missing dependencies but does **not** install them automatically.

## Reproducibility

The following measures were implemented to support computational reproducibility:

- Random seeds are explicitly defined for stochastic analyses.
- Model directions derived from the training cohort are applied unchanged to the validation and external cohorts.
- Scaling parameters derived from the training cohort are applied unchanged to the validation and external cohorts.
- Decision thresholds derived from the training cohort are applied unchanged to the validation and external cohorts.
- Validation and external cohorts are not used for model fitting.
- Multiple-testing correction is performed using the **Benjamini–Hochberg procedure** where specified.
- Sensitivity analyses write their results to separate output directories.
- Sensitivity analyses do not overwrite the primary analysis results.

## Data Availability

The scripts require the following input files and objects:

- Gene-expression matrices
- Phenotype tables
- Cell-annotation files
- Gene-order file
- Saved intermediate RDS objects
- Dataset-specific metadata files

The required file locations are described in the path section of each script.

Raw and processed patient-level data are not distributed with this repository.

## Public Datasets

Public datasets used in the workflow are available through the following accession identifiers:

- **GSE127165**
- **GSE142083**
- **GSE130605**
- **GSE206332**

## Recommended Execution Order

For complete reproduction of the analysis, run the scripts in numerical order:

```text
01 → 02 → 03 → 04 → 05 → 06 → 07 → 08 → 09
```

The sensitivity and benchmark analyses can then be run as follows:

```text
10 → 11 → 12 → 13 → 14 → 15 → 16
```

Some scripts require intermediate RDS objects or output files generated by earlier stages of the workflow. Verify all input paths before execution.

## Citation

If you use this code or workflow, please cite the associated article.

Full citation details will be added upon publication.

## License

Please refer to the repository license for the applicable terms of reuse and distribution.



