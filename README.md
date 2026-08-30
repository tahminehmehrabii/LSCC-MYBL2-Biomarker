LSCC MYBL2 Reproducible Analysis

This repository contains the R code used for the reproducible analysis of MYBL2 in lung squamous cell carcinoma (LSCC). The workflow integrates single-cell RNA sequencing, bulk transcriptomics, machine learning, immune-microenvironment profiling, functional-programme analysis, cell–cell communication, trajectory inference, and sensitivity analyses.

Repository Structure

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

Main Analysis Workflow

Run the primary analyses in the following order:

01_scRNAseq_high_CNV_marker_analysis.R
Performs Seurat processing, cell-type annotation, inferCNV analysis, malignant-cell identification, and high-CNV marker detection.

02_bulk_RNAseq_ML_pipeline.R
Performs bulk RNA-seq preprocessing, differential-expression analysis, WGCNA, integration of single-cell, bulk, and WGCNA results, LASSO feature selection, internal ROC analysis, and external machine-learning validation.

03_immune_microenvironment_analysis.R
Characterizes the MYBL2-associated immune microenvironment using ssGSEA and ESTIMATE.

04_functional_programme_analysis.R
Scores functional programmes across high-CNV malignant-cell subclusters.

05_MYBL2_continuous_Hallmark_GSEA.R
Performs Hallmark GSEA using continuous MYBL2 expression in bulk tumor samples.

06_cluster_signature_ssGSEA_validation.R
Validates the single-cell signatures of malignant Clusters 0 and 8 in bulk transcriptomic cohorts.

07_CellChat_clusters_0_8.R
Analyzes outgoing cell–cell signaling from malignant Clusters 0 and 8 using CellChat.

08_prepare_Monocle2_input.R
Prepares balanced Monocle 2 input from the saved high-CNV Seurat objects.

09_Monocle2_pseudotime_analysis.R
Performs Monocle 2 DDRTree trajectory inference and evaluates MYBL2 expression across pseudotime.

Sensitivity and Benchmark Analyses

The following scripts assess the robustness, stability, and comparative performance of the main findings:

10_marker_threshold_sensitivity.R — Evaluates permissive, primary, and stringent single-cell marker definitions, followed by repeated intersections with bulk DEGs and WGCNA modules.

11_WGCNA_beta_sensitivity.R — Assesses WGCNA network stability across alternative soft-thresholding powers.

12_multiseed_stability_analysis.R — Evaluates the stability of feature selection and model performance across multiple random seeds.

13_multigene_panel_benchmark.R — Compares MYBL2 with other LASSO-selected genes and a five-gene logistic model.

14_RF_vs_SVM_DeLong_test.R — Performs a paired DeLong test comparing external random-forest and support-vector-machine AUROCs.

15_established_biomarker_benchmark.R — Benchmarks MYBL2 against the established biomarkers CDKN2A, EGFR, and MKI67.

16_purity_adjusted_immune_sensitivity.R — Performs tumor-purity- and dataset-adjusted immune sensitivity analyses.

Configuration

The original analyses were run on Windows using absolute paths under D:/LSCC and E:/LSCC. Before running the scripts:

Copy config/paths.example.R to config/paths.R.

Edit config/paths.R to define the project paths for your local system.

Update the path block at the beginning of each analysis script to match the configured paths.

Absolute paths are retained in the scripts to preserve the provenance of the reported analyses. No raw or processed patient-level data are included in this repository.

Software Requirements

The trajectory workflow was developed using R 4.4.3 and Monocle 2.34.0.

Other scripts should be run with the package versions recorded in their generated sessionInfo.txt files.

Run requirements.R to identify missing packages. The script reports missing dependencies but does not install them automatically.

Reproducibility

Random seeds are set explicitly for stochastic analyses.

Model directions, scaling parameters, and decision thresholds derived from the training cohort are applied unchanged to validation and external cohorts.

Validation and external cohorts are not used for model fitting.

Multiple-testing correction uses the Benjamini–Hochberg procedure where specified.

Sensitivity analyses write to separate output directories and do not overwrite the primary results.

Data Availability

The scripts require the input matrices, phenotype tables, annotation files, gene-order file, and intermediate RDS objects specified in their path sections. These data are not distributed with this repository.

Public datasets used in the workflow are available under the following accession identifiers:

GSE127165

GSE142083

GSE130605

GSE206332

Citation

If you use this code, please cite the associated article. Full citation details will be added upon publication.

License

Please refer to the repository license for terms of reuse and distribution.



