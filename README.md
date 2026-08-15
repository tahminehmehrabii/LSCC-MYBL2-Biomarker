LSCC MYBL2 reproducible analysis

This repository contains the R code used for the single-cell, bulk-transcriptomic, machine-learning, immune, functional, trajectory, and reviewer-requested sensitivity analyses of MYBL2 in lung squamous cell carcinoma (LSCC).

Repository structure

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

Main workflow

Run the main analyses in the following order:

01_scRNAseq_high_CNV_marker_analysis.R: Seurat processing, cell annotation, inferCNV, malignant-cell identification, and high-CNV marker detection.

02_bulk_RNAseq_ML_pipeline.R: bulk RNA-seq preprocessing, differential expression, WGCNA, single-cell/bulk/WGCNA integration, LASSO selection, internal ROC analysis, and external machine-learning validation.

03_immune_microenvironment_analysis.R: ssGSEA and ESTIMATE analyses of the MYBL2-associated immune microenvironment.

04_functional_programme_analysis.R: functional-programme scoring of high-CNV malignant subclusters.

05_MYBL2_continuous_Hallmark_GSEA.R: MYBL2-continuous Hallmark GSEA in bulk tumor samples.

06_cluster_signature_ssGSEA_validation.R: bulk validation of Cluster 0 and Cluster 8 single-cell signatures.

07_CellChat_clusters_0_8.R: CellChat analysis of outgoing signaling from malignant Clusters 0 and 8.

08_prepare_Monocle2_input.R: preparation of balanced Monocle 2 input from the saved high-CNV Seurat objects.

09_Monocle2_pseudotime_analysis.R: Monocle 2 DDRTree trajectory and MYBL2 pseudotime analysis.

Reviewer-requested analyses

10_marker_threshold_sensitivity.R: permissive, primary, and stringent single-cell marker definitions followed by repeated bulk-DEG and WGCNA intersections.

11_WGCNA_beta_sensitivity.R: WGCNA network stability across alternative soft-threshold powers.

12_multiseed_stability_analysis.R: stability of feature selection and model performance across sampling seeds.

13_multigene_panel_benchmark.R: comparison of MYBL2, other LASSO-selected genes, and a five-gene logistic model.

14_RF_vs_SVM_DeLong_test.R: paired DeLong comparison of external RF and SVM AUROCs.

15_established_biomarker_benchmark.R: comparison of MYBL2 with CDKN2A, EGFR, and MKI67.

16_purity_adjusted_immune_sensitivity.R: tumor-purity and dataset-adjusted immune sensitivity analyses.

Configuration

The original analyses were run on Windows and use absolute paths under D:/LSCC or E:/LSCC. Before running the scripts:

Copy config/paths.example.R to config/paths.R.

Edit the project paths for the local system.

Update the path block at the beginning of each analysis script to match config/paths.R.

Absolute paths are retained in the scripts to preserve the exact provenance of the reported analysis. No raw or processed patient-level data are included in this code repository.

Software

The trajectory workflow was developed for R 4.4.3 and Monocle 2.34.0. Other scripts should be run with the package versions reported by their generated sessionInfo.txt files. Run requirements.R to identify missing packages; it does not install packages automatically.

Reproducibility notes

Random seeds are set explicitly in stochastic analyses.

Training-derived model directions, scaling parameters, and decision thresholds are applied unchanged to validation and external cohorts.

Validation and external cohorts are not used for model fitting.

Multiple-testing correction is performed with the Benjamini–Hochberg procedure where specified.

Reviewer sensitivity analyses write to separate output directories and do not overwrite the primary analysis.

Data availability

The scripts expect the input matrices, phenotype tables, annotation files, gene-order file, and saved intermediate RDS objects described in their path sections. Public accession identifiers used in the workflow include GSE127165, GSE142083, GSE130605, and GSE206332.

Citation

If the code is released with the associated article, cite the final published article and the repository release or DOI.
