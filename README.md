# LSCC-MYBL2-Biomarkers

## Overview

This github repository contains the data files and analysis code used for the scientific paper titled "Integrative bulk and single-cell transcriptomic analysis identifies MYBL2-associated molecular and immune features in laryngeal squamous cell carcinoma". The files are organised into three folders:

Data: which contains all the transcriptomic data required to perform the analyses.

Codes: contains the R scripts to reproduce all analyses.

Results: contains all the results produced by the R scripts.

## Reproducing the results

This repository contains all the code necessary to reproduce the results.

First, download the repository and place it in your project directory.

```bash
git clone https://github.com/your-username/LSCC-MYBL2-Biomarkers.git path/to/directory
```

In this command, "path/to/directory" refers to your project path.

Then, run the following commands in order.

Before running the code, make sure to set the project path using the setwd() command at the beginning of each code.

```r
setwd(project_path)
```

1. Run scRNAseqAnalysis.R to process the single-cell RNA-seq dataset, perform quality control, normalization, dimensionality reduction, Harmony integration, clustering, cell-type annotation, malignant epithelial cell identification, and inferCNV analysis.

2. Run bulkRNAseqML.R to process the bulk RNA-seq datasets, prepare expression matrices, perform differential expression analysis, identify MYBL2-associated genes, and evaluate diagnostic performance using machine learning models.

3. Run immuneMicroenvironmentAnalysis.R to assess immune cell infiltration in LSCC and evaluate the association between MYBL2 expression and immune microenvironment features.

4. Run functionalProgrammeAnalysis.R to investigate functional programme activity in High-CNV malignant epithelial subclusters and compare MYBL2-related biological programmes across malignant subtypes.

5. Run pseudotimeInputPreparation.R to prepare the required High-CNV malignant epithelial cell input object for Monocle 2 pseudotime analysis.

6. Run pseudotimeMonocleAnalysis.R to perform Monocle 2 DDRTree pseudotime analysis and evaluate MYBL2 expression dynamics along the malignant-cell trajectory.

7. Run clusterSignatureValidation.R to construct scRNA-seq-derived signatures for High-CNV malignant Cluster 0 and Cluster 8 and validate these signatures in bulk LSCC cohorts using ssGSEA.

8. Run MYBL2HallmarkGSEA.R to perform MYBL2-continuous Hallmark gene set enrichment analysis in bulk tumor samples.

9. Run CellChatAnalysis.R to assess ligand-receptor communication from High-CNV malignant Cluster 0 and Cluster 8 to selected tumor microenvironment cell populations.

## Required software

The scripts use core R functionality and several publicly available R packages listed below. Version numbers in brackets correspond to the approximate package versions that were used to develop, test, and debug these scripts.

R (4.4.3)

RStudio (2024.12.0+467): Optional, testing functions and running the code step-by-step.

Seurat (5.1.0)

SeuratObject (5.0.2)

Matrix (1.7.0)

data.table (1.16.0)

dplyr (1.1.4)

tidyr (1.3.1)

tibble (3.2.1)

stringr (1.5.1)

ggplot2 (3.5.1)

patchwork (1.2.0)

pheatmap (1.0.12)

ggrepel (0.9.5)

harmony (1.2.0)

infercnv (1.20.0)

limma (3.60.0)

edgeR (4.2.0)

WGCNA (1.72.5)

impute (1.78.0)

glmnet (4.1.8)

pROC (1.18.5)

randomForest (4.7.1.1)

e1071 (1.7.14)

clusterProfiler (4.12.0)

msigdbr (7.5.1)

fgsea (1.30.0)

GSVA (1.52.0)

IOBR (0.99.9)

readxl (1.4.3)

CellChat (2.1.2)

igraph (2.0.3)

Biobase (2.64.0)

Monocle (2.34.0)

DDRTree (0.1.5)

VGAM (1.1.11)

scales (1.3.0)

cowplot (1.1.3)

magick (2.8.5)

png (0.1.8)

grid (4.4.3)
