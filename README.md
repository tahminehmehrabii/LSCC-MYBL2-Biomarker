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

The scripts use core R functionality and several publicly available R packages listed below. 
limma (3.62.2)

WGCNA (1.74)

Seurat (5.5.0)

harmony (1.2.4)

glmnet (4.1.10)

pROC (1.19.0.1)

e1071 (1.7-17)

randomForest (4.7-1.2)

GSVA (2.0.7)

IOBR (2.2.3)

fgsea (1.32.4)

Monocle 2 (2.34.0)

CellChat (2.2.0.9001)
