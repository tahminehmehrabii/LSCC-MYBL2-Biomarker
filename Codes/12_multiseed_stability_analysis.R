###############################################################################
# MULTI-SEED STABILITY ANALYSIS
# LSCC MYBL2 reproducible analysis
#
# This script runs ONLY the random-seed stability analysis.
# It assumes the main preprocessing pipeline has already generated:
#   - train_discovery_merged_CPM_nonnegative_noZero.csv
#   - GSE130605_meta_CPM.csv
# and that the precomputed scRNA High-CNV marker file exists.
#
# Only the sample-selection seed changes across runs. Randomness in LASSO CV,
# WGCNA, and random forest fitting is held fixed at seed 123.
###############################################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
options(scipen = 100)
options(timeout = 3600)

# -----------------------------------------------------------------------------
# 01. PATHS
# -----------------------------------------------------------------------------

results_path <- "E:/LSCC/Results_LSCC/ML"
scrna_figdir <- "E:/LSCC/ScRNAseq_Results/GSE206332/Results/figures"
TARGET_GENE <- "MYBL2"

dir.create(results_path, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 02. PACKAGES
# -----------------------------------------------------------------------------

required_pkgs <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "limma", "WGCNA",
  "glmnet", "pROC", "randomForest", "e1071"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nPlease install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(limma)
  library(WGCNA)
  library(glmnet)
  library(pROC)
  library(randomForest)
  library(e1071)
})

# -----------------------------------------------------------------------------
# 03. SEED MODES
# -----------------------------------------------------------------------------
# Choose ONE mode below.
# QUICK_3   : fastest 3-seed sensitivity check.
# STANDARD_4: recommended balance between robustness and runtime.
# EXTENDED_9: strongest sensitivity check here, with WGCNA run 9 times.
#
# Seed 42 has been intentionally excluded from all sensitivity-analysis modes.
# All values are ordinary integer random seeds; they are not different
# statistical methods. They simply generate alternative random sample selections
# and 70:30 train-validation partitions.

SEED_MODE <- "STANDARD_4"

SEED_SETS <- list(
  QUICK_3 = c(1, 123, 999),
  STANDARD_4 = c(1, 123, 456, 999),
  EXTENDED_9 = c(1, 7, 21, 77, 123, 321, 456, 777, 999)
)

if (!SEED_MODE %in% names(SEED_SETS)) {
  stop(
    "Invalid SEED_MODE. Choose one of: ",
    paste(names(SEED_SETS), collapse = ", ")
  )
}

STABILITY_SEEDS <- SEED_SETS[[SEED_MODE]]
STABILITY_MODEL_SEED <- 123
STABILITY_TARGET_GENE <- TARGET_GENE

cat("\nSelected seed mode:", SEED_MODE, "\n")
cat("Sampling seeds:", paste(STABILITY_SEEDS, collapse = ", "), "\n")
cat("Fixed model-fitting seed:", STABILITY_MODEL_SEED, "\n")

# -----------------------------------------------------------------------------
# 04. ANALYSIS THRESHOLDS
# -----------------------------------------------------------------------------

DEG_FDR_THRESH <- 0.01
DEG_LOGFC_MIN <- 2
WGCNA_TOP_MAD <- 8000
WGCNA_TARGET_R2 <- 0.92
WGCNA_MIN_MEANK <- 10
WGCNA_POWERS <- 1:20
WGCNA_MM_MIN <- 0.55
WGCNA_GS_MIN <- 0.55
WGCNA_RULE <- "AND"
LASSO_NFOLDS <- 10
ROC_AUC_MIN <- 0.95
RF_NTREE <- 500
SVM_KERNEL <- "linear"

# -----------------------------------------------------------------------------
# 05. FIGURE SETTINGS
# -----------------------------------------------------------------------------

FONT_FAMILY <- "Arial"
FIG_DPI <- 600
FIG_BACKGROUND <- "white"
COL_NORMAL <- "#7470B2"
COL_TUMOR <- "#D81B60"

if (.Platform$OS.type == "windows") {
  try(
    grDevices::windowsFonts(
      Arial = grDevices::windowsFont("Arial")
    ),
    silent = TRUE
  )
}

theme_manuscript <- function(show_grid = FALSE,
                             legend_position = "right") {
  ggplot2::theme_bw(base_size = 10, base_family = FONT_FAMILY) +
    ggplot2::theme(
      text = ggplot2::element_text(family = FONT_FAMILY, color = "black"),
      plot.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = 11,
        hjust = 0.5,
        color = "black"
      ),
      axis.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = 10,
        color = "black"
      ),
      axis.text = ggplot2::element_text(
        family = FONT_FAMILY,
        size = 9,
        color = "black"
      ),
      legend.position = legend_position,
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.45),
      panel.grid.major = if (show_grid) ggplot2::element_line(color = "grey92") else ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

save_plot_all_formats <- function(plot_obj,
                                  filename_stem,
                                  dir_path,
                                  width = 8.8,
                                  height = 5.8) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  ggplot2::ggsave(
    file.path(dir_path, paste0(filename_stem, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = FIG_DPI,
    bg = FIG_BACKGROUND,
    limitsize = FALSE
  )
  
  ggplot2::ggsave(
    file.path(dir_path, paste0(filename_stem, ".tiff")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = FIG_DPI,
    compression = "lzw",
    bg = FIG_BACKGROUND,
    limitsize = FALSE
  )
  
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  
  ggplot2::ggsave(
    file.path(dir_path, paste0(filename_stem, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height,
    device = pdf_device,
    bg = FIG_BACKGROUND,
    limitsize = FALSE
  )
}

# -----------------------------------------------------------------------------
# 06. BASIC HELPERS
# -----------------------------------------------------------------------------

to_numeric_df <- function(df) {
  as.data.frame(
    lapply(df, function(x) as.numeric(as.character(x))),
    check.names = FALSE
  )
}

clean_gene_vector <- function(x) {
  x <- unique(trimws(as.character(x)))
  x[!is.na(x) & x != ""]
}

clean_id <- function(x) trimws(as.character(x))

coerce_group_12 <- function(g) {
  g0 <- as.character(g)
  
  if (all(g0 %in% c("1", "2"))) {
    return(as.integer(g0))
  }
  
  g_low <- tolower(g0)
  
  if (all(g_low %in% c("normal", "tumor", "non", "lscc", "cancer", "margin"))) {
    return(ifelse(g_low %in% c("tumor", "lscc", "cancer"), 2L, 1L))
  }
  
  u <- sort(unique(g0))
  
  if (length(u) != 2) {
    stop("group must have exactly 2 classes. Found: ", paste(u, collapse = ", "))
  }
  
  map <- setNames(c(1L, 2L), u)
  as.integer(map[g0])
}

make_group_factor <- function(g) {
  g12 <- coerce_group_12(g)
  factor(g12, levels = c(1, 2), labels = c("Normal", "Tumor"))
}

safe_read_gene_list <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(character(0))
  }
  
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    dat <- as.data.frame(data.table::fread(path), check.names = FALSE)
    gene_column <- intersect(
      c("Gene", "Genes", "gene", "Gene.symbol", "Symbol", "symbol"),
      colnames(dat)
    )
    
    if (length(gene_column) < 1) {
      return(character(0))
    }
    
    genes <- dat[[gene_column[1]]]
  } else {
    genes <- readLines(path, warn = FALSE)
  }
  
  clean_gene_vector(genes)
}

# -----------------------------------------------------------------------------
# 07. MULTI-SEED ANALYSIS
# -----------------------------------------------------------------------------
stability_dir <- file.path(results_path, paste0("Seed_Stability_Analysis_", SEED_MODE))
dir.create(stability_dir, recursive = TRUE, showWarnings = FALSE)

# Use the same thresholds as the main analysis pipeline.
if (!exists("DEG_FDR_THRESH", inherits = FALSE)) DEG_FDR_THRESH <- 0.01
if (!exists("DEG_LOGFC_MIN", inherits = FALSE))  DEG_LOGFC_MIN  <- 2
if (!exists("WGCNA_TOP_MAD", inherits = FALSE)) WGCNA_TOP_MAD <- 8000
if (!exists("WGCNA_TARGET_R2", inherits = FALSE)) WGCNA_TARGET_R2 <- 0.92
if (!exists("WGCNA_MIN_MEANK", inherits = FALSE)) WGCNA_MIN_MEANK <- 10
if (!exists("WGCNA_POWERS", inherits = FALSE)) WGCNA_POWERS <- 1:20
if (!exists("WGCNA_MM_MIN", inherits = FALSE)) WGCNA_MM_MIN <- 0.55
if (!exists("WGCNA_GS_MIN", inherits = FALSE)) WGCNA_GS_MIN <- 0.55
if (!exists("WGCNA_RULE", inherits = FALSE)) WGCNA_RULE <- "AND"
if (!exists("LASSO_NFOLDS", inherits = FALSE)) LASSO_NFOLDS <- 10
if (!exists("ROC_AUC_MIN", inherits = FALSE)) ROC_AUC_MIN <- 0.95
if (!exists("RF_NTREE", inherits = FALSE)) RF_NTREE <- 500
if (!exists("SVM_KERNEL", inherits = FALSE)) SVM_KERNEL <- "linear"

# ---------------------------------------------------------------------------
# 08. REQUIRED INPUTS
# ---------------------------------------------------------------------------

# The main pipeline creates this file before sample splitting.
if (!exists("s2_merged_nozero_csv", inherits = FALSE)) {
  s2_merged_nozero_csv <- file.path(
    results_path,
    "train_discovery_merged_CPM_nonnegative_noZero.csv"
  )
}

# External CPM with metadata created in Step 02.
if (!exists("s2_te_meta", inherits = FALSE)) {
  s2_te_meta <- file.path(results_path, "GSE130605_meta_CPM.csv")
}

required_stability_files <- c(
  s2_merged_nozero_csv,
  s2_te_meta
)

missing_stability_files <- required_stability_files[
  !file.exists(required_stability_files)
]

if (length(missing_stability_files) > 0) {
  stop(
    "This analysis requires outputs from script 02. Missing file(s):\n",
    paste(missing_stability_files, collapse = "\n")
  )
}

# Load the fixed scRNA-seq high-CNV malignant marker list.
if (!exists("sc_marker_file_candidates", inherits = FALSE)) {
  sc_marker_file_candidates <- c(
    file.path(scrna_figdir, "High_CNV_Malignant_marker_genes_after_lowCNV_removed_logFC1.txt"),
    file.path(scrna_figdir, "Final_High_CNV_Malignant_Gene_Names.txt"),
    file.path(scrna_figdir, "Final_High_CNV_Malignant_Gene_Names.csv"),
    file.path(scrna_figdir, "High_CNV_Malignant_FindAllMarkers_after_lowCNV_removed_logFC1.csv"),
    file.path(scrna_figdir, "Final_High_CNV_Malignant_Genes.csv"),
    file.path(scrna_figdir, "High_CNV_Malignant_specific_markers_FINAL.csv")
  )
}

sc_marker_file_stability <- sc_marker_file_candidates[
  file.exists(sc_marker_file_candidates)
][1]

if (is.na(sc_marker_file_stability) || !file.exists(sc_marker_file_stability)) {
  stop("Could not find the precomputed scRNA-seq high-CNV marker file.")
}

SC_MARKERS_STABILITY <- safe_read_gene_list(sc_marker_file_stability)

if (length(SC_MARKERS_STABILITY) < 1) {
  stop("The scRNA-seq marker list is empty.")
}

# ---------------------------------------------------------------------------
# 09. LOAD DEVELOPMENT AND EXTERNAL COHORTS
# ---------------------------------------------------------------------------

stability_samples <- as.data.frame(
  data.table::fread(s2_merged_nozero_csv),
  check.names = FALSE
)

if (!all(c("Sample", "group", "batch") %in% colnames(stability_samples))) {
  stop("The pre-split development cohort must contain Sample, group, and batch.")
}

stability_samples$Sample <- clean_id(stability_samples$Sample)
stability_samples$group <- coerce_group_12(stability_samples$group)
rownames(stability_samples) <- stability_samples$Sample

stability_gene_cols <- setdiff(
  colnames(stability_samples),
  c("Sample", "group", "batch")
)

external_stability_cpm <- as.data.frame(
  data.table::fread(s2_te_meta),
  check.names = FALSE
)

if (!all(c("Sample", "group", "batch") %in% colnames(external_stability_cpm))) {
  stop("The external CPM file must contain Sample, group, and batch.")
}

external_stability_cpm$Sample <- clean_id(external_stability_cpm$Sample)
external_stability_cpm$group <- coerce_group_12(external_stability_cpm$group)

common_stability_genes <- intersect(
  stability_gene_cols,
  setdiff(colnames(external_stability_cpm), c("Sample", "group", "batch"))
)

if (length(common_stability_genes) < 2) {
  stop("Too few common genes between the development and external cohorts.")
}

# Restrict once to the common gene universe used in every seed run.
stability_samples <- stability_samples[
  , c("Sample", "group", "batch", common_stability_genes),
  drop = FALSE
]

external_stability_cpm <- external_stability_cpm[
  , c("Sample", "group", "batch", common_stability_genes),
  drop = FALSE
]

# Sample counts before any balancing.
all_tumor_ids <- stability_samples$Sample[stability_samples$group == 2]
all_normal_ids <- stability_samples$Sample[stability_samples$group == 1]
n_balance_each <- min(length(all_tumor_ids), length(all_normal_ids))
n_train_each <- ceiling(0.7 * n_balance_each)

cat("\nPre-balance development cohort:\n")
cat("Tumor =", length(all_tumor_ids), "\n")
cat("Normal =", length(all_normal_ids), "\n")
cat("Balanced size per class =", n_balance_each, "\n")
cat("Training size per class =", n_train_each, "\n")
cat("Validation size per class =", n_balance_each - n_train_each, "\n")
cat("Sampling seeds =", paste(STABILITY_SEEDS, collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# 10. ANALYSIS HELPERS
# ---------------------------------------------------------------------------

stability_prepare_log2 <- function(df) {
  genes <- setdiff(colnames(df), c("Sample", "group", "batch"))
  M <- as.matrix(to_numeric_df(df[, genes, drop = FALSE]))
  rownames(M) <- df$Sample
  log2(M + 1)
}

stability_gene_auc <- function(expr, group) {
  group <- factor(group, levels = c("Normal", "Tumor"))
  expr <- as.numeric(expr)
  ok <- is.finite(expr) & !is.na(group)
  
  if (sum(ok) < 3 || length(unique(group[ok])) < 2) {
    return(NA_real_)
  }
  
  rr <- pROC::roc(
    response = group[ok],
    predictor = expr[ok],
    levels = c("Normal", "Tumor"),
    direction = "<",
    quiet = TRUE
  )
  
  as.numeric(pROC::auc(rr))
}

stability_scale_by_train <- function(train_mat, new_mat) {
  train_mat <- as.matrix(train_mat)
  new_mat <- as.matrix(new_mat)
  
  mu <- colMeans(train_mat, na.rm = TRUE)
  sdv <- apply(train_mat, 2, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  
  sweep(sweep(new_mat, 2, mu, "-"), 2, sdv, "/")
}

# Robust helper for per-seed gene records.
# If a stage returns zero genes, return a valid zero-row data frame instead of
# combining a length-1 Seed/Stage column with a length-0 Gene column.
stability_make_gene_record <- function(seed_value, stage_name, genes) {
  genes <- as.character(genes)
  genes <- genes[!is.na(genes) & genes != ""]
  
  if (length(genes) == 0) {
    return(data.frame(
      Seed = integer(0),
      Stage = character(0),
      Gene = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    Seed = rep(seed_value, length(genes)),
    Stage = rep(stage_name, length(genes)),
    Gene = genes,
    stringsAsFactors = FALSE
  )
}

stability_external_ml <- function(train_df_seed,
                                  external_df_seed,
                                  model_genes_seed) {
  
  model_genes_seed <- intersect(
    model_genes_seed,
    intersect(
      setdiff(colnames(train_df_seed), c("Sample", "group", "batch")),
      setdiff(colnames(external_df_seed), c("Sample", "group", "batch"))
    )
  )
  
  if (length(model_genes_seed) < 1) {
    return(c(RF_AUROC = NA_real_, SVM_AUROC = NA_real_))
  }
  
  trainX_seed <- as.matrix(
    to_numeric_df(train_df_seed[, model_genes_seed, drop = FALSE])
  )
  extX_seed <- as.matrix(
    to_numeric_df(external_df_seed[, model_genes_seed, drop = FALSE])
  )
  
  trainY_seed <- make_group_factor(train_df_seed$group)
  extY_seed <- make_group_factor(external_df_seed$group)
  
  train_scaled_seed <- stability_scale_by_train(trainX_seed, trainX_seed)
  ext_scaled_seed <- stability_scale_by_train(trainX_seed, extX_seed)
  
  train_ml_seed <- data.frame(
    group = trainY_seed,
    train_scaled_seed,
    check.names = FALSE
  )
  
  ext_ml_seed <- data.frame(
    group = extY_seed,
    ext_scaled_seed,
    check.names = FALSE
  )
  
  colnames(train_ml_seed) <- make.names(colnames(train_ml_seed), unique = TRUE)
  colnames(ext_ml_seed) <- make.names(colnames(ext_ml_seed), unique = TRUE)
  
  # Keep model randomness fixed so only sampling seed changes across runs.
  set.seed(STABILITY_MODEL_SEED)
  rf_seed <- randomForest::randomForest(
    group ~ .,
    data = train_ml_seed,
    ntree = RF_NTREE,
    importance = FALSE
  )
  
  rf_prob_seed <- predict(
    rf_seed,
    newdata = ext_ml_seed,
    type = "prob"
  )[, "Tumor"]
  
  rf_roc_seed <- pROC::roc(
    extY_seed,
    rf_prob_seed,
    levels = c("Normal", "Tumor"),
    direction = "<",
    quiet = TRUE
  )
  
  set.seed(STABILITY_MODEL_SEED)
  svm_seed <- e1071::svm(
    group ~ .,
    data = train_ml_seed,
    kernel = SVM_KERNEL,
    probability = TRUE
  )
  
  svm_pred_seed <- predict(
    svm_seed,
    newdata = ext_ml_seed,
    probability = TRUE
  )
  
  svm_prob_mat_seed <- attr(svm_pred_seed, "probabilities")
  
  svm_prob_seed <- if ("Tumor" %in% colnames(svm_prob_mat_seed)) {
    svm_prob_mat_seed[, "Tumor"]
  } else {
    1 - svm_prob_mat_seed[, "Normal"]
  }
  
  svm_roc_seed <- pROC::roc(
    extY_seed,
    svm_prob_seed,
    levels = c("Normal", "Tumor"),
    direction = "<",
    quiet = TRUE
  )
  
  c(
    RF_AUROC = as.numeric(pROC::auc(rf_roc_seed)),
    SVM_AUROC = as.numeric(pROC::auc(svm_roc_seed))
  )
}

# ---------------------------------------------------------------------------
# 11. SINGLE-SEED ANALYSIS
# ---------------------------------------------------------------------------

run_stability_seed <- function(seed_value) {
  
  cat("\n------------------------------------------------------------\n")
  cat("Sampling seed:", seed_value, "\n")
  cat("------------------------------------------------------------\n")
  
  tryCatch({
    
    # -----------------------------------------------------------------------
    # A. BALANCING + TRAIN/VALIDATION SPLIT
    # -----------------------------------------------------------------------
    
    set.seed(seed_value)
    
    tumor_bal_seed <- sample(all_tumor_ids, n_balance_each)
    normal_bal_seed <- sample(all_normal_ids, n_balance_each)
    
    training_tumor_seed <- sample(tumor_bal_seed, n_train_each)
    validation_tumor_seed <- setdiff(tumor_bal_seed, training_tumor_seed)
    
    training_normal_seed <- sample(normal_bal_seed, n_train_each)
    validation_normal_seed <- setdiff(normal_bal_seed, training_normal_seed)
    
    training_ids_seed <- c(training_tumor_seed, training_normal_seed)
    validation_ids_seed <- c(validation_tumor_seed, validation_normal_seed)
    
    excluded_balance_ids <- c(
      setdiff(all_tumor_ids, tumor_bal_seed),
      setdiff(all_normal_ids, normal_bal_seed)
    )
    
    train_seed_cpm <- stability_samples[
      match(training_ids_seed, stability_samples$Sample),
      , drop = FALSE
    ]
    
    valid_seed_cpm <- stability_samples[
      match(validation_ids_seed, stability_samples$Sample),
      , drop = FALSE
    ]
    
    # -----------------------------------------------------------------------
    # B. LOG2-CPM
    # -----------------------------------------------------------------------
    
    train_seed_log2 <- stability_prepare_log2(train_seed_cpm)
    valid_seed_log2 <- stability_prepare_log2(valid_seed_cpm)
    ext_seed_log2 <- stability_prepare_log2(external_stability_cpm)
    
    train_seed_df <- data.frame(
      Sample = rownames(train_seed_log2),
      group = train_seed_cpm$group,
      batch = train_seed_cpm$batch,
      as.data.frame(train_seed_log2, check.names = FALSE),
      check.names = FALSE
    )
    
    valid_seed_df <- data.frame(
      Sample = rownames(valid_seed_log2),
      group = valid_seed_cpm$group,
      batch = valid_seed_cpm$batch,
      as.data.frame(valid_seed_log2, check.names = FALSE),
      check.names = FALSE
    )
    
    ext_seed_df <- data.frame(
      Sample = rownames(ext_seed_log2),
      group = external_stability_cpm$group,
      batch = external_stability_cpm$batch,
      as.data.frame(ext_seed_log2, check.names = FALSE),
      check.names = FALSE
    )
    
    # -----------------------------------------------------------------------
    # C. DIFFERENTIAL EXPRESSION ON TRAINING SET
    # -----------------------------------------------------------------------
    
    seed_genes_bulk <- setdiff(
      colnames(train_seed_df),
      c("Sample", "group", "batch")
    )
    
    E_seed <- as.matrix(
      to_numeric_df(train_seed_df[, seed_genes_bulk, drop = FALSE])
    )
    rownames(E_seed) <- train_seed_df$Sample
    
    grp_seed <- factor(
      train_seed_df$group,
      levels = c(1, 2),
      labels = c("Normal", "Tumor")
    )
    grp_seed <- relevel(grp_seed, ref = "Normal")
    
    design_seed <- model.matrix(~ grp_seed)
    fit_seed <- limma::lmFit(t(E_seed), design_seed)
    fit_seed <- limma::eBayes(fit_seed, trend = TRUE)
    
    TT_seed <- limma::topTable(
      fit_seed,
      coef = 2,
      number = Inf,
      adjust.method = "BH",
      sort.by = "B"
    )
    
    TT_seed$Gene.symbol <- rownames(TT_seed)
    
    DEG_seed <- TT_seed[
      is.finite(TT_seed$adj.P.Val) &
        TT_seed$adj.P.Val < DEG_FDR_THRESH &
        abs(TT_seed$logFC) > DEG_LOGFC_MIN,
      , drop = FALSE
    ]
    
    Bulk_UP_seed <- unique(
      DEG_seed$Gene.symbol[DEG_seed$logFC > DEG_LOGFC_MIN]
    )
    
    # -----------------------------------------------------------------------
    # D. WGCNA ON THE SAME SEED-SPECIFIC TRAINING SET
    # -----------------------------------------------------------------------
    
    trait_seed <- data.frame(
      Tumor = as.integer(train_seed_df$group == 2),
      Normal = as.integer(train_seed_df$group == 1),
      row.names = train_seed_df$Sample
    )
    
    datExpr0_seed <- as.data.frame(
      lapply(
        train_seed_df[, seed_genes_bulk, drop = FALSE],
        function(x) as.numeric(as.character(x))
      ),
      check.names = FALSE
    )
    rownames(datExpr0_seed) <- train_seed_df$Sample
    
    gsg_seed <- WGCNA::goodSamplesGenes(datExpr0_seed, verbose = 0)
    
    if (!gsg_seed$allOK) {
      datExpr0_seed <- datExpr0_seed[
        gsg_seed$goodSamples,
        gsg_seed$goodGenes,
        drop = FALSE
      ]
      trait_seed <- trait_seed[rownames(datExpr0_seed), , drop = FALSE]
    }
    
    mad_seed <- apply(datExpr0_seed, 2, mad, na.rm = TRUE)
    mad_seed <- mad_seed[is.finite(mad_seed) & !is.na(mad_seed) & mad_seed > 0]
    
    topN_seed <- min(WGCNA_TOP_MAD, length(mad_seed))
    
    if (topN_seed < 2) {
      stop("Too few genes available for WGCNA in seed ", seed_value)
    }
    
    keep_wgcna_seed <- names(
      sort(mad_seed, decreasing = TRUE)
    )[seq_len(topN_seed)]
    
    datExpr_seed <- datExpr0_seed[, keep_wgcna_seed, drop = FALSE]
    
    sft_seed <- WGCNA::pickSoftThreshold(
      datExpr_seed,
      powerVector = WGCNA_POWERS,
      networkType = "signed",
      corFnc = "cor",
      verbose = 0
    )
    
    fitR2_seed <- -sign(sft_seed$fitIndices[, 3]) * sft_seed$fitIndices[, 2]
    meanK_seed <- sft_seed$fitIndices[, 5]
    
    eligible_seed <- which(
      fitR2_seed >= WGCNA_TARGET_R2 &
        meanK_seed >= WGCNA_MIN_MEANK
    )
    
    softPower_seed <- if (length(eligible_seed) > 0) {
      sft_seed$fitIndices[min(eligible_seed), 1]
    } else {
      sft_seed$fitIndices[which.max(fitR2_seed), 1]
    }
    
    # Keep any internal WGCNA randomness fixed across sampling seeds.
    set.seed(STABILITY_MODEL_SEED)
    
    net_seed <- WGCNA::blockwiseModules(
      datExpr_seed,
      power = softPower_seed,
      networkType = "signed",
      TOMType = "signed",
      corType = "pearson",
      minModuleSize = 60,
      mergeCutHeight = 0.18,
      numericLabels = FALSE,
      pamRespectsDendro = TRUE,
      maxBlockSize = ncol(datExpr_seed),
      verbose = 0
    )
    
    module_colors_seed <- net_seed$colors
    names(module_colors_seed) <- colnames(datExpr_seed)
    
    MEs_seed <- WGCNA::orderMEs(net_seed$MEs)
    
    modTraitCor_seed <- cor(
      MEs_seed,
      trait_seed,
      use = "p",
      method = "pearson"
    )
    
    valid_me_seed <- rownames(modTraitCor_seed)
    valid_me_seed <- valid_me_seed[valid_me_seed != "MEgrey"]
    
    if (length(valid_me_seed) < 1) {
      stop("No non-grey WGCNA module for seed ", seed_value)
    }
    
    tumor_cor_seed <- modTraitCor_seed[valid_me_seed, "Tumor"]
    tumor_cor_seed[tumor_cor_seed <= 0] <- NA_real_
    
    if (all(is.na(tumor_cor_seed))) {
      stop("No positively tumor-associated WGCNA module for seed ", seed_value)
    }
    
    bestME_seed <- valid_me_seed[which.max(tumor_cor_seed)]
    targetModule_seed <- sub("^ME", "", bestME_seed)
    
    targetGenes_seed <- names(module_colors_seed)[
      module_colors_seed == targetModule_seed
    ]
    
    GS_seed <- as.numeric(
      cor(datExpr_seed, trait_seed$Tumor, use = "p", method = "pearson")
    )
    
    MM_seed <- as.numeric(
      cor(datExpr_seed, MEs_seed[, bestME_seed], use = "p", method = "pearson")
    )
    
    names(GS_seed) <- colnames(datExpr_seed)
    names(MM_seed) <- colnames(datExpr_seed)
    
    wgcna_seed_table <- data.frame(
      Gene = targetGenes_seed,
      absGS = abs(GS_seed[targetGenes_seed]),
      absMM = abs(MM_seed[targetGenes_seed]),
      stringsAsFactors = FALSE
    )
    
    if (toupper(WGCNA_RULE) == "AND") {
      WGCNA_seed_genes <- wgcna_seed_table$Gene[
        wgcna_seed_table$absMM >= WGCNA_MM_MIN &
          wgcna_seed_table$absGS >= WGCNA_GS_MIN
      ]
    } else {
      WGCNA_seed_genes <- wgcna_seed_table$Gene[
        wgcna_seed_table$absMM >= WGCNA_MM_MIN |
          wgcna_seed_table$absGS >= WGCNA_GS_MIN
      ]
    }
    
    WGCNA_seed_genes <- unique(WGCNA_seed_genes)
    
    # -----------------------------------------------------------------------
    # E. TRIPLE OVERLAP
    # -----------------------------------------------------------------------
    
    triple_seed <- sort(
      unique(
        Reduce(
          intersect,
          list(
            Bulk_UP_seed,
            WGCNA_seed_genes,
            SC_MARKERS_STABILITY
          )
        )
      )
    )
    
    # MYBL2 direct AUC is calculated even if a particular seed does not carry
    # it through every screening step, which makes stability assessment clear.
    train_group_seed <- make_group_factor(train_seed_df$group)
    valid_group_seed <- make_group_factor(valid_seed_df$group)
    
    mybl2_train_auc_seed <- if (STABILITY_TARGET_GENE %in% colnames(train_seed_df)) {
      stability_gene_auc(
        train_seed_df[[STABILITY_TARGET_GENE]],
        train_group_seed
      )
    } else {
      NA_real_
    }
    
    mybl2_valid_auc_seed <- if (STABILITY_TARGET_GENE %in% colnames(valid_seed_df)) {
      stability_gene_auc(
        valid_seed_df[[STABILITY_TARGET_GENE]],
        valid_group_seed
      )
    } else {
      NA_real_
    }
    
    # -----------------------------------------------------------------------
    # F. LASSO + INTERNAL ROC
    # -----------------------------------------------------------------------
    
    common_model_seed <- Reduce(
      intersect,
      list(
        triple_seed,
        setdiff(colnames(train_seed_df), c("Sample", "group", "batch")),
        setdiff(colnames(valid_seed_df), c("Sample", "group", "batch")),
        setdiff(colnames(ext_seed_df), c("Sample", "group", "batch"))
      )
    )
    
    common_model_seed <- sort(unique(common_model_seed))
    
    lasso_genes_seed <- character(0)
    final_biomarkers_seed <- character(0)
    lambda_min_seed <- NA_real_
    
    if (length(common_model_seed) >= 2) {
      
      x_seed <- as.matrix(
        to_numeric_df(train_seed_df[, common_model_seed, drop = FALSE])
      )
      y_seed <- ifelse(train_group_seed == "Tumor", 1, 0)
      
      nfolds_seed <- min(LASSO_NFOLDS, length(y_seed))
      
      if (nfolds_seed >= 3) {
        
        # Fixed model-fitting seed: isolate the effect of SAMPLE-SELECTION seed.
        set.seed(STABILITY_MODEL_SEED)
        
        cv_seed <- glmnet::cv.glmnet(
          x = x_seed,
          y = y_seed,
          alpha = 1,
          family = "binomial",
          type.measure = "deviance",
          nfolds = nfolds_seed
        )
        
        lambda_min_seed <- cv_seed$lambda.min
        
        fit_lasso_seed <- glmnet::glmnet(
          x = x_seed,
          y = y_seed,
          alpha = 1,
          family = "binomial",
          lambda = lambda_min_seed
        )
        
        beta_seed <- as.vector(fit_lasso_seed$beta[, 1])
        names(beta_seed) <- rownames(fit_lasso_seed$beta)
        
        lasso_genes_seed <- names(beta_seed)[beta_seed != 0]
        lasso_genes_seed <- sort(unique(lasso_genes_seed))
        
        if (length(lasso_genes_seed) > 0) {
          
          seed_roc_table <- lapply(
            lasso_genes_seed,
            function(g) {
              data.frame(
                Gene = g,
                Train_AUC = stability_gene_auc(
                  train_seed_df[[g]],
                  train_group_seed
                ),
                Valid_AUC = stability_gene_auc(
                  valid_seed_df[[g]],
                  valid_group_seed
                ),
                stringsAsFactors = FALSE
              )
            }
          ) %>%
            dplyr::bind_rows() %>%
            dplyr::mutate(
              Pass_ROC_Cutoff = is.finite(Train_AUC) &
                is.finite(Valid_AUC) &
                Train_AUC >= ROC_AUC_MIN &
                Valid_AUC >= ROC_AUC_MIN
            )
          
          final_biomarkers_seed <- seed_roc_table$Gene[
            seed_roc_table$Pass_ROC_Cutoff
          ]
          
          final_biomarkers_seed <- sort(unique(final_biomarkers_seed))
        }
      }
    }
    
    # -----------------------------------------------------------------------
    # G. EXTERNAL RF/SVM USING THE SEED-SPECIFIC FINAL BIOMARKER SET
    # -----------------------------------------------------------------------
    
    external_ml_seed <- stability_external_ml(
      train_df_seed = train_seed_df,
      external_df_seed = ext_seed_df,
      model_genes_seed = final_biomarkers_seed
    )
    
    # -----------------------------------------------------------------------
    # H. SAVE COMPACT PER-SEED GENE RECORD
    # -----------------------------------------------------------------------
    
    gene_record_seed <- dplyr::bind_rows(
      stability_make_gene_record(seed_value, "Bulk_UP_DEG", Bulk_UP_seed),
      stability_make_gene_record(seed_value, "WGCNA_feature", WGCNA_seed_genes),
      stability_make_gene_record(seed_value, "Triple_overlap", triple_seed),
      stability_make_gene_record(seed_value, "LASSO", lasso_genes_seed),
      stability_make_gene_record(seed_value, "Final_biomarker", final_biomarkers_seed)
    )
    
    # -----------------------------------------------------------------------
    # I. ONE-ROW SUMMARY
    # -----------------------------------------------------------------------
    
    summary_seed <- data.frame(
      Seed = seed_value,
      Status = "SUCCESS",
      Prebalance_Tumor = length(all_tumor_ids),
      Prebalance_Normal = length(all_normal_ids),
      Balanced_per_class = n_balance_each,
      Train_Tumor = sum(train_seed_df$group == 2),
      Train_Normal = sum(train_seed_df$group == 1),
      Valid_Tumor = sum(valid_seed_df$group == 2),
      Valid_Normal = sum(valid_seed_df$group == 1),
      Balance_Excluded_Sample = ifelse(
        length(excluded_balance_ids) > 0,
        paste(excluded_balance_ids, collapse = ";"),
        "None"
      ),
      DEG_Up_Count = length(Bulk_UP_seed),
      WGCNA_Feature_Count = length(WGCNA_seed_genes),
      Triple_Overlap_Count = length(triple_seed),
      LASSO_Count = length(lasso_genes_seed),
      Final_Biomarker_Count = length(final_biomarkers_seed),
      MYBL2_in_Bulk_UP = STABILITY_TARGET_GENE %in% Bulk_UP_seed,
      MYBL2_in_WGCNA = STABILITY_TARGET_GENE %in% WGCNA_seed_genes,
      MYBL2_in_Triple = STABILITY_TARGET_GENE %in% triple_seed,
      MYBL2_in_LASSO = STABILITY_TARGET_GENE %in% lasso_genes_seed,
      MYBL2_Final = STABILITY_TARGET_GENE %in% final_biomarkers_seed,
      MYBL2_Train_AUC = mybl2_train_auc_seed,
      MYBL2_Valid_AUC = mybl2_valid_auc_seed,
      External_RF_AUROC = unname(external_ml_seed["RF_AUROC"]),
      External_SVM_AUROC = unname(external_ml_seed["SVM_AUROC"]),
      Lambda_min = lambda_min_seed,
      Final_Biomarkers = ifelse(
        length(final_biomarkers_seed) > 0,
        paste(final_biomarkers_seed, collapse = ";"),
        "None"
      ),
      Error_Message = "",
      stringsAsFactors = FALSE
    )
    
    cat(
      "Seed", seed_value,
      "| MYBL2 train AUC =", sprintf("%.4f", mybl2_train_auc_seed),
      "| validation AUC =", sprintf("%.4f", mybl2_valid_auc_seed),
      "| MYBL2 final =", STABILITY_TARGET_GENE %in% final_biomarkers_seed,
      "\n"
    )
    
    rm(
      datExpr0_seed,
      datExpr_seed,
      net_seed,
      MEs_seed,
      E_seed,
      fit_seed
    )
    gc(verbose = FALSE)
    
    list(
      summary = summary_seed,
      genes = gene_record_seed
    )
    
  }, error = function(e) {
    
    cat("Seed", seed_value, "FAILED:", conditionMessage(e), "\n")
    
    failed_summary <- data.frame(
      Seed = seed_value,
      Status = "FAILED",
      Prebalance_Tumor = length(all_tumor_ids),
      Prebalance_Normal = length(all_normal_ids),
      Balanced_per_class = n_balance_each,
      Train_Tumor = NA_integer_,
      Train_Normal = NA_integer_,
      Valid_Tumor = NA_integer_,
      Valid_Normal = NA_integer_,
      Balance_Excluded_Sample = NA_character_,
      DEG_Up_Count = NA_integer_,
      WGCNA_Feature_Count = NA_integer_,
      Triple_Overlap_Count = NA_integer_,
      LASSO_Count = NA_integer_,
      Final_Biomarker_Count = NA_integer_,
      MYBL2_in_Bulk_UP = NA,
      MYBL2_in_WGCNA = NA,
      MYBL2_in_Triple = NA,
      MYBL2_in_LASSO = NA,
      MYBL2_Final = NA,
      MYBL2_Train_AUC = NA_real_,
      MYBL2_Valid_AUC = NA_real_,
      External_RF_AUROC = NA_real_,
      External_SVM_AUROC = NA_real_,
      Lambda_min = NA_real_,
      Final_Biomarkers = NA_character_,
      Error_Message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
    
    list(
      summary = failed_summary,
      genes = data.frame(
        Seed = integer(0),
        Stage = character(0),
        Gene = character(0),
        stringsAsFactors = FALSE
      )
    )
  })
}

# ---------------------------------------------------------------------------
# 12. RUN ALL SAMPLING SEEDS
# ---------------------------------------------------------------------------

stability_runs <- lapply(STABILITY_SEEDS, run_stability_seed)

seed_stability_summary <- dplyr::bind_rows(
  lapply(stability_runs, `[[`, "summary")
)

seed_stability_genes <- dplyr::bind_rows(
  lapply(stability_runs, `[[`, "genes")
)

write.csv(
  seed_stability_summary,
  file.path(stability_dir, "Seed_stability_summary.csv"),
  row.names = FALSE
)

write.csv(
  seed_stability_genes,
  file.path(stability_dir, "Seed_stability_selected_genes_long.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# 13. SUMMARIZE AUC VALUES
# ---------------------------------------------------------------------------

successful_seed_summary <- seed_stability_summary %>%
  dplyr::filter(Status == "SUCCESS")

if (nrow(successful_seed_summary) < 2) {
  stop(
    "Fewer than two sampling seeds completed successfully; cannot create a stability figure."
  )
}

seed_auc_long <- successful_seed_summary %>%
  dplyr::select(
    Seed,
    MYBL2_Train_AUC,
    MYBL2_Valid_AUC,
    External_RF_AUROC,
    External_SVM_AUROC
  ) %>%
  tidyr::pivot_longer(
    cols = -Seed,
    names_to = "Metric",
    values_to = "AUROC"
  ) %>%
  dplyr::mutate(
    Metric = dplyr::recode(
      Metric,
      MYBL2_Train_AUC = "MYBL2 training AUC",
      MYBL2_Valid_AUC = "MYBL2 validation AUC",
      External_RF_AUROC = "External RF AUROC",
      External_SVM_AUROC = "External SVM AUROC"
    ),
    Seed = factor(Seed, levels = STABILITY_SEEDS)
  )

write.csv(
  seed_auc_long,
  file.path(stability_dir, "Seed_stability_AUC_long.csv"),
  row.names = FALSE
)

# Metric summary: mean, SD, min, max across successful sampling seeds.
seed_metric_summary <- seed_auc_long %>%
  dplyr::group_by(Metric) %>%
  dplyr::summarise(
    N_Seeds = sum(is.finite(AUROC)),
    Mean = mean(AUROC, na.rm = TRUE),
    SD = sd(AUROC, na.rm = TRUE),
    Min = min(AUROC, na.rm = TRUE),
    Max = max(AUROC, na.rm = TRUE),
    Range = Max - Min,
    .groups = "drop"
  )

write.csv(
  seed_metric_summary,
  file.path(stability_dir, "Seed_stability_metric_summary.csv"),
  row.names = FALSE
)

# MYBL2 selection frequency at key screening stages.
mybl2_final_n <- sum(
  successful_seed_summary$MYBL2_Final %in% TRUE,
  na.rm = TRUE
)

mybl2_lasso_n <- sum(
  successful_seed_summary$MYBL2_in_LASSO %in% TRUE,
  na.rm = TRUE
)

mybl2_triple_n <- sum(
  successful_seed_summary$MYBL2_in_Triple %in% TRUE,
  na.rm = TRUE
)

n_success_seeds <- nrow(successful_seed_summary)

# ---------------------------------------------------------------------------
# 14. STABILITY FIGURE
# ---------------------------------------------------------------------------

metric_order <- c(
  "MYBL2 training AUC",
  "MYBL2 validation AUC",
  "External RF AUROC",
  "External SVM AUROC"
)

seed_auc_long$Metric <- factor(seed_auc_long$Metric, levels = metric_order)

stability_colors <- c(
  "MYBL2 training AUC" = COL_NORMAL,
  "MYBL2 validation AUC" = COL_TUMOR,
  "External RF AUROC" = "#4C78A8",
  "External SVM AUROC" = "#2A9D8F"
)

finite_auc_values <- seed_auc_long$AUROC[is.finite(seed_auc_long$AUROC)]
plot_y_min <- if (length(finite_auc_values) > 0) {
  max(0.50, floor((min(finite_auc_values) - 0.03) * 20) / 20)
} else {
  0.50
}

p_seed_stability <- ggplot2::ggplot(
  seed_auc_long,
  ggplot2::aes(
    x = Seed,
    y = AUROC,
    group = Metric,
    color = Metric
  )
) +
  ggplot2::geom_hline(
    yintercept = ROC_AUC_MIN,
    linetype = "dashed",
    linewidth = 0.55,
    color = "grey45"
  ) +
  ggplot2::geom_line(
    linewidth = 0.9,
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    size = 2.6,
    na.rm = TRUE
  ) +
  ggplot2::scale_color_manual(
    values = stability_colors,
    drop = FALSE
  ) +
  ggplot2::coord_cartesian(
    ylim = c(plot_y_min, 1.00)
  ) +
  ggplot2::labs(
    title = "Stability across random sampling seeds",
    x = "Sampling seed",
    y = "AUC / AUROC",
    color = NULL
  ) +
  theme_manuscript(
    show_grid = FALSE,
    legend_position = "bottom"
  ) +
  ggplot2::theme(
    legend.direction = "horizontal",
    legend.box = "vertical"
  )

save_plot_all_formats(
  plot_obj = p_seed_stability,
  filename_stem = "Figure_Seed_Stability_Across_Sampling_Seeds",
  dir_path = stability_dir,
  width = 8.8,
  height = 5.8
)

# ---------------------------------------------------------------------------
# 15. TEXT SUMMARY
# ---------------------------------------------------------------------------

stability_summary_lines <- c(
  "MULTI-SEED STABILITY ANALYSIS",
  paste0("Seed mode: ", SEED_MODE),
  paste0("Sampling seeds tested: ", paste(STABILITY_SEEDS, collapse = ", ")),
  paste0("Successful seed runs: ", n_success_seeds, "/", length(STABILITY_SEEDS)),
  paste0(
    "MYBL2 present in triple overlap: ",
    mybl2_triple_n, "/", n_success_seeds
  ),
  paste0(
    "MYBL2 selected by LASSO: ",
    mybl2_lasso_n, "/", n_success_seeds
  ),
  paste0(
    "MYBL2 retained as final biomarker: ",
    mybl2_final_n, "/", n_success_seeds
  ),
  "",
  "AUC/AUROC summary across successful sampling seeds:",
  capture.output(print(seed_metric_summary, row.names = FALSE)),
  "",
  paste0(
    "Model-fitting randomness was held fixed at seed ",
    STABILITY_MODEL_SEED,
    "; only the sample-selection seed was varied."
  )
)

writeLines(
  stability_summary_lines,
  file.path(stability_dir, "Seed_stability_summary.txt")
)

# ---------------------------------------------------------------------------
# 16. CONSOLE SUMMARY
# ---------------------------------------------------------------------------

cat("\n============================================================\n")
cat("MULTI-SEED STABILITY ANALYSIS FINISHED\n")
cat("============================================================\n")
cat("Successful seeds:", n_success_seeds, "/", length(STABILITY_SEEDS), "\n")
cat("MYBL2 in triple overlap:", mybl2_triple_n, "/", n_success_seeds, "\n")
cat("MYBL2 selected by LASSO:", mybl2_lasso_n, "/", n_success_seeds, "\n")
cat("MYBL2 final biomarker:", mybl2_final_n, "/", n_success_seeds, "\n")
cat("\nSeed-level results:\n")
print(seed_stability_summary)
cat("\nMetric stability summary:\n")
print(seed_metric_summary)
cat("\nOutputs saved to:\n", stability_dir, "\n", sep = "")
cat("Stability figure:\n")
cat(
  file.path(
    stability_dir,
    "Figure_Seed_Stability_Across_Sampling_Seeds.png"
  ),
  "\n"
)

# End of script
