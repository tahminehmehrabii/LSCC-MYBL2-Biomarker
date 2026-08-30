# RANDOM FOREST VERSUS SVM COMPARISON
# LSCC MYBL2 reproducible analysis
#
# Purpose:
#   Statistically compare the external AUROCs of the Random Forest and
#   linear SVM models in GSE130605 using a paired DeLong test.
#
# Design:
#   - Same external samples are evaluated by both models -> paired comparison.
#   - Reuses rf_prob, svm_prob, and testY_ml if they already exist in memory.
#   - Otherwise reconstructs the RF/SVM analysis from saved pipeline files
#     using the primary model settings:
#       RF: 500 trees
#       SVM: linear kernel
#       model-fitting seed: 123
#
# Main outputs:
#   1) RF_vs_SVM_DeLong_test.csv
#   2) RF_vs_SVM_external_predictions.csv
#   3) RF_vs_SVM_DeLong_summary.txt
###############################################################################

cat("\n============================================================\n")
cat("PAIRED DELONG TEST: EXTERNAL RF vs SVM AUROC\n")
cat("============================================================\n")

# -----------------------------------------------------------------------------
# 1. SETTINGS
# -----------------------------------------------------------------------------

results_path <- "E:/LSCC/Results_LSCC/ML"
DELONG_OUTPUT_DIR <- file.path(results_path, "RF_vs_SVM_DeLong_Test")
dir.create(DELONG_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

RF_NTREE_DELONG <- 500
SVM_KERNEL_DELONG <- "linear"
MODEL_SEED_DELONG <- 123

# -----------------------------------------------------------------------------
# 2. REQUIRED PACKAGES
# -----------------------------------------------------------------------------

required_packages_delong <- c(
  "data.table",
  "pROC",
  "randomForest",
  "e1071"
)

missing_packages_delong <- required_packages_delong[
  !vapply(
    required_packages_delong,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages_delong) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages_delong, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# 3. LOCAL HELPERS
# -----------------------------------------------------------------------------

clean_gene_vector_delong <- function(x) {
  x <- unique(trimws(as.character(x)))
  x[!is.na(x) & x != ""]
}

to_numeric_df_delong <- function(df) {
  as.data.frame(
    lapply(df, function(x) as.numeric(as.character(x))),
    check.names = FALSE
  )
}

coerce_group_12_delong <- function(g) {
  g0 <- as.character(g)
  
  if (all(g0 %in% c("1", "2"))) {
    return(as.integer(g0))
  }
  
  g_low <- tolower(g0)
  
  if (all(g_low %in% c(
    "normal", "tumor", "non", "lscc", "cancer", "margin"
  ))) {
    return(
      ifelse(
        g_low %in% c("tumor", "lscc", "cancer"),
        2L,
        1L
      )
    )
  }
  
  u <- sort(unique(g0))
  
  if (length(u) != 2) {
    stop(
      "group must contain exactly two classes. Found: ",
      paste(u, collapse = ", ")
    )
  }
  
  map <- setNames(c(1L, 2L), u)
  as.integer(map[g0])
}

make_group_factor_delong <- function(g) {
  g12 <- coerce_group_12_delong(g)
  factor(
    g12,
    levels = c(1, 2),
    labels = c("Normal", "Tumor")
  )
}

scale_by_train_delong <- function(train_mat, new_mat) {
  train_mat <- as.matrix(train_mat)
  new_mat <- as.matrix(new_mat)
  
  mu <- colMeans(train_mat, na.rm = TRUE)
  sdv <- apply(train_mat, 2, stats::sd, na.rm = TRUE)
  
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  
  sweep(
    sweep(new_mat, 2, mu, "-"),
    2,
    sdv,
    "/"
  )
}

# -----------------------------------------------------------------------------
# 4. REUSE EXISTING PREDICTIONS OR RECONSTRUCT MODELS
# -----------------------------------------------------------------------------

can_reuse_predictions <- all(
  c("testY_ml", "rf_prob", "svm_prob") %in% ls(envir = .GlobalEnv)
)

if (can_reuse_predictions) {
  
  cat("\nUsing RF/SVM external probabilities already present in memory.\n")
  
  y_external <- factor(
    testY_ml,
    levels = c("Normal", "Tumor")
  )
  
  rf_probability <- as.numeric(rf_prob)
  svm_probability <- as.numeric(svm_prob)
  
  if (exists("test_df", inherits = FALSE) &&
      "Sample" %in% colnames(test_df) &&
      nrow(test_df) == length(y_external)) {
    external_sample_ids <- as.character(test_df$Sample)
  } else {
    external_sample_ids <- paste0("Sample_", seq_along(y_external))
  }
  
} else {
  
  cat("\nExisting predictions not found; reconstructing RF and SVM models...\n")
  
  train_log2_file_delong <- file.path(
    results_path,
    "train_log2CPM.csv"
  )
  
  external_log2_file_delong <- file.path(
    results_path,
    "external_GSE130605_log2CPM.csv"
  )
  
  final_biomarkers_file_delong <- file.path(
    results_path,
    "FINAL_BIOMARKERS.txt"
  )
  
  required_input_files_delong <- c(
    train_log2_file_delong,
    external_log2_file_delong,
    final_biomarkers_file_delong
  )
  
  missing_input_files_delong <- required_input_files_delong[
    !file.exists(required_input_files_delong)
  ]
  
  if (length(missing_input_files_delong) > 0) {
    stop(
      "Required pipeline file(s) not found:\n",
      paste(missing_input_files_delong, collapse = "\n"),
      "\n\nRun the main pipeline first, or keep ",
      "testY_ml/rf_prob/svm_prob in the current R session."
    )
  }
  
  train_df_delong <- as.data.frame(
    data.table::fread(train_log2_file_delong),
    check.names = FALSE
  )
  
  test_df_delong <- as.data.frame(
    data.table::fread(external_log2_file_delong),
    check.names = FALSE
  )
  
  final_biomarkers_delong <- clean_gene_vector_delong(
    readLines(
      final_biomarkers_file_delong,
      warn = FALSE
    )
  )
  
  model_genes_delong <- intersect(
    final_biomarkers_delong,
    intersect(
      setdiff(
        colnames(train_df_delong),
        c("Sample", "group", "batch")
      ),
      setdiff(
        colnames(test_df_delong),
        c("Sample", "group", "batch")
      )
    )
  )
  
  model_genes_delong <- sort(unique(model_genes_delong))
  
  if (length(model_genes_delong) < 1) {
    stop(
      "No final biomarker is available in both training and external matrices."
    )
  }
  
  cat(
    "Model gene(s):",
    paste(model_genes_delong, collapse = ", "),
    "\n"
  )
  
  trainX_delong <- as.matrix(
    to_numeric_df_delong(
      train_df_delong[, model_genes_delong, drop = FALSE]
    )
  )
  
  testX_delong <- as.matrix(
    to_numeric_df_delong(
      test_df_delong[, model_genes_delong, drop = FALSE]
    )
  )
  
  trainY_delong <- make_group_factor_delong(train_df_delong$group)
  y_external <- make_group_factor_delong(test_df_delong$group)
  
  train_scaled_delong <- scale_by_train_delong(
    trainX_delong,
    trainX_delong
  )
  
  test_scaled_delong <- scale_by_train_delong(
    trainX_delong,
    testX_delong
  )
  
  train_ml_delong <- data.frame(
    group = trainY_delong,
    train_scaled_delong,
    check.names = FALSE
  )
  
  test_ml_delong <- data.frame(
    group = y_external,
    test_scaled_delong,
    check.names = FALSE
  )
  
  colnames(train_ml_delong) <- make.names(
    colnames(train_ml_delong),
    unique = TRUE
  )
  
  colnames(test_ml_delong) <- make.names(
    colnames(test_ml_delong),
    unique = TRUE
  )
  
  # Random forest
  set.seed(MODEL_SEED_DELONG)
  
  rf_model_delong <- randomForest::randomForest(
    group ~ .,
    data = train_ml_delong,
    ntree = RF_NTREE_DELONG,
    importance = TRUE
  )
  
  rf_probability <- predict(
    rf_model_delong,
    newdata = test_ml_delong,
    type = "prob"
  )[, "Tumor"]
  
  # Linear SVM
  set.seed(MODEL_SEED_DELONG)
  
  svm_model_delong <- e1071::svm(
    group ~ .,
    data = train_ml_delong,
    kernel = SVM_KERNEL_DELONG,
    probability = TRUE
  )
  
  svm_pred_delong <- predict(
    svm_model_delong,
    newdata = test_ml_delong,
    probability = TRUE
  )
  
  svm_prob_mat_delong <- attr(
    svm_pred_delong,
    "probabilities"
  )
  
  if (is.null(svm_prob_mat_delong)) {
    stop(
      "The SVM model did not return probability estimates."
    )
  }
  
  svm_probability <- if (
    "Tumor" %in% colnames(svm_prob_mat_delong)
  ) {
    svm_prob_mat_delong[, "Tumor"]
  } else if (
    "Normal" %in% colnames(svm_prob_mat_delong)
  ) {
    1 - svm_prob_mat_delong[, "Normal"]
  } else {
    stop(
      "Could not identify Tumor/Normal probability columns in SVM output."
    )
  }
  
  external_sample_ids <- as.character(test_df_delong$Sample)
}

# -----------------------------------------------------------------------------
# 5. INPUT VALIDATION
# -----------------------------------------------------------------------------

if (
  length(y_external) != length(rf_probability) ||
  length(y_external) != length(svm_probability)
) {
  stop(
    "Outcome, RF probability, and SVM probability vectors have different lengths."
  )
}

complete_cases_delong <- (
  !is.na(y_external) &
    is.finite(rf_probability) &
    is.finite(svm_probability)
)

if (!all(complete_cases_delong)) {
  warning(
    sum(!complete_cases_delong),
    " external sample(s) had missing/non-finite values and were removed."
  )
  
  y_external <- y_external[complete_cases_delong]
  rf_probability <- rf_probability[complete_cases_delong]
  svm_probability <- svm_probability[complete_cases_delong]
  external_sample_ids <- external_sample_ids[complete_cases_delong]
}

if (length(unique(y_external)) != 2) {
  stop("Both Normal and Tumor classes are required for the DeLong test.")
}

# -----------------------------------------------------------------------------
# 6. EXTERNAL ROC CURVES
# -----------------------------------------------------------------------------

roc_rf_delong <- pROC::roc(
  response = y_external,
  predictor = rf_probability,
  levels = c("Normal", "Tumor"),
  direction = "<",
  quiet = TRUE
)

roc_svm_delong <- pROC::roc(
  response = y_external,
  predictor = svm_probability,
  levels = c("Normal", "Tumor"),
  direction = "<",
  quiet = TRUE
)

rf_auc_delong <- as.numeric(
  pROC::auc(roc_rf_delong)
)

svm_auc_delong <- as.numeric(
  pROC::auc(roc_svm_delong)
)

rf_ci_delong <- as.numeric(
  pROC::ci.auc(
    roc_rf_delong,
    method = "delong",
    conf.level = 0.95
  )
)

svm_ci_delong <- as.numeric(
  pROC::ci.auc(
    roc_svm_delong,
    method = "delong",
    conf.level = 0.95
  )
)

# -----------------------------------------------------------------------------
# 7. PAIRED DELONG TEST
# -----------------------------------------------------------------------------
# Paired is appropriate because RF and SVM were evaluated on the SAME
# GSE130605 samples.

delong_test <- pROC::roc.test(
  roc_rf_delong,
  roc_svm_delong,
  method = "delong",
  paired = TRUE,
  alternative = "two.sided",
  conf.level = 0.95
)

delong_statistic <- if (
  !is.null(delong_test$statistic) &&
  length(delong_test$statistic) > 0
) {
  as.numeric(delong_test$statistic[1])
} else {
  NA_real_
}

delong_statistic_name <- if (
  !is.null(names(delong_test$statistic)) &&
  length(names(delong_test$statistic)) > 0
) {
  names(delong_test$statistic)[1]
} else {
  "Statistic"
}

delong_p_value <- as.numeric(delong_test$p.value)

auc_difference <- svm_auc_delong - rf_auc_delong

# -----------------------------------------------------------------------------
# 8. SAVE EXTERNAL PREDICTIONS
# -----------------------------------------------------------------------------

prediction_table <- data.frame(
  Sample = external_sample_ids,
  True_Class = as.character(y_external),
  RF_Tumor_Probability = as.numeric(rf_probability),
  SVM_Tumor_Probability = as.numeric(svm_probability),
  stringsAsFactors = FALSE
)

write.csv(
  prediction_table,
  file.path(
    DELONG_OUTPUT_DIR,
    "RF_vs_SVM_external_predictions.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 9. SAVE DELONG RESULT TABLE
# -----------------------------------------------------------------------------

delong_result_table <- data.frame(
  External_Cohort = "GSE130605",
  N = length(y_external),
  Normal_N = sum(y_external == "Normal"),
  Tumor_N = sum(y_external == "Tumor"),
  
  RF_AUROC = rf_auc_delong,
  RF_AUROC_CI95_Lower = rf_ci_delong[1],
  RF_AUROC_CI95_Upper = rf_ci_delong[3],
  
  SVM_AUROC = svm_auc_delong,
  SVM_AUROC_CI95_Lower = svm_ci_delong[1],
  SVM_AUROC_CI95_Upper = svm_ci_delong[3],
  
  AUROC_Difference_SVM_minus_RF = auc_difference,
  
  DeLong_Test = "Paired, two-sided",
  DeLong_Statistic_Name = delong_statistic_name,
  DeLong_Statistic = delong_statistic,
  DeLong_P_Value = delong_p_value,
  Statistically_Significant_P_lt_0.05 = delong_p_value < 0.05,
  
  stringsAsFactors = FALSE
)

write.csv(
  delong_result_table,
  file.path(
    DELONG_OUTPUT_DIR,
    "RF_vs_SVM_DeLong_test.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 10. TEXT SUMMARY
# -----------------------------------------------------------------------------

interpretation_line <- if (delong_p_value < 0.05) {
  paste0(
    "The difference between the external RF and SVM AUROCs was statistically ",
    "significant by the paired DeLong test (P = ",
    format.pval(delong_p_value, digits = 4),
    ")."
  )
} else {
  paste0(
    "The difference between the external RF and SVM AUROCs was not statistically ",
    "significant by the paired DeLong test (P = ",
    format.pval(delong_p_value, digits = 4),
    ")."
  )
}

delong_summary <- c(
  "PAIRED DELONG TEST: RF vs SVM EXTERNAL AUROC",
  paste0("External cohort: GSE130605 (n = ", length(y_external), ")"),
  paste0(
    "RF AUROC = ",
    sprintf("%.4f", rf_auc_delong),
    " (95% CI ",
    sprintf("%.4f", rf_ci_delong[1]),
    "-",
    sprintf("%.4f", rf_ci_delong[3]),
    ")"
  ),
  paste0(
    "SVM AUROC = ",
    sprintf("%.4f", svm_auc_delong),
    " (95% CI ",
    sprintf("%.4f", svm_ci_delong[1]),
    "-",
    sprintf("%.4f", svm_ci_delong[3]),
    ")"
  ),
  paste0(
    "AUROC difference (SVM - RF) = ",
    sprintf("%.4f", auc_difference)
  ),
  paste0(
    "Paired DeLong ",
    delong_statistic_name,
    " = ",
    sprintf("%.4f", delong_statistic),
    ", P = ",
    format.pval(delong_p_value, digits = 4)
  ),
  interpretation_line
)

writeLines(
  delong_summary,
  file.path(
    DELONG_OUTPUT_DIR,
    "RF_vs_SVM_DeLong_summary.txt"
  )
)

# -----------------------------------------------------------------------------
# 11. CONSOLE OUTPUT
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("DELONG TEST FINISHED\n")
cat("============================================================\n")
cat("External cohort: GSE130605\n")
cat("N =", length(y_external), "\n")
cat("Normal =", sum(y_external == "Normal"), "\n")
cat("Tumor =", sum(y_external == "Tumor"), "\n\n")

cat(
  "RF AUROC =",
  sprintf("%.4f", rf_auc_delong),
  " | 95% CI:",
  sprintf("%.4f", rf_ci_delong[1]),
  "-",
  sprintf("%.4f", rf_ci_delong[3]),
  "\n"
)

cat(
  "SVM AUROC =",
  sprintf("%.4f", svm_auc_delong),
  " | 95% CI:",
  sprintf("%.4f", svm_ci_delong[1]),
  "-",
  sprintf("%.4f", svm_ci_delong[3]),
  "\n"
)

cat(
  "Difference (SVM - RF) =",
  sprintf("%.4f", auc_difference),
  "\n"
)

cat(
  "Paired DeLong",
  delong_statistic_name,
  "=",
  sprintf("%.4f", delong_statistic),
  "\n"
)

cat(
  "P-value =",
  format.pval(delong_p_value, digits = 6),
  "\n"
)

cat(
  "Significant at P < 0.05:",
  ifelse(delong_p_value < 0.05, "YES", "NO"),
  "\n"
)

cat("\n", interpretation_line, "\n", sep = "")

cat(
  "\nOutputs saved to:\n",
  DELONG_OUTPUT_DIR,
  "\n",
  sep = ""
)

cat("\nKey files:\n")
cat(" - RF_vs_SVM_DeLong_test.csv\n")
cat(" - RF_vs_SVM_external_predictions.csv\n")
cat(" - RF_vs_SVM_DeLong_summary.txt\n")

# End of script
