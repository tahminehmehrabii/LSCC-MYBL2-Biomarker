###############################################################################
# MULTIGENE-PANEL BENCHMARK
# LSCC MYBL2 reproducible analysis
#
# Compares MYBL2 and the other LASSO-selected genes with a five-gene
# logistic-regression panel in the training, validation, and external cohorts.
###############################################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(scipen = 999)
options(width = 200)

set.seed(123)

# =============================================================================
# 1. PATHS
# =============================================================================

RESULTS_PATH <- "E:/LSCC/Results_LSCC/ML"

OUT_DIR <- file.path(
  RESULTS_PATH,
  "Multigene_Panel_Benchmark"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

TRAIN_FILE <- file.path(
  RESULTS_PATH,
  "train_log2CPM.csv"
)

VALID_FILE <- file.path(
  RESULTS_PATH,
  "valid_log2CPM.csv"
)

EXTERNAL_FILE <- file.path(
  RESULTS_PATH,
  "external_GSE130605_log2CPM.csv"
)

# =============================================================================
# 2. GENES
# =============================================================================

GENES <- c(
  "MYBL2",
  "KRT17",
  "UCHL1",
  "PTHLH",
  "ODC1"
)

TARGET_GENE <- "MYBL2"

PANEL_KEY <- "FiveGenePanel"

PANEL_LABEL <- "Five-gene panel"

# =============================================================================
# 3. REQUIRED PACKAGES
# =============================================================================

required_packages <- c(
  "pROC",
  "ggplot2",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  
  stop(
    "Missing required package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\nPlease install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
  library(patchwork)
})

# =============================================================================
# 4. HELPER FUNCTIONS
# =============================================================================

safe_numeric <- function(x) {
  
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}


make_group_factor <- function(g) {
  
  g0 <- trimws(
    as.character(g)
  )
  
  # Original coding:
  # 1 = Normal
  # 2 = Tumor
  
  if (all(g0 %in% c("1", "2"))) {
    
    return(
      factor(
        g0,
        levels = c("1", "2"),
        labels = c("Normal", "Tumor")
      )
    )
  }
  
  gl <- tolower(g0)
  
  allowed <- c(
    "normal",
    "tumor",
    "non",
    "lscc",
    "cancer",
    "margin"
  )
  
  if (all(gl %in% allowed)) {
    
    label <- ifelse(
      gl %in% c(
        "tumor",
        "lscc",
        "cancer"
      ),
      "Tumor",
      "Normal"
    )
    
    return(
      factor(
        label,
        levels = c(
          "Normal",
          "Tumor"
        )
      )
    )
  }
  
  stop(
    "Could not recognize group labels: ",
    paste(
      unique(g0),
      collapse = ", "
    )
  )
}


make_roc <- function(
    response,
    predictor,
    direction) {
  
  ok <- (
    !is.na(response) &
      is.finite(predictor)
  )
  
  y <- response[ok]
  x <- predictor[ok]
  
  if (
    length(y) < 4 ||
    length(unique(y)) < 2
  ) {
    
    return(NULL)
  }
  
  pROC::roc(
    response = y,
    predictor = x,
    levels = c(
      "Normal",
      "Tumor"
    ),
    direction = direction,
    quiet = TRUE
  )
}


get_auc_ci <- function(
    cohort_name,
    model_key,
    model_label,
    model_type,
    response,
    predictor,
    direction) {
  
  ok <- (
    !is.na(response) &
      is.finite(predictor)
  )
  
  y <- response[ok]
  x <- predictor[ok]
  
  rr <- make_roc(
    response = y,
    predictor = x,
    direction = direction
  )
  
  if (is.null(rr)) {
    
    return(
      data.frame(
        Cohort = cohort_name,
        Model_Key = model_key,
        Model = model_label,
        Model_Type = model_type,
        N = length(y),
        Normal_N = sum(y == "Normal"),
        Tumor_N = sum(y == "Tumor"),
        ROC_Direction = direction,
        AUROC = NA_real_,
        AUROC_CI95_Lower = NA_real_,
        AUROC_CI95_Upper = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }
  
  auc_value <- as.numeric(
    pROC::auc(rr)
  )
  
  ci_value <- tryCatch(
    as.numeric(
      pROC::ci.auc(
        rr,
        method = "delong",
        conf.level = 0.95
      )
    ),
    error = function(e) {
      c(
        NA_real_,
        NA_real_,
        NA_real_
      )
    }
  )
  
  data.frame(
    Cohort = cohort_name,
    Model_Key = model_key,
    Model = model_label,
    Model_Type = model_type,
    N = length(y),
    Normal_N = sum(y == "Normal"),
    Tumor_N = sum(y == "Tumor"),
    ROC_Direction = direction,
    AUROC = auc_value,
    AUROC_CI95_Lower = ci_value[1],
    AUROC_CI95_Upper = ci_value[3],
    stringsAsFactors = FALSE
  )
}


apply_training_scaling <- function(
    df,
    means,
    sds) {
  
  out <- df
  
  for (gene in names(means)) {
    
    out[[gene]] <- (
      safe_numeric(
        out[[gene]]
      ) -
        means[[gene]]
    ) /
      sds[[gene]]
  }
  
  out
}


save_plot_all_formats <- function(
    plot_object,
    stem,
    width = 16.5,
    height = 6.3) {
  
  ggplot2::ggsave(
    filename = file.path(
      OUT_DIR,
      paste0(stem, ".png")
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
  
  ggplot2::ggsave(
    filename = file.path(
      OUT_DIR,
      paste0(stem, ".tiff")
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white",
    limitsize = FALSE
  )
  
  tryCatch(
    {
      
      if (capabilities("cairo")) {
        
        ggplot2::ggsave(
          filename = file.path(
            OUT_DIR,
            paste0(stem, ".pdf")
          ),
          plot = plot_object,
          width = width,
          height = height,
          units = "in",
          device = grDevices::cairo_pdf,
          bg = "white",
          limitsize = FALSE
        )
        
      } else {
        
        ggplot2::ggsave(
          filename = file.path(
            OUT_DIR,
            paste0(stem, ".pdf")
          ),
          plot = plot_object,
          width = width,
          height = height,
          units = "in",
          device = "pdf",
          bg = "white",
          limitsize = FALSE
        )
      }
      
    },
    error = function(e) {
      
      warning(
        "PDF export failed: ",
        conditionMessage(e)
      )
    }
  )
}

# =============================================================================
# 5. CHECK INPUT FILES
# =============================================================================

required_files <- c(
  TRAIN_FILE,
  VALID_FILE,
  EXTERNAL_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    "Required file(s) missing:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

# =============================================================================
# 6. LOAD COHORTS
# =============================================================================

train_df <- read.csv(
  TRAIN_FILE,
  check.names = FALSE
)

valid_df <- read.csv(
  VALID_FILE,
  check.names = FALSE
)

external_df <- read.csv(
  EXTERNAL_FILE,
  check.names = FALSE
)

cohorts <- list(
  Training = train_df,
  Validation = valid_df,
  External_GSE130605 = external_df
)

# =============================================================================
# 7. CHECK GROUPS AND REQUIRED GENES
# =============================================================================

for (cohort_name in names(cohorts)) {
  
  dat <- cohorts[[cohort_name]]
  
  if (!("group" %in% names(dat))) {
    
    stop(
      "'group' column is missing from ",
      cohort_name
    )
  }
  
  dat$Class <- make_group_factor(
    dat$group
  )
  
  missing_genes <- setdiff(
    GENES,
    names(dat)
  )
  
  if (length(missing_genes) > 0) {
    
    stop(
      "Missing gene(s) in ",
      cohort_name,
      ": ",
      paste(
        missing_genes,
        collapse = ", "
      )
    )
  }
  
  # Convert gene-expression columns to numeric
  
  for (gene in GENES) {
    
    dat[[gene]] <- safe_numeric(
      dat[[gene]]
    )
    
    if (any(!is.finite(dat[[gene]]))) {
      
      stop(
        "Non-finite expression values detected for ",
        gene,
        " in ",
        cohort_name
      )
    }
  }
  
  cohorts[[cohort_name]] <- dat
}

train_df <- cohorts$Training
valid_df <- cohorts$Validation
external_df <- cohorts$External_GSE130605

# =============================================================================
# 8. COHORT SAMPLE SUMMARY
# =============================================================================

cat("\n============================================================\n")
cat("COHORT SAMPLE SIZES\n")
cat("============================================================\n")

for (cohort_name in names(cohorts)) {
  
  dat <- cohorts[[cohort_name]]
  
  cat(
    cohort_name,
    ": N = ",
    nrow(dat),
    " | Normal = ",
    sum(dat$Class == "Normal"),
    " | Tumor = ",
    sum(dat$Class == "Tumor"),
    "\n",
    sep = ""
  )
}

# =============================================================================
# 9. LEARN SINGLE-GENE ROC DIRECTIONS FROM TRAINING ONLY
# =============================================================================

gene_directions <- setNames(
  rep(
    NA_character_,
    length(GENES)
  ),
  GENES
)

for (gene in GENES) {
  
  x <- train_df[[gene]]
  y <- train_df$Class
  
  normal_median <- median(
    x[
      y == "Normal"
    ],
    na.rm = TRUE
  )
  
  tumor_median <- median(
    x[
      y == "Tumor"
    ],
    na.rm = TRUE
  )
  
  gene_directions[[gene]] <- ifelse(
    tumor_median >= normal_median,
    "<",
    ">"
  )
}

cat("\n============================================================\n")
cat("ROC DIRECTIONS — LEARNED FROM TRAINING ONLY\n")
cat("============================================================\n")

print(
  data.frame(
    Gene = GENES,
    ROC_Direction = unname(
      gene_directions[GENES]
    )
  ),
  row.names = FALSE
)

# =============================================================================
# 10. TRAINING-ONLY STANDARDIZATION FOR FIVE-GENE PANEL
# =============================================================================

train_gene_matrix <- train_df[
  ,
  GENES,
  drop = FALSE
]

valid_gene_matrix <- valid_df[
  ,
  GENES,
  drop = FALSE
]

external_gene_matrix <- external_df[
  ,
  GENES,
  drop = FALSE
]

train_means <- vapply(
  train_gene_matrix,
  mean,
  numeric(1),
  na.rm = TRUE
)

train_sds <- vapply(
  train_gene_matrix,
  stats::sd,
  numeric(1),
  na.rm = TRUE
)

if (
  any(
    !is.finite(train_sds) |
    train_sds == 0
  )
) {
  
  bad_genes <- names(
    train_sds
  )[
    !is.finite(train_sds) |
      train_sds == 0
  ]
  
  stop(
    "Invalid training SD for: ",
    paste(
      bad_genes,
      collapse = ", "
    )
  )
}

train_z <- apply_training_scaling(
  train_gene_matrix,
  train_means,
  train_sds
)

valid_z <- apply_training_scaling(
  valid_gene_matrix,
  train_means,
  train_sds
)

external_z <- apply_training_scaling(
  external_gene_matrix,
  train_means,
  train_sds
)

# =============================================================================
# 11. FIT FIVE-GENE LOGISTIC PANEL — TRAINING ONLY
# =============================================================================

train_model_df <- data.frame(
  Outcome = ifelse(
    train_df$Class == "Tumor",
    1,
    0
  ),
  train_z,
  check.names = FALSE
)

panel_formula <- stats::as.formula(
  paste(
    "Outcome ~",
    paste(
      GENES,
      collapse = " + "
    )
  )
)

panel_fit <- suppressWarnings(
  stats::glm(
    formula = panel_formula,
    data = train_model_df,
    family = stats::binomial(
      link = "logit"
    )
  )
)

if (!isTRUE(panel_fit$converged)) {
  
  stop(
    "Five-gene logistic-regression panel did not converge."
  )
}

# Predictions:
# NO retraining in Validation or External

panel_prob_train <- as.numeric(
  stats::predict(
    panel_fit,
    newdata = train_z,
    type = "response"
  )
)

panel_prob_valid <- as.numeric(
  stats::predict(
    panel_fit,
    newdata = valid_z,
    type = "response"
  )
)

panel_prob_external <- as.numeric(
  stats::predict(
    panel_fit,
    newdata = external_z,
    type = "response"
  )
)

# =============================================================================
# 12. BUILD MODEL SCORES
# =============================================================================

scores <- list(
  Training = list(),
  Validation = list(),
  External_GSE130605 = list()
)

for (gene in GENES) {
  
  scores$Training[[gene]] <-
    train_df[[gene]]
  
  scores$Validation[[gene]] <-
    valid_df[[gene]]
  
  scores$External_GSE130605[[gene]] <-
    external_df[[gene]]
}

scores$Training[[PANEL_KEY]] <-
  panel_prob_train

scores$Validation[[PANEL_KEY]] <-
  panel_prob_valid

scores$External_GSE130605[[PANEL_KEY]] <-
  panel_prob_external

truths <- list(
  Training = train_df$Class,
  Validation = valid_df$Class,
  External_GSE130605 = external_df$Class
)

directions <- c(
  gene_directions,
  FiveGenePanel = "<"
)

model_labels <- c(
  setNames(
    GENES,
    GENES
  ),
  FiveGenePanel = PANEL_LABEL
)

model_types <- c(
  setNames(
    rep(
      "Single gene",
      length(GENES)
    ),
    GENES
  ),
  FiveGenePanel = "Multivariable logistic panel"
)

# =============================================================================
# 13. AUROC + 95% CI
#     ALL FIVE GENES + FIVE-GENE PANEL
# =============================================================================

auc_rows <- list()

row_index <- 1L

for (cohort_name in names(scores)) {
  
  for (model_key in names(scores[[cohort_name]])) {
    
    auc_rows[[row_index]] <- get_auc_ci(
      cohort_name = cohort_name,
      model_key = model_key,
      model_label = unname(
        model_labels[
          model_key
        ]
      ),
      model_type = unname(
        model_types[
          model_key
        ]
      ),
      response = truths[[cohort_name]],
      predictor = scores[[cohort_name]][[model_key]],
      direction = unname(
        directions[
          model_key
        ]
      )
    )
    
    row_index <- row_index + 1L
  }
}

auc_table <- do.call(
  rbind,
  auc_rows
)

# Cohort ordering

cohort_order <- c(
  "Training",
  "Validation",
  "External_GSE130605"
)

model_order <- c(
  "MYBL2",
  "KRT17",
  "UCHL1",
  "PTHLH",
  "ODC1",
  "FiveGenePanel"
)

auc_table$Cohort_Order <- match(
  auc_table$Cohort,
  cohort_order
)

auc_table$Model_Order <- match(
  auc_table$Model_Key,
  model_order
)

auc_table <- auc_table[
  order(
    auc_table$Cohort_Order,
    auc_table$Model_Order
  ),
  ,
  drop = FALSE
]

auc_table$Cohort_Order <- NULL
auc_table$Model_Order <- NULL

auc_table$AUROC_95CI <- ifelse(
  is.finite(auc_table$AUROC),
  sprintf(
    "%.4f (%.4f–%.4f)",
    auc_table$AUROC,
    auc_table$AUROC_CI95_Lower,
    auc_table$AUROC_CI95_Upper
  ),
  NA_character_
)

# =============================================================================
# 14. PAIRED TWO-SIDED DELONG
#     ONLY FIVE-GENE PANEL vs MYBL2
# =============================================================================

delong_rows <- list()

for (i in seq_along(cohort_order)) {
  
  cohort_name <- cohort_order[i]
  
  y_all <- truths[[cohort_name]]
  
  panel_score <- scores[[cohort_name]][[
    PANEL_KEY
  ]]
  
  mybl2_score <- scores[[cohort_name]][[
    TARGET_GENE
  ]]
  
  ok <- (
    !is.na(y_all) &
      is.finite(panel_score) &
      is.finite(mybl2_score)
  )
  
  y <- y_all[ok]
  panel_x <- panel_score[ok]
  mybl2_x <- mybl2_score[ok]
  
  roc_panel <- pROC::roc(
    response = y,
    predictor = panel_x,
    levels = c(
      "Normal",
      "Tumor"
    ),
    direction = "<",
    quiet = TRUE
  )
  
  roc_mybl2 <- pROC::roc(
    response = y,
    predictor = mybl2_x,
    levels = c(
      "Normal",
      "Tumor"
    ),
    direction = unname(
      gene_directions[
        TARGET_GENE
      ]
    ),
    quiet = TRUE
  )
  
  delong_test <- pROC::roc.test(
    roc_panel,
    roc_mybl2,
    method = "delong",
    paired = TRUE,
    alternative = "two.sided"
  )
  
  panel_auc <- as.numeric(
    pROC::auc(
      roc_panel
    )
  )
  
  mybl2_auc <- as.numeric(
    pROC::auc(
      roc_mybl2
    )
  )
  
  z_value <- if (
    !is.null(
      delong_test$statistic
    )
  ) {
    
    as.numeric(
      delong_test$statistic[1]
    )
    
  } else {
    
    NA_real_
  }
  
  p_value <- as.numeric(
    delong_test$p.value
  )
  
  delong_rows[[i]] <- data.frame(
    Cohort = cohort_name,
    N_Paired = length(y),
    Panel_AUROC = panel_auc,
    MYBL2_AUROC = mybl2_auc,
    Delta_AUROC_Panel_minus_MYBL2 =
      panel_auc - mybl2_auc,
    DeLong_Z = z_value,
    DeLong_P_Value = p_value,
    Significant_P_lt_0.05 =
      p_value < 0.05,
    stringsAsFactors = FALSE
  )
}

delong_table <- do.call(
  rbind,
  delong_rows
)

# =============================================================================
# 15. ADD DELONG RESULTS TO THE SUPPLEMENTARY TABLE
# =============================================================================

auc_table$Panel_vs_MYBL2_Delta_AUROC <- NA_real_
auc_table$Panel_vs_MYBL2_DeLong_Z <- NA_real_
auc_table$Panel_vs_MYBL2_DeLong_P <- NA_real_
auc_table$Panel_vs_MYBL2_Significant <- NA

for (i in seq_len(nrow(delong_table))) {
  
  cohort_name <- delong_table$Cohort[i]
  
  idx <- (
    auc_table$Cohort == cohort_name &
      auc_table$Model_Key == PANEL_KEY
  )
  
  auc_table$Panel_vs_MYBL2_Delta_AUROC[idx] <-
    delong_table$Delta_AUROC_Panel_minus_MYBL2[i]
  
  auc_table$Panel_vs_MYBL2_DeLong_Z[idx] <-
    delong_table$DeLong_Z[i]
  
  auc_table$Panel_vs_MYBL2_DeLong_P[idx] <-
    delong_table$DeLong_P_Value[i]
  
  auc_table$Panel_vs_MYBL2_Significant[idx] <-
    delong_table$Significant_P_lt_0.05[i]
}

# =============================================================================
# 16. SAVE BENCHMARK TABLE
# =============================================================================

SUPPLEMENTARY_TABLE_FILE <- file.path(
  OUT_DIR,
  "Supplementary_Table_S5_AUROC_and_DeLong.csv"
)

write.csv(
  auc_table,
  SUPPLEMENTARY_TABLE_FILE,
  row.names = FALSE
)

# =============================================================================
# 17. ROC FIGURE
# =============================================================================
#
# Final layout:
#   A = Training
#   B = Validation
#   C = External GSE130605
# =============================================================================

model_colors <- c(
  FiveGenePanel = "#111111",
  MYBL2         = "#D81B60",
  KRT17         = "#4C78A8",
  PTHLH         = "#F4A261",
  UCHL1         = "#2A9D8F",
  ODC1          = "#8E6C8A"
)

plot_model_order <- c(
  "FiveGenePanel",
  "MYBL2",
  "KRT17",
  "PTHLH",
  "UCHL1",
  "ODC1"
)

legend_labels <- c(
  FiveGenePanel = "Five-gene panel",
  MYBL2         = "MYBL2",
  KRT17         = "KRT17",
  PTHLH         = "PTHLH",
  UCHL1         = "UCHL1",
  ODC1          = "ODC1"
)


# =============================================================================
# 18. ROC PLOT FUNCTION
# =============================================================================

make_roc_plot_publication <- function(
    cohort_name,
    title_text,
    panel_tag) {
  
  plot_rows <- list()
  
  auc_values <- setNames(
    rep(
      NA_real_,
      length(plot_model_order)
    ),
    plot_model_order
  )
  
  row_index <- 1L
  
  
  # ---------------------------------------------------------------------------
  # Build ROC curves
  # ---------------------------------------------------------------------------
  
  for (model_key in plot_model_order) {
    
    rr <- make_roc(
      response = truths[[cohort_name]],
      predictor = scores[[cohort_name]][[model_key]],
      direction = unname(
        directions[model_key]
      )
    )
    
    if (is.null(rr)) {
      next
    }
    
    auc_values[model_key] <- as.numeric(
      pROC::auc(rr)
    )
    
    plot_rows[[row_index]] <- data.frame(
      
      False_Positive_Rate =
        1 - rev(
          rr$specificities
        ),
      
      Sensitivity =
        rev(
          rr$sensitivities
        ),
      
      Model_Key = model_key,
      
      stringsAsFactors = FALSE
    )
    
    row_index <- row_index + 1L
  }
  
  
  if (length(plot_rows) == 0) {
    
    stop(
      "No valid ROC curve could be generated for cohort: ",
      cohort_name
    )
  }
  
  
  plot_df <- do.call(
    rbind,
    plot_rows
  )
  
  plot_df$Model_Key <- factor(
    plot_df$Model_Key,
    levels = plot_model_order
  )
  
  
  # ---------------------------------------------------------------------------
  # AUROC annotation
  # ---------------------------------------------------------------------------
  
  available_models <- plot_model_order[
    is.finite(
      auc_values[
        plot_model_order
      ]
    )
  ]
  
  auc_annotation <- paste(
    
    paste0(
      legend_labels[available_models],
      ": ",
      sprintf(
        "%.3f",
        auc_values[available_models]
      )
    ),
    
    collapse = "\n"
  )
  
  
  # ---------------------------------------------------------------------------
  # Draw panel
  # ---------------------------------------------------------------------------
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = False_Positive_Rate,
      y = Sensitivity,
      color = Model_Key,
      group = Model_Key
    )
  ) +
    
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      linewidth = 0.50,
      color = "grey65"
    ) +
    
    ggplot2::geom_line(
      linewidth = 1.35,
      lineend = "round"
    ) +
    
    ggplot2::scale_color_manual(
      values = model_colors,
      breaks = plot_model_order,
      labels = legend_labels[plot_model_order],
      drop = FALSE
    ) +
    
    ggplot2::scale_x_continuous(
      breaks = seq(
        0,
        1,
        by = 0.25
      ),
      labels = sprintf(
        "%.2f",
        seq(
          0,
          1,
          by = 0.25
        )
      ),
      expand = c(0, 0)
    ) +
    
    ggplot2::scale_y_continuous(
      breaks = seq(
        0,
        1,
        by = 0.25
      ),
      labels = sprintf(
        "%.2f",
        seq(
          0,
          1,
          by = 0.25
        )
      ),
      expand = c(0, 0)
    ) +
    
    # A true square ROC panel.
    ggplot2::coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    
    # AUROC box placed in the relatively empty bottom-right region.
    ggplot2::annotate(
      geom = "label",
      x = 0.965,
      y = 0.035,
      label = paste0(
        "AUROC\n",
        auc_annotation
      ),
      hjust = 1,
      vjust = 0,
      size = 3.15,
      lineheight = 1.05,
      label.size = 0.25,
      label.padding = grid::unit(
        0.14,
        "lines"
      ),
      fill = "white",
      alpha = 0.96,
      color = "black"
    ) +
    
    ggplot2::labs(
      title = title_text,
      tag = panel_tag,
      x = "1 - Specificity",
      y = "Sensitivity",
      color = NULL
    ) +
    
    ggplot2::theme_bw(
      base_size = 13
    ) +
    
    ggplot2::theme(
      
      panel.grid = ggplot2::element_blank(),
      
      panel.border = ggplot2::element_rect(
        color = "black",
        linewidth = 0.75,
        fill = NA
      ),
      
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 16,
        hjust = 0.5,
        margin = ggplot2::margin(
          b = 5
        )
      ),
      
      # Put A/B/C at the upper-left of each panel.
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 16,
        hjust = 0,
        vjust = 1
      ),
      
      plot.tag.position = c(
        0.015,
        0.985
      ),
      
      axis.title.x = ggplot2::element_text(
        size = 13.5,
        margin = ggplot2::margin(
          t = 5
        )
      ),
      
      axis.title.y = ggplot2::element_text(
        size = 13.5,
        margin = ggplot2::margin(
          r = 5
        )
      ),
      
      axis.text.x = ggplot2::element_text(
        size = 11
      ),
      
      axis.text.y = ggplot2::element_text(
        size = 11
      ),
      
      axis.ticks = ggplot2::element_line(
        linewidth = 0.50
      ),
      
      legend.position = "bottom",
      
      legend.text = ggplot2::element_text(
        size = 12
      ),
      
      legend.key.width = grid::unit(
        1.20,
        "cm"
      ),
      
      legend.key.height = grid::unit(
        0.35,
        "cm"
      ),
      
      legend.spacing.x = grid::unit(
        0.25,
        "cm"
      ),
      
      # Very small margins make the ROC area occupy most of each panel.
      plot.margin = ggplot2::margin(
        t = 4,
        r = 4,
        b = 4,
        l = 4
      )
    ) +
    
    ggplot2::guides(
      color = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          linewidth = 2.0
        )
      )
    )
}


# =============================================================================
# 19. CREATE THE THREE ROC PANELS
# =============================================================================

p_train <- make_roc_plot_publication(
  cohort_name = "Training",
  title_text = "Training",
  panel_tag = "A"
)

p_valid <- make_roc_plot_publication(
  cohort_name = "Validation",
  title_text = "Validation",
  panel_tag = "B"
)

p_external <- make_roc_plot_publication(
  cohort_name = "External_GSE130605",
  title_text = "External GSE130605",
  panel_tag = "C"
)


# =============================================================================
# 20. HORIZONTAL THREE-PANEL LAYOUT
# =============================================================================
#
combined_plot <- patchwork::wrap_plots(
  p_train,
  p_valid,
  p_external,
  ncol = 3,
  guides = "collect"
) &
  ggplot2::theme(
    legend.position = "bottom"
  )


# =============================================================================
# 21. DISPLAY FIGURE
# =============================================================================

print(
  combined_plot
)


# =============================================================================
# 22. SAVE FIGURE
# =============================================================================

FIGURE_STEM <-
  "Supplementary_Figure_S4_LASSO_Genes_and_FiveGene_Panel_HORIZONTAL"

FIG_WIDTH  <- 18.5
FIG_HEIGHT <- 6.6


ggplot2::ggsave(
  filename = file.path(
    OUT_DIR,
    paste0(
      FIGURE_STEM,
      ".png"
    )
  ),
  plot = combined_plot,
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)


ggplot2::ggsave(
  filename = file.path(
    OUT_DIR,
    paste0(
      FIGURE_STEM,
      ".tiff"
    )
  ),
  plot = combined_plot,
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white",
  limitsize = FALSE
)


tryCatch(
  {
    
    if (capabilities("cairo")) {
      
      ggplot2::ggsave(
        filename = file.path(
          OUT_DIR,
          paste0(
            FIGURE_STEM,
            ".pdf"
          )
        ),
        plot = combined_plot,
        width = FIG_WIDTH,
        height = FIG_HEIGHT,
        units = "in",
        device = grDevices::cairo_pdf,
        bg = "white",
        limitsize = FALSE
      )
      
    } else {
      
      ggplot2::ggsave(
        filename = file.path(
          OUT_DIR,
          paste0(
            FIGURE_STEM,
            ".pdf"
          )
        ),
        plot = combined_plot,
        width = FIG_WIDTH,
        height = FIG_HEIGHT,
        units = "in",
        device = "pdf",
        bg = "white",
        limitsize = FALSE
      )
    }
    
  },
  error = function(e) {
    
    warning(
      "PDF export failed: ",
      conditionMessage(e)
    )
  }
)


# =============================================================================
# 23. CONSOLE RESULTS
# =============================================================================

cat("\n\n")
cat("============================================================\n")
cat("SUPPLEMENTARY TABLE S5 — AUROC + 95% CI\n")
cat("============================================================\n\n")

print(
  auc_table[
    ,
    c(
      "Cohort",
      "Model",
      "Model_Type",
      "N",
      "Normal_N",
      "Tumor_N",
      "AUROC",
      "AUROC_CI95_Lower",
      "AUROC_CI95_Upper",
      "AUROC_95CI"
    )
  ],
  row.names = FALSE
)


cat("\n\n")
cat("============================================================\n")
cat("FIVE-GENE PANEL vs MYBL2 — PAIRED TWO-SIDED DELONG\n")
cat("============================================================\n\n")

print(
  delong_table,
  row.names = FALSE
)


cat("\n\n")
cat("============================================================\n")
cat("KEY INTERPRETATION\n")
cat("============================================================\n\n")

for (i in seq_len(nrow(delong_table))) {
  
  x <- delong_table[i, ]
  
  cat(
    x$Cohort,
    ": Panel AUC = ",
    sprintf(
      "%.4f",
      x$Panel_AUROC
    ),
    "; MYBL2 AUC = ",
    sprintf(
      "%.4f",
      x$MYBL2_AUROC
    ),
    "; Delta = ",
    sprintf(
      "%.4f",
      x$Delta_AUROC_Panel_minus_MYBL2
    ),
    "; DeLong P = ",
    format(
      x$DeLong_P_Value,
      digits = 4
    ),
    "; Significant = ",
    ifelse(
      x$Significant_P_lt_0.05,
      "YES",
      "NO"
    ),
    "\n",
    sep = ""
  )
}


cat("\n")
cat("============================================================\n")
cat("ANALYSIS FINISHED\n")
cat("============================================================\n")

cat(
  "\nOutput folder:\n",
  OUT_DIR,
  "\n",
  sep = ""
)

cat(
  "\nRequired supplementary table:\n",
  SUPPLEMENTARY_TABLE_FILE,
  "\n",
  sep = ""
)

cat(
  "\nFinal horizontal supplementary figure:\n",
  file.path(
    OUT_DIR,
    paste0(
      FIGURE_STEM,
      ".png/.tiff/.pdf"
    )
  ),
  "\n",
  sep = ""
)

cat(
  "\nFigure dimensions: ",
  FIG_WIDTH,
  " x ",
  FIG_HEIGHT,
  " inches\n",
  sep = ""
)

cat("============================================================\n")
