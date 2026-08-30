###############################################################################
# ESTABLISHED-BIOMARKER BENCHMARK
# LSCC MYBL2 reproducible analysis
#
# Purpose:
#   Benchmark MYBL2 against established/reported markers:
#     - CDKN2A (p16)
#     - EGFR
#     - MKI67 (Ki-67)
#
# Cohorts:
#   - Training
#   - Internal validation
#   - External GSE130605
#
# Analyses:
#   1) Single-gene ROC/AUC with 95% CI in each cohort
#   2) Paired two-sided DeLong tests: MYBL2 vs each comparator
#   3) Tumor-vs-normal Wilcoxon tests for expression
#   4) ROC overlay figures for each cohort
#   5) Combined 3-panel ROC figure
#   6) Summary tables
#
#   - ROC direction is learned ONLY from the training cohort for each gene and
#     then kept fixed in validation/external cohorts to avoid data leakage.
###############################################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(scipen = 999)
set.seed(123)

# =============================================================================
# 1. SETTINGS
# =============================================================================

RESULTS_PATH <- "E:/LSCC/Results_LSCC/ML"

RUN_NAME <- "Established_Biomarker_Benchmark"

OUT_DIR <- file.path(
  RESULTS_PATH,
  RUN_NAME
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

MARKERS <- c(
  "MYBL2",
  "CDKN2A",
  "EGFR",
  "MKI67"
)

MARKER_LABELS <- c(
  MYBL2 = "MYBL2",
  CDKN2A = "p16 (CDKN2A)",
  EGFR = "EGFR",
  MKI67 = "Ki-67 (MKI67)"
)

TARGET_GENE <- "MYBL2"

FIG_DPI <- 600

# =============================================================================
# 2. REQUIRED PACKAGES
# =============================================================================

required_packages <- c(
  "data.table",
  "dplyr",
  "ggplot2",
  "pROC",
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
  
  cat(
    "\nInstalling missing CRAN package(s): ",
    paste(missing_packages, collapse = ", "),
    "\n",
    sep = ""
  )
  
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}

still_missing <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(still_missing) > 0) {
  stop(
    "Could not install/load required package(s): ",
    paste(still_missing, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(pROC)
  library(patchwork)
})

# =============================================================================
# 3. HELPERS
# =============================================================================

safe_numeric <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}

make_group_factor <- function(g) {
  
  g0 <- as.character(g)
  
  if (all(g0 %in% c("1", "2"))) {
    return(
      factor(
        g0,
        levels = c("1", "2"),
        labels = c("Normal", "Tumor")
      )
    )
  }
  
  gl <- tolower(trimws(g0))
  
  if (all(gl %in% c(
    "normal", "tumor", "non", "lscc",
    "cancer", "margin"
  ))) {
    
    lab <- ifelse(
      gl %in% c("tumor", "lscc", "cancer"),
      "Tumor",
      "Normal"
    )
    
    return(
      factor(
        lab,
        levels = c("Normal", "Tumor")
      )
    )
  }
  
  stop(
    "Could not convert group labels. Found: ",
    paste(unique(g0), collapse = ", ")
  )
}

format_p <- function(p) {
  
  if (is.na(p)) {
    return("NA")
  }
  
  if (p < 0.0001) {
    return(
      format(
        p,
        scientific = TRUE,
        digits = 3
      )
    )
  }
  
  sprintf("%.4f", p)
}

save_plot_all_formats <- function(
    plot_obj,
    stem,
    width = 7.2,
    height = 5.4) {
  
  png_file <- file.path(
    OUT_DIR,
    paste0(stem, ".png")
  )
  
  tiff_file <- file.path(
    OUT_DIR,
    paste0(stem, ".tiff")
  )
  
  pdf_file <- file.path(
    OUT_DIR,
    paste0(stem, ".pdf")
  )
  
  ggplot2::ggsave(
    filename = png_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = FIG_DPI,
    bg = "white",
    limitsize = FALSE
  )
  
  ggplot2::ggsave(
    filename = tiff_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = FIG_DPI,
    compression = "lzw",
    bg = "white",
    limitsize = FALSE
  )
  
  tryCatch(
    {
      if (capabilities("cairo")) {
        
        ggplot2::ggsave(
          filename = pdf_file,
          plot = plot_obj +
            ggplot2::theme(
              text = ggplot2::element_text(
                family = "sans"
              )
            ),
          width = width,
          height = height,
          units = "in",
          device = grDevices::cairo_pdf,
          bg = "white",
          limitsize = FALSE
        )
        
      } else {
        
        ggplot2::ggsave(
          filename = pdf_file,
          plot = plot_obj +
            ggplot2::theme(
              text = ggplot2::element_text(
                family = "Helvetica"
              )
            ),
          width = width,
          height = height,
          units = "in",
          device = "pdf",
          family = "Helvetica",
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
# 4. CHECK INPUT FILES
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
    "Required input file(s) not found:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

cat("\n============================================================\n")
cat("LSCC BIOMARKER BENCHMARKING\n")
cat("============================================================\n")
cat("Output folder:\n", OUT_DIR, "\n", sep = "")

# =============================================================================
# 5. LOAD COHORTS
# =============================================================================

train_df <- as.data.frame(
  data.table::fread(
    TRAIN_FILE
  ),
  check.names = FALSE
)

valid_df <- as.data.frame(
  data.table::fread(
    VALID_FILE
  ),
  check.names = FALSE
)

external_df <- as.data.frame(
  data.table::fread(
    EXTERNAL_FILE
  ),
  check.names = FALSE
)

cohorts <- list(
  Training = train_df,
  Validation = valid_df,
  External_GSE130605 = external_df
)

for (nm in names(cohorts)) {
  
  dat <- cohorts[[nm]]
  
  if (!("group" %in% colnames(dat))) {
    stop(
      "group column not found in cohort: ",
      nm
    )
  }
  
  cohorts[[nm]]$Class <- make_group_factor(
    dat$group
  )
}

# =============================================================================
# 6. CHECK MARKER AVAILABILITY
# =============================================================================

availability_table <- do.call(
  rbind,
  lapply(
    names(cohorts),
    function(cohort_name) {
      
      dat <- cohorts[[cohort_name]]
      
      data.frame(
        Cohort = cohort_name,
        Gene = MARKERS,
        Marker = unname(
          MARKER_LABELS[MARKERS]
        ),
        Available = MARKERS %in%
          colnames(dat),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  availability_table,
  file.path(
    OUT_DIR,
    "Biomarker_Availability_by_Cohort.csv"
  ),
  row.names = FALSE
)

cat("\nMarker availability:\n")
print(
  availability_table,
  row.names = FALSE
)

if (!all(
  availability_table$Available[
    availability_table$Gene == TARGET_GENE
  ]
)) {
  stop(
    TARGET_GENE,
    " is not available in all three cohorts."
  )
}

# Only benchmark markers available in ALL cohorts.
available_all <- availability_table %>%
  dplyr::group_by(
    Gene
  ) %>%
  dplyr::summarise(
    Available_All = all(Available),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    Available_All
  ) %>%
  dplyr::pull(
    Gene
  )

if (length(available_all) < 2) {
  stop(
    "Fewer than two benchmark markers are available in all cohorts."
  )
}

cat(
  "\nMarkers available in all cohorts: ",
  paste(
    available_all,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

# =============================================================================
# 7. LEARN ROC DIRECTION FROM TRAINING ONLY
# =============================================================================

# pROC direction:
# "<" means controls tend to have lower predictor values than cases.
# ">" means controls tend to have higher predictor values than cases.
#
# Direction is selected once in Training and then fixed in Validation/External.

direction_table <- data.frame(
  Gene = available_all,
  Marker = unname(
    MARKER_LABELS[available_all]
  ),
  Training_Normal_Median = NA_real_,
  Training_Tumor_Median = NA_real_,
  ROC_Direction = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(available_all)) {
  
  gene <- available_all[i]
  
  x <- safe_numeric(
    train_df[[gene]]
  )
  
  y <- cohorts$Training$Class
  
  med_normal <- median(
    x[y == "Normal"],
    na.rm = TRUE
  )
  
  med_tumor <- median(
    x[y == "Tumor"],
    na.rm = TRUE
  )
  
  direction_table$Training_Normal_Median[i] <- med_normal
  direction_table$Training_Tumor_Median[i] <- med_tumor
  
  direction_table$ROC_Direction[i] <- ifelse(
    med_tumor >= med_normal,
    "<",
    ">"
  )
}

write.csv(
  direction_table,
  file.path(
    OUT_DIR,
    "ROC_Direction_Learned_From_Training.csv"
  ),
  row.names = FALSE
)

cat("\nROC directions learned from Training:\n")
print(
  direction_table,
  row.names = FALSE
)

# =============================================================================
# 8. ROC/AUC + 95% CI FOR EACH MARKER AND COHORT
# =============================================================================

roc_objects <- list()
auc_rows <- list()

row_index <- 1L

for (cohort_name in names(cohorts)) {
  
  dat <- cohorts[[cohort_name]]
  
  roc_objects[[cohort_name]] <- list()
  
  for (gene in available_all) {
    
    y <- dat$Class
    
    x <- safe_numeric(
      dat[[gene]]
    )
    
    ok <- (
      !is.na(y) &
        is.finite(x)
    )
    
    y2 <- y[ok]
    x2 <- x[ok]
    
    direction_now <- direction_table$ROC_Direction[
      direction_table$Gene == gene
    ][1]
    
    if (
      length(unique(y2)) < 2 ||
      length(x2) < 4
    ) {
      
      auc_rows[[row_index]] <- data.frame(
        Cohort = cohort_name,
        Gene = gene,
        Marker = unname(
          MARKER_LABELS[gene]
        ),
        N = length(x2),
        Normal_N = sum(y2 == "Normal"),
        Tumor_N = sum(y2 == "Tumor"),
        ROC_Direction = direction_now,
        AUROC = NA_real_,
        AUROC_CI95_Lower = NA_real_,
        AUROC_CI95_Upper = NA_real_,
        stringsAsFactors = FALSE
      )
      
      row_index <- row_index + 1L
      next
    }
    
    roc_obj <- pROC::roc(
      response = y2,
      predictor = x2,
      levels = c(
        "Normal",
        "Tumor"
      ),
      direction = direction_now,
      quiet = TRUE
    )
    
    roc_objects[[cohort_name]][[gene]] <- roc_obj
    
    ci_auc <- as.numeric(
      pROC::ci.auc(
        roc_obj,
        method = "delong",
        conf.level = 0.95
      )
    )
    
    auc_rows[[row_index]] <- data.frame(
      Cohort = cohort_name,
      Gene = gene,
      Marker = unname(
        MARKER_LABELS[gene]
      ),
      N = length(y2),
      Normal_N = sum(
        y2 == "Normal"
      ),
      Tumor_N = sum(
        y2 == "Tumor"
      ),
      ROC_Direction = direction_now,
      AUROC = as.numeric(
        pROC::auc(
          roc_obj
        )
      ),
      AUROC_CI95_Lower = ci_auc[1],
      AUROC_CI95_Upper = ci_auc[3],
      stringsAsFactors = FALSE
    )
    
    row_index <- row_index + 1L
  }
}

auc_table <- dplyr::bind_rows(
  auc_rows
) %>%
  dplyr::arrange(
    factor(
      Cohort,
      levels = names(cohorts)
    ),
    dplyr::desc(
      AUROC
    )
  )

write.csv(
  auc_table,
  file.path(
    OUT_DIR,
    "Biomarker_Benchmark_AUROC_95CI.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 9. PAIRED DELONG: MYBL2 VS EACH COMPARATOR
# =============================================================================

delong_rows <- list()
k <- 1L

for (cohort_name in names(cohorts)) {
  
  dat <- cohorts[[cohort_name]]
  
  comparators <- setdiff(
    available_all,
    TARGET_GENE
  )
  
  for (gene in comparators) {
    
    y <- dat$Class
    x_target <- safe_numeric(dat[[TARGET_GENE]])
    x_comp <- safe_numeric(dat[[gene]])
    
    # Paired DeLong requires both markers to be evaluated on exactly
    # the same subjects. Therefore use pairwise complete observations.
    ok_pair <- (
      !is.na(y) &
        is.finite(x_target) &
        is.finite(x_comp)
    )
    
    y_pair <- y[ok_pair]
    target_pair <- x_target[ok_pair]
    comp_pair <- x_comp[ok_pair]
    
    if (
      length(y_pair) < 4 ||
      length(unique(y_pair)) < 2
    ) {
      delong_rows[[k]] <- data.frame(
        Cohort = cohort_name,
        N_Paired = length(y_pair),
        Reference_Gene = TARGET_GENE,
        Reference_Marker = unname(MARKER_LABELS[TARGET_GENE]),
        Reference_AUROC = NA_real_,
        Comparator_Gene = gene,
        Comparator_Marker = unname(MARKER_LABELS[gene]),
        Comparator_AUROC = NA_real_,
        AUROC_Difference_MYBL2_minus_Comparator = NA_real_,
        DeLong_Test = "Paired, two-sided",
        DeLong_Z = NA_real_,
        DeLong_P_Value = NA_real_,
        Significant_P_lt_0.05 = NA,
        stringsAsFactors = FALSE
      )
      
      k <- k + 1L
      next
    }
    
    target_direction <- direction_table$ROC_Direction[
      direction_table$Gene == TARGET_GENE
    ][1]
    
    comparator_direction <- direction_table$ROC_Direction[
      direction_table$Gene == gene
    ][1]
    
    target_roc_pair <- pROC::roc(
      response = y_pair,
      predictor = target_pair,
      levels = c("Normal", "Tumor"),
      direction = target_direction,
      quiet = TRUE
    )
    
    comparator_roc_pair <- pROC::roc(
      response = y_pair,
      predictor = comp_pair,
      levels = c("Normal", "Tumor"),
      direction = comparator_direction,
      quiet = TRUE
    )
    
    target_auc <- as.numeric(
      pROC::auc(target_roc_pair)
    )
    
    comparator_auc <- as.numeric(
      pROC::auc(comparator_roc_pair)
    )
    
    dt <- tryCatch(
      pROC::roc.test(
        target_roc_pair,
        comparator_roc_pair,
        method = "delong",
        paired = TRUE,
        alternative = "two.sided"
      ),
      error = function(e) {
        warning(
          "DeLong test failed for ",
          cohort_name,
          " | MYBL2 vs ",
          gene,
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
    
    if (is.null(dt)) {
      zval <- NA_real_
      pval <- NA_real_
    } else {
      zval <- if (!is.null(dt$statistic)) {
        as.numeric(dt$statistic[1])
      } else {
        NA_real_
      }
      
      pval <- as.numeric(dt$p.value)
    }
    
    delong_rows[[k]] <- data.frame(
      Cohort = cohort_name,
      N_Paired = length(y_pair),
      Reference_Gene = TARGET_GENE,
      Reference_Marker = unname(
        MARKER_LABELS[TARGET_GENE]
      ),
      Reference_AUROC = target_auc,
      Comparator_Gene = gene,
      Comparator_Marker = unname(
        MARKER_LABELS[gene]
      ),
      Comparator_AUROC = comparator_auc,
      AUROC_Difference_MYBL2_minus_Comparator =
        target_auc - comparator_auc,
      DeLong_Test = "Paired, two-sided",
      DeLong_Z = zval,
      DeLong_P_Value = pval,
      Significant_P_lt_0.05 =
        ifelse(
          is.na(pval),
          NA,
          pval < 0.05
        ),
      stringsAsFactors = FALSE
    )
    
    k <- k + 1L
  }
}

delong_table <- dplyr::bind_rows(
  delong_rows
)

write.csv(
  delong_table,
  file.path(
    OUT_DIR,
    "Biomarker_Benchmark_DeLong_MYBL2_vs_Comparators.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 10. TUMOR-vs-NORMAL EXPRESSION COMPARISON
# =============================================================================

wilcox_rows <- list()
k <- 1L

for (cohort_name in names(cohorts)) {
  
  dat <- cohorts[[cohort_name]]
  
  for (gene in available_all) {
    
    x <- safe_numeric(
      dat[[gene]]
    )
    
    y <- dat$Class
    
    ok <- (
      !is.na(y) &
        is.finite(x)
    )
    
    x2 <- x[ok]
    y2 <- y[ok]
    
    tumor_vals <- x2[
      y2 == "Tumor"
    ]
    
    normal_vals <- x2[
      y2 == "Normal"
    ]
    
    if (
      length(tumor_vals) < 2 ||
      length(normal_vals) < 2
    ) {
      
      wilcox_rows[[k]] <- data.frame(
        Cohort = cohort_name,
        Gene = gene,
        Marker = unname(
          MARKER_LABELS[gene]
        ),
        Tumor_N = length(
          tumor_vals
        ),
        Normal_N = length(
          normal_vals
        ),
        Tumor_Median = NA_real_,
        Normal_Median = NA_real_,
        Median_Difference_Tumor_minus_Normal = NA_real_,
        Wilcoxon_W = NA_real_,
        P_Value = NA_real_,
        stringsAsFactors = FALSE
      )
      
      k <- k + 1L
      next
    }
    
    wt <- stats::wilcox.test(
      tumor_vals,
      normal_vals,
      alternative = "two.sided",
      exact = FALSE
    )
    
    med_t <- median(
      tumor_vals,
      na.rm = TRUE
    )
    
    med_n <- median(
      normal_vals,
      na.rm = TRUE
    )
    
    wilcox_rows[[k]] <- data.frame(
      Cohort = cohort_name,
      Gene = gene,
      Marker = unname(
        MARKER_LABELS[gene]
      ),
      Tumor_N = length(
        tumor_vals
      ),
      Normal_N = length(
        normal_vals
      ),
      Tumor_Median = med_t,
      Normal_Median = med_n,
      Median_Difference_Tumor_minus_Normal =
        med_t - med_n,
      Wilcoxon_W = as.numeric(
        wt$statistic
      ),
      P_Value = as.numeric(
        wt$p.value
      ),
      stringsAsFactors = FALSE
    )
    
    k <- k + 1L
  }
}

wilcox_table <- dplyr::bind_rows(
  wilcox_rows
)

write.csv(
  wilcox_table,
  file.path(
    OUT_DIR,
    "Biomarker_Benchmark_Expression_Wilcoxon.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 11. ROC FIGURES — FINAL HORIZONTAL 3-PANEL PUBLICATION VERSION
# =============================================================================
#
# Figure layout matched to the requested reference image:
#   A = Training
#   B = Internal validation
#   C = External GSE130605
#
# Main figure-format settings:
#   - All three ROC panels placed SIDE BY SIDE
#   - True 1:1 ROC aspect ratio in every panel
#   - Wide landscape canvas (18.5 x 6.6 inches)
#   - Minimal outer margins so panels appear large / zoomed
#   - One shared legend below all three panels
#   - Cohort-specific AUROC values shown in a compact box at bottom-right
#   - Consistent colors across cohorts
#   - 600-dpi PNG/TIFF + vector PDF output
# =============================================================================

# Fixed colors for readability and consistency across all cohorts.
marker_colors <- c(
  MYBL2  = "#D81B60",
  CDKN2A = "#4C78A8",
  EGFR   = "#2A9D8F",
  MKI67  = "#F4A261"
)

# Keep the same plotting/legend order in every cohort.
plot_marker_order <- c(
  "MYBL2",
  "CDKN2A",
  "EGFR",
  "MKI67"
)

plot_marker_order <- plot_marker_order[
  plot_marker_order %in% available_all
]

make_roc_plot <- function(
    cohort_name,
    title_text,
    panel_tag) {
  
  rocs_now <- roc_objects[[cohort_name]]
  
  rocs_now <- rocs_now[
    names(rocs_now) %in% plot_marker_order
  ]
  
  if (length(rocs_now) < 1) {
    return(NULL)
  }
  
  plot_rows <- list()
  auc_values <- setNames(
    rep(NA_real_, length(plot_marker_order)),
    plot_marker_order
  )
  
  j <- 1L
  
  for (gene in plot_marker_order) {
    
    rr <- rocs_now[[gene]]
    
    if (is.null(rr)) {
      next
    }
    
    auc_now <- as.numeric(
      pROC::auc(rr)
    )
    
    auc_values[gene] <- auc_now
    
    coords_df <- data.frame(
      False_Positive_Rate = 1 - rev(
        rr$specificities
      ),
      Sensitivity = rev(
        rr$sensitivities
      ),
      Gene = gene,
      stringsAsFactors = FALSE
    )
    
    plot_rows[[j]] <- coords_df
    j <- j + 1L
  }
  
  plot_df <- dplyr::bind_rows(
    plot_rows
  )
  
  if (nrow(plot_df) == 0) {
    return(NULL)
  }
  
  # Lock color/legend order.
  plot_df$Gene <- factor(
    plot_df$Gene,
    levels = plot_marker_order
  )
  
  available_plot_markers <- plot_marker_order[
    is.finite(
      auc_values[plot_marker_order]
    )
  ]
  
  # Compact AUROC annotation inside the lower-right region of each ROC panel.
  auc_annotation <- paste(
    paste0(
      unname(
        MARKER_LABELS[available_plot_markers]
      ),
      ": ",
      sprintf(
        "%.3f",
        auc_values[available_plot_markers]
      )
    ),
    collapse = "\n"
  )
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = False_Positive_Rate,
      y = Sensitivity,
      color = Gene,
      group = Gene
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
      values = marker_colors,
      breaks = plot_marker_order,
      labels = unname(
        MARKER_LABELS[plot_marker_order]
      ),
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
    # Keep each ROC panel perfectly square.
    ggplot2::coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
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
      
      # A / B / C at the upper-left, matching the reference layout.
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
      
      legend.margin = ggplot2::margin(
        t = 2,
        r = 0,
        b = 0,
        l = 0
      ),
      
      # Minimal margins maximize the visible ROC plotting area.
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

# -----------------------------------------------------------------------------
# Create the three labeled panels
# -----------------------------------------------------------------------------

p_train <- make_roc_plot(
  cohort_name = "Training",
  title_text = "Training",
  panel_tag = "A"
)

p_valid <- make_roc_plot(
  cohort_name = "Validation",
  title_text = "Validation",
  panel_tag = "B"
)

p_external <- make_roc_plot(
  cohort_name = "External_GSE130605",
  title_text = "External GSE130605",
  panel_tag = "C"
)

# -----------------------------------------------------------------------------
# Save individual cohort figures
# -----------------------------------------------------------------------------

if (!is.null(p_train)) {
  save_plot_all_formats(
    p_train,
    "Figure_Benchmark_ROC_Training",
    width = 7.2,
    height = 6.4
  )
}

if (!is.null(p_valid)) {
  save_plot_all_formats(
    p_valid,
    "Figure_Benchmark_ROC_Validation",
    width = 7.2,
    height = 6.4
  )
}

if (!is.null(p_external)) {
  save_plot_all_formats(
    p_external,
    "Figure_Benchmark_ROC_External_GSE130605",
    width = 7.2,
    height = 6.4
  )
}

# -----------------------------------------------------------------------------
# FINAL HORIZONTAL 3-PANEL LAYOUT
# -----------------------------------------------------------------------------
# ncol = 3 places Training, Validation, and External side by side.
# guides = "collect" creates one shared legend beneath the full figure.
# -----------------------------------------------------------------------------

if (
  !is.null(p_train) &&
  !is.null(p_valid) &&
  !is.null(p_external)
) {
  
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
  
  # Display combined figure in the R graphics device.
  print(
    combined_plot
  )
  
  # Wide landscape dimensions matched to the reference image proportions.
  FIGURE_BENCHMARK_WIDTH  <- 18.5
  FIGURE_BENCHMARK_HEIGHT <- 6.6
  
  save_plot_all_formats(
    combined_plot,
    "Figure_Benchmark_ROC_All_Cohorts_HORIZONTAL",
    width = FIGURE_BENCHMARK_WIDTH,
    height = FIGURE_BENCHMARK_HEIGHT
  )
}


# =============================================================================
# 12. RANKING TABLE
# =============================================================================

ranking_table <- auc_table %>%
  dplyr::group_by(
    Cohort
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      AUROC
    ),
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    AUROC_Rank = dplyr::row_number()
  ) %>%
  dplyr::ungroup()

write.csv(
  ranking_table,
  file.path(
    OUT_DIR,
    "Biomarker_Benchmark_Ranking.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 13. TEXT SUMMARY
# =============================================================================

summary_lines <- c(
  "LSCC BIOMARKER BENCHMARKING",
  "",
  paste0(
    "Markers evaluated: ",
    paste(
      unname(
        MARKER_LABELS[
          available_all
        ]
      ),
      collapse = ", "
    )
  ),
  "",
  "ROC directions were determined from the training cohort only and then fixed for validation and external testing.",
  "",
  "AUROC RESULTS:"
)

for (cohort_name in names(cohorts)) {
  
  sub_auc <- auc_table %>%
    dplyr::filter(
      Cohort == cohort_name
    )
  
  summary_lines <- c(
    summary_lines,
    paste0(
      "\n",
      cohort_name,
      ":"
    )
  )
  
  for (i in seq_len(
    nrow(sub_auc)
  )) {
    
    summary_lines <- c(
      summary_lines,
      paste0(
        "  ",
        sub_auc$Marker[i],
        " AUROC = ",
        ifelse(
          is.na(
            sub_auc$AUROC[i]
          ),
          "NA",
          sprintf(
            "%.4f",
            sub_auc$AUROC[i]
          )
        ),
        " (95% CI ",
        ifelse(
          is.na(
            sub_auc$AUROC_CI95_Lower[i]
          ),
          "NA",
          sprintf(
            "%.4f",
            sub_auc$AUROC_CI95_Lower[i]
          )
        ),
        "-",
        ifelse(
          is.na(
            sub_auc$AUROC_CI95_Upper[i]
          ),
          "NA",
          sprintf(
            "%.4f",
            sub_auc$AUROC_CI95_Upper[i]
          )
        ),
        ")"
      )
    )
  }
}

summary_lines <- c(
  summary_lines,
  "",
  "PAIRED DELONG: MYBL2 VS COMPARATORS"
)

if (nrow(delong_table) > 0) {
  
  for (i in seq_len(
    nrow(delong_table)
  )) {
    
    summary_lines <- c(
      summary_lines,
      paste0(
        delong_table$Cohort[i],
        ": MYBL2 vs ",
        delong_table$Comparator_Marker[i],
        " | AUC difference = ",
        sprintf(
          "%.4f",
          delong_table$AUROC_Difference_MYBL2_minus_Comparator[i]
        ),
        " | Z = ",
        ifelse(
          is.na(
            delong_table$DeLong_Z[i]
          ),
          "NA",
          sprintf(
            "%.4f",
            delong_table$DeLong_Z[i]
          )
        ),
        " | P = ",
        format_p(
          delong_table$DeLong_P_Value[i]
        )
      )
    )
  }
}

writeLines(
  summary_lines,
  file.path(
    OUT_DIR,
    "Biomarker_Benchmarking_Summary.txt"
  )
)

# =============================================================================
# 14. CONSOLE SUMMARY
# =============================================================================

cat("\n============================================================\n")
cat("BENCHMARKING FINISHED\n")
cat("============================================================\n")

cat("\nAUROC TABLE:\n")
print(
  auc_table,
  row.names = FALSE
)

cat("\nPAIRED DELONG RESULTS:\n")

if (nrow(delong_table) > 0) {
  print(
    delong_table,
    row.names = FALSE
  )
} else {
  cat("No DeLong comparisons were available.\n")
}

cat("\nWILCOXON EXPRESSION RESULTS:\n")
print(
  wilcox_table,
  row.names = FALSE
)

cat("\nRANKING:\n")
print(
  ranking_table[
    ,
    c(
      "Cohort",
      "Marker",
      "AUROC",
      "AUROC_Rank"
    )
  ],
  row.names = FALSE
)

figure_files <- list.files(
  OUT_DIR,
  pattern = "\\.(png|tiff|pdf)$",
  ignore.case = TRUE
)

cat("\nFigure files created:", length(figure_files), "\n")

if (length(figure_files) > 0) {
  cat(
    paste0(
      " - ",
      figure_files
    ),
    sep = "\n"
  )
  cat("\n")
}

cat("\nOutputs saved to:\n")
cat(
  OUT_DIR,
  "\n"
)

cat("\nKey files:\n")
cat(" - Biomarker_Benchmark_AUROC_95CI.csv\n")
cat(" - Biomarker_Benchmark_DeLong_MYBL2_vs_Comparators.csv\n")
cat(" - Biomarker_Benchmark_Expression_Wilcoxon.csv\n")
cat(" - Biomarker_Benchmark_Ranking.csv\n")
cat(" - Figure_Benchmark_ROC_All_Cohorts.png/tiff/pdf\n")
cat(" - Biomarker_Benchmarking_Summary.txt\n")

cat("\n============================================================\n")
cat("END\n")
cat("============================================================\n")

# End of script
