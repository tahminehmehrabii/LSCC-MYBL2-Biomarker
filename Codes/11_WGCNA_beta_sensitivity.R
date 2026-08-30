
# LSCC WGCNA soft-threshold sensitivity analysis
# Evaluates network stability across alternative beta values.

required_packages <- c("WGCNA", "ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

cat("\n============================================================\n")
cat("WGCNA BETA SENSITIVITY ANALYSIS\n")
cat("============================================================\n")

# Step 1: settings

BETA_SENSITIVITY_POWERS <- c(12, 14, 16, 18, 20)
BETA_SENSITIVITY_REFERENCE <- 20
BETA_SENSITIVITY_SEED <- 123
BETA_TARGET_GENE <- "MYBL2"

# Keep the exact settings of the main WGCNA analysis.
if (!exists("WGCNA_TOP_MAD", inherits = FALSE)) WGCNA_TOP_MAD <- 8000
if (!exists("WGCNA_MM_MIN", inherits = FALSE)) WGCNA_MM_MIN <- 0.55
if (!exists("WGCNA_GS_MIN", inherits = FALSE)) WGCNA_GS_MIN <- 0.55
if (!exists("WGCNA_RULE", inherits = FALSE)) WGCNA_RULE <- "AND"
if (!exists("WGCNA_MIN_MEANK", inherits = FALSE)) WGCNA_MIN_MEANK <- 10

# Main module construction settings from Step 04.
BETA_MIN_MODULE_SIZE <- 60
BETA_MERGE_CUT_HEIGHT <- 0.18
BETA_NETWORK_TYPE <- "signed"
BETA_TOM_TYPE <- "signed"
BETA_COR_TYPE <- "pearson"
BETA_PAM_RESPECTS_DENDRO <- TRUE

# Output directory.
if (!exists("step04_dir", inherits = FALSE)) {
  if (exists("results_path", inherits = FALSE)) {
    step04_dir <- results_path
  } else {
    stop("step04_dir is not available. Run Step 04 first or define step04_dir.")
  }
}

beta_sensitivity_dir <- file.path(step04_dir, "Beta_Sensitivity_Analysis")
dir.create(beta_sensitivity_dir, recursive = TRUE, showWarnings = FALSE)

# Step 2: load the exact wgcna input used in the main analysis

# Prefer in-memory objects from Step 04. If they are absent, load the saved RDS.
if (!exists("datExpr", inherits = FALSE) || !exists("traitDF", inherits = FALSE)) {
  wgcna_rds <- file.path(step04_dir, "WGCNA_analysis_objects.rds")
  if (!file.exists(wgcna_rds)) {
    stop(
      "Could not find datExpr/traitDF in memory and the saved WGCNA_analysis_objects.rds is missing: ",
      wgcna_rds
    )
  }
  wgcna_main_obj <- readRDS(wgcna_rds)
  datExpr <- wgcna_main_obj$datExpr
  traitDF <- wgcna_main_obj$traitDF
  if (is.null(datExpr) || is.null(traitDF)) {
    stop("The saved WGCNA_analysis_objects.rds does not contain datExpr and traitDF.")
  }
}

# Ensure sample order is identical.
traitDF <- traitDF[rownames(datExpr), , drop = FALSE]
if (!"Tumor" %in% colnames(traitDF)) {
  stop("traitDF must contain the Tumor trait.")
}

cat("Samples used:", nrow(datExpr), "\n")
cat("Genes used:", ncol(datExpr), "\n")
cat("Powers tested:", paste(BETA_SENSITIVITY_POWERS, collapse = ", "), "\n")
cat("Reference power:", BETA_SENSITIVITY_REFERENCE, "\n")
cat("Target gene:", BETA_TARGET_GENE, "\n")

# Step 3: scale-free topology information for the tested powers

# Reuse main Step 04 soft-threshold results if available; otherwise recompute.
if (exists("sft_table", inherits = FALSE) &&
    all(c("Power", "SignedR2", "MeanK") %in% colnames(sft_table))) {
  beta_sft_table <- as.data.frame(sft_table)
} else {
  cat("\nRecomputing pickSoftThreshold for powers 1:20...\n")
  sft_beta <- WGCNA::pickSoftThreshold(
    datExpr,
    powerVector = 1:20,
    networkType = BETA_NETWORK_TYPE,
    corFnc = "cor",
    verbose = 0
  )
  beta_fitR2 <- -sign(sft_beta$fitIndices[, 3]) * sft_beta$fitIndices[, 2]
  beta_meanK <- sft_beta$fitIndices[, 5]
  beta_sft_table <- data.frame(
    Power = sft_beta$fitIndices[, 1],
    SignedR2 = beta_fitR2,
    MeanK = beta_meanK,
    stringsAsFactors = FALSE
  )
}

beta_sft_tested <- beta_sft_table[
  beta_sft_table$Power %in% BETA_SENSITIVITY_POWERS,
  c("Power", "SignedR2", "MeanK"),
  drop = FALSE
]

write.csv(
  beta_sft_tested,
  file.path(beta_sensitivity_dir, "WGCNA_beta_scale_free_topology.csv"),
  row.names = FALSE
)

# Step 4: helpers

safe_jaccard <- function(a, b) {
  a <- unique(as.character(a))
  b <- unique(as.character(b))
  a <- a[nzchar(a) & !is.na(a)]
  b <- b[nzchar(b) & !is.na(b)]
  u <- union(a, b)
  if (length(u) == 0) return(NA_real_)
  length(intersect(a, b)) / length(u)
}

safe_overlap_pct_ref <- function(current, ref) {
  current <- unique(as.character(current))
  ref <- unique(as.character(ref))
  current <- current[nzchar(current) & !is.na(current)]
  ref <- ref[nzchar(ref) & !is.na(ref)]
  if (length(ref) == 0) return(NA_real_)
  100 * length(intersect(current, ref)) / length(ref)
}

# Create a minimal plotting theme if the global manuscript theme is unavailable.
beta_theme <- function() {
  if (exists("theme_manuscript", mode = "function", inherits = TRUE)) {
    return(theme_manuscript(show_grid = FALSE, legend_position = "bottom"))
  }
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold")
    )
}

# Step 5: run one wgcna analysis at a forced beta

run_beta_wgcna <- function(beta_value) {

  cat("\n------------------------------------------------------------\n")
  cat("Running WGCNA sensitivity at beta =", beta_value, "\n")
  cat("------------------------------------------------------------\n")

  set.seed(BETA_SENSITIVITY_SEED)

  net_b <- WGCNA::blockwiseModules(
    datExpr,
    power = beta_value,
    networkType = BETA_NETWORK_TYPE,
    TOMType = BETA_TOM_TYPE,
    corType = BETA_COR_TYPE,
    minModuleSize = BETA_MIN_MODULE_SIZE,
    mergeCutHeight = BETA_MERGE_CUT_HEIGHT,
    numericLabels = FALSE,
    pamRespectsDendro = BETA_PAM_RESPECTS_DENDRO,
    maxBlockSize = ncol(datExpr),
    verbose = 0
  )

  module_colors_b <- net_b$colors
  names(module_colors_b) <- colnames(datExpr)

  MEs_b <- WGCNA::orderMEs(net_b$MEs)

  modTraitCor_b <- cor(
    MEs_b,
    traitDF,
    use = "p",
    method = "pearson"
  )

  modTraitP_b <- WGCNA::corPvalueStudent(
    modTraitCor_b,
    nrow(datExpr)
  )

  # All module sizes, including grey for completeness.
  module_size_b <- as.data.frame(table(module_colors_b), stringsAsFactors = FALSE)
  colnames(module_size_b) <- c("Module", "Module_Size")
  module_size_b$Module_Size <- as.integer(module_size_b$Module_Size)

  # Long table of ALL module-trait associations.
  module_trait_long_b <- do.call(
    rbind,
    lapply(rownames(modTraitCor_b), function(me_name) {
      module_name <- sub("^ME", "", me_name)
      size_now <- module_size_b$Module_Size[match(module_name, module_size_b$Module)]
      do.call(
        rbind,
        lapply(colnames(modTraitCor_b), function(trait_name) {
          data.frame(
            Beta = beta_value,
            Module = module_name,
            ME = me_name,
            Trait = trait_name,
            Correlation = as.numeric(modTraitCor_b[me_name, trait_name]),
            P_value = as.numeric(modTraitP_b[me_name, trait_name]),
            Module_Size = ifelse(length(size_now) == 1, size_now, NA_integer_),
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )

  # Find strongest positively tumor-associated non-grey module.
  valid_me_b <- rownames(modTraitCor_b)
  valid_me_b <- valid_me_b[valid_me_b != "MEgrey"]

  if (length(valid_me_b) < 1) {
    stop("No non-grey modules were detected at beta = ", beta_value)
  }

  tumor_cor_b <- modTraitCor_b[valid_me_b, "Tumor"]
  tumor_cor_pos_b <- tumor_cor_b
  tumor_cor_pos_b[tumor_cor_pos_b <= 0] <- NA_real_

  if (all(is.na(tumor_cor_pos_b))) {
    stop("No positively tumor-associated module at beta = ", beta_value)
  }

  bestME_b <- valid_me_b[which.max(tumor_cor_pos_b)]
  targetModule_b <- sub("^ME", "", bestME_b)
  targetGenes_b <- names(module_colors_b)[module_colors_b == targetModule_b]

  best_tumor_cor_b <- as.numeric(modTraitCor_b[bestME_b, "Tumor"])
  best_tumor_p_b <- as.numeric(modTraitP_b[bestME_b, "Tumor"])

  # Gene significance and module membership for the selected tumor module.
  GS_b <- as.numeric(cor(datExpr, traitDF$Tumor, use = "p", method = "pearson"))
  MM_b <- as.numeric(cor(datExpr, MEs_b[, bestME_b], use = "p", method = "pearson"))
  names(GS_b) <- colnames(datExpr)
  names(MM_b) <- colnames(datExpr)

  feature_table_b <- data.frame(
    Gene = targetGenes_b,
    GS_Tumor = GS_b[targetGenes_b],
    MM_Module = MM_b[targetGenes_b],
    stringsAsFactors = FALSE
  )
  feature_table_b$absGS <- abs(feature_table_b$GS_Tumor)
  feature_table_b$absMM <- abs(feature_table_b$MM_Module)

  if (toupper(WGCNA_RULE) == "AND") {
    selected_features_b <- feature_table_b$Gene[
      feature_table_b$absMM >= WGCNA_MM_MIN &
        feature_table_b$absGS >= WGCNA_GS_MIN
    ]
  } else {
    selected_features_b <- feature_table_b$Gene[
      feature_table_b$absMM >= WGCNA_MM_MIN |
        feature_table_b$absGS >= WGCNA_GS_MIN
    ]
  }
  selected_features_b <- sort(unique(selected_features_b))

  # Module-size metrics for over-aggregation concern.
  nongrey_sizes_b <- module_size_b[module_size_b$Module != "grey", , drop = FALSE]
  n_modules_b <- nrow(nongrey_sizes_b)
  largest_module_size_b <- if (n_modules_b > 0) max(nongrey_sizes_b$Module_Size) else NA_integer_
  median_module_size_b <- if (n_modules_b > 0) stats::median(nongrey_sizes_b$Module_Size) else NA_real_
  target_module_size_b <- length(targetGenes_b)

  # MYBL2-specific metrics.
  mybl2_in_input_b <- BETA_TARGET_GENE %in% colnames(datExpr)
  mybl2_module_b <- if (mybl2_in_input_b) unname(module_colors_b[BETA_TARGET_GENE]) else NA_character_
  mybl2_module_size_b <- if (mybl2_in_input_b && !is.na(mybl2_module_b)) {
    as.integer(sum(module_colors_b == mybl2_module_b))
  } else {
    NA_integer_
  }
  mybl2_gs_b <- if (mybl2_in_input_b) unname(GS_b[BETA_TARGET_GENE]) else NA_real_
  mybl2_mm_target_b <- if (mybl2_in_input_b) unname(MM_b[BETA_TARGET_GENE]) else NA_real_
  mybl2_in_target_b <- mybl2_in_input_b && (BETA_TARGET_GENE %in% targetGenes_b)
  mybl2_selected_b <- mybl2_in_input_b && (BETA_TARGET_GENE %in% selected_features_b)

  # Scale-free topology values for this beta.
  sft_row_b <- beta_sft_table[beta_sft_table$Power == beta_value, , drop = FALSE]
  signed_r2_b <- if (nrow(sft_row_b) > 0) sft_row_b$SignedR2[1] else NA_real_
  mean_k_b <- if (nrow(sft_row_b) > 0) sft_row_b$MeanK[1] else NA_real_

  summary_b <- data.frame(
    Beta = beta_value,
    Signed_R2 = signed_r2_b,
    Mean_Connectivity = mean_k_b,
    N_Modules_NonGrey = n_modules_b,
    Largest_Module_Size = largest_module_size_b,
    Median_Module_Size = median_module_size_b,
    Best_Tumor_Module = targetModule_b,
    Best_Tumor_Correlation = best_tumor_cor_b,
    Best_Tumor_P = best_tumor_p_b,
    Best_Tumor_Module_Size = target_module_size_b,
    Selected_Feature_Count = length(selected_features_b),
    MYBL2_in_WGCNA_Input = mybl2_in_input_b,
    MYBL2_Module = mybl2_module_b,
    MYBL2_Module_Size = mybl2_module_size_b,
    MYBL2_GS_Tumor = mybl2_gs_b,
    MYBL2_MM_to_Best_Tumor_Module = mybl2_mm_target_b,
    MYBL2_in_Best_Tumor_Module = mybl2_in_target_b,
    MYBL2_Selected_Feature = mybl2_selected_b,
    stringsAsFactors = FALSE
  )

  # Save per-beta detailed tables.
  write.csv(
    module_trait_long_b,
    file.path(beta_sensitivity_dir, paste0("beta_", beta_value, "_all_module_trait_associations.csv")),
    row.names = FALSE
  )

  write.csv(
    module_size_b,
    file.path(beta_sensitivity_dir, paste0("beta_", beta_value, "_module_sizes.csv")),
    row.names = FALSE
  )

  write.csv(
    feature_table_b,
    file.path(beta_sensitivity_dir, paste0("beta_", beta_value, "_best_tumor_module_GS_MM.csv")),
    row.names = FALSE
  )

  write.csv(
    data.frame(Gene = selected_features_b, stringsAsFactors = FALSE),
    file.path(beta_sensitivity_dir, paste0("beta_", beta_value, "_selected_WGCNA_features.csv")),
    row.names = FALSE
  )

  # Return compact objects. Avoid retaining full TOM/net objects to control RAM.
  out <- list(
    summary = summary_b,
    target_genes = targetGenes_b,
    selected_features = selected_features_b,
    module_trait = module_trait_long_b,
    module_sizes = module_size_b
  )

  rm(net_b, MEs_b, modTraitCor_b, modTraitP_b, feature_table_b)
  gc(verbose = FALSE)

  out
}

# Step 6: run all powers sequentially

beta_runs <- lapply(BETA_SENSITIVITY_POWERS, function(b) {
  tryCatch(
    run_beta_wgcna(b),
    error = function(e) {
      cat("beta", b, "FAILED:", conditionMessage(e), "\n")
      list(
        summary = data.frame(
          Beta = b,
          Signed_R2 = NA_real_,
          Mean_Connectivity = NA_real_,
          N_Modules_NonGrey = NA_integer_,
          Largest_Module_Size = NA_integer_,
          Median_Module_Size = NA_real_,
          Best_Tumor_Module = NA_character_,
          Best_Tumor_Correlation = NA_real_,
          Best_Tumor_P = NA_real_,
          Best_Tumor_Module_Size = NA_integer_,
          Selected_Feature_Count = NA_integer_,
          MYBL2_in_WGCNA_Input = BETA_TARGET_GENE %in% colnames(datExpr),
          MYBL2_Module = NA_character_,
          MYBL2_Module_Size = NA_integer_,
          MYBL2_GS_Tumor = NA_real_,
          MYBL2_MM_to_Best_Tumor_Module = NA_real_,
          MYBL2_in_Best_Tumor_Module = NA,
          MYBL2_Selected_Feature = NA,
          Error_Message = conditionMessage(e),
          stringsAsFactors = FALSE
        ),
        target_genes = character(0),
        selected_features = character(0),
        module_trait = data.frame(),
        module_sizes = data.frame()
      )
    }
  )
})

beta_summary <- dplyr::bind_rows(lapply(beta_runs, `[[`, "summary"))
beta_module_trait_all <- dplyr::bind_rows(lapply(beta_runs, `[[`, "module_trait"))

# Step 7: compare each power with beta = 20

ref_index <- match(BETA_SENSITIVITY_REFERENCE, BETA_SENSITIVITY_POWERS)
if (is.na(ref_index)) {
  stop("Reference beta must be included in BETA_SENSITIVITY_POWERS.")
}

ref_target_genes <- beta_runs[[ref_index]]$target_genes
ref_selected_features <- beta_runs[[ref_index]]$selected_features

beta_overlap <- do.call(
  rbind,
  lapply(seq_along(BETA_SENSITIVITY_POWERS), function(i) {
    b <- BETA_SENSITIVITY_POWERS[i]
    current_target <- beta_runs[[i]]$target_genes
    current_features <- beta_runs[[i]]$selected_features

    data.frame(
      Beta = b,
      Target_Module_Overlap_Count_vs_Beta20 = length(intersect(current_target, ref_target_genes)),
      Target_Module_Jaccard_vs_Beta20 = safe_jaccard(current_target, ref_target_genes),
      Target_Module_RefRecovery_pct = safe_overlap_pct_ref(current_target, ref_target_genes),
      Feature_Overlap_Count_vs_Beta20 = length(intersect(current_features, ref_selected_features)),
      Feature_Jaccard_vs_Beta20 = safe_jaccard(current_features, ref_selected_features),
      Feature_RefRecovery_pct = safe_overlap_pct_ref(current_features, ref_selected_features),
      MYBL2_in_Target_Module = BETA_TARGET_GENE %in% current_target,
      MYBL2_in_Selected_Features = BETA_TARGET_GENE %in% current_features,
      stringsAsFactors = FALSE
    )
  })
)

beta_summary_final <- dplyr::left_join(beta_summary, beta_overlap, by = "Beta")

write.csv(
  beta_summary_final,
  file.path(beta_sensitivity_dir, "WGCNA_beta_sensitivity_summary.csv"),
  row.names = FALSE
)

write.csv(
  beta_module_trait_all,
  file.path(beta_sensitivity_dir, "WGCNA_beta_all_module_trait_associations_long.csv"),
  row.names = FALSE
)

write.csv(
  beta_overlap,
  file.path(beta_sensitivity_dir, "WGCNA_beta_overlap_vs_beta20.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    settings = list(
      powers = BETA_SENSITIVITY_POWERS,
      reference_beta = BETA_SENSITIVITY_REFERENCE,
      seed = BETA_SENSITIVITY_SEED,
      minModuleSize = BETA_MIN_MODULE_SIZE,
      mergeCutHeight = BETA_MERGE_CUT_HEIGHT,
      networkType = BETA_NETWORK_TYPE,
      TOMType = BETA_TOM_TYPE,
      corType = BETA_COR_TYPE,
      MM_threshold = WGCNA_MM_MIN,
      GS_threshold = WGCNA_GS_MIN,
      selection_rule = WGCNA_RULE
    ),
    summary = beta_summary_final,
    runs = beta_runs
  ),
  file.path(beta_sensitivity_dir, "WGCNA_beta_sensitivity_objects.rds")
)

# Step 8: sensitivity figure

plot_df <- beta_summary_final

p_beta_cor <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = Beta, y = Best_Tumor_Correlation)
) +
  ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
  ggplot2::geom_point(size = 2.7, na.rm = TRUE) +
  ggplot2::scale_x_continuous(breaks = BETA_SENSITIVITY_POWERS) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Tumor-module association",
    x = expression(beta),
    y = "Strongest positive tumor correlation"
  ) +
  beta_theme()

p_beta_size <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = Beta, y = Best_Tumor_Module_Size)
) +
  ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
  ggplot2::geom_point(size = 2.7, na.rm = TRUE) +
  ggplot2::scale_x_continuous(breaks = BETA_SENSITIVITY_POWERS) +
  ggplot2::labs(
    title = "Size of strongest tumor-associated module",
    x = expression(beta),
    y = "Number of genes"
  ) +
  beta_theme()

p_beta_jaccard <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = Beta, y = Feature_Jaccard_vs_Beta20)
) +
  ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
  ggplot2::geom_point(size = 2.7, na.rm = TRUE) +
  ggplot2::scale_x_continuous(breaks = BETA_SENSITIVITY_POWERS) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Feature-gene overlap with beta = 20",
    x = expression(beta),
    y = "Jaccard index"
  ) +
  beta_theme()

mybl2_plot_df <- data.frame(
  Beta = plot_df$Beta,
  Status = ifelse(
    plot_df$MYBL2_Selected_Feature %in% TRUE,
    "Selected WGCNA feature",
    ifelse(
      plot_df$MYBL2_in_Best_Tumor_Module %in% TRUE,
      "In best tumor module",
      "Not in best tumor module"
    )
  ),
  stringsAsFactors = FALSE
)

mybl2_plot_df$Status <- factor(
  mybl2_plot_df$Status,
  levels = c(
    "Selected WGCNA feature",
    "In best tumor module",
    "Not in best tumor module"
  )
)

p_beta_mybl2 <- ggplot2::ggplot(
  mybl2_plot_df,
  ggplot2::aes(x = Beta, y = 1, shape = Status)
) +
  ggplot2::geom_point(size = 4, na.rm = TRUE) +
  ggplot2::scale_x_continuous(breaks = BETA_SENSITIVITY_POWERS) +
  ggplot2::scale_y_continuous(NULL, breaks = NULL) +
  ggplot2::labs(
    title = paste0(BETA_TARGET_GENE, " robustness"),
    x = expression(beta),
    shape = NULL
  ) +
  beta_theme()

if (requireNamespace("patchwork", quietly = TRUE)) {
  p_beta_sensitivity <- (
    p_beta_cor + p_beta_size + p_beta_jaccard + p_beta_mybl2
  ) + patchwork::plot_layout(ncol = 2)
} else {
  warning("Package 'patchwork' is not installed. Saving the four panels separately.")
  p_beta_sensitivity <- NULL
}

# FIXED EXPORT FUNCTION: PNG + TIFF + PDF
# Fixes: "failed to find or load PDF CID font"

save_ggplot_three_formats <- function(plot_obj, stem, width, height) {

  if (is.null(plot_obj)) {
    return(invisible(NULL))
  }

  dpi_use <- if (exists("FIG_DPI", inherits = TRUE)) {
    FIG_DPI
  } else {
    600
  }

  # -------------------------
  # PNG
  # -------------------------
  ggplot2::ggsave(
    filename = file.path(
      beta_sensitivity_dir,
      paste0(stem, ".png")
    ),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi_use,
    bg = "white",
    limitsize = FALSE
  )

  # -------------------------
  # TIFF
  # -------------------------
  ggplot2::ggsave(
    filename = file.path(
      beta_sensitivity_dir,
      paste0(stem, ".tiff")
    ),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi_use,
    compression = "lzw",
    bg = "white",
    limitsize = FALSE
  )

  # -------------------------
  # PDF
  # Use safe generic font for PDF only
  # -------------------------

  pdf_plot <- plot_obj

  safe_pdf_theme <- ggplot2::theme(
    text = ggplot2::element_text(
      family = "sans"
    ),
    plot.title = ggplot2::element_text(
      family = "sans",
      face = "bold"
    ),
    axis.title = ggplot2::element_text(
      family = "sans"
    ),
    axis.text = ggplot2::element_text(
      family = "sans"
    ),
    legend.title = ggplot2::element_text(
      family = "sans"
    ),
    legend.text = ggplot2::element_text(
      family = "sans"
    ),
    strip.text = ggplot2::element_text(
      family = "sans"
    )
  )

  # Apply font correction to all panels if plot is patchwork
  if (inherits(pdf_plot, "patchwork")) {
    pdf_plot <- pdf_plot & safe_pdf_theme
  } else {
    pdf_plot <- pdf_plot + safe_pdf_theme
  }

  pdf_file <- file.path(
    beta_sensitivity_dir,
    paste0(stem, ".pdf")
  )

  tryCatch({

    if (capabilities("cairo")) {

      ggplot2::ggsave(
        filename = pdf_file,
        plot = pdf_plot,
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
        plot = pdf_plot,
        width = width,
        height = height,
        units = "in",
        device = "pdf",
        family = "Helvetica",
        bg = "white",
        limitsize = FALSE
      )
    }

    cat("PDF successfully saved:\n", pdf_file, "\n")

  }, error = function(e) {

    warning(
      "PDF export failed, but PNG/TIFF were saved successfully: ",
      conditionMessage(e)
    )
  })

  invisible(
    list(
      PNG = file.path(
        beta_sensitivity_dir,
        paste0(stem, ".png")
      ),
      TIFF = file.path(
        beta_sensitivity_dir,
        paste0(stem, ".tiff")
      ),
      PDF = pdf_file
    )
  )
}

# Save sensitivity figure

save_ggplot_three_formats(
  plot_obj = p_beta_sensitivity,
  stem = "Figure_WGCNA_Beta_Sensitivity",
  width = 10,
  height = 8
)

cat("\nWGCNA beta sensitivity figure export finished.\n")

if (is.null(p_beta_sensitivity)) {
  save_ggplot_three_formats(p_beta_cor, "Figure_WGCNA_Beta_Sensitivity_A_TumorCorrelation", 5.5, 4.5)
  save_ggplot_three_formats(p_beta_size, "Figure_WGCNA_Beta_Sensitivity_B_ModuleSize", 5.5, 4.5)
  save_ggplot_three_formats(p_beta_jaccard, "Figure_WGCNA_Beta_Sensitivity_C_FeatureJaccard", 5.5, 4.5)
  save_ggplot_three_formats(p_beta_mybl2, "Figure_WGCNA_Beta_Sensitivity_D_MYBL2", 5.5, 4.5)
}

# Step 9: analysis summary

successful_beta_df <- beta_summary_final[is.finite(beta_summary_final$Best_Tumor_Correlation), , drop = FALSE]

beta_lines <- c(
  "WGCNA BETA SENSITIVITY ANALYSIS",
  paste0("Powers tested: ", paste(BETA_SENSITIVITY_POWERS, collapse = ", ")),
  paste0("Reference power: beta = ", BETA_SENSITIVITY_REFERENCE),
  paste0("Fixed random seed for WGCNA runs: ", BETA_SENSITIVITY_SEED),
  "",
  "All other WGCNA parameters were kept identical to the primary analysis:",
  paste0("networkType = ", BETA_NETWORK_TYPE),
  paste0("TOMType = ", BETA_TOM_TYPE),
  paste0("minModuleSize = ", BETA_MIN_MODULE_SIZE),
  paste0("mergeCutHeight = ", BETA_MERGE_CUT_HEIGHT),
  paste0("MM threshold = ", WGCNA_MM_MIN),
  paste0("GS threshold = ", WGCNA_GS_MIN),
  paste0("MM/GS rule = ", WGCNA_RULE),
  "",
  "Summary table:",
  capture.output(print(beta_summary_final, row.names = FALSE)),
  ""
)

if (nrow(successful_beta_df) > 0) {
  beta_lines <- c(
    beta_lines,
    paste0(
      "Strongest positive tumor-module correlations ranged from ",
      sprintf("%.3f", min(successful_beta_df$Best_Tumor_Correlation, na.rm = TRUE)),
      " to ",
      sprintf("%.3f", max(successful_beta_df$Best_Tumor_Correlation, na.rm = TRUE)),
      " across the tested powers."
    ),
    paste0(
      BETA_TARGET_GENE,
      " was retained as a selected WGCNA feature in ",
      sum(successful_beta_df$MYBL2_Selected_Feature %in% TRUE, na.rm = TRUE),
      "/",
      nrow(successful_beta_df),
      " successful beta runs."
    )
  )
}

writeLines(
  beta_lines,
  file.path(beta_sensitivity_dir, "WGCNA_beta_sensitivity_summary.txt")
)

# Step 10: console summary

cat("\n============================================================\n")
cat("WGCNA BETA SENSITIVITY FINISHED\n")
cat("============================================================\n")
print(beta_summary_final)
cat("\nOutputs saved to:\n", beta_sensitivity_dir, "\n", sep = "")
cat("\nKey output files:\n")
cat(" - WGCNA_beta_sensitivity_summary.csv\n")
cat(" - WGCNA_beta_all_module_trait_associations_long.csv\n")
cat(" - WGCNA_beta_overlap_vs_beta20.csv\n")
cat(" - Figure_WGCNA_Beta_Sensitivity.png/.tiff/.pdf\n")
cat(" - WGCNA_beta_sensitivity_summary.txt\n")

# END
