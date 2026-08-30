# LSCC bulk validation of Cluster 0 and Cluster 8 signatures
# Builds strict scRNA-seq-derived signatures and scores them with ssGSEA.

rm(list = ls(all.names = TRUE))
gc()

set.seed(123)

options(
  stringsAsFactors = FALSE,
  scipen = 100,
  timeout = 7200,
  width = 160
)

# Step 1: version check, paths, and analytical settings

if (as.character(getRversion()) != "4.4.3") {
  stop(
    "Run this script only in R 4.4.3. Current R version: ",
    as.character(getRversion())
  )
}

BASE_RESULTS_DIR <- "D:/LSCC/Results_LSCC"
BULK_ML_DIR <- file.path(BASE_RESULTS_DIR, "ML")

SC_PROJECT_DIR <- "D:/LSCC/ScRNAseq_Results/GSE206332"
SC_RESULTS_DIR <- file.path(SC_PROJECT_DIR, "Results")
SC_FIG_DIR <- file.path(SC_RESULTS_DIR, "figures")
SC_RDS_DIR <- file.path(SC_RESULTS_DIR, "rds")

OUTPUT_DIR <- file.path(
  BASE_RESULTS_DIR,
  "LSCC_ssGSEA_Clusters0_8"
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

TARGET_GENE <- "MYBL2"

# Only the two MYBL2-associated malignant subclusters are included.
SELECTED_SUBCLUSTERS <- c(
  "Cluster 0",
  "Cluster 8"
)

# Strict marker criteria, aligned with the High-CNV malignant marker pipeline.
MARKER_PADJ_CUTOFF <- 0.05
MARKER_LOGFC_CUTOFF <- 1.00
MARKER_PCT_CUTOFF <- 0.25

# Each cluster signature contains its top strict, unique markers after excluding
# MYBL2, mitochondrial genes, ribosomal genes, and markers shared by Clusters 0 and 8.
TOP_SIGNATURE_GENES <- 50L
MIN_SIGNATURE_GENES <- 10L

REMOVE_MITOCHONDRIAL_GENES <- TRUE
REMOVE_RIBOSOMAL_GENES <- TRUE

# Existing bulk matrices. These must be non-negative TMM-CPM matrices.
DEVELOPMENT_CPM_CANDIDATES <- c(
  file.path(
    BULK_ML_DIR,
    "train_discovery_merged_CPM_nonnegative.csv"
  ),
  file.path(
    BULK_ML_DIR,
    "train_discovery_merged_CPM_nonnegative_noZero.csv"
  )
)

EXTERNAL_CPM_CANDIDATES <- c(
  file.path(
    BULK_ML_DIR,
    "external_GSE130605_CPM.csv"
  )
)

MARKER_CSV_CANDIDATES <- c(
  file.path(
    SC_FIG_DIR,
    "High_CNV_Malignant_Cluster_FindAllMarkers_after_lowCNV_removed_logFC1_all.csv"
  ),
  file.path(
    SC_FIG_DIR,
    "High_CNV_Malignant_FindAllMarkers_after_lowCNV_removed_logFC1.csv"
  ),
  file.path(
    SC_FIG_DIR,
    "Final_High_CNV_Malignant_Genes.csv"
  ),
  file.path(
    SC_RDS_DIR,
    "High_CNV_Malignant_Cluster_FindAllMarkers_after_lowCNV_removed_logFC1_all.csv"
  )
)

# Output files
SIGNATURE_MARKERS_CSV <- file.path(
  OUTPUT_DIR,
  "Clusters0_8_ssGSEA_signature_genes.csv"
)

SIGNATURE_OVERLAP_CSV <- file.path(
  OUTPUT_DIR,
  "Clusters0_8_signature_gene_overlap.csv"
)

SIGNATURE_GMT <- file.path(
  OUTPUT_DIR,
  "Clusters0_8_ssGSEA_signatures.gmt"
)

SIGNATURE_STATUS_CSV <- file.path(
  OUTPUT_DIR,
  "Clusters0_8_signature_gene_overlap_in_bulk.csv"
)

DISCOVERY_SCORES_CSV <- file.path(
  OUTPUT_DIR,
  "Discovery_ssGSEA_scores_Clusters0_8.csv"
)

EXTERNAL_SCORES_CSV <- file.path(
  OUTPUT_DIR,
  "External_GSE130605_ssGSEA_scores_Clusters0_8.csv"
)

DISCOVERY_GROUP_TESTS_CSV <- file.path(
  OUTPUT_DIR,
  "Discovery_tumor_vs_normal_ssGSEA_tests_Clusters0_8.csv"
)

EXTERNAL_GROUP_TESTS_CSV <- file.path(
  OUTPUT_DIR,
  "External_GSE130605_tumor_vs_normal_ssGSEA_tests_Clusters0_8.csv"
)

DISCOVERY_CORRELATION_CSV <- file.path(
  OUTPUT_DIR,
  "Discovery_tumor_MYBL2_ssGSEA_Spearman_correlations_Clusters0_8.csv"
)

PROVENANCE_TXT <- file.path(
  OUTPUT_DIR,
  "Clusters0_8_ssGSEA_input_provenance.txt"
)

SUMMARY_TXT <- file.path(
  OUTPUT_DIR,
  "Clusters0_8_ssGSEA_summary.txt"
)

STATUS_TXT <- file.path(
  OUTPUT_DIR,
  "Clusters0_8_ssGSEA_status.txt"
)

SESSION_TXT <- file.path(
  OUTPUT_DIR,
  "Clusters0_8_ssGSEA_sessionInfo.txt"
)

RUN_LOG <- file.path(
  OUTPUT_DIR,
  "Clusters0_8_ssGSEA_run.log"
)

RESULT_RDS <- file.path(
  OUTPUT_DIR,
  "Clusters0_8_ssGSEA_analysis_objects.rds"
)

# The only figure produced by this script, exported as PNG, TIFF and PDF.
FINAL_FIGURE_STEM <- "Figure_08_cluster_signature_ssGSEA_validation"

# Step 2: manuscript figure settings

FONT_FAMILY <- "Arial"
FIG_DPI <- 600L
FIG_BACKGROUND <- "white"

# Three-row, two-column panel layout:
# a) discovery Tumor vs Normal
# b) external Tumor vs Normal
# c) discovery-tumor MYBL2 correlation
FIGURE_W <- 12.80
FIGURE_H <- 11.20

BASE_TEXT_PT <- 9.5
AXIS_TEXT_PT <- 8.5
AXIS_TITLE_PT <- 10
LEGEND_TEXT_PT <- 8.5
LEGEND_TITLE_PT <- 9
PLOT_TITLE_PT <- 10
PANEL_TAG_PT <- 14

PANEL_BORDER_LWD <- 0.45
AXIS_LWD <- 0.40
GEOM_LWD <- 0.70
GRID_LWD <- 0.30

# Chosen colorblind-friendly palette:
# Normal = teal, Tumor = burnt orange, Cluster 8 scatter = plum.
# Panel b uses a subtle blue-tinted background and blue border to denote
# independent external validation, while Normal/Tumor colours remain unchanged.
COL_NORMAL <- "#2A9D8F"
COL_TUMOR <- "#E76F51"
COL_CLUSTER0 <- "#287A78"
COL_CLUSTER8 <- "#825B9B"
COL_FIT <- "#303030"
COL_BAND <- "#D9D9D9"

COL_EXTERNAL_PANEL_BG <- "#F3F7FB"
COL_EXTERNAL_PANEL_BORDER <- "#4C78A8"

DISCOVERY_COHORT_HEADER <- "Discovery cohort (GSE127165 + GSE142083)"
EXTERNAL_COHORT_HEADER <- "Independent validation cohort (GSE130605)"

SIGNATURE_COLOURS <- c(
  "Cluster 0 signature" = COL_CLUSTER0,
  "Cluster 8 signature" = COL_CLUSTER8
)

MANUSCRIPT_MARGIN <- ggplot2::margin(
  t = 5.5,
  r = 6,
  b = 5.5,
  l = 6,
  unit = "pt"
)

if (.Platform$OS.type == "windows") {
  try(
    grDevices::windowsFonts(
      Arial = grDevices::windowsFont("Arial")
    ),
    silent = TRUE
  )
}

# Step 3: package check and loading

required_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "patchwork",
  "GSVA",
  "scales"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "These required packages are missing:\n",
    paste(missing_packages, collapse = ", "),
    "\n\nInstall them in R 4.4.3, restart RGui, and run this script again."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# Step 4: general helper functions

log_message <- function(...) {
  txt <- paste0(..., collapse = "")
  message(txt)

  cat(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ",
    txt,
    "\n",
    file = RUN_LOG,
    append = TRUE,
    sep = ""
  )
}

theme_manuscript <- function(show_grid = FALSE,
                             legend_position = "right",
                             panel_fill = FIG_BACKGROUND,
                             panel_border_colour = "black") {
  ggplot2::theme_bw(
    base_size = BASE_TEXT_PT,
    base_family = FONT_FAMILY
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = FONT_FAMILY,
        color = "black"
      ),
      plot.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PLOT_TITLE_PT,
        hjust = 0.5,
        color = "black"
      ),
      plot.tag = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PANEL_TAG_PT,
        color = "black"
      ),
      plot.tag.position = c(0.01, 0.99),
      axis.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = AXIS_TITLE_PT,
        color = "black"
      ),
      axis.text = ggplot2::element_text(
        family = FONT_FAMILY,
        size = AXIS_TEXT_PT,
        color = "black"
      ),
      axis.line = ggplot2::element_line(
        color = "black",
        linewidth = AXIS_LWD
      ),
      axis.ticks = ggplot2::element_line(
        color = "black",
        linewidth = AXIS_LWD
      ),
      axis.ticks.length = grid::unit(2, "pt"),
      legend.position = legend_position,
      legend.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = LEGEND_TITLE_PT
      ),
      legend.text = ggplot2::element_text(
        family = FONT_FAMILY,
        size = LEGEND_TEXT_PT
      ),
      legend.background = ggplot2::element_rect(
        fill = FIG_BACKGROUND,
        color = NA
      ),
      legend.key = ggplot2::element_rect(
        fill = FIG_BACKGROUND,
        color = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = panel_fill,
        color = NA
      ),
      plot.background = ggplot2::element_rect(
        fill = FIG_BACKGROUND,
        color = NA
      ),
      panel.border = ggplot2::element_rect(
        color = panel_border_colour,
        fill = NA,
        linewidth = PANEL_BORDER_LWD
      ),
      panel.grid.major = if (show_grid) {
        ggplot2::element_line(
          color = "grey92",
          linewidth = GRID_LWD
        )
      } else {
        ggplot2::element_blank()
      },
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = MANUSCRIPT_MARGIN
    )
}

save_plot_all_formats <- function(plot_obj,
                                  filename_stem,
                                  width = FIGURE_W,
                                  height = FIGURE_H) {
  filename_stem <- sub(
    "\\.(png|tiff|tif|pdf)$",
    "",
    filename_stem,
    ignore.case = TRUE
  )

  png_file <- file.path(OUTPUT_DIR, paste0(filename_stem, ".png"))
  tiff_file <- file.path(OUTPUT_DIR, paste0(filename_stem, ".tiff"))
  pdf_file <- file.path(OUTPUT_DIR, paste0(filename_stem, ".pdf"))

  ggplot2::ggsave(
    filename = png_file,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = FIG_DPI,
    bg = FIG_BACKGROUND,
    limitsize = FALSE
  )

  ggplot2::ggsave(
    filename = tiff_file,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = FIG_DPI,
    compression = "lzw",
    bg = FIG_BACKGROUND,
    limitsize = FALSE
  )

  pdf_device <- if (capabilities("cairo")) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }

  ggplot2::ggsave(
    filename = pdf_file,
    plot = plot_obj,
    width = width,
    height = height,
    device = pdf_device,
    bg = FIG_BACKGROUND,
    limitsize = FALSE
  )

  invisible(
    list(PNG = png_file, TIFF = tiff_file, PDF = pdf_file)
  )
}

find_first_existing <- function(candidates,
                                label) {
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0L) {
    stop(
      label,
      " was not found. Checked:\n",
      paste(candidates, collapse = "\n")
    )
  }

  normalizePath(existing[1], winslash = "/", mustWork = TRUE)
}

cluster_number <- function(x) {
  x <- trimws(as.character(x))

  value <- suppressWarnings(
    as.integer(sub("^.*?([0-9]+)\\s*$", "\\1", x))
  )

  value[!grepl("[0-9]+\\s*$", x)] <- NA_integer_
  value
}

normalise_cluster_label <- function(x) {
  x <- trimws(as.character(x))
  cluster_id <- cluster_number(x)

  x[!is.na(cluster_id)] <- paste0(
    "Cluster ",
    cluster_id[!is.na(cluster_id)]
  )

  x
}

make_group_factor <- function(x) {
  raw <- as.character(x)
  raw_lower <- tolower(trimws(raw))

  if (all(raw %in% c("1", "2"))) {
    return(
      factor(
        raw,
        levels = c("1", "2"),
        labels = c("Normal", "Tumor")
      )
    )
  }

  if (all(
    raw_lower %in% c(
      "normal", "tumor", "non", "lscc", "cancer", "margin"
    )
  )) {
    mapped <- ifelse(
      raw_lower %in% c("tumor", "lscc", "cancer"),
      "Tumor",
      "Normal"
    )

    return(
      factor(
        mapped,
        levels = c("Normal", "Tumor")
      )
    )
  }

  values <- sort(unique(raw))

  if (length(values) != 2L) {
    stop(
      "The bulk group column must contain exactly two classes. Found: ",
      paste(values, collapse = ", ")
    )
  }

  factor(
    ifelse(raw == values[1], "Normal", "Tumor"),
    levels = c("Normal", "Tumor")
  )
}

read_bulk_cpm <- function(path,
                          cohort_label) {
  dat <- as.data.frame(
    data.table::fread(path, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required_columns <- c("Sample", "group")
  missing_columns <- setdiff(required_columns, colnames(dat))

  if (length(missing_columns) > 0L) {
    stop(
      cohort_label,
      " input is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  dat$Sample <- trimws(as.character(dat$Sample))
  dat$Group <- make_group_factor(dat$group)

  if ("batch" %in% colnames(dat)) {
    dat$Batch <- as.character(dat$batch)
  } else {
    dat$Batch <- cohort_label
  }

  metadata_columns <- c("Sample", "group", "batch", "Group", "Batch")
  gene_columns <- setdiff(colnames(dat), metadata_columns)

  if (length(gene_columns) < 100L) {
    stop(
      cohort_label,
      " has too few gene columns after metadata removal."
    )
  }

  expression_samples_by_genes <- as.data.frame(
    lapply(
      dat[, gene_columns, drop = FALSE],
      function(x) suppressWarnings(as.numeric(as.character(x)))
    ),
    check.names = FALSE
  )

  colnames(expression_samples_by_genes) <- toupper(
    colnames(expression_samples_by_genes)
  )

  expression_samples_by_genes <- expression_samples_by_genes[
    ,
    !duplicated(colnames(expression_samples_by_genes)),
    drop = FALSE
  ]

  expression_matrix <- as.matrix(expression_samples_by_genes)
  storage.mode(expression_matrix) <- "numeric"
  expression_matrix[!is.finite(expression_matrix)] <- 0

  if (any(expression_matrix < 0, na.rm = TRUE)) {
    stop(
      cohort_label,
      " contains negative values; expected a non-negative TMM-CPM matrix."
    )
  }

  expression_gene_by_sample <- t(
    log2(expression_matrix + 1)
  )

  rownames(expression_gene_by_sample) <- colnames(
    expression_samples_by_genes
  )

  colnames(expression_gene_by_sample) <- dat$Sample

  list(
    metadata = data.frame(
      Sample = dat$Sample,
      Group = dat$Group,
      Batch = dat$Batch,
      stringsAsFactors = FALSE
    ),
    expression = expression_gene_by_sample,
    source = path,
    cohort = cohort_label
  )
}

run_ssgsea_safe <- function(expression_gene_by_sample,
                            gene_sets,
                            min_size = MIN_SIGNATURE_GENES) {
  expression_gene_by_sample <- as.matrix(expression_gene_by_sample)
  storage.mode(expression_gene_by_sample) <- "numeric"
  expression_gene_by_sample[!is.finite(expression_gene_by_sample)] <- 0

  result <- tryCatch(
    {
      parameter_object <- GSVA::ssgseaParam(
        exprData = expression_gene_by_sample,
        geneSets = gene_sets,
        minSize = min_size,
        maxSize = 500,
        normalize = TRUE
      )

      GSVA::gsva(
        parameter_object,
        verbose = FALSE
      )
    },
    error = function(e_new) {
      message(
        "New GSVA ssGSEA API failed; trying legacy API. Reason: ",
        e_new$message
      )

      tryCatch(
        {
          GSVA::gsva(
            expr = expression_gene_by_sample,
            gset.idx.list = gene_sets,
            method = "ssgsea",
            kcdf = "Gaussian",
            min.sz = min_size,
            max.sz = 500,
            ssgsea.norm = TRUE,
            verbose = FALSE
          )
        },
        error = function(e_old) {
          stop(
            "ssGSEA failed with both GSVA APIs.\n",
            "New API: ", e_new$message,
            "\nLegacy API: ", e_old$message
          )
        }
      )
    }
  )

  as.matrix(result)
}

safe_wilcox <- function(x,
                        group) {
  keep <- is.finite(x) & !is.na(group)
  x <- x[keep]
  group <- droplevels(factor(group[keep]))

  if (length(levels(group)) < 2L || min(table(group)) < 2L) {
    return(NA_real_)
  }

  tryCatch(
    stats::wilcox.test(
      x ~ group,
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}

safe_spearman <- function(x,
                          y) {
  keep <- is.finite(x) & is.finite(y)

  if (sum(keep) < 5L) {
    return(
      data.frame(
        rho = NA_real_,
        p_value = NA_real_,
        n = sum(keep)
      )
    )
  }

  output <- suppressWarnings(
    stats::cor.test(
      x[keep],
      y[keep],
      method = "spearman",
      exact = FALSE
    )
  )

  data.frame(
    rho = as.numeric(output$estimate),
    p_value = as.numeric(output$p.value),
    n = sum(keep)
  )
}

format_fdr <- function(x) {
  if (!is.finite(x)) return("FDR = NA")
  if (x < 0.001) return("FDR < 0.001")
  paste0("FDR = ", sprintf("%.3f", x))
}

safe_annotation_y <- function(y) {
  y <- y[is.finite(y)]

  if (length(y) == 0L) {
    return(1)
  }

  span <- diff(range(y))

  if (!is.finite(span) || span <= 0) {
    span <- max(abs(y), 1) * 0.15
  }

  max(y) + 0.10 * span
}

make_cohort_row_header <- function(label,
                                   text_colour = "black") {
  patchwork::wrap_elements(
    full = grid::textGrob(
      label = label,
      x = 0.5,
      y = 0.5,
      just = "centre",
      gp = grid::gpar(
        fontfamily = FONT_FAMILY,
        fontface = "bold",
        fontsize = 11,
        col = text_colour
      )
    )
  )
}

make_signature_gmt <- function(gene_sets,
                               file_path) {
  gmt_lines <- vapply(
    names(gene_sets),
    function(signature_name) {
      paste(
        c(
          signature_name,
          "Strict_HighCNV_marker_signature_without_MYBL2",
          gene_sets[[signature_name]]
        ),
        collapse = "\t"
      )
    },
    character(1)
  )

  writeLines(gmt_lines, file_path)
  invisible(file_path)
}

# Step 5: signature construction from strict scrna-seq markers

build_two_cluster_signatures <- function(marker_csv,
                                         expression_gene_names) {
  marker_table <- read.csv(
    marker_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_columns <- c(
    "Gene",
    "cluster",
    "p_val_adj",
    "pct.1"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(marker_table)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "The High-CNV marker CSV is incompatible. Missing:\n",
      paste(missing_columns, collapse = ", ")
    )
  }

  fold_change_column <- if ("avg_log2FC" %in% colnames(marker_table)) {
    "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(marker_table)) {
    "avg_logFC"
  } else {
    stop(
      "The High-CNV marker CSV has neither avg_log2FC nor avg_logFC."
    )
  }

  markers <- data.frame(
    Gene = toupper(trimws(as.character(marker_table$Gene))),
    Subcluster = normalise_cluster_label(marker_table$cluster),
    p_val_adj = suppressWarnings(as.numeric(marker_table$p_val_adj)),
    avg_logFC = suppressWarnings(
      as.numeric(marker_table[[fold_change_column]])
    ),
    pct.1 = suppressWarnings(as.numeric(marker_table$pct.1)),
    stringsAsFactors = FALSE
  )

  technical_pattern <- "^MT-|^RPL|^RPS"

  markers <- markers %>%
    dplyr::filter(
      Subcluster %in% SELECTED_SUBCLUSTERS,
      !is.na(Gene),
      nzchar(Gene),
      Gene != TARGET_GENE,
      is.finite(p_val_adj),
      p_val_adj < MARKER_PADJ_CUTOFF,
      is.finite(avg_logFC),
      avg_logFC >= MARKER_LOGFC_CUTOFF,
      is.finite(pct.1),
      pct.1 >= MARKER_PCT_CUTOFF
    )

  if (REMOVE_MITOCHONDRIAL_GENES && REMOVE_RIBOSOMAL_GENES) {
    markers <- markers %>%
      dplyr::filter(!grepl(technical_pattern, Gene))
  } else if (REMOVE_MITOCHONDRIAL_GENES) {
    markers <- markers %>%
      dplyr::filter(!grepl("^MT-", Gene))
  } else if (REMOVE_RIBOSOMAL_GENES) {
    markers <- markers %>%
      dplyr::filter(!grepl("^RPL|^RPS", Gene))
  }

  markers <- markers %>%
    dplyr::filter(Gene %in% expression_gene_names) %>%
    dplyr::group_by(Subcluster, Gene) %>%
    dplyr::slice_min(
      order_by = p_val_adj,
      n = 1L,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup()

  # Keep only genes that are unique to one of the two selected High-CNV
  # subclusters. This prevents shared malignant/cell-cycle markers from
  # contributing to both signatures and makes each signature cluster-specific.
  gene_cluster_counts <- markers %>%
    dplyr::distinct(Subcluster, Gene) %>%
    dplyr::count(Gene, name = "n_selected_clusters")

  markers <- markers %>%
    dplyr::left_join(gene_cluster_counts, by = "Gene") %>%
    dplyr::filter(n_selected_clusters == 1L) %>%
    dplyr::select(-n_selected_clusters) %>%
    dplyr::arrange(
      Subcluster,
      p_val_adj,
      dplyr::desc(avg_logFC),
      dplyr::desc(pct.1),
      Gene
    )

  signature_table <- markers %>%
    dplyr::group_by(Subcluster) %>%
    dplyr::slice_head(n = TOP_SIGNATURE_GENES) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Subcluster) %>%
    dplyr::mutate(
      Signature_rank = dplyr::row_number()
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Signature = paste0(Subcluster, " signature")
    ) %>%
    dplyr::select(
      Signature,
      Subcluster,
      Signature_rank,
      Gene,
      p_val_adj,
      avg_logFC,
      pct.1
    )

  counts <- signature_table %>%
    dplyr::group_by(Signature, Subcluster) %>%
    dplyr::summarise(
      Genes_retained = dplyr::n_distinct(Gene),
      .groups = "drop"
    )

  missing_subclusters <- setdiff(
    SELECTED_SUBCLUSTERS,
    counts$Subcluster
  )

  undersized_subclusters <- counts$Subcluster[
    counts$Genes_retained < MIN_SIGNATURE_GENES
  ]

  if (
    length(missing_subclusters) > 0L ||
      length(undersized_subclusters) > 0L
  ) {
    stop(
      "Too few eligible signature genes after MYBL2/technical-gene removal.\n",
      "Missing signatures: ",
      ifelse(
        length(missing_subclusters) > 0L,
        paste(missing_subclusters, collapse = ", "),
        "none"
      ),
      "\nSignatures below minimum size: ",
      ifelse(
        length(undersized_subclusters) > 0L,
        paste(undersized_subclusters, collapse = ", "),
        "none"
      )
    )
  }

  gene_sets <- split(
    signature_table$Gene,
    signature_table$Signature
  )

  gene_sets <- lapply(gene_sets, unique)

  ordered_signatures <- paste0(
    SELECTED_SUBCLUSTERS,
    " signature"
  )

  gene_sets <- gene_sets[ordered_signatures]

  overlap_matrix <- sapply(
    gene_sets,
    function(x) {
      sapply(
        gene_sets,
        function(y) length(intersect(x, y))
      )
    }
  )

  overlap_table <- as.data.frame(
    overlap_matrix,
    check.names = FALSE
  )

  overlap_table <- tibble::rownames_to_column(
    overlap_table,
    var = "Signature"
  )

  list(
    signature_table = signature_table,
    gene_sets = gene_sets,
    overlap_table = overlap_table
  )
}

# Step 6: statistical and figure helpers

make_group_comparison_table <- function(scores_long,
                                        cohort_label) {
  scores_long %>%
    dplyr::group_by(Signature) %>%
    dplyr::group_modify(
      ~ {
        current <- .x

        normal_values <- current$Score[current$Group == "Normal"]
        tumor_values <- current$Score[current$Group == "Tumor"]

        p_value <- safe_wilcox(
          current$Score,
          current$Group
        )

        data.frame(
          Cohort = cohort_label,
          n_Normal = sum(current$Group == "Normal"),
          n_Tumor = sum(current$Group == "Tumor"),
          Median_Normal = stats::median(normal_values, na.rm = TRUE),
          Median_Tumor = stats::median(tumor_values, na.rm = TRUE),
          Median_difference_Tumor_minus_Normal =
            stats::median(tumor_values, na.rm = TRUE) -
            stats::median(normal_values, na.rm = TRUE),
          Wilcoxon_P_value = p_value,
          stringsAsFactors = FALSE
        )
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      BH_FDR = stats::p.adjust(
        Wilcoxon_P_value,
        method = "BH"
      )
    )
}

make_correlation_table <- function(scores_long) {
  tumor_scores <- scores_long %>%
    dplyr::filter(Group == "Tumor")

  tumor_scores %>%
    dplyr::group_by(Signature) %>%
    dplyr::group_modify(
      ~ {
        current <- .x

        correlation <- safe_spearman(
          current$Score,
          current$MYBL2_expression
        )

        data.frame(
          n_Tumor = correlation$n,
          Spearman_rho = correlation$rho,
          P_value = correlation$p_value,
          stringsAsFactors = FALSE
        )
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      BH_FDR = stats::p.adjust(
        P_value,
        method = "BH"
      )
    )
}

make_single_boxplot <- function(scores_long,
                                group_test_table,
                                signature_name,
                                panel_tag = NULL,
                                y_axis_title = "ssGSEA score",
                                panel_fill = FIG_BACKGROUND,
                                panel_border_colour = "black") {
  data_one <- scores_long %>%
    dplyr::filter(Signature == signature_name)

  test_one <- group_test_table %>%
    dplyr::filter(Signature == signature_name)

  if (nrow(data_one) == 0L || nrow(test_one) != 1L) {
    stop(
      "Cannot create box plot for ",
      signature_name,
      "."
    )
  }

  y_top <- safe_annotation_y(data_one$Score)

  annotation <- paste0(
    "\u0394median = ",
    sprintf(
      "%.2f",
      test_one$Median_difference_Tumor_minus_Normal
    ),
    "\n",
    format_fdr(test_one$BH_FDR)
  )

  ggplot2::ggplot(
    data_one,
    ggplot2::aes(
      x = Group,
      y = Score,
      colour = Group
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.48,
      fill = "white",
      outlier.shape = NA,
      linewidth = 0.85
    ) +
    ggplot2::geom_jitter(
      width = 0.10,
      height = 0,
      size = 1.55,
      alpha = 0.68
    ) +
    ggplot2::annotate(
      "text",
      x = 1.5,
      y = y_top,
      label = annotation,
      family = FONT_FAMILY,
      size = 2.9,
      lineheight = 0.95
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Normal" = COL_NORMAL,
        "Tumor" = COL_TUMOR
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(0.04, 0.18)
      )
    ) +
    ggplot2::labs(
      title = signature_name,
      tag = panel_tag,
      x = NULL,
      y = y_axis_title,
      colour = NULL
    ) +
    theme_manuscript(
      show_grid = FALSE,
      legend_position = "none",
      panel_fill = panel_fill,
      panel_border_colour = panel_border_colour
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PLOT_TITLE_PT,
        hjust = 0.5
      )
    )
}

make_single_correlation_plot <- function(scores_long,
                                         correlation_table,
                                         signature_name,
                                         panel_tag = NULL,
                                         point_colour = COL_CLUSTER0,
                                         y_axis_title = "ssGSEA score") {
  data_one <- scores_long %>%
    dplyr::filter(
      Group == "Tumor",
      Signature == signature_name
    )

  test_one <- correlation_table %>%
    dplyr::filter(Signature == signature_name)

  if (nrow(data_one) == 0L || nrow(test_one) != 1L) {
    stop(
      "Cannot create correlation plot for ",
      signature_name,
      "."
    )
  }

  x_range <- range(data_one$MYBL2_expression, na.rm = TRUE)
  y_range <- range(data_one$Score, na.rm = TRUE)

  x_span <- diff(x_range)
  y_span <- diff(y_range)

  if (!is.finite(x_span) || x_span <= 0) x_span <- 1
  if (!is.finite(y_span) || y_span <= 0) y_span <- 1

  annotation <- paste0(
    "\u03c1 = ",
    sprintf("%.2f", test_one$Spearman_rho),
    "\n",
    format_fdr(test_one$BH_FDR)
  )

  ggplot2::ggplot(
    data_one,
    ggplot2::aes(
      x = MYBL2_expression,
      y = Score
    )
  ) +
    ggplot2::geom_point(
      color = point_colour,
      size = 1.70,
      alpha = 0.78
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 0.78,
      color = COL_FIT,
      fill = COL_BAND,
      alpha = 0.65
    ) +
    ggplot2::annotate(
      "text",
      x = x_range[1] + 0.05 * x_span,
      y = y_range[2] - 0.04 * y_span,
      label = annotation,
      hjust = 0,
      vjust = 1,
      family = FONT_FAMILY,
      size = 2.95,
      lineheight = 0.95
    ) +
    ggplot2::labs(
      title = signature_name,
      tag = panel_tag,
      x = NULL,
      y = y_axis_title
    ) +
    theme_manuscript(
      show_grid = TRUE,
      legend_position = "none"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PLOT_TITLE_PT,
        hjust = 0.5
      )
    )
}

# Step 7: run the complete two-cluster ssgsea workflow

run_two_cluster_ssgsea <- function() {
  writeLines(character(0), RUN_LOG)

  status <- list(
    completed = FALSE,
    stage = "started",
    error = NA_character_
  )

  write_status <- function() {
    writeLines(
      c(
        "Bulk ssGSEA validation of High-CNV malignant Clusters 0 and 8",
        paste0("Completed: ", status$completed),
        paste0("Final stage: ", status$stage),
        paste0(
          "Error: ",
          ifelse(is.na(status$error), "none", status$error)
        ),
        paste0("Output directory: ", OUTPUT_DIR)
      ),
      STATUS_TXT
    )

    writeLines(
      capture.output(sessionInfo()),
      SESSION_TXT
    )
  }

  tryCatch(
    {
      status$stage <- "locating_inputs"
      log_message(
        "STAGE: locating bulk expression matrices and High-CNV marker CSV"
      )

      discovery_cpm_file <- find_first_existing(
        DEVELOPMENT_CPM_CANDIDATES,
        "Discovery bulk CPM matrix"
      )

      external_cpm_file <- find_first_existing(
        EXTERNAL_CPM_CANDIDATES,
        "External GSE130605 bulk CPM matrix"
      )

      marker_csv <- find_first_existing(
        MARKER_CSV_CANDIDATES,
        "High-CNV malignant marker CSV"
      )

      status$stage <- "reading_bulk_data"
      log_message(
        "STAGE: reading discovery and external bulk cohorts"
      )

      discovery <- read_bulk_cpm(
        discovery_cpm_file,
        cohort_label = "Discovery"
      )

      external <- read_bulk_cpm(
        external_cpm_file,
        cohort_label = "External_GSE130605"
      )

      if (!(TARGET_GENE %in% rownames(discovery$expression))) {
        stop("MYBL2 is absent from the discovery expression matrix.")
      }

      if (!(TARGET_GENE %in% rownames(external$expression))) {
        stop("MYBL2 is absent from the external expression matrix.")
      }

      status$stage <- "building_cluster0_cluster8_signatures"
      log_message(
        "STAGE: building strict Cluster 0 and Cluster 8 signatures without MYBL2"
      )

      shared_gene_universe <- intersect(
        rownames(discovery$expression),
        rownames(external$expression)
      )

      signature_object <- build_two_cluster_signatures(
        marker_csv = marker_csv,
        expression_gene_names = shared_gene_universe
      )

      signature_table <- signature_object$signature_table
      signature_gene_sets <- signature_object$gene_sets
      signature_overlap <- signature_object$overlap_table

      write.csv(
        signature_table,
        SIGNATURE_MARKERS_CSV,
        row.names = FALSE
      )

      write.csv(
        signature_overlap,
        SIGNATURE_OVERLAP_CSV,
        row.names = FALSE
      )

      make_signature_gmt(
        signature_gene_sets,
        SIGNATURE_GMT
      )

      signature_status <- signature_table %>%
        dplyr::group_by(Signature, Subcluster) %>%
        dplyr::summarise(
          Genes_in_signature = dplyr::n_distinct(Gene),
          Genes_in_discovery =
            sum(Gene %in% rownames(discovery$expression)),
          Genes_in_external =
            sum(Gene %in% rownames(external$expression)),
          MYBL2_present_in_signature =
            any(Gene == TARGET_GENE),
          .groups = "drop"
        )

      if (any(signature_status$MYBL2_present_in_signature)) {
        stop(
          "MYBL2 was unexpectedly retained in a signature; analysis halted."
        )
      }

      if (any(
        signature_status$Genes_in_discovery < MIN_SIGNATURE_GENES |
          signature_status$Genes_in_external < MIN_SIGNATURE_GENES
      )) {
        stop(
          "At least one signature has fewer than ",
          MIN_SIGNATURE_GENES,
          " genes represented in one bulk cohort."
        )
      }

      write.csv(
        signature_status,
        SIGNATURE_STATUS_CSV,
        row.names = FALSE
      )

      writeLines(
        c(
          "Input provenance",
          paste0("Discovery CPM matrix: ", discovery$source),
          paste0("External CPM matrix: ", external$source),
          paste0("High-CNV marker CSV: ", marker_csv),
          paste0("Output directory: ", OUTPUT_DIR),
          "Analysis scope: existing strict High-CNV Cluster 0 and Cluster 8 markers only; no re-clustering or single-cell preprocessing was rerun.",
          "MYBL2 was removed from both signatures before ssGSEA scoring.",
          "Mitochondrial, ribosomal, and cross-cluster shared marker genes were excluded before signature construction.",
          paste0(
            "Strict marker filter: adjusted P < ",
            MARKER_PADJ_CUTOFF,
            "; avg logFC >= ",
            MARKER_LOGFC_CUTOFF,
            "; pct.1 >= ",
            MARKER_PCT_CUTOFF,
            "."
          ),
          paste0(
            "Top strict markers retained per signature: ",
            TOP_SIGNATURE_GENES,
            "."
          ),
          "ssGSEA input: log2(TMM-CPM + 1).",
          "No batch-effect correction was applied to score comparisons; each cohort is tested separately.",
          "Color palette: teal (Normal), burnt orange (Tumor), plum (Cluster 8 correlation).",
          "External-validation boxplots (panel b) use a subtle blue-tinted panel background and blue border; Normal/Tumor colours are unchanged.",
          "Full-width row headers identify the discovery cohort and independent validation cohort, including GEO accession numbers; the figure uses one explicit patchwork layout.",
          "The x-axis title was intentionally removed from the two MYBL2 correlation panels."
        ),
        PROVENANCE_TXT
      )

      status$stage <- "running_ssgsea"
      log_message(
        "STAGE: running ssGSEA for Cluster 0 and Cluster 8 signatures"
      )

      discovery_gene_sets <- lapply(
        signature_gene_sets,
        function(genes) intersect(
          genes,
          rownames(discovery$expression)
        )
      )

      external_gene_sets <- lapply(
        signature_gene_sets,
        function(genes) intersect(
          genes,
          rownames(external$expression)
        )
      )

      discovery_score_matrix <- run_ssgsea_safe(
        discovery$expression,
        discovery_gene_sets,
        min_size = MIN_SIGNATURE_GENES
      )

      external_score_matrix <- run_ssgsea_safe(
        external$expression,
        external_gene_sets,
        min_size = MIN_SIGNATURE_GENES
      )

      discovery_scores <- as.data.frame(
        t(discovery_score_matrix),
        check.names = FALSE
      )

      external_scores <- as.data.frame(
        t(external_score_matrix),
        check.names = FALSE
      )

      discovery_scores$Sample <- rownames(discovery_scores)
      external_scores$Sample <- rownames(external_scores)

      discovery_scores <- dplyr::left_join(
        discovery$metadata,
        discovery_scores,
        by = "Sample"
      )

      external_scores <- dplyr::left_join(
        external$metadata,
        external_scores,
        by = "Sample"
      )

      discovery_scores$MYBL2_expression <- as.numeric(
        discovery$expression[
          TARGET_GENE,
          discovery_scores$Sample
        ]
      )

      external_scores$MYBL2_expression <- as.numeric(
        external$expression[
          TARGET_GENE,
          external_scores$Sample
        ]
      )

      signature_columns <- names(signature_gene_sets)

      discovery_scores <- discovery_scores[
        ,
        c(
          "Sample",
          "Group",
          "Batch",
          "MYBL2_expression",
          signature_columns
        ),
        drop = FALSE
      ]

      external_scores <- external_scores[
        ,
        c(
          "Sample",
          "Group",
          "Batch",
          "MYBL2_expression",
          signature_columns
        ),
        drop = FALSE
      ]

      write.csv(
        discovery_scores,
        DISCOVERY_SCORES_CSV,
        row.names = FALSE
      )

      write.csv(
        external_scores,
        EXTERNAL_SCORES_CSV,
        row.names = FALSE
      )

      status$stage <- "testing_groups_and_correlations"
      log_message(
        "STAGE: testing Tumor/Normal differences and tumor-only MYBL2 correlations"
      )

      discovery_long <- discovery_scores %>%
        tidyr::pivot_longer(
          cols = dplyr::all_of(signature_columns),
          names_to = "Signature",
          values_to = "Score"
        ) %>%
        dplyr::mutate(
          Cohort = "Discovery",
          Group = factor(
            Group,
            levels = c("Normal", "Tumor")
          )
        )

      external_long <- external_scores %>%
        tidyr::pivot_longer(
          cols = dplyr::all_of(signature_columns),
          names_to = "Signature",
          values_to = "Score"
        ) %>%
        dplyr::mutate(
          Cohort = "External GSE130605",
          Group = factor(
            Group,
            levels = c("Normal", "Tumor")
          )
        )

      discovery_group_tests <- make_group_comparison_table(
        discovery_long,
        cohort_label = "Discovery"
      )

      external_group_tests <- make_group_comparison_table(
        external_long,
        cohort_label = "External_GSE130605"
      )

      discovery_correlation_tests <- make_correlation_table(
        discovery_long
      )

      write.csv(
        discovery_group_tests,
        DISCOVERY_GROUP_TESTS_CSV,
        row.names = FALSE
      )

      write.csv(
        external_group_tests,
        EXTERNAL_GROUP_TESTS_CSV,
        row.names = FALSE
      )

      write.csv(
        discovery_correlation_tests,
        DISCOVERY_CORRELATION_CSV,
        row.names = FALSE
      )

      status$stage <- "creating_one_combined_figure"
      log_message(
        "STAGE: creating final two-cluster ssGSEA manuscript figure"
      )

      cluster0_signature <- "Cluster 0 signature"
      cluster8_signature <- "Cluster 8 signature"

      # Panel a: pooled discovery cohort.
      p_a_cluster0 <- make_single_boxplot(
        discovery_long,
        discovery_group_tests,
        signature_name = cluster0_signature,
        panel_tag = "a"
      )

      p_a_cluster8 <- make_single_boxplot(
        discovery_long,
        discovery_group_tests,
        signature_name = cluster8_signature,
        panel_tag = NULL
      )

      # Panel b: independent external cohort GSE130605.
      # A subtle blue background and border distinguish external validation
      # while preserving the Normal/Tumor colour mapping used in panel a.
      p_b_cluster0 <- make_single_boxplot(
        external_long,
        external_group_tests,
        signature_name = cluster0_signature,
        panel_tag = "b",
        panel_fill = COL_EXTERNAL_PANEL_BG,
        panel_border_colour = COL_EXTERNAL_PANEL_BORDER
      )

      p_b_cluster8 <- make_single_boxplot(
        external_long,
        external_group_tests,
        signature_name = cluster8_signature,
        panel_tag = NULL,
        panel_fill = COL_EXTERNAL_PANEL_BG,
        panel_border_colour = COL_EXTERNAL_PANEL_BORDER
      )

      # Panel c: discovery-tumor MYBL2 associations.
      p_c_cluster0 <- make_single_correlation_plot(
        discovery_long,
        discovery_correlation_tests,
        signature_name = cluster0_signature,
        panel_tag = "c",
        point_colour = COL_CLUSTER0
      )

      p_c_cluster8 <- make_single_correlation_plot(
        discovery_long,
        discovery_correlation_tests,
        signature_name = cluster8_signature,
        panel_tag = NULL,
        point_colour = COL_CLUSTER8
      )

      # Full-width cohort labels are placed above rows a and b.
      # A single explicit patchwork design avoids nested-layout conflicts.
      # Panel c contains tumor-only discovery samples and is defined in the caption.
      header_a <- make_cohort_row_header(
        DISCOVERY_COHORT_HEADER,
        text_colour = "black"
      )

      header_b <- make_cohort_row_header(
        EXTERNAL_COHORT_HEADER,
        text_colour = COL_EXTERNAL_PANEL_BORDER
      )

      figure_design <- "
        AA
        BC
        DD
        EF
        GH
      "

      final_figure <- (
        header_a +
          p_a_cluster0 +
          p_a_cluster8 +
          header_b +
          p_b_cluster0 +
          p_b_cluster8 +
          p_c_cluster0 +
          p_c_cluster8
      ) +
        patchwork::plot_layout(
          design = figure_design,
          heights = c(0.10, 1, 0.10, 1, 1.08)
        )

      save_plot_all_formats(
        final_figure,
        FINAL_FIGURE_STEM,
        width = FIGURE_W,
        height = FIGURE_H
      )

      saveRDS(
        list(
          signatures = signature_gene_sets,
          signature_table = signature_table,
          signature_overlap = signature_overlap,
          discovery_scores = discovery_scores,
          external_scores = external_scores,
          discovery_group_tests = discovery_group_tests,
          external_group_tests = external_group_tests,
          discovery_mybl2_correlations = discovery_correlation_tests,
          discovery_source = discovery$source,
          external_source = external$source,
          marker_source = marker_csv,
          palette = list(
            Normal = COL_NORMAL,
            Tumor = COL_TUMOR,
            Cluster0_correlation = COL_CLUSTER0,
            Cluster8_correlation = COL_CLUSTER8,
            External_panel_background = COL_EXTERNAL_PANEL_BG,
            External_panel_border = COL_EXTERNAL_PANEL_BORDER
          )
        ),
        RESULT_RDS
      )

      writeLines(
        c(
          "Bulk ssGSEA validation of High-CNV malignant Clusters 0 and 8 completed.",
          "",
          "Signatures:",
          "  Cluster 0 and Cluster 8 strict-marker signatures only.",
          "  MYBL2 was excluded from both signature gene lists before ssGSEA.",
          "",
          "Cohorts:",
          "  Discovery: pooled GSE127165 + GSE142083.",
          "  External validation: GSE130605.",
          "",
          "Statistics:",
          "  Tumor-Normal differences: Wilcoxon rank-sum test with BH correction.",
          "  MYBL2-score associations: discovery-tumor Spearman correlation with BH correction.",
          "",
          "Final figure:",
          paste0(
            "  ",
            file.path(OUTPUT_DIR, FINAL_FIGURE_STEM),
            ".png/.tiff/.pdf"
          ),
          "",
          "Figure cohorts:",
          "  Discovery: GSE127165 + GSE142083.",
          "  Independent validation: GSE130605."
        ),
        SUMMARY_TXT
      )

      status$completed <- TRUE
      status$stage <- "finished"
      write_status()

      cat(
        "\n============================================================\n",
        "TWO-CLUSTER BULK ssGSEA ANALYSIS FINISHED SUCCESSFULLY\n",
        "============================================================\n",
        "All outputs were saved in:\n",
        OUTPUT_DIR,
        "\n\nFinal figure:\n",
        file.path(
          OUTPUT_DIR,
          paste0(FINAL_FIGURE_STEM, ".png")
        ),
        "\n============================================================\n",
        sep = ""
      )
    },
    error = function(e) {
      status$error <- conditionMessage(e)
      status$stage <- "failed"
      write_status()

      writeLines(
        c(
          "Two-cluster bulk ssGSEA analysis failed.",
          "",
          conditionMessage(e)
        ),
        file.path(
          OUTPUT_DIR,
          "Clusters0_8_ssGSEA_error.txt"
        )
      )

      stop(e)
    }
  )
}

# Step 8: execute the complete analysis

run_two_cluster_ssgsea()
