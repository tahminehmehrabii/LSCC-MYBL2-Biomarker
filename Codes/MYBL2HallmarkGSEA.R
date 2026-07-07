
# STEP 01: Setup.
rm(list = ls(all.names = TRUE))
invisible(gc())

set.seed(123)

options(
  stringsAsFactors = FALSE,
  scipen = 100,
  timeout = 7200,
  width = 160
)

# STEP 02: Set paths and parameters.
if (as.character(getRversion()) != "4.4.3") {
  stop(
    "Run this script only in R 4.4.3. Current R version: ",
    as.character(getRversion())
  )
}

BASE_RESULTS_DIR <- "D:/LSCC/Results_LSCC"
BULK_ML_DIR <- file.path(BASE_RESULTS_DIR, "ML")

OUTPUT_DIR <- file.path(
  BASE_RESULTS_DIR,
  "MYBL2_Continuous_GSEA"
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

TARGET_GENE <- "MYBL2"

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

MIN_GENE_SET_SIZE <- 10L
MAX_GENE_SET_SIZE <- 500L
FDR_CUTOFF <- 0.05
TOP_PATHWAYS_PER_DIRECTION <- 8L

FOCUS_HALLMARK_PATHWAYS <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_MYC_TARGETS_V1"
)

MYBL2_ASSOCIATION_CSV <- file.path(
  OUTPUT_DIR,
  "mybl2_assoc.csv"
)

HALLMARK_GMT <- file.path(
  OUTPUT_DIR,
  "hallmark.gmt"
)

HALLMARK_GENESET_CSV <- file.path(
  OUTPUT_DIR,
  "hallmark_genes.csv"
)

GSEA_RESULTS_CSV <- file.path(
  OUTPUT_DIR,
  "gsea_results.csv"
)

GSEA_LEADING_EDGE_CSV <- file.path(
  OUTPUT_DIR,
  "leading_edge.csv"
)

FOCUS_PATHWAYS_CSV <- file.path(
  OUTPUT_DIR,
  "focus_pathways.csv"
)

TUMOR_SAMPLE_METADATA_CSV <- file.path(
  OUTPUT_DIR,
  "tumor_meta.csv"
)

RESULT_RDS <- file.path(
  OUTPUT_DIR,
  "gsea_obj.rds"
)

GSEA_FIGURE_STEM <- (
  "fig_gsea_overview"
)

CURVE_FIGURE_STEM <- (
  "fig_gsea_curves"
)

# STEP 03: Set figure style.
FONT_FAMILY <- "Arial"
FIG_DPI <- 600L
FIG_BACKGROUND <- "white"

GSEA_FIGURE_W <- 10.60
GSEA_FIGURE_H <- 8.20

CURVE_FIGURE_W <- 12.20
CURVE_FIGURE_H <- 9.10

BASE_TEXT_PT <- 9.5
AXIS_TEXT_PT <- 8.5
AXIS_TITLE_PT <- 10
LEGEND_TEXT_PT <- 8.5
LEGEND_TITLE_PT <- 9
PLOT_TITLE_PT <- 10
PANEL_TAG_PT <- 14

PANEL_BORDER_LWD <- 0.45
AXIS_LWD <- 0.40
GRID_LWD <- 0.30

COL_NEGATIVE <- "#355C7D"
COL_POSITIVE <- "#C95B74"

CURVE_PATHWAY_COLORS <- c(
  "HALLMARK_E2F_TARGETS" = "#1F5AA6",
  "HALLMARK_G2M_CHECKPOINT" = "#C62828",
  "HALLMARK_DNA_REPAIR" = "#167C5A",
  "HALLMARK_MYC_TARGETS_V1" = "#B26A00"
)

COL_CURVE_ZERO <- "#A8B0B9"
COL_CURVE_TEXT <- "#222831"
COL_CURVE_GRID <- "#EEF1F4"
COL_CURVE_BORDER <- "#D8DEE5"
COL_ANNOTATION_FILL <- "#FFFFFF"
COL_ANNOTATION_BORDER <- "#AEB7C2"

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

# STEP 04: Load packages.
required_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "ggplot2",
  "patchwork",
  "limma",
  "msigdbr",
  "fgsea",
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
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(limma)
})

# STEP 05: Define helper functions.
theme_manuscript <- function(show_grid = FALSE,
                             legend_position = "right") {
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
        size = PLOT_TITLE_PT + 0.5,
        hjust = 0.5,
        color = "black"
      ),
      plot.subtitle = ggplot2::element_text(
        family = FONT_FAMILY,
        size = 8.2,
        hjust = 0.5,
        color = "black"
      ),
      plot.tag = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PANEL_TAG_PT
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
        fill = FIG_BACKGROUND,
        color = NA
      ),
      plot.background = ggplot2::element_rect(
        fill = FIG_BACKGROUND,
        color = NA
      ),
      panel.border = ggplot2::element_rect(
        color = "black",
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
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        color = "black"
      ),
      plot.margin = MANUSCRIPT_MARGIN
    )
}

make_custom_enrichment_plot <- function(pathway_genes,
                                        ranked_stats,
                                        pathway_name,
                                        result_row,
                                        line_color,
                                        panel_tag = NULL) {
  gene_names <- names(ranked_stats)
  ranked_stats <- as.numeric(ranked_stats)
  names(ranked_stats) <- gene_names
  
  if (length(ranked_stats) < 2L ||
      is.null(gene_names) ||
      length(gene_names) != length(ranked_stats) ||
      any(is.na(gene_names)) ||
      any(!nzchar(gene_names))) {
    stop("The ranked statistic must be a named numeric vector.")
  }
  
  pathway_genes <- unique(toupper(as.character(pathway_genes)))
  hit_index <- which(gene_names %in% pathway_genes)
  
  if (length(hit_index) < MIN_GENE_SET_SIZE) {
    stop(
      "Too few ranked genes overlap ",
      pathway_name,
      " to draw the enrichment curve."
    )
  }
  
  hit_weight <- abs(ranked_stats[hit_index])
  hit_weight[!is.finite(hit_weight)] <- 0
  
  if (sum(hit_weight) <= 0) {
    hit_weight <- rep(1, length(hit_index))
  }
  
  increment <- numeric(length(ranked_stats))
  increment[hit_index] <- hit_weight / sum(hit_weight)
  
  miss_count <- length(ranked_stats) - length(hit_index)
  if (miss_count > 0L) {
    increment[-hit_index] <- -1 / miss_count
  }
  
  running_score <- cumsum(increment)
  
  curve_data <- data.frame(
    Rank = seq_along(ranked_stats),
    Running_score = running_score,
    stringsAsFactors = FALSE
  )
  
  y_range <- range(running_score, finite = TRUE)
  y_span <- diff(y_range)
  
  if (!is.finite(y_span) || y_span <= 0) {
    y_span <- max(abs(y_range), 1) * 0.25
  }
  
  rug_bottom <- y_range[1] - 0.105 * y_span
  rug_top <- y_range[1] - 0.030 * y_span
  
  peak_rank <- if (as.numeric(result_row$NES[1]) >= 0) {
    which.max(running_score)
  } else {
    which.min(running_score)
  }
  
  annotation <- paste0(
    "NES = ",
    sprintf("%.2f", as.numeric(result_row$NES[1])),
    "\n",
    format_fdr(as.numeric(result_row$padj[1]))
  )
  
  annotation_x <- max(curve_data$Rank) * 0.055
  annotation_y <- y_range[2] - 0.055 * y_span
  
  ggplot2::ggplot(curve_data, ggplot2::aes(x = Rank, y = Running_score)) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = COL_CURVE_ZERO,
      linewidth = 0.34,
      linetype = "dashed"
    ) +
    ggplot2::geom_segment(
      data = data.frame(Rank = hit_index),
      ggplot2::aes(
        x = Rank,
        xend = Rank,
        y = rug_bottom,
        yend = rug_top
      ),
      inherit.aes = FALSE,
      color = line_color,
      linewidth = 0.24,
      alpha = 0.64,
      lineend = "butt"
    ) +
    ggplot2::geom_line(
      color = line_color,
      linewidth = 1.22,
      lineend = "round",
      linejoin = "round"
    ) +
    ggplot2::geom_point(
      data = curve_data[peak_rank, , drop = FALSE],
      color = line_color,
      fill = "white",
      shape = 21,
      stroke = 0.82,
      size = 2.65
    ) +
    ggplot2::annotate(
      "label",
      x = annotation_x,
      y = annotation_y,
      label = annotation,
      hjust = 0,
      vjust = 1,
      family = FONT_FAMILY,
      size = 3.0,
      lineheight = 0.98,
      color = COL_CURVE_TEXT,
      fill = scales::alpha(COL_ANNOTATION_FILL, 0.86),
      label.size = 0.28,
      label.r = grid::unit(0.12, "lines"),
      label.padding = grid::unit(0.28, "lines")
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.025, 0.055))
    ) +
    ggplot2::labs(
      title = pretty_hallmark_name(pathway_name),
      tag = panel_tag,
      x = "Gene Rank",
      y = "Running enrichment score"
    ) +
    theme_manuscript(
      show_grid = FALSE,
      legend_position = "none"
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(
        fill = "#FFFFFF",
        color = NA
      ),
      panel.border = ggplot2::element_rect(
        color = COL_CURVE_BORDER,
        fill = NA,
        linewidth = 0.34
      ),
      axis.line = ggplot2::element_line(
        color = "#7F8994",
        linewidth = 0.32
      ),
      axis.ticks = ggplot2::element_line(
        color = "#7F8994",
        linewidth = 0.30
      ),
      axis.text = ggplot2::element_text(color = "#3C4650"),
      plot.title = ggplot2::element_text(
        margin = ggplot2::margin(b = 4),
        color = COL_CURVE_TEXT
      ),
      plot.tag = ggplot2::element_text(
        size = PANEL_TAG_PT,
        face = "bold",
        color = COL_CURVE_TEXT
      )
    )
}

save_plot_all_formats <- function(plot_object,
                                  filename_stem,
                                  width,
                                  height) {
  filename_stem <- sub(
    "\\.(png|tiff|tif|pdf)$",
    "",
    filename_stem,
    ignore.case = TRUE
  )
  
  png_file <- file.path(
    OUTPUT_DIR,
    paste0(filename_stem, ".png")
  )
  
  tiff_file <- file.path(
    OUTPUT_DIR,
    paste0(filename_stem, ".tiff")
  )
  
  pdf_file <- file.path(
    OUTPUT_DIR,
    paste0(filename_stem, ".pdf")
  )
  
  ggplot2::ggsave(
    filename = png_file,
    plot = plot_object,
    width = width,
    height = height,
    dpi = FIG_DPI,
    bg = FIG_BACKGROUND,
    limitsize = FALSE
  )
  
  ggplot2::ggsave(
    filename = tiff_file,
    plot = plot_object,
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
    plot = plot_object,
    width = width,
    height = height,
    device = pdf_device,
    bg = FIG_BACKGROUND,
    limitsize = FALSE
  )
  
  invisible(
    list(
      PNG = png_file,
      TIFF = tiff_file,
      PDF = pdf_file
    )
  )
}

find_first_existing <- function(candidates,
                                label) {
  candidates <- unique(
    candidates[
      !is.na(candidates) & nzchar(candidates)
    ]
  )
  
  existing <- candidates[file.exists(candidates)]
  
  if (length(existing) == 0L) {
    stop(
      label,
      " was not found. Checked:\n",
      paste(candidates, collapse = "\n")
    )
  }
  
  normalizePath(
    existing[1],
    winslash = "/",
    mustWork = TRUE
  )
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
  
  known_labels <- c(
    "normal",
    "tumor",
    "non",
    "lscc",
    "cancer",
    "margin"
  )
  
  if (all(raw_lower %in% known_labels)) {
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
  
  observed_values <- sort(unique(raw))
  
  if (length(observed_values) != 2L) {
    stop(
      "The bulk group column must have exactly two classes. Found: ",
      paste(observed_values, collapse = ", ")
    )
  }
  
  factor(
    ifelse(raw == observed_values[1], "Normal", "Tumor"),
    levels = c("Normal", "Tumor")
  )
}

read_discovery_cpm <- function(path) {
  data <- as.data.frame(
    data.table::fread(
      path,
      check.names = FALSE
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  required_columns <- c("Sample", "group")
  
  missing_columns <- setdiff(
    required_columns,
    colnames(data)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Discovery CPM input is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  data$Sample <- trimws(as.character(data$Sample))
  data$Group <- make_group_factor(data$group)
  
  if ("batch" %in% colnames(data)) {
    data$Batch <- as.character(data$batch)
  } else {
    data$Batch <- "1"
  }
  
  metadata_columns <- c(
    "Sample",
    "group",
    "batch",
    "Group",
    "Batch"
  )
  
  gene_columns <- setdiff(
    colnames(data),
    metadata_columns
  )
  
  if (length(gene_columns) < 100L) {
    stop("Discovery CPM input has too few gene columns.")
  }
  
  expression_samples_by_genes <- as.data.frame(
    lapply(
      data[, gene_columns, drop = FALSE],
      function(x) {
        suppressWarnings(as.numeric(as.character(x)))
      }
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
      "Discovery CPM matrix contains negative values; expected non-negative CPM."
    )
  }
  
  expression_gene_by_sample <- t(
    log2(expression_matrix + 1)
  )
  
  rownames(expression_gene_by_sample) <- colnames(
    expression_samples_by_genes
  )
  
  colnames(expression_gene_by_sample) <- data$Sample
  
  list(
    metadata = data.frame(
      Sample = data$Sample,
      Group = data$Group,
      Batch = data$Batch,
      stringsAsFactors = FALSE
    ),
    expression = expression_gene_by_sample,
    source = path
  )
}

get_hallmark_gene_sets <- function() {
  hallmark_table <- tryCatch(
    msigdbr::msigdbr(
      species = "Homo sapiens",
      collection = "H"
    ),
    error = function(new_api_error) {
      tryCatch(
        msigdbr::msigdbr(
          species = "Homo sapiens",
          category = "H"
        ),
        error = function(old_api_error) {
          stop(
            "Could not retrieve Hallmark gene sets through msigdbr.\n",
            "New API: ",
            new_api_error$message,
            "\nLegacy API: ",
            old_api_error$message
          )
        }
      )
    }
  )
  
  required_columns <- c("gs_name", "gene_symbol")
  
  missing_columns <- setdiff(
    required_columns,
    colnames(hallmark_table)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "msigdbr Hallmark table is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  hallmark_table <- hallmark_table %>%
    dplyr::transmute(
      Pathway = as.character(gs_name),
      Gene = toupper(trimws(as.character(gene_symbol)))
    ) %>%
    dplyr::filter(
      !is.na(Gene),
      nzchar(Gene)
    ) %>%
    dplyr::distinct()
  
  pathways <- split(
    hallmark_table$Gene,
    hallmark_table$Pathway
  )
  
  pathways <- lapply(pathways, unique)
  
  list(
    pathways = pathways,
    table = hallmark_table
  )
}

write_gmt <- function(gene_sets,
                      file_path) {
  gmt_lines <- vapply(
    names(gene_sets),
    function(pathway_name) {
      paste(
        c(
          pathway_name,
          "MSigDB_Hallmark_Homo_sapiens",
          gene_sets[[pathway_name]]
        ),
        collapse = "\t"
      )
    },
    character(1)
  )
  
  writeLines(gmt_lines, file_path)
  
  invisible(file_path)
}

run_fgsea_safe <- function(pathways,
                           stats_vector) {
  if (!is.numeric(stats_vector) ||
      !is.atomic(stats_vector) ||
      is.list(stats_vector) ||
      length(stats_vector) == 0L ||
      length(names(stats_vector)) != length(stats_vector)) {
    stop(
      "The GSEA ranked statistic must be a named atomic numeric vector."
    )
  }
  
  output <- tryCatch(
    fgsea::fgseaMultilevel(
      pathways = pathways,
      stats = stats_vector,
      minSize = MIN_GENE_SET_SIZE,
      maxSize = MAX_GENE_SET_SIZE,
      eps = 0
    ),
    error = function(multilevel_error) {
      
      tryCatch(
        fgsea::fgseaSimple(
          pathways = pathways,
          stats = stats_vector,
          minSize = MIN_GENE_SET_SIZE,
          maxSize = MAX_GENE_SET_SIZE,
          nperm = 10000
        ),
        error = function(simple_error) {
          stop(
            "GSEA failed with both fgsea methods.\n",
            "fgseaMultilevel: ",
            multilevel_error$message,
            "\nfgseaSimple: ",
            simple_error$message
          )
        }
      )
    }
  )
  
  as.data.frame(output)
}

format_fdr <- function(x) {
  if (!is.finite(x)) {
    return("FDR = NA")
  }
  
  if (x < 0.001) {
    return("FDR < 0.001")
  }
  
  paste0("FDR = ", sprintf("%.3f", x))
}

pretty_hallmark_name <- function(x) {
  x <- as.character(x)
  
  labels <- c(
    "HALLMARK_E2F_TARGETS" = "E2F Targets",
    "HALLMARK_G2M_CHECKPOINT" = "G2/M Checkpoint",
    "HALLMARK_DNA_REPAIR" = "DNA Repair",
    "HALLMARK_MYC_TARGETS_V1" = "MYC Targets V1",
    "HALLMARK_MYC_TARGETS_V2" = "MYC Targets V2",
    "HALLMARK_MTORC1_SIGNALING" = "mTORC1 Signaling",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING" = "IL6/JAK/STAT3 Signaling",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = "TNFα Signaling via NF-κB",
    "HALLMARK_KRAS_SIGNALING_UP" = "KRAS Signaling Up"
  )
  
  output <- unname(labels[x])
  missing <- is.na(output)
  
  if (any(missing)) {
    fallback <- x[missing]
    fallback <- gsub("^HALLMARK_", "", fallback)
    fallback <- gsub("_", " ", fallback)
    output[missing] <- tools::toTitleCase(tolower(fallback))
  }
  
  output
}

wrap_hallmark_name <- function(x,
                               width = 28L) {
  vapply(
    pretty_hallmark_name(x),
    function(label) {
      paste(
        strwrap(label, width = width),
        collapse = "\n"
      )
    },
    character(1)
  )
}

make_fgsea_csv_safe <- function(data) {
  data <- as.data.frame(data)
  
  if (!("leadingEdge" %in% colnames(data))) {
    return(data)
  }
  
  data$leadingEdge_collapsed <- vapply(
    data$leadingEdge,
    function(genes) {
      paste(as.character(genes), collapse = ";")
    },
    character(1)
  )
  
  data$leadingEdge <- NULL
  
  data
}

make_leading_edge_table <- function(gsea_results) {
  if (!("leadingEdge" %in% colnames(gsea_results))) {
    stop(
      "The GSEA result table does not contain the leadingEdge column."
    )
  }
  
  output_rows <- lapply(
    seq_len(nrow(gsea_results)),
    function(index) {
      genes <- as.character(
        gsea_results$leadingEdge[[index]]
      )
      
      if (length(genes) == 0L) {
        return(NULL)
      }
      
      data.frame(
        pathway = rep(
          as.character(gsea_results$pathway[index]),
          length(genes)
        ),
        NES = rep(
          as.numeric(gsea_results$NES[index]),
          length(genes)
        ),
        padj = rep(
          as.numeric(gsea_results$padj[index]),
          length(genes)
        ),
        Direction = rep(
          as.character(gsea_results$Direction[index]),
          length(genes)
        ),
        Leading_edge_gene = genes,
        stringsAsFactors = FALSE
      )
    }
  )
  
  dplyr::bind_rows(output_rows)
}

make_ranked_limma_stats <- function(association_table,
                                    target_gene,
                                    epsilon = 1e-10) {
  ranked_table <- association_table %>%
    dplyr::transmute(
      Gene = toupper(trimws(as.character(Gene))),
      Statistic = as.numeric(t)
    ) %>%
    dplyr::filter(
      is.finite(Statistic),
      !is.na(Gene),
      nzchar(Gene),
      Gene != target_gene
    ) %>%
    dplyr::group_by(Gene) %>%
    dplyr::slice_max(
      order_by = abs(Statistic),
      n = 1L,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(
      dplyr::desc(Statistic),
      Gene
    )
  
  if (nrow(ranked_table) < 1000L) {
    stop(
      "Too few genes remain after removing MYBL2 and duplicated gene symbols."
    )
  }
  
  statistic_values <- as.numeric(ranked_table$Statistic)
  gene_names <- as.character(ranked_table$Gene)
  
  run_lengths <- rle(statistic_values)$lengths
  run_start <- 1L
  
  for (run_length in run_lengths) {
    if (run_length > 1L) {
      positions <- seq.int(
        from = run_start,
        length.out = run_length
      )
      
      offsets <- seq(
        from = run_length,
        to = 1L,
        by = -1L
      ) * epsilon
      
      statistic_values[positions] <- (
        statistic_values[positions] + offsets
      )
    }
    
    run_start <- run_start + run_length
  }
  
  final_order <- order(
    -statistic_values,
    gene_names,
    method = "radix"
  )
  
  ranked_stats <- as.numeric(
    statistic_values[final_order]
  )
  
  names(ranked_stats) <- gene_names[final_order]
  storage.mode(ranked_stats) <- "double"
  
  if (!is.atomic(ranked_stats) ||
      is.list(ranked_stats) ||
      !is.numeric(ranked_stats) ||
      any(!is.finite(ranked_stats))) {
    stop(
      "The final GSEA ranking could not be constructed as an atomic numeric vector."
    )
  }
  
  ranked_stats
}

# STEP 06: Define plotting functions.
make_gsea_overview_plot <- function(gsea_results) {
  usable_results <- gsea_results %>%
    dplyr::filter(
      is.finite(NES),
      is.finite(padj)
    )
  
  selected <- usable_results %>%
    dplyr::filter(padj < FDR_CUTOFF) %>%
    dplyr::mutate(
      Direction_short = ifelse(
        NES >= 0,
        "Positive",
        "Negative"
      )
    ) %>%
    dplyr::group_by(Direction_short) %>%
    dplyr::arrange(
      dplyr::desc(abs(NES)),
      padj,
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(
      n = TOP_PATHWAYS_PER_DIRECTION
    ) %>%
    dplyr::ungroup()
  
  if (nrow(selected) == 0L) {
    selected <- usable_results %>%
      dplyr::arrange(
        dplyr::desc(abs(NES)),
        padj
      ) %>%
      dplyr::slice_head(
        n = 2L * TOP_PATHWAYS_PER_DIRECTION
      ) %>%
      dplyr::mutate(
        Direction_short = ifelse(
          NES >= 0,
          "Positive",
          "Negative"
        )
      )
  }
  
  selected <- selected %>%
    dplyr::mutate(
      Direction_short = factor(
        Direction_short,
        levels = c("Negative", "Positive")
      ),
      Label = wrap_hallmark_name(pathway),
      neglog10FDR = -log10(pmax(padj, 1e-300))
    ) %>%
    dplyr::arrange(NES)
  
  selected$Label <- factor(
    selected$Label,
    levels = selected$Label
  )
  
  ggplot2::ggplot(
    selected,
    ggplot2::aes(
      x = NES,
      y = Label
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      color = "grey40",
      linewidth = AXIS_LWD
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = NES,
        y = Label,
        yend = Label,
        color = Direction_short
      ),
      linewidth = 0.85,
      alpha = 0.72
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = neglog10FDR,
        color = Direction_short
      ),
      alpha = 0.95
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Negative" = COL_NEGATIVE,
        "Positive" = COL_POSITIVE
      ),
      name = "Association direction"
    ) +
    ggplot2::scale_size_continuous(
      range = c(2.8, 7.2),
      name = expression(-log[10]~FDR)
    ) +
    ggplot2::labs(
      tag = "a",
      x = "Normalized enrichment score",
      y = NULL
    ) +
    theme_manuscript(
      show_grid = FALSE,
      legend_position = "right"
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        size = 8.0,
        lineheight = 0.92
      )
    )
}

make_focus_enrichment_plots <- function(pathways,
                                        ranked_stats,
                                        gsea_results) {
  focus_existing <- FOCUS_HALLMARK_PATHWAYS[
    FOCUS_HALLMARK_PATHWAYS %in% names(pathways)
  ]
  
  if (length(focus_existing) != length(FOCUS_HALLMARK_PATHWAYS)) {
    missing_focus <- setdiff(
      FOCUS_HALLMARK_PATHWAYS,
      focus_existing
    )
    
    stop(
      "The following requested Hallmark pathways were unavailable after ",
      "gene-overlap filtering: ",
      paste(missing_focus, collapse = ", ")
    )
  }
  
  plots <- lapply(
    focus_existing,
    function(pathway_name) {
      current_result <- gsea_results[
        gsea_results$pathway == pathway_name,
        ,
        drop = FALSE
      ]
      
      current_color <- unname(CURVE_PATHWAY_COLORS[[pathway_name]])
      if (is.null(current_color) || is.na(current_color)) {
        current_color <- "#1F5AA6"
      }
      
      if (nrow(current_result) != 1L) {
        stop("No unique GSEA result was available for ", pathway_name, ".")
      }
      
      make_custom_enrichment_plot(
        pathway_genes = pathways[[pathway_name]],
        ranked_stats = ranked_stats,
        pathway_name = pathway_name,
        result_row = current_result,
        line_color = current_color
      )
    }
  )
  
  names(plots) <- focus_existing
  
  plots <- Map(
    function(plot_object, panel_tag) {
      plot_object + ggplot2::labs(tag = panel_tag)
    },
    plots,
    letters[seq_along(plots)]
  )
  
  list(
    plots = plots,
    data = gsea_results %>%
      dplyr::filter(
        pathway %in% focus_existing
      ) %>%
      dplyr::arrange(
        match(pathway, focus_existing)
      )
  )
}

# STEP 07: Run analysis.
run_mybl2_limma_batch_adjusted_gsea <- function() {
  tryCatch(
    {
      
      discovery_cpm_file <- find_first_existing(
        DEVELOPMENT_CPM_CANDIDATES,
        "Discovery bulk CPM matrix"
      )
      
      input <- read_discovery_cpm(
        discovery_cpm_file
      )
      
      if (!(TARGET_GENE %in% rownames(input$expression))) {
        stop(
          "MYBL2 is absent from the discovery expression matrix."
        )
      }
      
      tumor_metadata <- input$metadata %>%
        dplyr::filter(Group == "Tumor") %>%
        dplyr::arrange(Sample)
      
      if (nrow(tumor_metadata) < 20L) {
        stop(
          "Fewer than 20 tumor samples are available for MYBL2-continuous GSEA."
        )
      }
      
      tumor_expression <- input$expression[
        ,
        tumor_metadata$Sample,
        drop = FALSE
      ]
      
      tumor_metadata$MYBL2_expression <- as.numeric(
        tumor_expression[
          TARGET_GENE,
          tumor_metadata$Sample
        ]
      )
      
      if (
        sum(is.finite(tumor_metadata$MYBL2_expression)) < 20L ||
        stats::sd(
          tumor_metadata$MYBL2_expression,
          na.rm = TRUE
        ) == 0
      ) {
        stop(
          "MYBL2 expression has insufficient finite variation among tumor samples."
        )
      }
      
      tumor_metadata$MYBL2_z <- as.numeric(
        scale(tumor_metadata$MYBL2_expression)
      )
      
      tumor_metadata$Batch <- factor(
        as.character(tumor_metadata$Batch)
      )
      
      write.csv(
        tumor_metadata,
        TUMOR_SAMPLE_METADATA_CSV,
        row.names = FALSE
      )
      
      gene_variances <- apply(
        tumor_expression,
        1,
        stats::var,
        na.rm = TRUE
      )
      
      keep_variable_genes <- (
        is.finite(gene_variances) &
          gene_variances > 0
      )
      
      if (sum(keep_variable_genes) < 1000L) {
        stop(
          "Too few variable genes remain after zero-variance filtering."
        )
      }
      
      tumor_expression <- tumor_expression[
        keep_variable_genes,
        ,
        drop = FALSE
      ]
      
      if (nlevels(tumor_metadata$Batch) >= 2L) {
        design <- stats::model.matrix(
          ~ MYBL2_z + Batch,
          data = tumor_metadata
        )
        batch_adjustment_used <- TRUE
      } else {
        design <- stats::model.matrix(
          ~ MYBL2_z,
          data = tumor_metadata
        )
        batch_adjustment_used <- FALSE
      }
      
      if (!("MYBL2_z" %in% colnames(design))) {
        stop(
          "The MYBL2 continuous coefficient was not found in the limma design."
        )
      }
      
      if (qr(design)$rank < ncol(design)) {
        stop(
          "The limma design matrix is rank deficient. MYBL2 and batch may be ",
          "perfectly confounded; batch-adjusted analysis cannot be estimated."
        )
      }
      
      fit <- limma::lmFit(
        tumor_expression,
        design
      )
      
      fit <- limma::eBayes(
        fit,
        trend = TRUE
      )
      
      association_table <- limma::topTable(
        fit,
        coef = "MYBL2_z",
        number = Inf,
        adjust.method = "BH",
        sort.by = "t"
      )
      
      association_table$Gene <- rownames(association_table)
      
      association_table <- association_table %>%
        dplyr::transmute(
          Gene = as.character(Gene),
          logFC = as.numeric(logFC),
          t = as.numeric(t),
          P.Value = as.numeric(P.Value),
          adj.P.Val = as.numeric(adj.P.Val),
          B = as.numeric(B)
        ) %>%
        dplyr::arrange(
          dplyr::desc(t),
          P.Value,
          Gene
        )
      
      if (nrow(association_table) < 1000L) {
        stop(
          "Too few genes were retained by limma for GSEA."
        )
      }
      
      write.csv(
        association_table,
        MYBL2_ASSOCIATION_CSV,
        row.names = FALSE
      )
      
      ranked_stats <- make_ranked_limma_stats(
        association_table = association_table,
        target_gene = TARGET_GENE,
        epsilon = 1e-10
      )
      
      if (length(ranked_stats) < 1000L) {
        stop(
          "Too few genes remain in the MYBL2-associated ranked list for GSEA."
        )
      }
      
      hallmark <- get_hallmark_gene_sets()
      
      hallmark_pathways <- lapply(
        hallmark$pathways,
        function(genes) {
          intersect(
            genes,
            names(ranked_stats)
          )
        }
      )
      
      hallmark_pathways <- hallmark_pathways[
        vapply(
          hallmark_pathways,
          length,
          integer(1)
        ) >= MIN_GENE_SET_SIZE
      ]
      
      if (length(hallmark_pathways) < 10L) {
        stop(
          "Too few Hallmark pathways overlap the MYBL2-ranked gene list."
        )
      }
      
      hallmark_table_used <- hallmark$table %>%
        dplyr::filter(
          Pathway %in% names(hallmark_pathways),
          Gene %in% names(ranked_stats)
        ) %>%
        dplyr::arrange(
          Pathway,
          Gene
        )
      
      write.csv(
        hallmark_table_used,
        HALLMARK_GENESET_CSV,
        row.names = FALSE
      )
      
      write_gmt(
        hallmark_pathways,
        HALLMARK_GMT
      )
      
      gsea_results <- run_fgsea_safe(
        pathways = hallmark_pathways,
        stats_vector = ranked_stats
      )
      
      if (nrow(gsea_results) == 0L) {
        stop("Hallmark GSEA returned no results.")
      }
      
      gsea_results <- gsea_results %>%
        dplyr::mutate(
          leadingEdge_collapsed = vapply(
            leadingEdge,
            function(genes) {
              paste(as.character(genes), collapse = ";")
            },
            character(1)
          ),
          Direction = ifelse(
            NES >= 0,
            "Positive association with MYBL2",
            "Negative association with MYBL2"
          )
        ) %>%
        dplyr::arrange(
          padj,
          dplyr::desc(abs(NES))
        )
      
      gsea_results_export <- make_fgsea_csv_safe(
        gsea_results
      ) %>%
        dplyr::select(
          pathway,
          pval,
          padj,
          ES,
          NES,
          size,
          Direction,
          leadingEdge_collapsed
        )
      
      write.csv(
        gsea_results_export,
        GSEA_RESULTS_CSV,
        row.names = FALSE
      )
      
      leading_edge_table <- make_leading_edge_table(
        gsea_results
      )
      
      write.csv(
        leading_edge_table,
        GSEA_LEADING_EDGE_CSV,
        row.names = FALSE
      )
      
      gsea_overview_plot <- make_gsea_overview_plot(
        gsea_results
      )
      
      focus <- make_focus_enrichment_plots(
        pathways = hallmark_pathways,
        ranked_stats = ranked_stats,
        gsea_results = gsea_results
      )
      
      focus_export <- make_fgsea_csv_safe(
        focus$data
      ) %>%
        dplyr::select(
          pathway,
          pval,
          padj,
          ES,
          NES,
          size,
          Direction,
          leadingEdge_collapsed
        )
      
      write.csv(
        focus_export,
        FOCUS_PATHWAYS_CSV,
        row.names = FALSE
      )
      
      if (length(focus$plots) != 4L) {
        stop(
          "Exactly four focus pathways are required to create the 2 x 2 curve figure."
        )
      }
      
      curve_plot <- (
        focus$plots[[1]] |
          focus$plots[[2]]
      ) / (
        focus$plots[[3]] |
          focus$plots[[4]]
      )
      
      save_plot_all_formats(
        gsea_overview_plot,
        GSEA_FIGURE_STEM,
        width = GSEA_FIGURE_W,
        height = GSEA_FIGURE_H
      )
      
      save_plot_all_formats(
        curve_plot,
        CURVE_FIGURE_STEM,
        width = CURVE_FIGURE_W,
        height = CURVE_FIGURE_H
      )
      
      saveRDS(
        list(
          input_source = input$source,
          tumor_metadata = tumor_metadata,
          batch_adjustment_used = batch_adjustment_used,
          design = design,
          association_table = association_table,
          ranked_stats = ranked_stats,
          hallmark_pathways = hallmark_pathways,
          gsea_results = gsea_results,
          focus_pathways = focus$data,
          focus_pathways_export = focus_export
        ),
        RESULT_RDS
      )
    },
    error = function(error_object) {
      stop(error_object)
    }
  )
}

# STEP 08: Execute.
run_mybl2_limma_batch_adjusted_gsea()
