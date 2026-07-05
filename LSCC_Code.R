# LSCC MYBL2 Bioinformatics Pipeline — GitHub Single-Script Version

# MODULE 01 — Bulk RNA-seq / ML pipeline and Figure 2 outputs

# STEP 01. SETUP, PACKAGES, PATHS, AND GLOBAL HELPERS

rm(list = ls())
gc()

set.seed(123)
options(stringsAsFactors = FALSE)
options(scipen = 100)
options(timeout = 3600)

project_path <- "D:/LSCC/Bulk_RNASeq_ML_Results"
results_path <- "D:/LSCC/Results_LSCC/ML"
manuscript_figures_path <- results_path
FIGURE_DIR <- results_path

dir.create(results_path, recursive = TRUE, showWarnings = FALSE)

scrna_figdir <- "D:/LSCC/ScRNAseq_Results/GSE206332/Results/figures"

TARGET_GENE <- "MYBL2"

# 03. STANDARD MANUSCRIPT FIGURE SETTINGS

if (.Platform$OS.type == "windows") {
  try(
    grDevices::windowsFonts(
      Arial = grDevices::windowsFont("Arial")
    ),
    silent = TRUE
  )
}

FONT_FAMILY <- "Arial"
FIG_DPI <- 600
FIG_BACKGROUND <- "white"

FIG_SINGLE_W <- 3.50
FIG_SINGLE_H <- 4.20

FIG_DOUBLE_W <- 7.20
FIG_DOUBLE_H <- 5.40

FIGURE2_COMBINED_W <- 14.00
FIGURE2_COMBINED_H <- 6.40

BASE_TEXT_PT <- 9.5

AXIS_TEXT_PT <- 8.5
AXIS_TITLE_PT <- 10

LEGEND_TEXT_PT <- 8.5
LEGEND_TITLE_PT <- 9

PLOT_TITLE_PT <- 10
PANEL_TAG_PT <- 14

PANEL_TAG_POSITION <- c(0.015, 0.985)
PANEL_TAG_RASTER_SIZE <- 135
PANEL_TAG_RASTER_OFFSET_X <- 80
PANEL_TAG_RASTER_OFFSET_Y <- 65
PANEL_BORDER_RASTER_PX <- 3
PANEL_CONTENT_INSET_PX <- 260
PANEL_CONTENT_VERTICAL_INSET_PX <- 220

TICK_LENGTH_PT <- 2
PANEL_SPACING_PT <- 6
BASE_AXIS_CEX <- 0.90
BASE_AXIS_TITLE_CEX <- 1.00
BASE_ANNOTATION_CEX <- 0.82

PANEL_BORDER_LWD <- 0.45
AXIS_LWD <- 0.40
GEOM_LWD <- 0.55
GRID_LWD <- 0.30

HEATMAP_FONT_PT <- 10
HEATMAP_ROW_FONT_PT <- 8
HEATMAP_COL_FONT_PT <- 8
HEATMAP_LEGEND_FONT_PT <- 10

COL_NORMAL <- "#7470B2"
COL_TUMOR <- "#D81B60"

COL_MYBL2_LOW <- COL_NORMAL
COL_MYBL2_HIGH <- COL_TUMOR

COL_UP <- COL_TUMOR
COL_DOWN <- COL_NORMAL
COL_NS <- "#B8B8B8"

CELLCHAT_PALETTE <- grDevices::colorRampPalette(
  c("#E8F6F4", "#B8E1DA", "#73C6BE", "#2D9D8C", "#006E63")
)(100)

MANUSCRIPT_MARGIN <- ggplot2::margin(
  t = 5.5,
  r = 6,
  b = 5.5,
  l = 6,
  unit = "pt"
)

# 04. STANDARD GGPLOT THEME

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
        size = PLOT_TITLE_PT,
        hjust = 0.5,
        vjust = 0.5,
        color = "black"
      ),
      plot.title.position = "panel",

      plot.tag = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PANEL_TAG_PT,
        hjust = 0,
        vjust = 1,
        color = "black"
      ),
      plot.tag.position = PANEL_TAG_POSITION,

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
      axis.ticks.length = grid::unit(TICK_LENGTH_PT, "pt"),

      legend.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = LEGEND_TITLE_PT,
        color = "black"
      ),
      legend.text = ggplot2::element_text(
        family = FONT_FAMILY,
        size = LEGEND_TEXT_PT,
        color = "black"
      ),
      legend.position = legend_position,
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

      strip.background = ggplot2::element_rect(
        fill = "grey95",
        color = "black",
        linewidth = PANEL_BORDER_LWD
      ),
      strip.text = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = AXIS_TEXT_PT,
        color = "black"
      ),

      panel.spacing = grid::unit(PANEL_SPACING_PT, "pt"),
      plot.margin = MANUSCRIPT_MARGIN
    )
}

# 05. EXPORT FUNCTIONS: PNG + TIFF + PDF

save_plot_all_formats <- function(plot_obj,
                                  filename_stem,
                                  dir_path,
                                  width = FIG_DOUBLE_W,
                                  height = FIG_DOUBLE_H) {

  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  png_file <- file.path(dir_path, paste0(filename_stem, ".png"))
  tiff_file <- file.path(dir_path, paste0(filename_stem, ".tiff"))
  pdf_file <- file.path(dir_path, paste0(filename_stem, ".pdf"))

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
    bg = FIG_BACKGROUND,
    compression = "lzw",
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
    list(
      PNG = png_file,
      TIFF = tiff_file,
      PDF = pdf_file
    )
  )
}

set_grob_font_family <- function(grob_object,
                                 font_family = FONT_FAMILY) {

  if (is.null(grob_object)) {
    return(grob_object)
  }

  if (!is.null(grob_object$gp)) {
    grob_object$gp$fontfamily <- font_family
  }

  if (!is.null(grob_object$children)) {
    for (i in seq_along(grob_object$children)) {
      grob_object$children[[i]] <- set_grob_font_family(
        grob_object$children[[i]],
        font_family
      )
    }
  }

  if (!is.null(grob_object$grobs)) {
    for (i in seq_along(grob_object$grobs)) {
      grob_object$grobs[[i]] <- set_grob_font_family(
        grob_object$grobs[[i]],
        font_family
      )
    }
  }

  grob_object
}

save_grob_all_formats <- function(grob_obj,
                                  filename_stem,
                                  dir_path,
                                  width = FIG_DOUBLE_W,
                                  height = FIG_DOUBLE_H) {

  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  png_file <- file.path(dir_path, paste0(filename_stem, ".png"))
  tiff_file <- file.path(dir_path, paste0(filename_stem, ".tiff"))
  pdf_file <- file.path(dir_path, paste0(filename_stem, ".pdf"))

  grob_obj <- set_grob_font_family(grob_obj, FONT_FAMILY)

  draw_grob <- function() {
    grid::grid.newpage()
    grid::grid.rect(
      gp = grid::gpar(
        fill = FIG_BACKGROUND,
        col = NA
      )
    )
    grid::grid.draw(grob_obj)
  }

  grDevices::png(
    filename = png_file,
    width = width,
    height = height,
    units = "in",
    res = FIG_DPI,
    bg = FIG_BACKGROUND
  )
  draw_grob()
  grDevices::dev.off()

  grDevices::tiff(
    filename = tiff_file,
    width = width,
    height = height,
    units = "in",
    res = FIG_DPI,
    compression = "lzw",
    bg = FIG_BACKGROUND
  )
  draw_grob()
  grDevices::dev.off()

  if (capabilities("cairo")) {
    grDevices::cairo_pdf(
      filename = pdf_file,
      width = width,
      height = height,
      family = FONT_FAMILY
    )
  } else {
    grDevices::pdf(
      file = pdf_file,
      width = width,
      height = height,
      family = FONT_FAMILY
    )
  }

  draw_grob()
  grDevices::dev.off()

  invisible(
    list(
      PNG = png_file,
      TIFF = tiff_file,
      PDF = pdf_file
    )
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(limma)
  library(pheatmap)
  library(patchwork)
  library(grid)
  library(impute)
  library(WGCNA)
  library(magick)
  library(ggvenn)
  library(glmnet)
  library(pROC)
  library(cowplot)
  library(randomForest)
  library(e1071)
  library(msigdbr)
  library(stringr)
  library(tibble)
  library(clusterProfiler)
  library(GSVA)
})

HAS_GBM <- requireNamespace("gbm", quietly = TRUE)

make_stage_dir <- function(step_number, step_name) {
  dir.create(results_path, recursive = TRUE, showWarnings = FALSE)
  results_path
}

save_fig <- function(plot_obj, filename, dir_path, width = 8, height = 6, dpi = FIG_DPI) {
  ggplot2::ggsave(
    filename = file.path(dir_path, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    bg = FIG_BACKGROUND,
    limitsize = FALSE
  )
}

save_manuscript_copy <- function(from_file, figure_name) {
  if (!file.exists(from_file)) return(invisible(NULL))

  to_file <- file.path(manuscript_figures_path, figure_name)

  if (identical(
    normalizePath(from_file, winslash = "/", mustWork = FALSE),
    normalizePath(to_file, winslash = "/", mustWork = FALSE)
  )) {
    return(invisible(to_file))
  }

  file.copy(from = from_file, to = to_file, overwrite = TRUE)
  invisible(to_file)
}

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

guess_sample_col <- function(ph) {
  cand <- c("Sample", "sample", "GSM", "gsm", "geo_accession", "geo", "run", "Run")
  hit <- cand[cand %in% names(ph)]
  if (length(hit) > 0) return(hit[1])
  names(ph)[1]
}

guess_group_col <- function(ph) {
  cand <- c(
    "group", "Group", "label", "Label", "condition", "Condition",
    "status", "Status", "phenotype", "Phenotype", "class", "Class"
  )
  hit <- cand[cand %in% names(ph)]
  if (length(hit) > 0) return(hit[1])
  stop("Could not find group/label column in pheno. Columns are: ",
       paste(names(ph), collapse = ", "))
}

first_existing <- function(paths) {
  paths <- paths[file.exists(paths)]
  if (length(paths) < 1) return(NA_character_)
  paths[1]
}

add_panel_label_img <- function(img,
                                label,
                                size = PANEL_TAG_RASTER_SIZE,
                                offset_x = PANEL_TAG_RASTER_OFFSET_X,
                                offset_y = PANEL_TAG_RASTER_OFFSET_Y) {
  magick::image_annotate(
    img,
    text = label,
    size = size,
    font = FONT_FAMILY,
    gravity = "northwest",
    location = paste0("+", offset_x, "+", offset_y),
    weight = 700,
    color = "black"
  )
}

make_clean_panel_img <- function(path,
                                 label,
                                 panel_width = 3600,
                                 panel_height = 2850,
                                 label_size = PANEL_TAG_RASTER_SIZE) {
  img <- magick::image_read(path)
  img <- magick::image_background(img, FIG_BACKGROUND, flatten = TRUE)
  img <- magick::image_trim(img)

  img <- magick::image_resize(
    img,
    geometry = paste0(
      panel_width - PANEL_CONTENT_INSET_PX, "x",
      panel_height - PANEL_CONTENT_VERTICAL_INSET_PX, ">"
    )
  )

  img <- magick::image_extent(
    img,
    geometry = paste0(panel_width, "x", panel_height),
    gravity = "center",
    color = FIG_BACKGROUND
  )

  img <- magick::image_border(
    img,
    color = "black",
    geometry = paste0(PANEL_BORDER_RASTER_PX, "x", PANEL_BORDER_RASTER_PX)
  )

  add_panel_label_img(
    img = img,
    label = label,
    size = label_size
  )
}

combine_panels_grid <- function(panel_files, panel_labels, output_file,
                                ncol = 2, panel_width = 3600, panel_height = 2700,
                                label_size = 130) {
  if (length(panel_files) != length(panel_labels)) {
    stop("panel_files and panel_labels must have the same length.")
  }

  missing <- panel_files[!file.exists(panel_files)]
  if (length(missing) > 0) {
    stop("Missing panel files:\n", paste(missing, collapse = "\n"))
  }

  imgs <- Map(
    function(path, lab) {
      make_clean_panel_img(
        path = path,
        label = lab,
        panel_width = panel_width,
        panel_height = panel_height,
        label_size = label_size
      )
    },
    panel_files,
    panel_labels
  )

  rows <- list()
  row_index <- 1
  for (i in seq(1, length(imgs), by = ncol)) {
    row_imgs <- imgs[i:min(i + ncol - 1, length(imgs))]
    if (length(row_imgs) < ncol) {
      blank <- magick::image_blank(panel_width, panel_height, color = "white")
      row_imgs <- c(row_imgs, rep(list(blank), ncol - length(row_imgs)))
    }
    rows[[row_index]] <- magick::image_append(do.call(c, row_imgs), stack = FALSE)
    row_index <- row_index + 1
  }

  combined <- magick::image_append(do.call(c, rows), stack = TRUE)
  combined <- magick::image_background(combined, "white", flatten = TRUE)

  magick::image_write(
    combined,
    path = output_file,
    format = "png",
    density = paste0(FIG_DPI, "x", FIG_DPI)
  )

  output_file
}

common_theme <- theme_manuscript(show_grid = FALSE, legend_position = "right")

cat("\nSTEP 01 finished: setup completed.\n")
cat("\nFigure-style standard applied: Arial, consistent panel tags, titles, axes, borders, margins, and legends.\n")

# STEP 02. RAW BULK PREPROCESSING, TMM-CPM, SYMBOL AGGREGATION, SPLIT SETS

step02_dir <- make_stage_dir(2, "Bulk_Preprocessing_Train_Validation_External")

train1_expr  <- file.path(project_path, "GSE127165_raw_counts_GRCh38.p13_NCBI.csv")
train2_expr  <- file.path(project_path, "GSE142083_raw_counts_GRCh38.p13_NCBI.csv")
test_expr    <- file.path(project_path, "GSE130605_raw_counts_GRCh38.p13_NCBI.csv")

train1_annot <- file.path(project_path, "GSE127165_annot.csv")
train2_annot <- file.path(project_path, "GSE142083_annot.csv")
test_annot   <- file.path(project_path, "GSE130605_annot.csv")

pheno_tr1    <- file.path(project_path, "Pheno_Data_GSE127165.csv")
pheno_tr2    <- file.path(project_path, "Pheno_Data_GSE142083.csv")
pheno_te     <- file.path(project_path, "Pheno_Data_GSE130605.csv")

all_input_files <- c(
  train1_expr, train2_expr, test_expr,
  train1_annot, train2_annot, test_annot,
  pheno_tr1, pheno_tr2, pheno_te
)

missing_input_files <- all_input_files[!file.exists(all_input_files)]
if (length(missing_input_files) > 0) {
  stop("Missing raw bulk input files:\n", paste(missing_input_files, collapse = "\n"))
}

calc_cpm_tmm_nonnegative <- function(counts_csv, out_csv) {
  X_raw <- data.table::fread(counts_csv, colClasses = c(GeneID = "character"))
  data.table::setnames(X_raw, 1, "GeneID")

  cts <- as.matrix(X_raw[, -1, with = FALSE])
  rownames(cts) <- X_raw$GeneID
  mode(cts) <- "numeric"

  y <- edgeR::DGEList(counts = cts)
  y <- edgeR::calcNormFactors(y, method = "TMM")

  CPM <- edgeR::cpm(y, log = FALSE)

  out_dt <- data.table::data.table(GeneID = rownames(CPM), CPM)
  data.table::fwrite(out_dt, out_csv)
  out_csv
}

aggregate_geneid_to_symbol <- function(annot_csv, cpm_csv, out_csv) {
  ann <- data.table::fread(
    annot_csv,
    colClasses = c(GeneID = "character", Symbol = "character")
  )[, .(GeneID, Symbol)]

  dat <- data.table::fread(cpm_csv, colClasses = c(GeneID = "character"))
  M <- merge(ann, dat, by = "GeneID", all.y = TRUE)
  M <- M[!is.na(Symbol) & Symbol != ""]

  agg <- M[, lapply(.SD, mean, na.rm = TRUE),
           by = Symbol,
           .SDcols = setdiff(names(M), c("GeneID", "Symbol"))]

  data.table::setorder(agg, Symbol)
  data.table::fwrite(agg, out_csv)
  out_csv
}

transpose_symbol_matrix <- function(agg_csv, out_csv) {
  A <- data.table::fread(agg_csv)
  data.table::setorder(A, Symbol)

  genes <- A$Symbol
  mat <- as.data.frame(t(A[, -1, with = FALSE]))
  colnames(mat) <- genes
  mat$Sample <- rownames(mat)
  rownames(mat) <- NULL
  mat <- mat[, c("Sample", setdiff(names(mat), "Sample"))]

  data.table::fwrite(as.data.table(mat), out_csv)
  out_csv
}

add_metadata_and_batch <- function(trp_csv, pheno_csv, out_csv, batch_value) {
  X  <- data.table::fread(trp_csv)
  ph <- data.table::fread(pheno_csv)

  if (!("Sample" %in% names(X))) {
    stop("Transposed file is missing Sample column: ", trp_csv)
  }

  X[, Sample := clean_id(Sample)]

  smp_col <- guess_sample_col(ph)
  grp_col <- guess_group_col(ph)

  data.table::setnames(ph, smp_col, "Sample")
  data.table::setnames(ph, grp_col, "group")

  ph[, Sample := clean_id(Sample)]
  ph[, group  := coerce_group_12(group)]
  ph[, batch  := as.integer(batch_value)]

  M <- merge(
    ph[, .(Sample, group, batch)],
    X,
    by = "Sample",
    all.y = TRUE
  )

  M <- M[
    !is.na(Sample) & Sample != "" &
      !is.na(group) & group %in% c(1, 2),
  ]

  M <- M[, c(
    "Sample", "group", "batch",
    setdiff(names(M), c("Sample", "group", "batch"))
  ), with = FALSE]

  data.table::fwrite(M, out_csv)
  out_csv
}

filter_zero_genes_by_rate <- function(input_file, output_file, zero_cutoff = 0) {
  df <- data.table::fread(input_file)

  meta_cols <- c("Sample", "group", "batch")
  gene_cols <- setdiff(names(df), meta_cols)

  mat <- as.matrix(df[, gene_cols, with = FALSE])
  mode(mat) <- "numeric"

  zero_rate <- colMeans(mat == 0, na.rm = TRUE)
  keep_genes <- names(zero_rate)[zero_rate <= zero_cutoff]
  removed_genes <- names(zero_rate)[zero_rate > zero_cutoff]

  if (length(keep_genes) == 0) {
    stop("No genes left after zero filtering. zero_cutoff = ", zero_cutoff)
  }

  df_filtered <- df[, c(meta_cols, keep_genes), with = FALSE]
  data.table::fwrite(df_filtered, output_file)

  report <- data.frame(
    Gene = names(zero_rate),
    ZeroRate = as.numeric(zero_rate),
    Status = ifelse(names(zero_rate) %in% keep_genes, "Kept", "Removed")
  )

  report_file <- sub("\\.csv$", "_zero_filter_report.csv", output_file)
  data.table::fwrite(report, report_file)

  cat("\nZero filtering:", input_file, "\n")
  cat("Original genes:", length(gene_cols), "\n")
  cat("Kept genes:", length(keep_genes), "\n")
  cat("Removed genes:", length(removed_genes), "\n")
  cat("Output:", output_file, "\n")

  output_file
}

s2_tr1_cpm <- file.path(step02_dir, "GSE127165_TMM_CPM_nonnegative.csv")
s2_tr2_cpm <- file.path(step02_dir, "GSE142083_TMM_CPM_nonnegative.csv")
s2_te_cpm  <- file.path(step02_dir, "GSE130605_TMM_CPM_nonnegative.csv")

calc_cpm_tmm_nonnegative(train1_expr, s2_tr1_cpm)
calc_cpm_tmm_nonnegative(train2_expr, s2_tr2_cpm)
calc_cpm_tmm_nonnegative(test_expr,   s2_te_cpm)

s2_tr1_sym <- file.path(step02_dir, "GSE127165_symbol_CPM.csv")
s2_tr2_sym <- file.path(step02_dir, "GSE142083_symbol_CPM.csv")
s2_te_sym  <- file.path(step02_dir, "GSE130605_symbol_CPM.csv")

aggregate_geneid_to_symbol(train1_annot, s2_tr1_cpm, s2_tr1_sym)
aggregate_geneid_to_symbol(train2_annot, s2_tr2_cpm, s2_tr2_sym)
aggregate_geneid_to_symbol(test_annot,   s2_te_cpm,  s2_te_sym)

s2_tr1_trp <- file.path(step02_dir, "GSE127165_transposed_CPM.csv")
s2_tr2_trp <- file.path(step02_dir, "GSE142083_transposed_CPM.csv")
s2_te_trp  <- file.path(step02_dir, "GSE130605_transposed_CPM.csv")

transpose_symbol_matrix(s2_tr1_sym, s2_tr1_trp)
transpose_symbol_matrix(s2_tr2_sym, s2_tr2_trp)
transpose_symbol_matrix(s2_te_sym,  s2_te_trp)

s2_tr1_meta <- file.path(step02_dir, "GSE127165_meta_CPM.csv")
s2_tr2_meta <- file.path(step02_dir, "GSE142083_meta_CPM.csv")
s2_te_meta  <- file.path(step02_dir, "GSE130605_meta_CPM.csv")

add_metadata_and_batch(s2_tr1_trp, pheno_tr1, s2_tr1_meta, batch_value = 1)
add_metadata_and_batch(s2_tr2_trp, pheno_tr2, s2_tr2_meta, batch_value = 2)
add_metadata_and_batch(s2_te_trp,  pheno_te,  s2_te_meta,  batch_value = 3)

S_tr1 <- data.table::fread(s2_tr1_meta)
S_tr2 <- data.table::fread(s2_tr2_meta)

req_cols <- c("Sample", "group", "batch")
common_train_genes <- intersect(
  setdiff(names(S_tr1), req_cols),
  setdiff(names(S_tr2), req_cols)
)

S_tr_merged <- data.table::rbindlist(
  list(
    S_tr1[, c(req_cols, common_train_genes), with = FALSE],
    S_tr2[, c(req_cols, common_train_genes), with = FALSE]
  ),
  use.names = TRUE,
  fill = FALSE
)

S_tr_merged <- S_tr_merged[
  !is.na(Sample) & Sample != "" &
    group %in% c(1, 2) &
    batch %in% c(1, 2)
]

s2_merged_csv <- file.path(step02_dir, "train_discovery_merged_CPM_nonnegative.csv")
data.table::fwrite(S_tr_merged, s2_merged_csv)

s2_merged_nozero_csv <- file.path(step02_dir, "train_discovery_merged_CPM_nonnegative_noZero.csv")

filter_zero_genes_by_rate(
  input_file = s2_merged_csv,
  output_file = s2_merged_nozero_csv,
  zero_cutoff = 0
)

samples <- as.data.frame(data.table::fread(s2_merged_nozero_csv))
rownames(samples) <- samples$Sample
samples <- samples[, setdiff(colnames(samples), "Sample"), drop = FALSE]

cancer_samples <- rownames(samples)[samples$group == 2]
normal_samples <- rownames(samples)[samples$group == 1]

set.seed(123)
n_bal <- min(length(cancer_samples), length(normal_samples))

cancer_bal <- sample(cancer_samples, n_bal)
normal_bal <- sample(normal_samples, n_bal)

n_train <- ceiling(0.7 * n_bal)

training_cancer_samples <- sample(cancer_bal, n_train)
validation_cancer_samples <- setdiff(cancer_bal, training_cancer_samples)

training_normal_samples <- sample(normal_bal, n_train)
validation_normal_samples <- setdiff(normal_bal, training_normal_samples)

training_samples <- c(training_cancer_samples, training_normal_samples)
validation_samples <- c(validation_cancer_samples, validation_normal_samples)

training_set_cpm <- samples[training_samples, , drop = FALSE]
validation_set_cpm <- samples[validation_samples, , drop = FALSE]

training_set_cpm <- data.frame(
  Sample = rownames(training_set_cpm),
  training_set_cpm,
  row.names = NULL,
  check.names = FALSE
)

validation_set_cpm <- data.frame(
  Sample = rownames(validation_set_cpm),
  validation_set_cpm,
  row.names = NULL,
  check.names = FALSE
)

training_file_cpm <- file.path(step02_dir, "training_set_CPM.csv")
validation_file_cpm <- file.path(step02_dir, "validation_set_CPM.csv")

write.csv(training_set_cpm, training_file_cpm, row.names = FALSE)
write.csv(validation_set_cpm, validation_file_cpm, row.names = FALSE)

external_meta <- data.table::fread(s2_te_meta)
kept_discovery_genes <- setdiff(colnames(training_set_cpm), c("Sample", "group", "batch"))
external_common_genes <- intersect(kept_discovery_genes, setdiff(names(external_meta), req_cols))

if (length(external_common_genes) < 2) {
  stop("Too few common genes between discovery train set and external GSE130605.")
}

external_set_cpm <- external_meta[, c(req_cols, external_common_genes), with = FALSE]
external_file_cpm <- file.path(step02_dir, "external_GSE130605_CPM.csv")
data.table::fwrite(external_set_cpm, external_file_cpm)

split_summary <- data.frame(
  Set = c("Training", "Validation", "External_GSE130605"),
  Samples = c(nrow(training_set_cpm), nrow(validation_set_cpm), nrow(external_set_cpm)),
  Normal = c(
    sum(training_set_cpm$group == 1),
    sum(validation_set_cpm$group == 1),
    sum(external_set_cpm$group == 1)
  ),
  Tumor = c(
    sum(training_set_cpm$group == 2),
    sum(validation_set_cpm$group == 2),
    sum(external_set_cpm$group == 2)
  ),
  Genes = c(
    length(setdiff(colnames(training_set_cpm), req_cols)),
    length(setdiff(colnames(validation_set_cpm), req_cols)),
    length(setdiff(colnames(external_set_cpm), req_cols))
  )
)

write.csv(
  split_summary,
  file.path(step02_dir, "Step_02_train_validation_external_summary.csv"),
  row.names = FALSE
)

cat("\nSTEP 02 finished: bulk preprocessing and train/validation/external CPM files generated.\n")
print(split_summary)

# STEP 03. LOG2-CPM, DIFFERENTIAL EXPRESSION, AND MANUSCRIPT FIGURE 2

step03_dir <- make_stage_dir(3, "Log2CPM_DEG_Figure2_no_batch_correction")

DEG_FDR_THRESH <- 0.01
DEG_LOGFC_MIN  <- 2

prepare_log2_matrix <- function(df) {
  genes <- setdiff(colnames(df), c("Sample", "group", "batch"))
  mat <- to_numeric_df(df[, genes, drop = FALSE])
  rownames(mat) <- df$Sample
  log2(as.matrix(mat) + 1)
}

training_set_cpm <- as.data.frame(data.table::fread(training_file_cpm))
validation_set_cpm <- as.data.frame(data.table::fread(validation_file_cpm))
external_set_cpm <- as.data.frame(data.table::fread(external_file_cpm))

common_step03_genes <- Reduce(
  intersect,
  list(
    setdiff(colnames(training_set_cpm), c("Sample", "group", "batch")),
    setdiff(colnames(validation_set_cpm), c("Sample", "group", "batch")),
    setdiff(colnames(external_set_cpm), c("Sample", "group", "batch"))
  )
)

training_set_cpm <- training_set_cpm[
  ,
  c("Sample", "group", "batch", common_step03_genes),
  drop = FALSE
]

validation_set_cpm <- validation_set_cpm[
  ,
  c("Sample", "group", "batch", common_step03_genes),
  drop = FALSE
]

external_set_cpm <- external_set_cpm[
  ,
  c("Sample", "group", "batch", common_step03_genes),
  drop = FALSE
]

train_log2 <- prepare_log2_matrix(training_set_cpm)
valid_log2 <- prepare_log2_matrix(validation_set_cpm)
test_log2  <- prepare_log2_matrix(external_set_cpm)

train_group <- as.integer(training_set_cpm$group)
valid_group <- as.integer(validation_set_cpm$group)
test_group  <- as.integer(external_set_cpm$group)

train_batch <- as.integer(training_set_cpm$batch)
valid_batch <- as.integer(validation_set_cpm$batch)
test_batch  <- as.integer(external_set_cpm$batch)

train_log2_df <- data.frame(
  Sample = rownames(train_log2),
  group = train_group,
  batch = train_batch,
  as.data.frame(train_log2, check.names = FALSE),
  check.names = FALSE
)

valid_log2_df <- data.frame(
  Sample = rownames(valid_log2),
  group = valid_group,
  batch = valid_batch,
  as.data.frame(valid_log2, check.names = FALSE),
  check.names = FALSE
)

external_log2_df <- data.frame(
  Sample = rownames(test_log2),
  group = test_group,
  batch = test_batch,
  as.data.frame(test_log2, check.names = FALSE),
  check.names = FALSE
)

train_log2_file <- file.path(step03_dir, "train_log2CPM.csv")
valid_log2_file <- file.path(step03_dir, "valid_log2CPM.csv")
external_log2_file <- file.path(step03_dir, "external_GSE130605_log2CPM.csv")

write.csv(train_log2_df, train_log2_file, row.names = FALSE)
write.csv(valid_log2_df, valid_log2_file, row.names = FALSE)
write.csv(external_log2_df, external_log2_file, row.names = FALSE)

Tr <- data.table::fread(train_log2_file)

genes_bulk <- setdiff(names(Tr), c("Sample", "group", "batch"))
E <- as.matrix(to_numeric_df(Tr[, ..genes_bulk]))
rownames(E) <- Tr$Sample

grp <- factor(
  Tr$group,
  levels = c(1, 2),
  labels = c("Normal", "Tumor")
)
grp <- relevel(grp, ref = "Normal")

design <- model.matrix(~ grp)

fit <- limma::lmFit(t(E), design)
fit <- limma::eBayes(fit, trend = TRUE)

TT_full <- limma::topTable(
  fit,
  coef = 2,
  number = Inf,
  adjust.method = "BH",
  sort.by = "B"
)

TT_full$Gene.symbol <- rownames(TT_full)
TT_full <- TT_full[, c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")]

TT_deg <- TT_full %>%
  dplyr::filter(
    adj.P.Val < DEG_FDR_THRESH,
    abs(logFC) > DEG_LOGFC_MIN
  )

Bulk_UP <- TT_deg %>%
  dplyr::filter(logFC > DEG_LOGFC_MIN) %>%
  dplyr::pull(Gene.symbol) %>%
  unique()

Bulk_DOWN <- TT_deg %>%
  dplyr::filter(logFC < -DEG_LOGFC_MIN) %>%
  dplyr::pull(Gene.symbol) %>%
  unique()

write.csv(TT_full, file.path(step03_dir, "DE_topTable_full.csv"), row.names = FALSE)
write.csv(TT_deg, file.path(step03_dir, "DEG_logFC2_adjP0.01.csv"), row.names = FALSE)
writeLines(Bulk_UP, file.path(step03_dir, "DEG_up_logFC2_adjP0.01.txt"))
writeLines(Bulk_DOWN, file.path(step03_dir, "DEG_down_logFC2_adjP0.01.txt"))

# FIGURE 2A: VOLCANO PLOT

TT_plot <- TT_full %>%
  dplyr::mutate(
    neglog10FDR = -log10(
      pmax(adj.P.Val, 1e-300)
    ),

    Status = dplyr::case_when(
      adj.P.Val < DEG_FDR_THRESH &
        logFC > DEG_LOGFC_MIN ~ "Up",

      adj.P.Val < DEG_FDR_THRESH &
        logFC < -DEG_LOGFC_MIN ~ "Down",

      TRUE ~ "NS"
    ),

    Status = factor(
      Status,
      levels = c("Up", "Down", "NS")
    )

  )

p_fig2a <- ggplot(
  TT_plot,
  aes(
    x = logFC,
    y = neglog10FDR,
    color = Status
  )
) +
  geom_point(
    size = 1.15,
    alpha = 0.84
  ) +
  geom_vline(
    xintercept = c(
      -DEG_LOGFC_MIN,
      DEG_LOGFC_MIN
    ),
    linetype = "dashed",
    linewidth = GEOM_LWD,
    color = "black"
  ) +
  geom_hline(
    yintercept = -log10(DEG_FDR_THRESH),
    linetype = "dashed",
    linewidth = GEOM_LWD,
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      Up = COL_UP,
      Down = COL_DOWN,
      NS = COL_NS
    )
  ) +
  labs(
    title = NULL,
    tag = "a",
    x = "log2 fold change",
    y = "-log10 adjusted P-value",
    color = NULL
  ) +
  theme_manuscript(
    show_grid = FALSE,
    legend_position = "bottom"
  ) +
  theme(
    legend.direction = "horizontal"
  )

save_plot_all_formats(
  plot_obj = p_fig2a,
  filename_stem = "Figure_02A_DEG_volcano",
  dir_path = FIGURE_DIR,
  width = FIG_SINGLE_W,
  height = FIG_SINGLE_H
)

# FIGURE 2B: ALL-SIGNIFICANT-DEG HEATMAP

all_deg_gene_pool <- TT_deg %>%
  dplyr::filter(
    Gene.symbol %in% colnames(E)
  ) %>%
  dplyr::mutate(
    Regulation = ifelse(
      logFC > 0,
      "Up",
      "Down"
    )
  ) %>%
  dplyr::arrange(
    Regulation,
    adj.P.Val,
    dplyr::desc(abs(logFC))
  )

all_deg_genes <- unique(
  all_deg_gene_pool$Gene.symbol
)

if (length(all_deg_genes) < 2) {
  stop("There are not enough significant DEGs to construct the heatmap.")
}

hm_mat <- t(
  E[
    ,
    all_deg_genes,
    drop = FALSE
  ]
)

hm_scaled <- t(
  scale(
    t(hm_mat)
  )
)

hm_scaled[!is.finite(hm_scaled)] <- 0

sample_order <- order(grp)

hm_scaled <- hm_scaled[
  ,
  sample_order,
  drop = FALSE
]

grp_ordered <- grp[sample_order]

ann_col <- data.frame(
  Group = factor(
    grp_ordered,
    levels = c("Normal", "Tumor")
  )
)

rownames(ann_col) <- colnames(hm_scaled)

ann_row <- data.frame(
  Regulation = ifelse(
    rownames(hm_scaled) %in% Bulk_UP,
    "Up",
    "Down"
  )
)

rownames(ann_row) <- rownames(hm_scaled)

normal_count <- sum(
  ann_col$Group == "Normal"
)

gap_cols <- NULL

if (
  normal_count > 0 &&
  sum(ann_col$Group == "Tumor") > 0
) {
  gap_cols <- normal_count
}

heatmap_colors <- grDevices::colorRampPalette(
  c(
    COL_DOWN,
    "white",
    COL_UP
  )
)(101)

HEATMAP_GROUP_NORMAL <- "#4C78A8"
HEATMAP_GROUP_TUMOR <- "#D64F7A"
HEATMAP_REGULATION_DOWN <- "#2A9D8F"
HEATMAP_REGULATION_UP <- "#F4A261"

heatmap_annotation_colors <- list(
  Group = c(
    Normal = HEATMAP_GROUP_NORMAL,
    Tumor = HEATMAP_GROUP_TUMOR
  ),
  Regulation = c(
    Down = HEATMAP_REGULATION_DOWN,
    Up = HEATMAP_REGULATION_UP
  )
)

p_heatmap_obj <- pheatmap::pheatmap(
  mat = hm_scaled,
  color = heatmap_colors,
  border_color = NA,

  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = gap_cols,

  show_rownames = FALSE,
  show_colnames = FALSE,

  annotation_col = ann_col,
  annotation_row = ann_row,
  annotation_colors = heatmap_annotation_colors,

  annotation_names_row = FALSE,
  annotation_names_col = FALSE,

  fontsize = HEATMAP_FONT_PT,
  fontsize_row = HEATMAP_ROW_FONT_PT,
  fontsize_col = HEATMAP_COL_FONT_PT,

  main = "",
  silent = TRUE
)

heatmap_grob <- set_grob_font_family(
  p_heatmap_obj$gtable,
  FONT_FAMILY
)

heatmap_grob_with_tag <- grid::grobTree(
  heatmap_grob,
  grid::rectGrob(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(0.5, "npc"),
    width = grid::unit(0.995, "npc"),
    height = grid::unit(0.995, "npc"),
    gp = grid::gpar(
      fill = NA,
      col = "black",
      lwd = PANEL_BORDER_LWD
    )
  ),
  grid::textGrob(
    label = "b",
    x = grid::unit(PANEL_TAG_POSITION[1], "npc"),
    y = grid::unit(PANEL_TAG_POSITION[2], "npc"),
    just = c("left", "top"),
    gp = grid::gpar(
      fontfamily = FONT_FAMILY,
      fontface = "bold",
      fontsize = PANEL_TAG_PT,
      col = "black"
    )
  )
)

save_grob_all_formats(
  grob_obj = heatmap_grob_with_tag,
  filename_stem = "Figure_02B_DEG_heatmap_all_DEGs",
  dir_path = FIGURE_DIR,
  width = 10.80,
  height = 8.20
)

p_fig2b <- patchwork::wrap_elements(
  full = heatmap_grob_with_tag
)

# MANUSCRIPT FIGURE 2

fig2_combined <- patchwork::wrap_plots(
  list(
    p_fig2a,
    p_fig2b
  ),
  ncol = 2,
  widths = c(0.62, 1.38)
)

save_plot_all_formats(
  plot_obj = fig2_combined,
  filename_stem = "Figure_02_DEG_analysis_LavenderPurple_Fuchsia",
  dir_path = FIGURE_DIR,
  width = FIGURE2_COMBINED_W,
  height = FIGURE2_COMBINED_H
)

cat("\nSTEP 03 finished: log2-CPM, DEG analysis, and standardized Figure 2 generated.\n")
cat("No batch-effect correction was applied.\n")
cat("Total DEGs:", nrow(TT_deg), "\n")
cat("Upregulated DEGs:", length(Bulk_UP), "\n")
cat("Downregulated DEGs:", length(Bulk_DOWN), "\n")

# STEP 04. WGCNA AND MANUSCRIPT FIGURE 3

step04_dir <- make_stage_dir(4, "WGCNA_Figure3")

WGCNA_TOP_MAD <- 8000
WGCNA_TARGET_R2 <- 0.92
WGCNA_MIN_MEANK <- 10
WGCNA_POWERS <- 1:20
WGCNA_MM_MIN <- 0.55
WGCNA_GS_MIN <- 0.55
WGCNA_RULE <- "AND"

TrW <- data.table::fread(train_log2_file)

traitDF <- data.frame(
  Tumor  = as.integer(TrW$group == 2),
  Normal = as.integer(TrW$group == 1),
  row.names = TrW$Sample
)

gene_cols_w <- setdiff(names(TrW), c("Sample", "group", "batch"))

datExpr0 <- as.data.frame(
  lapply(TrW[, ..gene_cols_w], function(x) as.numeric(as.character(x))),
  check.names = FALSE
)
rownames(datExpr0) <- TrW$Sample

gsg <- WGCNA::goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  traitDF  <- traitDF[rownames(datExpr0), , drop = FALSE]
}

sampleTree <- hclust(dist(datExpr0), method = "average")

mad_vals <- apply(datExpr0, 2, mad, na.rm = TRUE)
mad_vals <- mad_vals[is.finite(mad_vals) & !is.na(mad_vals) & mad_vals > 0]

topN_mad_final <- min(WGCNA_TOP_MAD, length(mad_vals))
keep_genes_w <- names(sort(mad_vals, decreasing = TRUE))[1:topN_mad_final]
datExpr <- datExpr0[, keep_genes_w, drop = FALSE]

sft <- WGCNA::pickSoftThreshold(
  datExpr,
  powerVector = WGCNA_POWERS,
  networkType = "signed",
  corFnc = "cor",
  verbose = 3
)

fitR2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
meanK <- sft$fitIndices[, 5]

eligible <- which(fitR2 >= WGCNA_TARGET_R2 & meanK >= WGCNA_MIN_MEANK)
softPower <- if (length(eligible) > 0) {
  sft$fitIndices[min(eligible), 1]
} else {
  sft$fitIndices[which.max(fitR2), 1]
}

sft_table <- data.frame(
  Power = sft$fitIndices[, 1],
  SignedR2 = fitR2,
  MeanK = meanK
)

write.csv(sft_table, file.path(step04_dir, "pickSoftThreshold_full_table.csv"), row.names = FALSE)
writeLines(paste0("Chosen softPower = ", softPower), file.path(step04_dir, "softPower.txt"))

net <- WGCNA::blockwiseModules(
  datExpr,
  power = softPower,
  networkType = "signed",
  TOMType = "signed",
  corType = "pearson",
  minModuleSize = 60,
  mergeCutHeight = 0.18,
  numericLabels = FALSE,
  pamRespectsDendro = TRUE,
  maxBlockSize = ncol(datExpr),
  verbose = 3
)

moduleColors_all <- net$colors
names(moduleColors_all) <- colnames(datExpr)

MEs <- WGCNA::orderMEs(net$MEs)

modTraitCor <- cor(MEs, traitDF, use = "p", method = "pearson")
modTraitPval <- WGCNA::corPvalueStudent(modTraitCor, nrow(datExpr))

write.csv(modTraitCor, file.path(step04_dir, "module_trait_cor.csv"))
write.csv(modTraitPval, file.path(step04_dir, "module_trait_pvalue.csv"))

valid_me <- rownames(modTraitCor)
valid_me <- valid_me[valid_me != "MEgrey"]

tumor_cor <- modTraitCor[valid_me, "Tumor"]
tumor_cor[tumor_cor <= 0] <- NA

if (all(is.na(tumor_cor))) {
  stop("No positively tumor-associated module was found.")
}

bestME <- valid_me[which.max(tumor_cor)]
targetModule <- sub("^ME", "", bestME)
targetGenes <- names(moduleColors_all)[moduleColors_all == targetModule]

GS_all <- as.numeric(cor(datExpr, traitDF$Tumor, use = "p", method = "pearson"))
MM_all <- as.numeric(cor(datExpr, MEs[, bestME], use = "p", method = "pearson"))

names(GS_all) <- colnames(datExpr)
names(MM_all) <- colnames(datExpr)

mod_table <- data.frame(
  Gene = targetGenes,
  GS_Tumor = GS_all[targetGenes],
  MM_Module = MM_all[targetGenes],
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    absGS = abs(GS_Tumor),
    absMM = abs(MM_Module)
  ) %>%
  dplyr::arrange(dplyr::desc(absMM), dplyr::desc(absGS))

if (toupper(WGCNA_RULE) == "AND") {
  WGCNA_Strong_Tumor_Module_Genes <- mod_table %>%
    dplyr::filter(absMM >= WGCNA_MM_MIN, absGS >= WGCNA_GS_MIN) %>%
    dplyr::pull(Gene) %>%
    unique()
} else {
  WGCNA_Strong_Tumor_Module_Genes <- mod_table %>%
    dplyr::filter(absMM >= WGCNA_MM_MIN | absGS >= WGCNA_GS_MIN) %>%
    dplyr::pull(Gene) %>%
    unique()
}

write.csv(mod_table, file.path(step04_dir, paste0(targetModule, "_GS_MM_full.csv")), row.names = FALSE)
write.csv(data.frame(Gene = WGCNA_Strong_Tumor_Module_Genes), file.path(step04_dir, "WGCNA_feature_genes.csv"), row.names = FALSE)

saveRDS(
  list(
    datExpr0 = datExpr0,
    datExpr = datExpr,
    traitDF = traitDF,
    sampleTree = sampleTree,
    sft = sft,
    sft_table = sft_table,
    softPower = softPower,
    net = net,
    moduleColors_all = moduleColors_all,
    MEs = MEs,
    modTraitCor = modTraitCor,
    modTraitPval = modTraitPval,
    bestME = bestME,
    targetModule = targetModule,
    mod_table = mod_table,
    WGCNA_Strong_Tumor_Module_Genes = WGCNA_Strong_Tumor_Module_Genes
  ),
  file.path(step04_dir, "WGCNA_analysis_objects.rds")
)

save_base_plot_all_formats <- function(draw_fun,
                                       filename_stem,
                                       dir_path,
                                       width,
                                       height) {

  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  png_file <- file.path(dir_path, paste0(filename_stem, ".png"))
  tiff_file <- file.path(dir_path, paste0(filename_stem, ".tiff"))
  pdf_file <- file.path(dir_path, paste0(filename_stem, ".pdf"))

  grDevices::png(
    filename = png_file,
    width = width,
    height = height,
    units = "in",
    res = FIG_DPI,
    bg = FIG_BACKGROUND
  )
  draw_fun()
  grDevices::dev.off()

  grDevices::tiff(
    filename = tiff_file,
    width = width,
    height = height,
    units = "in",
    res = FIG_DPI,
    compression = "lzw",
    bg = FIG_BACKGROUND
  )
  draw_fun()
  grDevices::dev.off()

  if (capabilities("cairo")) {
    grDevices::cairo_pdf(
      filename = pdf_file,
      width = width,
      height = height,
      family = FONT_FAMILY
    )
  } else {
    grDevices::pdf(
      file = pdf_file,
      width = width,
      height = height,
      family = FONT_FAMILY
    )
  }

  draw_fun()
  grDevices::dev.off()

  invisible(
    list(
      PNG = png_file,
      TIFF = tiff_file,
      PDF = pdf_file
    )
  )
}

traitColors <- data.frame(
  Tumor = ifelse(traitDF$Tumor == 1, COL_TUMOR, "white"),
  Normal = ifelse(traitDF$Normal == 1, COL_NORMAL, "white"),
  row.names = rownames(traitDF)
)

draw_fig3a <- function() {

  par(
    family = FONT_FAMILY,
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black",
    lwd = AXIS_LWD,
    cex.axis = BASE_AXIS_CEX,
    cex.lab = BASE_AXIS_TITLE_CEX,
    cex.main = BASE_AXIS_TITLE_CEX,
    las = 1,
    mar = c(1.2, 4.5, 1.0, 1.0)
  )

  WGCNA::plotDendroAndColors(
    sampleTree,
    colors = traitColors,
    groupLabels = colnames(traitColors),
    main = "",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = FALSE,
    cex.colorLabels = BASE_AXIS_CEX,
    cex.dendroLabels = 0.60,
    colorHeight = 0.18
  )
}

save_base_plot_all_formats(
  draw_fun = draw_fig3a,
  filename_stem = "Figure_03A_WGCNA_sample_dendrogram_trait_heatmap",
  dir_path = step04_dir,
  width = 13.50,
  height = 5.80
)

fig3a_file <- file.path(
  step04_dir,
  "Figure_03A_WGCNA_sample_dendrogram_trait_heatmap.png"
)

sft_r2_y <- max(sft_table$SignedR2, na.rm = TRUE)
sft_k_y  <- max(sft_table$MeanK, na.rm = TRUE)

p_sft1 <- ggplot(sft_table, aes(Power, SignedR2)) +
  geom_line(
    color = COL_TUMOR,
    linewidth = GEOM_LWD * 1.8
  ) +
  geom_point(
    color = COL_TUMOR,
    size = 2.7
  ) +
  geom_vline(
    xintercept = softPower,
    linetype = "dashed",
    linewidth = GEOM_LWD * 1.15,
    color = "black"
  ) +
  annotate(
    "text",
    x = softPower,
    y = min(sft_r2_y, 0.96),
    label = paste0("β = ", softPower),
    family = FONT_FAMILY,
    size = 4.0,
    fontface = "bold",
    vjust = -0.55,
    hjust = 1.12,
    color = "black"
  ) +
  labs(
    x = "Soft-threshold power",
    y = "Signed R²"
  ) +
  theme_manuscript(
    show_grid = FALSE,
    legend_position = "none"
  )

p_sft2 <- ggplot(sft_table, aes(Power, MeanK)) +
  geom_line(
    color = COL_TUMOR,
    linewidth = GEOM_LWD * 1.8
  ) +
  geom_point(
    color = COL_TUMOR,
    size = 2.7
  ) +
  geom_hline(
    yintercept = WGCNA_MIN_MEANK,
    linetype = "dashed",
    linewidth = GEOM_LWD * 1.15,
    color = "black"
  ) +
  geom_vline(
    xintercept = softPower,
    linetype = "dashed",
    linewidth = GEOM_LWD * 1.15,
    color = "black"
  ) +
  annotate(
    "text",
    x = softPower,
    y = sft_k_y * 0.92,
    label = paste0("β = ", softPower),
    family = FONT_FAMILY,
    size = 4.0,
    fontface = "bold",
    hjust = 1.12,
    color = "black"
  ) +
  annotate(
    "text",
    x = min(sft_table$Power) + 1.2,
    y = WGCNA_MIN_MEANK + (sft_k_y * 0.06),
    label = paste0("Mean connectivity = ", WGCNA_MIN_MEANK),
    family = FONT_FAMILY,
    size = 3.45,
    fontface = "bold",
    hjust = 0,
    color = "black"
  ) +
  labs(
    x = "Soft-threshold power",
    y = "Mean connectivity"
  ) +
  theme_manuscript(
    show_grid = FALSE,
    legend_position = "none"
  )

fig3b_plot <- p_sft1 + p_sft2 +
  patchwork::plot_layout(ncol = 2)

save_plot_all_formats(
  plot_obj = fig3b_plot,
  filename_stem = "Figure_03B_WGCNA_soft_threshold",
  dir_path = step04_dir,
  width = 9.60,
  height = 4.70
)

fig3b_file <- file.path(
  step04_dir,
  "Figure_03B_WGCNA_soft_threshold.png"
)

geneTree <- net$dendrograms[[1]]
moduleColors_block <- moduleColors_all[net$blockGenes[[1]]]

draw_fig3c <- function() {

  par(
    family = FONT_FAMILY,
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black",
    lwd = AXIS_LWD,
    cex.axis = BASE_AXIS_CEX,
    cex.lab = BASE_AXIS_TITLE_CEX,
    cex.main = BASE_AXIS_TITLE_CEX,
    las = 1,
    mar = c(1.2, 4.5, 0.8, 1.0)
  )

  WGCNA::plotDendroAndColors(
    geneTree,
    colors = moduleColors_block,
    groupLabels = "Module colors",
    main = "",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = FALSE,
    guideHang = 0.05,
    cex.colorLabels = BASE_AXIS_CEX,
    colorHeight = 0.16
  )
}

save_base_plot_all_formats(
  draw_fun = draw_fig3c,
  filename_stem = "Figure_03C_WGCNA_gene_dendrogram_module_colors",
  dir_path = step04_dir,
  width = 8.20,
  height = 5.80
)

fig3c_file <- file.path(
  step04_dir,
  "Figure_03C_WGCNA_gene_dendrogram_module_colors.png"
)

textMatrix <- paste0(
  formatC(modTraitCor, format = "f", digits = 2),
  "\n(",
  formatC(modTraitPval, format = "e", digits = 1),
  ")"
)

dim(textMatrix) <- dim(modTraitCor)

wgcna_heatmap_colors <- grDevices::colorRampPalette(
  c(COL_NORMAL, "white", COL_TUMOR)
)(50)

draw_fig3d <- function() {

  par(
    family = FONT_FAMILY,
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black",
    lwd = AXIS_LWD,
    cex.axis = BASE_AXIS_CEX,
    cex.lab = BASE_AXIS_TITLE_CEX,
    cex.main = BASE_AXIS_TITLE_CEX,
    las = 1,
    mar = c(5.8, 8.3, 1.4, 3.2)
  )

  WGCNA::labeledHeatmap(
    Matrix = modTraitCor,
    xLabels = colnames(traitDF),
    yLabels = rownames(modTraitCor),
    ySymbols = rownames(modTraitCor),
    colorLabels = FALSE,
    colors = wgcna_heatmap_colors,
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.80,
    cex.lab.x = 1.25,
    cex.lab.y = 1.04,
    zlim = c(-1, 1),
    main = ""
  )
}

save_base_plot_all_formats(
  draw_fun = draw_fig3d,
  filename_stem = "Figure_03D_WGCNA_module_trait_heatmap",
  dir_path = step04_dir,
  width = 8.50,
  height = 10.20
)

fig3d_file <- file.path(
  step04_dir,
  "Figure_03D_WGCNA_module_trait_heatmap.png"
)

mod_table_plot <- mod_table %>%
  dplyr::mutate(
    Gene_set = ifelse(
      Gene %in% WGCNA_Strong_Tumor_Module_Genes,
      "Selected feature genes",
      "Other module genes"
    ),
    Gene_set = factor(
      Gene_set,
      levels = c("Other module genes", "Selected feature genes")
    )
  )

p_fig3e <- ggplot(
  mod_table_plot,
  aes(absMM, absGS, color = Gene_set)
) +
  geom_point(
    size = 1.8,
    alpha = 0.80
  ) +
  geom_vline(
    xintercept = WGCNA_MM_MIN,
    linetype = "dashed",
    linewidth = GEOM_LWD,
    color = "black"
  ) +
  geom_hline(
    yintercept = WGCNA_GS_MIN,
    linetype = "dashed",
    linewidth = GEOM_LWD,
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "Other module genes" = COL_NS,
      "Selected feature genes" = COL_TUMOR
    )
  ) +
  labs(
    x = paste0("|Module membership in ", targetModule, "|"),
    y = "|Gene significance for tumor|",
    color = NULL
  ) +
  theme_manuscript(
    show_grid = FALSE,
    legend_position = "top"
  ) +
  theme(
    legend.direction = "horizontal"
  )

save_plot_all_formats(
  plot_obj = p_fig3e,
  filename_stem = "Figure_03E_WGCNA_MM_GS_scatter",
  dir_path = step04_dir,
  width = 6.20,
  height = 5.80
)

fig3e_file <- file.path(
  step04_dir,
  "Figure_03E_WGCNA_MM_GS_scatter.png"
)

add_wgcna_panel_label <- function(img,
                                  label,
                                  label_size = PANEL_TAG_RASTER_SIZE,
                                  offset_x = PANEL_TAG_RASTER_OFFSET_X,
                                  offset_y = PANEL_TAG_RASTER_OFFSET_Y) {
  add_panel_label_img(
    img = img,
    label = label,
    size = label_size,
    offset_x = offset_x,
    offset_y = offset_y
  )
}

make_wgcna_panel_img <- function(path,
                                 label,
                                 panel_width,
                                 panel_height,
                                 label_size = PANEL_TAG_RASTER_SIZE) {
  img <- magick::image_read(path)
  img <- magick::image_background(img, FIG_BACKGROUND, flatten = TRUE)
  img <- magick::image_trim(img)

  img <- magick::image_resize(
    img,
    geometry = paste0(
      panel_width - PANEL_CONTENT_INSET_PX, "x",
      panel_height - PANEL_CONTENT_VERTICAL_INSET_PX, ">"
    )
  )

  img <- magick::image_extent(
    img,
    geometry = paste0(panel_width, "x", panel_height),
    gravity = "center",
    color = FIG_BACKGROUND
  )

  img <- magick::image_border(
    img,
    color = "black",
    geometry = paste0(PANEL_BORDER_RASTER_PX, "x", PANEL_BORDER_RASTER_PX)
  )

  add_wgcna_panel_label(
    img = img,
    label = label,
    label_size = label_size
  )
}

img_a <- make_wgcna_panel_img(
  path = fig3a_file,
  label = "a",
  panel_width = 7000,
  panel_height = 2550
)

img_b <- make_wgcna_panel_img(
  path = fig3b_file,
  label = "b",
  panel_width = 3500,
  panel_height = 2450
)

img_c <- make_wgcna_panel_img(
  path = fig3c_file,
  label = "c",
  panel_width = 3500,
  panel_height = 2450
)

img_d <- make_wgcna_panel_img(
  path = fig3d_file,
  label = "d",
  panel_width = 3500,
  panel_height = 4300
)

img_e <- make_wgcna_panel_img(
  path = fig3e_file,
  label = "e",
  panel_width = 3500,
  panel_height = 4300
)

row1 <- img_a
row2 <- magick::image_append(c(img_b, img_c), stack = FALSE)
row3 <- magick::image_append(c(img_d, img_e), stack = FALSE)

fig3_combined_img <- magick::image_append(
  c(row1, row2, row3),
  stack = TRUE
)

fig3_combined_img <- magick::image_background(
  fig3_combined_img,
  FIG_BACKGROUND,
  flatten = TRUE
)

FIG3_FINAL_HEIGHT_CM <- 19.0
FIG3_FINAL_HEIGHT_PX <- round((FIG3_FINAL_HEIGHT_CM / 2.54) * FIG_DPI)

fig3_combined_img <- magick::image_resize(
  fig3_combined_img,
  geometry = paste0("x", FIG3_FINAL_HEIGHT_PX)
)

fig3_png_file <- file.path(
  step04_dir,
  "Figure_03_WGCNA_combined.png"
)

fig3_tiff_file <- file.path(
  step04_dir,
  "Figure_03_WGCNA_combined.tiff"
)

fig3_pdf_file <- file.path(
  step04_dir,
  "Figure_03_WGCNA_combined.pdf"
)

magick::image_write(
  fig3_combined_img,
  path = fig3_png_file,
  format = "png",
  density = paste0(FIG_DPI, "x", FIG_DPI)
)

magick::image_write(
  fig3_combined_img,
  path = fig3_tiff_file,
  format = "tiff",
  density = paste0(FIG_DPI, "x", FIG_DPI)
)

try(
  magick::image_write(
    fig3_combined_img,
    path = fig3_pdf_file,
    format = "pdf",
    density = paste0(FIG_DPI, "x", FIG_DPI)
  ),
  silent = TRUE
)

save_manuscript_copy(
  fig3_png_file,
  "Figure_03_WGCNA_combined.png"
)

cat("\nSTEP 04 finished: WGCNA and standardized Figure 3 generated.\n")
cat("Target module:", targetModule, "\n")
cat("WGCNA feature genes:", length(WGCNA_Strong_Tumor_Module_Genes), "\n")

# STEP 05. TRIPLE OVERLAP USING BULK-UP, WGCNA, AND PRECOMPUTED scRNA MARKERS

if (!exists("results_path", inherits = FALSE)) {
  results_path <- "D:/LSCC/Results_LSCC/ML"
}

if (!exists("scrna_figdir", inherits = FALSE)) {
  scrna_figdir <- "D:/LSCC/ScRNAseq_Results/GSE206332/Results/figures"
}

if (!exists("FIG_DPI", inherits = FALSE)) {
  FIG_DPI <- 600
}

dir.create(results_path, recursive = TRUE, showWarnings = FALSE)

if (!exists("make_stage_dir", mode = "function", inherits = FALSE)) {
  make_stage_dir <- function(step_number, step_name) {
    dir.create(results_path, recursive = TRUE, showWarnings = FALSE)
    results_path
  }
}

if (!exists("clean_gene_vector", mode = "function", inherits = FALSE)) {
  clean_gene_vector <- function(x) {
    x <- unique(trimws(as.character(x)))
    x[!is.na(x) & x != ""]
  }
}

if (!exists("step03_dir", inherits = FALSE)) {
  step03_dir <- results_path
}

if (!exists("step04_dir", inherits = FALSE)) {
  step04_dir <- results_path
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required for Step 05.")
}

if (!requireNamespace("ggvenn", quietly = TRUE)) {
  stop("Package 'ggvenn' is required for Step 05.")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required for Step 05.")
}

if (.Platform$OS.type == "windows") {
  try(
    grDevices::windowsFonts(
      Arial = grDevices::windowsFont("Arial")
    ),
    silent = TRUE
  )
}

if (!exists("FONT_FAMILY", inherits = FALSE)) {
  FONT_FAMILY <- "Arial"
}

if (!exists("FIG_BACKGROUND", inherits = FALSE)) {
  FIG_BACKGROUND <- "white"
}

if (!exists("manuscript_figures_path", inherits = FALSE)) {
  manuscript_figures_path <- results_path
}

if (!exists("FIGURE_DIR", inherits = FALSE)) {
  FIGURE_DIR <- results_path
}

if (!exists("HAS_GBM", inherits = FALSE)) {
  HAS_GBM <- requireNamespace("gbm", quietly = TRUE)
}

downstream_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggvenn",
  "glmnet",
  "pROC",
  "cowplot",
  "magick",
  "randomForest",
  "e1071"
)

missing_downstream_packages <- downstream_packages[
  !vapply(
    downstream_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_downstream_packages) > 0) {
  stop(
    "The following packages are required for Steps 05-09:
",
    paste(missing_downstream_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggvenn)
  library(glmnet)
  library(pROC)
  library(cowplot)
  library(magick)
  library(randomForest)
  library(e1071)
})

if (!exists("save_fig", mode = "function", inherits = FALSE)) {
  save_fig <- function(plot_obj,
                       filename,
                       dir_path,
                       width = 8,
                       height = 6,
                       dpi = FIG_DPI) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

    ggplot2::ggsave(
      filename = file.path(dir_path, filename),
      plot = plot_obj,
      width = width,
      height = height,
      dpi = dpi,
      bg = FIG_BACKGROUND,
      limitsize = FALSE
    )
  }
}

if (!exists("save_manuscript_copy", mode = "function", inherits = FALSE)) {
  save_manuscript_copy <- function(from_file, figure_name) {
    if (!file.exists(from_file)) {
      return(invisible(NULL))
    }

    to_file <- file.path(manuscript_figures_path, figure_name)

    if (identical(
      normalizePath(from_file, winslash = "/", mustWork = FALSE),
      normalizePath(to_file, winslash = "/", mustWork = FALSE)
    )) {
      return(invisible(to_file))
    }

    file.copy(from = from_file, to = to_file, overwrite = TRUE)
    invisible(to_file)
  }
}

if (!exists("add_panel_label_img", mode = "function", inherits = FALSE)) {
  add_panel_label_img <- function(img,
                                  label,
                                  size = PANEL_TAG_RASTER_SIZE,
                                  offset_x = PANEL_TAG_RASTER_OFFSET_X,
                                  offset_y = PANEL_TAG_RASTER_OFFSET_Y) {
    magick::image_annotate(
      img,
      text = label,
      size = size,
      font = FONT_FAMILY,
      gravity = "northwest",
      location = paste0("+", offset_x, "+", offset_y),
      weight = 700,
      color = "black"
    )
  }
}

if (!exists("make_clean_panel_img", mode = "function", inherits = FALSE)) {
  make_clean_panel_img <- function(path,
                                   label,
                                   panel_width = 3600,
                                   panel_height = 2850,
                                   label_size = PANEL_TAG_RASTER_SIZE) {
    img <- magick::image_read(path)
    img <- magick::image_background(img, FIG_BACKGROUND, flatten = TRUE)
    img <- magick::image_trim(img)

    img <- magick::image_resize(
      img,
      geometry = paste0(
        panel_width - PANEL_CONTENT_INSET_PX, "x",
        panel_height - PANEL_CONTENT_VERTICAL_INSET_PX, ">"
      )
    )

    img <- magick::image_extent(
      img,
      geometry = paste0(panel_width, "x", panel_height),
      gravity = "center",
      color = FIG_BACKGROUND
    )

    img <- magick::image_border(
      img,
      color = "black",
      geometry = paste0(PANEL_BORDER_RASTER_PX, "x", PANEL_BORDER_RASTER_PX)
    )

    add_panel_label_img(
      img = img,
      label = label,
      size = label_size
    )
  }
}

if (!exists("combine_panels_grid", mode = "function", inherits = FALSE)) {
  combine_panels_grid <- function(panel_files,
                                  panel_labels,
                                  output_file,
                                  ncol = 2,
                                  panel_width = 3600,
                                  panel_height = 2700,
                                  label_size = 130) {
    if (length(panel_files) != length(panel_labels)) {
      stop("panel_files and panel_labels must have the same length.")
    }

    missing <- panel_files[!file.exists(panel_files)]

    if (length(missing) > 0) {
      stop("Missing panel files:
", paste(missing, collapse = "
"))
    }

    imgs <- Map(
      function(path, lab) {
        make_clean_panel_img(
          path = path,
          label = lab,
          panel_width = panel_width,
          panel_height = panel_height,
          label_size = label_size
        )
      },
      panel_files,
      panel_labels
    )

    rows <- list()
    row_index <- 1

    for (i in seq(1, length(imgs), by = ncol)) {
      row_imgs <- imgs[i:min(i + ncol - 1, length(imgs))]

      if (length(row_imgs) < ncol) {
        blank <- magick::image_blank(
          panel_width,
          panel_height,
          color = "white"
        )
        row_imgs <- c(row_imgs, rep(list(blank), ncol - length(row_imgs)))
      }

      rows[[row_index]] <- magick::image_append(
        do.call(c, row_imgs),
        stack = FALSE
      )

      row_index <- row_index + 1
    }

    combined <- magick::image_append(do.call(c, rows), stack = TRUE)
    combined <- magick::image_background(combined, "white", flatten = TRUE)

    magick::image_write(
      combined,
      path = output_file,
      format = "png",
      density = paste0(FIG_DPI, "x", FIG_DPI)
    )

    output_file
  }
}

step05_dir <- make_stage_dir(5, "Triple_Overlap_Figure5A")

bulk_up_file <- file.path(step03_dir, "DEG_up_logFC2_adjP0.01.txt")
wgcna_file   <- file.path(step04_dir, "WGCNA_feature_genes.csv")

sc_marker_file_candidates <- c(
  file.path(scrna_figdir, "High_CNV_Malignant_marker_genes_after_lowCNV_removed_logFC1.txt"),
  file.path(scrna_figdir, "Final_High_CNV_Malignant_Gene_Names.txt"),
  file.path(scrna_figdir, "Final_High_CNV_Malignant_Gene_Names.csv"),
  file.path(scrna_figdir, "High_CNV_Malignant_FindAllMarkers_after_lowCNV_removed_logFC1.csv"),
  file.path(scrna_figdir, "Final_High_CNV_Malignant_Genes.csv"),
  file.path(scrna_figdir, "High_CNV_Malignant_specific_markers_FINAL.csv")
)

safe_read_gene_list <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(character(0))
  }

  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    dt <- data.table::fread(path)
    gene_col <- intersect(
      c("Gene", "Genes", "gene", "Gene.symbol", "Symbol", "symbol"),
      colnames(dt)
    )
    if (length(gene_col) < 1) {
      return(character(0))
    }
    x <- dt[[gene_col[1]]]
  } else {
    x <- readLines(path, warn = FALSE)
  }

  clean_gene_vector(x)
}

if (!file.exists(bulk_up_file)) stop("Bulk upregulated DEG file not found: ", bulk_up_file)
if (!file.exists(wgcna_file)) stop("WGCNA feature gene file not found: ", wgcna_file)

Bulk_UP <- clean_gene_vector(readLines(bulk_up_file, warn = FALSE))

wgcna_df <- as.data.frame(data.table::fread(wgcna_file))
if (!"Gene" %in% colnames(wgcna_df)) {
  stop("WGCNA_feature_genes.csv must contain a column named Gene.")
}

WGCNA_Strong_Tumor_Module_Genes <- clean_gene_vector(wgcna_df$Gene)

sc_marker_file <- sc_marker_file_candidates[file.exists(sc_marker_file_candidates)][1]
if (is.na(sc_marker_file) || !file.exists(sc_marker_file)) {
  stop(
    "No precomputed scRNA marker file found. This script does not run single-cell analysis.\n",
    "Check these paths:\n",
    paste(sc_marker_file_candidates, collapse = "\n")
  )
}

SC_HighCNV_Malignant_Specific_Genes <- safe_read_gene_list(sc_marker_file)

if (length(Bulk_UP) < 1) stop("Bulk_UP gene list is empty.")
if (length(WGCNA_Strong_Tumor_Module_Genes) < 1) stop("WGCNA gene list is empty.")
if (length(SC_HighCNV_Malignant_Specific_Genes) < 1) stop("scRNA marker gene list is empty.")

Triple_Overlap_Genes <- Reduce(
  intersect,
  list(
    Bulk_UP,
    WGCNA_Strong_Tumor_Module_Genes,
    SC_HighCNV_Malignant_Specific_Genes
  )
)

Triple_Overlap_Genes <- sort(unique(Triple_Overlap_Genes))

write.csv(
  data.frame(Gene = Triple_Overlap_Genes),
  file.path(step05_dir, "TripleOverlap_BulkUP_WGCNA_HighCNV_Malignant.csv"),
  row.names = FALSE
)

writeLines(Triple_Overlap_Genes, file.path(step05_dir, "TripleOverlap_genes.txt"))

if (length(Triple_Overlap_Genes) < 2) {
  stop("Too few triple-overlap genes for LASSO. Check DEG/WGCNA/scRNA thresholds.")
}

cat("\nSTEP 05 finished: triple overlap generated.\n")
cat("Triple-overlap genes:", length(Triple_Overlap_Genes), "\n")

# STEP 06. LASSO FEATURE SELECTION AND FIGURE 5B

if (!exists("results_path", inherits = FALSE)) {
  results_path <- "D:/LSCC/Results_LSCC/ML"
}

if (!exists("FIG_DPI", inherits = FALSE)) {
  FIG_DPI <- 600
}

dir.create(results_path, recursive = TRUE, showWarnings = FALSE)

if (!exists("make_stage_dir", mode = "function", inherits = FALSE)) {
  make_stage_dir <- function(step_number, step_name) {
    dir.create(results_path, recursive = TRUE, showWarnings = FALSE)
    results_path
  }
}

if (!exists("clean_gene_vector", mode = "function", inherits = FALSE)) {
  clean_gene_vector <- function(x) {
    x <- unique(trimws(as.character(x)))
    x[!is.na(x) & x != ""]
  }
}

if (!exists("to_numeric_df", mode = "function", inherits = FALSE)) {
  to_numeric_df <- function(df) {
    as.data.frame(
      lapply(df, function(x) as.numeric(as.character(x))),
      check.names = FALSE
    )
  }
}

if (!exists("coerce_group_12", mode = "function", inherits = FALSE)) {
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
}

if (!exists("make_group_factor", mode = "function", inherits = FALSE)) {
  make_group_factor <- function(g) {
    g12 <- coerce_group_12(g)
    factor(g12, levels = c(1, 2), labels = c("Normal", "Tumor"))
  }
}

if (!exists("train_log2_file", inherits = FALSE)) {
  train_log2_file <- file.path(results_path, "train_log2CPM.csv")
}

if (!exists("valid_log2_file", inherits = FALSE)) {
  valid_log2_file <- file.path(results_path, "valid_log2CPM.csv")
}

if (!exists("external_log2_file", inherits = FALSE)) {
  external_log2_file <- file.path(
    results_path,
    "external_GSE130605_log2CPM.csv"
  )
}

required_step06_files <- c(
  train_log2_file,
  valid_log2_file,
  external_log2_file
)

missing_step06_files <- required_step06_files[!file.exists(required_step06_files)]

if (length(missing_step06_files) > 0) {
  stop(
    "Step 06 requires the log2-CPM files generated in Step 03. Missing file(s):
",
    paste(missing_step06_files, collapse = "
"),
    "

Run the complete script from Step 01, or run Step 03 before Step 06."
  )
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required for Step 06.")
}

if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("Package 'glmnet' is required for Step 06.")
}

if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required for Step 06.")
}

if (!("package:dplyr" %in% search())) {
  suppressPackageStartupMessages(library(dplyr))
}

step06_dir <- make_stage_dir(6, "LASSO_Figure5B")

LASSO_NFOLDS <- 10

triple_file <- file.path(step05_dir, "TripleOverlap_genes.txt")
Triple_Overlap_Genes <- clean_gene_vector(readLines(triple_file, warn = FALSE))

if (length(Triple_Overlap_Genes) < 2) {
  stop("Too few triple-overlap genes for LASSO.")
}

train_df <- as.data.frame(data.table::fread(train_log2_file))
valid_df <- as.data.frame(data.table::fread(valid_log2_file))
test_df  <- as.data.frame(data.table::fread(external_log2_file))

common_model_genes <- Reduce(
  intersect,
  list(
    Triple_Overlap_Genes,
    setdiff(colnames(train_df), c("Sample", "group", "batch")),
    setdiff(colnames(valid_df), c("Sample", "group", "batch")),
    setdiff(colnames(test_df),  c("Sample", "group", "batch"))
  )
)

common_model_genes <- sort(unique(common_model_genes))

if (length(common_model_genes) < 2) {
  stop("Too few triple-overlap genes available in train/valid/test matrices for LASSO.")
}

write.csv(
  data.frame(Gene = common_model_genes),
  file.path(step06_dir, "Genes_entering_LASSO_from_triple_overlap.csv"),
  row.names = FALSE
)
writeLines(common_model_genes, file.path(step06_dir, "Genes_entering_LASSO_from_triple_overlap.txt"))

trainX <- as.matrix(to_numeric_df(train_df[, common_model_genes, drop = FALSE]))
validX <- as.matrix(to_numeric_df(valid_df[, common_model_genes, drop = FALSE]))
testX  <- as.matrix(to_numeric_df(test_df[,  common_model_genes, drop = FALSE]))

rownames(trainX) <- train_df$Sample
rownames(validX) <- valid_df$Sample
rownames(testX)  <- test_df$Sample

trainY <- make_group_factor(train_df$group)
validY <- make_group_factor(valid_df$group)
testY  <- make_group_factor(test_df$group)

x <- trainX
y <- ifelse(trainY == "Tumor", 1, 0)

nfolds_use <- min(LASSO_NFOLDS, length(y))
if (nfolds_use < 3) stop("Too few samples for LASSO cross-validation.")

set.seed(123)

cv_lasso <- glmnet::cv.glmnet(
  x = x,
  y = y,
  alpha = 1,
  family = "binomial",
  type.measure = "deviance",
  nfolds = nfolds_use
)

fit_lasso <- glmnet::glmnet(
  x = x,
  y = y,
  alpha = 1,
  family = "binomial",
  lambda = cv_lasso$lambda.min
)

beta <- as.vector(fit_lasso$beta[, 1])
names(beta) <- rownames(fit_lasso$beta)

lasso_coef_table <- data.frame(
  Gene = names(beta),
  Coefficient = beta,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(Coefficient != 0) %>%
  dplyr::arrange(dplyr::desc(abs(Coefficient)))

if (nrow(lasso_coef_table) < 1) {
  stop("No genes selected by LASSO.")
}

write.csv(
  lasso_coef_table,
  file.path(step06_dir, "Final_LASSO_coefficients_lambda_min.csv"),
  row.names = FALSE
)

writeLines(lasso_coef_table$Gene, file.path(step06_dir, "Final_LASSO_genes.txt"))

lasso_genes <- lasso_coef_table$Gene

lasso_summary <- data.frame(
  Input_gene_count = length(common_model_genes),
  Selected_gene_count = length(lasso_genes),
  Lambda_min = cv_lasso$lambda.min,
  Lambda_1se = cv_lasso$lambda.1se,
  Nfolds = nfolds_use
)

write.csv(lasso_summary, file.path(step06_dir, "LASSO_summary.csv"), row.names = FALSE)

cat("\nSTEP 06 finished: LASSO feature selection generated.\n")
print(lasso_coef_table)

# STEP 07. INTERNAL ROC ANALYSIS AND FIGURE 5C

step07_dir <- make_stage_dir(7, "Internal_ROC_Figure5C")

ROC_AUC_MIN <- 0.95

calc_gene_auc <- function(expr, group_factor) {
  group_factor <- factor(group_factor, levels = c("Normal", "Tumor"))
  expr <- as.numeric(expr)

  ok <- is.finite(expr) & !is.na(group_factor)

  if (sum(ok) < 3 || length(unique(group_factor[ok])) < 2) {
    return(NA_real_)
  }

  roc_obj <- pROC::roc(
    response = group_factor[ok],
    predictor = expr[ok],
    levels = c("Normal", "Tumor"),
    direction = "<",
    quiet = TRUE
  )

  as.numeric(roc_obj$auc)
}

lasso_genes <- clean_gene_vector(readLines(file.path(step06_dir, "Final_LASSO_genes.txt"), warn = FALSE))

common_model_genes <- clean_gene_vector(readLines(file.path(step06_dir, "Genes_entering_LASSO_from_triple_overlap.txt"), warn = FALSE))
common_model_genes <- Reduce(
  intersect,
  list(
    common_model_genes,
    setdiff(colnames(train_df), c("Sample", "group", "batch")),
    setdiff(colnames(valid_df), c("Sample", "group", "batch")),
    setdiff(colnames(test_df),  c("Sample", "group", "batch"))
  )
)
common_model_genes <- sort(unique(common_model_genes))

trainX <- as.matrix(to_numeric_df(train_df[, common_model_genes, drop = FALSE]))
validX <- as.matrix(to_numeric_df(valid_df[, common_model_genes, drop = FALSE]))
testX  <- as.matrix(to_numeric_df(test_df[,  common_model_genes, drop = FALSE]))

rownames(trainX) <- train_df$Sample
rownames(validX) <- valid_df$Sample
rownames(testX)  <- test_df$Sample

trainY <- make_group_factor(train_df$group)
validY <- make_group_factor(valid_df$group)
testY  <- make_group_factor(test_df$group)

lasso_genes <- Reduce(intersect, list(lasso_genes, colnames(trainX), colnames(validX)))

if (length(lasso_genes) < 1) {
  stop("None of the LASSO genes are present in trainX and validX.")
}

roc_table <- lapply(lasso_genes, function(g) {
  data.frame(
    Gene = g,
    Train_AUC = calc_gene_auc(trainX[, g], trainY),
    Valid_AUC = calc_gene_auc(validX[, g], validY),
    stringsAsFactors = FALSE
  )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    Pass_ROC_Cutoff = is.finite(Train_AUC) &
      is.finite(Valid_AUC) &
      Train_AUC >= ROC_AUC_MIN &
      Valid_AUC >= ROC_AUC_MIN,
    Minimum_AUC = pmin(Train_AUC, Valid_AUC, na.rm = TRUE)
  ) %>%
  dplyr::arrange(dplyr::desc(Pass_ROC_Cutoff), dplyr::desc(Minimum_AUC))

write.csv(roc_table, file.path(step07_dir, "ROC_AUC_table_LASSO_genes.csv"), row.names = FALSE)

final_biomarkers <- roc_table %>%
  dplyr::filter(Pass_ROC_Cutoff) %>%
  dplyr::pull(Gene)

if (length(final_biomarkers) < 1) {
  write.csv(
    roc_table,
    file.path(step07_dir, "ROC_AUC_table_LASSO_genes_NO_GENE_PASSED_0.95.csv"),
    row.names = FALSE
  )
  stop(
    "No LASSO gene passed the strict ROC cutoff. Required: Train_AUC >= ",
    ROC_AUC_MIN,
    " AND Valid_AUC >= ",
    ROC_AUC_MIN
  )
}

writeLines(final_biomarkers, file.path(step07_dir, "FINAL_BIOMARKERS.txt"))

write.csv(
  roc_table %>% dplyr::filter(Gene %in% final_biomarkers),
  file.path(step07_dir, "FINAL_BIOMARKERS_AUC_table.csv"),
  row.names = FALSE
)

cat("\nSTEP 07 finished: internal ROC statistics generated.\n")
print(roc_table)

# STEP 08. EXPRESSION VALIDATION BOXPLOTS AND MANUSCRIPT FIGURE 5

step08_dir <- make_stage_dir(8, "Expression_Boxplots_Figure5D_and_Combined_Figure5")

NORMAL_FILL <- "#E7E5F2"
NORMAL_COL  <- "#7A74A8"
TUMOR_FILL  <- "#FFD6E3"
TUMOR_COL   <- "#E75480"

format_p_value <- function(p) {
  if (!is.finite(p)) return("p = NA")
  paste0("p = ", formatC(p, format = "e", digits = 1))
}

build_long_df <- function(df, dataset_name, genes) {
  genes_present <- intersect(
    genes,
    setdiff(colnames(df), c("Sample", "group", "batch"))
  )

  if (length(genes_present) < 1) {
    return(data.frame())
  }

  tmp <- df[, c("Sample", "group", genes_present), drop = FALSE]
  tmp$Group <- make_group_factor(tmp$group)
  tmp$Dataset <- dataset_name
  tmp$group <- NULL

  long_df <- tmp %>%
    tidyr::pivot_longer(
      cols = all_of(genes_present),
      names_to = "Gene",
      values_to = "Expression"
    ) %>%
    dplyr::mutate(
      Expression = as.numeric(Expression),
      Dataset = factor(Dataset, levels = c("Training", "Validation")),
      Group = factor(Group, levels = c("Normal", "Tumor"))
    ) %>%
    dplyr::filter(is.finite(Expression))

  as.data.frame(long_df)
}

final_biomarkers <- clean_gene_vector(readLines(file.path(step07_dir, "FINAL_BIOMARKERS.txt"), warn = FALSE))

writeLines(final_biomarkers, file.path(step08_dir, "Final_biomarkers_used_for_expression_boxplot.txt"))

expr_long <- dplyr::bind_rows(
  build_long_df(train_df, "Training", final_biomarkers),
  build_long_df(valid_df, "Validation", final_biomarkers)
)

if (nrow(expr_long) < 1) {
  stop("Expression table is empty. Final biomarker(s) were not found in train/validation matrices.")
}

expr_long$Gene <- factor(expr_long$Gene, levels = final_biomarkers)

write.csv(expr_long, file.path(step08_dir, "Final_biomarker_expression_long_table.csv"), row.names = FALSE)

pval_df <- expr_long %>%
  dplyr::group_by(Gene, Dataset) %>%
  dplyr::summarise(
    Wilcoxon_p = tryCatch(wilcox.test(Expression ~ Group)$p.value, error = function(e) NA_real_),
    y_max = max(Expression, na.rm = TRUE),
    y_min = min(Expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    y_range = ifelse(is.finite(y_max - y_min) & (y_max - y_min) > 0, y_max - y_min, 1),
    y_line = y_max + 0.10 * y_range,
    y_text = y_max + 0.17 * y_range,
    p_label = vapply(Wilcoxon_p, format_p_value, character(1))
  )

expr_summary <- expr_long %>%
  dplyr::group_by(Gene, Dataset, Group) %>%
  dplyr::summarise(
    Median = median(Expression, na.rm = TRUE),
    Mean = mean(Expression, na.rm = TRUE),
    SD = sd(Expression, na.rm = TRUE),
    N = dplyr::n(),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(names_from = Group, values_from = c(Median, Mean, SD, N)) %>%
  dplyr::left_join(
    pval_df %>% dplyr::select(Gene, Dataset, Wilcoxon_p),
    by = c("Gene", "Dataset")
  )

write.csv(expr_summary, file.path(step08_dir, "Final_biomarker_expression_validation_summary.csv"), row.names = FALSE)

cat("\nSTEP 08 finished: expression-validation tables generated.\n")
packageVersion("glmnet")
packageVersion("pROC")

options(stringsAsFactors = FALSE, scipen = 100, width = 140)
set.seed(123)

results_path <- "D:/LSCC/Results_LSCC/ML"
FIGURE_DIR <- results_path

scrna_figdir <- "D:/LSCC/ScRNAseq_Results/GSE206332/Results/figures"

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

if (.Platform$OS.type == "windows") {
  try(
    grDevices::windowsFonts(
      Arial = grDevices::windowsFont("Arial")
    ),
    silent = TRUE
  )
}

FONT_FAMILY <- "Arial"
FIG_DPI <- 600
FIG_BACKGROUND <- "white"

FIG_SINGLE_W <- 3.50
FIG_SINGLE_H <- 4.20

FIG_DOUBLE_W <- 7.20
FIG_DOUBLE_H <- 5.40

FIGURE2_COMBINED_W <- 14.00
FIGURE2_COMBINED_H <- 6.40

BASE_TEXT_PT <- 9.5

AXIS_TEXT_PT <- 8.5
AXIS_TITLE_PT <- 10

LEGEND_TEXT_PT <- 8.5
LEGEND_TITLE_PT <- 9

PLOT_TITLE_PT <- 10
PANEL_TAG_PT <- 14

PANEL_BORDER_LWD <- 0.45
AXIS_LWD <- 0.40
GEOM_LWD <- 0.55
GRID_LWD <- 0.30

HEATMAP_FONT_PT <- 8
HEATMAP_ROW_FONT_PT <- 7
HEATMAP_COL_FONT_PT <- 7
HEATMAP_LEGEND_FONT_PT <- 8

COL_NORMAL <- "#7470B2"

COL_TUMOR <- "#D81B60"

COL_MYBL2_LOW <- COL_NORMAL
COL_MYBL2_HIGH <- COL_TUMOR

COL_UP <- COL_TUMOR
COL_DOWN <- COL_NORMAL
COL_NS <- "#B8B8B8"

CELLCHAT_PALETTE <- grDevices::colorRampPalette(
  c("#E8F6F4", "#B8E1DA", "#73C6BE", "#2D9D8C", "#006E63")
)(100)

COL_WGCNA <- "#2D9D8C"

NORMAL_FILL <- "#DDD9EE"
TUMOR_FILL <- "#F5C6D8"

MANUSCRIPT_MARGIN <- ggplot2::margin(
  t = 5.5,
  r = 6,
  b = 5.5,
  l = 6,
  unit = "pt"
)

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
        size = PLOT_TITLE_PT,
        hjust = 0.5,
        vjust = 0.5,
        color = "black"
      ),
      plot.title.position = "panel",

      plot.tag = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PANEL_TAG_PT,
        hjust = 0,
        vjust = 1,
        color = "black"
      ),
      plot.tag.position = PANEL_TAG_POSITION,

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
      axis.ticks.length = grid::unit(TICK_LENGTH_PT, "pt"),

      legend.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = LEGEND_TITLE_PT,
        color = "black"
      ),
      legend.text = ggplot2::element_text(
        family = FONT_FAMILY,
        size = LEGEND_TEXT_PT,
        color = "black"
      ),
      legend.position = legend_position,
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

      strip.background = ggplot2::element_rect(
        fill = "grey95",
        color = "black",
        linewidth = PANEL_BORDER_LWD
      ),
      strip.text = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = AXIS_TEXT_PT,
        color = "black"
      ),

      panel.spacing = grid::unit(PANEL_SPACING_PT, "pt"),
      plot.margin = MANUSCRIPT_MARGIN
    )
}

save_plot_all_formats <- function(plot_obj,
                                  filename_stem,
                                  dir_path,
                                  width = FIG_DOUBLE_W,
                                  height = FIG_DOUBLE_H) {

  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  png_file <- file.path(dir_path, paste0(filename_stem, ".png"))
  tiff_file <- file.path(dir_path, paste0(filename_stem, ".tiff"))
  pdf_file <- file.path(dir_path, paste0(filename_stem, ".pdf"))

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
    bg = FIG_BACKGROUND,
    compression = "lzw",
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
    list(
      PNG = png_file,
      TIFF = tiff_file,
      PDF = pdf_file
    )
  )
}

save_base_plot_all_formats <- function(draw_fun,
                                       filename_stem,
                                       dir_path,
                                       width = FIG_DOUBLE_W,
                                       height = FIG_DOUBLE_H) {

  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  png_file <- file.path(dir_path, paste0(filename_stem, ".png"))
  tiff_file <- file.path(dir_path, paste0(filename_stem, ".tiff"))
  pdf_file <- file.path(dir_path, paste0(filename_stem, ".pdf"))

  grDevices::png(
    filename = png_file,
    width = width,
    height = height,
    units = "in",
    res = FIG_DPI,
    bg = FIG_BACKGROUND
  )
  draw_fun()
  grDevices::dev.off()

  grDevices::tiff(
    filename = tiff_file,
    width = width,
    height = height,
    units = "in",
    res = FIG_DPI,
    compression = "lzw",
    bg = FIG_BACKGROUND
  )
  draw_fun()
  grDevices::dev.off()

  if (capabilities("cairo")) {
    grDevices::cairo_pdf(
      filename = pdf_file,
      width = width,
      height = height,
      family = FONT_FAMILY
    )
  } else {
    grDevices::pdf(
      file = pdf_file,
      width = width,
      height = height,
      family = FONT_FAMILY
    )
  }
  draw_fun()
  grDevices::dev.off()

  invisible(
    list(
      PNG = png_file,
      TIFF = tiff_file,
      PDF = pdf_file
    )
  )
}

save_magick_all_formats <- function(image_object,
                                    filename_stem,
                                    dir_path) {

  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  png_file <- file.path(dir_path, paste0(filename_stem, ".png"))
  tiff_file <- file.path(dir_path, paste0(filename_stem, ".tiff"))
  pdf_file <- file.path(dir_path, paste0(filename_stem, ".pdf"))

  magick::image_write(
    image_object,
    path = png_file,
    format = "png",
    density = paste0(FIG_DPI, "x", FIG_DPI)
  )

  magick::image_write(
    image_object,
    path = tiff_file,
    format = "tiff",
    density = paste0(FIG_DPI, "x", FIG_DPI)
  )

  try(
    magick::image_write(
      image_object,
      path = pdf_file,
      format = "pdf",
      density = paste0(FIG_DPI, "x", FIG_DPI)
    ),
    silent = TRUE
  )

  invisible(
    list(
      PNG = png_file,
      TIFF = tiff_file,
      PDF = pdf_file
    )
  )
}

required_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggvenn",
  "glmnet",
  "pROC",
  "cowplot",
  "magick"
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
    "These required packages are missing:\n",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggvenn)
  library(glmnet)
  library(pROC)
  library(cowplot)
  library(magick)
})

clean_gene_vector <- function(x) {
  x <- unique(trimws(as.character(x)))
  x[!is.na(x) & x != ""]
}

to_numeric_df <- function(df) {
  as.data.frame(
    lapply(df, function(x) as.numeric(as.character(x))),
    check.names = FALSE
  )
}

coerce_group_12 <- function(g) {
  g0 <- as.character(g)

  if (all(g0 %in% c("1", "2"))) {
    return(as.integer(g0))
  }

  g_low <- tolower(g0)

  if (all(g_low %in% c("normal", "tumor", "non", "lscc", "cancer", "margin"))) {
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
      "group must have exactly 2 classes. Found: ",
      paste(u, collapse = ", ")
    )
  }

  map <- setNames(c(1L, 2L), u)
  as.integer(map[g0])
}

make_group_factor <- function(g) {
  g12 <- coerce_group_12(g)
  factor(g12, levels = c(1, 2), labels = c("Normal", "Tumor"))
}

first_existing <- function(paths) {
  paths <- unique(paths)
  paths <- paths[file.exists(paths)]

  if (length(paths) < 1) {
    return(NA_character_)
  }

  paths[1]
}

safe_read_gene_list <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(character(0))
  }

  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    dat <- as.data.frame(
      data.table::fread(path),
      check.names = FALSE
    )

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

save_fig <- function(plot_obj,
                     filename,
                     dir_path,
                     width = FIG_DOUBLE_W,
                     height = FIG_DOUBLE_H) {

  save_plot_all_formats(
    plot_obj = plot_obj,
    filename_stem = tools::file_path_sans_ext(filename),
    dir_path = dir_path,
    width = width,
    height = height
  )
}

add_panel_label_img <- function(img,
                                label,
                                size = PANEL_TAG_RASTER_SIZE,
                                offset_x = PANEL_TAG_RASTER_OFFSET_X,
                                offset_y = PANEL_TAG_RASTER_OFFSET_Y) {
  magick::image_annotate(
    img,
    text = label,
    size = size,
    font = FONT_FAMILY,
    gravity = "northwest",
    location = paste0("+", offset_x, "+", offset_y),
    weight = 700,
    color = "black"
  )
}

make_clean_panel_img <- function(path,
                                 label,
                                 panel_width = 3600,
                                 panel_height = 2850,
                                 label_size = PANEL_TAG_RASTER_SIZE) {
  img <- magick::image_read(path)
  img <- magick::image_background(img, FIG_BACKGROUND, flatten = TRUE)
  img <- magick::image_trim(img)

  img <- magick::image_resize(
    img,
    geometry = paste0(
      panel_width - PANEL_CONTENT_INSET_PX, "x",
      panel_height - PANEL_CONTENT_VERTICAL_INSET_PX, ">"
    )
  )

  img <- magick::image_extent(
    img,
    geometry = paste0(panel_width, "x", panel_height),
    gravity = "center",
    color = FIG_BACKGROUND
  )

  img <- magick::image_border(
    img,
    color = "black",
    geometry = paste0(PANEL_BORDER_RASTER_PX, "x", PANEL_BORDER_RASTER_PX)
  )

  add_panel_label_img(
    img = img,
    label = label,
    size = label_size
  )
}

bulk_up_file <- file.path(
  results_path,
  "DEG_up_logFC2_adjP0.01.txt"
)

wgcna_file <- file.path(
  results_path,
  "WGCNA_feature_genes.csv"
)

triple_file <- file.path(
  results_path,
  "TripleOverlap_genes.txt"
)

lasso_genes_file <- file.path(
  results_path,
  "Final_LASSO_genes.txt"
)

final_biomarkers_file <- file.path(
  results_path,
  "FINAL_BIOMARKERS.txt"
)

train_log2_file <- file.path(
  results_path,
  "train_log2CPM.csv"
)

valid_log2_file <- file.path(
  results_path,
  "valid_log2CPM.csv"
)

external_log2_file <- file.path(
  results_path,
  "external_GSE130605_log2CPM.csv"
)

required_files <- c(
  bulk_up_file,
  wgcna_file,
  triple_file,
  lasso_genes_file,
  final_biomarkers_file,
  train_log2_file,
  valid_log2_file,
  external_log2_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "These existing pipeline output files are required to regenerate Figure 5:\n",
    paste(missing_files, collapse = "\n")
  )
}

sc_marker_file <- first_existing(
  c(
    file.path(
      scrna_figdir,
      "High_CNV_Malignant_marker_genes_after_lowCNV_removed_logFC1.txt"
    ),
    file.path(
      scrna_figdir,
      "Final_High_CNV_Malignant_Gene_Names.txt"
    ),
    file.path(
      scrna_figdir,
      "Final_High_CNV_Malignant_Gene_Names.csv"
    ),
    file.path(
      scrna_figdir,
      "High_CNV_Malignant_FindAllMarkers_after_lowCNV_removed_logFC1.csv"
    ),
    file.path(
      scrna_figdir,
      "Final_High_CNV_Malignant_Genes.csv"
    ),
    file.path(
      scrna_figdir,
      "High_CNV_Malignant_specific_markers_FINAL.csv"
    )
  )
)

if (is.na(sc_marker_file)) {
  stop(
    "No precomputed high-CNV malignant scRNA marker file was found in:\n",
    scrna_figdir
  )
}

Bulk_UP <- clean_gene_vector(
  readLines(bulk_up_file, warn = FALSE)
)

wgcna_df <- as.data.frame(
  data.table::fread(wgcna_file),
  check.names = FALSE
)

if (!("Gene" %in% colnames(wgcna_df))) {
  stop("WGCNA_feature_genes.csv must contain a column named 'Gene'.")
}

WGCNA_Strong_Tumor_Module_Genes <- clean_gene_vector(
  wgcna_df$Gene
)

SC_HighCNV_Malignant_Specific_Genes <- safe_read_gene_list(
  sc_marker_file
)

Triple_Overlap_Genes <- clean_gene_vector(
  readLines(triple_file, warn = FALSE)
)

lasso_genes <- clean_gene_vector(
  readLines(lasso_genes_file, warn = FALSE)
)

final_biomarkers <- clean_gene_vector(
  readLines(final_biomarkers_file, warn = FALSE)
)

train_df <- as.data.frame(
  data.table::fread(train_log2_file),
  check.names = FALSE
)

valid_df <- as.data.frame(
  data.table::fread(valid_log2_file),
  check.names = FALSE
)

test_df <- as.data.frame(
  data.table::fread(external_log2_file),
  check.names = FALSE
)

# FIGURE 5A. TRIPLE-OVERLAP VENN DIAGRAM

sets3 <- list(
  "Upregulated\nDEGs" = Bulk_UP,
  "Tumor-related\nWGCNA genes" = WGCNA_Strong_Tumor_Module_Genes,
  "High-CNV malignant\nscRNA genes" = SC_HighCNV_Malignant_Specific_Genes
)

p_fig5a <- ggvenn::ggvenn(
  sets3,
  fill_color = c(COL_UP, COL_WGCNA, COL_NORMAL),
  fill_alpha = 0.76,
  stroke_size = PANEL_BORDER_LWD,
  text_size = 4.6,
  set_name_size = 3.25,
  show_percentage = FALSE
) +
  ggplot2::theme_void(
    base_family = FONT_FAMILY,
    base_size = BASE_TEXT_PT
  ) +
  ggplot2::theme(
    text = ggplot2::element_text(
      family = FONT_FAMILY,
      color = "black"
    ),
    plot.background = ggplot2::element_rect(
      fill = FIG_BACKGROUND,
      color = NA
    ),
    panel.background = ggplot2::element_rect(
      fill = FIG_BACKGROUND,
      color = NA
    ),
    panel.border = ggplot2::element_rect(
      color = "black",
      fill = NA,
      linewidth = PANEL_BORDER_LWD
    ),
    plot.margin = MANUSCRIPT_MARGIN
  )

save_plot_all_formats(
  plot_obj = p_fig5a,
  filename_stem = "Figure_05A_Triple_overlap_Venn",
  dir_path = FIGURE_DIR,
  width = FIG_DOUBLE_W,
  height = FIG_DOUBLE_H
)

fig5a_file <- file.path(
  FIGURE_DIR,
  "Figure_05A_Triple_overlap_Venn.png"
)

# FIGURE 5B. LASSO CROSS-VALIDATION CURVE

common_model_genes <- Reduce(
  intersect,
  list(
    Triple_Overlap_Genes,
    setdiff(colnames(train_df), c("Sample", "group", "batch")),
    setdiff(colnames(valid_df), c("Sample", "group", "batch")),
    setdiff(colnames(test_df), c("Sample", "group", "batch"))
  )
)

common_model_genes <- sort(unique(common_model_genes))

if (length(common_model_genes) < 2) {
  stop(
    "Too few triple-overlap genes are available for reconstructing Figure 5B."
  )
}

trainX <- as.matrix(
  to_numeric_df(
    train_df[, common_model_genes, drop = FALSE]
  )
)

rownames(trainX) <- train_df$Sample

trainY <- make_group_factor(train_df$group)

set.seed(123)

cv_lasso <- glmnet::cv.glmnet(
  x = trainX,
  y = ifelse(trainY == "Tumor", 1, 0),
  alpha = 1,
  family = "binomial",
  type.measure = "deviance",
  nfolds = min(10, length(trainY))
)

draw_fig5b <- function() {

  par(
    family = FONT_FAMILY,
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black",
    lwd = AXIS_LWD,
    cex.axis = BASE_AXIS_CEX,
    cex.lab = BASE_AXIS_TITLE_CEX,
    cex.main = BASE_AXIS_TITLE_CEX,
    las = 1,
    mar = c(4.2, 4.5, 1.1, 1.0)
  )

  plot(
    cv_lasso,
    xlab = "log(λ)",
    ylab = "Binomial deviance",
    main = ""
  )

  plot_usr <- par("usr")

  lambda_1se_x <- -log(cv_lasso$lambda.1se)
  lambda_min_x <- -log(cv_lasso$lambda.min)

  label_y_low <- plot_usr[3] + 0.075 * diff(plot_usr[3:4])
  label_y_high <- plot_usr[4] - 0.060 * diff(plot_usr[3:4])

  text(
    x = lambda_1se_x,
    y = label_y_high,
    labels = "lambda.1se",
    pos = 2,
    offset = 0.30,
    cex = BASE_ANNOTATION_CEX,
    col = "black"
  )

  text(
    x = lambda_min_x,
    y = label_y_low,
    labels = "lambda.min",
    pos = 4,
    offset = 0.30,
    cex = BASE_ANNOTATION_CEX,
    col = "black"
  )
}

save_base_plot_all_formats(
  draw_fun = draw_fig5b,
  filename_stem = "Figure_05B_LASSO_cv_curve",
  dir_path = FIGURE_DIR,
  width = FIG_DOUBLE_W,
  height = FIG_DOUBLE_H
)

fig5b_file <- file.path(
  FIGURE_DIR,
  "Figure_05B_LASSO_cv_curve.png"
)

# FIGURE 5C. INTERNAL ROC CURVES

calc_roc_object <- function(expression,
                            group_factor) {

  pROC::roc(
    response = group_factor,
    predictor = as.numeric(expression),
    levels = c("Normal", "Tumor"),
    direction = "<",
    quiet = TRUE
  )
}

make_internal_roc_plot <- function(gene) {

  roc_train <- calc_roc_object(
    expression = train_df[[gene]],
    group_factor = make_group_factor(train_df$group)
  )

  roc_valid <- calc_roc_object(
    expression = valid_df[[gene]],
    group_factor = make_group_factor(valid_df$group)
  )

  roc_list <- list(
    Training = roc_train,
    Validation = roc_valid
  )

  pROC::ggroc(
    roc_list,
    legacy.axes = TRUE,
    linewidth = GEOM_LWD * 2
  ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey55",
      linewidth = GRID_LWD + 0.10
    ) +
    ggplot2::scale_color_manual(
      values = c(
        Training = COL_TUMOR,
        Validation = COL_NORMAL
      ),
      labels = c(
        Training = paste0(
          "Training AUC = ",
          sprintf("%.2f", as.numeric(roc_train$auc))
        ),
        Validation = paste0(
          "Validation AUC = ",
          sprintf("%.2f", as.numeric(roc_valid$auc))
        )
      )
    ) +
    ggplot2::coord_fixed(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    ggplot2::labs(
      title = gene,
      x = "1 - Specificity",
      y = "Sensitivity",
      color = "Dataset"
    ) +
    theme_manuscript(
      show_grid = FALSE,
      legend_position = c(0.66, 0.18)
    ) +
    ggplot2::theme(
      legend.background = ggplot2::element_rect(
        fill = FIG_BACKGROUND,
        color = "black",
        linewidth = PANEL_BORDER_LWD
      ),
      legend.key.width = grid::unit(0.40, "cm"),
      legend.spacing.y = grid::unit(0.04, "cm")
    )
}

roc_genes <- intersect(
  final_biomarkers,
  Reduce(
    intersect,
    list(
      lasso_genes,
      colnames(train_df),
      colnames(valid_df)
    )
  )
)

if (length(roc_genes) < 1) {
  stop(
    "No final biomarker was found in the training and validation matrices."
  )
}

roc_plots <- lapply(roc_genes, make_internal_roc_plot)
names(roc_plots) <- roc_genes

for (gene in names(roc_plots)) {
  safe_gene_name <- gsub("[^A-Za-z0-9_\\-]", "_", gene)

  save_plot_all_formats(
    plot_obj = roc_plots[[gene]],
    filename_stem = paste0("Figure_05C_ROC_", safe_gene_name),
    dir_path = FIGURE_DIR,
    width = FIG_SINGLE_W,
    height = FIG_SINGLE_H
  )
}

fig5c_combined <- cowplot::plot_grid(
  plotlist = roc_plots,
  ncol = min(3, length(roc_plots)),
  align = "hv"
)

save_plot_all_formats(
  plot_obj = fig5c_combined,
  filename_stem = "Figure_05C_Internal_ROC_final_biomarker",
  dir_path = FIGURE_DIR,
  width = max(FIG_DOUBLE_W, min(3, length(roc_plots)) * FIG_DOUBLE_W),
  height = ceiling(length(roc_plots) / 3) * FIG_DOUBLE_H
)

fig5c_file <- file.path(
  FIGURE_DIR,
  "Figure_05C_Internal_ROC_final_biomarker.png"
)

# FIGURE 5D. EXPRESSION BOXPLOTS IN TRAINING AND VALIDATION COHORTS

format_p_value <- function(p) {
  if (!is.finite(p)) {
    return("P = NA")
  }

  paste0(
    "P = ",
    formatC(p, format = "e", digits = 1)
  )
}

build_long_df <- function(df,
                          dataset_name,
                          genes) {

  genes_present <- intersect(
    genes,
    setdiff(colnames(df), c("Sample", "group", "batch"))
  )

  if (length(genes_present) < 1) {
    return(data.frame())
  }

  tmp <- df[, c("Sample", "group", genes_present), drop = FALSE]

  tmp$Group <- make_group_factor(tmp$group)
  tmp$Dataset <- dataset_name
  tmp$group <- NULL

  long_df <- tmp %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(genes_present),
      names_to = "Gene",
      values_to = "Expression"
    ) %>%
    dplyr::mutate(
      Expression = as.numeric(Expression),
      Dataset = factor(
        Dataset,
        levels = c("Training", "Validation")
      ),
      Group = factor(
        Group,
        levels = c("Normal", "Tumor")
      )
    ) %>%
    dplyr::filter(is.finite(Expression))

  as.data.frame(long_df)
}

expr_long <- dplyr::bind_rows(
  build_long_df(
    df = train_df,
    dataset_name = "Training",
    genes = final_biomarkers
  ),
  build_long_df(
    df = valid_df,
    dataset_name = "Validation",
    genes = final_biomarkers
  )
)

if (nrow(expr_long) < 1) {
  stop(
    "Expression data for the final biomarker(s) were not found in the existing matrices."
  )
}

expr_long$Gene <- factor(
  expr_long$Gene,
  levels = final_biomarkers
)

pval_df <- expr_long %>%
  dplyr::group_by(Gene, Dataset) %>%
  dplyr::summarise(
    Wilcoxon_p = tryCatch(
      wilcox.test(Expression ~ Group)$p.value,
      error = function(e) NA_real_
    ),
    y_max = max(Expression, na.rm = TRUE),
    y_min = min(Expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    y_range = ifelse(
      is.finite(y_max - y_min) & (y_max - y_min) > 0,
      y_max - y_min,
      1
    ),
    y_line = y_max + 0.10 * y_range,
    y_text = y_max + 0.17 * y_range,
    p_label = vapply(
      Wilcoxon_p,
      format_p_value,
      character(1)
    )
  )

make_gene_boxplot <- function(gene_name) {

  df_gene <- expr_long %>%
    dplyr::filter(Gene == gene_name)

  stat_gene <- pval_df %>%
    dplyr::filter(Gene == gene_name)

  ggplot2::ggplot(
    df_gene,
    ggplot2::aes(
      x = Group,
      y = Expression,
      fill = Group,
      color = Group
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.52,
      outlier.shape = NA,
      alpha = 0.45,
      linewidth = GEOM_LWD
    ) +
    ggplot2::geom_jitter(
      width = 0.11,
      height = 0,
      size = 1.45,
      alpha = 0.78
    ) +
    ggplot2::geom_segment(
      data = stat_gene,
      ggplot2::aes(
        x = 1,
        xend = 2,
        y = y_line,
        yend = y_line
      ),
      inherit.aes = FALSE,
      linewidth = GEOM_LWD,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = stat_gene,
      ggplot2::aes(
        x = 1.5,
        y = y_text,
        label = p_label
      ),
      inherit.aes = FALSE,
      family = FONT_FAMILY,
      size = 3.0,
      color = "black"
    ) +
    ggplot2::facet_wrap(
      ~ Dataset,
      nrow = 1,
      scales = "free_y"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        Normal = NORMAL_FILL,
        Tumor = TUMOR_FILL
      )
    ) +
    ggplot2::scale_color_manual(
      values = c(
        Normal = COL_NORMAL,
        Tumor = COL_TUMOR
      )
    ) +
    ggplot2::labs(
      title = gene_name,
      x = NULL,
      y = paste0(gene_name, " expression")
    ) +
    theme_manuscript(
      show_grid = TRUE,
      legend_position = "none"
    ) +
    ggplot2::theme() +
    ggplot2::coord_cartesian(clip = "off")
}

gene_plots <- lapply(
  as.character(final_biomarkers),
  make_gene_boxplot
)

names(gene_plots) <- as.character(final_biomarkers)

for (gene in names(gene_plots)) {
  safe_gene_name <- gsub("[^A-Za-z0-9_\\-]", "_", gene)

  save_plot_all_formats(
    plot_obj = gene_plots[[gene]],
    filename_stem = paste0(
      "Figure_05D_Expression_",
      safe_gene_name
    ),
    dir_path = FIGURE_DIR,
    width = FIG_DOUBLE_W,
    height = FIG_SINGLE_H
  )
}

fig5d_combined <- cowplot::plot_grid(
  plotlist = gene_plots,
  ncol = ifelse(length(gene_plots) == 1, 1, 2),
  align = "hv"
)

save_plot_all_formats(
  plot_obj = fig5d_combined,
  filename_stem = "Figure_05D_Final_biomarker_expression_boxplots",
  dir_path = FIGURE_DIR,
  width = ifelse(
    length(gene_plots) == 1,
    FIG_DOUBLE_W,
    FIGURE2_COMBINED_W
  ),
  height = ifelse(
    length(gene_plots) == 1,
    FIG_SINGLE_H,
    ceiling(length(gene_plots) / 2) * FIG_SINGLE_H
  )
)

fig5d_file <- file.path(
  FIGURE_DIR,
  "Figure_05D_Final_biomarker_expression_boxplots.png"
)

# FINAL FIGURE 5: REBUILD PANEL C WITH 4-DECIMAL AUC VALUES + COMBINE A-D

rm(list = intersect(
  c("train_df", "valid_df", "roc_train", "roc_valid", "p_fig5c",
    "fig5c_img", "fig5_combined_img"),
  ls()
))

options(stringsAsFactors = FALSE, scipen = 100)

# 1) PATHS AND REQUIRED PACKAGES

results_path <- "D:/LSCC/Results_LSCC/ML"
FIGURE_DIR <- results_path
FONT_FAMILY <- "Arial"
FIG_DPI <- 600
FIG_BACKGROUND <- "white"

required_packages <- c("pROC", "ggplot2", "magick")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install these package(s) first:\n",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
  library(magick)
})

if (.Platform$OS.type == "windows") {
  try(
    grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial")),
    silent = TRUE
  )
}

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# 2) HELPER FUNCTIONS

make_group_factor <- function(g) {
  factor(
    ifelse(
      as.character(g) %in% c("2", "Tumor", "tumor", "LSCC", "lscc", "Cancer", "cancer"),
      "Tumor",
      "Normal"
    ),
    levels = c("Normal", "Tumor")
  )
}

save_plot_all_formats <- function(plot_obj, filename_stem, dir_path,
                                  width, height) {

  png_file <- file.path(dir_path, paste0(filename_stem, ".png"))
  tiff_file <- file.path(dir_path, paste0(filename_stem, ".tiff"))
  pdf_file <- file.path(dir_path, paste0(filename_stem, ".pdf"))

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

  ggplot2::ggsave(
    filename = pdf_file,
    plot = plot_obj,
    width = width,
    height = height,
    device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf,
    bg = FIG_BACKGROUND,
    limitsize = FALSE
  )

  invisible(list(PNG = png_file, TIFF = tiff_file, PDF = pdf_file))
}

add_panel_label_img <- function(img,
                                label,
                                size = PANEL_TAG_RASTER_SIZE,
                                offset_x = PANEL_TAG_RASTER_OFFSET_X,
                                offset_y = PANEL_TAG_RASTER_OFFSET_Y) {
  magick::image_annotate(
    img,
    text = label,
    size = size,
    font = FONT_FAMILY,
    gravity = "northwest",
    location = paste0("+", offset_x, "+", offset_y),
    weight = 700,
    color = "black"
  )
}

make_clean_panel_img <- function(path,
                                 label,
                                 panel_width = 3600,
                                 panel_height = 2850,
                                 label_size = PANEL_TAG_RASTER_SIZE) {
  img <- magick::image_read(path)
  img <- magick::image_background(img, FIG_BACKGROUND, flatten = TRUE)
  img <- magick::image_trim(img)

  img <- magick::image_resize(
    img,
    geometry = paste0(
      panel_width - PANEL_CONTENT_INSET_PX, "x",
      panel_height - PANEL_CONTENT_VERTICAL_INSET_PX, ">"
    )
  )

  img <- magick::image_extent(
    img,
    geometry = paste0(panel_width, "x", panel_height),
    gravity = "center",
    color = FIG_BACKGROUND
  )

  img <- magick::image_border(
    img,
    color = "black",
    geometry = paste0(PANEL_BORDER_RASTER_PX, "x", PANEL_BORDER_RASTER_PX)
  )

  add_panel_label_img(
    img = img,
    label = label,
    size = label_size
  )
}

save_magick_all_formats <- function(image_object, filename_stem, dir_path) {

  png_file <- file.path(dir_path, paste0(filename_stem, ".png"))
  tiff_file <- file.path(dir_path, paste0(filename_stem, ".tiff"))
  pdf_file <- file.path(dir_path, paste0(filename_stem, ".pdf"))

  magick::image_write(
    image_object,
    path = png_file,
    format = "png",
    density = paste0(FIG_DPI, "x", FIG_DPI)
  )

  magick::image_write(
    image_object,
    path = tiff_file,
    format = "tiff",
    density = paste0(FIG_DPI, "x", FIG_DPI)
  )

  try(
    magick::image_write(
      image_object,
      path = pdf_file,
      format = "pdf",
      density = paste0(FIG_DPI, "x", FIG_DPI)
    ),
    silent = TRUE
  )

  invisible(list(PNG = png_file, TIFF = tiff_file, PDF = pdf_file))
}

# 3) READ TRAINING AND VALIDATION DATA

gene <- "MYBL2"

train_file <- file.path(results_path, "train_log2CPM.csv")
valid_file <- file.path(results_path, "valid_log2CPM.csv")

required_data_files <- c(train_file, valid_file)
missing_data_files <- required_data_files[!file.exists(required_data_files)]

if (length(missing_data_files) > 0) {
  stop(
    "Missing input file(s):\n",
    paste(missing_data_files, collapse = "\n")
  )
}

train_df <- read.csv(
  train_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

valid_df <- read.csv(
  valid_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c("Sample", "group", gene)

if (!all(required_columns %in% colnames(train_df))) {
  stop("Training file does not contain Sample, group, and MYBL2.")
}

if (!all(required_columns %in% colnames(valid_df))) {
  stop("Validation file does not contain Sample, group, and MYBL2.")
}

train_group <- make_group_factor(train_df$group)
valid_group <- make_group_factor(valid_df$group)

# 4) REBUILD PANEL C DIRECTLY FROM THE ORIGINAL DATA

roc_train <- pROC::roc(
  response = train_group,
  predictor = as.numeric(train_df[[gene]]),
  levels = c("Normal", "Tumor"),
  direction = "<",
  quiet = TRUE
)

roc_valid <- pROC::roc(
  response = valid_group,
  predictor = as.numeric(valid_df[[gene]]),
  levels = c("Normal", "Tumor"),
  direction = "<",
  quiet = TRUE
)

training_auc <- as.numeric(roc_train$auc)
validation_auc <- as.numeric(roc_valid$auc)

training_auc_label <- sprintf("Training AUC = %.4f", training_auc)
validation_auc_label <- sprintf("Validation AUC = %.4f", validation_auc)

cat("\nPanel C AUC values calculated from the data:\n")
cat(training_auc_label, "\n")
cat(validation_auc_label, "\n")

if (
  abs(training_auc - 0.9712) > 0.00005 ||
  abs(validation_auc - 0.9658) > 0.00005
) {
  warning(
    "The current input data do not round to 0.9712 and 0.9658. ",
    "The panel will show the true recalculated values instead."
  )
}

roc_list <- list(
  Training = roc_train,
  Validation = roc_valid
)

p_fig5c <- pROC::ggroc(
  roc_list,
  legacy.axes = TRUE,
  linewidth = 1.25
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey55",
    linewidth = 0.55
  ) +
  ggplot2::scale_color_manual(
    values = c(
      Training = COL_TUMOR,
      Validation = COL_NORMAL
    ),
    labels = c(
      Training = training_auc_label,
      Validation = validation_auc_label
    )
  ) +
  ggplot2::coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  ggplot2::labs(
    title = "MYBL2",
    x = "1 - Specificity",
    y = "Sensitivity",
    color = "Dataset"
  ) +
  theme_manuscript(
    show_grid = FALSE,
    legend_position = c(0.67, 0.17)
  ) +
  ggplot2::theme(
    legend.justification = c(0.5, 0.5),
    legend.background = ggplot2::element_rect(
      fill = FIG_BACKGROUND,
      color = "black",
      linewidth = PANEL_BORDER_LWD
    ),
    legend.key.width = grid::unit(0.55, "cm"),
    legend.key.height = grid::unit(0.45, "cm")
  )

save_plot_all_formats(
  plot_obj = p_fig5c,
  filename_stem = "Figure_05C_Internal_ROC_final_biomarker",
  dir_path = FIGURE_DIR,
  width = 7.2,
  height = 6.5
)

fig5c_file <- file.path(
  FIGURE_DIR,
  "Figure_05C_Internal_ROC_final_biomarker.png"
)

# 5) REQUIRED EXISTING PANELS A, B, AND D

fig5a_file <- file.path(
  FIGURE_DIR,
  "Figure_05A_Triple_overlap_Venn.png"
)

fig5b_file <- file.path(
  FIGURE_DIR,
  "Figure_05B_LASSO_cv_curve.png"
)

fig5d_file <- file.path(
  FIGURE_DIR,
  "Figure_05D_Final_biomarker_expression_boxplots.png"
)

final_panel_files <- c(
  fig5a_file,
  fig5b_file,
  fig5c_file,
  fig5d_file
)

missing_panel_files <- final_panel_files[!file.exists(final_panel_files)]

if (length(missing_panel_files) > 0) {
  stop(
    "The following figure panels are missing:\n",
    paste(missing_panel_files, collapse = "\n")
  )
}

# 6) FINAL COMBINED FIGURE 5

panel_a_width <- 3000
panel_b_width <- 4200

panel_c_width <- 2800
panel_d_width <- 4400

panel_top_height <- 2850
panel_bottom_height <- 2800

label_size <- 135

img_a <- make_clean_panel_img(
  path = fig5a_file,
  label = "a",
  panel_width = panel_a_width,
  panel_height = panel_top_height,
  label_size = label_size
)

img_b <- make_clean_panel_img(
  path = fig5b_file,
  label = "b",
  panel_width = panel_b_width,
  panel_height = panel_top_height,
  label_size = label_size
)

img_c <- make_clean_panel_img(
  path = fig5c_file,
  label = "c",
  panel_width = panel_c_width,
  panel_height = panel_bottom_height,
  label_size = label_size
)

img_d <- make_clean_panel_img(
  path = fig5d_file,
  label = "d",
  panel_width = panel_d_width,
  panel_height = panel_bottom_height,
  label_size = label_size
)

row_1 <- magick::image_append(
  c(img_a, img_b),
  stack = FALSE
)

row_2 <- magick::image_append(
  c(img_c, img_d),
  stack = FALSE
)

fig5_combined_img <- magick::image_append(
  c(row_1, row_2),
  stack = TRUE
)

fig5_combined_img <- magick::image_background(
  fig5_combined_img,
  FIG_BACKGROUND,
  flatten = TRUE
)

save_magick_all_formats(
  image_object = fig5_combined_img,
  filename_stem = "Figure_05_Screening_validation_MYBL2_combined",
  dir_path = FIGURE_DIR
)

cat("\nFINAL FIGURE 5 FINISHED.\n")
cat("Panel C legend:\n")
cat(training_auc_label, "\n")
cat(validation_auc_label, "\n")
cat("\nOutput folder:\n", FIGURE_DIR, "\n", sep = "")

# STEP 09. EXTERNAL MACHINE-LEARNING VALIDATION AND FIGURE 6A

if (!exists("results_path", inherits = FALSE)) {
  results_path <- "D:/LSCC/Results_LSCC/ML"
}

if (!exists("FIGURE_DIR", inherits = FALSE)) {
  FIGURE_DIR <- results_path
}

if (!exists("FIG_DPI", inherits = FALSE)) {
  FIG_DPI <- 600
}

if (!exists("FIG_BACKGROUND", inherits = FALSE)) {
  FIG_BACKGROUND <- "white"
}

if (.Platform$OS.type == "windows") {
  try(
    grDevices::windowsFonts(
      Arial = grDevices::windowsFont("Arial")
    ),
    silent = TRUE
  )
}

if (!exists("FONT_FAMILY", inherits = FALSE)) {
  FONT_FAMILY <- "Arial"
}

dir.create(results_path, recursive = TRUE, showWarnings = FALSE)

if (!exists("make_stage_dir", mode = "function", inherits = FALSE)) {
  make_stage_dir <- function(step_number, step_name) {
    dir.create(results_path, recursive = TRUE, showWarnings = FALSE)
    results_path
  }
}

if (!exists("clean_gene_vector", mode = "function", inherits = FALSE)) {
  clean_gene_vector <- function(x) {
    x <- unique(trimws(as.character(x)))
    x[!is.na(x) & x != ""]
  }
}

if (!exists("to_numeric_df", mode = "function", inherits = FALSE)) {
  to_numeric_df <- function(df) {
    as.data.frame(
      lapply(df, function(x) as.numeric(as.character(x))),
      check.names = FALSE
    )
  }
}

if (!exists("coerce_group_12", mode = "function", inherits = FALSE)) {
  coerce_group_12 <- function(g) {
    g0 <- as.character(g)

    if (all(g0 %in% c("1", "2"))) {
      return(as.integer(g0))
    }

    g_low <- tolower(g0)

    if (all(g_low %in% c("normal", "tumor", "non", "lscc", "cancer", "margin"))) {
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
        "group must have exactly 2 classes. Found: ",
        paste(u, collapse = ", ")
      )
    }

    map <- setNames(c(1L, 2L), u)
    as.integer(map[g0])
  }
}

if (!exists("make_group_factor", mode = "function", inherits = FALSE)) {
  make_group_factor <- function(g) {
    g12 <- coerce_group_12(g)
    factor(g12, levels = c(1, 2), labels = c("Normal", "Tumor"))
  }
}

required_packages <- c(
  "data.table",
  "dplyr",
  "ggplot2",
  "pROC",
  "randomForest",
  "e1071"
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
    "The following packages are required for Step 09:\n",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(pROC)
  library(randomForest)
  library(e1071)
})

if (!exists("step07_dir", inherits = FALSE)) {
  step07_dir <- results_path
}

if (!exists("train_log2_file", inherits = FALSE)) {
  train_log2_file <- file.path(results_path, "train_log2CPM.csv")
}

if (!exists("external_log2_file", inherits = FALSE)) {
  external_log2_file <- file.path(
    results_path,
    "external_GSE130605_log2CPM.csv"
  )
}

final_biomarkers_file <- file.path(
  step07_dir,
  "FINAL_BIOMARKERS.txt"
)

required_input_files <- c(
  train_log2_file,
  external_log2_file,
  final_biomarkers_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    "Step 09 requires output files from Steps 03 and 07. Missing file(s):\n",
    paste(missing_input_files, collapse = "\n"),
    "\n\nRun the complete pipeline through Step 07 before running Step 09."
  )
}

if (!exists("train_df", inherits = FALSE)) {
  train_df <- as.data.frame(
    data.table::fread(train_log2_file),
    check.names = FALSE
  )
}

if (!exists("test_df", inherits = FALSE)) {
  test_df <- as.data.frame(
    data.table::fread(external_log2_file),
    check.names = FALSE
  )
}

# STEP 09. EXTERNAL MACHINE-LEARNING VALIDATION AND FIGURE 6A

step09_dir <- make_stage_dir(9, "External_ML_Figure6A")

RF_NTREE <- 500
SVM_KERNEL <- "linear"

scale_by_train <- function(train_mat, new_mat) {
  train_mat <- as.matrix(train_mat)
  new_mat <- as.matrix(new_mat)

  mu <- colMeans(train_mat, na.rm = TRUE)
  sdv <- apply(train_mat, 2, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1

  sweep(sweep(new_mat, 2, mu, "-"), 2, sdv, "/")
}

metrics_binary <- function(y_true, y_pred, y_prob) {
  y_true <- factor(y_true, levels = c("Normal", "Tumor"))
  y_pred <- factor(y_pred, levels = c("Normal", "Tumor"))

  cm <- table(True = y_true, Pred = y_pred)

  TP <- ifelse(
    "Tumor" %in% rownames(cm) && "Tumor" %in% colnames(cm),
    cm["Tumor", "Tumor"],
    0
  )

  TN <- ifelse(
    "Normal" %in% rownames(cm) && "Normal" %in% colnames(cm),
    cm["Normal", "Normal"],
    0
  )

  FP <- ifelse(
    "Normal" %in% rownames(cm) && "Tumor" %in% colnames(cm),
    cm["Normal", "Tumor"],
    0
  )

  FN <- ifelse(
    "Tumor" %in% rownames(cm) && "Normal" %in% colnames(cm),
    cm["Tumor", "Normal"],
    0
  )

  accuracy <- (TP + TN) / sum(cm)
  sensitivity <- ifelse((TP + FN) == 0, NA_real_, TP / (TP + FN))
  specificity <- ifelse((TN + FP) == 0, NA_real_, TN / (TN + FP))
  precision <- ifelse((TP + FP) == 0, NA_real_, TP / (TP + FP))
  f1 <- ifelse(
    is.na(precision) | is.na(sensitivity) | (precision + sensitivity) == 0,
    NA_real_,
    2 * precision * sensitivity / (precision + sensitivity)
  )

  roc_obj <- pROC::roc(
    response = y_true,
    predictor = as.numeric(y_prob),
    levels = c("Normal", "Tumor"),
    direction = "<",
    quiet = TRUE
  )

  data.frame(
    AUROC = as.numeric(roc_obj$auc),
    Accuracy = as.numeric(accuracy),
    Sensitivity = as.numeric(sensitivity),
    Specificity = as.numeric(specificity),
    Precision = as.numeric(precision),
    F1 = as.numeric(f1)
  )
}

final_biomarkers <- clean_gene_vector(
  readLines(final_biomarkers_file, warn = FALSE)
)

model_genes <- intersect(
  final_biomarkers,
  Reduce(
    intersect,
    list(
      setdiff(colnames(train_df), c("Sample", "group", "batch")),
      setdiff(colnames(test_df), c("Sample", "group", "batch"))
    )
  )
)

model_genes <- sort(unique(model_genes))

if (length(model_genes) < 1) {
  stop("No final biomarker is available in both train and external matrices.")
}

write.csv(
  data.frame(Gene = model_genes),
  file.path(step09_dir, "Genes_used_for_external_ML.csv"),
  row.names = FALSE
)

trainX_ml <- as.matrix(
  to_numeric_df(
    train_df[, model_genes, drop = FALSE]
  )
)

testX_ml <- as.matrix(
  to_numeric_df(
    test_df[, model_genes, drop = FALSE]
  )
)

rownames(trainX_ml) <- train_df$Sample
rownames(testX_ml) <- test_df$Sample

trainY_ml <- make_group_factor(train_df$group)
testY_ml <- make_group_factor(test_df$group)

train_ml_scaled <- scale_by_train(trainX_ml, trainX_ml)
test_ml_scaled <- scale_by_train(trainX_ml, testX_ml)

train_ml_df <- data.frame(
  group = trainY_ml,
  train_ml_scaled,
  check.names = FALSE
)

test_ml_df <- data.frame(
  group = testY_ml,
  test_ml_scaled,
  check.names = FALSE
)

colnames(train_ml_df) <- make.names(
  colnames(train_ml_df),
  unique = TRUE
)

colnames(test_ml_df) <- make.names(
  colnames(test_ml_df),
  unique = TRUE
)

set.seed(123)

rf_model <- randomForest::randomForest(
  group ~ .,
  data = train_ml_df,
  ntree = RF_NTREE,
  importance = TRUE
)

rf_prob <- predict(
  rf_model,
  newdata = test_ml_df,
  type = "prob"
)[, "Tumor"]

rf_pred <- predict(
  rf_model,
  newdata = test_ml_df,
  type = "response"
)

rf_metrics <- metrics_binary(testY_ml, rf_pred, rf_prob)
rf_metrics$Model <- "RF"

set.seed(123)

svm_model <- e1071::svm(
  group ~ .,
  data = train_ml_df,
  kernel = SVM_KERNEL,
  probability = TRUE
)

svm_pred <- predict(
  svm_model,
  newdata = test_ml_df,
  probability = TRUE
)

svm_prob_mat <- attr(svm_pred, "probabilities")

svm_prob <- if ("Tumor" %in% colnames(svm_prob_mat)) {
  svm_prob_mat[, "Tumor"]
} else {
  1 - svm_prob_mat[, "Normal"]
}

svm_metrics <- metrics_binary(testY_ml, svm_pred, svm_prob)
svm_metrics$Model <- "SVM"

ml_metrics <- dplyr::bind_rows(
  rf_metrics,
  svm_metrics
) %>%
  dplyr::select(Model, everything())

write.csv(
  ml_metrics,
  file.path(step09_dir, "External_GSE130605_ML_metrics.csv"),
  row.names = FALSE
)

roc_models <- list(
  RF = pROC::roc(
    testY_ml,
    rf_prob,
    levels = c("Normal", "Tumor"),
    direction = "<",
    quiet = TRUE
  ),
  SVM = pROC::roc(
    testY_ml,
    svm_prob,
    levels = c("Normal", "Tumor"),
    direction = "<",
    quiet = TRUE
  )
)

legend_labels <- sapply(
  names(roc_models),
  function(nm) {
    paste0(
      nm,
      " (AUC = ",
      sprintf("%.2f%%", 100 * as.numeric(pROC::auc(roc_models[[nm]]))),
      ")"
    )
  }
)

roc_colors <- c(
  RF = "#5B7FC7",
  SVM = "#E6007E"
)

p_fig6a <- pROC::ggroc(
  roc_models,
  legacy.axes = TRUE,
  linewidth = 1.35
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.8
  ) +
  ggplot2::scale_color_manual(
    values = roc_colors,
    labels = legend_labels
  ) +
  ggplot2::coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  ggplot2::labs(
    tag = "a",
    x = "1 - Specificity",
    y = "Sensitivity",
    color = NULL,
    title = "External validation in GSE130605"
  ) +
  theme_manuscript(
    show_grid = FALSE,
    legend_position = c(0.67, 0.18)
  ) +
  ggplot2::theme(
    legend.background = ggplot2::element_rect(
      fill = FIG_BACKGROUND,
      color = "black",
      linewidth = PANEL_BORDER_LWD
    )
  )

fig6a_file <- file.path(
  step09_dir,
  "Figure_06A_External_validation_RF_SVM_ROC.png"
)

ggplot2::ggsave(
  filename = fig6a_file,
  plot = p_fig6a,
  width = 7.2,
  height = 6.5,
  dpi = FIG_DPI,
  bg = "white",
  limitsize = FALSE
)

cat("\nSTEP 09 finished: external RF and SVM validation and Figure 6A generated.\n")
print(ml_metrics)

# MODULE 02 — Immune microenvironment analysis and Figure 6 outputs

# FINAL FIGURE 6 PANEL ORDER:

# STEP 12.1. CLEAN ENVIRONMENT, SETTINGS, PACKAGES, AND PATHS

rm(list = ls())
gc()

set.seed(123)
options(stringsAsFactors = FALSE)
options(scipen = 100)
options(timeout = 7200)

TARGET_GENE <- "MYBL2"

FONT_FAMILY <- "Arial"
FIG_DPI <- 600
FIG_BACKGROUND <- "white"

FIG_SINGLE_W <- 3.50
FIG_SINGLE_H <- 4.20
FIG_DOUBLE_W <- 7.20
FIG_DOUBLE_H <- 5.40
FIGURE2_COMBINED_W <- 14.00
FIGURE2_COMBINED_H <- 6.40

BASE_TEXT_PT <- 9.5
AXIS_TEXT_PT <- 8.5
AXIS_TITLE_PT <- 10
LEGEND_TEXT_PT <- 8.5
LEGEND_TITLE_PT <- 9
PLOT_TITLE_PT <- 10
PANEL_TAG_PT <- 14

PANEL_BORDER_LWD <- 0.45
AXIS_LWD <- 0.40
TICK_LWD <- 0.40
GEOM_LWD <- 0.55
GRID_LWD <- 0.30
TICK_LENGTH_PT <- 2
PANEL_SPACING_PT <- 6

RASTER_PANEL_TAG_SIZE <- 145
RASTER_PANEL_TAG_OFFSET_X <- 100
RASTER_PANEL_TAG_OFFSET_Y <- 82
RASTER_PANEL_BORDER_PX <- 3
RASTER_PANEL_INSET_X <- 280
RASTER_PANEL_INSET_Y <- 80

HEATMAP_FONT_PT <- 8
HEATMAP_ROW_FONT_PT <- 7
HEATMAP_COL_FONT_PT <- 7
HEATMAP_LEGEND_FONT_PT <- 8

COL_NORMAL <- "#5B6FB2"
COL_TUMOR <- "#CE5C8A"
COL_MYBL2_LOW <- COL_NORMAL
COL_MYBL2_HIGH <- COL_TUMOR
COL_UP <- COL_TUMOR
COL_DOWN <- COL_NORMAL
COL_NS <- "#D1D3DA"

CELLCHAT_PALETTE <- grDevices::colorRampPalette(
  c("#F4F1F8", "#DDD5EA", "#BDA9D5", "#8B70B4", "#573D7E")
)(100)

MANUSCRIPT_MARGIN <- ggplot2::margin(
  t = 7,
  r = 8,
  b = 7,
  l = 8,
  unit = "pt"
)

LEGEND_BOX_MARGIN <- ggplot2::margin(
  t = 0,
  r = 0,
  b = 0,
  l = 0,
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

FDR_CUTOFF <- 0.05
CORRELATION_RHO_CUTOFF <- 0.30

PANEL_B_N_SIGNATURES <- 15L
PANEL_B_X_TEXT_SIZE <- 8.0
PANEL_B_X_TEXT_ANGLE <- 45
PANEL_B_LABEL_WRAP_WIDTH <- 19
PANEL_B_BOX_WIDTH <- 0.42
PANEL_B_DODGE_WIDTH <- 0.52
PANEL_B_POINT_SIZE <- 0.72

CHAROENTONG_MMC3_FILE <- "D:/LSCC/Results_LSCC/mmc3.xlsx"
CHAROENTONG_SHEET <- "Sheet1"
CHAROENTONG_REFERENCE_CITATION <- paste0(
  "Charoentong et al. 2017, Cell Reports 18:248-262; ",
  "Table S6 (Pancancer immune metagenes), local mmc3.xlsx"
)
MIN_SIGNATURE_GENE_OVERLAP <- 3

ml_input_dir <- "D:/LSCC/Results_LSCC/ML"

immune_output_dir <- "D:/LSCC/Results_LSCC/Immune Microenvironment"
dir.create(immune_output_dir, recursive = TRUE, showWarnings = FALSE)

step12_dir <- immune_output_dir
step12_figdir <- immune_output_dir
step12_tabledir <- immune_output_dir
step12_rdsdir <- immune_output_dir
manuscript_figures_path <- immune_output_dir

required_packages <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "patchwork",
  "tibble", "stringr", "GSVA", "IOBR", "readxl", "magick",
  "randomForest", "e1071", "pROC"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install these package(s) before running this script:\n",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(stringr)
  library(GSVA)
  library(IOBR)
  library(readxl)
  library(magick)
  library(randomForest)
  library(e1071)
  library(pROC)
})

# STEP 12.2. HELPER FUNCTIONS

first_existing <- function(paths) {
  paths <- paths[!is.na(paths) & file.exists(paths)]
  if (length(paths) == 0) return(NA_character_)
  paths[1]
}

to_numeric_df <- function(df) {
  as.data.frame(
    lapply(df, function(x) suppressWarnings(as.numeric(as.character(x)))),
    check.names = FALSE
  )
}

make_group_factor <- function(g) {
  g0 <- as.character(g)
  g_low <- tolower(g0)

  if (all(g0 %in% c("1", "2"))) {
    return(factor(g0, levels = c("1", "2"), labels = c("Normal", "Tumor")))
  }

  if (all(g_low %in% c("normal", "tumor", "non", "lscc", "cancer", "margin"))) {
    out <- ifelse(g_low %in% c("tumor", "lscc", "cancer"), "Tumor", "Normal")
    return(factor(out, levels = c("Normal", "Tumor")))
  }

  u <- sort(unique(g0))
  if (length(u) != 2) {
    stop("The bulk group column must contain exactly two classes. Found: ",
         paste(u, collapse = ", "))
  }

  map <- setNames(c("Normal", "Tumor"), u)
  factor(map[g0], levels = c("Normal", "Tumor"))
}

make_mybl2_high_low <- function(x) {
  med <- median(as.numeric(x), na.rm = TRUE)
  factor(
    ifelse(as.numeric(x) >= med, "MYBL2-high", "MYBL2-low"),
    levels = c("MYBL2-low", "MYBL2-high")
  )
}

format_p <- function(p) {
  if (!is.finite(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", sprintf("%.3f", p))
}

format_fdr <- function(fdr) {
  if (!is.finite(fdr)) return("FDR = NA")
  if (fdr < 0.001) return("FDR < 0.001")
  paste0("FDR = ", sprintf("%.3f", fdr))
}

fdr_star <- function(fdr) {
  if (!is.finite(fdr)) return("")
  if (fdr < 0.001) return("***")
  if (fdr < 0.010) return("**")
  if (fdr < 0.050) return("*")
  ""
}

full_panel_b_label <- function(x) {
  stringr::str_wrap(
    as.character(x),
    width = PANEL_B_LABEL_WRAP_WIDTH
  )
}

safe_wilcox <- function(x, group) {
  keep <- is.finite(x) & !is.na(group)
  x <- x[keep]
  group <- droplevels(factor(group[keep]))

  if (length(unique(group)) < 2 || min(table(group)) < 2) return(NA_real_)

  tryCatch(
    wilcox.test(x ~ group, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5) {
    return(data.frame(rho = NA_real_, p_value = NA_real_, n = sum(keep)))
  }

  out <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  data.frame(
    rho = as.numeric(out$estimate),
    p_value = as.numeric(out$p.value),
    n = sum(keep)
  )
}

save_plot <- function(plot_obj, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path(step12_figdir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = FIG_DPI,
    bg = "white",
    limitsize = FALSE
  )
}

save_manuscript_copy <- function(from_file, figure_name) {
  if (!file.exists(from_file)) return(invisible(NULL))

  to_file <- file.path(manuscript_figures_path, figure_name)

  if (identical(
    normalizePath(from_file, winslash = "/", mustWork = FALSE),
    normalizePath(to_file, winslash = "/", mustWork = FALSE)
  )) {
    return(invisible(to_file))
  }

  file.copy(from = from_file, to = to_file, overwrite = TRUE)
  invisible(to_file)
}

add_panel_label_img <- function(img,
                                label,
                                size = RASTER_PANEL_TAG_SIZE,
                                offset_x = RASTER_PANEL_TAG_OFFSET_X,
                                offset_y = RASTER_PANEL_TAG_OFFSET_Y) {
  magick::image_annotate(
    img,
    text = label,
    size = size,
    font = FONT_FAMILY,
    gravity = "northwest",
    location = paste0("+", offset_x, "+", offset_y),
    weight = 700,
    color = "black"
  )
}

make_clean_panel_img <- function(path,
                                 label,
                                 panel_width,
                                 panel_height,
                                 label_size = RASTER_PANEL_TAG_SIZE) {
  img <- magick::image_read(path)
  img <- magick::image_background(img, FIG_BACKGROUND, flatten = TRUE)
  img <- magick::image_trim(img)

  img <- magick::image_resize(
    img,
    geometry = paste0(
      panel_width - RASTER_PANEL_INSET_X,
      "x",
      panel_height - RASTER_PANEL_INSET_Y,
      ">"
    )
  )

  img <- magick::image_extent(
    img,
    geometry = paste0(panel_width, "x", panel_height),
    gravity = "center",
    color = FIG_BACKGROUND
  )

  img <- magick::image_border(
    img,
    color = "black",
    geometry = paste0(
      RASTER_PANEL_BORDER_PX,
      "x",
      RASTER_PANEL_BORDER_PX
    )
  )

  add_panel_label_img(
    img = img,
    label = label,
    size = label_size
  )
}

save_magick_all_formats <- function(image_object,
                                    filename_stem,
                                    dir_path) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  png_file <- file.path(dir_path, paste0(filename_stem, ".png"))
  tiff_file <- file.path(dir_path, paste0(filename_stem, ".tiff"))
  pdf_file <- file.path(dir_path, paste0(filename_stem, ".pdf"))

  magick::image_write(
    image_object,
    path = png_file,
    format = "png",
    density = paste0(FIG_DPI, "x", FIG_DPI)
  )

  magick::image_write(
    image_object,
    path = tiff_file,
    format = "tiff",
    density = paste0(FIG_DPI, "x", FIG_DPI)
  )

  try(
    magick::image_write(
      image_object,
      path = pdf_file,
      format = "pdf",
      density = paste0(FIG_DPI, "x", FIG_DPI)
    ),
    silent = TRUE
  )

  invisible(
    list(
      PNG = png_file,
      TIFF = tiff_file,
      PDF = pdf_file
    )
  )
}

make_clean_theme <- function(base_size = BASE_TEXT_PT,
                             legend_position = "right",
                             show_grid = FALSE) {
  ggplot2::theme_bw(
    base_size = base_size,
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
        vjust = 0.5,
        color = "black"
      ),
      plot.title.position = "panel",

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
        linewidth = TICK_LWD
      ),
      axis.ticks.length = grid::unit(TICK_LENGTH_PT, "pt"),

      strip.background = ggplot2::element_rect(
        fill = FIG_BACKGROUND,
        color = "black",
        linewidth = PANEL_BORDER_LWD
      ),
      strip.text = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PLOT_TITLE_PT,
        color = "black"
      ),

      legend.position = legend_position,
      legend.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = LEGEND_TITLE_PT,
        color = "black"
      ),
      legend.text = ggplot2::element_text(
        family = FONT_FAMILY,
        size = LEGEND_TEXT_PT,
        color = "black"
      ),
      legend.background = ggplot2::element_rect(
        fill = FIG_BACKGROUND,
        color = NA
      ),
      legend.key = ggplot2::element_rect(
        fill = FIG_BACKGROUND,
        color = NA
      ),
      legend.key.width = grid::unit(0.38, "cm"),
      legend.key.height = grid::unit(0.38, "cm"),
      legend.box.margin = LEGEND_BOX_MARGIN,
      legend.box.spacing = grid::unit(PANEL_SPACING_PT, "pt"),

      panel.background = ggplot2::element_rect(
        fill = FIG_BACKGROUND,
        color = NA
      ),
      plot.background = ggplot2::element_rect(
        fill = FIG_BACKGROUND,
        color = NA
      ),
      panel.border = ggplot2::element_rect(
        colour = "black",
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

      panel.spacing = grid::unit(PANEL_SPACING_PT, "pt"),
      plot.margin = MANUSCRIPT_MARGIN
    )
}

normalize_signature_name <- function(x) {
  x <- tolower(as.character(x))
  gsub("[^a-z0-9]+", "", x)
}

immune28_labels <- c(
  "Activated B cell",
  "Activated CD4 T cell",
  "Activated CD8 T cell",
  "Activated dendritic cell",
  "CD56bright natural killer cell",
  "CD56dim natural killer cell",
  "Central memory CD4 T cell",
  "Central memory CD8 T cell",
  "Effector memory CD4 T cell",
  "Effector memory CD8 T cell",
  "Eosinophil",
  "Gamma delta T cell",
  "Immature B cell",
  "Immature dendritic cell",
  "Macrophage",
  "Mast cell",
  "MDSC",
  "Memory B cell",
  "Monocyte",
  "Natural killer cell",
  "Natural killer T cell",
  "Neutrophil",
  "Plasmacytoid dendritic cell",
  "Regulatory T cell",
  "T follicular helper cell",
  "Type 1 T helper cell",
  "Type 17 T helper cell",
  "Type 2 T helper cell"
)

immune28_aliases <- list(
  "Activated B cell" = c("activatedbcell", "activatedbcells"),
  "Activated CD4 T cell" = c("activatedcd4tcell", "activatedcd4tcells"),
  "Activated CD8 T cell" = c("activatedcd8tcell", "activatedcd8tcells"),
  "Activated dendritic cell" = c("activateddendriticcell", "activateddendriticcells", "activateddc"),
  "CD56bright natural killer cell" = c("cd56brightnaturalkillercell", "cd56brightnaturalkillercells", "cd56brightnkcell", "cd56brightnkcells"),
  "CD56dim natural killer cell" = c("cd56dimnaturalkillercell", "cd56dimnaturalkillercells", "cd56dimnkcell", "cd56dimnkcells"),
  "Central memory CD4 T cell" = c("centralmemorycd4tcell", "centralmemorycd4tcells", "cd4tcm"),
  "Central memory CD8 T cell" = c("centralmemorycd8tcell", "centralmemorycd8tcells", "cd8tcm"),
  "Effector memory CD4 T cell" = c("effectormemorycd4tcell", "effectormemorycd4tcells", "effectormemeorycd4tcell", "effectormemeorycd4tcells", "cd4tem"),
  "Effector memory CD8 T cell" = c("effectormemorycd8tcell", "effectormemorycd8tcells", "effectormemeorycd8tcell", "effectormemeorycd8tcells", "cd8tem"),
  "Eosinophil" = c("eosinophil", "eosinophils"),
  "Gamma delta T cell" = c("gammadeltatcell", "gammadeltatcells", "gammadeltat"),
  "Immature B cell" = c("immaturebcell", "immaturebcells"),
  "Immature dendritic cell" = c("immaturedendriticcell", "immaturedendriticcells", "immaturedc"),
  "Macrophage" = c("macrophage", "macrophages"),
  "Mast cell" = c("mastcell", "mastcells"),
  "MDSC" = c("mdsc", "mdscs", "myeloidderivedsuppressorcell", "myeloidderivedsuppressorcells"),
  "Memory B cell" = c("memorybcell", "memorybcells"),
  "Monocyte" = c("monocyte", "monocytes"),
  "Natural killer cell" = c("naturalkillercell", "naturalkillercells", "nkcell", "nkcells"),
  "Natural killer T cell" = c("naturalkillertcell", "naturalkillertcells", "nktcell", "nktcells"),
  "Neutrophil" = c("neutrophil", "neutrophils"),
  "Plasmacytoid dendritic cell" = c("plasmacytoiddendriticcell", "plasmacytoiddendriticcells", "pdc", "pdcs"),
  "Regulatory T cell" = c("regulatorytcell", "regulatorytcells", "treg", "tregs"),
  "T follicular helper cell" = c("tfollicularhelpercell", "tfollicularhelpercells", "tfh", "tfhs"),
  "Type 1 T helper cell" = c("type1thelpercell", "type1thelpercells", "th1"),
  "Type 17 T helper cell" = c("type17thelpercell", "type17thelpercells", "th17"),
  "Type 2 T helper cell" = c("type2thelpercell", "type2thelpercells", "th2")
)

canonicalize_immune28_label <- function(x) {
  x_norm <- normalize_signature_name(x)
  output <- rep(NA_character_, length(x_norm))

  for (label in immune28_labels) {
    hit <- x_norm %in% immune28_aliases[[label]]
    output[is.na(output) & hit] <- label
  }

  output
}

if (!identical(
  canonicalize_immune28_label(c(
    "Effector memeory CD4 T cell",
    "Effector memeory CD8 T cell"
  )),
  c("Effector memory CD4 T cell", "Effector memory CD8 T cell")
)) {
  stop("Internal alias check failed for the known 'memeory' spelling in mmc3.xlsx.")
}

write_gmt <- function(gene_sets, gmt_file) {
  gmt_lines <- vapply(
    names(gene_sets),
    function(label) {
      paste(
        c(
          label,
          "Charoentong_et_al_2017_Cell_Reports_Table_S6_mmc3",
          unique(gene_sets[[label]])
        ),
        collapse = "\t"
      )
    },
    character(1)
  )

  writeLines(gmt_lines, con = gmt_file, useBytes = TRUE)
  invisible(gmt_file)
}

read_local_charoentong_immune28 <- function() {
  if (!file.exists(CHAROENTONG_MMC3_FILE)) {
    stop(
      "Local Charoentong source file was not found:\\n",
      CHAROENTONG_MMC3_FILE,
      "\\n\\nPlace mmc3.xlsx in D:/LSCC/Results_LSCC and run again."
    )
  }

  available_sheets <- readxl::excel_sheets(CHAROENTONG_MMC3_FILE)

  if (!(CHAROENTONG_SHEET %in% available_sheets)) {
    stop(
      "Expected sheet '", CHAROENTONG_SHEET, "' was not found in mmc3.xlsx.\\n",
      "Available sheets: ", paste(available_sheets, collapse = ", ")
    )
  }

  raw_tbl <- suppressMessages(
    readxl::read_excel(
      path = CHAROENTONG_MMC3_FILE,
      sheet = CHAROENTONG_SHEET,
      col_names = FALSE,
      col_types = "text",
      .name_repair = "minimal"
    )
  )

  raw_df <- as.data.frame(raw_tbl, stringsAsFactors = FALSE, check.names = FALSE)

  if (ncol(raw_df) < 2 || nrow(raw_df) < 4) {
    stop(
      "The local mmc3.xlsx Sheet1 does not contain the expected Table S6 structure."
    )
  }

  first_col_norm <- normalize_signature_name(raw_df[[1]])
  second_col_norm <- normalize_signature_name(raw_df[[2]])

  header_row <- which(
    first_col_norm %in% c("metagene", "gene", "genesymbol", "symbol") &
      second_col_norm %in% c("celltype", "immunecell", "immunecelltype")
  )[1]

  if (is.na(header_row) || header_row >= nrow(raw_df)) {
    stop(
      "Could not locate the official 'Metagene | Cell type' header row in Sheet1."
    )
  }

  data_rows <- seq.int(header_row + 1L, nrow(raw_df))

  reference_df <- data.frame(
    Metagene = toupper(trimws(as.character(raw_df[[1]][data_rows]))),
    Source_cell_label = trimws(as.character(raw_df[[2]][data_rows])),
    Immunity_class = if (ncol(raw_df) >= 3) {
      trimws(as.character(raw_df[[3]][data_rows]))
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(
      !is.na(Metagene),
      nzchar(Metagene),
      !is.na(Source_cell_label),
      nzchar(Source_cell_label)
    ) %>%
    dplyr::mutate(
      Immune_cell = canonicalize_immune28_label(Source_cell_label)
    ) %>%
    dplyr::filter(!is.na(Immune_cell)) %>%
    dplyr::distinct(Immune_cell, Metagene, .keep_all = TRUE)

  detected_labels <- sort(unique(as.character(reference_df$Immune_cell)))
  missing_labels <- setdiff(immune28_labels, detected_labels)

  if (length(missing_labels) > 0) {
    original_labels <- sort(unique(trimws(as.character(raw_df[[2]][data_rows]))))
    original_labels <- original_labels[!is.na(original_labels) & nzchar(original_labels)]

    stop(
      "The local Charoentong Table S6 did not yield all 28 required immune-cell signatures.\\n",
      "Missing: ", paste(missing_labels, collapse = ", "),
      "\\n\\nDetected original cell-type labels in Sheet1:\\n",
      paste(original_labels, collapse = " | ")
    )
  }

  reference_df <- reference_df %>%
    dplyr::mutate(
      Immune_cell = factor(Immune_cell, levels = immune28_labels)
    ) %>%
    dplyr::arrange(Immune_cell, Metagene)

  gene_sets <- lapply(immune28_labels, function(label) {
    unique(
      as.character(
        reference_df$Metagene[
          as.character(reference_df$Immune_cell) == label
        ]
      )
    )
  })
  names(gene_sets) <- immune28_labels

  source_names <- vapply(
    immune28_labels,
    function(label) {
      values <- unique(
        as.character(
          reference_df$Source_cell_label[
            as.character(reference_df$Immune_cell) == label
          ]
        )
      )
      values[1]
    },
    character(1)
  )
  names(source_names) <- immune28_labels

  list(
    gene_sets = gene_sets,
    reference_table = reference_df,
    source_names = source_names,
    source_file = CHAROENTONG_MMC3_FILE,
    source_sheet = CHAROENTONG_SHEET,
    header_row = header_row
  )
}

run_ssgsea_safe <- function(expr_mat_gene_by_sample, gene_sets,
                            min_size = MIN_SIGNATURE_GENE_OVERLAP) {
  expr_mat_gene_by_sample <- as.matrix(expr_mat_gene_by_sample)
  storage.mode(expr_mat_gene_by_sample) <- "numeric"
  expr_mat_gene_by_sample[!is.finite(expr_mat_gene_by_sample)] <- 0

  out <- tryCatch(
    {
      param <- GSVA::ssgseaParam(
        exprData = expr_mat_gene_by_sample,
        geneSets = gene_sets,
        minSize = min_size,
        maxSize = 500,
        normalize = TRUE
      )
      GSVA::gsva(param, verbose = FALSE)
    },
    error = function(e_new) {
      message("New GSVA ssGSEA API failed. Trying legacy API. Reason: ", e_new$message)
      tryCatch(
        {
          GSVA::gsva(
            expr = expr_mat_gene_by_sample,
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
            "ssGSEA failed with both GSVA APIs.\\n",
            "New API: ", e_new$message, "\\n",
            "Legacy API: ", e_old$message
          )
        }
      )
    }
  )

  as.matrix(out)
}

run_estimate_safe <- function(expr_gene_by_sample) {
  expr_gene_by_sample <- as.matrix(expr_gene_by_sample)
  storage.mode(expr_gene_by_sample) <- "numeric"
  expr_gene_by_sample[!is.finite(expr_gene_by_sample)] <- 0
  expr_gene_by_sample[expr_gene_by_sample < 0] <- 0

  if (nrow(expr_gene_by_sample) < 100 || ncol(expr_gene_by_sample) < 5) {
    stop("ESTIMATE input matrix has insufficient genes or samples.")
  }

  estimate_res <- tryCatch(
    IOBR::deconvo_tme(
      eset = expr_gene_by_sample,
      method = "estimate"
    ),
    error = function(e) {
      stop("ESTIMATE calculation through IOBR failed: ", e$message)
    }
  )

  estimate_res <- as.data.frame(estimate_res, check.names = FALSE)
  if (!"ID" %in% colnames(estimate_res)) {
    estimate_res$ID <- rownames(estimate_res)
  }
  estimate_res
}

# STEP 12.3. POOLED TUMOR COHORT WITHOUT BATCH-EFFECT CORRECTION

discovery_file <- file.path(
  ml_input_dir,
  "train_discovery_merged_CPM_nonnegative.csv"
)

if (!file.exists(discovery_file)) {
  stop(
    "The immune-analysis input file was not found:\n",
    discovery_file,
    "\n\nRun the second bulk RNA-seq/ML pipeline through Step 02 first."
  )
}

discovery_df <- as.data.frame(
  data.table::fread(discovery_file),
  check.names = FALSE
)

required_metadata <- c("Sample", "group", "batch")
if (!all(required_metadata %in% colnames(discovery_df))) {
  stop(
    "The pooled discovery CPM matrix must contain Sample, group, and batch columns."
  )
}

discovery_df$Sample <- trimws(as.character(discovery_df$Sample))
discovery_df$Tumor_Normal_Group <- make_group_factor(discovery_df$group)
discovery_df$batch <- suppressWarnings(as.integer(as.character(discovery_df$batch)))

if (any(!discovery_df$batch %in% c(1L, 2L))) {
  stop(
    "This immune workflow expects original discovery batches 1 and 2 only ",
    "(GSE127165 and GSE142083)."
  )
}

dataset_map <- c(
  "1" = "GSE127165",
  "2" = "GSE142083"
)

tumor_raw_df <- discovery_df %>%
  dplyr::filter(Tumor_Normal_Group == "Tumor") %>%
  dplyr::mutate(
    Dataset = factor(
      unname(dataset_map[as.character(batch)]),
      levels = c("GSE127165", "GSE142083")
    )
  )

if (nrow(tumor_raw_df) < 10) {
  stop("Too few tumor samples were retained for immune-microenvironment analysis.")
}

if (length(unique(tumor_raw_df$batch)) != 2) {
  stop(
    "Both GSE127165 and GSE142083 must contribute tumor samples to the pooled analysis."
  )
}

metadata_columns <- c(
  "Sample", "group", "batch", "Tumor_Normal_Group", "Dataset"
)

common_genes <- setdiff(colnames(tumor_raw_df), metadata_columns)

if (!(TARGET_GENE %in% common_genes)) {
  stop(TARGET_GENE, " was not found in the pooled discovery expression matrix.")
}

expr_cpm_sample_by_gene <- as.matrix(
  to_numeric_df(tumor_raw_df[, common_genes, drop = FALSE])
)

rownames(expr_cpm_sample_by_gene) <- tumor_raw_df$Sample
storage.mode(expr_cpm_sample_by_gene) <- "numeric"
expr_cpm_sample_by_gene[!is.finite(expr_cpm_sample_by_gene)] <- 0
expr_cpm_sample_by_gene[expr_cpm_sample_by_gene < 0] <- 0

if (anyDuplicated(rownames(expr_cpm_sample_by_gene)) > 0) {
  stop("Duplicated Sample IDs were detected in the pooled tumor cohort.")
}

expr_sample_by_gene <- log2(expr_cpm_sample_by_gene + 1)
expr_sample_by_gene <- as.matrix(expr_sample_by_gene)
storage.mode(expr_sample_by_gene) <- "numeric"

gene_sd <- apply(expr_sample_by_gene, 2, stats::sd, na.rm = TRUE)
keep_genes <- names(gene_sd)[is.finite(gene_sd) & gene_sd > 0]

if (!(TARGET_GENE %in% keep_genes)) {
  stop(TARGET_GENE, " is absent or has zero variance in the pooled tumor cohort.")
}

expr_sample_by_gene <- expr_sample_by_gene[, keep_genes, drop = FALSE]
estimate_cpm_sample_by_gene <- expr_cpm_sample_by_gene[, keep_genes, drop = FALSE]

mybl2_expression <- as.numeric(expr_sample_by_gene[, TARGET_GENE])
names(mybl2_expression) <- rownames(expr_sample_by_gene)

global_mybl2_median <- median(mybl2_expression, na.rm = TRUE)

mybl2_group <- factor(
  ifelse(
    mybl2_expression >= global_mybl2_median,
    "MYBL2-high",
    "MYBL2-low"
  ),
  levels = c("MYBL2-low", "MYBL2-high")
)

sample_metadata <- data.frame(
  Sample = rownames(expr_sample_by_gene),
  Dataset = factor(
    as.character(tumor_raw_df$Dataset),
    levels = c("GSE127165", "GSE142083")
  ),
  Batch = tumor_raw_df$batch,
  MYBL2_expression = mybl2_expression,
  MYBL2_group = mybl2_group,
  stringsAsFactors = FALSE
)

if (length(unique(sample_metadata$MYBL2_group)) != 2) {
  stop("MYBL2-high and MYBL2-low groups could not both be formed.")
}

write.csv(
  sample_metadata,
  file.path(
    step12_tabledir,
    "Pooled_tumor_samples_no_batch_correction_MYBL2_high_low_metadata.csv"
  ),
  row.names = FALSE
)

writeLines(
  c(
    "Pooled tumor cohort and MYBL2 split",
    paste0("Input TMM-CPM matrix: ", discovery_file),
    "Tumor samples from GSE127165 and GSE142083 were pooled.",
    "No batch-effect correction was applied.",
    "Batch 1 = GSE127165; Batch 2 = GSE142083.",
    paste0("Global MYBL2 median on log2(TMM-CPM + 1): ", format(global_mybl2_median, digits = 8)),
    paste0("MYBL2-low tumors: ", sum(sample_metadata$MYBL2_group == "MYBL2-low")),
    paste0("MYBL2-high tumors: ", sum(sample_metadata$MYBL2_group == "MYBL2-high"))
  ),
  file.path(step12_tabledir, "Pooled_no_batch_correction_global_MYBL2_median_cutoff.txt")
)

write.csv(
  estimate_cpm_sample_by_gene,
  file.path(step12_tabledir, "ESTIMATE_input_nonlog_TMM_CPM_no_batch_correction.csv")
)

saveRDS(
  list(
    expression_log2_TMM_CPM_no_batch_correction = expr_sample_by_gene,
    ESTIMATE_input_nonlog_TMM_CPM_no_batch_correction = estimate_cpm_sample_by_gene,
    sample_metadata = sample_metadata,
    global_MYBL2_median = global_mybl2_median,
    input_TMM_CPM_file = discovery_file
  ),
  file.path(step12_rdsdir, "Pooled_tumor_no_batch_correction_expression_objects.rds")
)

# STEP 12.4. ssGSEA OF 28 IMMUNE CELL TYPES

official_immune28 <- read_local_charoentong_immune28()
immune28_sets <- official_immune28$gene_sets
immune28_source <- CHAROENTONG_REFERENCE_CITATION
immune28_source_names <- official_immune28$source_names

official_gmt_file <- file.path(
  step12_tabledir,
  "Charoentong_2017_28_immune_cell_metagenes_official_source.gmt"
)
write_gmt(immune28_sets, official_gmt_file)

local_workbook_copy <- file.path(
  step12_tabledir,
  "Charoentong_2017_Cell_Reports_Table_S6_mmc3.xlsx"
)

if (!identical(
  normalizePath(CHAROENTONG_MMC3_FILE, winslash = "/", mustWork = TRUE),
  normalizePath(local_workbook_copy, winslash = "/", mustWork = FALSE)
)) {
  file.copy(
    from = CHAROENTONG_MMC3_FILE,
    to = local_workbook_copy,
    overwrite = TRUE
  )
}

write.csv(
  official_immune28$reference_table,
  file.path(
    step12_tabledir,
    "Charoentong_2017_28_immune_cell_metagenes_official_source.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    Immune_cell = names(immune28_sets),
    Signature_source = immune28_source,
    Source_signature_name = unname(immune28_source_names[names(immune28_sets)]),
    Official_workbook = basename(official_immune28$source_file),
    Official_sheet = official_immune28$source_sheet,
    Official_header_row = official_immune28$header_row,
    Exported_GMT = basename(official_gmt_file),
    stringsAsFactors = FALSE
  ),
  file.path(step12_tabledir, "Immune28_signature_provenance.csv"),
  row.names = FALSE
)

immune28_sets_used <- lapply(immune28_sets, function(gs) {
  intersect(unique(gs), colnames(expr_sample_by_gene))
})

overlap_table <- data.frame(
  Immune_cell = names(immune28_sets),
  Source_signature_name = unname(immune28_source_names[names(immune28_sets)]),
  Signature_gene_count = vapply(immune28_sets, length, integer(1)),
  Genes_present_in_LSCC = vapply(immune28_sets_used, length, integer(1)),
  stringsAsFactors = FALSE
)

write.csv(
  overlap_table,
  file.path(step12_tabledir, "Immune28_signature_gene_overlap.csv"),
  row.names = FALSE
)

insufficient <- overlap_table$Immune_cell[
  overlap_table$Genes_present_in_LSCC < MIN_SIGNATURE_GENE_OVERLAP
]
if (length(insufficient) > 0) {
  stop(
    "Fewer than ", MIN_SIGNATURE_GENE_OVERLAP,
    " genes were available for these immune-cell signatures:\n",
    paste(insufficient, collapse = ", "),
    "\n\nInspect: ", file.path(step12_tabledir, "Immune28_signature_gene_overlap.csv")
  )
}

expr_gene_by_sample <- t(expr_sample_by_gene)
ssgsea_scores <- run_ssgsea_safe(expr_gene_by_sample, immune28_sets_used)

ssgsea_scores <- ssgsea_scores[immune28_labels, sample_metadata$Sample, drop = FALSE]

write.csv(
  ssgsea_scores,
  file.path(step12_tabledir, "ssGSEA_28_immune_cell_scores_tumor_samples.csv")
)

saveRDS(
  list(
    immune28_sets = immune28_sets,
    immune28_sets_used = immune28_sets_used,
    ssGSEA_scores = ssgsea_scores,
    sample_metadata = sample_metadata,
    source = immune28_source
  ),
  file.path(step12_rdsdir, "ssGSEA28_immune_microenvironment_objects.rds")
)

# STEP 12.5. MYBL2-HIGH VS MYBL2-LOW IMMUNE COMPARISON, BH-FDR, AND BOXPLOTS

immune_long <- as.data.frame(t(ssgsea_scores)) %>%
  tibble::rownames_to_column("Sample") %>%
  tidyr::pivot_longer(
    cols = -Sample,
    names_to = "Immune_cell",
    values_to = "ssGSEA_score"
  ) %>%
  dplyr::left_join(sample_metadata, by = "Sample") %>%
  dplyr::mutate(
    Immune_cell = factor(Immune_cell, levels = immune28_labels),
    MYBL2_group = factor(MYBL2_group, levels = c("MYBL2-low", "MYBL2-high")),
    Dataset = factor(Dataset, levels = c("GSE127165", "GSE142083"))
  )

immune_wilcox <- immune_long %>%
  dplyr::group_by(Immune_cell) %>%
  dplyr::summarise(
    MYBL2_low_median = median(ssGSEA_score[MYBL2_group == "MYBL2-low"], na.rm = TRUE),
    MYBL2_high_median = median(ssGSEA_score[MYBL2_group == "MYBL2-high"], na.rm = TRUE),
    Median_difference_high_minus_low = MYBL2_high_median - MYBL2_low_median,
    Wilcoxon_p = safe_wilcox(ssGSEA_score, MYBL2_group),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    BH_FDR = p.adjust(Wilcoxon_p, method = "BH"),
    Significant_FDR = is.finite(BH_FDR) & BH_FDR < FDR_CUTOFF,
    Direction = dplyr::case_when(
      Median_difference_high_minus_low > 0 ~ "Higher in MYBL2-high",
      Median_difference_high_minus_low < 0 ~ "Higher in MYBL2-low",
      TRUE ~ "No median difference"
    )
  ) %>%
  dplyr::arrange(BH_FDR, dplyr::desc(abs(Median_difference_high_minus_low)))

write.csv(
  immune_wilcox,
  file.path(step12_tabledir, "ssGSEA28_MYBL2_high_low_Wilcoxon_BH_FDR_results.csv"),
  row.names = FALSE
)

write.csv(
  immune_wilcox %>% dplyr::filter(Significant_FDR),
  file.path(step12_tabledir, "ssGSEA28_MYBL2_high_low_FDR_significant_results.csv"),
  row.names = FALSE
)

immune_wilcox_ranked <- immune_wilcox %>%
  dplyr::mutate(
    BH_FDR_sort = dplyr::if_else(is.finite(BH_FDR), BH_FDR, Inf)
  ) %>%
  dplyr::arrange(
    BH_FDR_sort,
    dplyr::desc(abs(Median_difference_high_minus_low))
  )

significant_boxplot_cells <- immune_wilcox_ranked %>%
  dplyr::filter(Significant_FDR) %>%
  dplyr::pull(Immune_cell) %>%
  as.character()

if (length(significant_boxplot_cells) >= PANEL_B_N_SIGNATURES) {
  boxplot_cells <- head(significant_boxplot_cells, PANEL_B_N_SIGNATURES)

  boxplot_panel_note <- paste0(
    "Final Figure 6 panel b displays the ",
    PANEL_B_N_SIGNATURES,
    " most significant immune-cell signatures (BH-FDR < 0.05)."
  )
} else {
  boxplot_cells <- immune_wilcox_ranked %>%
    dplyr::slice_head(n = PANEL_B_N_SIGNATURES) %>%
    dplyr::pull(Immune_cell) %>%
    as.character()

  boxplot_panel_note <- paste0(
    "Final Figure 6 panel b displays the top ",
    PANEL_B_N_SIGNATURES,
    " BH-FDR-ranked immune-cell signatures. Only ",
    length(significant_boxplot_cells),
    " met BH-FDR < 0.05 and therefore receive significance stars."
  )
}

writeLines(
  boxplot_panel_note,
  file.path(step12_tabledir, "Figure_06B_boxplot_panel_note.txt")
)

box_df <- immune_long %>%
  dplyr::filter(as.character(Immune_cell) %in% boxplot_cells) %>%
  dplyr::mutate(
    Immune_cell_label = factor(
      as.character(Immune_cell),
      levels = boxplot_cells
    )
  )

box_stats <- box_df %>%
  dplyr::group_by(Immune_cell_label) %>%
  dplyr::summarise(
    y_max = max(ssGSEA_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    immune_wilcox %>%
      dplyr::mutate(
        Immune_cell_label = as.character(Immune_cell)
      ) %>%
      dplyr::select(Immune_cell_label, BH_FDR, Significant_FDR),
    by = "Immune_cell_label"
  ) %>%
  dplyr::mutate(
    significance = vapply(BH_FDR, fdr_star, character(1))
  )

box_y_min <- min(box_df$ssGSEA_score, na.rm = TRUE)
box_y_max <- max(box_df$ssGSEA_score, na.rm = TRUE)
box_y_range <- max(box_y_max - box_y_min, 0.05)

box_stats <- box_stats %>%
  dplyr::mutate(
    y_text = y_max + 0.055 * box_y_range
  )

p_immune_box <- ggplot(
  box_df,
  aes(x = Immune_cell_label, y = ssGSEA_score, fill = MYBL2_group)
) +
  geom_boxplot(
    position = position_dodge(width = PANEL_B_DODGE_WIDTH),
    width = PANEL_B_BOX_WIDTH,
    outlier.shape = NA,
    alpha = 0.82,
    linewidth = 0.40
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = PANEL_B_DODGE_WIDTH,
      seed = 123
    ),
    size = PANEL_B_POINT_SIZE,
    alpha = 0.58,
    color = "black"
  ) +
  geom_text(
    data = dplyr::filter(box_stats, significance != ""),
    aes(x = Immune_cell_label, y = y_text, label = significance),
    inherit.aes = FALSE,
    size = 3.20,
    vjust = 0,
    color = "black"
  ) +
  scale_x_discrete(labels = full_panel_b_label) +
  scale_fill_manual(
    values = c("MYBL2-low" = COL_MYBL2_LOW, "MYBL2-high" = COL_MYBL2_HIGH)
  ) +
  labs(
    x = NULL,
    y = "ssGSEA enrichment score",
    fill = "MYBL2 group",
    title = NULL
  ) +
  make_clean_theme(base_size = BASE_TEXT_PT, legend_position = "right") +
  theme(
    axis.text.x = element_text(
      family = FONT_FAMILY,
      size = AXIS_TEXT_PT,
      angle = PANEL_B_X_TEXT_ANGLE,
      hjust = 1,
      vjust = 1,
      lineheight = 0.92,
      color = "black"
    ),
    panel.grid.major.x = element_blank()
  ) +
  coord_cartesian(
    ylim = c(
      box_y_min - 0.04 * box_y_range,
      box_y_max + 0.18 * box_y_range
    ),
    clip = "off"
  )

fig_06b_file <- file.path(
  step12_figdir,
  "Figure_06B_Top15_BH_FDR_ssGSEA_immune_cell_boxplot.png"
)

save_plot(
  p_immune_box,
  basename(fig_06b_file),
  width = 16.8,
  height = 10.9
)

# STEP 12.6. MYBL2–IMMUNE-CELL SPEARMAN CORRELATION WITH BH-FDR

immune_cor <- immune_long %>%
  dplyr::group_by(Immune_cell) %>%
  dplyr::summarise(
    tmp = list(safe_spearman(MYBL2_expression, ssGSEA_score)),
    .groups = "drop"
  ) %>%
  tidyr::unnest(tmp) %>%
  dplyr::mutate(
    BH_FDR = p.adjust(p_value, method = "BH"),
    Significant_FDR_rule = is.finite(rho) &
      is.finite(BH_FDR) &
      abs(rho) > CORRELATION_RHO_CUTOFF &
      BH_FDR < FDR_CUTOFF,
    Direction = dplyr::case_when(
      rho > 0 ~ "Positive",
      rho < 0 ~ "Negative",
      TRUE ~ "Zero"
    )
  ) %>%
  dplyr::arrange(dplyr::desc(abs(rho)), BH_FDR)

write.csv(
  immune_cor,
  file.path(step12_tabledir, "MYBL2_ssGSEA28_Spearman_BH_FDR_correlation_results.csv"),
  row.names = FALSE
)

write.csv(
  immune_cor %>% dplyr::filter(Significant_FDR_rule),
  file.path(step12_tabledir, "MYBL2_ssGSEA28_Spearman_FDR_significant_correlations.csv"),
  row.names = FALSE
)

immune_cor_plot <- immune_cor %>%
  dplyr::filter(is.finite(rho)) %>%
  dplyr::mutate(
    Immune_cell_label = factor(
      as.character(Immune_cell),
      levels = rev(immune28_labels)
    ),
    Correlation_class = dplyr::case_when(
      Significant_FDR_rule & rho > 0 ~ "Significant positive",
      Significant_FDR_rule & rho < 0 ~ "Significant negative",
      TRUE ~ "Not significant"
    ),
    rho_label = ifelse(
      Significant_FDR_rule,
      sprintf("%.2f", rho),
      ""
    )
  )

p_immune_cor <- ggplot(
  immune_cor_plot,
  aes(x = "MYBL2", y = Immune_cell_label)
) +
  geom_point(
    aes(size = abs(rho), color = Correlation_class),
    alpha = 0.95
  ) +
  geom_text(
    data = dplyr::filter(immune_cor_plot, Significant_FDR_rule),
    aes(label = rho_label),
    size = 3.10,
    color = "black",
    fontface = "plain"
  ) +
  scale_color_manual(
    values = c(
      "Significant positive" = COL_TUMOR,
      "Significant negative" = COL_NORMAL,
      "Not significant" = "grey78"
    ),
    breaks = c(
      "Significant positive",
      "Significant negative",
      "Not significant"
    ),
    name = "Correlation result"
  ) +
  scale_size_continuous(
    range = c(2.7, 9.5),
    breaks = c(0.10, 0.30, 0.50),
    name = "|Spearman rho|"
  ) +
  labs(x = NULL, y = NULL, title = NULL) +
  make_clean_theme(base_size = BASE_TEXT_PT, legend_position = "right") +
  theme(
    panel.grid.major.x = element_blank(),
    legend.box = "vertical"
  )

fig_06c_file <- file.path(
  step12_figdir,
  "Figure_06C_MYBL2_ssGSEA28_BH_FDR_correlation_bubbleplot.png"
)

save_plot(
  p_immune_cor,
  basename(fig_06c_file),
  width = 9.4,
  height = 16.9
)

# STEP 12.7. ESTIMATE SCORES, GROUP COMPARISONS, AND CORRELATIONS

estimate_gene_by_sample <- t(estimate_cpm_sample_by_gene)
message("Running ESTIMATE on non-log TMM-CPM values without batch-effect correction.")
estimate_df <- run_estimate_safe(estimate_gene_by_sample)
colnames(estimate_df) <- gsub(" ", "_", colnames(estimate_df), fixed = TRUE)

estimate_df <- estimate_df %>%
  dplyr::rename(Sample = ID) %>%
  dplyr::filter(Sample %in% sample_metadata$Sample) %>%
  dplyr::left_join(sample_metadata, by = "Sample")

score_name_candidates <- list(
  Stromal_score = c("StromalScore_estimate", "StromalScore", "Stromal_score"),
  Immune_score = c("ImmuneScore_estimate", "ImmuneScore", "Immune_score"),
  ESTIMATE_score = c("ESTIMATEScore_estimate", "ESTIMATEScore", "EstimateScore", "ESTIMATE_score")
)

estimate_score_columns <- vapply(score_name_candidates, function(candidates) {
  hit <- intersect(candidates, colnames(estimate_df))
  if (length(hit) == 0) NA_character_ else hit[1]
}, character(1))

if (anyNA(estimate_score_columns)) {
  stop(
    "Could not identify all three ESTIMATE score columns. Available columns:\n",
    paste(colnames(estimate_df), collapse = ", ")
  )
}

estimate_scores <- estimate_df %>%
  dplyr::select(
    Sample,
    Dataset,
    Batch,
    MYBL2_expression,
    MYBL2_group,
    dplyr::all_of(unname(estimate_score_columns))
  )

colnames(estimate_scores)[
  match(unname(estimate_score_columns), colnames(estimate_scores))
] <- names(estimate_score_columns)

write.csv(
  estimate_scores,
  file.path(step12_tabledir, "ESTIMATE_scores_nonlog_TMM_CPM_no_batch_correction_tumor_samples.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "ESTIMATE input information",
    "Input scale: original non-log TMM-normalized CPM values.",
    "Source: train_discovery_merged_CPM_nonnegative.csv from the bulk RNA-seq/ML pipeline.",
    "No batch-effect correction was applied."
  ),
  file.path(step12_tabledir, "ESTIMATE_input_scale_no_batch_correction.txt")
)

estimate_long <- estimate_scores %>%
  tidyr::pivot_longer(
    cols = c("Stromal_score", "Immune_score", "ESTIMATE_score"),
    names_to = "TME_score",
    values_to = "Score"
  ) %>%
  dplyr::mutate(
    TME_score = factor(
      TME_score,
      levels = c("Stromal_score", "Immune_score", "ESTIMATE_score"),
      labels = c("Stromal score", "Immune score", "ESTIMATE score")
    )
  )

estimate_wilcox <- estimate_long %>%
  dplyr::group_by(TME_score) %>%
  dplyr::summarise(
    MYBL2_low_median = median(Score[MYBL2_group == "MYBL2-low"], na.rm = TRUE),
    MYBL2_high_median = median(Score[MYBL2_group == "MYBL2-high"], na.rm = TRUE),
    Median_difference_high_minus_low = MYBL2_high_median - MYBL2_low_median,
    Wilcoxon_p = safe_wilcox(Score, MYBL2_group),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    BH_FDR = p.adjust(Wilcoxon_p, method = "BH"),
    Significant_FDR = is.finite(BH_FDR) & BH_FDR < FDR_CUTOFF
  )

write.csv(
  estimate_wilcox,
  file.path(step12_tabledir, "ESTIMATE_scores_MYBL2_high_low_Wilcoxon_BH_FDR_results.csv"),
  row.names = FALSE
)

estimate_panel_stats <- estimate_long %>%
  dplyr::group_by(TME_score) %>%
  dplyr::summarise(
    y_min = min(Score, na.rm = TRUE),
    y_max = max(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    estimate_wilcox %>%
      dplyr::select(TME_score, Wilcoxon_p, BH_FDR, Significant_FDR),
    by = "TME_score"
  ) %>%
  dplyr::mutate(
    y_range = pmax(y_max - y_min, 1),
    y_line = y_max + 0.07 * y_range,
    y_text = y_max + 0.13 * y_range,
    fdr_label = vapply(BH_FDR, format_fdr, character(1))
  )

p_estimate_group <- ggplot(
  estimate_long,
  aes(x = MYBL2_group, y = Score, fill = MYBL2_group)
) +
  geom_violin(trim = FALSE, alpha = 0.82, color = NA) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.55
  ) +
  geom_jitter(width = 0.11, size = 1.10, alpha = 0.70, color = "black") +
  geom_segment(
    data = estimate_panel_stats,
    aes(x = 1, xend = 2, y = y_line, yend = y_line),
    inherit.aes = FALSE,
    linewidth = 0.38,
    color = "black"
  ) +
  geom_text(
    data = estimate_panel_stats,
    aes(x = 1.5, y = y_text, label = fdr_label),
    inherit.aes = FALSE,
    size = 3.0,
    color = "black"
  ) +
  facet_wrap(~ TME_score, scales = "free_y", nrow = 1) +
  scale_fill_manual(
    values = c("MYBL2-low" = COL_MYBL2_LOW, "MYBL2-high" = COL_MYBL2_HIGH)
  ) +
  labs(x = NULL, y = "ESTIMATE-derived score", title = NULL) +
  make_clean_theme(base_size = BASE_TEXT_PT, legend_position = "none") +
  coord_cartesian(clip = "off")

fig_06d_file <- file.path(
  step12_figdir,
  "Figure_06D_ESTIMATE_score_MYBL2_high_low_BH_FDR_comparison.png"
)

save_plot(
  p_estimate_group,
  basename(fig_06d_file),
  width = 13.4,
  height = 6.9
)

estimate_cor <- estimate_long %>%
  dplyr::group_by(TME_score) %>%
  dplyr::summarise(
    tmp = list(safe_spearman(MYBL2_expression, Score)),
    .groups = "drop"
  ) %>%
  tidyr::unnest(tmp) %>%
  dplyr::mutate(
    BH_FDR = p.adjust(p_value, method = "BH"),
    Significant_FDR_rule = is.finite(rho) &
      is.finite(BH_FDR) &
      abs(rho) > CORRELATION_RHO_CUTOFF &
      BH_FDR < FDR_CUTOFF,
    Label = paste0(
      "rho = ", sprintf("%.2f", rho),
      "\n", vapply(BH_FDR, format_fdr, character(1))
    )
  )

write.csv(
  estimate_cor,
  file.path(step12_tabledir, "MYBL2_ESTIMATE_score_Spearman_BH_FDR_correlation_results.csv"),
  row.names = FALSE
)

p_estimate_cor <- ggplot(
  estimate_long,
  aes(x = MYBL2_expression, y = Score, color = MYBL2_group, shape = Dataset)
) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "grey60",
    alpha = 0.16,
    linewidth = 0.65,
    show.legend = FALSE
  ) +
  geom_point(size = 2.45, alpha = 0.82, stroke = 0.2) +
  geom_label(
    data = estimate_cor,
    aes(x = -Inf, y = Inf, label = Label),
    inherit.aes = FALSE,
    hjust = -0.06,
    vjust = 1.15,
    size = 3.40,
    linewidth = 0.24,
    fill = "white",
    alpha = 0.92,
    color = "black"
  ) +
  facet_wrap(~ TME_score, scales = "free_y", nrow = 1) +
  scale_color_manual(
    values = c("MYBL2-low" = COL_MYBL2_LOW, "MYBL2-high" = COL_MYBL2_HIGH)
  ) +
  scale_shape_manual(
    values = c("GSE127165" = 16, "GSE142083" = 17)
  ) +
  labs(
    x = "MYBL2 expression",
    y = "ESTIMATE-derived score",
    color = "MYBL2 group",
    shape = "Dataset",
    title = NULL
  ) +
  make_clean_theme(base_size = BASE_TEXT_PT, legend_position = "bottom") +
  theme(
    legend.box = "vertical"
  ) +
  guides(
    color = guide_legend(order = 1, nrow = 1),
    shape = guide_legend(order = 2, nrow = 1)
  )

fig_06e_file <- file.path(
  step12_figdir,
  "Figure_06E_MYBL2_ESTIMATE_score_BH_FDR_correlations.png"
)

save_plot(
  p_estimate_cor,
  basename(fig_06e_file),
  width = 14.0,
  height = 7.2
)

# STEP 12.8. FINAL COMBINED MANUSCRIPT FIGURE 6

COL_RF_ROC  <- "#274C77"
COL_SVM_ROC <- "#D06B43"

ml_dir <- ml_input_dir
immune_dir <- immune_output_dir
FIGURE_DIR <- immune_output_dir

train_ml_file <- file.path(ml_dir, "train_log2CPM.csv")
external_ml_file <- file.path(ml_dir, "external_GSE130605_log2CPM.csv")

if (!file.exists(train_ml_file) || !file.exists(external_ml_file)) {
  stop(
    "Panel a cannot be rebuilt because one or both required ML matrices are missing:\n",
    train_ml_file, "\n", external_ml_file
  )
}

train_ml <- as.data.frame(data.table::fread(train_ml_file), check.names = FALSE)
external_ml <- as.data.frame(data.table::fread(external_ml_file), check.names = FALSE)

if (!all(c("group", TARGET_GENE) %in% colnames(train_ml))) {
  stop("The training ML matrix must contain 'group' and ", TARGET_GENE, ".")
}
if (!all(c("group", TARGET_GENE) %in% colnames(external_ml))) {
  stop("The external ML matrix must contain 'group' and ", TARGET_GENE, ".")
}

make_binary_ml_group <- function(x) {
  z <- trimws(tolower(as.character(x)))
  out <- ifelse(
    z %in% c("2", "tumor", "lscc", "cancer", "case"), "Tumor",
    ifelse(z %in% c("1", "normal", "non", "control", "margin"), "Normal", NA_character_)
  )
  if (anyNA(out) || length(unique(out)) != 2) {
    stop("The ML group column could not be mapped unambiguously to Normal and Tumor.")
  }
  factor(out, levels = c("Normal", "Tumor"))
}

train_y <- make_binary_ml_group(train_ml$group)
external_y <- make_binary_ml_group(external_ml$group)

train_x_raw <- suppressWarnings(as.numeric(train_ml[[TARGET_GENE]]))
external_x_raw <- suppressWarnings(as.numeric(external_ml[[TARGET_GENE]]))

if (any(!is.finite(train_x_raw)) || any(!is.finite(external_x_raw))) {
  stop("Non-finite MYBL2 values were detected in the ML matrices.")
}

train_center <- mean(train_x_raw)
train_scale <- stats::sd(train_x_raw)
if (!is.finite(train_scale) || train_scale == 0) {
  stop("MYBL2 has zero or non-finite variance in the training dataset.")
}

train_x <- data.frame(MYBL2 = (train_x_raw - train_center) / train_scale)
external_x <- data.frame(MYBL2 = (external_x_raw - train_center) / train_scale)

set.seed(123)
rf_fit <- randomForest::randomForest(
  x = train_x,
  y = train_y,
  ntree = 500
)

svm_fit <- e1071::svm(
  x = train_x,
  y = train_y,
  kernel = "linear",
  probability = TRUE,
  scale = FALSE
)

rf_prob <- predict(rf_fit, newdata = external_x, type = "prob")[, "Tumor"]
svm_pred <- predict(svm_fit, newdata = external_x, probability = TRUE)
svm_prob <- attr(svm_pred, "probabilities")[, "Tumor"]

roc_rf <- pROC::roc(
  response = external_y,
  predictor = rf_prob,
  levels = c("Normal", "Tumor"),
  direction = "<",
  quiet = TRUE
)
roc_svm <- pROC::roc(
  response = external_y,
  predictor = svm_prob,
  levels = c("Normal", "Tumor"),
  direction = "<",
  quiet = TRUE
)

auc_rf <- as.numeric(pROC::auc(roc_rf))
auc_svm <- as.numeric(pROC::auc(roc_svm))

roc_df <- dplyr::bind_rows(
  data.frame(
    FPR = 1 - roc_rf$specificities,
    Sensitivity = roc_rf$sensitivities,
    Model = sprintf("RF (AUC = %.2f%%)", 100 * auc_rf)
  ),
  data.frame(
    FPR = 1 - roc_svm$specificities,
    Sensitivity = roc_svm$sensitivities,
    Model = sprintf("SVM (AUC = %.2f%%)", 100 * auc_svm)
  )
) %>%
  dplyr::mutate(
    Model = factor(
      Model,
      levels = c(
        sprintf("RF (AUC = %.2f%%)", 100 * auc_rf),
        sprintf("SVM (AUC = %.2f%%)", 100 * auc_svm)
      )
    )
  )

p_external_roc <- ggplot(
  roc_df,
  aes(x = FPR, y = Sensitivity, colour = Model)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey55",
    linewidth = 0.45
  ) +
  geom_step(linewidth = 0.85) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_colour_manual(
    values = setNames(
      c(COL_RF_ROC, COL_SVM_ROC),
      levels(roc_df$Model)
    ),
    name = NULL
  ) +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity",
    title = NULL
  ) +
  make_clean_theme(
    base_size = BASE_TEXT_PT,
    legend_position = c(0.70, 0.19),
    show_grid = TRUE
  ) +
  theme(
    legend.background = element_rect(
      fill = FIG_BACKGROUND,
      colour = "black",
      linewidth = PANEL_BORDER_LWD
    )
  )

fig6a_external_file <- file.path(
  step12_figdir,
  "Figure_06A_External_validation_RF_SVM_ROC_navy_copper_no_title.png"
)

save_plot(
  p_external_roc,
  basename(fig6a_external_file),
  width = 6.0,
  height = 8.0
)

write.csv(
  data.frame(
    Model = c("RF", "SVM"),
    AUROC = c(auc_rf, auc_svm),
    stringsAsFactors = FALSE
  ),
  file.path(step12_tabledir, "Figure_06A_navy_copper_external_ROC_AUCs.csv"),
  row.names = FALSE
)

fig6b_boxplot_file <- fig_06b_file
fig6c_bubble_file <- fig_06c_file
fig6d_est_group_file <- fig_06d_file
fig6e_est_cor_file <- fig_06e_file

required_panel_files <- c(
  fig6a_external_file,
  fig6b_boxplot_file,
  fig6c_bubble_file,
  fig6d_est_group_file,
  fig6e_est_cor_file
)

missing_panel_files <- required_panel_files[!file.exists(required_panel_files)]

if (length(missing_panel_files) > 0) {
  stop(
    "The following panel file(s) are missing:\n",
    paste(missing_panel_files, collapse = "\n")
  )
}

panel_a_width <- 3600
panel_b_width <- 7400
panel_top_height <- 4800

panel_c_width <- 4000
panel_right_width <- 7000
panel_bottom_each_height <- 3600
panel_c_height <- panel_bottom_each_height * 2

label_size <- RASTER_PANEL_TAG_SIZE

img_a <- make_clean_panel_img(
  path = fig6a_external_file,
  label = "a",
  panel_width = panel_a_width,
  panel_height = panel_top_height,
  label_size = label_size
)

img_b <- make_clean_panel_img(
  path = fig6b_boxplot_file,
  label = "b",
  panel_width = panel_b_width,
  panel_height = panel_top_height,
  label_size = label_size
)

img_c <- make_clean_panel_img(
  path = fig6c_bubble_file,
  label = "c",
  panel_width = panel_c_width,
  panel_height = panel_c_height,
  label_size = label_size
)

img_d <- make_clean_panel_img(
  path = fig6d_est_group_file,
  label = "d",
  panel_width = panel_right_width,
  panel_height = panel_bottom_each_height,
  label_size = label_size
)

img_e <- make_clean_panel_img(
  path = fig6e_est_cor_file,
  label = "e",
  panel_width = panel_right_width,
  panel_height = panel_bottom_each_height,
  label_size = label_size
)

row_1 <- magick::image_append(c(img_a, img_b), stack = FALSE)
right_stack <- magick::image_append(c(img_d, img_e), stack = TRUE)
row_2 <- magick::image_append(c(img_c, right_stack), stack = FALSE)

fig6_combined_img <- magick::image_append(
  c(row_1, row_2),
  stack = TRUE
)

fig6_combined_img <- magick::image_background(
  fig6_combined_img,
  FIG_BACKGROUND,
  flatten = TRUE
)

figure_06_file <- file.path(
  FIGURE_DIR,
  "Figure_06_NavyRoseImmune_MYBL2_combined_larger_panels_less_whitespace.png"
)

save_magick_all_formats(
  image_object = fig6_combined_img,
  filename_stem = "Figure_06_NavyRoseImmune_MYBL2_combined_larger_panels_less_whitespace",
  dir_path = FIGURE_DIR
)

cat("\nFIGURE 6 RECOMBINATION FINISHED.\n")
cat("Panel order:\n")
cat("a = External validation ROC (RF + SVM)\n")
cat("b = Top 15 BH-FDR-ranked ssGSEA boxplots (complete signature names)\n")
cat("c = MYBL2–immune correlation bubble plot\n")
cat("d = ESTIMATE group comparison\n")
cat("e = ESTIMATE correlation plots\n\n")
cat("The descriptive ssGSEA heatmap was intentionally omitted.\n")
cat("All combined outputs were saved to:\n")
cat(FIGURE_DIR, "\n")

# STEP 12.9. ANALYSIS SUMMARY AND SESSION INFORMATION

summary_lines <- c(
  "MYBL2-associated immune microenvironment analysis completed.",
  "Final Figure 6 contains five panels; panel b shows exactly 15 BH-FDR-ranked ssGSEA signatures and the descriptive ssGSEA heatmap was intentionally omitted.",
  "Panels b-e use a indigo-blue / rose palette; panel a uses a navy / copper ROC palette without a title banner.",
  paste0("Pooled tumor samples analysed: ", nrow(sample_metadata)),
  paste0("GSE127165 tumor samples: ", sum(sample_metadata$Dataset == "GSE127165")),
  paste0("GSE142083 tumor samples: ", sum(sample_metadata$Dataset == "GSE142083")),
  "Batch correction: none; pooled TMM-CPM values were analysed without correction.",
  "ESTIMATE input: original non-log TMM-CPM values from the pooled discovery matrix.",
  paste0("Global MYBL2 median on log2(TMM-CPM + 1): ", format(global_mybl2_median, digits = 8)),
  paste0("MYBL2-low samples: ", sum(sample_metadata$MYBL2_group == "MYBL2-low")),
  paste0("MYBL2-high samples: ", sum(sample_metadata$MYBL2_group == "MYBL2-high")),
  paste0("28 official Charoentong immune-signature source: ", immune28_source),
  paste0(
    "Panel b signatures displayed: ", length(boxplot_cells),
    " (target = ", PANEL_B_N_SIGNATURES, ")."
  ),
  paste0(
    "Immune signatures with BH-FDR < ", FDR_CUTOFF, ": ",
    sum(immune_wilcox$Significant_FDR, na.rm = TRUE)
  ),
  paste0(
    "MYBL2–immune correlations meeting |rho| > ", CORRELATION_RHO_CUTOFF,
    " and BH-FDR < ", FDR_CUTOFF, ": ",
    sum(immune_cor$Significant_FDR_rule, na.rm = TRUE)
  ),
  paste0(
    "ESTIMATE group comparisons with BH-FDR < ", FDR_CUTOFF, ": ",
    sum(estimate_wilcox$Significant_FDR, na.rm = TRUE)
  ),
  paste0(
    "ESTIMATE correlations meeting |rho| > ", CORRELATION_RHO_CUTOFF,
    " and BH-FDR < ", FDR_CUTOFF, ": ",
    sum(estimate_cor$Significant_FDR_rule, na.rm = TRUE)
  ),
  "",
  "Interpretation: ssGSEA values are relative immune-signature enrichment scores,",
  "not direct immune-cell fractions. ESTIMATE values are computational TME scores.",
  "Associations do not establish a causal effect of MYBL2 on immune regulation."
)

writeLines(
  summary_lines,
  file.path(step12_dir, "Figure_06_analysis_summary.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(step12_dir, "Figure_06_sessionInfo.txt")
)

cat("\n============================================================\n")
cat("MYBL2 IMMUNE MICROENVIRONMENT ANALYSIS FINISHED SUCCESSFULLY\n")
cat("============================================================\n")
cat("Pooled input CPM file:\n", discovery_file, "\n", sep = "")
cat("Output directory (all files):\n", immune_output_dir, "\n", sep = "")
cat("Combined Figure 6 (navy/copper external ROC as panel a):\n", figure_06_file, "\n", sep = "")
cat("============================================================\n")

# MODULE 03 — Functional programme analysis and Figure 7 outputs

# STEP 01.1 — SETUP, PATHS, AND MANUSCRIPT FIGURE SETTINGS

RESET_WORKSPACE_AT_START <- TRUE

if (RESET_WORKSPACE_AT_START) {
  rm(list = ls())
}

gc()

completed_steps <- character(0)

mark_step_done <- function(step_id) {
  completed_steps <<- unique(c(completed_steps, step_id))
  message("\n✓ STEP ", step_id, " completed. Move to the next numbered step.\n")
}

require_previous_step <- function(required_step, current_step) {
  if (!(required_step %in% completed_steps)) {
    stop(
      "STEP ", current_step,
      " cannot run yet. First run the complete STEP ",
      required_step,
      " section."
    )
  }
}

set.seed(123)

options(
  stringsAsFactors = FALSE,
  scipen = 100,
  timeout = 7200
)

TARGET_GENE <- "MYBL2"

FONT_FAMILY <- "Arial"
FIG_DPI <- 600
FIG_BACKGROUND <- "white"

if (.Platform$OS.type == "windows") {
  try(
    grDevices::windowsFonts(
      Arial = grDevices::windowsFont("Arial")
    ),
    silent = TRUE
  )
}

FIG_SINGLE_W <- 3.50
FIG_SINGLE_H <- 4.20

FIG_DOUBLE_W <- 7.20
FIG_DOUBLE_H <- 5.40

FIGURE7_COMBINED_W <- 16.50
FIGURE7_COMBINED_H <- 12.20

FIGURE7_TOP_ROW_HEIGHT <- 1.15
FIGURE7_BOTTOM_ROW_HEIGHT <- 0.82

BASE_TEXT_PT <- 9.5
AXIS_TEXT_PT <- 8.5
AXIS_TITLE_PT <- 10
LEGEND_TEXT_PT <- 8.5
LEGEND_TITLE_PT <- 9
PLOT_TITLE_PT <- 10
PANEL_TAG_PT <- 14

PANEL_BORDER_LWD <- 0.45
PANEL_SECTION_BORDER_LWD <- 0.45
AXIS_LWD <- 0.40
TICK_LWD <- 0.40
GEOM_LWD <- 0.55
GRID_LWD <- 0.30
TICK_LENGTH_PT <- 2
PANEL_SPACING_PT <- 3
PANEL_TAG_POSITION <- c(0.045, 0.955)

VERTICAL_DOTPLOT_X_LIMITS <- c(0.55, 1.45)
VERTICAL_DOTPLOT_POINT_RANGE <- c(3.0, 14.0)

PANEL_B_TAG_POSITION <- c(0.018, 0.988)
PANEL_B_TAG_MARGIN <- ggplot2::margin(
  t = 13,
  r = 4,
  b = 3,
  l = 19,
  unit = "pt"
)

PANEL_A_TAG_POSITION <- c(0.018, 0.988)
PANEL_A_TAG_MARGIN <- ggplot2::margin(
  t = 13,
  r = 4,
  b = 3,
  l = 19,
  unit = "pt"
)

HEATMAP_FONT_PT <- 8
HEATMAP_ROW_FONT_PT <- AXIS_TEXT_PT
HEATMAP_COL_FONT_PT <- AXIS_TEXT_PT
HEATMAP_LEGEND_FONT_PT <- LEGEND_TEXT_PT

COL_NORMAL <- "#7470B2"
COL_TUMOR <- "#D81B60"

COL_MYBL2_LOW <- COL_NORMAL
COL_MYBL2_HIGH <- COL_TUMOR
COL_UP <- COL_TUMOR
COL_DOWN <- COL_NORMAL
COL_NS <- "#B8B8B8"

MANUSCRIPT_MARGIN <- ggplot2::margin(
  t = 3,
  r = 4,
  b = 3,
  l = 4,
  unit = "pt"
)

LEGEND_BOX_MARGIN <- ggplot2::margin(
  t = 0,
  r = 0,
  b = 0,
  l = 0,
  unit = "pt"
)

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
        size = PLOT_TITLE_PT,
        hjust = 0.5,
        vjust = 0.5,
        color = "black"
      ),
      plot.title.position = "panel",

      plot.tag = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PANEL_TAG_PT,
        hjust = 0,
        vjust = 1,
        color = "black"
      ),
      plot.tag.position = PANEL_TAG_POSITION,

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
        linewidth = TICK_LWD
      ),
      axis.ticks.length = grid::unit(TICK_LENGTH_PT, "pt"),

      legend.position = legend_position,
      legend.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = LEGEND_TITLE_PT,
        color = "black"
      ),
      legend.text = ggplot2::element_text(
        family = FONT_FAMILY,
        size = LEGEND_TEXT_PT,
        color = "black"
      ),
      legend.background = ggplot2::element_rect(
        fill = FIG_BACKGROUND,
        color = NA
      ),
      legend.key = ggplot2::element_rect(
        fill = FIG_BACKGROUND,
        color = NA
      ),
      legend.key.width = grid::unit(0.38, "cm"),
      legend.key.height = grid::unit(0.38, "cm"),
      legend.box.margin = LEGEND_BOX_MARGIN,
      legend.box.spacing = grid::unit(PANEL_SPACING_PT, "pt"),

      panel.background = ggplot2::element_rect(
        fill = FIG_BACKGROUND,
        color = NA
      ),
      plot.background = ggplot2::element_rect(
        fill = FIG_BACKGROUND,
        color = "black",
        linewidth = PANEL_SECTION_BORDER_LWD
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

      panel.spacing = grid::unit(PANEL_SPACING_PT, "pt"),
      plot.margin = MANUSCRIPT_MARGIN
    )
}

TSNE_PANEL_WIDTH <- 1.20
MYBL2_PANEL_WIDTH <- 0.80

PANEL_C_LEFT_SPACER_WIDTH <- 0.35
PANEL_C_HEATMAP_WIDTH <- 1.30
PANEL_C_RIGHT_SPACER_WIDTH <- 0.35

PANEL_C_TAG_POSITION <- c(0.018, 0.988)
PANEL_C_TAG_MARGIN <- ggplot2::margin(
  t = 13,
  r = 4,
  b = 3,
  l = 19,
  unit = "pt"
)

RESOLUTION_MALIGNANT_FROM_PREVIOUS_PIPELINE <- 0.5

FIXED_CLUSTER_ORDER <- c(
  "Cluster 0", "Cluster 1", "Cluster 2", "Cluster 3",
  "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 8",
  "Cluster 10", "Cluster 11"
)

PANEL_C_CLUSTER_ORDER <- c("Cluster 0", "Cluster 8")

PROGRAMME_DISPLAY_ORDER <- c(
  "DNA repair",
  "E2F targets",
  "EMT",
  "Epithelial differentiation",
  "G2/M checkpoint",
  "Hypoxia",
  "MYC targets",
  "Senescence-like",
  "Stress response"
)

SC_PROJECT_DIR <- "D:/LSCC/ScRNAseq_Results/GSE206332"
BASE_RESULTS_DIR <- "D:/LSCC/Results_LSCC"

OUTPUT_DIR <- file.path(
  BASE_RESULTS_DIR,
  "Functional"
)

PROGRAMME_INTERPRETATION_NOTE <- file.path(
  OUTPUT_DIR,
  "Functional_programme_interpretation_note.txt"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

CLEAN_OLD_FUNCTIONAL_FIGURES <- TRUE

if (CLEAN_OLD_FUNCTIONAL_FIGURES) {
  old_functional_figure_files <- list.files(
    OUTPUT_DIR,
    pattern = "^Figure_7_HighCNV_subtypes_MYBL2_functional_programmes_FixedAllPanelLabels[.](png|tiff|pdf)$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(old_functional_figure_files) > 0) {
    unlink(old_functional_figure_files, force = TRUE)
  }
}

HIGHCNV_OBJECT_RDS_OVERRIDE <- NA_character_
HIGHCNV_DE_OBJECT_RDS_OVERRIDE <- NA_character_

mark_step_done("01.1")

# STEP 01.2 — DETECT COMPLETED PRIOR-PIPELINE RDS FILES

require_previous_step("01.1", "01.2")

candidate_rds_dirs <- c(
  file.path(
    SC_PROJECT_DIR,
    "Results_Modular",
    "Step_20_Final_Exports_and_Objects",
    "rds"
  ),
  file.path(
    SC_PROJECT_DIR,
    "Results",
    "rds"
  )
)

detect_file <- function(override, candidates, label) {
  if (!is.na(override) && nzchar(override)) {
    if (!file.exists(override)) {
      stop(
        label,
        " override does not exist:\n",
        override
      )
    }

    return(
      normalizePath(
        override,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  found <- candidates[file.exists(candidates)]

  if (length(found) == 0) {
    stop(
      label,
      " was not found. Checked:\n",
      paste(candidates, collapse = "\n")
    )
  }

  normalizePath(
    found[1],
    winslash = "/",
    mustWork = TRUE
  )
}

HIGHCNV_OBJECT_RDS <- detect_file(
  HIGHCNV_OBJECT_RDS_OVERRIDE,
  file.path(
    candidate_rds_dirs,
    "HIGH_CNV_MALIGNANT_object_after_lowCNV_removal.rds"
  ),
  "HIGH_CNV_MALIGNANT_object_after_lowCNV_removal.rds"
)

HIGHCNV_DE_OBJECT_RDS <- detect_file(
  HIGHCNV_DE_OBJECT_RDS_OVERRIDE,
  file.path(
    candidate_rds_dirs,
    "HIGH_CNV_MALIGNANT_object_DE_fixed.rds"
  ),
  "HIGH_CNV_MALIGNANT_object_DE_fixed.rds"
)

writeLines(
  c(
    "Input provenance",
    paste0(
      "Existing malignant-clustering resolution retained: ",
      RESOLUTION_MALIGNANT_FROM_PREVIOUS_PIPELINE
    ),
    paste0("High-CNV object: ", HIGHCNV_OBJECT_RDS),
    paste0("High-CNV DE object: ", HIGHCNV_DE_OBJECT_RDS),
    paste0("Output directory: ", OUTPUT_DIR),
    "No upstream single-cell processing was rerun."
  ),
  file.path(OUTPUT_DIR, "Functional_programme_input_provenance.txt")
)

mark_step_done("01.2")

# STEP 01.3 — PACKAGE CHECK AND LOADING

require_previous_step("01.2", "01.3")

REQUIRED_PACKAGES <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "dplyr",
  "tidyr",
  "ggplot2",
  "patchwork",
  "tibble",
  "msigdbr"
)

MISSING_PACKAGES <- REQUIRED_PACKAGES[
  !vapply(
    REQUIRED_PACKAGES,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(MISSING_PACKAGES) > 0) {
  stop(
    "The following required packages are not installed: ",
    paste(MISSING_PACKAGES, collapse = ", "),
    ".\nPlease install them manually, restart RStudio, and rerun the script."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(msigdbr)
})

writeLines(
  c(
    "Required packages were found and loaded successfully.",
    paste0("R: ", R.version.string),
    paste0(
      "Packages: ",
      paste(
        paste0(
          REQUIRED_PACKAGES,
          " (",
          vapply(
            REQUIRED_PACKAGES,
            function(pkg) as.character(utils::packageVersion(pkg)),
            character(1)
          ),
          ")"
        ),
        collapse = "; "
      )
    )
  ),
  file.path(OUTPUT_DIR, "Functional_programme_package_versions.txt")
)

mark_step_done("01.3")

# STEP 01.4 — HELPER FUNCTIONS

require_previous_step("01.3", "01.4")

save_plot_all_formats <- function(plot_obj,
                                  filename_stem,
                                  dir_path = OUTPUT_DIR,
                                  width = FIG_DOUBLE_W,
                                  height = FIG_DOUBLE_H) {
  is_full_path <- grepl("^[A-Za-z]:", filename_stem) ||
    startsWith(filename_stem, "/") ||
    startsWith(filename_stem, "~")

  if (is_full_path) {
    dir_path <- dirname(filename_stem)
    filename_stem <- basename(filename_stem)
  }

  filename_stem <- sub("\\.(png|tiff|pdf)$", "", filename_stem, ignore.case = TRUE)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(dir_path)) {
    stop("Could not create the output directory:\n", dir_path)
  }

  png_file <- file.path(dir_path, paste0(filename_stem, ".png"))
  tiff_file <- file.path(dir_path, paste0(filename_stem, ".tiff"))
  pdf_file <- file.path(dir_path, paste0(filename_stem, ".pdf"))

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
    bg = FIG_BACKGROUND,
    compression = "lzw",
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
    list(
      PNG = png_file,
      TIFF = tiff_file,
      PDF = pdf_file
    )
  )
}

meta_with_cell_id <- function(meta) {
  x <- as.data.frame(
    meta,
    stringsAsFactors = FALSE
  )

  if ("Cell" %in% colnames(x)) {
    x$Cell <- NULL
  }

  tibble::rownames_to_column(
    x,
    var = "Cell"
  )
}

get_existing_layer <- function(object, assay, layer = "data") {
  if (!(assay %in% names(object@assays))) {
    stop(
      "Assay not found: ",
      assay
    )
  }

  DefaultAssay(object) <- assay

  assay_obj <- object[[assay]]

  layer_names <- tryCatch(
    SeuratObject::Layers(assay_obj),
    error = function(e) character(0)
  )

  matching_layers <- layer_names[
    grepl(
      paste0("^", layer, "(\\.|$)"),
      layer_names
    )
  ]

  if (length(matching_layers) == 0) {
    legacy_layer <- tryCatch(
      Seurat::GetAssayData(
        object,
        assay = assay,
        slot = layer
      ),
      error = function(e) NULL
    )

    if (!is.null(legacy_layer) && ncol(legacy_layer) > 0) {
      return(as(legacy_layer, "dgCMatrix"))
    }

    stop(
      "No usable ",
      layer,
      " layer was found in assay ",
      assay
    )
  }

  layer_matrices <- list()

  for (layer_name in matching_layers) {
    current_matrix <- tryCatch(
      SeuratObject::LayerData(
        object,
        assay = assay,
        layer = layer_name,
        fast = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(current_matrix) || ncol(current_matrix) == 0) {
      next
    }

    if (
      length(intersect(
        rownames(current_matrix),
        colnames(object)
      )) >
      length(intersect(
        colnames(current_matrix),
        colnames(object)
      ))
    ) {
      current_matrix <- Matrix::t(current_matrix)
    }

    retained_cells <- intersect(
      colnames(current_matrix),
      colnames(object)
    )

    if (length(retained_cells) > 0) {
      layer_matrices[[layer_name]] <- current_matrix[
        ,
        retained_cells,
        drop = FALSE
      ]
    }
  }

  if (length(layer_matrices) == 0) {
    stop(
      "No valid ",
      layer,
      " layers could be extracted."
    )
  }

  common_genes <- Reduce(
    intersect,
    lapply(layer_matrices, rownames)
  )

  if (length(common_genes) < 2) {
    stop(
      "Too few shared genes across ",
      layer,
      " layers."
    )
  }

  layer_matrices <- lapply(
    layer_matrices,
    function(x) {
      x[
        common_genes,
        ,
        drop = FALSE
      ]
    }
  )

  joined_matrix <- do.call(
    cbind,
    layer_matrices
  )

  joined_matrix <- joined_matrix[
    ,
    !duplicated(colnames(joined_matrix)),
    drop = FALSE
  ]

  final_cells <- intersect(
    colnames(object),
    colnames(joined_matrix)
  )

  joined_matrix <- joined_matrix[
    ,
    final_cells,
    drop = FALSE
  ]

  if (ncol(joined_matrix) < 2) {
    stop(
      "Fewer than two cells remained after joining ",
      layer,
      " layers."
    )
  }

  as(joined_matrix, "dgCMatrix")
}

get_reduction_name <- function(object) {
  available_reductions <- Reductions(object)

  if ("tsne" %in% available_reductions) {
    return("tsne")
  }

  if ("umap" %in% available_reductions) {
    return("umap")
  }

  if ("pca" %in% available_reductions) {
    return("pca")
  }

  stop(
    "No t-SNE, UMAP, or PCA reduction was found in the saved High-CNV object."
  )
}

program_z_score_table <- function(
    data_mat,
    gene_sets,
    min_genes = 5
) {
  result <- lapply(gene_sets, function(gene_set) {
    retained_genes <- intersect(
      unique(gene_set),
      rownames(data_mat)
    )

    if (length(retained_genes) < min_genes) {
      return(rep(NA_real_, ncol(data_mat)))
    }

    x <- as.matrix(
      data_mat[
        retained_genes,
        ,
        drop = FALSE
      ]
    )

    gene_means <- rowMeans(
      x,
      na.rm = TRUE
    )

    gene_sds <- apply(
      x,
      1,
      stats::sd,
      na.rm = TRUE
    )

    valid_genes <- is.finite(gene_means) &
      is.finite(gene_sds) &
      gene_sds > 0

    if (sum(valid_genes) < min_genes) {
      return(rep(NA_real_, ncol(data_mat)))
    }

    z <- sweep(
      x[valid_genes, , drop = FALSE],
      1,
      gene_means[valid_genes],
      "-"
    )

    z <- sweep(
      z,
      1,
      gene_sds[valid_genes],
      "/"
    )

    colMeans(
      z,
      na.rm = TRUE
    )
  })

  result <- as.data.frame(
    result,
    check.names = FALSE
  )

  rownames(result) <- colnames(data_mat)

  result
}

cluster_number <- function(x) {
  x <- as.character(x)

  output <- suppressWarnings(
    as.integer(
      sub(
        "^.*?([0-9]+)\\s*$",
        "\\1",
        x
      )
    )
  )

  output[
    !grepl(
      "[0-9]+\\s*$",
      x
    )
  ] <- NA_integer_

  output
}

normalise_cluster_label <- function(x) {
  x <- trimws(as.character(x))
  numbers <- cluster_number(x)
  output <- x
  output[!is.na(numbers)] <- paste0("Cluster ", numbers[!is.na(numbers)])
  output
}

mark_step_done("01.4")

# STEP 01.5 — LOAD AND VALIDATE EXISTING HIGH-CNV OBJECTS

require_previous_step("01.4", "01.5")

highcnv_obj <- readRDS(HIGHCNV_OBJECT_RDS)
highcnv_de_obj <- readRDS(HIGHCNV_DE_OBJECT_RDS)

if (!inherits(highcnv_obj, "Seurat")) {
  stop(
    "highcnv_obj is not a Seurat object."
  )
}

if (!inherits(highcnv_de_obj, "Seurat")) {
  stop(
    "highcnv_de_obj is not a Seurat object."
  )
}

required_metadata <- c(
  "sample",
  "Cluster_label",
  "CNV_class"
)

if (!all(required_metadata %in% colnames(highcnv_obj@meta.data))) {
  stop(
    "The High-CNV object needs this metadata:\n",
    paste(required_metadata, collapse = ", ")
  )
}

if (!(TARGET_GENE %in% rownames(highcnv_obj))) {
  stop(
    TARGET_GENE,
    " is absent from the High-CNV object."
  )
}

shared_highcnv_de_cells <- intersect(
  colnames(highcnv_obj),
  colnames(highcnv_de_obj)
)

if (length(shared_highcnv_de_cells) < 2) {
  stop(
    "High-CNV cell IDs do not sufficiently overlap the prerequisite DE object."
  )
}

highcnv_analysis_obj <- subset(
  highcnv_de_obj,
  cells = shared_highcnv_de_cells
)

highcnv_analysis_obj$Cluster_label <- highcnv_obj$Cluster_label[
  match(
    colnames(highcnv_analysis_obj),
    colnames(highcnv_obj)
  )
]

highcnv_analysis_obj$sample <- highcnv_obj$sample[
  match(
    colnames(highcnv_analysis_obj),
    colnames(highcnv_obj)
  )
]

highcnv_analysis_obj$CNV_class <- highcnv_obj$CNV_class[
  match(
    colnames(highcnv_analysis_obj),
    colnames(highcnv_obj)
  )
]

cnv_labels <- as.character(highcnv_analysis_obj$CNV_class)

if (
  any(
    is.na(cnv_labels) |
    cnv_labels != "High-CNV malignant"
  )
) {
  stop(
    "Only non-missing High-CNV malignant cells are allowed in this analysis."
  )
}

normalised_subclusters <- normalise_cluster_label(
  highcnv_analysis_obj$Cluster_label
)

observed_subclusters <- sort(unique(normalised_subclusters))

if (!setequal(observed_subclusters, FIXED_CLUSTER_ORDER)) {
  stop(
    "The retained High-CNV cells do not match the predefined ten-subcluster set. Observed: ",
    paste(observed_subclusters, collapse = ", ")
  )
}

highcnv_analysis_obj$Malignant_subcluster <- factor(
  normalised_subclusters,
  levels = FIXED_CLUSTER_ORDER
)

if (any(is.na(highcnv_analysis_obj$Malignant_subcluster))) {
  stop("Missing or unrecognised Cluster_label values were detected in the retained High-CNV cells.")
}

MALIGNANT_ASSAY <- if (
  "RNA_DE_HIGH_CNV" %in%
  names(highcnv_analysis_obj@assays)
) {
  "RNA_DE_HIGH_CNV"
} else if (
  "RNA" %in%
  names(highcnv_analysis_obj@assays)
) {
  "RNA"
} else {
  stop(
    "No RNA_DE_HIGH_CNV or RNA assay was found."
  )
}

malignant_data <- get_existing_layer(
  highcnv_analysis_obj,
  MALIGNANT_ASSAY,
  "data"
)

valid_data_cells <- intersect(
  colnames(highcnv_analysis_obj),
  colnames(malignant_data)
)

if (length(valid_data_cells) < 100) {
  stop(
    "Fewer than 100 cells remained after checking the expression layer."
  )
}

highcnv_analysis_obj <- subset(
  highcnv_analysis_obj,
  cells = valid_data_cells
)

remaining_subclusters <- sort(unique(
  as.character(highcnv_analysis_obj$Malignant_subcluster)
))

if (!setequal(remaining_subclusters, FIXED_CLUSTER_ORDER)) {
  stop(
    "After expression-layer filtering, all ten predefined High-CNV subclusters must remain. Observed: ",
    paste(remaining_subclusters, collapse = ", ")
  )
}

highcnv_analysis_obj$Malignant_subcluster <- factor(
  as.character(highcnv_analysis_obj$Malignant_subcluster),
  levels = FIXED_CLUSTER_ORDER
)

malignant_data <- malignant_data[
  ,
  colnames(highcnv_analysis_obj),
  drop = FALSE
]

if (!(TARGET_GENE %in% rownames(malignant_data))) {
  stop(
    TARGET_GENE,
    " is absent from the expression layer used for programme scoring and DotPlot."
  )
}

highcnv_embedding_obj <- subset(
  highcnv_obj,
  cells = colnames(highcnv_analysis_obj)
)

highcnv_embedding_obj$Malignant_subcluster <- factor(
  as.character(
    highcnv_analysis_obj$Malignant_subcluster[
      match(
        colnames(highcnv_embedding_obj),
        colnames(highcnv_analysis_obj)
      )
    ]
  ),
  levels = FIXED_CLUSTER_ORDER
)

if (any(is.na(highcnv_embedding_obj$Malignant_subcluster))) {
  stop("Embedding metadata could not be aligned to the fixed High-CNV subcluster order.")
}

writeLines(
  c(
    "Existing-cluster continuity passed.",
    paste0(
      "Retained High-CNV cells: ",
      ncol(highcnv_analysis_obj)
    ),
    paste0(
      "Existing High-CNV subclusters: ",
      paste(
        levels(
          highcnv_analysis_obj$Malignant_subcluster
        ),
        collapse = ", "
      )
    ),
    paste0(
      "Expression assay used: ",
      MALIGNANT_ASSAY
    ),
    "No FindClusters or FindAllMarkers operation was run in this script."
  ),
  file.path(
    OUTPUT_DIR,
    "Existing_cluster_continuity.txt"
  )
)

mark_step_done("01.5")

# STEP 01.6 — FUNCTIONAL PROGRAMMES AND FIGURE 7

require_previous_step("01.5", "01.6")

hallmark <- tryCatch(
  msigdbr::msigdbr(
    species = "Homo sapiens",
    collection = "H"
  ),
  error = function(e) {
    msigdbr::msigdbr(
      species = "Homo sapiens",
      category = "H"
    )
  }
)

hallmark_sets <- split(
  hallmark$gene_symbol,
  hallmark$gs_name
)

program_gene_sets <- list(
  E2F_targets = hallmark_sets[["HALLMARK_E2F_TARGETS"]],
  G2M_checkpoint = hallmark_sets[["HALLMARK_G2M_CHECKPOINT"]],
  DNA_repair = hallmark_sets[["HALLMARK_DNA_REPAIR"]],
  MYC_targets = hallmark_sets[["HALLMARK_MYC_TARGETS_V1"]],
  EMT = hallmark_sets[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
  Hypoxia = hallmark_sets[["HALLMARK_HYPOXIA"]],
  Stress_response = c(
    "ATF3",
    "DDIT3",
    "HSPA1A",
    "HSPA1B",
    "JUN",
    "FOS",
    "HSP90AA1",
    "DNAJB1",
    "XBP1",
    "HMOX1"
  ),
  Senescence_like = c(
    "CDKN1A",
    "CDKN2A",
    "GDF15",
    "SERPINE1",
    "IL6",
    "CXCL8",
    "IGFBP7",
    "GLB1",
    "MMP3",
    "MMP9"
  ),
  Epithelial_differentiation = c(
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT15",
    "KRT17",
    "CLDN1",
    "CLDN4"
  )
)

program_gene_set_sizes <- vapply(
  program_gene_sets,
  function(gene_set) {
    length(intersect(unique(gene_set), rownames(malignant_data)))
  },
  integer(1)
)

if (any(program_gene_set_sizes < 5)) {
  stop(
    "The following programme(s) have fewer than five represented genes in the expression matrix: ",
    paste(
      names(program_gene_set_sizes)[program_gene_set_sizes < 5],
      collapse = ", "
    )
  )
}

program_gene_set_export <- do.call(
  rbind,
  lapply(names(program_gene_sets), function(program_name) {
    data.frame(
      Programme = program_name,
      Gene = intersect(
        unique(program_gene_sets[[program_name]]),
        rownames(malignant_data)
      ),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  program_gene_set_export,
  file.path(OUTPUT_DIR, "Functional_programme_gene_sets_used.csv"),
  row.names = FALSE
)

program_scores <- program_z_score_table(
  malignant_data,
  program_gene_sets
)

colnames(program_scores) <- paste0(
  "ProgZ_",
  colnames(program_scores)
)

mybl2_expression <- as.numeric(
  malignant_data[
    TARGET_GENE,
    colnames(highcnv_analysis_obj)
  ]
)

names(mybl2_expression) <- colnames(highcnv_analysis_obj)

highcnv_meta <- meta_with_cell_id(
  highcnv_analysis_obj@meta.data
) %>%
  dplyr::mutate(
    MYBL2_logNorm = mybl2_expression[Cell],
    sample = as.character(sample),
    Malignant_subcluster = factor(Malignant_subcluster)
  ) %>%
  dplyr::left_join(
    program_scores %>%
      tibble::rownames_to_column("Cell"),
    by = "Cell"
  )

program_cols <- grep(
  "^ProgZ_",
  colnames(highcnv_meta),
  value = TRUE
)

program_summary <- highcnv_meta %>%
  dplyr::group_by(Malignant_subcluster) %>%
  dplyr::summarise(
    Cells = dplyr::n(),
    Samples_represented = dplyr::n_distinct(sample),
    Mean_MYBL2 = mean(MYBL2_logNorm, na.rm = TRUE),
    dplyr::across(
      dplyr::all_of(program_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Malignant_subcluster = factor(
      as.character(Malignant_subcluster),
      levels = FIXED_CLUSTER_ORDER
    )
  ) %>%
  dplyr::arrange(Malignant_subcluster)

write.csv(
  highcnv_meta,
  file.path(OUTPUT_DIR, "HighCNV_cell_metadata_MYBL2_programZ.csv"),
  row.names = FALSE
)

write.csv(
  program_summary,
  file.path(OUTPUT_DIR, "HighCNV_subcluster_programZ_summary.csv"),
  row.names = FALSE
)

existing_reduction <- get_reduction_name(highcnv_embedding_obj)

reduction_x <- ifelse(
  existing_reduction == "tsne",
  "t-SNE 1",
  ifelse(existing_reduction == "umap", "UMAP 1", "PC 1")
)

reduction_y <- ifelse(
  existing_reduction == "tsne",
  "t-SNE 2",
  ifelse(existing_reduction == "umap", "UMAP 2", "PC 2")
)

p_subclusters <- DimPlot(
  highcnv_embedding_obj,
  reduction = existing_reduction,
  group.by = "Malignant_subcluster",
  label = TRUE,
  repel = TRUE
) +
  labs(
    title = NULL,
    x = reduction_x,
    y = reduction_y,
    color = "High-CNV subtype"
  ) +
  theme_manuscript(show_grid = FALSE) +
  theme(
    plot.title = element_blank()
  )

mybl2_dot_raw <- Seurat::DotPlot(
  highcnv_analysis_obj,
  features = TARGET_GENE,
  assay = MALIGNANT_ASSAY,
  group.by = "Malignant_subcluster"
)

mybl2_dot_data <- mybl2_dot_raw$data %>%
  dplyr::mutate(
    cluster_label = normalise_cluster_label(id),
    avg.exp = as.numeric(avg.exp),
    avg.exp.scaled = as.numeric(avg.exp.scaled),
    pct.exp = as.numeric(pct.exp)
  ) %>%
  dplyr::left_join(
    program_summary %>%
      dplyr::transmute(
        cluster_label = as.character(Malignant_subcluster),
        Mean_MYBL2
      ),
    by = "cluster_label"
  ) %>%
  dplyr::mutate(
    cluster_label = factor(
      cluster_label,
      levels = FIXED_CLUSTER_ORDER
    )
  ) %>%
  dplyr::arrange(cluster_label)

if (nrow(mybl2_dot_data) != length(FIXED_CLUSTER_ORDER)) {
  stop(
    "The MYBL2 DotPlot did not contain exactly one record for each High-CNV subtype."
  )
}

write.csv(
  mybl2_dot_data,
  file.path(OUTPUT_DIR, "Figure_7B_MYBL2_DotPlot_data.csv"),
  row.names = FALSE
)

p_mybl2_only <- ggplot(
  mybl2_dot_data,
  aes(
    x = 1,
    y = cluster_label
  )
) +
  geom_point(
    aes(
      size = pct.exp,
      color = avg.exp.scaled
    ),
    alpha = 0.95
  ) +
  scale_x_continuous(
    breaks = 1,
    labels = TARGET_GENE,
    limits = VERTICAL_DOTPLOT_X_LIMITS,
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(
    limits = rev(FIXED_CLUSTER_ORDER),
    drop = FALSE
  ) +
  scale_size_continuous(
    range = VERTICAL_DOTPLOT_POINT_RANGE,
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    name = "Percent expressed"
  ) +
  scale_color_gradientn(
    colours = c("#EEEAF7", COL_NORMAL, COL_TUMOR),
    name = "Average expression",
    na.value = "grey85"
  ) +
  guides(
    size = guide_legend(
      order = 1,
      direction = "vertical",
      override.aes = list(
        size = c(2.5, 4.8, 7.2, 9.6, 12.0)
      )
    ),
    color = guide_colorbar(
      order = 2,
      direction = "vertical",
      barwidth = grid::unit(0.30, "cm"),
      barheight = grid::unit(1.65, "cm")
    )
  ) +
  labs(
    x = NULL,
    y = "High-CNV subtype"
  ) +
  theme_manuscript(show_grid = FALSE) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(
      family = FONT_FAMILY,
      face = "bold",
      size = AXIS_TITLE_PT,
      color = "black",
      margin = ggplot2::margin(t = 2, unit = "pt")
    ),
    axis.text.y = element_text(
      family = FONT_FAMILY,
      size = AXIS_TEXT_PT,
      color = "black"
    ),
    axis.title.y = element_text(
      margin = ggplot2::margin(r = 3, unit = "pt")
    ),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.justification = "center",
    legend.key.width = grid::unit(0.32, "cm"),
    legend.key.height = grid::unit(0.32, "cm"),
    legend.box.spacing = grid::unit(2, "pt"),
    legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt")
  )

summary_long_all <- program_summary %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(program_cols),
    names_to = "Program",
    values_to = "Mean_program_Z"
  ) %>%
  dplyr::mutate(
    Program = sub("^ProgZ_", "", Program),
    Program = dplyr::recode(
      Program,
      E2F_targets = "E2F targets",
      G2M_checkpoint = "G2/M checkpoint",
      DNA_repair = "DNA repair",
      MYC_targets = "MYC targets",
      EMT = "EMT",
      Hypoxia = "Hypoxia",
      Stress_response = "Stress response",
      Senescence_like = "Senescence-like",
      Epithelial_differentiation = "Epithelial differentiation"
    ),
    Program = factor(
      Program,
      levels = PROGRAMME_DISPLAY_ORDER
    ),
    Malignant_subcluster = factor(
      as.character(Malignant_subcluster),
      levels = FIXED_CLUSTER_ORDER
    )
  )

programme_scale_limit <- max(
  abs(summary_long_all$Mean_program_Z),
  na.rm = TRUE
)

if (!is.finite(programme_scale_limit) || programme_scale_limit <= 0) {
  programme_scale_limit <- 1
}

make_programme_heatmap <- function(heatmap_data, x_title) {
  ggplot(
    heatmap_data,
    aes(
      x = Malignant_subcluster,
      y = Program,
      fill = Mean_program_Z
    )
  ) +
    geom_tile(
      color = "white",
      linewidth = GRID_LWD
    ) +
    geom_text(
      aes(
        label = sprintf("%.2f", Mean_program_Z)
      ),
      size = HEATMAP_FONT_PT / 3
    ) +
    scale_fill_gradient2(
      low = COL_NORMAL,
      mid = "white",
      high = COL_TUMOR,
      midpoint = 0,
      limits = c(-programme_scale_limit, programme_scale_limit),
      name = "Relative programme\nactivity (Z)",
      na.value = "grey92"
    ) +
    labs(
      x = x_title,
      y = NULL
    ) +
    theme_manuscript(show_grid = FALSE) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(
        family = FONT_FAMILY,
        size = HEATMAP_ROW_FONT_PT,
        color = "black"
      ),
      axis.text.x = element_text(
        family = FONT_FAMILY,
        angle = 35,
        hjust = 1,
        vjust = 1,
        size = HEATMAP_COL_FONT_PT,
        color = "black"
      ),
      legend.title = element_text(
        family = FONT_FAMILY,
        size = HEATMAP_LEGEND_FONT_PT,
        face = "bold",
        color = "black"
      ),
      legend.text = element_text(
        family = FONT_FAMILY,
        size = HEATMAP_LEGEND_FONT_PT,
        color = "black"
      )
    )
}

summary_long_panel_c <- summary_long_all %>%
  dplyr::filter(
    as.character(Malignant_subcluster) %in% PANEL_C_CLUSTER_ORDER
  ) %>%
  dplyr::mutate(
    Malignant_subcluster = factor(
      as.character(Malignant_subcluster),
      levels = PANEL_C_CLUSTER_ORDER
    )
  )

if (!setequal(
  unique(as.character(summary_long_panel_c$Malignant_subcluster)),
  PANEL_C_CLUSTER_ORDER
)) {
  stop(
    "Panel c could not retain exactly the requested subclusters: ",
    paste(PANEL_C_CLUSTER_ORDER, collapse = ", ")
  )
}

write.csv(
  summary_long_panel_c,
  file.path(
    OUTPUT_DIR,
    "Figure_7C_functional_programme_heatmap_Cluster0_Cluster8_data.csv"
  ),
  row.names = FALSE
)

p_programs <- make_programme_heatmap(
  summary_long_panel_c,
  x_title = "High-CNV malignant subtype"
)

writeLines(
  c(
    "Functional programme interpretation note",
    "Programme heatmaps show descriptive, relative gene-wise Z-score averages at the subcluster level.",
    "They are not inferential statistical tests; do not describe programme differences as statistically significant without a separate sample-level or pseudobulk analysis.",
    "The Stress_response, Senescence_like, and Epithelial_differentiation programmes are curated custom gene sets.",
    "The exact genes used for every programme are exported in Functional_programme_gene_sets_used.csv."
  ),
  PROGRAMME_INTERPRETATION_NOTE
)

p_subclusters_tagged <- p_subclusters +
  ggplot2::labs(tag = "a") +
  ggplot2::theme(
    plot.tag.position = PANEL_A_TAG_POSITION,
    plot.margin = PANEL_A_TAG_MARGIN
  )

p_mybl2_only_tagged <- p_mybl2_only +
  ggplot2::labs(tag = "b") +
  ggplot2::theme(
    plot.tag.position = PANEL_B_TAG_POSITION,
    plot.margin = PANEL_B_TAG_MARGIN
  )

p_programs_tagged <- p_programs +
  ggplot2::labs(tag = "c") +
  ggplot2::theme(
    plot.tag.position = PANEL_C_TAG_POSITION,
    plot.margin = PANEL_C_TAG_MARGIN
  )

top_row_7 <- (
  p_subclusters_tagged |
    p_mybl2_only_tagged
) +
  patchwork::plot_layout(
    widths = c(TSNE_PANEL_WIDTH, MYBL2_PANEL_WIDTH)
  )

bottom_row_7 <- (
  patchwork::plot_spacer() |
    p_programs_tagged |
    patchwork::plot_spacer()
) +
  patchwork::plot_layout(
    widths = c(
      PANEL_C_LEFT_SPACER_WIDTH,
      PANEL_C_HEATMAP_WIDTH,
      PANEL_C_RIGHT_SPACER_WIDTH
    )
  )

figure_7 <- (
  top_row_7 /
    bottom_row_7
) +
  patchwork::plot_layout(
    heights = c(FIGURE7_TOP_ROW_HEIGHT, FIGURE7_BOTTOM_ROW_HEIGHT)
  )

save_plot_all_formats(
  plot_obj = figure_7,
  filename_stem = "Figure_7_HighCNV_subtypes_MYBL2_functional_programmes_FixedAllPanelLabels",
  dir_path = OUTPUT_DIR,
  width = FIGURE7_COMBINED_W,
  height = FIGURE7_COMBINED_H
)

mark_step_done("01.6")

# STEP 01.7 — EXPORT OBJECTS AND FINAL STATUS

require_previous_step("01.6", "01.7")

saveRDS(
  list(
    highcnv_analysis_object = highcnv_analysis_obj,
    highcnv_metadata = highcnv_meta,
    programZ_summary = program_summary,
    mybl2_dotplot_data = mybl2_dot_data,
    programme_gene_sets = program_gene_set_export
  ),
  file.path(
    OUTPUT_DIR,
    "Functional_Programme_Objects.rds"
  )
)

summary_lines <- c(
  "Functional programme analysis completed.",
  paste0(
    "Existing clustering resolution retained: ",
    RESOLUTION_MALIGNANT_FROM_PREVIOUS_PIPELINE
  ),
  paste0(
    "High-CNV cells analysed: ",
    ncol(highcnv_analysis_obj)
  ),
  paste0(
    "Existing High-CNV subclusters: ",
    paste(
      levels(highcnv_analysis_obj$Malignant_subcluster),
      collapse = ", "
    )
  ),
  "Figure 7a: existing High-CNV subtype embedding with the technical in-image title removed.",
  "Figure 7b: MYBL2 DotPlot with a slightly wider panel for improved readability; colour denotes average expression.",
  "Figure 7c: descriptive functional-programme heatmap for Cluster 0 and Cluster 8 only.",
  "Functional programmes were gene-wise Z-scored before averaging and are descriptive, not inferential statistical comparisons.",
  "Custom programme gene sets: Stress_response, Senescence_like, and Epithelial_differentiation; exact genes are exported separately.",
  "Publication Figure 7: all High-CNV subtypes in panels a-b and functional programmes for Cluster 0 and Cluster 8 in panel c.",
  "Figure 7 section frames: each complete a, b, and c section is enclosed by the same thin black outer border for clear visual division.",
  "Figure 7 rendering size: enlarged export canvas, larger panel-row allocations, reduced shared margins, and reduced DotPlot vertical white space.",
  "Panel b: vertical ten-subtype MYBL2 DotPlot with enlarged symbols, wider panel allocation, and bottom legends to use the full panel height.",
  "Final layout revision: panel b uses compact right-side legends; panel c is narrowed and centered with plot_spacer() on both sides.",
  "Panel c refinement: shortened bottom-row height and panel-c-only top/left tag clearance prevent overlap between tag c and the Stress response label.",
  "Panel b refinement: panel-b-only top/left tag clearance prevents overlap between tag b and the Cluster 0 label.",
  "Panel a refinement: panel-a-only top/left tag clearance positions tag a away from the t-SNE embedding and cluster labels."
)

writeLines(
  summary_lines,
  file.path(OUTPUT_DIR, "Functional_Programme_Analysis_summary.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "Functional_Programme_sessionInfo.txt")
)

cat("\n============================================================\n")
cat("FUNCTIONAL PROGRAMME ANALYSIS FINISHED\n")
cat("============================================================\n")
cat("All outputs: ", OUTPUT_DIR, "\n", sep = "")
cat("Figure 7: Figure_7_HighCNV_subtypes_MYBL2_functional_programmes_FixedAllPanelLabels.{png,tiff,pdf}\n")
cat("============================================================\n")

mark_step_done("01.7")

message(
  "All numbered steps finished. Outputs are in: ",
  OUTPUT_DIR,
  "\nFigure 7 includes all ten High-CNV malignant subtypes in panels a-b; panel c includes Cluster 0 and Cluster 8 only."
)

# MODULE 04 — CellChat communication analysis and Figure 13 outputs

rm(list = ls(all.names = TRUE))
gc()

set.seed(123)
options(
  stringsAsFactors = FALSE,
  scipen = 100,
  timeout = 7200,
  width = 160
)

# 1. PATHS AND ANALYSIS SETTINGS

BASE_RESULTS_DIR <- "D:/LSCC/Results_LSCC"
SC_PROJECT_DIR <- "D:/LSCC/ScRNAseq_Results/GSE206332"
OUTPUT_DIR <- file.path(BASE_RESULTS_DIR, "CellChat")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

FULL_OBJECT_RDS_OVERRIDE <- NA_character_
HIGHCNV_OBJECT_RDS_OVERRIDE <- NA_character_

ASSAY_TO_USE <- "RNA"
MIN_CELLS_PER_GROUP <- 15L
MAX_CELLS_PER_SAMPLE_GROUP <- 350L
TOP_LR_TO_SHOW <- 5L
PVAL_CUTOFF <- 0.05

FORCE_RECOMPUTE_CELLCHAT <- FALSE

SELECTED_SUBCLUSTERS <- c("Cluster 0", "Cluster 8")

TARGET_POPULATIONS <- c(
  "T cell", "NK cell", "B cell", "Myeloid", "Fibroblast",
  "Endothelial", "Epithelial"
)

# 2. STANDARD MANUSCRIPT FIGURE SETTINGS

FONT_FAMILY <- "Arial"
FIG_DPI <- 600L
FIG_BACKGROUND <- "white"

FIG_SINGLE_W <- 3.50
FIG_SINGLE_H <- 4.20
FIG_DOUBLE_W <- 7.20
FIG_DOUBLE_H <- 5.40

FIGURE13_COMBINED_W <- 15.80
FIGURE13_COMBINED_H <- 8.80

BASE_TEXT_PT <- 9.5
AXIS_TEXT_PT <- 8.5
AXIS_TITLE_PT <- 10
LEGEND_TEXT_PT <- 8.5
LEGEND_TITLE_PT <- 9
PLOT_TITLE_PT <- 10
PANEL_TAG_PT <- 14

PANEL_BORDER_LWD <- 0.45
AXIS_LWD <- 0.40
GEOM_LWD <- 0.55
GRID_LWD <- 0.30

HEATMAP_FONT_PT <- 8
HEATMAP_ROW_FONT_PT <- 7
HEATMAP_COL_FONT_PT <- 7
HEATMAP_LEGEND_FONT_PT <- 8

COL_NORMAL <- "#7470B2"
COL_TUMOR <- "#D81B60"
COL_MYBL2_LOW <- COL_NORMAL
COL_MYBL2_HIGH <- COL_TUMOR
COL_UP <- COL_TUMOR
COL_DOWN <- COL_NORMAL
COL_NS <- "#B8B8B8"

CELLCHAT_PALETTE <- grDevices::colorRampPalette(
  c("#365F86", "#8DBBD5", "#D7D6C7", "#E7B58F", "#B96B63")
)(100)

CELLCHAT_NODE_COLOURS <- c(
  "Source cluster" = "#365F86",
  "T cell" = "#A9CBE0",
  "NK cell" = "#E6A15A",
  "B cell" = "#E6EDF2",
  "Myeloid" = "#E5A96E",
  "Fibroblast" = "#E5B59F",
  "Endothelial" = "#D7D3A8",
  "Epithelial" = "#B7D2E3"
)

NETWORK_EDGE_COLOUR <- "#A98274"

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
      legend.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = LEGEND_TITLE_PT,
        color = "black"
      ),
      legend.text = ggplot2::element_text(
        family = FONT_FAMILY,
        size = LEGEND_TEXT_PT,
        color = "black"
      ),
      legend.position = legend_position,
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
          color = "grey94",
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
                                  dir_path,
                                  width = FIG_DOUBLE_W,
                                  height = FIG_DOUBLE_H) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  filename_stem <- sub("\\.(png|tiff|tif|pdf)$", "", filename_stem, ignore.case = TRUE)

  png_file <- file.path(dir_path, paste0(filename_stem, ".png"))
  tiff_file <- file.path(dir_path, paste0(filename_stem, ".tiff"))
  pdf_file <- file.path(dir_path, paste0(filename_stem, ".pdf"))

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

# 3. REQUIRED LIBRARIES

required_pkgs <- c(
  "Seurat", "SeuratObject", "Matrix", "dplyr", "ggplot2",
  "patchwork", "grid", "future", "CellChat"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_pkgs) > 0L) {
  stop(
    "These packages are missing from the active R library:\n",
    paste(missing_pkgs, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(future)
  library(CellChat)
})

future::plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)

# 4. GENERAL HELPER FUNCTIONS

placeholder_plot <- function(text_string) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = text_string,
      size = 4.2,
      family = FONT_FAMILY
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_void(base_family = FONT_FAMILY)
}

cluster_number <- function(x) {
  x <- as.character(x)
  value <- suppressWarnings(as.integer(sub("^.*?([0-9]+)\\s*$", "\\1", x)))
  value[!grepl("[0-9]+\\s*$", x)] <- NA_integer_
  value
}

normalise_cluster_label <- function(x) {
  x <- trimws(as.character(x))
  n <- cluster_number(x)
  out <- x
  out[!is.na(n)] <- paste0("Cluster ", n[!is.na(n)])
  out
}

sort_cluster_labels <- function(x) {
  x[order(cluster_number(x), x, na.last = TRUE)]
}

first_existing_file <- function(paths, label) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  found <- paths[file.exists(paths)]

  if (length(found) == 0L) {
    stop(label, " was not found. Checked:\n", paste(paths, collapse = "\n"))
  }

  normalizePath(found[1], winslash = "/", mustWork = TRUE)
}

safe_read_rds <- function(path) {
  tryCatch(readRDS(path), error = function(e) NULL)
}

get_existing_layer <- function(object, assay, layer = "data") {
  if (!(assay %in% names(object@assays))) {
    stop("Assay not found: ", assay)
  }

  DefaultAssay(object) <- assay
  assay_obj <- object[[assay]]

  layer_names <- tryCatch(
    SeuratObject::Layers(assay_obj),
    error = function(e) character(0)
  )

  matching_layers <- layer_names[
    grepl(paste0("^", layer, "(\\.|$)"), layer_names)
  ]

  if (length(matching_layers) == 0L) {
    old_matrix <- tryCatch(
      Seurat::GetAssayData(object, assay = assay, slot = layer),
      error = function(e) NULL
    )

    if (!is.null(old_matrix) && ncol(old_matrix) > 0L) {
      return(as(old_matrix, "dgCMatrix"))
    }

    stop("No usable ", layer, " layer was found in assay ", assay)
  }

  mats <- list()

  for (layer_name in matching_layers) {
    mat <- tryCatch(
      SeuratObject::LayerData(
        object,
        assay = assay,
        layer = layer_name,
        fast = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(mat) || ncol(mat) == 0L) next

    if (
      length(intersect(rownames(mat), colnames(object))) >
      length(intersect(colnames(mat), colnames(object)))
    ) {
      mat <- Matrix::t(mat)
    }

    usable_cells <- intersect(colnames(mat), colnames(object))

    if (length(usable_cells) > 0L) {
      mats[[layer_name]] <- mat[, usable_cells, drop = FALSE]
    }
  }

  if (length(mats) == 0L) {
    stop("No valid ", layer, " layers could be extracted.")
  }

  common_genes <- Reduce(intersect, lapply(mats, rownames))

  if (length(common_genes) < 2L) {
    stop("Too few shared genes across ", layer, " layers.")
  }

  mats <- lapply(mats, function(m) m[common_genes, , drop = FALSE])
  joined <- do.call(cbind, mats)
  joined <- joined[, !duplicated(colnames(joined)), drop = FALSE]

  final_cells <- intersect(colnames(object), colnames(joined))
  joined <- joined[, final_cells, drop = FALSE]

  if (ncol(joined) < 2L) {
    stop("Fewer than two cells remained after joining expression layers.")
  }

  as(joined, "dgCMatrix")
}

get_sample_column <- function(meta) {
  candidates <- c(
    "sample", "Sample", "orig.ident", "orig_ident", "patient", "Patient"
  )
  found <- candidates[candidates %in% colnames(meta)]

  if (length(found) == 0L) {
    stop(
      "No sample metadata column was found. Checked: ",
      paste(candidates, collapse = ", ")
    )
  }

  found[1]
}

balance_cells_by_sample_group <- function(meta, sample_col, group_col, max_n) {
  split_ids <- split(
    rownames(meta),
    paste(meta[[sample_col]], meta[[group_col]], sep = "__")
  )

  kept <- unlist(
    lapply(split_ids, function(ids) {
      if (length(ids) <= max_n) return(ids)
      sample(ids, max_n)
    }),
    use.names = FALSE
  )

  unique(kept)
}

find_source_group <- function(group_names, short_cluster) {
  pattern <- paste0(
    "(^High-CNV malignant:[[:space:]]*)?",
    gsub(" ", "[[:space:]]*", short_cluster),
    "$"
  )

  hit <- group_names[grepl(pattern, group_names, ignore.case = TRUE)]

  if (length(hit) == 0L) return(NA_character_)
  hit[1]
}

is_compatible_cellchat <- function(x) {
  if (!inherits(x, "CellChat")) return(FALSE)
  if (is.null(x@net$count) || is.null(x@net$weight)) return(FALSE)

  group_names <- tryCatch(levels(x@idents), error = function(e) character(0))

  selected_found <- vapply(
    SELECTED_SUBCLUSTERS,
    function(cl) !is.na(find_source_group(group_names, cl)),
    logical(1)
  )

  all(selected_found)
}

get_interaction_label <- function(d) {
  if ("interaction_name_2" %in% colnames(d)) {
    return(as.character(d$interaction_name_2))
  }

  if ("interaction_name" %in% colnames(d)) {
    return(as.character(d$interaction_name))
  }

  if (all(c("ligand", "receptor") %in% colnames(d))) {
    return(paste(d$ligand, d$receptor, sep = " - "))
  }

  rep("Ligand-receptor pair", nrow(d))
}

safe_limits <- function(x, padding_fraction = 0.06) {
  x <- x[is.finite(x)]

  if (length(x) == 0L) return(c(0, 1))

  limits <- range(x)

  if (limits[1] == limits[2]) {
    pad <- max(abs(limits[1]) * padding_fraction, 0.01)
    return(c(limits[1] - pad, limits[2] + pad))
  }

  limits
}

make_network_layout <- function(cluster_id, available_targets) {
  all_target_positions <- data.frame(
    cell = TARGET_POPULATIONS,
    x = c(-1.28, -0.98, -0.52, 0.02, 0.56, 1.02, 1.30),
    y = c(0.40, 1.03, 1.36, 1.48, 1.27, 0.88, 0.30),
    node_type = TARGET_POPULATIONS,
    label_vjust = rep(-1.05, length(TARGET_POPULATIONS)),
    stringsAsFactors = FALSE
  )

  source_node <- data.frame(
    cell = paste0("Cluster ", cluster_id),
    x = 0,
    y = -1.18,
    node_type = "Source cluster",
    label_vjust = 1.65,
    stringsAsFactors = FALSE
  )

  target_nodes <- all_target_positions[
    all_target_positions$cell %in% available_targets,
    ,
    drop = FALSE
  ]

  rbind(source_node, target_nodes)
}

# 5. DETECT A COMPLETED CELLCCHAT OBJECT OR BUILD A NEW ONE

cellchat_rds_candidates <- unique(c(
  file.path(BASE_RESULTS_DIR, "CellChat", "CellChat_SelectedHighCNV_Clusters0_8_TME.rds"),
  file.path(BASE_RESULTS_DIR, "CellChat", "CellChat_object.rds"),
  file.path(BASE_RESULTS_DIR, "CellChat", "LSCC_Stage02_selected_HighCNV_CellChat.rds"),
  file.path(BASE_RESULTS_DIR, "LSCC_02_HighCNV_CellChat", "LSCC_Stage02_selected_HighCNV_CellChat.rds"),
  file.path(BASE_RESULTS_DIR, "CellChat", "GSE206332_all_annotated_cells_existing_resolution0.5_HighCNV_CellChat.rds"),
  file.path(BASE_RESULTS_DIR, "CellChat", "GSE206332_all_annotated_cells_existing_resolution0.5_HighCNV_CellChat_core.rds"),
  list.files(BASE_RESULTS_DIR, pattern = "CellChat.*\\.rds$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
))

cellchat_rds_candidates <- cellchat_rds_candidates[file.exists(cellchat_rds_candidates)]

cellchat <- NULL
cellchat_input_description <- NA_character_

if (!FORCE_RECOMPUTE_CELLCHAT && length(cellchat_rds_candidates) > 0L) {
  for (candidate in cellchat_rds_candidates) {
    candidate_obj <- safe_read_rds(candidate)

    if (is_compatible_cellchat(candidate_obj)) {
      cellchat <- candidate_obj
      cellchat_input_description <- paste0(
        "Reused compatible CellChat RDS: ",
        normalizePath(candidate, winslash = "/", mustWork = TRUE)
      )
      break
    }
  }
}

if (is.null(cellchat)) {
  candidate_rds_dirs <- c(
    file.path(SC_PROJECT_DIR, "Results_Modular", "Step_20_Final_Exports_and_Objects", "rds"),
    file.path(SC_PROJECT_DIR, "Results", "rds")
  )

  full_rds_candidates <- c(
    FULL_OBJECT_RDS_OVERRIDE,
    file.path(candidate_rds_dirs, "FINAL_seurat_object_with_malignant_annotation.rds")
  )

  highcnv_rds_candidates <- c(
    HIGHCNV_OBJECT_RDS_OVERRIDE,
    file.path(candidate_rds_dirs, "HIGH_CNV_MALIGNANT_object_after_lowCNV_removal.rds")
  )

  FULL_OBJECT_RDS <- first_existing_file(
    full_rds_candidates,
    "FINAL_seurat_object_with_malignant_annotation.rds"
  )

  HIGHCNV_OBJECT_RDS <- first_existing_file(
    highcnv_rds_candidates,
    "HIGH_CNV_MALIGNANT_object_after_lowCNV_removal.rds"
  )

  full_obj <- readRDS(FULL_OBJECT_RDS)
  highcnv_obj <- readRDS(HIGHCNV_OBJECT_RDS)

  if (!inherits(full_obj, "Seurat") || !inherits(highcnv_obj, "Seurat")) {
    stop("The prerequisite RDS files are not valid Seurat objects.")
  }

  if (!"final_celltype" %in% colnames(full_obj@meta.data)) {
    stop("The full Seurat object does not have a final_celltype metadata column.")
  }

  highcnv_cluster_column <- c("Cluster_label", "Malignant_subcluster")
  highcnv_cluster_column <- highcnv_cluster_column[
    highcnv_cluster_column %in% colnames(highcnv_obj@meta.data)
  ]

  if (length(highcnv_cluster_column) == 0L) {
    stop("The High-CNV object has neither Cluster_label nor Malignant_subcluster metadata.")
  }

  highcnv_cluster_column <- highcnv_cluster_column[1]
  high_meta <- highcnv_obj@meta.data

  if ("CNV_class" %in% colnames(high_meta)) {
    high_meta <- high_meta[
      is.na(high_meta$CNV_class) | as.character(high_meta$CNV_class) == "High-CNV malignant",
      , drop = FALSE
    ]
  }

  highcnv_map <- data.frame(
    Cell = rownames(high_meta),
    Cluster = normalise_cluster_label(high_meta[[highcnv_cluster_column]]),
    stringsAsFactors = FALSE
  )

  full_meta <- full_obj@meta.data
  full_meta$.__cell__ <- rownames(full_meta)
  sample_col <- get_sample_column(full_meta)

  meta <- full_meta %>%
    dplyr::left_join(highcnv_map, by = c(".__cell__" = "Cell")) %>%
    dplyr::mutate(
      final_celltype = as.character(final_celltype),
      CellChat_group = dplyr::case_when(
        !is.na(Cluster) ~ paste0("High-CNV malignant: ", Cluster),
        final_celltype == "Tcell" ~ "T cell",
        final_celltype == "NKcell" ~ "NK cell",
        final_celltype == "Bcell" ~ "B cell",
        final_celltype == "Myeloid" ~ "Myeloid",
        final_celltype == "Fibroblast" ~ "Fibroblast",
        final_celltype == "Endothelial" ~ "Endothelial",
        final_celltype == "Epithelial" ~ "Epithelial",
        TRUE ~ NA_character_
      ),
      .__sample__ = as.character(.data[[sample_col]])
    )

  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  rownames(meta) <- meta$.__cell__

  required_source_groups <- paste0("High-CNV malignant: ", SELECTED_SUBCLUSTERS)
  missing_sources <- setdiff(required_source_groups, unique(meta$CellChat_group))

  if (length(missing_sources) > 0L) {
    stop(
      "These selected High-CNV subtypes are absent from the input objects:\n",
      paste(missing_sources, collapse = ", ")
    )
  }

  meta <- meta[!is.na(meta$CellChat_group), , drop = FALSE]

  pre_counts <- as.data.frame(table(meta$CellChat_group), stringsAsFactors = FALSE)
  colnames(pre_counts) <- c("CellChat_group", "Cells_before_filtering")
  valid_groups <- pre_counts$CellChat_group[
    pre_counts$Cells_before_filtering >= MIN_CELLS_PER_GROUP
  ]
  meta <- meta[meta$CellChat_group %in% valid_groups, , drop = FALSE]

  missing_after_filter <- setdiff(required_source_groups, unique(meta$CellChat_group))
  if (length(missing_after_filter) > 0L) {
    stop(
      "At least one selected cluster has fewer than MIN_CELLS_PER_GROUP cells:\n",
      paste(missing_after_filter, collapse = ", ")
    )
  }

  kept_cells <- balance_cells_by_sample_group(
    meta = meta,
    sample_col = ".__sample__",
    group_col = "CellChat_group",
    max_n = MAX_CELLS_PER_SAMPLE_GROUP
  )

  meta <- meta[rownames(meta) %in% kept_cells, , drop = FALSE]

  if (!(ASSAY_TO_USE %in% names(full_obj@assays))) {
    ASSAY_TO_USE <- DefaultAssay(full_obj)
  }

  expression_data <- get_existing_layer(full_obj, ASSAY_TO_USE, layer = "data")
  retained_cells <- intersect(rownames(meta), colnames(expression_data))
  meta <- meta[retained_cells, , drop = FALSE]
  expression_data <- expression_data[, retained_cells, drop = FALSE]

  if (!identical(rownames(meta), colnames(expression_data))) {
    meta <- meta[colnames(expression_data), , drop = FALSE]
  }

  if (!identical(rownames(meta), colnames(expression_data))) {
    stop("CellChat expression matrix and metadata are not aligned after filtering.")
  }

  highcnv_groups <- sort_cluster_labels(
    gsub("^High-CNV malignant: ", "", unique(grep("^High-CNV malignant:", meta$CellChat_group, value = TRUE)))
  )
  highcnv_groups <- paste0("High-CNV malignant: ", highcnv_groups)
  tme_groups <- intersect(TARGET_POPULATIONS, unique(meta$CellChat_group))
  final_group_order <- c(highcnv_groups, tme_groups)

  cellchat_meta <- data.frame(
    labels = factor(meta$CellChat_group, levels = final_group_order),
    samples = factor(meta$.__sample__),
    row.names = rownames(meta),
    stringsAsFactors = FALSE
  )

  write.csv(
    data.frame(
      CellChat_group = names(table(cellchat_meta$labels)),
      Cells = as.integer(table(cellchat_meta$labels))
    ),
    file.path(OUTPUT_DIR, "CellChat_group_cell_counts.csv"),
    row.names = FALSE
  )

  cellchat <- CellChat::createCellChat(
    object = expression_data,
    meta = cellchat_meta,
    group.by = "labels",
    do.sparse = TRUE
  )

  data(CellChatDB.human)
  cellchat@DB <- CellChatDB.human

  cellchat <- CellChat::subsetData(cellchat)

  if (nrow(cellchat@data.signaling) < 10L) {
    stop("Fewer than 10 signalling genes overlap the CellChat human database.")
  }

  cellchat <- CellChat::identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, type = "triMean", population.size = FALSE)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = MIN_CELLS_PER_GROUP)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)

  cellchat_input_description <- paste0(
    "Fresh CellChat inference from: ", FULL_OBJECT_RDS,
    " | High-CNV source: ", HIGHCNV_OBJECT_RDS
  )
}

STABLE_CELLCHAT_RDS <- file.path(OUTPUT_DIR, "CellChat_SelectedHighCNV_Clusters0_8_TME.rds")
saveRDS(cellchat, STABLE_CELLCHAT_RDS)

# 6. VALIDATE CELLCCHAT CONTENT AND EXPORT THE COMMUNICATION TABLES

if (is.null(cellchat@net$count) || is.null(cellchat@net$weight)) {
  stop("The CellChat object does not contain aggregated network results.")
}

group_names <- levels(cellchat@idents)
weight_matrix <- cellchat@net$weight

group_size <- as.numeric(table(cellchat@idents))
names(group_size) <- names(table(cellchat@idents))
group_size <- group_size[group_names]

source_groups <- setNames(
  vapply(
    SELECTED_SUBCLUSTERS,
    function(x) find_source_group(group_names, x),
    character(1)
  ),
  SELECTED_SUBCLUSTERS
)

if (any(is.na(source_groups))) {
  stop(
    "Could not find the selected source group(s): ",
    paste(names(source_groups)[is.na(source_groups)], collapse = ", "),
    "\nAvailable groups:\n",
    paste(group_names, collapse = " | ")
  )
}

available_targets <- intersect(TARGET_POPULATIONS, group_names)

if (length(available_targets) < 2L) {
  stop(
    "Too few selected target populations were found in the CellChat object. ",
    "Available groups: ", paste(group_names, collapse = " | ")
  )
}

write.csv(
  data.frame(
    CellChat_group = group_names,
    Cells = as.integer(group_size),
    stringsAsFactors = FALSE
  ),
  file.path(OUTPUT_DIR, "SelectedClusters_CellChat_group_counts.csv"),
  row.names = FALSE
)

all_communications <- tryCatch(
  CellChat::subsetCommunication(cellchat),
  error = function(e) data.frame()
)

if (nrow(all_communications) > 0L) {
  write.csv(
    all_communications,
    file.path(OUTPUT_DIR, "CellChat_all_inferred_communications.csv"),
    row.names = FALSE
  )
}

# 7. PREPARE NETWORK AND TOP-FIVE LIGAND–RECEPTOR DATA

network_table <- dplyr::bind_rows(
  lapply(names(source_groups), function(cluster_label) {
    source_name <- source_groups[[cluster_label]]

    data.frame(
      Cluster = cluster_label,
      Target = available_targets,
      Weight = as.numeric(
        weight_matrix[source_name, available_targets, drop = TRUE]
      ),
      stringsAsFactors = FALSE
    )
  })
) %>%
  dplyr::mutate(
    Target = factor(Target, levels = TARGET_POPULATIONS),
    Weight = ifelse(is.finite(Weight), Weight, 0)
  ) %>%
  dplyr::filter(Weight > 0)

if (nrow(network_table) == 0L) {
  stop(
    "No non-zero outgoing Cluster 0/8-to-TME communication weights were found."
  )
}

bubble_table <- dplyr::bind_rows(
  lapply(names(source_groups), function(cluster_label) {
    source_name <- source_groups[[cluster_label]]
    source_index <- match(source_name, group_names)
    target_indices <- match(available_targets, group_names)

    communication_data <- tryCatch(
      CellChat::subsetCommunication(
        object = cellchat,
        sources.use = source_index,
        targets.use = target_indices
      ),
      error = function(e) NULL
    )

    if (is.null(communication_data) || nrow(communication_data) == 0L) {
      return(NULL)
    }

    d <- as.data.frame(communication_data, stringsAsFactors = FALSE)

    if (!"prob" %in% colnames(d)) return(NULL)
    if (!"pval" %in% colnames(d)) d$pval <- 1

    d$prob <- as.numeric(d$prob)
    d$pval <- as.numeric(d$pval)
    d$Target <- as.character(d$target)
    d$LR_pair <- get_interaction_label(d)

    d <- d %>%
      dplyr::filter(
        Target %in% available_targets,
        is.finite(prob),
        prob > 0,
        is.finite(pval),
        pval < PVAL_CUTOFF
      )

    if (nrow(d) == 0L) return(NULL)

    top_pairs <- d %>%
      dplyr::group_by(LR_pair) %>%
      dplyr::summarise(Total_probability = sum(prob), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(Total_probability)) %>%
      dplyr::slice_head(n = TOP_LR_TO_SHOW) %>%
      dplyr::pull(LR_pair)

    d %>%
      dplyr::filter(LR_pair %in% top_pairs) %>%
      dplyr::mutate(
        Cluster = cluster_label,
        Target = factor(Target, levels = TARGET_POPULATIONS),
        Minus_log10_P = pmin(-log10(pmax(pval, 1e-16)), 16)
      ) %>%
      dplyr::select(
        Cluster, Target, LR_pair, prob, Minus_log10_P, pval
      )
  })
)

if (is.null(bubble_table) || nrow(bubble_table) == 0L) {
  stop(
    "No significant ligand–receptor interactions were retained for the selected clusters."
  )
}

write.csv(
  network_table,
  file.path(OUTPUT_DIR, "Figure13_Clusters0_8_outgoing_network_weights.csv"),
  row.names = FALSE
)

write.csv(
  bubble_table,
  file.path(OUTPUT_DIR, "Figure13_Clusters0_8_top5_LR_pairs.csv"),
  row.names = FALSE
)

probability_limits <- safe_limits(bubble_table$prob)
significance_limits <- safe_limits(bubble_table$Minus_log10_P)
edge_limits <- safe_limits(network_table$Weight)

# 8. MANUSCRIPT-STYLE PLOTTING FUNCTIONS

make_network_plot <- function(cluster_id, panel_tag) {
  source_name <- paste0("Cluster ", cluster_id)
  source_group <- source_groups[[source_name]]

  nodes <- make_network_layout(
    cluster_id = cluster_id,
    available_targets = available_targets
  )

  node_cell_counts <- c(
    setNames(group_size[source_group], source_name),
    group_size[available_targets]
  )

  nodes$Cells <- as.numeric(node_cell_counts[nodes$cell])
  nodes$Cells[!is.finite(nodes$Cells) | nodes$Cells <= 0] <- 1

  edges <- network_table %>%
    dplyr::filter(Cluster == source_name) %>%
    dplyr::mutate(
      Source = source_name,
      Target_label = as.character(Target)
    ) %>%
    dplyr::left_join(
      nodes %>% dplyr::select(cell, x, y),
      by = c("Source" = "cell")
    ) %>%
    dplyr::rename(x_start = x, y_start = y) %>%
    dplyr::left_join(
      nodes %>% dplyr::select(cell, x, y),
      by = c("Target_label" = "cell")
    ) %>%
    dplyr::rename(x_end = x, y_end = y)

  if (nrow(edges) == 0L) {
    return(
      placeholder_plot(
        paste0("No outgoing interactions for ", source_name)
      )
    )
  }

  ggplot2::ggplot() +
    ggplot2::geom_curve(
      data = edges,
      ggplot2::aes(
        x = x_start,
        y = y_start,
        xend = x_end,
        yend = y_end,
        linewidth = Weight
      ),
      curvature = 0.10,
      colour = NETWORK_EDGE_COLOUR,
      alpha = 0.78
    ) +
    ggplot2::scale_linewidth_continuous(
      limits = edge_limits,
      range = c(0.65, 2.85),
      guide = "none"
    ) +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = x, y = y, size = Cells, fill = node_type),
      shape = 21,
      stroke = 0.75,
      colour = "white"
    ) +
    ggplot2::scale_size_continuous(
      range = c(6.7, 13.1),
      guide = "none"
    ) +
    ggplot2::scale_fill_manual(
      values = CELLCHAT_NODE_COLOURS,
      guide = "none"
    ) +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(
        x = x,
        y = y,
        label = cell,
        vjust = label_vjust
      ),
      family = FONT_FAMILY,
      size = 3.55,
      colour = "black"
    ) +
    ggplot2::coord_equal(
      xlim = c(-1.70, 1.70),
      ylim = c(-1.62, 1.86),
      clip = "off"
    ) +
    ggplot2::labs(
      title = paste0(source_name, " outgoing signaling"),
      tag = panel_tag
    ) +
    ggplot2::theme_void(base_family = FONT_FAMILY) +
    ggplot2::theme(
      text = ggplot2::element_text(family = FONT_FAMILY, color = "black"),
      plot.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PLOT_TITLE_PT,
        hjust = 0.5
      ),
      plot.tag = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PANEL_TAG_PT
      ),
      plot.tag.position = c(0.01, 0.99),
      plot.margin = MANUSCRIPT_MARGIN
    )
}

make_bubble_plot <- function(cluster_id) {
  d <- bubble_table %>%
    dplyr::filter(Cluster == paste0("Cluster ", cluster_id))

  if (nrow(d) == 0L) {
    return(
      placeholder_plot(
        paste0("No significant ligand–receptor pairs for Cluster ", cluster_id)
      )
    )
  }

  pair_order <- d %>%
    dplyr::group_by(LR_pair) %>%
    dplyr::summarise(Total_probability = sum(prob), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(Total_probability)) %>%
    dplyr::pull(LR_pair)

  d <- d %>%
    dplyr::mutate(
      LR_pair = factor(LR_pair, levels = rev(pair_order))
    )

  ggplot2::ggplot(
    d,
    ggplot2::aes(x = Target, y = LR_pair)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = Minus_log10_P,
        colour = prob
      ),
      alpha = 0.96
    ) +
    ggplot2::scale_colour_gradientn(
      colours = CELLCHAT_PALETTE,
      limits = probability_limits,
      name = "Communication\nprobability"
    ) +
    ggplot2::scale_size_continuous(
      limits = significance_limits,
      range = c(2.6, 6.8),
      breaks = pretty(significance_limits, n = 4),
      name = "Significance, −log10(P value)"
    ) +
    ggplot2::labs(
      title = paste0(
        "Cluster ", cluster_id,
        ": ligand–receptor pairs"
      ),
      x = "Target cell type",
      y = NULL
    ) +
    theme_manuscript(
      show_grid = FALSE,
      legend_position = "right"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = PLOT_TITLE_PT
      ),
      axis.text.x = ggplot2::element_text(
        angle = 42,
        hjust = 1,
        vjust = 1,
        size = AXIS_TEXT_PT + 1.5
      ),
      axis.text.y = ggplot2::element_text(
        size = AXIS_TEXT_PT
      ),
      axis.title.x = ggplot2::element_text(
        size = AXIS_TITLE_PT
      ),
      legend.title = ggplot2::element_text(
        size = LEGEND_TITLE_PT
      ),
      legend.text = ggplot2::element_text(
        size = LEGEND_TEXT_PT
      ),
      panel.grid.major.x = ggplot2::element_line(
        color = "grey94",
        linewidth = GRID_LWD
      )
    )
}

# 9. BUILD AND EXPORT FIGURE 13

network_0 <- make_network_plot("0", "a")
network_8 <- make_network_plot("8", "b")

bubble_0 <- make_bubble_plot("0")
bubble_8 <- make_bubble_plot("8")

row_0 <- network_0 + bubble_0 + patchwork::plot_layout(widths = c(2.15, 1.85))
row_8 <- network_8 + bubble_8 + patchwork::plot_layout(widths = c(2.15, 1.85))

figure_13_combined <- (row_0 / row_8) +
  patchwork::plot_layout(
    heights = c(1, 1),
    guides = "collect"
  ) &
  ggplot2::theme(
    legend.position = "right"
  )

print(figure_13_combined)

save_plot_all_formats(
  plot_obj = figure_13_combined,
  filename_stem = "Figure13_CellChat_Clusters0_8_EnhancedNetwork_LRpairs",
  dir_path = OUTPUT_DIR,
  width = FIGURE13_COMBINED_W,
  height = FIGURE13_COMBINED_H
)

# 10. RUN SUMMARY

final_png <- file.path(
  OUTPUT_DIR,
  "Figure13_CellChat_Clusters0_8_EnhancedNetwork_LRpairs.png"
)
final_tiff <- file.path(
  OUTPUT_DIR,
  "Figure13_CellChat_Clusters0_8_EnhancedNetwork_LRpairs.tiff"
)
final_pdf <- file.path(
  OUTPUT_DIR,
  "Figure13_CellChat_Clusters0_8_EnhancedNetwork_LRpairs.pdf"
)

summary_lines <- c(
  "CellChat selected-cluster Figure 13 run completed.",
  paste0("CellChat source: ", cellchat_input_description),
  paste0("Stable CellChat RDS: ", STABLE_CELLCHAT_RDS),
  paste0("Output folder: ", OUTPUT_DIR),
  "Final figure includes outgoing communication from Clusters 0 and 8.",
  "Each source row contains a widened radial outgoing network and the selected ligand–receptor pairs.",
  "All probability colour scales, point-size scales, and edge-width scales are shared across clusters.",
  "Bubble-plot x-axis labels were enlarged and the size legend is titled: Significance, −log10(P value).",
  "Figure settings: radial fan layout; muted pastel palette; Arial; 600 dpi; white background.",
  paste0("Final PNG: ", final_png),
  paste0("Final TIFF: ", final_tiff),
  paste0("Final PDF: ", final_pdf)
)

writeLines(
  summary_lines,
  file.path(OUTPUT_DIR, "Figure13_CellChat_Clusters0_8_run_summary.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "Figure13_CellChat_Clusters0_8_sessionInfo.txt")
)

cat(paste(summary_lines, collapse = "\n"), "\n")

# MODULE 05 — MYBL2-continuous Hallmark GSEA analysis

rm(list = ls(all.names = TRUE))
gc()

set.seed(123)

options(
  stringsAsFactors = FALSE,
  scipen = 100,
  timeout = 7200,
  width = 160
)

# 1. VERSION CHECK, PATHS, AND ANALYTICAL SETTINGS

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
  "MYBL2_limma_batch_adjusted_gene_association.csv"
)

HALLMARK_GMT <- file.path(
  OUTPUT_DIR,
  "Hallmark_gene_sets_used_LimmaBatchAdjusted.gmt"
)

HALLMARK_GENESET_CSV <- file.path(
  OUTPUT_DIR,
  "Hallmark_gene_sets_used_LimmaBatchAdjusted.csv"
)

GSEA_RESULTS_CSV <- file.path(
  OUTPUT_DIR,
  "MYBL2_limma_batch_adjusted_Hallmark_GSEA_results.csv"
)

GSEA_LEADING_EDGE_CSV <- file.path(
  OUTPUT_DIR,
  "MYBL2_limma_batch_adjusted_Hallmark_GSEA_leading_edge_genes.csv"
)

FOCUS_PATHWAYS_CSV <- file.path(
  OUTPUT_DIR,
  "MYBL2_limma_batch_adjusted_Hallmark_focus_pathways.csv"
)

TUMOR_SAMPLE_METADATA_CSV <- file.path(
  OUTPUT_DIR,
  "MYBL2_limma_batch_adjusted_tumor_samples.csv"
)

PROVENANCE_TXT <- file.path(
  OUTPUT_DIR,
  "MYBL2_limma_batch_adjusted_input_provenance.txt"
)

SUMMARY_TXT <- file.path(
  OUTPUT_DIR,
  "MYBL2_limma_batch_adjusted_summary.txt"
)

STATUS_TXT <- file.path(
  OUTPUT_DIR,
  "MYBL2_limma_batch_adjusted_status.txt"
)

SESSION_TXT <- file.path(
  OUTPUT_DIR,
  "MYBL2_limma_batch_adjusted_sessionInfo.txt"
)

RUN_LOG <- file.path(
  OUTPUT_DIR,
  "MYBL2_limma_batch_adjusted_run.log"
)

RESULT_RDS <- file.path(
  OUTPUT_DIR,
  "MYBL2_limma_batch_adjusted_analysis_objects.rds"
)

SUMMARY_FIGURE_STEM <- (
  "Figure_MYBL2_LimmaBatchAdjusted_Hallmark_GSEA_Summary_GeneRankStyled"
)

CURVE_FIGURE_STEM <- (
  "Figure_MYBL2_LimmaBatchAdjusted_Hallmark_GSEA_EnrichmentCurves_GeneRankStyled"
)

# 2. FIGURE SETTINGS

FONT_FAMILY <- "Arial"
FIG_DPI <- 600L
FIG_BACKGROUND <- "white"

SUMMARY_FIGURE_W <- 10.60
SUMMARY_FIGURE_H <- 8.20

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

# 3. PACKAGE CHECK AND LOADING

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

# 4. GENERAL HELPER FUNCTIONS

log_message <- function(...) {
  text <- paste0(..., collapse = "")
  message(text)

  cat(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ",
    text,
    "\n",
    file = RUN_LOG,
    append = TRUE,
    sep = ""
  )
}

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
      message(
        "fgseaMultilevel failed; trying fgseaSimple. Reason: ",
        multilevel_error$message
      )

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

# 5. FIGURE HELPERS

make_gsea_summary_plot <- function(gsea_results) {
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

# 6. MAIN ANALYSIS

run_mybl2_limma_batch_adjusted_gsea <- function() {
  writeLines(character(0), RUN_LOG)

  status <- list(
    completed = FALSE,
    stage = "started",
    error = NA_character_
  )

  write_status <- function() {
    writeLines(
      c(
        "MYBL2-continuous Hallmark GSEA status — limma batch-adjusted version",
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
      status$stage <- "locating_input"
      log_message("STAGE: locating pooled discovery CPM input")

      discovery_cpm_file <- find_first_existing(
        DEVELOPMENT_CPM_CANDIDATES,
        "Discovery bulk CPM matrix"
      )

      status$stage <- "reading_tumor_expression"
      log_message("STAGE: reading pooled discovery cohort and retaining tumors")

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

      if (nlevels(tumor_metadata$Batch) < 2L) {
        log_message(
          "NOTICE: only one batch is available; the model will include MYBL2 only."
        )
      } else {
        log_message(
          "Batch-adjusted model enabled with ",
          nlevels(tumor_metadata$Batch),
          " batch levels."
        )
      }

      write.csv(
        tumor_metadata,
        TUMOR_SAMPLE_METADATA_CSV,
        row.names = FALSE
      )

      status$stage <- "filtering_zero_variance_genes"
      log_message("STAGE: removing zero-variance genes before limma")

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

      status$stage <- "limma_continuous_association"
      log_message(
        "STAGE: limma model of continuous MYBL2 expression adjusted for batch"
      )

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

      status$stage <- "retrieving_hallmark_gene_sets"
      log_message("STAGE: retrieving MSigDB Hallmark gene sets")

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

      status$stage <- "running_preranked_gsea"
      log_message("STAGE: pre-ranked Hallmark GSEA")

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

      status$stage <- "creating_figures"
      log_message("STAGE: creating Hallmark GSEA manuscript figures")

      summary_plot <- make_gsea_summary_plot(
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
        summary_plot,
        SUMMARY_FIGURE_STEM,
        width = SUMMARY_FIGURE_W,
        height = SUMMARY_FIGURE_H
      )

      save_plot_all_formats(
        curve_plot,
        CURVE_FIGURE_STEM,
        width = CURVE_FIGURE_W,
        height = CURVE_FIGURE_H
      )

      status$stage <- "writing_reproducibility_outputs"
      log_message("STAGE: writing summary, provenance, and RDS outputs")

      writeLines(
        c(
          "Input provenance",
          paste0("Discovery CPM matrix: ", input$source),
          paste0("Output directory: ", OUTPUT_DIR),
          "Analysis samples: tumor samples from pooled GSE127165 + GSE142083 only.",
          "Expression scale: log2(TMM-CPM + 1).",
          paste0(
            "Batch adjustment used: ",
            ifelse(batch_adjustment_used, "yes", "no; only one batch level")
          ),
          "Model: limma continuous MYBL2 z-score association plus batch covariates when available.",
          "Ranking statistic: moderated limma t-statistic for MYBL2_z.",
          "MYBL2 itself was removed from the ranked list before GSEA.",
          "Genes with zero variance across tumors were removed before limma fitting.",
          "Exact ties in the limma t-statistic were given deterministic negligible offsets only for fgsea ranking.",
          "Pathways: MSigDB Hallmark gene sets retrieved using msigdbr.",
          "GSEA: fgsea pre-ranked analysis.",
          "All fgsea list columns were converted to semicolon-separated strings before CSV export.",
          "Figure labels use canonical pathway capitalization (G2/M, DNA, and MYC). Enrichment panels use a refined schematic style with a new four-colour palette (cobalt blue, crimson red, forest green, and amber), thicker smooth curves, short subtle hit-rank marks, semi-transparent NES/FDR annotation boxes, a compact lower strip without a Leading edge label, and a white background with very light panel borders; original fgsea NES and FDR values are retained. The custom renderer preserves ranked-gene names during numeric conversion."
        ),
        PROVENANCE_TXT
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

      significant_pathway_count <- sum(
        is.finite(gsea_results$padj) &
        gsea_results$padj < FDR_CUTOFF
      )

      writeLines(
        c(
          "MYBL2-continuous Hallmark GSEA completed successfully.",
          "",
          "Analysis population:",
          paste0(
            "  ",
            nrow(tumor_metadata),
            " tumor samples from pooled GSE127165 + GSE142083."
          ),
          "",
          "Continuous association model:",
          paste0(
            "  Gene expression ~ standardized MYBL2 expression",
            ifelse(
              batch_adjustment_used,
              " + batch.",
              "."
            )
          ),
          "  Ranking statistic: moderated limma t-statistic.",
          "  MYBL2 was removed from the ranked gene list before GSEA.",
          "",
          "GSEA:",
          "  Database: MSigDB Hallmark gene sets.",
          paste0(
            "  Significant pathways at BH-FDR < ",
            FDR_CUTOFF,
            ": ",
            significant_pathway_count,
            "."
          ),
          "",
          "Figures:",
          paste0(
            "  ",
            file.path(
              OUTPUT_DIR,
              SUMMARY_FIGURE_STEM
            ),
            ".png/.tiff/.pdf"
          ),
          paste0(
            "  ",
            file.path(
              OUTPUT_DIR,
              CURVE_FIGURE_STEM
            ),
            ".png/.tiff/.pdf"
          )
        ),
        SUMMARY_TXT
      )

      status$completed <- TRUE
      status$stage <- "finished"
      write_status()

      cat(
        "\n============================================================\n",
        "MYBL2-CONTINUOUS LIMMA BATCH-ADJUSTED HALLMARK GSEA FINISHED\n",
        "============================================================\n",
        "All outputs were saved in:\n",
        OUTPUT_DIR,
        "\n============================================================\n",
        sep = ""
      )
    },
    error = function(error_object) {
      status$error <- conditionMessage(error_object)
      status$stage <- "failed"

      write_status()

      writeLines(
        c(
          "MYBL2-continuous limma batch-adjusted Hallmark GSEA failed.",
          "",
          conditionMessage(error_object)
        ),
        file.path(
          OUTPUT_DIR,
          "MYBL2_limma_batch_adjusted_GSEA_error.txt"
        )
      )

      stop(error_object)
    }
  )
}

# 7. EXECUTE

run_mybl2_limma_batch_adjusted_gsea()

# MODULE 06 — scRNA-seq preprocessing, QC, clustering, and inferCNV workflow

# STEP 1. CLEAN ENVIRONMENT

rm(list = ls())
gc()

set.seed(123)
options(stringsAsFactors = FALSE)
options(scipen = 100)

# STEP 2. LIBRARIES

required_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "data.table",
  "dplyr",
  "ggplot2",
  "harmony",
  "infercnv",
  "patchwork",
  "pheatmap",
  "grid",
  "png"
)

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Install these packages first:\n",
    paste(missing_pkgs, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(harmony)
  library(infercnv)
  library(patchwork)
  library(pheatmap)
  library(grid)
  library(png)
})

# STEP 3. PATHS

scrna_root_dir   <- "D:/LSCC/ScRNAseq_Results"
scrna_project_id <- "GSE206332"

scrna_base_dir <- file.path(scrna_root_dir, scrna_project_id)
scrna_rawdir   <- file.path(scrna_base_dir, "raw_extracted")
scrna_out_dir  <- file.path(scrna_base_dir, "Results")

scrna_figdir <- file.path(scrna_out_dir, "figures")
scrna_rdsdir <- file.path(scrna_out_dir, "rds")
infercnv_out <- file.path(scrna_out_dir, "infercnv")

dir.create(scrna_rawdir, recursive = TRUE, showWarnings = FALSE)
dir.create(scrna_figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(scrna_rdsdir, recursive = TRUE, showWarnings = FALSE)
dir.create(infercnv_out, recursive = TRUE, showWarnings = FALSE)

gene_order_file <- "D:/LSCC/ScRNAseq_Results/hg38_gencode_v27.txt"

if (!file.exists(gene_order_file)) {
  stop("gene_order_file not found.")
}

# STEP 4. PARAMETERS

QC_MIN_FEATURES <- 200
QC_MAX_MT <- 20

N_HVG <- 2000
NPCS <- 40
HARMONY_DIMS <- 1:30
RESOLUTION_GLOBAL <- 0.8
RESOLUTION_MALIGNANT <- 0.5

FINDALL_LOGFC <- 0.10
FINDALL_MINPCT <- 0.10
GLOBAL_MARKER_SIG_LOGFC <- 0.25
GLOBAL_MARKER_SIG_PADJ <- 0.05

HIGH_CNV_CLUSTER_MARKER_LOGFC <- 1
HIGH_CNV_CLUSTER_MARKER_MINPCT <- 0.25
HIGH_CNV_CLUSTER_MARKER_PADJ <- 0.05

LOW_CNV_AUTO_QUANTILE <- 0.25
LOW_CNV_MIN_CLUSTERS <- 1
LOW_CNV_MAX_FRACTION <- 0.35

USE_MANUAL_LOW_CNV_CLUSTERS <- FALSE
MANUAL_LOW_CNV_CLUSTERS <- c(
)

INITIAL_INFERCNV_CUTOFF <- 0.1
REFINED_INFERCNV_CUTOFF <- 0.1
MIN_LOW_CNV_REFERENCE_CELLS <- 10

FIG_DPI <- 600
FONT_FAMILY <- "Arial"

# STEP 5. HELPER FUNCTIONS

join_layers_safe <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj2 <- tryCatch(
    Seurat::JoinLayers(object = obj, assay = "RNA"),
    error = function(e) obj
  )
  return(obj2)
}

get_counts_safe <- function(obj) {

  DefaultAssay(obj) <- "RNA"
  assay_obj <- obj[["RNA"]]

  layer_names <- tryCatch(
    Layers(assay_obj),
    error = function(e) character(0)
  )

  count_layers <- layer_names[grepl("^counts", layer_names)]

  if (length(count_layers) == 0) {
    counts_old <- tryCatch(
      GetAssayData(obj, assay = "RNA", slot = "counts"),
      error = function(e) NULL
    )

    if (!is.null(counts_old) && ncol(counts_old) > 0) {
      return(as(counts_old, "dgCMatrix"))
    }

    stop("No counts layer found.")
  }

  count_list <- list()

  for (ly in count_layers) {

    mat <- tryCatch(
      LayerData(obj, assay = "RNA", layer = ly, fast = FALSE),
      error = function(e) NULL
    )

    if (is.null(mat)) next

    cell_overlap_col <- length(intersect(colnames(mat), colnames(obj)))
    cell_overlap_row <- length(intersect(rownames(mat), colnames(obj)))

    if (cell_overlap_row > cell_overlap_col) {
      mat <- t(mat)
    }

    keep_cells <- intersect(colnames(mat), colnames(obj))

    if (length(keep_cells) > 0) {
      mat <- mat[, keep_cells, drop = FALSE]
      count_list[[ly]] <- mat
    }
  }

  if (length(count_list) == 0) {
    stop("No usable count layer found.")
  }

  common_genes <- Reduce(intersect, lapply(count_list, rownames))

  count_list <- lapply(count_list, function(m) {
    m[common_genes, , drop = FALSE]
  })

  counts <- do.call(cbind, count_list)
  counts <- counts[, !duplicated(colnames(counts)), drop = FALSE]

  wanted_cells <- intersect(colnames(obj), colnames(counts))
  counts <- counts[, wanted_cells, drop = FALSE]

  return(as(counts, "dgCMatrix"))
}

make_DE_assay <- function(obj, assay_name = "RNA_DE") {

  DefaultAssay(obj) <- "RNA"

  counts <- get_counts_safe(obj)

  if (assay_name %in% names(obj@assays)) {
    obj[[assay_name]] <- NULL
  }

  obj[[assay_name]] <- CreateAssayObject(counts = counts)
  DefaultAssay(obj) <- assay_name

  obj <- NormalizeData(
    obj,
    assay = assay_name,
    normalization.method = "LogNormalize",
    verbose = FALSE
  )

  obj <- tryCatch(
    JoinLayers(obj, assay = assay_name),
    error = function(e) obj
  )

  DefaultAssay(obj) <- assay_name

  return(obj)
}

save_plot <- function(plot_obj, filename, width = 8, height = 6) {
  ggplot2::ggsave(
    filename = file.path(scrna_figdir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = FIG_DPI,
    bg = "white",
    limitsize = FALSE
  )
}

make_clean_theme <- function() {
  theme_bw(base_size = 12, base_family = FONT_FAMILY) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.45),
      panel.grid.minor = element_blank()
    )
}

extract_cluster_num <- function(x) {
  as.integer(gsub("[^0-9]+", "", as.character(x)))
}

make_cluster_factor <- function(x) {
  nums <- sort(unique(extract_cluster_num(x)))
  factor(as.character(x), levels = paste0("Cluster ", nums))
}

# STEP 6. SAMPLE INFO

sample_info <- data.frame(
  gsm = c("GSM6251294", "GSM6251297", "GSM6251300"),
  sample = c("LSCC1", "LSCC2", "LSCC3"),
  stringsAsFactors = FALSE
)

# STEP 7. LOAD RAW FILES

if (length(list.files(scrna_rawdir, recursive = TRUE, full.names = TRUE)) == 0) {
  tar_file <- file.choose()
  untar(tar_file, exdir = scrna_rawdir)
}

raw_files <- list.files(scrna_rawdir, recursive = TRUE, full.names = TRUE)

# STEP 8. CREATE SEURAT OBJECTS AND APPLY BASIC QC FILTERS

seu_list <- list()
qc_report <- list()

for (i in 1:nrow(sample_info)) {

  gsm <- sample_info$gsm[i]

  mtx <- raw_files[grepl(paste0("^", gsm, ".*_matrix"), basename(raw_files))]
  bc  <- raw_files[grepl(paste0("^", gsm, ".*_barcodes"), basename(raw_files))]
  fea <- raw_files[grepl(paste0("^", gsm, ".*_features"), basename(raw_files))]

  if (length(mtx) == 0 | length(bc) == 0 | length(fea) == 0) {
    stop("Missing 10x files for sample: ", gsm)
  }

  mat <- Matrix::readMM(gzfile(mtx[1]))
  mat <- as(mat, "CsparseMatrix")

  features <- read.delim(gzfile(fea[1]), header = FALSE)
  barcodes <- read.delim(gzfile(bc[1]), header = FALSE)

  rownames(mat) <- make.unique(features[[2]])
  colnames(mat) <- paste0(gsm, "_", barcodes[[1]])

  seu <- CreateSeuratObject(counts = mat)
  DefaultAssay(seu) <- "RNA"

  seu$sample <- sample_info$sample[i]
  seu$gsm <- gsm
  seu$percent.mt <- PercentageFeatureSet(seu, pattern = "^MT-")

  cells_before <- ncol(seu)

  seu <- subset(
    seu,
    subset = nFeature_RNA >= QC_MIN_FEATURES &
      percent.mt < QC_MAX_MT
  )

  cells_after <- ncol(seu)

  qc_report[[gsm]] <- data.frame(
    GSM = gsm,
    Sample = sample_info$sample[i],
    Cells_before_QC = cells_before,
    Cells_after_QC = cells_after
  )

  seu_list[[gsm]] <- seu
}

qc_report_df <- bind_rows(qc_report)

write.csv(
  qc_report_df,
  file.path(scrna_figdir, "Simple_scRNA_QC_Report.csv"),
  row.names = FALSE
)

cat("\nQC finished.\n")
print(qc_report_df)

# STEP 9. MERGE, NORMALIZE, HVG, PCA, HARMONY, CLUSTERING, t-SNE

seurat_obj <- Reduce(merge, seu_list)
DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- join_layers_safe(seurat_obj)

seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize"
)

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = N_HVG
)

seurat_obj <- ScaleData(seurat_obj)

seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(seurat_obj),
  npcs = NPCS
)

seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars = "sample",
  dims.use = HARMONY_DIMS
)

seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "harmony",
  dims = HARMONY_DIMS
)

seurat_obj <- FindClusters(
  seurat_obj,
  resolution = RESOLUTION_GLOBAL
)

seurat_obj <- RunTSNE(
  seurat_obj,
  reduction = "harmony",
  dims = HARMONY_DIMS
)

p_tsne_clusters <- DimPlot(
  seurat_obj,
  reduction = "tsne",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  make_clean_theme() +
  labs(x = "t-SNE 1", y = "t-SNE 2", title = NULL)

save_plot(
  p_tsne_clusters,
  "Figure_SC_01A_tSNE_global_clusters.png",
  width = 8,
  height = 6
)

# STEP 10. FINDALLMARKERS FOR GLOBAL CLUSTERS

seurat_de_obj <- make_DE_assay(
  seurat_obj,
  assay_name = "RNA_DE_GLOBAL"
)

Idents(seurat_de_obj) <- "seurat_clusters"

global_cluster_markers_all <- FindAllMarkers(
  seurat_de_obj,
  assay = "RNA_DE_GLOBAL",
  logfc.threshold = FINDALL_LOGFC,
  min.pct = FINDALL_MINPCT,
  return.thresh = 1,
  only.pos = FALSE
)

if (nrow(global_cluster_markers_all) > 0) {

  global_cluster_markers_all$Gene <- rownames(global_cluster_markers_all)

  fc_col_global <- if ("avg_log2FC" %in% colnames(global_cluster_markers_all)) {
    "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(global_cluster_markers_all)) {
    "avg_logFC"
  } else {
    NA_character_
  }

  if (!is.na(fc_col_global)) {
    global_cluster_markers_all <- global_cluster_markers_all %>%
      dplyr::arrange(
        cluster,
        p_val_adj,
        dplyr::desc(abs(.data[[fc_col_global]]))
      )
  }

} else {

  warning("Global FindAllMarkers returned zero rows even after creating RNA_DE_GLOBAL.")
  global_cluster_markers_all$Gene <- character(0)
  fc_col_global <- NA_character_
}

write.csv(
  global_cluster_markers_all,
  file.path(scrna_figdir, "Global_Cluster_FindAllMarkers_all_logFC0.10_minpct0.10.csv"),
  row.names = FALSE
)

write.csv(
  global_cluster_markers_all,
  file.path(scrna_figdir, "Global_Cluster_FindAllMarkers_logFC0.25_minpct0.25.csv"),
  row.names = FALSE
)

if (nrow(global_cluster_markers_all) > 0 && !is.na(fc_col_global)) {

  global_cluster_markers_sig <- global_cluster_markers_all %>%
    dplyr::filter(
      abs(.data[[fc_col_global]]) > GLOBAL_MARKER_SIG_LOGFC,
      p_val_adj < GLOBAL_MARKER_SIG_PADJ
    )

} else {

  global_cluster_markers_sig <- global_cluster_markers_all[0, , drop = FALSE]
}

write.csv(
  global_cluster_markers_sig,
  file.path(scrna_figdir, "Global_Cluster_FindAllMarkers_significant.csv"),
  row.names = FALSE
)

cat("\nGlobal cluster marker detection finished.\n")
cat("Global all marker rows:", nrow(global_cluster_markers_all), "\n")
cat("Global significant marker rows:", nrow(global_cluster_markers_sig), "\n")

rm(seurat_de_obj)
gc()

# STEP 11. SIMPLE CELL TYPE ANNOTATION USING CANONICAL MARKERS

marker_sets <- list(
  Bcell       = c("CD79A", "CD19", "CD79B"),
  Endothelial = c("KDR", "FLT1", "TEK", "ICAM1"),
  Epithelial  = c("EPCAM", "KRT15", "KRT18", "KRT19"),
  Fibroblast  = c("ACTA2", "FAP", "S100A4"),
  Myeloid     = c("CD68", "CD33", "CD1E", "LYZ", "LAMP3"),
  NKcell      = c("NCAM1", "FCGR3A", "NCR1", "NCR3"),
  Tcell       = c("CD2", "CD3D", "CD3E", "CD3G")
)

marker_sets <- lapply(marker_sets, function(x) {
  intersect(x, rownames(seurat_obj))
})

for (nm in names(marker_sets)) {
  if (length(marker_sets[[nm]]) > 0) {
    seurat_obj <- AddModuleScore(
      seurat_obj,
      features = list(marker_sets[[nm]]),
      name = paste0("MS_", nm, "_")
    )
  }
}

score_cols <- grep("^MS_.*_1$", colnames(seurat_obj@meta.data), value = TRUE)
score_mat <- seurat_obj@meta.data[, score_cols, drop = FALSE]
labels <- sub("^MS_(.*)_1$", "\\1", score_cols)

seurat_obj$celltype <- apply(score_mat, 1, function(x) {
  if (all(is.na(x))) return("Unknown")
  labels[which.max(x)]
})

write.csv(
  data.frame(
    Cell = colnames(seurat_obj),
    Sample = seurat_obj$sample,
    Cluster = seurat_obj$seurat_clusters,
    Celltype = seurat_obj$celltype
  ),
  file.path(scrna_figdir, "Simple_Celltype_Annotation.csv"),
  row.names = FALSE
)

p_tsne_celltype <- DimPlot(
  seurat_obj,
  reduction = "tsne",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE
) +
  make_clean_theme() +
  labs(x = "t-SNE 1", y = "t-SNE 2", title = NULL)

save_plot(
  p_tsne_celltype,
  "Figure_SC_01B_tSNE_initial_celltypes.png",
  width = 8,
  height = 6
)

cat("\nCell type annotation finished.\n")
print(table(seurat_obj$celltype))

# STEP 12. INITIAL inferCNV AND CNV SCORE FOR MALIGNANT CELL DETECTION

samples_use <- unique(seurat_obj$sample)

seurat_obj$CNV_score <- NA_real_
seurat_obj$CNV_cutoff <- NA_real_
seurat_obj$malignant <- "non_malignant"

cnv_report <- list()

for (s in samples_use) {

  obj_s <- subset(
    seurat_obj,
    subset = sample == s & celltype %in% c("Epithelial", "Myeloid")
  )

  obj_s <- join_layers_safe(obj_s)

  if (ncol(obj_s) == 0) {
    warning("No epithelial/myeloid cells found for sample: ", s)
    next
  }

  counts <- get_counts_safe(obj_s)

  anno <- data.frame(
    cell = colnames(obj_s),
    type = obj_s$celltype
  )

  anno <- anno[anno$cell %in% colnames(counts), , drop = FALSE]
  anno <- anno[match(colnames(counts), anno$cell), , drop = FALSE]

  if (!"Myeloid" %in% anno$type) {
    warning("No Myeloid reference cells for sample: ", s)
    next
  }

  if (!"Epithelial" %in% anno$type) {
    warning("No Epithelial cells for sample: ", s)
    next
  }

  dir_s <- file.path(infercnv_out, paste0("Initial_", s))
  dir.create(dir_s, recursive = TRUE, showWarnings = FALSE)

  anno_file <- file.path(dir_s, "anno.txt")

  write.table(
    anno,
    file = anno_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  infer_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts,
    annotations_file = anno_file,
    gene_order_file = gene_order_file,
    ref_group_names = c("Myeloid")
  )

  infer_obj <- infercnv::run(
    infer_obj,
    cutoff = INITIAL_INFERCNV_CUTOFF,
    out_dir = dir_s,
    denoise = TRUE,
    HMM = FALSE,
    cluster_by_groups = TRUE,
    plot_steps = FALSE,
    output_format = "png"
  )

  rds_path <- file.path(dir_s, "run.final.infercnv_obj")

  if (!file.exists(rds_path)) {
    warning("inferCNV final object not found for sample: ", s)
    next
  }

  infer_final <- readRDS(rds_path)
  cnv_expr <- infer_final@expr.data

  cnv_score_s <- apply(
    cnv_expr,
    2,
    function(x) mean(abs(x - 1), na.rm = TRUE)
  )

  ref_cells <- anno$cell[anno$type == "Myeloid"]
  ref_cells <- intersect(ref_cells, names(cnv_score_s))

  epithelial_cells <- anno$cell[anno$type == "Epithelial"]
  epithelial_cells <- intersect(epithelial_cells, names(cnv_score_s))

  cnv_cutoff_s <- as.numeric(
    quantile(cnv_score_s[ref_cells], 0.95, na.rm = TRUE)
  )

  malignant_cells_s <- epithelial_cells[
    cnv_score_s[epithelial_cells] > cnv_cutoff_s
  ]

  common_cells <- intersect(names(cnv_score_s), colnames(seurat_obj))

  seurat_obj$CNV_score[common_cells] <- cnv_score_s[common_cells]
  seurat_obj$CNV_cutoff[common_cells] <- cnv_cutoff_s

  seurat_obj$malignant[
    colnames(seurat_obj) %in% malignant_cells_s
  ] <- "malignant"

  cnv_report[[s]] <- data.frame(
    Sample = s,
    Reference_Myeloid_cells = length(ref_cells),
    Epithelial_cells = length(epithelial_cells),
    Malignant_cells = length(malignant_cells_s),
    CNV_cutoff = cnv_cutoff_s
  )
}

cnv_report_df <- bind_rows(cnv_report)

write.csv(
  cnv_report_df,
  file.path(scrna_figdir, "Initial_inferCNV_CNV_score_report.csv"),
  row.names = FALSE
)

cat("\nInitial inferCNV finished.\n")
print(cnv_report_df)

# STEP 13. FINAL CELL TYPE LABEL WITH MALIGNANT CELLS

seurat_obj$final_celltype <- seurat_obj$celltype
seurat_obj$final_celltype[seurat_obj$malignant == "malignant"] <- "Malignant"

celltype_order <- c(
  "Bcell",
  "Endothelial",
  "Epithelial",
  "Fibroblast",
  "Malignant",
  "Myeloid",
  "NKcell",
  "Tcell",
  "Unknown"
)

seurat_obj$final_celltype <- factor(
  seurat_obj$final_celltype,
  levels = intersect(celltype_order, unique(as.character(seurat_obj$final_celltype)))
)

p_tsne_final <- DimPlot(
  seurat_obj,
  reduction = "tsne",
  group.by = "final_celltype",
  label = TRUE,
  repel = TRUE
) +
  make_clean_theme() +
  labs(x = "t-SNE 1", y = "t-SNE 2", title = NULL)

save_plot(
  p_tsne_final,
  "Figure_SC_01B_tSNE_final_celltypes_with_malignant.png",
  width = 8,
  height = 6
)

cat("\nFinal cell type counts:\n")
print(table(seurat_obj$final_celltype))

# STEP 14. DOTPLOT OF CANONICAL CELL-TYPE MARKERS

dotplot_markers <- c(
  "KRT18", "EPCAM", "KRT15",
  "FLT1", "FCGR3A", "CD68", "LYZ",
  "CD3E", "CD3D", "CD2",
  "S100A4", "CD79B", "CD19", "CD79A"
)

dotplot_markers <- unique(dotplot_markers)
dotplot_markers <- intersect(dotplot_markers, rownames(seurat_obj))

Idents(seurat_obj) <- "final_celltype"

if (length(dotplot_markers) >= 2) {

  p_dotplot <- DotPlot(
    seurat_obj,
    features = dotplot_markers,
    group.by = "final_celltype"
  ) +
    coord_flip() +
    make_clean_theme() +
    labs(x = NULL, y = NULL, title = NULL) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      legend.position = "right"
    )

  save_plot(
    p_dotplot,
    "Figure_SC_01C_DotPlot_celltype_markers.png",
    width = 12,
    height = 6
  )

} else {

  p_dotplot <- ggplot() +
    annotate("text", x = 0, y = 0, label = "Not enough marker genes for DotPlot") +
    theme_void()
}

# STEP 15. RE-CLUSTERING OF MALIGNANT EPITHELIAL CELLS

malignant_obj <- subset(seurat_obj, subset = malignant == "malignant")

if (ncol(malignant_obj) < 2) {
  stop("Not enough malignant cells found.")
}

malignant_obj <- join_layers_safe(malignant_obj)

malignant_obj <- NormalizeData(
  malignant_obj,
  normalization.method = "LogNormalize"
)

malignant_obj <- FindVariableFeatures(
  malignant_obj,
  selection.method = "vst",
  nfeatures = N_HVG
)

malignant_obj <- ScaleData(malignant_obj)

malignant_obj <- RunPCA(
  malignant_obj,
  features = VariableFeatures(malignant_obj),
  npcs = NPCS
)

dims_use_malig <- 1:min(30, ncol(Embeddings(malignant_obj, "pca")))

malignant_obj <- FindNeighbors(
  malignant_obj,
  dims = dims_use_malig
)

malignant_obj <- FindClusters(
  malignant_obj,
  resolution = RESOLUTION_MALIGNANT
)

malignant_obj <- RunTSNE(
  malignant_obj,
  dims = dims_use_malig
)

malignant_obj$Cluster_label <- make_cluster_factor(
  paste0("Cluster ", as.character(malignant_obj$seurat_clusters))
)

cluster_levels_num <- levels(malignant_obj$Cluster_label)

p_malig_sub <- DimPlot(
  malignant_obj,
  reduction = "tsne",
  group.by = "Cluster_label",
  label = TRUE,
  repel = TRUE
) +
  make_clean_theme() +
  labs(
    x = "t-SNE 1",
    y = "t-SNE 2",
    title = NULL,
    color = NULL
  )

save_plot(
  p_malig_sub,
  "Figure_SC_01D_tSNE_malignant_subclusters.png",
  width = 8,
  height = 6
)

cnv_cluster_summary_initial <- malignant_obj@meta.data %>%
  group_by(Cluster_label) %>%
  summarise(
    Cells = n(),
    Mean_initial_CNV_score = mean(CNV_score, na.rm = TRUE),
    Median_initial_CNV_score = median(CNV_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Cluster_num = extract_cluster_num(Cluster_label)) %>%
  arrange(Cluster_num)

write.csv(
  cnv_cluster_summary_initial,
  file.path(scrna_figdir, "Initial_Malignant_subcluster_CNV_score_summary.csv"),
  row.names = FALSE
)

cat("\nInitial malignant subcluster CNV score summary:\n")
print(cnv_cluster_summary_initial)

# STEP 16. INITIAL CNV VIOLIN + WEAK LOW-CNV REFERENCE SELECTION + REFINED inferCNV

# STEP 17. INITIAL CNV SCORE VIOLIN PLOT FOR REFERENCE SELECTION

p_initial_cnv_violin <- ggplot(
  malignant_obj@meta.data,
  aes(
    x = Cluster_label,
    y = CNV_score,
    fill = Cluster_label
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    alpha = 0.85,
    linewidth = 0.25
  ) +
  geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.85,
    linewidth = 0.25
  ) +
  make_clean_theme() +
  labs(
    x = "Malignant subclusters",
    y = "Initial CNV score",
    title = NULL,
    fill = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "right"
  )

save_plot(
  p_initial_cnv_violin,
  "Figure_SC_INITIAL_CNV_score_violin_for_reference_selection.png",
  width = 8,
  height = 6
)

save_plot(
  p_initial_cnv_violin,
  "Figure_SC_01E_Initial_CNV_score_violin_for_reference_selection.png",
  width = 8,
  height = 6
)

# STEP 18. SELECT WEAK / LOW-CNV CLUSTERS FROM INITIAL CNV SCORE

cnv_cluster_summary_initial_ranked <- cnv_cluster_summary_initial %>%
  dplyr::arrange(Median_initial_CNV_score, Cluster_num)

n_malignant_clusters <- nrow(cnv_cluster_summary_initial_ranked)

if (n_malignant_clusters < 2) {
  stop("At least two malignant subclusters are required for low-CNV reference selection.")
}

low_cnv_auto_cutoff <- as.numeric(
  quantile(
    cnv_cluster_summary_initial_ranked$Median_initial_CNV_score,
    probs = LOW_CNV_AUTO_QUANTILE,
    na.rm = TRUE
  )
)

max_low_cnv_clusters <- floor(n_malignant_clusters * LOW_CNV_MAX_FRACTION)
max_low_cnv_clusters <- max(LOW_CNV_MIN_CLUSTERS, max_low_cnv_clusters)
max_low_cnv_clusters <- min(max_low_cnv_clusters, n_malignant_clusters - 1)

low_cnv_selection_table <- cnv_cluster_summary_initial_ranked %>%
  dplyr::mutate(
    Selection_basis = "Initial CNV score before refined inferCNV",
    Auto_low_CNV_cutoff = low_cnv_auto_cutoff,
    Initially_selected_as_weak_low_CNV =
      Median_initial_CNV_score <= low_cnv_auto_cutoff
  ) %>%
  dplyr::arrange(Median_initial_CNV_score)

if (USE_MANUAL_LOW_CNV_CLUSTERS) {

  if (length(MANUAL_LOW_CNV_CLUSTERS) < 1) {
    stop(
      "USE_MANUAL_LOW_CNV_CLUSTERS is TRUE, but MANUAL_LOW_CNV_CLUSTERS is empty. ",
      "Add cluster labels such as: c('Cluster 0', 'Cluster 2')."
    )
  }

  missing_manual_clusters <- setdiff(
    MANUAL_LOW_CNV_CLUSTERS,
    as.character(cnv_cluster_summary_initial_ranked$Cluster_label)
  )

  if (length(missing_manual_clusters) > 0) {
    stop(
      "These manual low-CNV clusters were not found in malignant clusters: ",
      paste(missing_manual_clusters, collapse = ", ")
    )
  }

  low_cnv_clusters <- MANUAL_LOW_CNV_CLUSTERS
  selection_mode_used <- "Manual selection after inspecting the initial CNV violin plot"

} else {

  low_cnv_clusters <- low_cnv_selection_table %>%
    dplyr::filter(Initially_selected_as_weak_low_CNV) %>%
    dplyr::slice_head(n = max_low_cnv_clusters) %>%
    dplyr::pull(Cluster_label) %>%
    as.character() %>%
    unique()

  if (length(low_cnv_clusters) < LOW_CNV_MIN_CLUSTERS) {
    low_cnv_clusters <- cnv_cluster_summary_initial_ranked %>%
      dplyr::slice_head(n = LOW_CNV_MIN_CLUSTERS) %>%
      dplyr::pull(Cluster_label) %>%
      as.character() %>%
      unique()
  }

  selection_mode_used <- "Automatic selection from the initial CNV-score table"
}

low_cnv_clusters <- cluster_levels_num[cluster_levels_num %in% low_cnv_clusters]

low_cnv_reference_cells <- colnames(malignant_obj)[
  malignant_obj$Cluster_label %in% low_cnv_clusters
]

if (length(low_cnv_reference_cells) < MIN_LOW_CNV_REFERENCE_CELLS) {

  ordered_clusters <- as.character(cnv_cluster_summary_initial_ranked$Cluster_label)

  for (cl in ordered_clusters) {

    if (!(cl %in% low_cnv_clusters)) {
      low_cnv_clusters <- unique(c(low_cnv_clusters, cl))
    }

    low_cnv_reference_cells <- colnames(malignant_obj)[
      malignant_obj$Cluster_label %in% low_cnv_clusters
    ]

    if (length(low_cnv_reference_cells) >= MIN_LOW_CNV_REFERENCE_CELLS) {
      break
    }

    if (!USE_MANUAL_LOW_CNV_CLUSTERS && length(low_cnv_clusters) >= max_low_cnv_clusters) {
      break
    }
  }
}

low_cnv_selection_table <- low_cnv_selection_table %>%
  dplyr::mutate(
    Selection_mode_used = selection_mode_used,
    Final_selected_as_weak_low_CNV_reference =
      as.character(Cluster_label) %in% low_cnv_clusters
  )

weak_low_cnv_reference_clusters_df <- low_cnv_selection_table %>%
  dplyr::filter(Final_selected_as_weak_low_CNV_reference) %>%
  dplyr::select(
    Cluster_label,
    Cells,
    Mean_initial_CNV_score,
    Median_initial_CNV_score,
    Selection_basis,
    Selection_mode_used
  )

write.csv(
  low_cnv_selection_table,
  file.path(scrna_figdir, "Initial_CNV_weak_low_CNV_reference_selection_table.csv"),
  row.names = FALSE
)

write.csv(
  weak_low_cnv_reference_clusters_df,
  file.path(scrna_figdir, "Initial_CNV_weak_low_CNV_reference_clusters.csv"),
  row.names = FALSE
)

write.csv(
  low_cnv_selection_table,
  file.path(scrna_figdir, "AUTO_low_CNV_cluster_selection_table.csv"),
  row.names = FALSE
)

if (length(low_cnv_reference_cells) < MIN_LOW_CNV_REFERENCE_CELLS) {
  stop(
    "Too few weak/low-CNV reference cells after selection. Found: ",
    length(low_cnv_reference_cells)
  )
}

cat("\n###############################################################################\n")
cat("WEAK / LOW-CNV REFERENCE SELECTION FROM INITIAL CNV SCORE\n")
cat("###############################################################################\n")
cat("Selection basis: Initial CNV score before refined inferCNV\n")
cat("Selection mode: ", selection_mode_used, "\n", sep = "")
cat("Initial CNV violin plot saved in:\n")
cat(file.path(scrna_figdir, "Figure_SC_INITIAL_CNV_score_violin_for_reference_selection.png"), "\n")
cat("\nSelected weak/low-CNV malignant clusters used as refined inferCNV reference:\n")
print(low_cnv_clusters)
cat("\nNumber of weak/low-CNV reference cells:\n")
cat(length(low_cnv_reference_cells), "\n")
cat("\nReference selection CSV:\n")
cat(file.path(scrna_figdir, "Initial_CNV_weak_low_CNV_reference_clusters.csv"), "\n")
cat("###############################################################################\n")

# STEP 19. RUN REFINED inferCNV USING INITIAL WEAK LOW-CNV CLUSTERS AS REFERENCE

malignant_counts <- get_counts_safe(malignant_obj)

anno_refined <- data.frame(
  cell = colnames(malignant_obj),
  type = as.character(malignant_obj$Cluster_label),
  stringsAsFactors = FALSE
)

anno_refined <- anno_refined[
  anno_refined$cell %in% colnames(malignant_counts),
  ,
  drop = FALSE
]

anno_refined <- anno_refined[
  match(colnames(malignant_counts), anno_refined$cell),
  ,
  drop = FALSE
]

anno_refined$type <- factor(anno_refined$type, levels = cluster_levels_num)
anno_refined <- anno_refined %>%
  dplyr::arrange(type)
anno_refined$type <- as.character(anno_refined$type)

refined_infercnv_dir <- file.path(infercnv_out, "Refined_malignant_initial_lowCNV_reference")
dir.create(refined_infercnv_dir, recursive = TRUE, showWarnings = FALSE)

anno_refined_file <- file.path(refined_infercnv_dir, "anno_refined_malignant.txt")

write.table(
  anno_refined,
  file = anno_refined_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

refined_infer_obj <- CreateInfercnvObject(
  raw_counts_matrix = malignant_counts,
  annotations_file = anno_refined_file,
  gene_order_file = gene_order_file,
  ref_group_names = low_cnv_clusters
)

refined_infer_obj <- infercnv::run(
  refined_infer_obj,
  cutoff = REFINED_INFERCNV_CUTOFF,
  out_dir = refined_infercnv_dir,
  denoise = TRUE,
  HMM = FALSE,
  cluster_by_groups = TRUE,
  plot_steps = FALSE,
  output_format = "png"
)

refined_rds_path <- file.path(refined_infercnv_dir, "run.final.infercnv_obj")

if (!file.exists(refined_rds_path)) {
  stop("Refined inferCNV final object was not found.")
}

refined_final <- readRDS(refined_rds_path)
refined_cnv_expr <- refined_final@expr.data

refined_cnv_score <- apply(
  refined_cnv_expr,
  2,
  function(x) mean(abs(x - 1), na.rm = TRUE)
)

malignant_obj$Refined_CNV_score <- NA_real_

common_refined_cells <- intersect(
  names(refined_cnv_score),
  colnames(malignant_obj)
)

malignant_obj$Refined_CNV_score[common_refined_cells] <- refined_cnv_score[common_refined_cells]

malignant_obj$CNV_class <- ifelse(
  malignant_obj$Cluster_label %in% low_cnv_clusters,
  "Low-CNV malignant",
  "High-CNV malignant"
)

malignant_obj$CNV_class <- factor(
  malignant_obj$CNV_class,
  levels = c("Low-CNV malignant", "High-CNV malignant")
)

refined_cnv_cluster_summary <- malignant_obj@meta.data %>%
  group_by(Cluster_label, CNV_class) %>%
  summarise(
    Cells = n(),
    Mean_refined_CNV_score = mean(Refined_CNV_score, na.rm = TRUE),
    Median_refined_CNV_score = median(Refined_CNV_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Cluster_num = extract_cluster_num(Cluster_label)) %>%
  arrange(Cluster_num)

write.csv(
  refined_cnv_cluster_summary,
  file.path(scrna_figdir, "Refined_Malignant_subcluster_CNV_score_summary.csv"),
  row.names = FALSE
)

write.csv(
  malignant_obj@meta.data,
  file.path(scrna_rdsdir, "Refined_all_malignant_metadata_with_CNV_class.csv"),
  row.names = TRUE
)

cat("\nRefined malignant subcluster CNV score summary:\n")
print(refined_cnv_cluster_summary)

# STEP 20. VIOLIN PLOT OF REFINED CNV SCORES

p_refined_cnv_violin <- ggplot(
  malignant_obj@meta.data,
  aes(
    x = Cluster_label,
    y = Refined_CNV_score,
    fill = Cluster_label
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    alpha = 0.85,
    linewidth = 0.25
  ) +
  geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.85,
    linewidth = 0.25
  ) +
  make_clean_theme() +
  labs(
    x = "Malignant subclusters",
    y = "Refined CNV score",
    title = NULL,
    fill = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "right"
  )

save_plot(
  p_refined_cnv_violin,
  "Figure_SC_01E_Refined_CNV_score_violin_malignant_subclusters.png",
  width = 8,
  height = 6
)

save_plot(
  p_refined_cnv_violin,
  "Figure_SC_01G_Refined_CNV_score_violin_malignant_subclusters.png",
  width = 8,
  height = 6
)

# STEP 21. USE ORIGINAL inferCNV OUTPUT IMAGE AS PANEL F

old_custom_heatmap_file <- file.path(
  scrna_figdir,
  "Figure_SC_01F_inferCNV.png"
)

if (file.exists(old_custom_heatmap_file)) {
  unlink(old_custom_heatmap_file, force = TRUE)
}

infercnv_png_candidates <- c(
  file.path(refined_infercnv_dir, "infercnv.png"),
  list.files(
    refined_infercnv_dir,
    pattern = "infercnv.*\\.png$",
    full.names = TRUE,
    recursive = TRUE
  ),
  list.files(
    refined_infercnv_dir,
    pattern = "\\.png$",
    full.names = TRUE,
    recursive = TRUE
  )
)

infercnv_png_candidates <- unique(infercnv_png_candidates)
infercnv_png_candidates <- infercnv_png_candidates[file.exists(infercnv_png_candidates)]

if (length(infercnv_png_candidates) < 1) {
  stop(
    "No inferCNV PNG output was found in refined inferCNV folder:\n",
    refined_infercnv_dir
  )
}

infercnv_png_source <- infercnv_png_candidates[1]

infercnv_png_final <- file.path(
  scrna_figdir,
  "Figure_SC_01F_inferCNV.png"
)

file.copy(
  from = infercnv_png_source,
  to = infercnv_png_final,
  overwrite = TRUE
)

infercnv_img <- png::readPNG(infercnv_png_source)

p_infercnv_panel <- ggplot() +
  annotation_custom(
    grob = grid::rasterGrob(
      infercnv_img,
      width = grid::unit(1, "npc"),
      height = grid::unit(1, "npc"),
      interpolate = TRUE
    ),
    xmin = 0,
    xmax = 1,
    ymin = 0,
    ymax = 1
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  labs(title = NULL)

cat("\nOriginal inferCNV image used as panel F:\n")
cat(infercnv_png_source, "\n")

cat("\nCopied inferCNV panel file:\n")
cat(infercnv_png_final, "\n")

# STEP 22. REMOVE LOW-CNV MALIGNANT CLUSTERS

high_cnv_malignant_obj <- subset(
  malignant_obj,
  subset = CNV_class == "High-CNV malignant"
)

write.csv(
  malignant_obj@meta.data,
  file.path(scrna_rdsdir, "All_malignant_metadata_before_lowCNV_removal.csv"),
  row.names = TRUE
)

write.csv(
  high_cnv_malignant_obj@meta.data,
  file.path(scrna_rdsdir, "HighCNV_malignant_metadata_after_lowCNV_removal.csv"),
  row.names = TRUE
)

cat("\nLow-CNV clusters removed:\n")
print(low_cnv_clusters)

cat("\nCells before low-CNV removal:", ncol(malignant_obj), "\n")
cat("Cells after low-CNV removal:", ncol(high_cnv_malignant_obj), "\n")

# STEP 23. HIGH-CNV CLUSTER-SPECIFIC MARKER DETECTION AFTER LOW-CNV REMOVAL

needed_objects <- c(
  "high_cnv_malignant_obj",
  "low_cnv_clusters",
  "scrna_figdir",
  "scrna_rdsdir"
)

missing_objects <- needed_objects[!sapply(needed_objects, exists)]

if (length(missing_objects) > 0) {
  stop(
    "These objects are missing from the R environment:\n",
    paste(missing_objects, collapse = ", "),
    "\n\nRun the single-cell code until the end of STEP 22 first."
  )
}

# STEP 24. Prepare High-CNV malignant object for DE analysis

high_cnv_malignant_obj$Cluster_label <- droplevels(
  high_cnv_malignant_obj$Cluster_label
)

high_cnv_de_obj <- make_DE_assay(
  high_cnv_malignant_obj,
  assay_name = "RNA_DE_HIGH_CNV"
)

Idents(high_cnv_de_obj) <- "Cluster_label"

cat("\nHigh-CNV malignant clusters used for FindAllMarkers:\n")
print(table(Idents(high_cnv_de_obj)))

if (length(levels(Idents(high_cnv_de_obj))) < 2) {
  stop(
    "FindAllMarkers needs at least two High-CNV malignant clusters after low-CNV removal. ",
    "Current number of High-CNV clusters: ",
    length(levels(Idents(high_cnv_de_obj)))
  )
}

# STEP 25. Run FindAllMarkers among remaining High-CNV malignant clusters

high_cnv_cluster_markers_all <- FindAllMarkers(
  high_cnv_de_obj,
  assay = "RNA_DE_HIGH_CNV",
  logfc.threshold = HIGH_CNV_CLUSTER_MARKER_LOGFC,
  min.pct = HIGH_CNV_CLUSTER_MARKER_MINPCT,
  return.thresh = 1,
  only.pos = TRUE
)

if (nrow(high_cnv_cluster_markers_all) > 0) {

  high_cnv_cluster_markers_all$Gene <- rownames(high_cnv_cluster_markers_all)
  high_cnv_cluster_markers_all$Marker_strategy <-
    "FindAllMarkers among High-CNV malignant clusters after Low-CNV removal"

  fc_col_high_cnv <- if ("avg_log2FC" %in% colnames(high_cnv_cluster_markers_all)) {
    "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(high_cnv_cluster_markers_all)) {
    "avg_logFC"
  } else {
    stop("No avg_log2FC or avg_logFC column was found in FindAllMarkers output.")
  }

  high_cnv_cluster_markers_all <- high_cnv_cluster_markers_all %>%
    dplyr::arrange(
      p_val_adj,
      cluster,
      dplyr::desc(.data[[fc_col_high_cnv]])
    )

} else {

  warning(
    "FindAllMarkers returned zero rows among High-CNV malignant clusters. ",
    "The thresholds may be too strict."
  )

  high_cnv_cluster_markers_all$Gene <- character(0)
  fc_col_high_cnv <- NA_character_
}

write.csv(
  high_cnv_cluster_markers_all,
  file.path(
    scrna_figdir,
    "High_CNV_Malignant_Cluster_FindAllMarkers_after_lowCNV_removed_logFC1_all.csv"
  ),
  row.names = FALSE
)

# STEP 26. Final filtering and top-gene selection

if (nrow(high_cnv_cluster_markers_all) > 0 && !is.na(fc_col_high_cnv)) {

  final_high_cnv_markers <- high_cnv_cluster_markers_all %>%
    dplyr::filter(
      .data[[fc_col_high_cnv]] >= HIGH_CNV_CLUSTER_MARKER_LOGFC,
      p_val_adj < HIGH_CNV_CLUSTER_MARKER_PADJ
    ) %>%
    dplyr::arrange(
      p_val_adj,
      cluster,
      dplyr::desc(.data[[fc_col_high_cnv]])
    )

} else {

  final_high_cnv_markers <- high_cnv_cluster_markers_all[0, , drop = FALSE]
}

if (nrow(final_high_cnv_markers) < 1 && nrow(high_cnv_cluster_markers_all) > 0) {

  warning(
    "No genes passed p_val_adj < HIGH_CNV_CLUSTER_MARKER_PADJ. ",
    "Saving top ranked High-CNV cluster markers instead."
  )

  final_high_cnv_markers <- high_cnv_cluster_markers_all %>%
    dplyr::arrange(
      p_val_adj,
      cluster,
      dplyr::desc(.data[[fc_col_high_cnv]])
    )
}

final_high_cnv_markers <- final_high_cnv_markers %>%
  dplyr::filter(!is.na(Gene), Gene != "") %>%
  dplyr::group_by(Gene) %>%
  dplyr::slice_min(order_by = p_val_adj, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    p_val_adj,
    dplyr::desc(.data[[fc_col_high_cnv]])
  )

final_high_cnv_gene_names <- data.frame(
  Gene = unique(final_high_cnv_markers$Gene),
  stringsAsFactors = FALSE
)

# STEP 27. Save final marker files

write.csv(
  high_cnv_cluster_markers_all,
  file.path(
    scrna_figdir,
    "High_CNV_Malignant_Cluster_FindAllMarkers_after_lowCNV_removed_logFC1_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  final_high_cnv_markers,
  file.path(scrna_figdir, "Final_High_CNV_Malignant_Genes.csv"),
  row.names = FALSE
)

write.csv(
  final_high_cnv_gene_names,
  file.path(scrna_figdir, "Final_High_CNV_Malignant_Gene_Names.csv"),
  row.names = FALSE
)

writeLines(
  final_high_cnv_gene_names$Gene,
  file.path(scrna_figdir, "Final_High_CNV_Malignant_Gene_Names.txt")
)

write.csv(
  final_high_cnv_markers,
  file.path(scrna_figdir, "High_CNV_Malignant_FindAllMarkers_after_lowCNV_removed_logFC1.csv"),
  row.names = FALSE
)

writeLines(
  final_high_cnv_gene_names$Gene,
  file.path(scrna_figdir, "High_CNV_Malignant_marker_genes_after_lowCNV_removed_logFC1.txt")
)

write.csv(
  final_high_cnv_markers,
  file.path(scrna_figdir, "High_CNV_Malignant_specific_markers_FINAL.csv"),
  row.names = FALSE
)

writeLines(
  final_high_cnv_gene_names$Gene,
  file.path(scrna_figdir, "High_CNV_Malignant_marker_genes_FINAL.txt")
)

# STEP 28. Save corrected objects

saveRDS(
  high_cnv_de_obj,
  file.path(scrna_rdsdir, "HIGH_CNV_MALIGNANT_object_DE_fixed.rds")
)

saveRDS(
  high_cnv_malignant_obj,
  file.path(scrna_rdsdir, "HIGH_CNV_MALIGNANT_object_after_lowCNV_removal.rds")
)

cat("\n###############################################################################\n")
cat("HIGH-CNV CLUSTER-SPECIFIC MARKER DETECTION FINISHED SUCCESSFULLY\n")
cat("###############################################################################\n")
cat("Marker strategy: FindAllMarkers among High-CNV malignant clusters after Low-CNV removal\n")
cat("Low-CNV clusters removed before marker detection:\n")
print(low_cnv_clusters)
cat("\nHigh-CNV cluster marker logFC threshold:", HIGH_CNV_CLUSTER_MARKER_LOGFC, "\n")
cat("High-CNV cluster marker min.pct:", HIGH_CNV_CLUSTER_MARKER_MINPCT, "\n")
cat("All High-CNV cluster marker rows:", nrow(high_cnv_cluster_markers_all), "\n")
cat("Final marker rows:", nrow(final_high_cnv_markers), "\n")
cat("Final unique genes:", nrow(final_high_cnv_gene_names), "\n")
cat("\nMain full marker CSV:\n")
cat(file.path(scrna_figdir, "High_CNV_Malignant_Cluster_FindAllMarkers_after_lowCNV_removed_logFC1_all.csv"), "\n")
cat("\nMain final detailed CSV:\n")
cat(file.path(scrna_figdir, "Final_High_CNV_Malignant_Genes.csv"), "\n")
cat("\nFinal gene names CSV:\n")
cat(file.path(scrna_figdir, "Final_High_CNV_Malignant_Gene_Names.csv"), "\n")
cat("\nDownstream gene list TXT:\n")
cat(file.path(scrna_figdir, "High_CNV_Malignant_marker_genes_after_lowCNV_removed_logFC1.txt"), "\n")
cat("###############################################################################\n")

# STEP 29. SAVE OBJECTS

saveRDS(
  seurat_obj,
  file.path(scrna_rdsdir, "FINAL_seurat_object_with_malignant_annotation.rds")
)

saveRDS(
  malignant_obj,
  file.path(scrna_rdsdir, "MALIGNANT_all_object_with_refined_CNV.rds")
)

saveRDS(
  high_cnv_malignant_obj,
  file.path(scrna_rdsdir, "HIGH_CNV_MALIGNANT_object_after_lowCNV_removal.rds")
)

# STEP 30. COMBINED MANUSCRIPT FIGURE

p_tsne_clusters_combined <- p_tsne_clusters + theme(legend.position = "none")
p_tsne_final_combined    <- p_tsne_final + theme(legend.position = "none")
p_dotplot_combined       <- p_dotplot
p_malig_sub_combined     <- p_malig_sub + theme(legend.position = "none")
p_initial_cnv_combined   <- p_initial_cnv_violin + theme(legend.position = "none")
p_infercnv_combined      <- p_infercnv_panel
p_refined_cnv_combined   <- p_refined_cnv_violin + theme(legend.position = "none")

fig_sc_combined <- (
  (p_tsne_clusters_combined | p_tsne_final_combined) /
    p_dotplot_combined /
    (p_malig_sub_combined | p_initial_cnv_combined) /
    p_infercnv_combined /
    p_refined_cnv_combined
) +
  patchwork::plot_layout(
    heights = c(1.0, 1.05, 1.0, 1.85, 1.0),
    widths = c(1, 1)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 18,
      family = FONT_FAMILY,
      colour = "black"
    ),
    plot.tag.position = c(0.01, 0.98)
  )

ggsave(
  filename = file.path(scrna_figdir, "Figure_SC_01_combined_single_cell_malignant_refined_inferCNV_AG.png"),
  plot = fig_sc_combined,
  width = 13.5,
  height = 22.5,
  dpi = FIG_DPI,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(scrna_figdir, "Figure_SC_01_combined_single_cell_malignant_refined_inferCNV.png"),
  plot = fig_sc_combined,
  width = 13.5,
  height = 22.5,
  dpi = FIG_DPI,
  bg = "white",
  limitsize = FALSE
)

# STEP 31. FINAL CONSOLE SUMMARY

cat("\n###############################################################################\n")
cat("LSCC SINGLE-CELL PIPELINE WITH REFINED inferCNV FINISHED\n")
cat("###############################################################################\n\n")

cat("Dataset: GSE206332\n")
cat("Samples: GSM6251294, GSM6251297, GSM6251300\n\n")

cat("QC filters used:\n")
cat("nFeature_RNA >= ", QC_MIN_FEATURES, "\n", sep = "")
cat("percent.mt < ", QC_MAX_MT, "\n\n", sep = "")

cat("Highly variable genes:", N_HVG, "\n")
cat("Harmony batch correction: applied by sample\n")
cat("Global FindAllMarkers: logfc.threshold = ", FINDALL_LOGFC, ", min.pct = ", FINDALL_MINPCT, "\n", sep = "")
cat("High-CNV malignant markers: FindAllMarkers among High-CNV clusters after Low-CNV removal\n")
cat("High-CNV cluster marker logFC threshold = ", HIGH_CNV_CLUSTER_MARKER_LOGFC, "\n", sep = "")
cat("\n")

cat("Final cell counts:\n")
print(table(seurat_obj$final_celltype))

cat("\nMalignant vs non-malignant:\n")
print(table(seurat_obj$malignant))

cat("\nWeak/low-CNV clusters selected from the initial CNV score and used as refined inferCNV reference and removed:\n")
print(low_cnv_clusters)

cat("\nInitial CNV weak/low-CNV reference selection table saved in:\n")
cat(file.path(scrna_figdir, "Initial_CNV_weak_low_CNV_reference_selection_table.csv"), "\n")

cat("\nRefined CNV cluster summary:\n")
print(refined_cnv_cluster_summary)

cat("\nFinal marker gene list:\n")
cat(file.path(scrna_figdir, "Final_High_CNV_Malignant_Gene_Names.txt"), "\n")

cat("\nMain combined figure with panels A-G:\n")
cat(file.path(scrna_figdir, "Figure_SC_01_combined_single_cell_malignant_refined_inferCNV_AG.png"), "\n")

cat("\nIndividual figure files:\n")
cat(file.path(scrna_figdir, "Figure_SC_01A_tSNE_global_clusters.png"), "\n")
cat(file.path(scrna_figdir, "Figure_SC_01B_tSNE_final_celltypes_with_malignant.png"), "\n")
cat(file.path(scrna_figdir, "Figure_SC_01C_DotPlot_celltype_markers.png"), "\n")
cat(file.path(scrna_figdir, "Figure_SC_01D_tSNE_malignant_subclusters.png"), "\n")
cat(file.path(scrna_figdir, "Figure_SC_01E_Initial_CNV_score_violin_for_reference_selection.png"), "\n")
cat(file.path(scrna_figdir, "Figure_SC_01F_inferCNV.png"), "\n")
cat(file.path(scrna_figdir, "Figure_SC_01G_Refined_CNV_score_violin_malignant_subclusters.png"), "\n")

cat("\n###############################################################################\n")

cat("\n###############################################################################\n")
cat("FILES FOR THE NEXT DOWNSTREAM ANALYSIS\n")
cat("###############################################################################\n")
cat("Use this CSV if the next analysis needs only the final gene names:\n")
cat(file.path(scrna_figdir, "Final_High_CNV_Malignant_Gene_Names.csv"), "\n")
cat("\nUse this CSV if the next analysis needs the full marker table with logFC and adjusted p-values:\n")
cat(file.path(scrna_figdir, "Final_High_CNV_Malignant_Genes.csv"), "\n")
cat("\nUse this TXT file if the downstream overlap code expects one gene symbol per line:\n")
cat(file.path(scrna_figdir, "High_CNV_Malignant_marker_genes_after_lowCNV_removed_logFC1.txt"), "\n")
cat("\nUse this CSV if you want the full FindAllMarkers output among the remaining High-CNV clusters:\n")
cat(file.path(scrna_figdir, "High_CNV_Malignant_Cluster_FindAllMarkers_after_lowCNV_removed_logFC1_all.csv"), "\n")
cat("\nDo NOT use any High-CNV vs Low-CNV DE table as the final marker list in this article-like workflow.\n")
cat("###############################################################################\n")

rm(list = ls())
gc()

options(stringsAsFactors = FALSE, scipen = 100, width = 150)

# STEP 1. REQUIRED PACKAGES

required_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "ggplot2",
  "patchwork",
  "grid",
  "png"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "These required packages are missing:\n",
    paste(missing_pkgs, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(png)
})

# STEP 2. PATHS

scrna_root_dir   <- "D:/LSCC/ScRNAseq_Results"
scrna_project_id <- "GSE206332"

scrna_base_dir <- file.path(scrna_root_dir, scrna_project_id)
scrna_out_dir  <- file.path(scrna_base_dir, "Results")

scrna_figdir <- file.path(scrna_out_dir, "figures")
scrna_rdsdir <- file.path(scrna_out_dir, "rds")
infercnv_out <- file.path(scrna_out_dir, "infercnv")

figure_out_dir <- file.path(
  scrna_figdir,
  "Standardized_Arial600_Reconstructed_NoC_NoG_CroppedInferCNV_LegendsABC"
)

dir.create(figure_out_dir, recursive = TRUE, showWarnings = FALSE)

# STEP 3. STANDARD MANUSCRIPT FIGURE SETTINGS

if (.Platform$OS.type == "windows") {
  try(
    grDevices::windowsFonts(
      Arial = grDevices::windowsFont("Arial")
    ),
    silent = TRUE
  )
}

FONT_FAMILY <- "Arial"
FIG_DPI <- 600
FIG_BACKGROUND <- "white"

FIG_SINGLE_W <- 3.50
FIG_SINGLE_H <- 4.20

FIG_DOUBLE_W <- 7.20
FIG_DOUBLE_H <- 5.40

FIGURE2_COMBINED_W <- 14.00
FIGURE2_COMBINED_H <- 6.40

BASE_TEXT_PT <- 9.5

AXIS_TEXT_PT <- 8.5
AXIS_TITLE_PT <- 10

LEGEND_TEXT_PT <- 8.5
LEGEND_TITLE_PT <- 9

PLOT_TITLE_PT <- 10
PANEL_TAG_PT <- 14

PANEL_BORDER_LWD <- 0.45
AXIS_LWD <- 0.40
GEOM_LWD <- 0.55
GRID_LWD <- 0.30

HEATMAP_FONT_PT <- 8
HEATMAP_ROW_FONT_PT <- 7
HEATMAP_COL_FONT_PT <- 7
HEATMAP_LEGEND_FONT_PT <- 8

COL_NORMAL <- "#7470B2"

COL_TUMOR <- "#D81B60"

COL_MYBL2_LOW <- COL_NORMAL
COL_MYBL2_HIGH <- COL_TUMOR

COL_UP <- COL_TUMOR
COL_DOWN <- COL_NORMAL
COL_NS <- "#B8B8B8"

CELLCHAT_PALETTE <- grDevices::colorRampPalette(
  c("#E8F6F4", "#B8E1DA", "#73C6BE", "#2D9D8C", "#006E63")
)(100)

MANUSCRIPT_MARGIN <- ggplot2::margin(
  t = 5.5,
  r = 6,
  b = 5.5,
  l = 6,
  unit = "pt"
)

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

      legend.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = LEGEND_TITLE_PT,
        color = "black"
      ),

      legend.text = ggplot2::element_text(
        family = FONT_FAMILY,
        size = LEGEND_TEXT_PT,
        color = "black"
      ),

      legend.position = legend_position,

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

      plot.margin = MANUSCRIPT_MARGIN
    )
}

save_plot_all_formats <- function(plot_obj,
                                  filename_stem,
                                  dir_path,
                                  width = FIG_DOUBLE_W,
                                  height = FIG_DOUBLE_H) {

  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  png_file <- file.path(dir_path, paste0(filename_stem, ".png"))
  tiff_file <- file.path(dir_path, paste0(filename_stem, ".tiff"))
  pdf_file <- file.path(dir_path, paste0(filename_stem, ".pdf"))

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
    bg = FIG_BACKGROUND,
    compression = "lzw",
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
    list(
      PNG = png_file,
      TIFF = tiff_file,
      PDF = pdf_file
    )
  )
}

# STEP 4. HELPER FUNCTIONS

first_existing <- function(paths) {
  paths <- unique(paths)
  paths <- paths[file.exists(paths)]

  if (length(paths) < 1) {
    return(NA_character_)
  }

  paths[1]
}

find_first_file <- function(directory,
                            pattern,
                            recursive = TRUE) {

  if (!dir.exists(directory)) {
    return(NA_character_)
  }

  hits <- list.files(
    directory,
    pattern = pattern,
    full.names = TRUE,
    recursive = recursive,
    ignore.case = TRUE
  )

  if (length(hits) < 1) {
    return(NA_character_)
  }

  hits[1]
}

join_layers_safe <- function(obj) {

  DefaultAssay(obj) <- "RNA"

  tryCatch(
    Seurat::JoinLayers(
      object = obj,
      assay = "RNA"
    ),
    error = function(e) obj
  )
}

extract_cluster_num <- function(x) {
  out <- suppressWarnings(
    as.integer(gsub("[^0-9]+", "", as.character(x)))
  )

  out[is.na(out)] <- seq_len(sum(is.na(out)))
  out
}

order_cluster_labels <- function(x) {

  x <- as.character(x)
  nums <- extract_cluster_num(x)

  unique(x[order(nums, x)])
}

make_discrete_palette <- function(levels_vec) {

  levels_vec <- as.character(levels_vec)
  n <- length(levels_vec)

  base_palette <- c(
    "#7470B2", "#D81B60", "#6AA6A5", "#C5A0D8",
    "#E0A458", "#5C7CBA", "#8AB17D", "#B56576",
    "#4F7C82", "#A78BBA", "#D68C45", "#7E8AA2",
    "#C66B8C", "#5A9B92", "#9381FF", "#D4A373",
    "#7393B3", "#B8B8B8", "#8F9E7A", "#B07AA1"
  )

  if (n <= length(base_palette)) {
    values <- base_palette[seq_len(n)]
  } else {
    values <- c(
      base_palette,
      grDevices::hcl.colors(
        n - length(base_palette),
        palette = "Dynamic"
      )
    )
  }

  setNames(values, levels_vec)
}

assert_metadata_column <- function(obj,
                                   column_name,
                                   object_name) {

  if (!(column_name %in% colnames(obj@meta.data))) {
    stop(
      "The metadata column '", column_name,
      "' is missing from ", object_name, "."
    )
  }
}

assert_reduction <- function(obj,
                             reduction_name,
                             object_name) {

  if (!(reduction_name %in% names(obj@reductions))) {
    stop(
      "The reduction '", reduction_name,
      "' is missing from ", object_name, "."
    )
  }
}

make_cluster_label_if_missing <- function(obj) {

  if (!("Cluster_label" %in% colnames(obj@meta.data))) {

    if (!("seurat_clusters" %in% colnames(obj@meta.data))) {
      stop(
        "Neither Cluster_label nor seurat_clusters was found in malignant_obj."
      )
    }

    obj$Cluster_label <- paste0(
      "Cluster ",
      as.character(obj$seurat_clusters)
    )
  }

  cluster_levels <- order_cluster_labels(obj$Cluster_label)

  obj$Cluster_label <- factor(
    as.character(obj$Cluster_label),
    levels = cluster_levels
  )

  obj
}

# STEP 5. LOAD EXISTING OBJECTS ONLY

seurat_rds_file <- first_existing(
  c(
    file.path(
      scrna_rdsdir,
      "FINAL_seurat_object_with_malignant_annotation.rds"
    ),
    find_first_file(
      scrna_rdsdir,
      pattern = "FINAL_seurat_object.*\\.rds$"
    )
  )
)

malignant_rds_file <- first_existing(
  c(
    file.path(
      scrna_rdsdir,
      "MALIGNANT_all_object_with_refined_CNV.rds"
    ),
    find_first_file(
      scrna_rdsdir,
      pattern = "MALIGNANT_all_object.*\\.rds$"
    )
  )
)

if (is.na(seurat_rds_file)) {
  stop(
    "Could not find FINAL_seurat_object_with_malignant_annotation.rds in:\n",
    scrna_rdsdir
  )
}

if (is.na(malignant_rds_file)) {
  stop(
    "Could not find MALIGNANT_all_object_with_refined_CNV.rds in:\n",
    scrna_rdsdir
  )
}

cat("\nLoading existing objects only; no analysis will be rerun.\n")
cat("Seurat object:\n", seurat_rds_file, "\n", sep = "")
cat("Malignant object:\n", malignant_rds_file, "\n", sep = "")

seurat_obj <- readRDS(seurat_rds_file)
malignant_obj <- readRDS(malignant_rds_file)

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- join_layers_safe(seurat_obj)

malignant_obj <- make_cluster_label_if_missing(malignant_obj)

assert_reduction(seurat_obj, "tsne", "seurat_obj")
assert_reduction(malignant_obj, "tsne", "malignant_obj")

assert_metadata_column(seurat_obj, "seurat_clusters", "seurat_obj")
assert_metadata_column(seurat_obj, "final_celltype", "seurat_obj")
assert_metadata_column(malignant_obj, "CNV_score", "malignant_obj")
assert_metadata_column(malignant_obj, "Refined_CNV_score", "malignant_obj")

# STEP 6. PANEL A — GLOBAL t-SNE CLUSTERS

global_cluster_levels <- order_cluster_labels(seurat_obj$seurat_clusters)

seurat_obj$seurat_clusters <- factor(
  as.character(seurat_obj$seurat_clusters),
  levels = global_cluster_levels
)

global_cluster_palette <- make_discrete_palette(global_cluster_levels)

p_sc_a <- DimPlot(
  seurat_obj,
  reduction = "tsne",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  cols = unname(global_cluster_palette)
) +
  scale_color_manual(
    values = global_cluster_palette,
    drop = FALSE
  ) +
  theme_manuscript(
    show_grid = FALSE,
    legend_position = "right"
  ) +
  labs(
    x = "t-SNE 1",
    y = "t-SNE 2",
    title = NULL,
    color = "Cluster"
  )

save_plot_all_formats(
  plot_obj = p_sc_a,
  filename_stem = "Figure_SC_01A_tSNE_global_clusters_standardized",
  dir_path = figure_out_dir,
  width = FIG_DOUBLE_W,
  height = 5.60
)

# STEP 7. PANEL B — FINAL CELL-TYPE ANNOTATION

celltype_preferred_order <- c(
  "Bcell",
  "Endothelial",
  "Epithelial",
  "Fibroblast",
  "Malignant",
  "Myeloid",
  "NKcell",
  "Tcell",
  "Unknown"
)

available_celltypes <- unique(as.character(seurat_obj$final_celltype))

celltype_levels <- c(
  intersect(celltype_preferred_order, available_celltypes),
  setdiff(sort(available_celltypes), celltype_preferred_order)
)

seurat_obj$final_celltype <- factor(
  as.character(seurat_obj$final_celltype),
  levels = celltype_levels
)

celltype_palette <- c(
  Bcell       = "#A78BBA",
  Endothelial = "#6AA6A5",
  Epithelial  = "#B7B3E0",
  Fibroblast  = "#8E9AAF",
  Malignant   = COL_TUMOR,
  Myeloid     = "#5C7CBA",
  NKcell      = "#7AA095",
  Tcell       = "#7470B2",
  Unknown     = COL_NS
)

celltype_palette <- celltype_palette[
  intersect(names(celltype_palette), celltype_levels)
]

additional_celltypes <- setdiff(celltype_levels, names(celltype_palette))

if (length(additional_celltypes) > 0) {
  celltype_palette <- c(
    celltype_palette,
    make_discrete_palette(additional_celltypes)
  )
}

p_sc_b <- DimPlot(
  seurat_obj,
  reduction = "tsne",
  group.by = "final_celltype",
  label = TRUE,
  repel = TRUE,
  cols = unname(celltype_palette)
) +
  scale_color_manual(
    values = celltype_palette,
    drop = FALSE
  ) +
  theme_manuscript(
    show_grid = FALSE,
    legend_position = "right"
  ) +
  labs(
    x = "t-SNE 1",
    y = "t-SNE 2",
    title = NULL,
    color = "Cell type"
  )

save_plot_all_formats(
  plot_obj = p_sc_b,
  filename_stem = "Figure_SC_01B_tSNE_final_celltypes_standardized",
  dir_path = figure_out_dir,
  width = FIG_DOUBLE_W,
  height = 5.60
)

# STEP 8. PANEL C — MALIGNANT t-SNE SUBCLUSTERS

malignant_cluster_levels <- levels(malignant_obj$Cluster_label)
malignant_cluster_palette <- make_discrete_palette(malignant_cluster_levels)

p_sc_d <- DimPlot(
  malignant_obj,
  reduction = "tsne",
  group.by = "Cluster_label",
  label = TRUE,
  repel = TRUE,
  cols = unname(malignant_cluster_palette)
) +
  scale_color_manual(
    values = malignant_cluster_palette,
    drop = FALSE
  ) +
  theme_manuscript(
    show_grid = FALSE,
    legend_position = "right"
  ) +
  labs(
    x = "t-SNE 1",
    y = "t-SNE 2",
    title = NULL,
    color = "Malignant subcluster"
  )

save_plot_all_formats(
  plot_obj = p_sc_d,
  filename_stem = "Figure_SC_01C_tSNE_malignant_subclusters_standardized",
  dir_path = figure_out_dir,
  width = FIG_DOUBLE_W,
  height = 5.60
)

# STEP 9. PANEL D — INITIAL CNV SCORE VIOLIN PLOT

p_sc_e <- ggplot(
  malignant_obj@meta.data,
  aes(
    x = Cluster_label,
    y = CNV_score,
    fill = Cluster_label
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    alpha = 0.85,
    linewidth = GEOM_LWD
  ) +
  geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = FIG_BACKGROUND,
    alpha = 0.90,
    linewidth = GEOM_LWD
  ) +
  scale_fill_manual(
    values = malignant_cluster_palette,
    drop = FALSE
  ) +
  theme_manuscript(
    show_grid = FALSE,
    legend_position = "right"
  ) +
  labs(
    x = "Malignant subclusters",
    y = "Initial CNV score",
    title = NULL,
    fill = "Cluster"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = AXIS_TEXT_PT
    )
  )

save_plot_all_formats(
  plot_obj = p_sc_e,
  filename_stem = "Figure_SC_01D_Initial_CNV_score_violin_standardized",
  dir_path = figure_out_dir,
  width = FIG_DOUBLE_W,
  height = 5.60
)

# STEP 10. PANEL E — EXISTING REFINED inferCNV IMAGE

infercnv_png_source <- first_existing(
  c(
    file.path(scrna_figdir, "Figure_SC_01F_inferCNV.png"),
    find_first_file(
      infercnv_out,
      pattern = "infercnv.*\\.png$"
    ),
    find_first_file(
      infercnv_out,
      pattern = "\\.png$"
    )
  )
)

if (is.na(infercnv_png_source)) {
  stop(
    "No existing inferCNV PNG image was found in:\n",
    scrna_figdir,
    "\nor\n",
    infercnv_out
  )
}

infercnv_img <- png::readPNG(infercnv_png_source)

INFERCNV_CROP_TOP_FRACTION <- 0.14

infercnv_nrow <- dim(infercnv_img)[1]
infercnv_ncol <- dim(infercnv_img)[2]

infercnv_row_start <- max(
  1L,
  floor(infercnv_nrow * INFERCNV_CROP_TOP_FRACTION) + 1L
)

if (length(dim(infercnv_img)) == 2L) {
  infercnv_img_cropped <- infercnv_img[
    infercnv_row_start:infercnv_nrow,
    1:infercnv_ncol,
    drop = FALSE
  ]
} else {
  infercnv_img_cropped <- infercnv_img[
    infercnv_row_start:infercnv_nrow,
    1:infercnv_ncol,
    ,
    drop = FALSE
  ]
}

p_sc_f <- ggplot() +
  annotation_custom(
    grob = grid::rasterGrob(
      infercnv_img_cropped,
      width = grid::unit(1, "npc"),
      height = grid::unit(1, "npc"),
      interpolate = TRUE
    ),
    xmin = 0,
    xmax = 1,
    ymin = 0,
    ymax = 1
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  coord_fixed(expand = FALSE) +
  theme_void() +
  theme(
    plot.background = element_rect(
      fill = FIG_BACKGROUND,
      color = NA
    ),
    panel.background = element_rect(
      fill = FIG_BACKGROUND,
      color = NA
    ),
    plot.margin = ggplot2::margin(
      t = 2,
      r = 2,
      b = 2,
      l = 2,
      unit = "pt"
    )
  )

save_plot_all_formats(
  plot_obj = p_sc_f,
  filename_stem = "Figure_SC_01E_existing_inferCNV_standardized",
  dir_path = figure_out_dir,
  width = 13.50,
  height = 10.50
)

# STEP 11. COMBINED SINGLE-CELL MANUSCRIPT FIGURE

p_sc_a_combined <- p_sc_a + theme(legend.position = "right")
p_sc_b_combined <- p_sc_b + theme(legend.position = "right")
p_sc_c_combined <- p_sc_d + theme(legend.position = "right")
p_sc_d_combined <- p_sc_e + theme(legend.position = "none")
p_sc_e_combined <- p_sc_f

fig_sc_combined <- (
  (p_sc_a_combined | p_sc_b_combined) /
    (p_sc_c_combined | p_sc_d_combined) /
    p_sc_e_combined
) +
  patchwork::plot_layout(
    heights = c(1.00, 1.00, 1.95),
    widths = c(1, 1)
  ) +
  patchwork::plot_annotation(
    tag_levels = "a"
  ) &
  theme(
    plot.tag = element_text(
      family = FONT_FAMILY,
      face = "bold",
      size = PANEL_TAG_PT,
      color = "black"
    ),
    plot.tag.position = c(0.01, 0.99)
  )

save_plot_all_formats(
  plot_obj = fig_sc_combined,
  filename_stem = "Figure_SC_01_combined_single_cell_noC_noG_croppedInferCNV_legendsABC_standardized",
  dir_path = figure_out_dir,
  width = 15.00,
  height = 16.40
)

# STEP 12. FINAL SUMMARY

cat("\n###############################################################################\n")
cat("SINGLE-CELL FIGURE RECONSTRUCTION FINISHED — PANELS C AND G REMOVED; INFERCNV CROPPED; LEGENDS RESTORED FOR PANELS A, B, AND C\n")
cat("###############################################################################\n")
cat("No QC, normalization, Harmony, clustering, inferCNV, or marker analysis was rerun.\n")
cat("\nExisting RDS objects used:\n")
cat("- ", seurat_rds_file, "\n", sep = "")
cat("- ", malignant_rds_file, "\n", sep = "")
cat("\nExisting inferCNV image used:\n")
cat("- ", infercnv_png_source, "\n", sep = "")
cat("\nNew standardized figures saved in:\n")
cat(figure_out_dir, "\n")
cat("\nMain combined figure:\n")
cat(
  file.path(
    figure_out_dir,
    "Figure_SC_01_combined_single_cell_noC_noG_croppedInferCNV_legendsABC_standardized.png"
  ),
  "\n"
)
cat("###############################################################################\n")

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 100,
  width = 150
)

# STEP 1. REQUIRED PACKAGES

required_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "ggplot2",
  "patchwork",
  "ggrepel",
  "grid",
  "png"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    "These packages are missing:\n",
    paste(missing_pkgs, collapse = ", "),
    "\n\nInstall them first using:\n",
    "install.packages(c('patchwork', 'ggrepel'))"
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(grid)
  library(png)
})

# STEP 2. PATHS

scrna_root_dir   <- "D:/LSCC/ScRNAseq_Results"
scrna_project_id <- "GSE206332"

scrna_base_dir <- file.path(
  scrna_root_dir,
  scrna_project_id
)

scrna_out_dir <- file.path(
  scrna_base_dir,
  "Results"
)

scrna_figdir <- file.path(
  scrna_out_dir,
  "figures"
)

scrna_rdsdir <- file.path(
  scrna_out_dir,
  "rds"
)

infercnv_out <- file.path(
  scrna_out_dir,
  "infercnv"
)

figure_out_dir <- file.path(
  scrna_figdir,
  "Manuscript_Ready_HighResolution_LargeInferCNV_NoOverlapLabels_FIXED"
)

dir.create(
  figure_out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# STEP 3. FIGURE SETTINGS

if (.Platform$OS.type == "windows") {

  try(
    grDevices::windowsFonts(
      Arial = grDevices::windowsFont("Arial")
    ),
    silent = TRUE
  )
}

FONT_FAMILY <- "Arial"

FIG_DPI <- 600
FIG_BACKGROUND <- "white"

BASE_TEXT_PT <- 10
AXIS_TEXT_PT <- 8.5
AXIS_TITLE_PT <- 10

LEGEND_TEXT_PT <- 8
LEGEND_TITLE_PT <- 9

PANEL_TAG_PT <- 15

POINT_SIZE_GLOBAL <- 0.32
POINT_SIZE_MALIGNANT <- 0.42

PANEL_BORDER_LWD <- 0.45
AXIS_LWD <- 0.40
GEOM_LWD <- 0.55

COL_MALIGNANT <- "#D81B60"
COL_UNKNOWN <- "#B8B8B8"

# STEP 4. HELPER FUNCTIONS

first_existing <- function(paths) {

  paths <- unique(paths)
  paths <- paths[file.exists(paths)]

  if (length(paths) < 1) {
    return(NA_character_)
  }

  return(paths[1])
}

find_first_file <- function(directory,
                            pattern,
                            recursive = TRUE) {

  if (!dir.exists(directory)) {
    return(NA_character_)
  }

  hits <- list.files(
    path = directory,
    pattern = pattern,
    full.names = TRUE,
    recursive = recursive,
    ignore.case = TRUE
  )

  if (length(hits) < 1) {
    return(NA_character_)
  }

  return(hits[1])
}

join_layers_safe <- function(obj) {

  DefaultAssay(obj) <- "RNA"

  obj_out <- tryCatch(
    Seurat::JoinLayers(
      object = obj,
      assay = "RNA"
    ),
    error = function(e) obj
  )

  return(obj_out)
}

extract_cluster_num <- function(x) {

  x <- as.character(x)

  out <- suppressWarnings(
    as.integer(
      gsub(
        pattern = "[^0-9]+",
        replacement = "",
        x = x
      )
    )
  )

  if (any(is.na(out))) {
    out[is.na(out)] <- seq_len(sum(is.na(out)))
  }

  return(out)
}

order_cluster_labels <- function(x) {

  x <- as.character(x)
  x <- unique(x)

  nums <- extract_cluster_num(x)

  return(x[order(nums, x)])
}

make_discrete_palette <- function(levels_vec) {

  levels_vec <- as.character(levels_vec)
  n_levels <- length(levels_vec)

  base_palette <- c(
    "#7470B2",
    "#D81B60",
    "#6AA6A5",
    "#C5A0D8",
    "#E0A458",
    "#5C7CBA",
    "#8AB17D",
    "#B56576",
    "#4F7C82",
    "#A78BBA",
    "#D68C45",
    "#7E8AA2",
    "#C66B8C",
    "#5A9B92",
    "#9381FF",
    "#D4A373",
    "#7393B3",
    "#B8B8B8",
    "#8F9E7A",
    "#B07AA1",
    "#AA6F73",
    "#6C9A8B",
    "#C38D9E",
    "#A8DADC",
    "#457B9D",
    "#E9C46A",
    "#F4A261",
    "#2A9D8F",
    "#9C89B8",
    "#F28482"
  )

  if (n_levels <= length(base_palette)) {

    palette_values <- base_palette[seq_len(n_levels)]

  } else {

    palette_values <- c(
      base_palette,
      grDevices::hcl.colors(
        n = n_levels - length(base_palette),
        palette = "Dynamic"
      )
    )
  }

  return(
    setNames(
      palette_values,
      levels_vec
    )
  )
}

theme_manuscript <- function(legend_position = "right") {

  ggplot2::theme_bw(
    base_size = BASE_TEXT_PT,
    base_family = FONT_FAMILY
  ) +
    ggplot2::theme(

      text = ggplot2::element_text(
        family = FONT_FAMILY,
        colour = "black"
      ),

      plot.tag = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PANEL_TAG_PT,
        colour = "black"
      ),

      plot.tag.position = c(0.01, 0.99),

      axis.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = AXIS_TITLE_PT,
        colour = "black"
      ),

      axis.text = ggplot2::element_text(
        family = FONT_FAMILY,
        size = AXIS_TEXT_PT,
        colour = "black"
      ),

      axis.line = ggplot2::element_line(
        colour = "black",
        linewidth = AXIS_LWD
      ),

      axis.ticks = ggplot2::element_line(
        colour = "black",
        linewidth = AXIS_LWD
      ),

      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = PANEL_BORDER_LWD
      ),

      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      legend.position = legend_position,

      legend.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = LEGEND_TITLE_PT,
        colour = "black"
      ),

      legend.text = ggplot2::element_text(
        family = FONT_FAMILY,
        size = LEGEND_TEXT_PT,
        colour = "black"
      ),

      legend.key.height = grid::unit(0.34, "cm"),
      legend.key.width = grid::unit(0.34, "cm"),

      plot.background = ggplot2::element_rect(
        fill = FIG_BACKGROUND,
        colour = NA
      ),

      panel.background = ggplot2::element_rect(
        fill = FIG_BACKGROUND,
        colour = NA
      ),

      plot.margin = ggplot2::margin(
        t = 6,
        r = 8,
        b = 6,
        l = 7,
        unit = "pt"
      )
    )
}

save_plot_all_formats <- function(plot_obj,
                                  filename_stem,
                                  dir_path,
                                  width,
                                  height) {

  dir.create(
    dir_path,
    recursive = TRUE,
    showWarnings = FALSE
  )

  png_file <- file.path(
    dir_path,
    paste0(filename_stem, ".png")
  )

  tiff_file <- file.path(
    dir_path,
    paste0(filename_stem, ".tiff")
  )

  pdf_file <- file.path(
    dir_path,
    paste0(filename_stem, ".pdf")
  )

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
    device = "tiff",
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

  return(
    invisible(
      list(
        PNG = png_file,
        TIFF = tiff_file,
        PDF = pdf_file
      )
    )
  )
}

get_tsne_dataframe <- function(obj,
                               group_column) {

  tsne_mat <- Seurat::Embeddings(
    object = obj,
    reduction = "tsne"
  )

  if (!(group_column %in% colnames(obj@meta.data))) {
    stop(
      "Metadata column not found: ",
      group_column
    )
  }

  plot_df <- data.frame(
    tSNE_1 = tsne_mat[, 1],
    tSNE_2 = tsne_mat[, 2],
    Group = as.character(
      obj@meta.data[[group_column]]
    ),
    stringsAsFactors = FALSE
  )

  return(plot_df)
}

get_group_centers <- function(plot_df) {

  center_df <- plot_df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(
      tSNE_1 = median(tSNE_1, na.rm = TRUE),
      tSNE_2 = median(tSNE_2, na.rm = TRUE),
      .groups = "drop"
    )

  return(center_df)
}

make_tsne_plot_with_nonoverlap_labels <- function(obj,
                                                  group_column,
                                                  palette,
                                                  legend_title,
                                                  point_size = 0.35,
                                                  label_size = 3.0,
                                                  label_box_padding = 0.35,
                                                  label_point_padding = 0.15,
                                                  legend_position = "right",
                                                  seed_value = 123) {

  plot_df <- get_tsne_dataframe(
    obj = obj,
    group_column = group_column
  )

  plot_df$Group <- factor(
    plot_df$Group,
    levels = names(palette)
  )

  center_df <- get_group_centers(plot_df)

  center_df$Group <- factor(
    center_df$Group,
    levels = names(palette)
  )

  set.seed(seed_value)

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = tSNE_1,
      y = tSNE_2,
      colour = Group
    )
  ) +
    ggplot2::geom_point(
      size = point_size,
      alpha = 0.85
    ) +
    ggrepel::geom_label_repel(
      data = center_df,
      ggplot2::aes(
        x = tSNE_1,
        y = tSNE_2,
        label = Group
      ),
      inherit.aes = FALSE,
      family = FONT_FAMILY,
      fontface = "plain",
      size = label_size,
      colour = "black",
      fill = scales::alpha("white", 0.84),
      label.size = 0.18,
      label.padding = grid::unit(0.14, "lines"),
      box.padding = grid::unit(label_box_padding, "lines"),
      point.padding = grid::unit(label_point_padding, "lines"),
      min.segment.length = 0,
      segment.color = "grey45",
      segment.size = 0.28,
      max.overlaps = Inf,
      force = 2.5,
      force_pull = 0.8,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = palette,
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(
      clip = "off"
    ) +
    theme_manuscript(
      legend_position = legend_position
    ) +
    ggplot2::labs(
      x = "t-SNE 1",
      y = "t-SNE 2",
      colour = legend_title
    )

  return(p)
}

# STEP 5. LOAD EXISTING OBJECTS ONLY

seurat_rds_file <- first_existing(
  c(
    file.path(
      scrna_rdsdir,
      "FINAL_seurat_object_with_malignant_annotation.rds"
    ),
    find_first_file(
      directory = scrna_rdsdir,
      pattern = "FINAL_seurat_object.*\\.rds$"
    )
  )
)

malignant_rds_file <- first_existing(
  c(
    file.path(
      scrna_rdsdir,
      "MALIGNANT_all_object_with_refined_CNV.rds"
    ),
    find_first_file(
      directory = scrna_rdsdir,
      pattern = "MALIGNANT_all_object.*\\.rds$"
    )
  )
)

if (is.na(seurat_rds_file)) {
  stop(
    "Could not find FINAL_seurat_object_with_malignant_annotation.rds in:\n",
    scrna_rdsdir
  )
}

if (is.na(malignant_rds_file)) {
  stop(
    "Could not find MALIGNANT_all_object_with_refined_CNV.rds in:\n",
    scrna_rdsdir
  )
}

cat("\nLoading existing RDS objects only.\n")
cat("No single-cell analysis will be rerun.\n")

seurat_obj <- readRDS(seurat_rds_file)
malignant_obj <- readRDS(malignant_rds_file)

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- join_layers_safe(seurat_obj)

if (!("Cluster_label" %in% colnames(malignant_obj@meta.data))) {

  if (!("seurat_clusters" %in% colnames(malignant_obj@meta.data))) {
    stop(
      "Neither Cluster_label nor seurat_clusters was found in malignant_obj."
    )
  }

  malignant_obj$Cluster_label <- paste0(
    "Cluster ",
    as.character(malignant_obj$seurat_clusters)
  )
}

if (!("tsne" %in% names(seurat_obj@reductions))) {
  stop(
    "t-SNE reduction is missing from the global Seurat object."
  )
}

if (!("tsne" %in% names(malignant_obj@reductions))) {
  stop(
    "t-SNE reduction is missing from the malignant Seurat object."
  )
}

if (!("final_celltype" %in% colnames(seurat_obj@meta.data))) {
  stop(
    "final_celltype is missing from the global Seurat object."
  )
}

if (!("CNV_score" %in% colnames(malignant_obj@meta.data))) {
  stop(
    "CNV_score is missing from the malignant Seurat object."
  )
}

# STEP 6. PANEL A — GLOBAL t-SNE CLUSTERS

global_cluster_levels <- order_cluster_labels(
  seurat_obj$seurat_clusters
)

seurat_obj$seurat_clusters <- factor(
  as.character(seurat_obj$seurat_clusters),
  levels = global_cluster_levels
)

global_cluster_palette <- make_discrete_palette(
  global_cluster_levels
)

p_sc_a <- make_tsne_plot_with_nonoverlap_labels(
  obj = seurat_obj,
  group_column = "seurat_clusters",
  palette = global_cluster_palette,
  legend_title = "Cluster",
  point_size = POINT_SIZE_GLOBAL,
  label_size = 2.55,
  label_box_padding = 0.20,
  label_point_padding = 0.08,
  legend_position = "right",
  seed_value = 100
)

save_plot_all_formats(
  plot_obj = p_sc_a,
  filename_stem = "Figure_SC_01A_tSNE_global_clusters_high_resolution",
  dir_path = figure_out_dir,
  width = 8.2,
  height = 6.3
)

# STEP 7. PANEL B — FINAL CELL-TYPE ANNOTATION

celltype_preferred_order <- c(
  "Bcell",
  "Endothelial",
  "Epithelial",
  "Fibroblast",
  "Malignant",
  "Myeloid",
  "NKcell",
  "Tcell",
  "Unknown"
)

available_celltypes <- unique(
  as.character(seurat_obj$final_celltype)
)

celltype_levels <- c(
  intersect(
    celltype_preferred_order,
    available_celltypes
  ),
  setdiff(
    sort(available_celltypes),
    celltype_preferred_order
  )
)

seurat_obj$final_celltype <- factor(
  as.character(seurat_obj$final_celltype),
  levels = celltype_levels
)

celltype_palette <- c(
  Bcell = "#A78BBA",
  Endothelial = "#6AA6A5",
  Epithelial = "#B7B3E0",
  Fibroblast = "#8E9AAF",
  Malignant = COL_MALIGNANT,
  Myeloid = "#5C7CBA",
  NKcell = "#7AA095",
  Tcell = "#7470B2",
  Unknown = COL_UNKNOWN
)

celltype_palette <- celltype_palette[
  intersect(
    names(celltype_palette),
    celltype_levels
  )
]

additional_celltypes <- setdiff(
  celltype_levels,
  names(celltype_palette)
)

if (length(additional_celltypes) > 0) {

  celltype_palette <- c(
    celltype_palette,
    make_discrete_palette(additional_celltypes)
  )
}

p_sc_b <- make_tsne_plot_with_nonoverlap_labels(
  obj = seurat_obj,
  group_column = "final_celltype",
  palette = celltype_palette,
  legend_title = "Cell type",
  point_size = POINT_SIZE_GLOBAL,
  label_size = 3.00,
  label_box_padding = 0.90,
  label_point_padding = 0.22,
  legend_position = "right",
  seed_value = 200
)

save_plot_all_formats(
  plot_obj = p_sc_b,
  filename_stem = "Figure_SC_01B_tSNE_final_celltypes_no_overlap_labels",
  dir_path = figure_out_dir,
  width = 8.2,
  height = 6.3
)

# STEP 8. PANEL C — MALIGNANT EPITHELIAL SUBCLUSTERS

malignant_cluster_levels <- order_cluster_labels(
  malignant_obj$Cluster_label
)

malignant_obj$Cluster_label <- factor(
  as.character(malignant_obj$Cluster_label),
  levels = malignant_cluster_levels
)

malignant_cluster_palette <- make_discrete_palette(
  malignant_cluster_levels
)

p_sc_c <- make_tsne_plot_with_nonoverlap_labels(
  obj = malignant_obj,
  group_column = "Cluster_label",
  palette = malignant_cluster_palette,
  legend_title = "Malignant subcluster",
  point_size = POINT_SIZE_MALIGNANT,
  label_size = 2.50,
  label_box_padding = 0.62,
  label_point_padding = 0.16,
  legend_position = "right",
  seed_value = 300
)

save_plot_all_formats(
  plot_obj = p_sc_c,
  filename_stem = "Figure_SC_01C_tSNE_malignant_subclusters_no_overlap_labels",
  dir_path = figure_out_dir,
  width = 8.2,
  height = 6.3
)

# STEP 9. PANEL D — INITIAL CNV SCORE VIOLIN PLOT

p_sc_d <- ggplot2::ggplot(
  malignant_obj@meta.data,
  ggplot2::aes(
    x = Cluster_label,
    y = CNV_score,
    fill = Cluster_label
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    alpha = 0.88,
    linewidth = GEOM_LWD
  ) +
  ggplot2::geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.92,
    linewidth = GEOM_LWD
  ) +
  ggplot2::scale_fill_manual(
    values = malignant_cluster_palette,
    drop = FALSE
  ) +
  theme_manuscript(
    legend_position = "none"
  ) +
  ggplot2::labs(
    x = "Malignant subclusters",
    y = "Initial CNV score",
    fill = NULL
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 55,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    plot.margin = ggplot2::margin(
      t = 6,
      r = 6,
      b = 22,
      l = 7,
      unit = "pt"
    )
  )

save_plot_all_formats(
  plot_obj = p_sc_d,
  filename_stem = "Figure_SC_01D_initial_CNV_score_violin_high_resolution",
  dir_path = figure_out_dir,
  width = 8.2,
  height = 6.3
)

# STEP 10. PANEL E — LARGE REFINED inferCNV IMAGE

infercnv_png_source <- first_existing(
  c(
    file.path(
      scrna_figdir,
      "Figure_SC_01F_inferCNV.png"
    ),
    find_first_file(
      directory = infercnv_out,
      pattern = "infercnv.*\\.png$"
    ),
    find_first_file(
      directory = infercnv_out,
      pattern = "\\.png$"
    )
  )
)

if (is.na(infercnv_png_source)) {
  stop(
    "No existing inferCNV PNG image was found in:\n",
    scrna_figdir,
    "\nor\n",
    infercnv_out
  )
}

infercnv_img <- png::readPNG(
  infercnv_png_source
)

INFERCNV_CROP_TOP_FRACTION <- 0.14

infercnv_nrow <- dim(infercnv_img)[1]
infercnv_ncol <- dim(infercnv_img)[2]

infercnv_row_start <- max(
  1L,
  floor(
    infercnv_nrow *
      INFERCNV_CROP_TOP_FRACTION
  ) + 1L
)

if (length(dim(infercnv_img)) == 2L) {

  infercnv_img_cropped <- infercnv_img[
    infercnv_row_start:infercnv_nrow,
    1:infercnv_ncol,
    drop = FALSE
  ]

} else {

  infercnv_img_cropped <- infercnv_img[
    infercnv_row_start:infercnv_nrow,
    1:infercnv_ncol,
    ,
    drop = FALSE
  ]
}

p_sc_e <- ggplot2::ggplot() +
  ggplot2::annotation_custom(
    grob = grid::rasterGrob(
      infercnv_img_cropped,
      width = grid::unit(1, "npc"),
      height = grid::unit(1, "npc"),
      interpolate = TRUE
    ),
    xmin = 0,
    xmax = 1,
    ymin = 0,
    ymax = 1
  ) +
  ggplot2::xlim(0, 1) +
  ggplot2::ylim(0, 1) +
  ggplot2::coord_fixed(
    expand = FALSE
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(
      fill = FIG_BACKGROUND,
      colour = NA
    ),
    panel.background = ggplot2::element_rect(
      fill = FIG_BACKGROUND,
      colour = NA
    ),
    plot.margin = ggplot2::margin(
      t = 2,
      r = 2,
      b = 2,
      l = 2,
      unit = "pt"
    )
  )

save_plot_all_formats(
  plot_obj = p_sc_e,
  filename_stem = "Figure_SC_01E_refined_inferCNV_large_high_resolution",
  dir_path = figure_out_dir,
  width = 14.5,
  height = 11.5
)

# STEP 11. COMBINED MANUSCRIPT FIGURE

p_sc_a_combined <- p_sc_a +
  ggplot2::theme(
    legend.position = "right",
    legend.text = ggplot2::element_text(
      size = 6.8
    ),
    legend.key.height = grid::unit(
      0.24,
      "cm"
    ),
    legend.title = ggplot2::element_text(
      size = 7.5
    )
  )

p_sc_b_combined <- p_sc_b +
  ggplot2::theme(
    legend.position = "right",
    legend.text = ggplot2::element_text(
      size = 7.5
    ),
    legend.title = ggplot2::element_text(
      size = 8
    )
  )

p_sc_c_combined <- p_sc_c +
  ggplot2::theme(
    legend.position = "right",
    legend.text = ggplot2::element_text(
      size = 6.8
    ),
    legend.key.height = grid::unit(
      0.24,
      "cm"
    ),
    legend.title = ggplot2::element_text(
      size = 7.5
    )
  )

p_sc_d_combined <- p_sc_d
p_sc_e_combined <- p_sc_e

fig_sc_combined <- (
  (p_sc_a_combined | p_sc_b_combined) /
    (p_sc_c_combined | p_sc_d_combined) /
    p_sc_e_combined
) +
  patchwork::plot_layout(
    widths = c(1, 1),
    heights = c(1.08, 1.08, 2.75)
  ) +
  patchwork::plot_annotation(
    tag_levels = "a"
  ) &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      family = FONT_FAMILY,
      face = "bold",
      size = PANEL_TAG_PT,
      colour = "black"
    ),
    plot.tag.position = c(0.01, 0.99)
  )

save_plot_all_formats(
  plot_obj = fig_sc_combined,
  filename_stem = "Figure_SC_01_combined_manuscript_ready_large_inferCNV_no_overlap_labels",
  dir_path = figure_out_dir,
  width = 16.5,
  height = 20.8
)

# STEP 12. FINAL SUMMARY

cat("\n###############################################################################\n")
cat("MANUSCRIPT-READY SINGLE-CELL FIGURE RECONSTRUCTION FINISHED\n")
cat("###############################################################################\n")

cat("\nNo QC, normalization, Harmony, clustering, inferCNV, or marker analysis was rerun.\n")

cat("\nUpdated features:\n")
cat("- 600 dpi PNG output\n")
cat("- 600 dpi TIFF output\n")
cat("- PDF output for manuscript editing\n")
cat("- Larger refined inferCNV panel\n")
cat("- Non-overlapping cell-type labels\n")
cat("- Non-overlapping malignant-subcluster labels\n")
cat("- Larger combined manuscript figure\n")

cat("\nOutput folder:\n")
cat(figure_out_dir, "\n")

cat("\nMain combined figure:\n")
cat(
  file.path(
    figure_out_dir,
    "Figure_SC_01_combined_manuscript_ready_large_inferCNV_no_overlap_labels.png"
  ),
  "\n"
)

cat("\nIndividual panels:\n")

cat(
  file.path(
    figure_out_dir,
    "Figure_SC_01A_tSNE_global_clusters_high_resolution.png"
  ),
  "\n"
)

cat(
  file.path(
    figure_out_dir,
    "Figure_SC_01B_tSNE_final_celltypes_no_overlap_labels.png"
  ),
  "\n"
)

cat(
  file.path(
    figure_out_dir,
    "Figure_SC_01C_tSNE_malignant_subclusters_no_overlap_labels.png"
  ),
  "\n"
)

cat(
  file.path(
    figure_out_dir,
    "Figure_SC_01D_initial_CNV_score_violin_high_resolution.png"
  ),
  "\n"
)

cat(
  file.path(
    figure_out_dir,
    "Figure_SC_01E_refined_inferCNV_large_high_resolution.png"
  ),
  "\n"
)

cat("\n###############################################################################\n")

# MODULE 07 — Bulk ssGSEA validation of High-CNV Clusters 0 and 8

rm(list = ls(all.names = TRUE))
gc()

set.seed(123)

options(
  stringsAsFactors = FALSE,
  scipen = 100,
  timeout = 7200,
  width = 160
)

# 1. VERSION CHECK, PATHS, AND ANALYTICAL SETTINGS

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

SELECTED_SUBCLUSTERS <- c(
  "Cluster 0",
  "Cluster 8"
)

MARKER_PADJ_CUTOFF <- 0.05
MARKER_LOGFC_CUTOFF <- 1.00
MARKER_PCT_CUTOFF <- 0.25

TOP_SIGNATURE_GENES <- 50L
MIN_SIGNATURE_GENES <- 10L

REMOVE_MITOCHONDRIAL_GENES <- TRUE
REMOVE_RIBOSOMAL_GENES <- TRUE

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

FINAL_FIGURE_STEM <- "Figure_08_Bulk_ssGSEA_Clusters0_8"

# 2. MANUSCRIPT FIGURE SETTINGS

FONT_FAMILY <- "Arial"
FIG_DPI <- 600L
FIG_BACKGROUND <- "white"

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

# 3. PACKAGE CHECK AND LOADING

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

# 4. GENERAL HELPER FUNCTIONS

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

# 5. SIGNATURE CONSTRUCTION FROM STRICT scRNA-seq MARKERS

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

# 6. STATISTICAL AND FIGURE HELPERS

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

# 7. RUN THE COMPLETE TWO-CLUSTER ssGSEA WORKFLOW

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
          "Figure palette:",
          "  Normal = teal; Tumor = burnt orange; Cluster 8 scatter = plum.",
          "  Panel b (external validation) = subtle blue-tinted background and blue border.",
          "  Full-width row headers: Discovery cohort (GSE127165 + GSE142083) and Independent validation cohort (GSE130605)."
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

# 8. EXECUTE THE COMPLETE ANALYSIS

run_two_cluster_ssgsea()

# MODULE 08 — Monocle 2 input preparation / Part 1

rm(list = ls(all.names = TRUE))
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 100,
  timeout = 7200,
  width = 160
)
set.seed(123)

# 1. VERSION CHECK, REAL INPUT PATHS, AND OUTPUT PATHS

if (as.character(getRversion()) != "4.4.3") {
  stop(
    "Run this script only in R 4.4.3. Current R version: ",
    as.character(getRversion())
  )
}

SC_PROJECT_DIR <- "D:/LSCC/ScRNAseq_Results/GSE206332"
SC_RESULTS_DIR <- file.path(SC_PROJECT_DIR, "Results")
SC_RDS_DIR <- file.path(SC_RESULTS_DIR, "rds")
SC_FIG_DIR <- file.path(SC_RESULTS_DIR, "figures")

PSEUDOTIME_MONOCLE_DIR <- "D:/LSCC/Results_LSCC/Pseudotime_Monocle"

INPUT_RDS <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_full_analysis_input.rds"
)

CELL_BALANCE_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_input_cell_balance.csv"
)

ORDERING_GENES_INPUT_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_input_ordering_genes.csv"
)

PREP_PROVENANCE_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_input_provenance.txt"
)

PREP_STATUS_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_input_preparation_status.txt"
)

PREP_SESSION_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_input_preparation_sessionInfo.txt"
)

dir.create(PSEUDOTIME_MONOCLE_DIR, recursive = TRUE, showWarnings = FALSE)

# 2. ANALYSIS SETTINGS — MATCH THE EXISTING PART 2 SCRIPT

TARGET_GENE <- "MYBL2"

MAX_CELLS_PER_SAMPLE_SUBCLUSTER <- 350L
MAX_ORDERING_GENES_PER_SUBCLUSTER <- 200L
MIN_ORDERING_GENES <- 100L

FIXED_CLUSTER_ORDER <- c(
  "Cluster 0", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4",
  "Cluster 5", "Cluster 6", "Cluster 8", "Cluster 10", "Cluster 11"
)

G2M_ROOT_GENE_SET <- c(
  "CDK1", "CCNB1", "CCNB2", "CDC20", "CDC25C", "CENPA", "CENPE",
  "CENPF", "CENPM", "CENPU", "CKS1B", "CKS2", "MKI67", "NUSAP1",
  "PLK1", "PRC1", "TOP2A", "TPX2", "UBE2C", "BIRC5"
)

# 3. REQUIRED PACKAGES

required_packages <- c("Seurat", "SeuratObject", "Matrix")

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
    "These packages must be installed for corrected Part 1:\n",
    paste(missing_packages, collapse = ", "),
    "\n\nUse the normal R 4.4.3 library that contains Seurat."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

# 4. HELPER FUNCTIONS

find_first_existing <- function(candidates, label) {
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0L) {
    stop(
      label, " was not found. Checked these paths:\n",
      paste(candidates, collapse = "\n")
    )
  }

  normalizePath(existing[1], winslash = "/", mustWork = TRUE)
}

cluster_number <- function(x) {
  x <- trimws(as.character(x))
  output <- suppressWarnings(
    as.integer(sub("^.*?([0-9]+)\\s*$", "\\1", x))
  )
  output[!grepl("[0-9]+\\s*$", x)] <- NA_integer_
  output
}

normalise_cluster <- function(x) {
  x <- trimws(as.character(x))
  n <- cluster_number(x)
  x[!is.na(n)] <- paste0("Cluster ", n[!is.na(n)])
  x
}

assert_raw_counts <- function(x) {
  if (!inherits(x, "Matrix")) {
    stop("The extracted count matrix is not a sparse Matrix object.")
  }

  if (length(x@x) == 0L) {
    stop("The extracted count matrix has no non-zero values.")
  }

  if (any(!is.finite(x@x)) || any(x@x < 0)) {
    stop("The count matrix contains invalid values.")
  }

  if (any(abs(x@x - round(x@x)) > 1e-8)) {
    stop(
      "The extracted matrix is not raw integer counts. ",
      "Monocle 2 requires raw counts."
    )
  }

  invisible(TRUE)
}

extract_joined_counts <- function(object, assay, layer = "counts") {
  if (!(assay %in% names(object@assays))) {
    stop("Assay not found: ", assay)
  }

  DefaultAssay(object) <- assay

  layer_names <- tryCatch(
    SeuratObject::Layers(object[[assay]]),
    error = function(e) character(0)
  )

  layer_names <- layer_names[
    grepl(paste0("^", layer, "(\\.|$)"), layer_names)
  ]

  if (length(layer_names) <= 1L) {
    single_layer <- if (length(layer_names) == 1L) {
      tryCatch(
        SeuratObject::LayerData(
          object,
          assay = assay,
          layer = layer_names[1],
          fast = FALSE
        ),
        error = function(e) NULL
      )
    } else {
      NULL
    }

    if (is.null(single_layer) || ncol(single_layer) == 0L) {
      single_layer <- tryCatch(
        Seurat::GetAssayData(
          object,
          assay = assay,
          slot = layer
        ),
        error = function(e) NULL
      )
    }

    if (!is.null(single_layer) && ncol(single_layer) > 0L) {
      if (
        length(intersect(rownames(single_layer), colnames(object))) >
          length(intersect(colnames(single_layer), colnames(object)))
      ) {
        single_layer <- Matrix::t(single_layer)
      }

      retained_cells <- intersect(colnames(object), colnames(single_layer))

      return(
        as(
          single_layer[, retained_cells, drop = FALSE],
          "dgCMatrix"
        )
      )
    }
  }

  matrices <- lapply(layer_names, function(layer_name) {
    x <- SeuratObject::LayerData(
      object,
      assay = assay,
      layer = layer_name,
      fast = FALSE
    )

    if (
      length(intersect(rownames(x), colnames(object))) >
        length(intersect(colnames(x), colnames(object)))
    ) {
      x <- Matrix::t(x)
    }

    retained_cells <- intersect(colnames(x), colnames(object))

    x[, retained_cells, drop = FALSE]
  })

  matrices <- matrices[vapply(matrices, ncol, integer(1)) > 0L]

  if (length(matrices) == 0L) {
    stop("No usable counts layer was found in assay ", assay)
  }

  common_genes <- Reduce(intersect, lapply(matrices, rownames))

  if (length(common_genes) < 2L) {
    stop("Too few shared genes across count layers.")
  }

  output <- do.call(
    cbind,
    lapply(matrices, function(x) x[common_genes, , drop = FALSE])
  )

  output <- output[, !duplicated(colnames(output)), drop = FALSE]

  retained_cells <- intersect(colnames(object), colnames(output))

  as(output[, retained_cells, drop = FALSE], "dgCMatrix")
}

g2m_score <- function(count_matrix) {
  genes <- intersect(G2M_ROOT_GENE_SET, rownames(count_matrix))

  if (length(genes) < 5L) {
    stop(
      "Fewer than five G2/M root-score genes are present in the raw counts. ",
      "Present genes: ", paste(genes, collapse = ", ")
    )
  }

  library_size <- Matrix::colSums(count_matrix)

  normalised <- log1p(
    sweep(
      as.matrix(count_matrix[genes, , drop = FALSE]),
      2,
      pmax(library_size, 1),
      "/"
    ) * 10000
  )

  colMeans(normalised)
}

get_marker_columns <- function(marker_table) {
  gene_col <- intersect(
    c("Gene", "gene", "Gene.symbol", "gene_name", "Symbol", "symbol"),
    colnames(marker_table)
  )[1]

  cluster_col <- intersect(
    c("cluster", "Cluster", "Cluster_label", "Malignant_subcluster"),
    colnames(marker_table)
  )[1]

  padj_col <- intersect(
    c("p_val_adj", "p_adj", "adj.P.Val", "FDR"),
    colnames(marker_table)
  )[1]

  fc_col <- intersect(
    c("avg_log2FC", "avg_logFC", "logFC"),
    colnames(marker_table)
  )[1]

  pct_col <- intersect(
    c("pct.1", "pct1", "Percent.1", "pct_in"),
    colnames(marker_table)
  )[1]

  required <- c(
    Gene = gene_col,
    cluster = cluster_col,
    p_val_adj = padj_col,
    avg_logFC = fc_col,
    pct.1 = pct_col
  )

  if (any(is.na(required) | !nzchar(required))) {
    stop(
      "The selected marker CSV lacks a required column.\n",
      "Available columns: ", paste(colnames(marker_table), collapse = ", "),
      "\nRequired concepts: Gene, cluster, p_val_adj, avg_log2FC/avg_logFC, pct.1."
    )
  }

  required
}

# 5. LOCATE THE REAL PRIOR-PIPELINE OUTPUTS

HIGHCNV_OBJECT_RDS <- find_first_existing(
  c(
    file.path(
      SC_PROJECT_DIR,
      "Results_Modular",
      "Step_20_Final_Exports_and_Objects",
      "rds",
      "HIGH_CNV_MALIGNANT_object_after_lowCNV_removal.rds"
    ),
    file.path(
      SC_RDS_DIR,
      "HIGH_CNV_MALIGNANT_object_after_lowCNV_removal.rds"
    )
  ),
  "High-CNV malignant metadata object"
)

HIGHCNV_DE_OBJECT_RDS <- find_first_existing(
  c(
    file.path(
      SC_PROJECT_DIR,
      "Results_Modular",
      "Step_20_Final_Exports_and_Objects",
      "rds",
      "HIGH_CNV_MALIGNANT_object_DE_fixed.rds"
    ),
    file.path(
      SC_RDS_DIR,
      "HIGH_CNV_MALIGNANT_object_DE_fixed.rds"
    )
  ),
  "High-CNV malignant DE object"
)

MARKER_CSV <- find_first_existing(
  c(
    file.path(
      SC_FIG_DIR,
      "High_CNV_Malignant_Cluster_FindAllMarkers_after_lowCNV_removed_logFC1_all.csv"
    ),
    file.path(
      SC_FIG_DIR,
      "Final_High_CNV_Malignant_Genes.csv"
    ),
    file.path(
      SC_FIG_DIR,
      "High_CNV_Malignant_FindAllMarkers_after_lowCNV_removed_logFC1.csv"
    ),
    file.path(
      SC_FIG_DIR,
      "High_CNV_Malignant_specific_markers_FINAL.csv"
    )
  ),
  "Precomputed High-CNV marker CSV"
)

# 6. LOAD AND VALIDATE THE REAL OBJECTS

cat(
  "\n============================================================\n",
  "CORRECTED PART 1: PREPARING MONOCLE INPUT\n",
  "============================================================\n",
  sep = ""
)

cat("High-CNV metadata object:\n", HIGHCNV_OBJECT_RDS, "\n\n", sep = "")
cat("High-CNV DE/count object:\n", HIGHCNV_DE_OBJECT_RDS, "\n\n", sep = "")
cat("Precomputed marker CSV:\n", MARKER_CSV, "\n\n", sep = "")

highcnv_obj <- readRDS(HIGHCNV_OBJECT_RDS)
highcnv_de_obj <- readRDS(HIGHCNV_DE_OBJECT_RDS)

if (!inherits(highcnv_obj, "Seurat")) {
  stop("The High-CNV metadata RDS is not a Seurat object.")
}

if (!inherits(highcnv_de_obj, "Seurat")) {
  stop("The High-CNV DE RDS is not a Seurat object.")
}

shared_cells <- colnames(highcnv_de_obj)[
  colnames(highcnv_de_obj) %in% colnames(highcnv_obj)
]

if (length(shared_cells) < 100L) {
  stop(
    "Fewer than 100 shared cells were found between the two saved High-CNV objects."
  )
}

analysis_obj <- subset(highcnv_de_obj, cells = shared_cells)

source_meta <- as.data.frame(
  highcnv_obj@meta.data[shared_cells, , drop = FALSE],
  stringsAsFactors = FALSE
)

required_metadata <- c("sample", "Cluster_label", "CNV_class")
missing_metadata <- setdiff(required_metadata, colnames(source_meta))

if (length(missing_metadata) > 0L) {
  stop(
    "The High-CNV metadata object lacks: ",
    paste(missing_metadata, collapse = ", ")
  )
}

analysis_obj$sample <- as.character(source_meta$sample)
analysis_obj$Cluster_label <- as.character(source_meta$Cluster_label)
analysis_obj$CNV_class <- as.character(source_meta$CNV_class)

if (any(
  is.na(analysis_obj$CNV_class) |
    analysis_obj$CNV_class != "High-CNV malignant"
)) {
  stop(
    "The selected source object includes cells that are not High-CNV malignant."
  )
}

analysis_obj$Malignant_subcluster <- factor(
  normalise_cluster(analysis_obj$Cluster_label),
  levels = FIXED_CLUSTER_ORDER
)

if (any(is.na(analysis_obj$Malignant_subcluster))) {
  bad <- unique(
    as.character(analysis_obj$Cluster_label[
      is.na(analysis_obj$Malignant_subcluster)
    ])
  )

  stop(
    "Unrecognised Cluster_label values were found:\n",
    paste(bad, collapse = ", ")
  )
}

observed_subclusters <- unique(as.character(analysis_obj$Malignant_subcluster))

if (!setequal(observed_subclusters, FIXED_CLUSTER_ORDER)) {
  stop(
    "The expected ten High-CNV subclusters are not all present.\n",
    "Observed: ", paste(sort(observed_subclusters), collapse = ", ")
  )
}

ASSAY_NAME <- if ("RNA_DE_HIGH_CNV" %in% names(analysis_obj@assays)) {
  "RNA_DE_HIGH_CNV"
} else if ("RNA" %in% names(analysis_obj@assays)) {
  "RNA"
} else {
  stop("Neither RNA_DE_HIGH_CNV nor RNA assay is available in the DE object.")
}

count_matrix <- extract_joined_counts(
  analysis_obj,
  assay = ASSAY_NAME,
  layer = "counts"
)

assert_raw_counts(count_matrix)

retained_cells <- intersect(colnames(count_matrix), colnames(analysis_obj))

if (length(retained_cells) < 100L) {
  stop("Fewer than 100 cells overlap between raw counts and metadata.")
}

count_matrix <- count_matrix[, retained_cells, drop = FALSE]

metadata <- as.data.frame(
  analysis_obj@meta.data[retained_cells, , drop = FALSE],
  stringsAsFactors = FALSE
)

rownames(metadata) <- retained_cells

if (!identical(colnames(count_matrix), rownames(metadata))) {
  stop("Raw count matrix and cell metadata could not be aligned.")
}

if (!(TARGET_GENE %in% rownames(count_matrix))) {
  stop(TARGET_GENE, " is absent from the raw count matrix.")
}

# 7. BALANCED CELL SAMPLING

metadata$sample <- as.character(metadata$sample)
metadata$Malignant_subcluster <- factor(
  as.character(metadata$Malignant_subcluster),
  levels = FIXED_CLUSTER_ORDER
)

metadata$sample_subcluster <- paste(
  metadata$sample,
  metadata$Malignant_subcluster,
  sep = " || "
)

before_balance <- as.data.frame(
  table(
    Sample = metadata$sample,
    Subcluster = metadata$Malignant_subcluster
  )
)
names(before_balance)[3] <- "Cells_before"

groups <- split(seq_len(nrow(metadata)), metadata$sample_subcluster)

keep_index <- sort(
  unlist(
    lapply(groups, function(index) {
      if (length(index) > MAX_CELLS_PER_SAMPLE_SUBCLUSTER) {
        sample(index, MAX_CELLS_PER_SAMPLE_SUBCLUSTER)
      } else {
        index
      }
    }),
    use.names = FALSE
  )
)

retained_cells <- rownames(metadata)[keep_index]

count_matrix <- count_matrix[, retained_cells, drop = FALSE]
metadata <- metadata[retained_cells, , drop = FALSE]

if (!identical(colnames(count_matrix), rownames(metadata))) {
  stop("Counts and metadata are not aligned after balanced sampling.")
}

if (!setequal(
  unique(as.character(metadata$Malignant_subcluster)),
  FIXED_CLUSTER_ORDER
)) {
  stop(
    "At least one expected High-CNV subcluster disappeared after balanced sampling."
  )
}

after_balance <- as.data.frame(
  table(
    Sample = metadata$sample,
    Subcluster = metadata$Malignant_subcluster
  )
)
names(after_balance)[3] <- "Cells_after"

balance_table <- merge(
  before_balance,
  after_balance,
  by = c("Sample", "Subcluster"),
  all = TRUE
)
balance_table[is.na(balance_table)] <- 0L

write.csv(balance_table, CELL_BALANCE_CSV, row.names = FALSE)

# 8. SELECT ORDERING GENES FROM THE REAL PRECOMPUTED MARKER TABLE

markers <- read.csv(
  MARKER_CSV,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(markers) == 0L) {
  stop("The selected marker CSV has no rows.")
}

marker_columns <- get_marker_columns(markers)

markers_clean <- data.frame(
  Gene = as.character(markers[[marker_columns["Gene"]]]),
  cluster = as.character(markers[[marker_columns["cluster"]]]),
  p_val_adj = suppressWarnings(
    as.numeric(markers[[marker_columns["p_val_adj"]]])
  ),
  avg_logFC = suppressWarnings(
    as.numeric(markers[[marker_columns["avg_logFC"]]])
  ),
  pct.1 = suppressWarnings(
    as.numeric(markers[[marker_columns["pct.1"]]])
  ),
  stringsAsFactors = FALSE
)

markers_clean$cluster_label <- normalise_cluster(markers_clean$cluster)

markers_clean <- markers_clean[
  !is.na(markers_clean$Gene) &
    nzchar(markers_clean$Gene) &
    is.finite(markers_clean$p_val_adj) &
    markers_clean$p_val_adj < 0.05 &
    is.finite(markers_clean$avg_logFC) &
    markers_clean$avg_logFC >= 1 &
    is.finite(markers_clean$pct.1) &
    markers_clean$pct.1 >= 0.25 &
    markers_clean$cluster_label %in% FIXED_CLUSTER_ORDER &
    markers_clean$Gene %in% rownames(count_matrix),
  ,
  drop = FALSE
]

if (nrow(markers_clean) == 0L) {
  stop(
    "No precomputed markers passed the required filters after matching to raw counts.\n",
    "Thresholds: adjusted P < 0.05, logFC >= 1, pct.1 >= 0.25."
  )
}

markers_clean$cluster_order <- match(
  markers_clean$cluster_label,
  FIXED_CLUSTER_ORDER
)

markers_clean <- markers_clean[
  order(
    markers_clean$cluster_order,
    markers_clean$p_val_adj,
    -markers_clean$avg_logFC,
    -markers_clean$pct.1,
    markers_clean$Gene
  ),
  ,
  drop = FALSE
]

ordering_by_cluster <- lapply(FIXED_CLUSTER_ORDER, function(cluster_label) {
  cluster_genes <- markers_clean$Gene[
    markers_clean$cluster_label == cluster_label
  ]

  unique(head(cluster_genes, MAX_ORDERING_GENES_PER_SUBCLUSTER))
})

ordering_genes <- unique(unlist(ordering_by_cluster, use.names = FALSE))

ordering_genes <- ordering_genes[
  !grepl("^(MT-|RPL|RPS)", ordering_genes)
]

if (length(ordering_genes) < MIN_ORDERING_GENES) {
  stop(
    "Too few validated ordering genes remained: ",
    length(ordering_genes),
    ". Required minimum: ", MIN_ORDERING_GENES
  )
}

write.csv(
  data.frame(Gene = ordering_genes),
  ORDERING_GENES_INPUT_CSV,
  row.names = FALSE
)

# 9. CALCULATE G2/M ROOT SCORE AND WRITE THE EXACT PART-2 INPUT RDS

metadata$G2M_root_score <- g2m_score(count_matrix)

if (any(!is.finite(metadata$G2M_root_score))) {
  stop("Non-finite G2/M root scores were produced.")
}

phenotype_data <- metadata[
  ,
  c("sample", "Malignant_subcluster", "G2M_root_score"),
  drop = FALSE
]

phenotype_data$sample <- factor(as.character(phenotype_data$sample))
phenotype_data$Malignant_subcluster <- factor(
  as.character(phenotype_data$Malignant_subcluster),
  levels = FIXED_CLUSTER_ORDER
)

rownames(phenotype_data) <- colnames(count_matrix)

if (!identical(rownames(phenotype_data), colnames(count_matrix))) {
  stop("Final Monocle phenotype metadata and count matrix are not aligned.")
}

saveRDS(
  list(
    schema_version = "LSCC_Monocle2_real_scRNA_input_v2",
    count_mat = as(count_matrix, "dgCMatrix"),
    pd_df = phenotype_data,
    ordering_genes = ordering_genes,
    target_gene = TARGET_GENE,
    root_column = "G2M_root_score",
    expected_subclusters = FIXED_CLUSTER_ORDER,
    max_cells_per_sample_subcluster = MAX_CELLS_PER_SAMPLE_SUBCLUSTER,
    root_rule = "Root = Monocle state with the lowest median G2M_root_score.",
    source_highcnv_object_rds = HIGHCNV_OBJECT_RDS,
    source_highcnv_de_object_rds = HIGHCNV_DE_OBJECT_RDS,
    source_marker_csv = MARKER_CSV,
    raw_count_assay = ASSAY_NAME
  ),
  INPUT_RDS
)

writeLines(
  c(
    "SUCCESS: Corrected Monocle input preparation completed in RGui 4.4.3.",
    "",
    "This script intentionally bypassed the nonexistent Functional_Enrichment_Objects.rds.",
    "",
    paste0("High-CNV metadata object used: ", HIGHCNV_OBJECT_RDS),
    paste0("High-CNV DE/count object used: ", HIGHCNV_DE_OBJECT_RDS),
    paste0("Precomputed marker CSV used: ", MARKER_CSV),
    paste0("Raw-count assay used: ", ASSAY_NAME),
    paste0("Prepared input RDS: ", INPUT_RDS),
    paste0("Cells retained after balancing: ", ncol(count_matrix)),
    paste0("Genes retained: ", nrow(count_matrix)),
    paste0("Ordering genes retained: ", length(ordering_genes)),
    paste0(
      "Balanced-sampling cap: ",
      MAX_CELLS_PER_SAMPLE_SUBCLUSTER,
      " cells per sample-subcluster."
    ),
    "",
    "NEXT STEP:",
    "Close RGui completely.",
    "Open a fresh RGui 4.4.3 session.",
    "Run your existing Part 2 script unchanged."
  ),
  PREP_PROVENANCE_TXT
)

writeLines(
  c(
    "SUCCESS: Corrected Part 1 Monocle input preparation finished.",
    paste0("Input RDS: ", INPUT_RDS),
    "Part 2 may now be run only in a fresh RGui session without Seurat loaded."
  ),
  PREP_STATUS_TXT
)

writeLines(capture.output(sessionInfo()), PREP_SESSION_TXT)

cat(
  "\n============================================================\n",
  "CORRECTED PART 1 FINISHED SUCCESSFULLY\n",
  "Prepared input:\n",
  INPUT_RDS,
  "\n\nNOW CLOSE RGui COMPLETELY.\n",
  "Then reopen RGui 4.4.3 and run the existing Part 2 script.\n",
  "============================================================\n",
  sep = ""
)

stopifnot(file.exists(INPUT_RDS))
print(file.info(INPUT_RDS)[, c("size", "mtime")])

# MODULE 09 — Monocle 2 DDRTree pseudotime / Part 2

rm(list = ls(all.names = TRUE))
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 100,
  timeout = 7200,
  width = 160
)
set.seed(123)

# 1. VERSION CHECK AND PATHS

if (as.character(getRversion()) != "4.4.3") {
  stop(
    "Run this script only in R 4.4.3. Current R version: ",
    as.character(getRversion())
  )
}

PSEUDOTIME_STEP1_DIR <- "D:/LSCC/Results_LSCC/Pseudotime_Step1"

PSEUDOTIME_MONOCLE_DIR <- "D:/LSCC/Results_LSCC/Pseudotime_Monocle"

MONOCLE_LIB <- "D:/LSCC/Monocle2_R44_library"

FUNCTIONAL_RDS <- file.path(
  PSEUDOTIME_STEP1_DIR,
  "Functional_Enrichment_Objects.rds"
)

INPUT_RDS <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_full_analysis_input.rds"
)

CELL_BALANCE_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_input_cell_balance.csv"
)
ORDERING_GENES_INPUT_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_input_ordering_genes.csv"
)
PREP_PROVENANCE_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_input_provenance.txt"
)
PREP_STATUS_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_input_preparation_status.txt"
)
PREP_SESSION_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_input_preparation_sessionInfo.txt"
)

STATUS_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_status.txt"
)
ERROR_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_error.txt"
)
RUN_LOG <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_run.log"
)
SESSION_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_sessionInfo.txt"
)
COMPATIBILITY_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_dplyr_compatibility_note.txt"
)
REDUCTION_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_reduction_status.txt"
)
CDS_RDS <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "HighCNV_Monocle2_CDS.rds"
)
TRAJECTORY_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_pseudotime_metadata.csv"
)
ORDERING_GENES_TRAJECTORY_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_ordering_genes_used.csv"
)
ROOT_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_root_state_median_G2M_scores.csv"
)
ROOT_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_root_rule.txt"
)
MYBL2_CELL_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "MYBL2_cell_level_pseudotime_expression.csv"
)
LOESS_CURVE_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "MYBL2_selected_clusters_LOESS_curve_data.csv"
)
BACKBONE_STATUS_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_DDRTree_backbone_status.txt"
)
BACKBONE_SEGMENTS_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_DDRTree_backbone_segments.csv"
)

FIGURE_3A_STEM <- "Figure_3A_HighCNV_DDRTree_pseudotime"
FIGURE_3B_STEM <- "Figure_3B_MYBL2_across_pseudotime"
FIGURE_8_STEM  <- "Figure_8_Combined_Pseudotime_MYBL2"

dir.create(
  PSEUDOTIME_MONOCLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. ANALYSIS SETTINGS

TARGET_GENE <- "MYBL2"

MAX_CELLS_PER_SAMPLE_SUBCLUSTER <- 350L

MAX_ORDERING_GENES_PER_SUBCLUSTER <- 200L
MIN_ORDERING_GENES <- 100L

FIXED_CLUSTER_ORDER <- c(
  "Cluster 0", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4",
  "Cluster 5", "Cluster 6", "Cluster 8", "Cluster 10", "Cluster 11"
)

IMPORTANT_CLUSTERS <- c(
  "Cluster 0",
  "Cluster 8"
)

G2M_ROOT_GENE_SET <- c(
  "CDK1", "CCNB1", "CCNB2", "CDC20", "CDC25C", "CENPA", "CENPE",
  "CENPF", "CENPM", "CENPU", "CKS1B", "CKS2", "MKI67", "NUSAP1",
  "PLK1", "PRC1", "TOP2A", "TPX2", "UBE2C", "BIRC5"
)

LOESS_SPAN <- c(
  "Cluster 0" = 0.90,
  "Cluster 8" = 1.00
)

LOESS_CI_LEVEL <- 0.95

LOESS_RIBBON_Q_LOWER <- 0.05
LOESS_RIBBON_Q_UPPER <- 0.95
LOESS_RIBBON_MIN_LOCAL_CELLS <- 20L
LOESS_RIBBON_WINDOW_FRACTION <- 0.12

# 3. STANDARD MANUSCRIPT FIGURE SETTINGS

FONT_FAMILY <- "Arial"
FIG_DPI <- 600L
FIG_BACKGROUND <- "white"

FIG_SINGLE_W <- 3.50
FIG_SINGLE_H <- 4.20
FIG_DOUBLE_W <- 7.20
FIG_DOUBLE_H <- 5.40

FIGURE3_COMBINED_W <- 14.00
FIGURE3_COMBINED_H <- 7.20

FIGURE8_COMBINED_W <- 14.00
FIGURE8_COMBINED_H <- 12.40

BASE_TEXT_PT <- 9.5
AXIS_TEXT_PT <- 8.5
AXIS_TITLE_PT <- 10
LEGEND_TEXT_PT <- 8.5
LEGEND_TITLE_PT <- 9
PLOT_TITLE_PT <- 10
PANEL_TAG_PT <- 14

PANEL_BORDER_LWD <- 0.45
AXIS_LWD <- 0.40
GEOM_LWD <- 0.55
GRID_LWD <- 0.30

HEATMAP_FONT_PT <- 8
HEATMAP_ROW_FONT_PT <- 7
HEATMAP_COL_FONT_PT <- 7
HEATMAP_LEGEND_FONT_PT <- 8

COL_NORMAL <- "#7470B2"
COL_TUMOR <- "#D81B60"
COL_MYBL2_LOW <- COL_NORMAL
COL_MYBL2_HIGH <- COL_TUMOR
COL_UP <- COL_TUMOR
COL_DOWN <- COL_NORMAL
COL_NS <- "#B8B8B8"

CELLCHAT_PALETTE <- grDevices::colorRampPalette(
  c("#E8F6F4", "#B8E1DA", "#73C6BE", "#2D9D8C", "#006E63")
)(100)

PSEUDOTIME_COLOURS <- c(
  "#132B43", "#2166AC", "#4EB3D3",
  "#A1DAB4", "#FEE08B", "#D73027"
)

MONOCLE_STATE_COLOURS <- c(
  "#0072B2", "#E69F00", "#009E73", "#CC79A7",
  "#56B4E9", "#D55E00", "#000000", "#999999"
)

IMPORTANT_CLUSTER_COLOURS <- c(
  "Cluster 0" = COL_TUMOR,
  "Cluster 8" = COL_NORMAL,
  "Cluster 11" = "#333333"
)

theme_manuscript <- function(
    show_grid = FALSE,
    legend_position = "bottom"
) {
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
      plot.tag.position = c(0.012, 0.988),
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
      legend.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = LEGEND_TITLE_PT,
        color = "black"
      ),
      legend.text = ggplot2::element_text(
        family = FONT_FAMILY,
        size = LEGEND_TEXT_PT,
        color = "black"
      ),
      legend.position = legend_position,
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
      plot.margin = ggplot2::margin(
        t = 5.5, r = 6, b = 5.5, l = 6, unit = "pt"
      )
    )
}

save_plot_all_formats <- function(
    plot_obj,
    filename_stem,
    dir_path,
    width = FIG_DOUBLE_W,
    height = FIG_DOUBLE_H
) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  png_file <- file.path(dir_path, paste0(filename_stem, ".png"))
  tiff_file <- file.path(dir_path, paste0(filename_stem, ".tiff"))
  pdf_file <- file.path(dir_path, paste0(filename_stem, ".pdf"))

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
    device = "tiff",
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

# 4. GENERAL HELPER FUNCTIONS

cluster_number <- function(x) {
  x <- as.character(x)
  output <- suppressWarnings(
    as.integer(sub("^.*?([0-9]+)\\s*$", "\\1", x))
  )
  output[!grepl("[0-9]+\\s*$", x)] <- NA_integer_
  output
}

normalise_cluster <- function(x) {
  x <- trimws(as.character(x))
  numbers <- cluster_number(x)
  x[!is.na(numbers)] <- paste0("Cluster ", numbers[!is.na(numbers)])
  x
}

assert_raw_counts <- function(x) {
  if (!inherits(x, "Matrix")) {
    stop("The count matrix is not a sparse Matrix object.")
  }
  if (any(!is.finite(x@x)) || any(x@x < 0)) {
    stop("The count matrix contains invalid values.")
  }
  if (any(abs(x@x - round(x@x)) > 1e-8)) {
    stop("The count matrix is not raw integer counts.")
  }
  invisible(TRUE)
}

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

# 6. MONOCLE / FIGURE HELPER FUNCTIONS

estimate_dispersions_blind_compat <- function(
    cds,
    min_cells_detected = 1L,
    remove_outliers = TRUE,
    verbose = FALSE
) {
  family_name <- cds@expressionFamily@vfamily

  if (!(family_name %in% c("negbinomial", "negbinomial.size"))) {
    stop(
      "Compatibility dispersion estimation requires a negative-binomial CellDataSet."
    )
  }

  dispersion_table <- monocle:::disp_calc_helper_NB(
    cds,
    cds@expressionFamily,
    min_cells_detected
  )

  dispersion_table <- as.data.frame(
    dispersion_table,
    stringsAsFactors = FALSE
  )

  required_columns <- c("gene_id", "mu", "disp")

  if (!all(required_columns %in% colnames(dispersion_table))) {
    stop("Unexpected Monocle dispersion-table structure.")
  }

  dispersion_table <- dispersion_table[
    is.finite(dispersion_table$mu) &
      is.finite(dispersion_table$disp) &
      dispersion_table$mu > 0 &
      dispersion_table$disp >= 0,
    ,
    drop = FALSE
  ]

  if (nrow(dispersion_table) < 20L) {
    stop("Too few genes remained for pooled dispersion fitting.")
  }

  rownames(dispersion_table) <- make.unique(
    as.character(dispersion_table$gene_id)
  )

  fit_result <- monocle:::parametricDispersionFit(
    dispersion_table,
    verbose = verbose
  )

  fit <- fit_result[[1]]
  coefficients <- fit_result[[2]]

  if (isTRUE(remove_outliers)) {
    cooks <- tryCatch(
      stats::cooks.distance(fit),
      error = function(e) numeric(0)
    )

    if (length(cooks) > 0L) {
      cutoff <- 4 / nrow(dispersion_table)

      outlier_rows <- union(
        names(cooks)[is.finite(cooks) & cooks > cutoff],
        setdiff(rownames(dispersion_table), names(cooks))
      )

      keep_rows <- !(rownames(dispersion_table) %in% outlier_rows)

      if (sum(keep_rows) >= 20L) {
        refit <- tryCatch(
          monocle:::parametricDispersionFit(
            dispersion_table[keep_rows, , drop = FALSE],
            verbose = verbose
          ),
          error = function(e) NULL
        )

        if (!is.null(refit)) {
          fit <- refit[[1]]
          coefficients <- refit[[2]]
        }
      }
    }
  }

  if (
    length(coefficients) != 2L ||
      any(!is.finite(coefficients)) ||
      coefficients[1] <= 0 ||
      coefficients[2] < 0
  ) {
    stop("Pooled dispersion fit produced invalid coefficients.")
  }

  names(coefficients) <- c("asymptDisp", "extraPois")

  dispersion_function <- function(q) {
    q <- pmax(as.numeric(q), .Machine$double.eps)
    coefficients["asymptDisp"] + coefficients["extraPois"] / q
  }

  attr(dispersion_function, "coefficients") <- coefficients

  cds@dispFitInfo[["blind"]] <- list(
    disp_table = dispersion_table,
    disp_func = dispersion_function
  )

  cds
}

make_trajectory_table <- function(cds_object) {
  metadata <- as.data.frame(
    Biobase::pData(cds_object),
    stringsAsFactors = FALSE
  )

  required_metadata <- c(
    "Pseudotime",
    "State",
    "Malignant_subcluster",
    "sample"
  )

  missing_metadata <- setdiff(required_metadata, colnames(metadata))

  if (length(missing_metadata) > 0L) {
    stop(
      "Trajectory metadata is missing: ",
      paste(missing_metadata, collapse = ", ")
    )
  }

  coordinates <- as.matrix(monocle::reducedDimS(cds_object))

  if (ncol(coordinates) == nrow(metadata)) {
    coordinates <- t(coordinates)
  }

  if (
    nrow(coordinates) != nrow(metadata) ||
      ncol(coordinates) < 2L
  ) {
    stop("Invalid DDRTree coordinate matrix.")
  }

  if (
    !is.null(rownames(coordinates)) &&
      all(rownames(metadata) %in% rownames(coordinates))
  ) {
    coordinates <- coordinates[rownames(metadata), , drop = FALSE]
  }

  state_values <- suppressWarnings(
    as.integer(as.character(metadata$State))
  )
  state_levels <- paste0(
    "State ",
    sort(unique(state_values[is.finite(state_values)]))
  )

  data.frame(
    Cell = rownames(metadata),
    DDRTree_1 = as.numeric(coordinates[, 1]),
    DDRTree_2 = as.numeric(coordinates[, 2]),
    Pseudotime = as.numeric(metadata$Pseudotime),
    Monocle_state = factor(
      paste0("State ", as.character(metadata$State)),
      levels = state_levels
    ),
    sample = as.character(metadata$sample),
    Malignant_subcluster = factor(
      normalise_cluster(metadata$Malignant_subcluster),
      levels = FIXED_CLUSTER_ORDER
    ),
    stringsAsFactors = FALSE
  )
}

make_loess_curve <- function(data,
                             cluster_name,
                             span_value,
                             ci_level = LOESS_CI_LEVEL) {
  current <- data[
    as.character(data$Malignant_subcluster) == cluster_name,
    ,
    drop = FALSE
  ]

  current <- current[
    is.finite(current$Pseudotime) &
      is.finite(current$log1p_size_factor_normalised_MYBL2),
    ,
    drop = FALSE
  ]
  current <- current[order(current$Pseudotime), , drop = FALSE]

  if (nrow(current) < 10L) {
    stop("Too few cells for MYBL2 LOESS fitting in ", cluster_name)
  }

  fit <- stats::loess(
    log1p_size_factor_normalised_MYBL2 ~ Pseudotime,
    data = current,
    span = span_value,
    degree = 1,
    family = "gaussian",
    control = stats::loess.control(surface = "direct")
  )

  pseudotime_grid <- seq(
    min(current$Pseudotime, na.rm = TRUE),
    max(current$Pseudotime, na.rm = TRUE),
    length.out = 150L
  )

  prediction <- tryCatch(
    stats::predict(
      fit,
      newdata = data.frame(Pseudotime = pseudotime_grid),
      se = TRUE
    ),
    error = function(e) NULL
  )

  if (
    is.list(prediction) &&
      all(c("fit", "se.fit") %in% names(prediction))
  ) {
    fitted_values <- as.numeric(prediction$fit)
    standard_errors <- as.numeric(prediction$se.fit)
  } else {
    fitted_values <- as.numeric(
      stats::predict(
        fit,
        newdata = data.frame(Pseudotime = pseudotime_grid)
      )
    )
    standard_errors <- rep(NA_real_, length(fitted_values))
  }

  critical_value <- stats::qnorm(
    1 - (1 - ci_level) / 2
  )

  lower_ci <- pmax(
    0,
    fitted_values - critical_value * standard_errors
  )
  upper_ci <- fitted_values + critical_value * standard_errors

  pt_quantiles <- stats::quantile(
    current$Pseudotime,
    probs = c(LOESS_RIBBON_Q_LOWER, LOESS_RIBBON_Q_UPPER),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )

  pt_range <- diff(range(current$Pseudotime, na.rm = TRUE))
  local_half_window <- max(
    pt_range * LOESS_RIBBON_WINDOW_FRACTION / 2,
    0.10
  )

  local_cell_count <- vapply(
    pseudotime_grid,
    function(current_pt) {
      sum(abs(current$Pseudotime - current_pt) <= local_half_window)
    },
    integer(1)
  )

  required_local_cells <- max(
    as.integer(LOESS_RIBBON_MIN_LOCAL_CELLS),
    as.integer(ceiling(0.05 * nrow(current)))
  )

  ribbon_supported <- (
    pseudotime_grid >= pt_quantiles[1] &
      pseudotime_grid <= pt_quantiles[2] &
      local_cell_count >= required_local_cells &
      is.finite(lower_ci) &
      is.finite(upper_ci)
  )

  ribbon_segment <- rep(NA_integer_, length(pseudotime_grid))

  if (any(ribbon_supported)) {
    transition <- c(
      TRUE,
      ribbon_supported[-1] != ribbon_supported[-length(ribbon_supported)]
    )
    segment_id <- cumsum(transition)
    ribbon_segment[ribbon_supported] <- segment_id[ribbon_supported]
  }

  data.frame(
    Malignant_subcluster = cluster_name,
    Pseudotime = pseudotime_grid,
    LOESS_smoothed_log1p_expression = fitted_values,
    LOESS_standard_error = standard_errors,
    LOESS_lower_CI = lower_ci,
    LOESS_upper_CI = upper_ci,
    LOESS_CI_level = ci_level,
    Local_cell_count = local_cell_count,
    Required_local_cells = required_local_cells,
    Ribbon_supported = ribbon_supported,
    Ribbon_segment = ribbon_segment,
    Ribbon_pseudotime_lower = pt_quantiles[1],
    Ribbon_pseudotime_upper = pt_quantiles[2],
    stringsAsFactors = FALSE
  )
}

extract_ddrtree_backbone_segments <- function(cds_object, trajectory_df) {
  empty_segments <- data.frame(
    x = numeric(0),
    y = numeric(0),
    xend = numeric(0),
    yend = numeric(0),
    Source = character(0),
    stringsAsFactors = FALSE
  )

  graph_to_segments <- function(graph, coordinates, source_name) {
    if (is.null(graph) || igraph::ecount(graph) < 1L) {
      return(empty_segments)
    }

    edges <- tryCatch(
      igraph::as_data_frame(graph, what = "edges"),
      error = function(e) NULL
    )

    if (
      is.null(edges) ||
        !all(c("from", "to") %in% colnames(edges))
    ) {
      return(empty_segments)
    }

    coordinates <- as.data.frame(
      coordinates,
      stringsAsFactors = FALSE
    )

    if (!all(c("Node", "x", "y") %in% colnames(coordinates))) {
      return(empty_segments)
    }

    coordinates$Node <- as.character(coordinates$Node)
    coordinates <- coordinates[
      !duplicated(coordinates$Node) &
        is.finite(coordinates$x) &
        is.finite(coordinates$y),
      ,
      drop = FALSE
    ]

    from_id <- as.character(edges$from)
    to_id <- as.character(edges$to)

    from_index <- match(from_id, coordinates$Node)
    to_index <- match(to_id, coordinates$Node)

    if (
      all(is.na(from_index)) &&
        all(is.na(to_index)) &&
        igraph::vcount(graph) == nrow(coordinates)
    ) {
      vertex_names <- igraph::V(graph)$name

      if (is.null(vertex_names) || length(vertex_names) != nrow(coordinates)) {
        vertex_names <- as.character(seq_len(nrow(coordinates)))
      }

      coordinates$Node <- as.character(vertex_names)
      from_index <- match(from_id, coordinates$Node)
      to_index <- match(to_id, coordinates$Node)
    }

    segments <- data.frame(
      x = coordinates$x[from_index],
      y = coordinates$y[from_index],
      xend = coordinates$x[to_index],
      yend = coordinates$y[to_index],
      Source = source_name,
      stringsAsFactors = FALSE
    )

    segments <- segments[
      is.finite(segments$x) &
        is.finite(segments$y) &
        is.finite(segments$xend) &
        is.finite(segments$yend),
      ,
      drop = FALSE
    ]

    segments
  }

  cell_coordinates <- data.frame(
    Node = as.character(trajectory_df$Cell),
    x = trajectory_df$DDRTree_1,
    y = trajectory_df$DDRTree_2,
    stringsAsFactors = FALSE
  )

  aux <- cds_object@auxOrderingData
  aux_names <- names(aux)

  ddr_keys <- unique(
    c(
      "DDRTree",
      aux_names[grepl("ddrtree", aux_names, ignore.case = TRUE)]
    )
  )
  ddr_keys <- ddr_keys[ddr_keys %in% aux_names]

  for (key in ddr_keys) {
    projected_graph <- tryCatch(
      aux[[key]]$pr_graph_cell_proj_tree,
      error = function(e) NULL
    )

    projected_segments <- graph_to_segments(
      projected_graph,
      cell_coordinates,
      source_name = paste0("projected_cell_graph:", key)
    )

    if (nrow(projected_segments) > 0L) {
      return(
        list(
          segments = projected_segments,
          status = paste0(
            "DDRTree backbone drawn from ",
            nrow(projected_segments),
            " projected-cell graph edges (",
            key,
            ")."
          )
        )
      )
    }
  }

  principal_coordinates <- tryCatch(
    as.matrix(monocle::reducedDimK(cds_object)),
    error = function(e) NULL
  )

  mst_graph <- tryCatch(
    monocle::minSpanningTree(cds_object),
    error = function(e) NULL
  )

  if (
    !is.null(principal_coordinates) &&
      length(principal_coordinates) > 0L &&
      !is.null(mst_graph)
  ) {
    if (
      nrow(principal_coordinates) == 2L &&
        ncol(principal_coordinates) != 2L
    ) {
      principal_coordinates <- t(principal_coordinates)
    }

    if (
      ncol(principal_coordinates) >= 2L &&
        nrow(principal_coordinates) >= 2L
    ) {
      node_names <- rownames(principal_coordinates)

      if (
        is.null(node_names) ||
          length(node_names) != nrow(principal_coordinates)
      ) {
        node_names <- as.character(seq_len(nrow(principal_coordinates)))
      }

      principal_df <- data.frame(
        Node = as.character(node_names),
        x = as.numeric(principal_coordinates[, 1]),
        y = as.numeric(principal_coordinates[, 2]),
        stringsAsFactors = FALSE
      )

      principal_segments <- graph_to_segments(
        mst_graph,
        principal_df,
        source_name = "principal_nodes:minSpanningTree"
      )

      if (nrow(principal_segments) > 0L) {
        return(
          list(
            segments = principal_segments,
            status = paste0(
              "DDRTree backbone drawn from ",
              nrow(principal_segments),
              " principal-node MST edges."
            )
          )
        )
      }
    }
  }

  cell_mst_segments <- graph_to_segments(
    mst_graph,
    cell_coordinates,
    source_name = "cell_coordinates:minSpanningTree"
  )

  if (nrow(cell_mst_segments) > 0L) {
    return(
      list(
        segments = cell_mst_segments,
        status = paste0(
          "DDRTree backbone drawn from ",
          nrow(cell_mst_segments),
          " cell-level MST edges."
        )
      )
    )
  }

  list(
    segments = empty_segments,
    status = paste0(
      "DDRTree backbone could not be drawn because no stored graph vertices ",
      "could be aligned with their DDRTree coordinates."
    )
  )
}

# 8. COMPLETE IGRAPH 2.3 COMPATIBILITY PATCH FOR MONOCLE 2

install_igraph_legacy_monocle_compatibility_patch <- function() {

  active_version <- as.character(utils::packageVersion("igraph"))

  if (utils::compareVersion(active_version, "1.3.0") < 0) {
    return(invisible(FALSE))
  }

  monocle_namespace <- asNamespace("monocle")

  replace_existing_namespace_binding <- function(
      namespace_environment,
      binding_name,
      value
  ) {
    if (!exists(
      binding_name,
      envir = namespace_environment,
      inherits = FALSE
    )) {
      stop(
        "Internal patch error: Monocle binding does not exist: ",
        binding_name
      )
    }

    was_locked <- bindingIsLocked(binding_name, namespace_environment)

    if (was_locked) {
      unlockBinding(binding_name, namespace_environment)
    }

    on.exit({
      if (was_locked && !bindingIsLocked(binding_name, namespace_environment)) {
        lockBinding(binding_name, namespace_environment)
      }
    }, add = TRUE)

    assign(binding_name, value, envir = namespace_environment)
    invisible(TRUE)
  }

  namespace_objects <- ls(monocle_namespace, all.names = TRUE)
  patched_functions <- character(0)

  for (object_name in namespace_objects) {
    object_value <- tryCatch(
      get(
        object_name,
        envir = monocle_namespace,
        inherits = FALSE
      ),
      error = function(e) NULL
    )

    if (!is.function(object_value)) {
      next
    }

    body_text <- paste(
      deparse(body(object_value), control = "all"),
      collapse = "\n"
    )

    has_legacy_nei <- grepl(
      "(?<![[:alnum:]_.:])nei\\s*\\(",
      body_text,
      perl = TRUE
    )

    has_legacy_neimode <- grepl(
      "\\bneimode\\s*=",
      body_text,
      perl = TRUE
    )

    has_bare_dfs <- grepl(
      "(?<![[:alnum:]_.:])dfs\\s*\\(",
      body_text,
      perl = TRUE
    )

    if (!has_legacy_nei && !has_legacy_neimode && !has_bare_dfs) {
      next
    }

    updated_body_text <- body_text

    updated_body_text <- gsub(
      "(?<![[:alnum:]_.:])nei\\s*\\(",
      "igraph:::.nei(",
      updated_body_text,
      perl = TRUE
    )

    updated_body_text <- gsub(
      "\\bneimode\\s*=",
      "mode =",
      updated_body_text,
      perl = TRUE
    )

    updated_body_text <- gsub(
      "(?<![[:alnum:]_.:])dfs\\s*\\(",
      "igraph::dfs(",
      updated_body_text,
      perl = TRUE
    )

    updated_body <- tryCatch(
      parse(text = updated_body_text)[[1]],
      error = function(e) {
        stop(
          "Could not patch legacy igraph calls in Monocle function '",
          object_name, "': ", conditionMessage(e)
        )
      }
    )

    body(object_value) <- updated_body

    replace_existing_namespace_binding(
      monocle_namespace,
      object_name,
      object_value
    )

    patched_functions <- c(patched_functions, object_name)
  }

  if (length(patched_functions) == 0L) {
    message(
      "No legacy igraph calls were found inside the active Monocle namespace; ",
      "no compatibility rewrite was required for igraph ", active_version, "."
    )
  } else {
    message(
      "Applied safe Monocle 2 / igraph ", active_version,
      " compatibility rewrite (existing Monocle functions only): ",
      paste(sort(unique(patched_functions)), collapse = ", "),
      "."
    )
  }

  invisible(length(patched_functions) > 0L)
}

# 9. LOAD MONOCLE PACKAGES IN A FRESH RGui SESSION

load_monocle_environment <- function() {
  if (!dir.exists(MONOCLE_LIB)) {
    stop("Monocle library folder is missing:\n", MONOCLE_LIB)
  }

  .libPaths(unique(c(MONOCLE_LIB, .libPaths())))

  required_packages <- c(
    "monocle", "Biobase", "BiocGenerics", "VGAM",
    "igraph", "DDRTree", "ggplot2", "patchwork", "Matrix"
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
      "Missing package(s) in the Monocle R library:\n",
      paste(missing_packages, collapse = ", ")
    )
  }

  suppressPackageStartupMessages({
    library(monocle)
    library(ggplot2)
    library(patchwork)
  })

  if (as.character(utils::packageVersion("monocle")) != "2.34.0") {
    stop(
      "Required monocle version: 2.34.0. Active version: ",
      utils::packageVersion("monocle")
    )
  }

  igraph_active_version <- as.character(utils::packageVersion("igraph"))

  if (utils::compareVersion(igraph_active_version, "2.0.0") < 0) {
    stop(
      "igraph version 2.0.0 or newer is required. Active version: ",
      igraph_active_version
    )
  }

  message("Using igraph version ", igraph_active_version, ".")

  install_igraph_legacy_monocle_compatibility_patch()

  invisible(TRUE)
}

# 9. FIGURE CREATION FROM A COMPLETED MONOCLE CellDataSet

create_all_figures <- function(cds, target_gene, root_state) {
  trajectory_df <- make_trajectory_table(cds)

  trajectory_df <- trajectory_df[
    is.finite(trajectory_df$DDRTree_1) &
      is.finite(trajectory_df$DDRTree_2) &
      is.finite(trajectory_df$Pseudotime) &
      !is.na(trajectory_df$Malignant_subcluster),
    ,
    drop = FALSE
  ]

  if (nrow(trajectory_df) < 100L) {
    stop("Too few valid cells remained for trajectory plotting.")
  }

  write.csv(trajectory_df, TRAJECTORY_CSV, row.names = FALSE)

  trajectory_x_limits <- range(trajectory_df$DDRTree_1, na.rm = TRUE)
  trajectory_y_limits <- range(trajectory_df$DDRTree_2, na.rm = TRUE)
  pseudotime_limits <- range(trajectory_df$Pseudotime, na.rm = TRUE)

  backbone_result <- extract_ddrtree_backbone_segments(cds, trajectory_df)
  backbone_segments <- backbone_result$segments

  write.csv(
    backbone_segments,
    BACKBONE_SEGMENTS_CSV,
    row.names = FALSE
  )

  writeLines(backbone_result$status, BACKBONE_STATUS_TXT)

  backbone_layer <- if (nrow(backbone_segments) > 0L) {
    ggplot2::geom_segment(
      data = backbone_segments,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      colour = "grey20",
      linewidth = 0.48,
      alpha = 0.80,
      lineend = "round"
    )
  } else {
    NULL
  }

  root_state_label <- paste0("State ", as.integer(root_state))
  root_state_cells <- trajectory_df[
    as.character(trajectory_df$Monocle_state) == root_state_label,
    ,
    drop = FALSE
  ]

  if (nrow(root_state_cells) < 1L) {
    stop(
      "The selected root state ",
      root_state_label,
      " was not found in the plotted trajectory."
    )
  }

  root_point <- root_state_cells[
    which.min(root_state_cells$Pseudotime),
    ,
    drop = FALSE
  ]

  x_span <- diff(trajectory_x_limits)
  y_span <- diff(trajectory_y_limits)
  x_mid <- mean(trajectory_x_limits)
  y_mid <- mean(trajectory_y_limits)

  x_direction <- if (root_point$DDRTree_1 >= x_mid) -1 else 1
  y_direction <- if (root_point$DDRTree_2 >= y_mid) -1 else 1

  root_label_x <- root_point$DDRTree_1 + x_direction * max(0.08 * x_span, 0.32)
  root_label_y <- root_point$DDRTree_2 + y_direction * max(0.07 * y_span, 0.28)

  root_label_x <- min(
    max(root_label_x, trajectory_x_limits[1] + 0.09 * x_span),
    trajectory_x_limits[2] - 0.09 * x_span
  )

  root_label_y <- min(
    max(root_label_y, trajectory_y_limits[1] + 0.09 * y_span),
    trajectory_y_limits[2] - 0.09 * y_span
  )

  root_label_data <- data.frame(
    x = root_label_x,
    y = root_label_y,
    label = paste0("Root\n(", root_state_label, ")"),
    stringsAsFactors = FALSE
  )

  root_connector_layer <- ggplot2::geom_segment(
    data = data.frame(
      x = root_point$DDRTree_1,
      y = root_point$DDRTree_2,
      xend = root_label_data$x,
      yend = root_label_data$y
    ),
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.34,
    alpha = 0.70,
    lineend = "round"
  )

  root_marker_layer <- ggplot2::geom_point(
    data = root_point,
    ggplot2::aes(x = DDRTree_1, y = DDRTree_2),
    inherit.aes = FALSE,
    shape = 8,
    size = 3.35,
    stroke = 0.90,
    colour = "black"
  )

  root_label_layer <- ggplot2::geom_label(
    data = root_label_data,
    ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    family = FONT_FAMILY,
    size = 2.75,
    fontface = "bold",
    label.size = 0.32,
    label.padding = grid::unit(0.12, "lines"),
    fill = "white",
    colour = "black",
    alpha = 0.98
  )

  common_trajectory_theme <- theme_manuscript(
    show_grid = FALSE,
    legend_position = "bottom"
  ) +
    ggplot2::theme(
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.justification = "center",
      legend.spacing.x = grid::unit(12, "pt"),
      legend.box.spacing = grid::unit(5, "pt")
    )

  p_pseudotime <- ggplot2::ggplot(
    trajectory_df,
    ggplot2::aes(
      x = DDRTree_1,
      y = DDRTree_2,
      colour = Pseudotime
    )
  ) +
    ggplot2::geom_point(size = 0.80, alpha = 0.88) +
    backbone_layer +
    root_connector_layer +
    root_marker_layer +
    root_label_layer +
    ggplot2::scale_colour_gradientn(
      colours = PSEUDOTIME_COLOURS,
      limits = pseudotime_limits,
      name = "Pseudotime"
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_colorbar(
        barwidth = grid::unit(3.1, "cm"),
        barheight = grid::unit(0.38, "cm"),
        title.position = "top"
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = trajectory_x_limits,
      ylim = trajectory_y_limits,
      expand = FALSE
    ) +
    ggplot2::labs(
      x = "DDRTree component 1",
      y = "DDRTree component 2"
    ) +
    common_trajectory_theme

  state_levels <- levels(trajectory_df$Monocle_state)

  if (length(state_levels) <= length(MONOCLE_STATE_COLOURS)) {
    state_colours <- MONOCLE_STATE_COLOURS[seq_along(state_levels)]
  } else {
    state_colours <- grDevices::hcl.colors(
      length(state_levels),
      palette = "Dynamic"
    )
  }
  names(state_colours) <- state_levels

  p_state <- ggplot2::ggplot(
    trajectory_df,
    ggplot2::aes(
      x = DDRTree_1,
      y = DDRTree_2,
      colour = Monocle_state
    )
  ) +
    ggplot2::geom_point(size = 0.80, alpha = 0.88) +
    backbone_layer +
    root_marker_layer +
    ggplot2::scale_colour_manual(
      values = state_colours,
      name = "States",
      drop = FALSE
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE,
        title.position = "top",
        override.aes = list(size = 2.1, alpha = 1)
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = trajectory_x_limits,
      ylim = trajectory_y_limits,
      expand = FALSE
    ) +
    ggplot2::labs(
      x = "DDRTree component 1",
      y = "DDRTree component 2"
    ) +
    common_trajectory_theme

  background_cells <- trajectory_df[
    !(as.character(trajectory_df$Malignant_subcluster) %in% IMPORTANT_CLUSTERS),
    ,
    drop = FALSE
  ]

  selected_cells <- trajectory_df[
    as.character(trajectory_df$Malignant_subcluster) %in% IMPORTANT_CLUSTERS,
    ,
    drop = FALSE
  ]

  p_selected <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = background_cells,
      ggplot2::aes(x = DDRTree_1, y = DDRTree_2),
      colour = COL_NS,
      size = 0.68,
      alpha = 0.35
    ) +
    ggplot2::geom_point(
      data = selected_cells,
      ggplot2::aes(
        x = DDRTree_1,
        y = DDRTree_2,
        colour = Malignant_subcluster
      ),
      size = 0.92,
      alpha = 0.95
    ) +
    backbone_layer +
    root_marker_layer +
    ggplot2::scale_colour_manual(
      values = IMPORTANT_CLUSTER_COLOURS,
      limits = IMPORTANT_CLUSTERS,
      name = "Subclusters",
      drop = FALSE
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE,
        title.position = "top",
        override.aes = list(size = 2.1, alpha = 1)
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = trajectory_x_limits,
      ylim = trajectory_y_limits,
      expand = FALSE
    ) +
    ggplot2::labs(
      x = "DDRTree component 1",
      y = "DDRTree component 2"
    ) +
    common_trajectory_theme

  top_trajectory_panels <- patchwork::wrap_plots(
    p_pseudotime,
    p_state,
    p_selected,
    ncol = 3,
    guides = "collect"
  ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.justification = "center",
      legend.spacing.x = grid::unit(12, "pt"),
      legend.box.spacing = grid::unit(5, "pt")
    )

  figure_3a <- top_trajectory_panels +
    patchwork::plot_annotation(tag_levels = "a")

  figure_3a <- figure_3a &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PANEL_TAG_PT,
        colour = "black"
      ),
      plot.tag.position = c(0.012, 0.988)
    )

  save_plot_all_formats(
    plot_obj = figure_3a,
    filename_stem = FIGURE_3A_STEM,
    dir_path = PSEUDOTIME_MONOCLE_DIR,
    width = FIGURE3_COMBINED_W,
    height = FIGURE3_COMBINED_H
  )

  expression_matrix <- Biobase::exprs(cds)

  if (!(target_gene %in% rownames(expression_matrix))) {
    stop(target_gene, " is absent from the Monocle expression matrix.")
  }

  size_factors <- BiocGenerics::sizeFactors(cds)

  if (is.null(size_factors)) {
    stop("Monocle size factors are missing from the CellDataSet.")
  }

  if (is.null(names(size_factors))) {
    names(size_factors) <- colnames(expression_matrix)
  }

  valid_cells <- trajectory_df$Cell[
    trajectory_df$Cell %in% colnames(expression_matrix) &
      is.finite(trajectory_df$Pseudotime)
  ]

  if (length(valid_cells) < 20L) {
    stop("Too few cells have finite pseudotime values.")
  }

  valid_cells <- valid_cells[
    order(
      trajectory_df$Pseudotime[
        match(valid_cells, trajectory_df$Cell)
      ]
    )
  ]

  valid_size_factors <- size_factors[valid_cells]

  if (
    any(!is.finite(valid_size_factors)) ||
      any(valid_size_factors <= 0)
  ) {
    stop("Invalid Monocle size factors were detected.")
  }

  mybl2_expression <- log1p(
    as.numeric(expression_matrix[target_gene, valid_cells]) /
      pmax(valid_size_factors, .Machine$double.eps)
  )

  mybl2_df <- data.frame(
    Cell = valid_cells,
    Pseudotime = trajectory_df$Pseudotime[
      match(valid_cells, trajectory_df$Cell)
    ],
    Malignant_subcluster = as.character(
      trajectory_df$Malignant_subcluster[
        match(valid_cells, trajectory_df$Cell)
      ]
    ),
    log1p_size_factor_normalised_MYBL2 = mybl2_expression,
    stringsAsFactors = FALSE
  )

  write.csv(mybl2_df, MYBL2_CELL_CSV, row.names = FALSE)

  mybl2_selected <- mybl2_df[
    mybl2_df$Malignant_subcluster %in% IMPORTANT_CLUSTERS &
      is.finite(mybl2_df$Pseudotime) &
      is.finite(mybl2_df$log1p_size_factor_normalised_MYBL2),
    ,
    drop = FALSE
  ]

  mybl2_selected$Malignant_subcluster <- factor(
    mybl2_selected$Malignant_subcluster,
    levels = IMPORTANT_CLUSTERS
  )

  if (
    !setequal(
      unique(as.character(mybl2_selected$Malignant_subcluster)),
      IMPORTANT_CLUSTERS
    )
  ) {
    stop("At least one selected subtype is missing from MYBL2 pseudotime data.")
  }

  loess_curve <- do.call(
    rbind,
    lapply(IMPORTANT_CLUSTERS, function(cluster_name) {
      make_loess_curve(
        data = mybl2_selected,
        cluster_name = cluster_name,
        span_value = LOESS_SPAN[[cluster_name]],
        ci_level = LOESS_CI_LEVEL
      )
    })
  )

  loess_curve$Malignant_subcluster <- factor(
    loess_curve$Malignant_subcluster,
    levels = IMPORTANT_CLUSTERS
  )

  write.csv(loess_curve, LOESS_CURVE_CSV, row.names = FALSE)

  supported_upper_ci <- loess_curve$LOESS_upper_CI[
    !is.na(loess_curve$Ribbon_supported) &
      loess_curve$Ribbon_supported
  ]

  mybl2_y_max <- max(
    mybl2_selected$log1p_size_factor_normalised_MYBL2,
    supported_upper_ci,
    na.rm = TRUE
  )

  if (!is.finite(mybl2_y_max) || mybl2_y_max <= 0) {
    mybl2_y_max <- max(
      mybl2_selected$log1p_size_factor_normalised_MYBL2,
      na.rm = TRUE
    )
  }

  if (!is.finite(mybl2_y_max) || mybl2_y_max <= 0) {
    mybl2_y_max <- 1
  }

  make_mybl2_panel <- function(cluster_name) {
    point_data <- mybl2_selected[
      as.character(mybl2_selected$Malignant_subcluster) == cluster_name,
      ,
      drop = FALSE
    ]

    line_data <- loess_curve[
      as.character(loess_curve$Malignant_subcluster) == cluster_name,
      ,
      drop = FALSE
    ]

    ribbon_data <- line_data[
      !is.na(line_data$Ribbon_supported) &
        line_data$Ribbon_supported &
        is.finite(line_data$LOESS_lower_CI) &
        is.finite(line_data$LOESS_upper_CI),
      ,
      drop = FALSE
    ]

    cluster_colour <- IMPORTANT_CLUSTER_COLOURS[[cluster_name]]
    panel_title <- paste0(
      cluster_name,
      " (n = ",
      nrow(point_data),
      ")"
    )

    panel_x_range <- range(point_data$Pseudotime, na.rm = TRUE)
    panel_x_span <- diff(panel_x_range)

    panel_x_padding <- if (is.finite(panel_x_span) && panel_x_span > 0) {
      max(0.025 * panel_x_span, 0.06)
    } else {
      0.50
    }

    panel_x_limits <- c(
      panel_x_range[1] - panel_x_padding,
      panel_x_range[2] + panel_x_padding
    )

    ggplot2::ggplot(
      point_data,
      ggplot2::aes(
        x = Pseudotime,
        y = log1p_size_factor_normalised_MYBL2
      )
    ) +
      ggplot2::geom_point(
        colour = cluster_colour,
        alpha = 0.20,
        size = 0.78
      ) +
      ggplot2::geom_ribbon(
        data = ribbon_data,
        ggplot2::aes(
          x = Pseudotime,
          ymin = LOESS_lower_CI,
          ymax = LOESS_upper_CI,
          group = Ribbon_segment
        ),
        inherit.aes = FALSE,
        fill = cluster_colour,
        alpha = 0.16,
        na.rm = TRUE
      ) +
      ggplot2::geom_line(
        data = line_data,
        ggplot2::aes(
          x = Pseudotime,
          y = LOESS_smoothed_log1p_expression
        ),
        inherit.aes = FALSE,
        colour = cluster_colour,
        linewidth = max(GEOM_LWD, 0.85),
        na.rm = TRUE
      ) +
      ggplot2::coord_cartesian(
        xlim = panel_x_limits,
        ylim = c(0, mybl2_y_max * 1.05),
        expand = FALSE
      ) +
      ggplot2::labs(
        title = panel_title,
        x = "Inferred pseudotime",
        y = "log1p(MYBL2 expression)"
      ) +
      theme_manuscript(
        show_grid = FALSE,
        legend_position = "none"
      )
  }

  p_cluster_0 <- make_mybl2_panel("Cluster 0")
  p_cluster_8 <- make_mybl2_panel("Cluster 8")

  figure_3b <- (
    p_cluster_0 |
      p_cluster_8
  ) +
    patchwork::plot_annotation(tag_levels = "a")

  figure_3b <- figure_3b &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PANEL_TAG_PT,
        colour = "black"
      ),
      plot.tag.position = c(0.012, 0.988)
    )

  save_plot_all_formats(
    plot_obj = figure_3b,
    filename_stem = FIGURE_3B_STEM,
    dir_path = PSEUDOTIME_MONOCLE_DIR,
    width = FIGURE3_COMBINED_W,
    height = FIGURE3_COMBINED_H
  )

  bottom_mybl2_panels <- (
    p_cluster_0 |
      p_cluster_8
  )

  figure_8 <- (
    top_trajectory_panels /
      bottom_mybl2_panels
  ) +
    patchwork::plot_layout(heights = c(1.00, 0.82)) +
    patchwork::plot_annotation(tag_levels = "a")

  figure_8 <- figure_8 &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PANEL_TAG_PT,
        colour = "black"
      ),
      plot.tag.position = c(0.012, 0.988)
    )

  save_plot_all_formats(
    plot_obj = figure_8,
    filename_stem = FIGURE_8_STEM,
    dir_path = PSEUDOTIME_MONOCLE_DIR,
    width = FIGURE8_COMBINED_W,
    height = FIGURE8_COMBINED_H
  )

  invisible(
    list(
      Figure_3A = figure_3a,
      Figure_3B = figure_3b,
      Figure_8 = figure_8
    )
  )
}

# 10. STEP 2: RUN DDRTREE MONOCLE AND CREATE ALL FIGURES

run_trajectory_and_figures <- function() {
  cat(
    "\n============================================================\n",
    "STEP 2 / 2: MONOCLE 2 DDRTREE + ALL FIGURES IN RGui 4.4.3\n",
    "============================================================\n",
    sep = ""
  )

  if (!file.exists(INPUT_RDS)) {
    stop(
      "Prepared Monocle input is missing:\n",
      INPUT_RDS,
      "\n\nRun PART 1 (Prepare Monocle Input) first."
    )
  }

  load_monocle_environment()

  old_outputs <- c(
    STATUS_TXT, ERROR_TXT, RUN_LOG, SESSION_TXT, COMPATIBILITY_TXT,
    REDUCTION_TXT, CDS_RDS, TRAJECTORY_CSV, ORDERING_GENES_TRAJECTORY_CSV,
    ROOT_CSV, ROOT_TXT, MYBL2_CELL_CSV, LOESS_CURVE_CSV, BACKBONE_STATUS_TXT,
    BACKBONE_SEGMENTS_CSV,
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_3A_STEM, ".png")),
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_3A_STEM, ".tiff")),
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_3A_STEM, ".pdf")),
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_3B_STEM, ".png")),
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_3B_STEM, ".tiff")),
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_3B_STEM, ".pdf")),
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_8_STEM, ".png")),
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_8_STEM, ".tiff")),
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_8_STEM, ".pdf"))
  )

  unlink(old_outputs[file.exists(old_outputs)], force = TRUE)

  writeLines(
    c(
      "This RGui runner restores a pooled Monocle-compatible blind dispersion",
      "model before DDRTree.",
      "The native Monocle estimateDispersions() function is not called because",
      "it relies on legacy dplyr calls that can fail with current dplyr versions."
    ),
    COMPATIBILITY_TXT
  )

  result <- list(
    trajectory_success = FALSE,
    figure_success = FALSE,
    reduction_fallback_used = FALSE,
    error = NA_character_,
    figure_error = NA_character_,
    stage = "started",
    root_state = NA_integer_,
    n_cells = 0L,
    n_ordering_genes = 0L,
    monocle_version = as.character(utils::packageVersion("monocle")),
    igraph_version = as.character(utils::packageVersion("igraph"))
  )

  write_status <- function() {
    writeLines(
      c(
        "LSCC Monocle 2 pseudotime status",
        "Run mode: fresh plain R 4.4.3 (RGui).",
        "Pooled compatibility dispersion model: enabled.",
        "Dynamic-gene test: intentionally omitted in this trajectory-and-figure runner.",
        paste0("Final stage: ", result$stage),
        paste0("Trajectory completed: ", result$trajectory_success),
        paste0("Figures completed: ", result$figure_success),
        paste0("Reduction fallback used: ", result$reduction_fallback_used),
        paste0("Cells used: ", result$n_cells),
        paste0("Ordering genes used: ", result$n_ordering_genes),
        paste0("Root state: ", result$root_state),
        paste0("Monocle version: ", result$monocle_version),
        paste0("igraph version: ", result$igraph_version),
        paste0(
          "Core error: ",
          ifelse(is.na(result$error), "none", result$error)
        ),
        paste0(
          "Figure error: ",
          ifelse(is.na(result$figure_error), "none", result$figure_error)
        )
      ),
      STATUS_TXT
    )

    writeLines(capture.output(sessionInfo()), SESSION_TXT)
  }

  cds <- NULL
  target_gene <- NULL

  tryCatch(
    {
      result$stage <- "loading_input"
      log_message("STAGE: loading prepared Monocle input")

      input <- readRDS(INPUT_RDS)

      required_items <- c(
        "count_mat", "pd_df", "ordering_genes",
        "target_gene", "root_column", "expected_subclusters"
      )

      if (!is.list(input) || !all(required_items %in% names(input))) {
        stop("Prepared Monocle input RDS is incomplete.")
      }

      counts <- as(input$count_mat, "dgCMatrix")
      assert_raw_counts(counts)

      pd <- as.data.frame(input$pd_df, stringsAsFactors = FALSE)
      ordering_genes <- intersect(
        as.character(input$ordering_genes),
        rownames(counts)
      )
      target_gene <- as.character(input$target_gene)
      root_column <- as.character(input$root_column)
      expected_subclusters <- as.character(input$expected_subclusters)

      if (!identical(rownames(pd), colnames(counts))) {
        stop("Input metadata and count matrix are not aligned.")
      }

      if (
        !all(
          c("sample", "Malignant_subcluster", root_column) %in% names(pd)
        )
      ) {
        stop("Required metadata columns are missing.")
      }

      if (!(target_gene %in% rownames(counts))) {
        stop(target_gene, " is absent from raw counts.")
      }

      if (length(ordering_genes) < MIN_ORDERING_GENES) {
        stop("Too few validated ordering genes: ", length(ordering_genes))
      }

      pd$sample <- factor(as.character(pd$sample))
      pd$Malignant_subcluster <- factor(
        normalise_cluster(pd$Malignant_subcluster),
        levels = expected_subclusters
      )

      if (
        any(is.na(pd$Malignant_subcluster)) ||
          !setequal(
            unique(as.character(pd$Malignant_subcluster)),
            expected_subclusters
          )
      ) {
        stop("The expected ten High-CNV subclusters are not all present.")
      }

      if (any(!is.finite(pd[[root_column]]))) {
        stop("The G2/M root score contains invalid values.")
      }

      result$n_cells <- ncol(counts)
      result$n_ordering_genes <- length(ordering_genes)

      result$stage <- "creating_CellDataSet"
      log_message("STAGE: creating CellDataSet")

      feature_data <- data.frame(
        gene_short_name = rownames(counts),
        row.names = rownames(counts),
        stringsAsFactors = FALSE
      )

      cds <- monocle::newCellDataSet(
        counts,
        phenoData = Biobase::AnnotatedDataFrame(pd),
        featureData = Biobase::AnnotatedDataFrame(feature_data),
        lowerDetectionLimit = 0.5,
        expressionFamily = VGAM::negbinomial.size()
      )

      result$stage <- "estimating_size_factors"
      log_message("STAGE: estimating size factors via BiocGenerics S4 generic")

      cds <- BiocGenerics::estimateSizeFactors(cds)

      result$stage <- "detecting_genes"
      log_message("STAGE: detecting expressed genes")

      cds <- monocle::detectGenes(cds, min_expr = 0.5)

      result$stage <- "compatibility_dispersion"
      log_message("STAGE: creating pooled Monocle-compatible dispersion model")

      cds <- estimate_dispersions_blind_compat(
        cds,
        min_cells_detected = 1L,
        remove_outliers = TRUE,
        verbose = FALSE
      )

      cds <- monocle::setOrderingFilter(cds, ordering_genes)

      result$stage <- "DDRTree_reduction"
      log_message("STAGE: running DDRTree with sample adjustment")

      adjusted_arguments <- list(
        cds = cds,
        max_components = 2,
        reduction_method = "DDRTree",
        norm_method = "log",
        residualModelFormulaStr = "~sample",
        verbose = FALSE
      )

      reduction_note <- "Sample-adjusted DDRTree reduction completed."

      cds <- tryCatch(
        do.call(monocle::reduceDimension, adjusted_arguments),
        error = function(e) {
          result$reduction_fallback_used <<- TRUE

          reduction_note <<- paste0(
            "Sample-adjusted DDRTree failed; fallback without residual model used: ",
            conditionMessage(e)
          )

          log_message(reduction_note)

          fallback_arguments <- list(
            cds = cds,
            max_components = 2,
            reduction_method = "DDRTree",
            norm_method = "log",
            verbose = FALSE
          )

          do.call(monocle::reduceDimension, fallback_arguments)
        }
      )

      writeLines(reduction_note, REDUCTION_TXT)

      result$stage <- "initial_ordering"
      log_message("STAGE: initial cell ordering")
      write_ordercells_geometry_diagnostics(cds, "before_initial_orderCells")

      cds <- monocle::orderCells(cds)

      result$stage <- "root_selection"
      log_message("STAGE: selecting root state using lowest median G2/M score")

      cds_metadata <- as.data.frame(
        Biobase::pData(cds),
        stringsAsFactors = FALSE
      )

      state_values <- suppressWarnings(
        as.integer(as.character(cds_metadata$State))
      )

      state_medians <- tapply(
        cds_metadata[[root_column]],
        state_values,
        median,
        na.rm = TRUE
      )
      state_medians <- state_medians[is.finite(state_medians)]

      if (length(state_medians) == 0L) {
        stop("No valid root-state G2/M medians were produced.")
      }

      root_state <- as.integer(
        names(state_medians)[which.min(state_medians)]
      )

      write.csv(
        data.frame(
          State = as.integer(names(state_medians)),
          Median_G2M_root_score = as.numeric(state_medians),
          Is_selected_root = as.integer(names(state_medians)) == root_state
        ),
        ROOT_CSV,
        row.names = FALSE
      )

      write_ordercells_geometry_diagnostics(cds, paste0("before_rooted_orderCells_root_state_", root_state))
      cds <- monocle::orderCells(cds, root_state = root_state)

      result$root_state <- root_state
      result$trajectory_success <- TRUE
      result$stage <- "trajectory_completed"

      saveRDS(cds, CDS_RDS)

      write.csv(
        data.frame(Gene = ordering_genes),
        ORDERING_GENES_TRAJECTORY_CSV,
        row.names = FALSE
      )

      writeLines(
        c(
          paste0("Root state: ", root_state),
          "Root rule: state with the lowest median G2M_root_score.",
          "Pseudotime is an inferred transcriptional continuum, not observed time.",
          reduction_note
        ),
        ROOT_TXT
      )

      log_message("TRAJECTORY SUCCESSFUL | Root state = ", root_state)
    },
    error = function(e) {
      result$error <<- conditionMessage(e)
      result$stage <<- paste0("failed_at_", result$stage)

      writeLines(
        c(
          "CORE ERROR",
          paste0("Stage: ", result$stage),
          paste0("Message: ", conditionMessage(e)),
          "",
          "Call stack:",
          capture.output(sys.calls())
        ),
        ERROR_TXT
      )

      log_message("CORE ERROR: ", conditionMessage(e))
    }
  )

  if (isTRUE(result$trajectory_success)) {
    tryCatch(
      {
        result$stage <- "creating_figures"
        log_message("STAGE: creating Figure 3A, Figure 3B and Figure 8")

        create_all_figures(
          cds = cds,
          target_gene = target_gene,
          root_state = result$root_state
        )

        result$figure_success <- TRUE
        result$stage <- "figures_completed"

        log_message("FIGURES SUCCESSFUL")
      },
      error = function(e) {
        result$figure_error <<- conditionMessage(e)
        log_message("FIGURE ERROR: ", conditionMessage(e))
      }
    )
  }

  write_status()

  cat(
    "\n============================================================\n",
    "MONOCLE 2 RUN FINISHED\n",
    "Trajectory completed: ", result$trajectory_success, "\n",
    "Figures completed: ", result$figure_success, "\n",
    "Output folder:\n", PSEUDOTIME_MONOCLE_DIR,
    "\n============================================================\n",
    sep = ""
  )

  if (!isTRUE(result$trajectory_success) || !isTRUE(result$figure_success)) {
    warning(
      "The run did not fully complete. Read these files inside Pseudotime_Monocle:\n",
      "Monocle2_error.txt\n",
      "Monocle2_status.txt\n",
      "Monocle2_run.log"
    )
  }
}

# EXECUTE PART 2

FIGURE_REFRESH_STATUS_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_figure_refresh_status.txt"
)

find_prepared_input_rds <- function() {
  search_roots <- unique(
    c(
      PSEUDOTIME_MONOCLE_DIR,
      "D:/LSCC/Results_LSCC",
      "D:/LSCC/ScRNAseq_Results/GSE206332/Results"
    )
  )

  direct_candidates <- c(
    INPUT_RDS,
    file.path(
      PSEUDOTIME_MONOCLE_DIR,
      "Monocle2_full_analysis_input.RDS"
    )
  )

  discovered_candidates <- unlist(
    lapply(search_roots, function(root_dir) {
      if (!dir.exists(root_dir)) {
        return(character(0))
      }

      list.files(
        root_dir,
        pattern = "^Monocle2_full_analysis_input\\.(rds|RDS)$",
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
      )
    }),
    use.names = FALSE
  )

  candidates <- unique(c(direct_candidates, discovered_candidates))
  candidates <- candidates[file.exists(candidates)]

  if (length(candidates) == 0L) {
    return(NA_character_)
  }

  normalizePath(candidates[1], winslash = "/", mustWork = TRUE)
}

infer_root_state_for_figure_refresh <- function(cds_object) {
  if (file.exists(ROOT_CSV)) {
    root_table <- tryCatch(
      read.csv(ROOT_CSV, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )

    if (
      !is.null(root_table) &&
      all(c("State", "Is_selected_root") %in% colnames(root_table))
    ) {
      selected <- root_table$State[
        !is.na(root_table$Is_selected_root) &
          as.logical(root_table$Is_selected_root)
      ]

      if (length(selected) >= 1L && is.finite(selected[1])) {
        return(as.integer(selected[1]))
      }
    }
  }

  metadata <- as.data.frame(
    Biobase::pData(cds_object),
    stringsAsFactors = FALSE
  )

  state_values <- suppressWarnings(
    as.integer(as.character(metadata$State))
  )

  if ("G2M_root_score" %in% colnames(metadata)) {
    state_medians <- tapply(
      suppressWarnings(as.numeric(metadata$G2M_root_score)),
      state_values,
      median,
      na.rm = TRUE
    )

    state_medians <- state_medians[is.finite(state_medians)]

    if (length(state_medians) > 0L) {
      return(as.integer(names(state_medians)[which.min(state_medians)]))
    }
  }

  pseudotime <- suppressWarnings(as.numeric(metadata$Pseudotime))

  if (any(is.finite(pseudotime) & is.finite(state_values))) {
    return(
      state_values[
        which.min(replace(pseudotime, !is.finite(pseudotime), Inf))
      ]
    )
  }

  stop(
    "The root state could not be recovered from the saved CellDataSet."
  )
}

refresh_figures_from_completed_cds <- function() {
  cat(
    "\n============================================================\n",
    "AUTO-RESUME MODE: REDRAWING FIGURES FROM COMPLETED CellDataSet\n",
    "============================================================\n",
    sep = ""
  )

  if (!file.exists(CDS_RDS)) {
    stop("Completed CellDataSet is missing:\n", CDS_RDS)
  }

  load_monocle_environment()

  cds <- readRDS(CDS_RDS)

  if (!inherits(cds, "CellDataSet")) {
    stop("Saved HighCNV_Monocle2_CDS.rds is not a valid Monocle CellDataSet.")
  }

  expression_matrix <- Biobase::exprs(cds)

  if (!(TARGET_GENE %in% rownames(expression_matrix))) {
    stop(
      TARGET_GENE,
      " is absent from the saved completed CellDataSet."
    )
  }

  metadata <- as.data.frame(
    Biobase::pData(cds),
    stringsAsFactors = FALSE
  )

  required_metadata <- c(
    "Pseudotime",
    "State",
    "sample",
    "Malignant_subcluster"
  )

  missing_metadata <- setdiff(required_metadata, colnames(metadata))

  if (length(missing_metadata) > 0L) {
    stop(
      "Saved completed CellDataSet is missing: ",
      paste(missing_metadata, collapse = ", ")
    )
  }

  recovered_root_state <- infer_root_state_for_figure_refresh(cds)

  log_message(
    "AUTO-RESUME: using saved completed CellDataSet; no trajectory is recomputed."
  )
  log_message(
    "AUTO-RESUME: recovered root state = ",
    recovered_root_state
  )

  create_all_figures(
    cds = cds,
    target_gene = TARGET_GENE,
    root_state = recovered_root_state
  )

  writeLines(
    c(
      "SUCCESS: Figures were refreshed from an existing completed Monocle CellDataSet.",
      paste0("Saved CellDataSet: ", normalizePath(CDS_RDS, winslash = "/", mustWork = TRUE)),
      paste0("Recovered root state: ", recovered_root_state),
      "Trajectory was not recomputed.",
      "Updated features: in-panel root-state label, DDRTree backbone overlay,",
      "support-trimmed LOESS 95% confidence ribbons, a compact root annotation, and per-cluster pseudotime x-limits.",
      paste0("Output folder: ", PSEUDOTIME_MONOCLE_DIR)
    ),
    FIGURE_REFRESH_STATUS_TXT
  )

  cat(
    "\n============================================================\n",
    "FIGURE REFRESH FINISHED SUCCESSFULLY\n",
    "Trajectory was reused; figures were regenerated.\n",
    "Output folder:\n",
    PSEUDOTIME_MONOCLE_DIR,
    "\n============================================================\n",
    sep = ""
  )

  invisible(TRUE)
}

# AUTO-REBUILD OF MISSING PART-1 INPUT

AUTO_REBUILD_STATUS_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_auto_rebuild_input_status.txt"
)

find_first_existing_rebuild <- function(candidates, label) {
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0L) {
    stop(
      label, " was not found. Checked:\n",
      paste(candidates, collapse = "\n")
    )
  }

  normalizePath(existing[1], winslash = "/", mustWork = TRUE)
}

rebuild_cluster_number <- function(x) {
  x <- trimws(as.character(x))
  out <- suppressWarnings(as.integer(sub("^.*?([0-9]+)\\s*$", "\\1", x)))
  out[!grepl("[0-9]+\\s*$", x)] <- NA_integer_
  out
}

rebuild_normalise_cluster <- function(x) {
  x <- trimws(as.character(x))
  n <- rebuild_cluster_number(x)
  x[!is.na(n)] <- paste0("Cluster ", n[!is.na(n)])
  x
}

rebuild_assert_raw_counts <- function(x) {
  if (!inherits(x, "Matrix")) {
    stop("Auto-rebuild extracted a non-sparse count matrix.")
  }

  if (length(x@x) == 0L) {
    stop("Auto-rebuild extracted an empty count matrix.")
  }

  if (any(!is.finite(x@x)) || any(x@x < 0)) {
    stop("Auto-rebuild count matrix contains invalid values.")
  }

  if (any(abs(x@x - round(x@x)) > 1e-8)) {
    stop(
      "Auto-rebuild did not obtain raw integer counts. ",
      "Monocle 2 requires raw count data."
    )
  }

  invisible(TRUE)
}

rebuild_extract_joined_counts <- function(object, assay, layer = "counts") {
  if (!(assay %in% names(object@assays))) {
    stop("Assay not found during auto-rebuild: ", assay)
  }

  layer_names <- tryCatch(
    SeuratObject::Layers(object[[assay]]),
    error = function(e) character(0)
  )

  layer_names <- layer_names[
    grepl(paste0("^", layer, "(\\.|$)"), layer_names)
  ]

  if (length(layer_names) <= 1L) {
    count_matrix <- NULL

    if (length(layer_names) == 1L) {
      count_matrix <- tryCatch(
        SeuratObject::LayerData(
          object,
          assay = assay,
          layer = layer_names[1],
          fast = FALSE
        ),
        error = function(e) NULL
      )
    }

    if (is.null(count_matrix) || ncol(count_matrix) == 0L) {
      count_matrix <- tryCatch(
        Seurat::GetAssayData(
          object,
          assay = assay,
          slot = layer
        ),
        error = function(e) NULL
      )
    }

    if (!is.null(count_matrix) && ncol(count_matrix) > 0L) {
      if (
        length(intersect(rownames(count_matrix), colnames(object))) >
          length(intersect(colnames(count_matrix), colnames(object)))
      ) {
        count_matrix <- Matrix::t(count_matrix)
      }

      retained <- intersect(colnames(object), colnames(count_matrix))

      return(
        as(
          count_matrix[, retained, drop = FALSE],
          "dgCMatrix"
        )
      )
    }
  }

  layer_matrices <- lapply(layer_names, function(layer_name) {
    x <- SeuratObject::LayerData(
      object,
      assay = assay,
      layer = layer_name,
      fast = FALSE
    )

    if (
      length(intersect(rownames(x), colnames(object))) >
        length(intersect(colnames(x), colnames(object)))
    ) {
      x <- Matrix::t(x)
    }

    retained <- intersect(colnames(x), colnames(object))
    x[, retained, drop = FALSE]
  })

  layer_matrices <- layer_matrices[
    vapply(layer_matrices, ncol, integer(1)) > 0L
  ]

  if (length(layer_matrices) == 0L) {
    stop("No usable counts layer was found in assay ", assay)
  }

  shared_genes <- Reduce(intersect, lapply(layer_matrices, rownames))

  if (length(shared_genes) < 2L) {
    stop("Too few shared genes across count layers during auto-rebuild.")
  }

  combined <- do.call(
    cbind,
    lapply(
      layer_matrices,
      function(x) x[shared_genes, , drop = FALSE]
    )
  )

  combined <- combined[, !duplicated(colnames(combined)), drop = FALSE]

  retained <- intersect(colnames(object), colnames(combined))

  as(combined[, retained, drop = FALSE], "dgCMatrix")
}

rebuild_marker_columns <- function(marker_table) {
  gene_col <- intersect(
    c("Gene", "gene", "Gene.symbol", "gene_name", "Symbol", "symbol"),
    colnames(marker_table)
  )[1]

  cluster_col <- intersect(
    c("cluster", "Cluster", "Cluster_label", "Malignant_subcluster"),
    colnames(marker_table)
  )[1]

  padj_col <- intersect(
    c("p_val_adj", "p_adj", "adj.P.Val", "FDR"),
    colnames(marker_table)
  )[1]

  fc_col <- intersect(
    c("avg_log2FC", "avg_logFC", "logFC"),
    colnames(marker_table)
  )[1]

  pct_col <- intersect(
    c("pct.1", "pct1", "Percent.1", "pct_in"),
    colnames(marker_table)
  )[1]

  output <- c(
    Gene = gene_col,
    cluster = cluster_col,
    p_val_adj = padj_col,
    avg_logFC = fc_col,
    pct.1 = pct_col
  )

  if (any(is.na(output) | !nzchar(output))) {
    stop(
      "The auto-rebuild marker CSV lacks a required column.\n",
      "Available columns: ", paste(colnames(marker_table), collapse = ", ")
    )
  }

  output
}

rebuild_g2m_score <- function(count_matrix) {
  genes <- intersect(G2M_ROOT_GENE_SET, rownames(count_matrix))

  if (length(genes) < 5L) {
    stop(
      "Fewer than five G2/M root-score genes were found in raw counts. ",
      "Detected: ", paste(genes, collapse = ", ")
    )
  }

  library_size <- Matrix::colSums(count_matrix)

  normalised <- log1p(
    sweep(
      as.matrix(count_matrix[genes, , drop = FALSE]),
      2,
      pmax(library_size, 1),
      "/"
    ) * 10000
  )

  colMeans(normalised)
}

rebuild_prepared_input_from_real_sources <- function() {
  cat(
    "\n============================================================\n",
    "AUTO-REBUILD MODE: RECREATING MISSING MONOCLE INPUT\n",
    "============================================================\n",
    sep = ""
  )

  required_rebuild_packages <- c("Seurat", "SeuratObject", "Matrix")

  missing_rebuild_packages <- required_rebuild_packages[
    !vapply(
      required_rebuild_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(missing_rebuild_packages) > 0L) {
    stop(
      "Auto-rebuild requires these packages in R 4.4.3:\n",
      paste(missing_rebuild_packages, collapse = ", ")
    )
  }

  sc_project_dir <- "D:/LSCC/ScRNAseq_Results/GSE206332"
  sc_results_dir <- file.path(sc_project_dir, "Results")
  sc_rds_dir <- file.path(sc_results_dir, "rds")
  sc_fig_dir <- file.path(sc_results_dir, "figures")

  highcnv_metadata_rds <- find_first_existing_rebuild(
    c(
      file.path(
        sc_project_dir,
        "Results_Modular",
        "Step_20_Final_Exports_and_Objects",
        "rds",
        "HIGH_CNV_MALIGNANT_object_after_lowCNV_removal.rds"
      ),
      file.path(
        sc_rds_dir,
        "HIGH_CNV_MALIGNANT_object_after_lowCNV_removal.rds"
      )
    ),
    "High-CNV malignant metadata Seurat object"
  )

  highcnv_de_rds <- find_first_existing_rebuild(
    c(
      file.path(
        sc_project_dir,
        "Results_Modular",
        "Step_20_Final_Exports_and_Objects",
        "rds",
        "HIGH_CNV_MALIGNANT_object_DE_fixed.rds"
      ),
      file.path(
        sc_rds_dir,
        "HIGH_CNV_MALIGNANT_object_DE_fixed.rds"
      )
    ),
    "High-CNV malignant DE Seurat object"
  )

  marker_csv <- find_first_existing_rebuild(
    c(
      file.path(
        sc_fig_dir,
        "High_CNV_Malignant_Cluster_FindAllMarkers_after_lowCNV_removed_logFC1_all.csv"
      ),
      file.path(
        sc_fig_dir,
        "High_CNV_Malignant_FindAllMarkers_after_lowCNV_removed_logFC1.csv"
      ),
      file.path(
        sc_fig_dir,
        "Final_High_CNV_Malignant_Genes.csv"
      ),
      file.path(
        sc_fig_dir,
        "High_CNV_Malignant_specific_markers_FINAL.csv"
      )
    ),
    "Precomputed High-CNV marker CSV"
  )

  log_message("AUTO-REBUILD: loading real High-CNV source objects")

  highcnv_obj <- readRDS(highcnv_metadata_rds)
  highcnv_de_obj <- readRDS(highcnv_de_rds)

  if (!inherits(highcnv_obj, "Seurat") || !inherits(highcnv_de_obj, "Seurat")) {
    stop("Auto-rebuild source RDS files are not valid Seurat objects.")
  }

  shared_cells <- intersect(
    colnames(highcnv_de_obj),
    colnames(highcnv_obj)
  )

  if (length(shared_cells) < 100L) {
    stop("Fewer than 100 shared cells were found between High-CNV source objects.")
  }

  assay_name <- if ("RNA_DE_HIGH_CNV" %in% names(highcnv_de_obj@assays)) {
    "RNA_DE_HIGH_CNV"
  } else if ("RNA" %in% names(highcnv_de_obj@assays)) {
    "RNA"
  } else {
    stop("Neither RNA_DE_HIGH_CNV nor RNA is available in the DE source object.")
  }

  counts_all <- rebuild_extract_joined_counts(
    highcnv_de_obj,
    assay = assay_name,
    layer = "counts"
  )

  rebuild_assert_raw_counts(counts_all)

  retained_cells <- intersect(shared_cells, colnames(counts_all))

  if (length(retained_cells) < 100L) {
    stop("Fewer than 100 cells overlap between raw counts and High-CNV metadata.")
  }

  counts_all <- counts_all[, retained_cells, drop = FALSE]

  metadata <- as.data.frame(
    highcnv_obj@meta.data[retained_cells, , drop = FALSE],
    stringsAsFactors = FALSE
  )

  required_metadata <- c("sample", "Cluster_label", "CNV_class")
  missing_metadata <- setdiff(required_metadata, colnames(metadata))

  if (length(missing_metadata) > 0L) {
    stop(
      "High-CNV metadata object lacks: ",
      paste(missing_metadata, collapse = ", ")
    )
  }

  if (any(
    is.na(metadata$CNV_class) |
      as.character(metadata$CNV_class) != "High-CNV malignant"
  )) {
    stop(
      "The selected source object includes cells not labelled as High-CNV malignant."
    )
  }

  metadata$sample <- as.character(metadata$sample)
  metadata$Malignant_subcluster <- factor(
    rebuild_normalise_cluster(metadata$Cluster_label),
    levels = FIXED_CLUSTER_ORDER
  )

  if (any(is.na(metadata$Malignant_subcluster))) {
    bad <- unique(
      as.character(metadata$Cluster_label[
        is.na(metadata$Malignant_subcluster)
      ])
    )

    stop(
      "Unrecognised Cluster_label values in auto-rebuild:\n",
      paste(bad, collapse = ", ")
    )
  }

  if (!setequal(
    unique(as.character(metadata$Malignant_subcluster)),
    FIXED_CLUSTER_ORDER
  )) {
    stop("The expected ten High-CNV subclusters are not all present.")
  }

  metadata$sample_subcluster <- paste(
    metadata$sample,
    metadata$Malignant_subcluster,
    sep = " || "
  )

  before_balance <- as.data.frame(
    table(
      Sample = metadata$sample,
      Subcluster = metadata$Malignant_subcluster
    )
  )
  names(before_balance)[3] <- "Cells_before"

  grouped_indices <- split(
    seq_len(nrow(metadata)),
    metadata$sample_subcluster
  )

  keep_indices <- sort(
    unlist(
      lapply(grouped_indices, function(index) {
        if (length(index) > MAX_CELLS_PER_SAMPLE_SUBCLUSTER) {
          sample(index, MAX_CELLS_PER_SAMPLE_SUBCLUSTER)
        } else {
          index
        }
      }),
      use.names = FALSE
    )
  )

  retained_cells <- rownames(metadata)[keep_indices]
  counts <- counts_all[, retained_cells, drop = FALSE]
  metadata <- metadata[retained_cells, , drop = FALSE]

  if (!identical(colnames(counts), rownames(metadata))) {
    stop("Counts and metadata could not be aligned after balanced sampling.")
  }

  after_balance <- as.data.frame(
    table(
      Sample = metadata$sample,
      Subcluster = metadata$Malignant_subcluster
    )
  )
  names(after_balance)[3] <- "Cells_after"

  balance_table <- merge(
    before_balance,
    after_balance,
    by = c("Sample", "Subcluster"),
    all = TRUE
  )
  balance_table[is.na(balance_table)] <- 0L
  write.csv(balance_table, CELL_BALANCE_CSV, row.names = FALSE)

  markers <- read.csv(
    marker_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (nrow(markers) == 0L) {
    stop("The precomputed High-CNV marker CSV has zero rows.")
  }

  marker_columns <- rebuild_marker_columns(markers)

  marker_data <- data.frame(
    Gene = as.character(markers[[marker_columns["Gene"]]]),
    cluster = as.character(markers[[marker_columns["cluster"]]]),
    p_val_adj = suppressWarnings(
      as.numeric(markers[[marker_columns["p_val_adj"]]])
    ),
    avg_logFC = suppressWarnings(
      as.numeric(markers[[marker_columns["avg_logFC"]]])
    ),
    pct.1 = suppressWarnings(
      as.numeric(markers[[marker_columns["pct.1"]]])
    ),
    stringsAsFactors = FALSE
  )

  marker_data$cluster_label <- rebuild_normalise_cluster(marker_data$cluster)

  marker_data <- marker_data[
    !is.na(marker_data$Gene) &
      nzchar(marker_data$Gene) &
      is.finite(marker_data$p_val_adj) &
      marker_data$p_val_adj < 0.05 &
      is.finite(marker_data$avg_logFC) &
      marker_data$avg_logFC >= 1 &
      is.finite(marker_data$pct.1) &
      marker_data$pct.1 >= 0.25 &
      marker_data$cluster_label %in% FIXED_CLUSTER_ORDER &
      marker_data$Gene %in% rownames(counts),
    ,
    drop = FALSE
  ]

  if (nrow(marker_data) == 0L) {
    stop(
      "No marker genes passed FDR < 0.05, logFC >= 1, pct.1 >= 0.25 ",
      "after matching to the raw count matrix."
    )
  }

  marker_data$cluster_order <- match(
    marker_data$cluster_label,
    FIXED_CLUSTER_ORDER
  )

  marker_data <- marker_data[
    order(
      marker_data$cluster_order,
      marker_data$p_val_adj,
      -marker_data$avg_logFC,
      -marker_data$pct.1,
      marker_data$Gene
    ),
    ,
    drop = FALSE
  ]

  ordering_by_cluster <- lapply(FIXED_CLUSTER_ORDER, function(cluster_label) {
    genes <- marker_data$Gene[
      marker_data$cluster_label == cluster_label
    ]
    unique(head(genes, MAX_ORDERING_GENES_PER_SUBCLUSTER))
  })

  ordering_genes <- unique(unlist(ordering_by_cluster, use.names = FALSE))
  ordering_genes <- ordering_genes[
    !grepl("^(MT-|RPL|RPS)", ordering_genes)
  ]

  if (length(ordering_genes) < MIN_ORDERING_GENES) {
    stop(
      "Too few ordering genes remained after filtering: ",
      length(ordering_genes)
    )
  }

  write.csv(
    data.frame(Gene = ordering_genes),
    ORDERING_GENES_INPUT_CSV,
    row.names = FALSE
  )

  if (!(TARGET_GENE %in% rownames(counts))) {
    stop(TARGET_GENE, " is absent from the rebuilt raw count matrix.")
  }

  metadata$G2M_root_score <- rebuild_g2m_score(counts)

  if (any(!is.finite(metadata$G2M_root_score))) {
    stop("Auto-rebuild produced invalid G2/M root scores.")
  }

  phenotype_data <- metadata[
    ,
    c("sample", "Malignant_subcluster", "G2M_root_score"),
    drop = FALSE
  ]

  phenotype_data$sample <- factor(as.character(phenotype_data$sample))
  phenotype_data$Malignant_subcluster <- factor(
    as.character(phenotype_data$Malignant_subcluster),
    levels = FIXED_CLUSTER_ORDER
  )
  rownames(phenotype_data) <- colnames(counts)

  if (!identical(rownames(phenotype_data), colnames(counts))) {
    stop("Final rebuilt phenotype metadata is not aligned with counts.")
  }

  saveRDS(
    list(
      schema_version = "LSCC_Monocle2_real_scRNA_input_auto_rebuilt_v1",
      count_mat = as(counts, "dgCMatrix"),
      pd_df = phenotype_data,
      ordering_genes = ordering_genes,
      target_gene = TARGET_GENE,
      root_column = "G2M_root_score",
      expected_subclusters = FIXED_CLUSTER_ORDER,
      max_cells_per_sample_subcluster = MAX_CELLS_PER_SAMPLE_SUBCLUSTER,
      root_rule = "Root = Monocle state with the lowest median G2M_root_score.",
      source_highcnv_object_rds = highcnv_metadata_rds,
      source_highcnv_de_object_rds = highcnv_de_rds,
      source_marker_csv = marker_csv,
      raw_count_assay = assay_name
    ),
    INPUT_RDS
  )

  writeLines(
    c(
      "SUCCESS: Missing Monocle Part-1 input was rebuilt automatically.",
      paste0("Prepared input: ", normalizePath(INPUT_RDS, winslash = "/", mustWork = TRUE)),
      paste0("High-CNV metadata source: ", highcnv_metadata_rds),
      paste0("High-CNV DE/count source: ", highcnv_de_rds),
      paste0("Marker CSV source: ", marker_csv),
      paste0("Raw count assay: ", assay_name),
      paste0("Cells retained after balancing: ", ncol(counts)),
      paste0("Ordering genes retained: ", length(ordering_genes)),
      "The same script will now continue into the complete Part-2 DDRTree run."
    ),
    AUTO_REBUILD_STATUS_TXT
  )

  rm(
    highcnv_obj, highcnv_de_obj, counts_all, counts, metadata,
    phenotype_data, markers, marker_data
  )
  gc()

  cat(
    "\nAUTO-REBUILD SUCCESSFUL\n",
    "Prepared input:\n",
    INPUT_RDS,
    "\n\nProceeding automatically to the complete Part-2 trajectory run.\n",
    sep = ""
  )

  invisible(INPUT_RDS)
}

original_run_trajectory_and_figures <- run_trajectory_and_figures

run_trajectory_and_figures <- function() {
  prepared_input_found <- find_prepared_input_rds()

  if (file.exists(CDS_RDS)) {
    return(refresh_figures_from_completed_cds())
  }

  if (!is.na(prepared_input_found)) {
    INPUT_RDS <<- prepared_input_found

    message(
      "\nAUTO-DISCOVERY: prepared Monocle input found at:\n",
      INPUT_RDS,
      "\nThe complete Part-2 trajectory run will now proceed.\n"
    )

    return(original_run_trajectory_and_figures())
  }

  message(
    "\nAUTO-REBUILD: neither prior Monocle RDS file was found.\n",
    "Recreating Monocle2_full_analysis_input.rds directly from the real ",
    "High-CNV source objects, then continuing automatically.\n"
  )

  rebuild_prepared_input_from_real_sources()

  if (!file.exists(INPUT_RDS)) {
    stop(
      "Auto-rebuild finished without creating the required input RDS:\n",
      INPUT_RDS
    )
  }

  original_run_trajectory_and_figures()
}

# EXECUTE SELF-HEALING PART 2

FIGURE3_COMBINED_W <- 17.40
FIGURE3A_COMBINED_H <- 6.60
FIGURE3B_COMBINED_H <- 4.45
FIGURE8_COMBINED_W <- 17.40
FIGURE8_COMBINED_H <- 11.20

extract_ddrtree_backbone_segments <- function(cds_object, trajectory_df) {
  empty_segments <- data.frame(
    x = numeric(0), y = numeric(0),
    xend = numeric(0), yend = numeric(0),
    Source = character(0),
    stringsAsFactors = FALSE
  )

  as_principal_coordinate_table <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NULL)

    x <- tryCatch(as.matrix(x), error = function(e) NULL)
    if (is.null(x) || nrow(x) < 2L || ncol(x) < 2L) return(NULL)

    if (nrow(x) == 2L && ncol(x) > 2L) x <- t(x)

    if (ncol(x) < 2L || nrow(x) < 2L) return(NULL)

    node_names <- rownames(x)
    if (is.null(node_names) || length(node_names) != nrow(x)) {
      node_names <- as.character(seq_len(nrow(x)))
    }

    out <- data.frame(
      Node = as.character(node_names),
      x = as.numeric(x[, 1]),
      y = as.numeric(x[, 2]),
      stringsAsFactors = FALSE
    )

    out <- out[
      is.finite(out$x) & is.finite(out$y) & !duplicated(out$Node),
      , drop = FALSE
    ]

    if (nrow(out) < 2L) return(NULL)
    out
  }

  graph_to_segments <- function(graph, coordinates, source_name) {
    if (is.null(graph) || !inherits(graph, "igraph") || igraph::ecount(graph) < 1L) {
      return(empty_segments)
    }

    edges <- tryCatch(igraph::as_data_frame(graph, what = "edges"), error = function(e) NULL)
    if (is.null(edges) || !all(c("from", "to") %in% colnames(edges))) {
      return(empty_segments)
    }

    coordinates <- coordinates[
      is.finite(coordinates$x) & is.finite(coordinates$y),
      , drop = FALSE
    ]

    if (nrow(coordinates) < 2L) return(empty_segments)

    from_index <- match(as.character(edges$from), coordinates$Node)
    to_index <- match(as.character(edges$to), coordinates$Node)

    if (all(is.na(from_index)) && all(is.na(to_index)) &&
        igraph::vcount(graph) == nrow(coordinates)) {
      vertex_names <- igraph::V(graph)$name
      if (is.null(vertex_names) || length(vertex_names) != nrow(coordinates)) {
        vertex_names <- as.character(seq_len(nrow(coordinates)))
      }
      coordinates$Node <- as.character(vertex_names)
      from_index <- match(as.character(edges$from), coordinates$Node)
      to_index <- match(as.character(edges$to), coordinates$Node)
    }

    segments <- data.frame(
      x = coordinates$x[from_index],
      y = coordinates$y[from_index],
      xend = coordinates$x[to_index],
      yend = coordinates$y[to_index],
      Source = source_name,
      stringsAsFactors = FALSE
    )

    segments <- segments[
      is.finite(segments$x) & is.finite(segments$y) &
        is.finite(segments$xend) & is.finite(segments$yend),
      , drop = FALSE
    ]

    if (nrow(segments) > 0L) {
      key <- paste(
        pmin(segments$x, segments$xend),
        pmin(segments$y, segments$yend),
        pmax(segments$x, segments$xend),
        pmax(segments$y, segments$yend),
        sep = "|"
      )
      segments <- segments[!duplicated(key), , drop = FALSE]
    }

    segments
  }

  principal_df <- tryCatch(
    as_principal_coordinate_table(monocle::reducedDimK(cds_object)),
    error = function(e) NULL
  )

  if (is.null(principal_df)) {
    return(list(
      segments = empty_segments,
      status = "No compact DDRTree principal-node coordinates were available; backbone intentionally omitted rather than drawing cell-to-cell edges."
    ))
  }

  if (nrow(principal_df) >= max(100L, floor(0.50 * nrow(trajectory_df)))) {
    return(list(
      segments = empty_segments,
      status = paste0(
        "The available reducedDimK matrix contains ", nrow(principal_df),
        " nodes, which is too close to the cell count to be a compact principal-node skeleton; ",
        "backbone intentionally omitted rather than drawing cell-to-cell edges."
      )
    ))
  }

  stored_tree <- tryCatch(monocle::minSpanningTree(cds_object), error = function(e) NULL)

  if (!is.null(stored_tree) && inherits(stored_tree, "igraph") &&
      igraph::vcount(stored_tree) == nrow(principal_df)) {
    exact_segments <- graph_to_segments(
      stored_tree,
      principal_df,
      source_name = "stored_principal_node_MST"
    )

    if (nrow(exact_segments) > 0L) {
      return(list(
        segments = exact_segments,
        status = paste0(
          "DDRTree backbone drawn from ", nrow(exact_segments),
          " compact stored principal-node tree edges."
        )
      ))
    }
  }

  distances <- as.matrix(stats::dist(principal_df[, c("x", "y"), drop = FALSE]))
  compact_graph <- igraph::graph_from_adjacency_matrix(
    distances,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  compact_tree <- igraph::mst(
    compact_graph,
    weights = igraph::E(compact_graph)$weight
  )

  igraph::V(compact_tree)$name <- principal_df$Node

  reconstructed_segments <- graph_to_segments(
    compact_tree,
    principal_df,
    source_name = "reconstructed_principal_node_MST"
  )

  if (nrow(reconstructed_segments) == 0L) {
    return(list(
      segments = empty_segments,
      status = "A compact principal-node backbone could not be assembled; backbone intentionally omitted."
    ))
  }

  list(
    segments = reconstructed_segments,
    status = paste0(
      "DDRTree backbone reconstructed from ", nrow(principal_df),
      " compact principal nodes and ", nrow(reconstructed_segments),
      " branch-level MST edges; no cell-to-cell edges were drawn."
    )
  )
}

create_all_figures <- function(cds, target_gene, root_state) {
  trajectory_df <- make_trajectory_table(cds)

  trajectory_df <- trajectory_df[
    is.finite(trajectory_df$DDRTree_1) &
      is.finite(trajectory_df$DDRTree_2) &
      is.finite(trajectory_df$Pseudotime) &
      !is.na(trajectory_df$Malignant_subcluster),
    ,
    drop = FALSE
  ]

  if (nrow(trajectory_df) < 100L) {
    stop("Too few valid cells remained for trajectory plotting.")
  }

  write.csv(trajectory_df, TRAJECTORY_CSV, row.names = FALSE)

  trajectory_x_limits <- range(trajectory_df$DDRTree_1, na.rm = TRUE)
  trajectory_y_limits <- range(trajectory_df$DDRTree_2, na.rm = TRUE)
  pseudotime_limits <- range(trajectory_df$Pseudotime, na.rm = TRUE)

  backbone_result <- extract_ddrtree_backbone_segments(cds, trajectory_df)
  backbone_segments <- backbone_result$segments

  write.csv(
    backbone_segments,
    BACKBONE_SEGMENTS_CSV,
    row.names = FALSE
  )

  writeLines(backbone_result$status, BACKBONE_STATUS_TXT)

  backbone_layer <- if (nrow(backbone_segments) > 0L) {
    ggplot2::geom_segment(
      data = backbone_segments,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      colour = "grey20",
      linewidth = 0.48,
      alpha = 0.80,
      lineend = "round"
    )
  } else {
    NULL
  }

  root_state_label <- paste0("State ", as.integer(root_state))
  root_state_cells <- trajectory_df[
    as.character(trajectory_df$Monocle_state) == root_state_label,
    ,
    drop = FALSE
  ]

  if (nrow(root_state_cells) < 1L) {
    stop(
      "The selected root state ",
      root_state_label,
      " was not found in the plotted trajectory."
    )
  }

  root_point <- root_state_cells[
    which.min(root_state_cells$Pseudotime),
    ,
    drop = FALSE
  ]

  x_span <- diff(trajectory_x_limits)
  y_span <- diff(trajectory_y_limits)
  x_mid <- mean(trajectory_x_limits)
  y_mid <- mean(trajectory_y_limits)

  x_direction <- if (root_point$DDRTree_1 >= x_mid) -1 else 1
  y_direction <- if (root_point$DDRTree_2 >= y_mid) -1 else 1

  root_label_x <- root_point$DDRTree_1 + x_direction * max(0.08 * x_span, 0.32)
  root_label_y <- root_point$DDRTree_2 + y_direction * max(0.07 * y_span, 0.28)

  root_label_x <- min(
    max(root_label_x, trajectory_x_limits[1] + 0.09 * x_span),
    trajectory_x_limits[2] - 0.09 * x_span
  )

  root_label_y <- min(
    max(root_label_y, trajectory_y_limits[1] + 0.09 * y_span),
    trajectory_y_limits[2] - 0.09 * y_span
  )

  root_label_data <- data.frame(
    x = root_label_x,
    y = root_label_y,
    label = paste0("Root\n(", root_state_label, ")"),
    stringsAsFactors = FALSE
  )

  root_connector_layer <- ggplot2::geom_segment(
    data = data.frame(
      x = root_point$DDRTree_1,
      y = root_point$DDRTree_2,
      xend = root_label_data$x,
      yend = root_label_data$y
    ),
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.34,
    alpha = 0.70,
    lineend = "round"
  )

  root_marker_layer <- ggplot2::geom_point(
    data = root_point,
    ggplot2::aes(x = DDRTree_1, y = DDRTree_2),
    inherit.aes = FALSE,
    shape = 8,
    size = 3.35,
    stroke = 0.90,
    colour = "black"
  )

  root_label_layer <- ggplot2::geom_label(
    data = root_label_data,
    ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    family = FONT_FAMILY,
    size = 2.75,
    fontface = "bold",
    label.size = 0.32,
    label.padding = grid::unit(0.12, "lines"),
    fill = "white",
    colour = "black",
    alpha = 0.98
  )

  common_trajectory_theme <- theme_manuscript(
    show_grid = FALSE,
    legend_position = "right"
  ) +
    ggplot2::theme(
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.box.just = "top",
      legend.justification = c(0, 1),
      legend.spacing.y = grid::unit(3, "pt"),
      legend.box.spacing = grid::unit(4, "pt"),
      legend.key.height = grid::unit(10, "pt"),
      legend.key.width = grid::unit(10, "pt")
    )

  p_pseudotime <- ggplot2::ggplot(
    trajectory_df,
    ggplot2::aes(
      x = DDRTree_1,
      y = DDRTree_2,
      colour = Pseudotime
    )
  ) +
    backbone_layer +
    ggplot2::geom_point(size = 0.80, alpha = 0.88) +
    root_connector_layer +
    root_marker_layer +
    root_label_layer +
    ggplot2::scale_colour_gradientn(
      colours = PSEUDOTIME_COLOURS,
      limits = pseudotime_limits,
      name = "Pseudotime"
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_colorbar(
        barwidth = grid::unit(0.38, "cm"),
        barheight = grid::unit(2.8, "cm"),
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = trajectory_x_limits,
      ylim = trajectory_y_limits,
      expand = FALSE
    ) +
    ggplot2::labs(
      x = "DDRTree component 1",
      y = "DDRTree component 2"
    ) +
    common_trajectory_theme

  state_levels <- levels(trajectory_df$Monocle_state)

  if (length(state_levels) <= length(MONOCLE_STATE_COLOURS)) {
    state_colours <- MONOCLE_STATE_COLOURS[seq_along(state_levels)]
  } else {
    state_colours <- grDevices::hcl.colors(
      length(state_levels),
      palette = "Dynamic"
    )
  }
  names(state_colours) <- state_levels

  p_state <- ggplot2::ggplot(
    trajectory_df,
    ggplot2::aes(
      x = DDRTree_1,
      y = DDRTree_2,
      colour = Monocle_state
    )
  ) +
    backbone_layer +
    ggplot2::geom_point(size = 0.80, alpha = 0.88) +
    ggplot2::scale_colour_manual(
      values = state_colours,
      name = "States",
      drop = FALSE
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol = 1,
        byrow = TRUE,
        title.position = "top",
        override.aes = list(size = 2.1, alpha = 1)
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = trajectory_x_limits,
      ylim = trajectory_y_limits,
      expand = FALSE
    ) +
    ggplot2::labs(
      x = "DDRTree component 1",
      y = "DDRTree component 2"
    ) +
    common_trajectory_theme

  background_cells <- trajectory_df[
    !(as.character(trajectory_df$Malignant_subcluster) %in% IMPORTANT_CLUSTERS),
    ,
    drop = FALSE
  ]

  selected_cells <- trajectory_df[
    as.character(trajectory_df$Malignant_subcluster) %in% IMPORTANT_CLUSTERS,
    ,
    drop = FALSE
  ]

  p_selected <- ggplot2::ggplot() +
    backbone_layer +
    ggplot2::geom_point(
      data = background_cells,
      ggplot2::aes(x = DDRTree_1, y = DDRTree_2),
      colour = COL_NS,
      size = 0.68,
      alpha = 0.35
    ) +
    ggplot2::geom_point(
      data = selected_cells,
      ggplot2::aes(
        x = DDRTree_1,
        y = DDRTree_2,
        colour = Malignant_subcluster
      ),
      size = 0.92,
      alpha = 0.95
    ) +
    ggplot2::scale_colour_manual(
      values = IMPORTANT_CLUSTER_COLOURS,
      limits = IMPORTANT_CLUSTERS,
      name = "Subclusters",
      drop = FALSE
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol = 1,
        byrow = TRUE,
        title.position = "top",
        override.aes = list(size = 2.1, alpha = 1)
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = trajectory_x_limits,
      ylim = trajectory_y_limits,
      expand = FALSE
    ) +
    ggplot2::labs(
      x = "DDRTree component 1",
      y = "DDRTree component 2"
    ) +
    common_trajectory_theme

  top_trajectory_panels <- patchwork::wrap_plots(
    p_pseudotime,
    p_state,
    p_selected,
    ncol = 3
  )

  figure_3a <- top_trajectory_panels +
    patchwork::plot_annotation(tag_levels = "a")

  figure_3a <- figure_3a &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PANEL_TAG_PT,
        colour = "black"
      ),
      plot.tag.position = c(0.012, 0.988)
    )

  save_plot_all_formats(
    plot_obj = figure_3a,
    filename_stem = FIGURE_3A_STEM,
    dir_path = PSEUDOTIME_MONOCLE_DIR,
    width = FIGURE3_COMBINED_W,
    height = FIGURE3A_COMBINED_H
  )

  expression_matrix <- Biobase::exprs(cds)

  if (!(target_gene %in% rownames(expression_matrix))) {
    stop(target_gene, " is absent from the Monocle expression matrix.")
  }

  size_factors <- BiocGenerics::sizeFactors(cds)

  if (is.null(size_factors)) {
    stop("Monocle size factors are missing from the CellDataSet.")
  }

  if (is.null(names(size_factors))) {
    names(size_factors) <- colnames(expression_matrix)
  }

  valid_cells <- trajectory_df$Cell[
    trajectory_df$Cell %in% colnames(expression_matrix) &
      is.finite(trajectory_df$Pseudotime)
  ]

  if (length(valid_cells) < 20L) {
    stop("Too few cells have finite pseudotime values.")
  }

  valid_cells <- valid_cells[
    order(
      trajectory_df$Pseudotime[
        match(valid_cells, trajectory_df$Cell)
      ]
    )
  ]

  valid_size_factors <- size_factors[valid_cells]

  if (
    any(!is.finite(valid_size_factors)) ||
      any(valid_size_factors <= 0)
  ) {
    stop("Invalid Monocle size factors were detected.")
  }

  mybl2_expression <- log1p(
    as.numeric(expression_matrix[target_gene, valid_cells]) /
      pmax(valid_size_factors, .Machine$double.eps)
  )

  mybl2_df <- data.frame(
    Cell = valid_cells,
    Pseudotime = trajectory_df$Pseudotime[
      match(valid_cells, trajectory_df$Cell)
    ],
    Malignant_subcluster = as.character(
      trajectory_df$Malignant_subcluster[
        match(valid_cells, trajectory_df$Cell)
      ]
    ),
    log1p_size_factor_normalised_MYBL2 = mybl2_expression,
    stringsAsFactors = FALSE
  )

  write.csv(mybl2_df, MYBL2_CELL_CSV, row.names = FALSE)

  mybl2_selected <- mybl2_df[
    mybl2_df$Malignant_subcluster %in% IMPORTANT_CLUSTERS &
      is.finite(mybl2_df$Pseudotime) &
      is.finite(mybl2_df$log1p_size_factor_normalised_MYBL2),
    ,
    drop = FALSE
  ]

  mybl2_selected$Malignant_subcluster <- factor(
    mybl2_selected$Malignant_subcluster,
    levels = IMPORTANT_CLUSTERS
  )

  if (
    !setequal(
      unique(as.character(mybl2_selected$Malignant_subcluster)),
      IMPORTANT_CLUSTERS
    )
  ) {
    stop("At least one selected subtype is missing from MYBL2 pseudotime data.")
  }

  loess_curve <- do.call(
    rbind,
    lapply(IMPORTANT_CLUSTERS, function(cluster_name) {
      make_loess_curve(
        data = mybl2_selected,
        cluster_name = cluster_name,
        span_value = LOESS_SPAN[[cluster_name]],
        ci_level = LOESS_CI_LEVEL
      )
    })
  )

  loess_curve$Malignant_subcluster <- factor(
    loess_curve$Malignant_subcluster,
    levels = IMPORTANT_CLUSTERS
  )

  write.csv(loess_curve, LOESS_CURVE_CSV, row.names = FALSE)

  supported_upper_ci <- loess_curve$LOESS_upper_CI[
    !is.na(loess_curve$Ribbon_supported) &
      loess_curve$Ribbon_supported
  ]

  mybl2_y_max <- max(
    mybl2_selected$log1p_size_factor_normalised_MYBL2,
    supported_upper_ci,
    na.rm = TRUE
  )

  if (!is.finite(mybl2_y_max) || mybl2_y_max <= 0) {
    mybl2_y_max <- max(
      mybl2_selected$log1p_size_factor_normalised_MYBL2,
      na.rm = TRUE
    )
  }

  if (!is.finite(mybl2_y_max) || mybl2_y_max <= 0) {
    mybl2_y_max <- 1
  }

  make_mybl2_panel <- function(cluster_name) {
    point_data <- mybl2_selected[
      as.character(mybl2_selected$Malignant_subcluster) == cluster_name,
      ,
      drop = FALSE
    ]

    line_data <- loess_curve[
      as.character(loess_curve$Malignant_subcluster) == cluster_name,
      ,
      drop = FALSE
    ]

    ribbon_data <- line_data[
      !is.na(line_data$Ribbon_supported) &
        line_data$Ribbon_supported &
        is.finite(line_data$LOESS_lower_CI) &
        is.finite(line_data$LOESS_upper_CI),
      ,
      drop = FALSE
    ]

    cluster_colour <- IMPORTANT_CLUSTER_COLOURS[[cluster_name]]
    panel_title <- paste0(
      cluster_name,
      " (n = ",
      nrow(point_data),
      ")"
    )

    panel_x_range <- range(point_data$Pseudotime, na.rm = TRUE)
    panel_x_span <- diff(panel_x_range)

    panel_x_padding <- if (is.finite(panel_x_span) && panel_x_span > 0) {
      max(0.025 * panel_x_span, 0.06)
    } else {
      0.50
    }

    panel_x_limits <- c(
      panel_x_range[1] - panel_x_padding,
      panel_x_range[2] + panel_x_padding
    )

    ggplot2::ggplot(
      point_data,
      ggplot2::aes(
        x = Pseudotime,
        y = log1p_size_factor_normalised_MYBL2
      )
    ) +
      ggplot2::geom_point(
        colour = cluster_colour,
        alpha = 0.20,
        size = 0.78
      ) +
      ggplot2::geom_ribbon(
        data = ribbon_data,
        ggplot2::aes(
          x = Pseudotime,
          ymin = LOESS_lower_CI,
          ymax = LOESS_upper_CI,
          group = Ribbon_segment
        ),
        inherit.aes = FALSE,
        fill = cluster_colour,
        alpha = 0.16,
        na.rm = TRUE
      ) +
      ggplot2::geom_line(
        data = line_data,
        ggplot2::aes(
          x = Pseudotime,
          y = LOESS_smoothed_log1p_expression
        ),
        inherit.aes = FALSE,
        colour = cluster_colour,
        linewidth = max(GEOM_LWD, 0.85),
        na.rm = TRUE
      ) +
      ggplot2::coord_cartesian(
        xlim = panel_x_limits,
        ylim = c(0, mybl2_y_max * 1.05),
        expand = FALSE
      ) +
      ggplot2::labs(
        title = panel_title,
        x = "Inferred pseudotime",
        y = "log1p(MYBL2 expression)"
      ) +
      theme_manuscript(
        show_grid = FALSE,
        legend_position = "none"
      )
  }

  p_cluster_0 <- make_mybl2_panel("Cluster 0")
  p_cluster_8 <- make_mybl2_panel("Cluster 8")

  figure_3b <- (
    p_cluster_0 |
      p_cluster_8
  ) +
    patchwork::plot_annotation(tag_levels = "a")

  figure_3b <- figure_3b &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PANEL_TAG_PT,
        colour = "black"
      ),
      plot.tag.position = c(0.012, 0.988)
    )

  save_plot_all_formats(
    plot_obj = figure_3b,
    filename_stem = FIGURE_3B_STEM,
    dir_path = PSEUDOTIME_MONOCLE_DIR,
    width = FIGURE3_COMBINED_W,
    height = FIGURE3B_COMBINED_H
  )

  bottom_mybl2_panels <- (
    p_cluster_0 |
      p_cluster_8
  )

  figure_8 <- (
    top_trajectory_panels /
      bottom_mybl2_panels
  ) +
    patchwork::plot_layout(heights = c(1.00, 0.82)) +
    patchwork::plot_annotation(tag_levels = "a")

  figure_8 <- figure_8 &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PANEL_TAG_PT,
        colour = "black"
      ),
      plot.tag.position = c(0.012, 0.988)
    )

  save_plot_all_formats(
    plot_obj = figure_8,
    filename_stem = FIGURE_8_STEM,
    dir_path = PSEUDOTIME_MONOCLE_DIR,
    width = FIGURE8_COMBINED_W,
    height = FIGURE8_COMBINED_H
  )

  invisible(
    list(
      Figure_3A = figure_3a,
      Figure_3B = figure_3b,
      Figure_8 = figure_8
    )
  )
}

ORDERCELLS_DIAGNOSTIC_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_orderCells_DDRTree_diagnostics.txt"
)

replace_existing_monocle_function <- function(binding_name, value) {
  ns <- asNamespace("monocle")

  if (!exists(binding_name, envir = ns, inherits = FALSE)) {
    stop("Cannot install the orderCells patch: missing Monocle function ", binding_name)
  }

  was_locked <- bindingIsLocked(binding_name, ns)
  if (was_locked) unlockBinding(binding_name, ns)

  on.exit({
    if (was_locked && !bindingIsLocked(binding_name, ns)) lockBinding(binding_name, ns)
  }, add = TRUE)

  assign(binding_name, value, envir = ns)
  invisible(TRUE)
}

write_ordercells_geometry_diagnostics <- function(cds, stage_label) {
  z <- tryCatch(as.matrix(monocle::reducedDimS(cds)), error = function(e) NULL)
  y <- tryCatch(as.matrix(monocle::reducedDimK(cds)), error = function(e) NULL)
  graph <- tryCatch(monocle::minSpanningTree(cds), error = function(e) NULL)

  report <- c(
    "Monocle 2 DDRTree / orderCells dimension diagnostics",
    paste0("Stage: ", stage_label),
    paste0("Cells: ", ncol(cds)),
    paste0("reducedDimS dimensions: ", if (is.null(z)) "unavailable" else paste(dim(z), collapse = " x ")),
    paste0("reducedDimK dimensions: ", if (is.null(y)) "unavailable" else paste(dim(y), collapse = " x ")),
    paste0("MST vertices: ", if (is.null(graph)) "unavailable" else igraph::vcount(graph)),
    paste0("MST edges: ", if (is.null(graph)) "unavailable" else igraph::ecount(graph)),
    paste0("dim_reduce_type: ", cds@dim_reduce_type),
    "Patch: dimension-safe project2MST() active."
  )

  writeLines(report, ORDERCELLS_DIAGNOSTIC_TXT)
  invisible(report)
}

install_dimension_safe_project2mst_patch <- function() {
  monocle_ns <- asNamespace("monocle")

  project2MST_dimension_safe <- function(cds, Projection_Method) {
    dp_mst <- minSpanningTree(cds)
    Z <- as.matrix(reducedDimS(cds))
    Y <- as.matrix(reducedDimK(cds))
    n_cells <- ncol(cds)

    if (ncol(Z) != n_cells && nrow(Z) == n_cells) Z <- t(Z)

    if (ncol(Z) != n_cells || nrow(Z) < 2L) {
      stop(
        "Dimension-safe project2MST: reducedDimS must be components x cells; observed ",
        paste(dim(Z), collapse = " x "), ", expected second dimension ", n_cells, "."
      )
    }

    if (nrow(Y) != nrow(Z) && ncol(Y) == nrow(Z)) Y <- t(Y)

    if (nrow(Y) != nrow(Z) || ncol(Y) < 2L) {
      stop(
        "Dimension-safe project2MST: reducedDimK must be components x principal nodes; observed ",
        paste(dim(Y), collapse = " x "), "."
      )
    }

    cell_names <- colnames(Z)
    if (is.null(cell_names) || length(cell_names) != n_cells) cell_names <- colnames(cds)
    if (is.null(cell_names) || length(cell_names) != n_cells) {
      cell_names <- as.character(seq_len(n_cells))
    }
    colnames(Z) <- cell_names

    principal_names <- colnames(Y)
    if (is.null(principal_names) || length(principal_names) != ncol(Y)) {
      principal_names <- paste0("Y_", seq_len(ncol(Y)))
    }
    colnames(Y) <- principal_names

    if (is.null(dp_mst) || !inherits(dp_mst, "igraph") ||
        igraph::vcount(dp_mst) != ncol(Y)) {
      principal_distance <- as.matrix(stats::dist(t(Y)))
      dp_graph <- igraph::graph_from_adjacency_matrix(
        principal_distance,
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE
      )
      dp_mst <- igraph::mst(dp_graph, weights = igraph::E(dp_graph)$weight)
    }
    igraph::V(dp_mst)$name <- principal_names

    z_norm <- colSums(Z * Z)
    y_norm <- colSums(Y * Y)
    squared_distances <- outer(y_norm, z_norm, "+") - 2 * crossprod(Y, Z)
    squared_distances[squared_distances < 0 & squared_distances > -1e-10] <- 0
    closest_vertex <- apply(squared_distances, 2L, which.min)

    P <- matrix(
      NA_real_,
      nrow = nrow(Z),
      ncol = n_cells,
      dimnames = list(rownames(Z), cell_names)
    )

    tip_leaves <- igraph::as_ids(
      igraph::V(dp_mst)[igraph::degree(dp_mst) == 1L]
    )

    project_to_segment <- function(point, segment) {
      a <- as.numeric(segment[, 1])
      b <- as.numeric(segment[, 2])
      point <- as.numeric(point)
      ab <- b - a
      denom <- sum(ab * ab)
      if (!is.finite(denom) || denom <= .Machine$double.eps) return(a)
      t_value <- sum((point - a) * ab) / denom
      t_value <- max(0, min(1, t_value))
      a + t_value * ab
    }

    for (i in seq_len(n_cells)) {
      nearest_index <- closest_vertex[i]
      nearest_name <- principal_names[nearest_index]
      neighbor_names <- igraph::as_ids(
        igraph::neighbors(dp_mst, v = nearest_name, mode = "all")
      )

      if (length(neighbor_names) == 0L) {
        P[, i] <- Y[, nearest_index]
        next
      }

      point <- Z[, i]
      candidates <- lapply(neighbor_names, function(neighbor_name) {
        neighbor_index <- match(neighbor_name, principal_names)
        if (is.na(neighbor_index)) return(NULL)
        segment <- Y[, c(nearest_index, neighbor_index), drop = FALSE]

        projected <- tryCatch({
          if (is.function(Projection_Method) && !(nearest_name %in% tip_leaves)) {
            as.numeric(Projection_Method(point, segment))
          } else {
            project_to_segment(point, segment)
          }
        }, error = function(e) project_to_segment(point, segment))

        if (length(projected) != nrow(Z) || any(!is.finite(projected))) {
          projected <- project_to_segment(point, segment)
        }
        projected
      })
      candidates <- Filter(Negate(is.null), candidates)

      if (length(candidates) == 0L) {
        P[, i] <- Y[, nearest_index]
      } else {
        candidate_matrix <- do.call(cbind, candidates)
        candidate_sq_dist <- colSums((candidate_matrix - point)^2)
        P[, i] <- candidate_matrix[, which.min(candidate_sq_dist)]
      }
    }

    cell_distance <- as.matrix(stats::dist(t(P)))
    if (nrow(cell_distance) != n_cells || ncol(cell_distance) != n_cells) {
      stop("Dimension-safe project2MST produced an invalid projected-cell distance matrix.")
    }

    nonzero_distance <- cell_distance[cell_distance > 0 & is.finite(cell_distance)]
    minimum_distance <- if (length(nonzero_distance) > 0L) {
      min(nonzero_distance)
    } else {
      1e-8
    }

    cell_distance <- cell_distance + minimum_distance
    diag(cell_distance) <- 0
    rownames(cell_distance) <- cell_names
    colnames(cell_distance) <- cell_names

    projected_graph <- igraph::graph_from_adjacency_matrix(
      cell_distance,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )
    igraph::V(projected_graph)$name <- cell_names
    projected_tree <- igraph::mst(
      projected_graph,
      weights = igraph::E(projected_graph)$weight
    )

    cellPairwiseDistances(cds) <- cell_distance
    cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- projected_tree
    cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P
    closest_matrix <- matrix(
      as.integer(closest_vertex),
      ncol = 1L,
      dimnames = list(cell_names, "closest_principal_vertex")
    )
    cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_matrix

    cds
  }

  environment(project2MST_dimension_safe) <- monocle_ns
  attr(project2MST_dimension_safe, "LSCC_dimension_safe_project2MST") <- TRUE

  replace_existing_monocle_function("project2MST", project2MST_dimension_safe)

  message(
    "Applied dimension-safe DDRTree project2MST patch for orderCells()."
  )
  invisible(TRUE)
}

load_monocle_environment_before_dimension_safe_patch <- load_monocle_environment
load_monocle_environment <- function() {
  result <- load_monocle_environment_before_dimension_safe_patch()
  install_dimension_safe_project2mst_patch()
  invisible(result)
}

original_run_trajectory_and_figures_before_dimension_safe_patch <- original_run_trajectory_and_figures
original_run_trajectory_and_figures <- function() {
  if (file.exists(ORDERCELLS_DIAGNOSTIC_TXT)) unlink(ORDERCELLS_DIAGNOSTIC_TXT, force = TRUE)
  original_run_trajectory_and_figures_before_dimension_safe_patch()
}

COMPAT_ORDERING_NODE_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_DDRTree_principal_node_assignments.csv"
)
COMPAT_ORDERING_METHOD_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_DDRTree_graph_pseudotime_method.txt"
)

safe_ddrtree_matrix <- function(x, expected_n, role) {
  x <- as.matrix(x)
  if (length(dim(x)) != 2L) {
    stop(role, " is not a two-dimensional matrix.")
  }
  if (ncol(x) != expected_n && nrow(x) == expected_n) {
    x <- t(x)
  }
  if (ncol(x) != expected_n || nrow(x) < 2L) {
    stop(
      role, " has incompatible dimensions ",
      paste(dim(x), collapse = " x "),
      "; expected components x ", expected_n, "."
    )
  }
  x
}

make_structural_principal_states <- function(tree, root_distances) {
  node_names <- igraph::V(tree)$name
  n_nodes <- igraph::vcount(tree)
  degrees <- igraph::degree(tree, mode = "all")
  branch_vertices <- which(degrees >= 3L)

  if (length(branch_vertices) == 0L) {
    return(rep.int(1L, n_nodes))
  }

  pruned <- igraph::delete_vertices(tree, branch_vertices)
  component_membership <- igraph::components(pruned)$membership
  component_names <- names(component_membership)
  state <- rep.int(NA_integer_, n_nodes)
  names(state) <- node_names

  if (length(component_membership) > 0L) {
    state[component_names] <- as.integer(component_membership)
  }

  next_state <- if (all(is.na(state))) 1L else max(state, na.rm = TRUE) + 1L

  for (vertex_index in branch_vertices) {
    vertex_name <- node_names[vertex_index]
    neighbours <- igraph::as_ids(igraph::neighbors(tree, vertex_name, mode = "all"))
    neighbour_states <- state[neighbours]
    usable <- !is.na(neighbour_states)

    if (any(usable)) {
      candidates <- neighbours[usable]
      chosen <- candidates[which.min(root_distances[candidates])]
      state[vertex_name] <- state[chosen]
    } else {
      state[vertex_name] <- next_state
      next_state <- next_state + 1L
    }
  }

  if (any(is.na(state))) {
    state[is.na(state)] <- next_state
  }

  as.integer(state[node_names])
}

run_ddrtree_graph_pseudotime <- function(cds, root_score_column) {
  n_cells <- ncol(cds)
  cell_names <- colnames(cds)
  metadata <- as.data.frame(Biobase::pData(cds), stringsAsFactors = FALSE)

  Z <- safe_ddrtree_matrix(
    monocle::reducedDimS(cds),
    expected_n = n_cells,
    role = "reducedDimS"
  )
  colnames(Z) <- cell_names

  Y <- as.matrix(monocle::reducedDimK(cds))
  if (length(dim(Y)) != 2L) {
    stop("reducedDimK is not a two-dimensional matrix.")
  }
  if (nrow(Y) != nrow(Z) && ncol(Y) == nrow(Z)) {
    Y <- t(Y)
  }
  if (nrow(Y) != nrow(Z) || ncol(Y) < 2L) {
    stop(
      "reducedDimK has incompatible dimensions ",
      paste(dim(Y), collapse = " x "),
      "; expected ", nrow(Z), " x K principal nodes."
    )
  }

  principal_names <- colnames(Y)
  if (is.null(principal_names) || length(principal_names) != ncol(Y)) {
    principal_names <- paste0("Y_", seq_len(ncol(Y)))
  }
  colnames(Y) <- principal_names

  tree <- tryCatch(monocle::minSpanningTree(cds), error = function(e) NULL)
  if (is.null(tree) || !inherits(tree, "igraph") || igraph::vcount(tree) != ncol(Y)) {
    principal_distance <- as.matrix(stats::dist(t(Y)))
    principal_graph <- igraph::graph_from_adjacency_matrix(
      principal_distance,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )
    tree <- igraph::mst(principal_graph, weights = igraph::E(principal_graph)$weight)
  }
  igraph::V(tree)$name <- principal_names

  z_norm <- colSums(Z * Z)
  y_norm <- colSums(Y * Y)
  squared_distance <- outer(y_norm, z_norm, "+") - 2 * crossprod(Y, Z)
  squared_distance[squared_distance < 0 & squared_distance > -1e-10] <- 0
  closest_node_index <- apply(squared_distance, 2L, which.min)
  closest_node_name <- principal_names[closest_node_index]
  cell_to_node_distance <- sqrt(pmax(
    squared_distance[cbind(closest_node_index, seq_len(n_cells))], 0
  ))

  temporary_root <- principal_names[1]
  temporary_distances <- as.numeric(igraph::distances(
    tree, v = temporary_root, to = igraph::V(tree),
    mode = "all", weights = igraph::E(tree)$weight
  ))
  names(temporary_distances) <- principal_names
  principal_state <- make_structural_principal_states(tree, temporary_distances)
  cell_state <- principal_state[closest_node_index]

  root_scores <- suppressWarnings(as.numeric(metadata[[root_score_column]]))
  if (length(root_scores) != n_cells || any(!is.finite(root_scores))) {
    stop("The root-score column is missing or contains non-finite values.")
  }

  state_medians <- tapply(root_scores, cell_state, stats::median, na.rm = TRUE)
  state_medians <- state_medians[is.finite(state_medians)]
  if (length(state_medians) == 0L) {
    stop("No finite median G2/M scores were available to select the root state.")
  }
  root_state <- as.integer(names(state_medians)[which.min(state_medians)])

  node_medians <- tapply(root_scores, closest_node_name, stats::median, na.rm = TRUE)
  node_medians <- node_medians[is.finite(node_medians)]
  candidate_nodes <- names(node_medians)[
    principal_state[match(names(node_medians), principal_names)] == root_state
  ]
  if (length(candidate_nodes) == 0L) {
    candidate_nodes <- names(node_medians)
  }
  root_node <- candidate_nodes[which.min(node_medians[candidate_nodes])]

  graph_distance <- as.numeric(igraph::distances(
    tree, v = root_node, to = igraph::V(tree),
    mode = "all", weights = igraph::E(tree)$weight
  ))
  names(graph_distance) <- principal_names
  if (any(!is.finite(graph_distance))) {
    stop("The DDRTree principal graph is disconnected; graph pseudotime cannot be computed.")
  }

  pseudotime <- graph_distance[closest_node_name] + cell_to_node_distance
  pseudotime <- pseudotime - min(pseudotime, na.rm = TRUE)

  principal_state <- make_structural_principal_states(tree, graph_distance)
  cell_state <- principal_state[closest_node_index]
  root_state <- principal_state[match(root_node, principal_names)]

  node_table <- data.frame(
    Principal_node = principal_names,
    DDRTree_1 = as.numeric(Y[1, ]),
    DDRTree_2 = as.numeric(Y[2, ]),
    Structural_state = principal_state,
    Distance_from_root = as.numeric(graph_distance[principal_names]),
    Median_G2M_root_score = as.numeric(node_medians[principal_names]),
    Is_root_node = principal_names == root_node,
    stringsAsFactors = FALSE
  )
  write.csv(node_table, COMPAT_ORDERING_NODE_CSV, row.names = FALSE)

  pd <- Biobase::pData(cds)
  pd$Pseudotime <- as.numeric(pseudotime)
  pd$State <- as.integer(cell_state)
  Biobase::pData(cds) <- pd

  monocle::minSpanningTree(cds) <- tree

  list(
    cds = cds,
    root_state = as.integer(root_state),
    root_node = root_node,
    state_medians = state_medians,
    node_table = node_table
  )
}

original_run_trajectory_and_figures <- function() {
  cat(
    "\n============================================================\n",
    "STEP 2 / 2: MONOCLE 2 DDRTREE + GRAPH PSEUDOTIME IN RGui 4.4.3\n",
    "============================================================\n",
    sep = ""
  )

  if (!file.exists(INPUT_RDS)) {
    stop("Prepared Monocle input is missing:\n", INPUT_RDS)
  }

  load_monocle_environment()

  old_outputs <- c(
    STATUS_TXT, ERROR_TXT, RUN_LOG, SESSION_TXT, COMPATIBILITY_TXT,
    REDUCTION_TXT, CDS_RDS, TRAJECTORY_CSV, ORDERING_GENES_TRAJECTORY_CSV,
    ROOT_CSV, ROOT_TXT, MYBL2_CELL_CSV, LOESS_CURVE_CSV, BACKBONE_STATUS_TXT,
    BACKBONE_SEGMENTS_CSV, COMPAT_ORDERING_NODE_CSV, COMPAT_ORDERING_METHOD_TXT,
    ORDERCELLS_DIAGNOSTIC_TXT,
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_3A_STEM, ".png")),
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_3A_STEM, ".tiff")),
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_3A_STEM, ".pdf")),
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_3B_STEM, ".png")),
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_3B_STEM, ".tiff")),
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_3B_STEM, ".pdf")),
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_8_STEM, ".png")),
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_8_STEM, ".tiff")),
    file.path(PSEUDOTIME_MONOCLE_DIR, paste0(FIGURE_8_STEM, ".pdf"))
  )
  unlink(old_outputs[file.exists(old_outputs)], force = TRUE)

  writeLines(
    c(
      "DDRTree embedding: monocle 2.34.0 reduceDimension(reduction_method = 'DDRTree').",
      "Pseudotime assignment: direct weighted shortest-path distance on the learned DDRTree principal MST.",
      "Compatibility reason: monocle::orderCells() is bypassed because its legacy internal graph indexing fails under igraph 2.0.3 with 'incorrect number of dimensions'.",
      "Root rule: structural state with the lowest median G2M_root_score; root node is the lowest-median G2M principal node within that state.",
      "This runner does not call monocle::orderCells(), extract_ddrtree_ordering(), or project2MST()."
    ),
    COMPATIBILITY_TXT
  )

  result <- list(
    trajectory_success = FALSE, figure_success = FALSE,
    reduction_fallback_used = FALSE, error = NA_character_,
    figure_error = NA_character_, stage = "started", root_state = NA_integer_,
    root_node = NA_character_, n_cells = 0L, n_ordering_genes = 0L,
    monocle_version = as.character(utils::packageVersion("monocle")),
    igraph_version = as.character(utils::packageVersion("igraph"))
  )

  write_status <- function() {
    writeLines(c(
      "LSCC DDRTree graph-pseudotime status",
      "Run mode: fresh R 4.4.3 / Monocle 2.34.0 / igraph 2.0.3.",
      "DDRTree embedding: completed by Monocle reduceDimension().",
      "Pseudotime engine: compatibility graph-distance implementation; monocle::orderCells() intentionally bypassed.",
      paste0("Final stage: ", result$stage),
      paste0("Trajectory completed: ", result$trajectory_success),
      paste0("Figures completed: ", result$figure_success),
      paste0("Reduction fallback used: ", result$reduction_fallback_used),
      paste0("Cells used: ", result$n_cells),
      paste0("Ordering genes used: ", result$n_ordering_genes),
      paste0("Root state: ", result$root_state),
      paste0("Root principal node: ", result$root_node),
      paste0("Monocle version: ", result$monocle_version),
      paste0("igraph version: ", result$igraph_version),
      paste0("Core error: ", ifelse(is.na(result$error), "none", result$error)),
      paste0("Figure error: ", ifelse(is.na(result$figure_error), "none", result$figure_error))
    ), STATUS_TXT)
    writeLines(capture.output(sessionInfo()), SESSION_TXT)
  }

  cds <- NULL
  target_gene <- NULL

  tryCatch({
    result$stage <- "loading_input"
    log_message("STAGE: loading prepared Monocle input")
    input <- readRDS(INPUT_RDS)
    required_items <- c("count_mat", "pd_df", "ordering_genes", "target_gene", "root_column", "expected_subclusters")
    if (!is.list(input) || !all(required_items %in% names(input))) stop("Prepared Monocle input RDS is incomplete.")

    counts <- as(input$count_mat, "dgCMatrix")
    assert_raw_counts(counts)
    pd <- as.data.frame(input$pd_df, stringsAsFactors = FALSE)
    ordering_genes <- intersect(as.character(input$ordering_genes), rownames(counts))
    target_gene <- as.character(input$target_gene)
    root_column <- as.character(input$root_column)
    expected_subclusters <- as.character(input$expected_subclusters)

    if (!identical(rownames(pd), colnames(counts))) stop("Input metadata and count matrix are not aligned.")
    if (!all(c("sample", "Malignant_subcluster", root_column) %in% names(pd))) stop("Required metadata columns are missing.")
    if (!(target_gene %in% rownames(counts))) stop(target_gene, " is absent from raw counts.")
    if (length(ordering_genes) < MIN_ORDERING_GENES) stop("Too few validated ordering genes: ", length(ordering_genes))

    pd$sample <- factor(as.character(pd$sample))
    pd$Malignant_subcluster <- factor(normalise_cluster(pd$Malignant_subcluster), levels = expected_subclusters)
    if (any(is.na(pd$Malignant_subcluster)) || !setequal(unique(as.character(pd$Malignant_subcluster)), expected_subclusters)) {
      stop("The expected ten High-CNV subclusters are not all present.")
    }
    if (any(!is.finite(pd[[root_column]]))) stop("The G2/M root score contains invalid values.")

    result$n_cells <- ncol(counts)
    result$n_ordering_genes <- length(ordering_genes)

    result$stage <- "creating_CellDataSet"
    log_message("STAGE: creating CellDataSet")
    feature_data <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts), stringsAsFactors = FALSE)
    cds <- monocle::newCellDataSet(
      counts,
      phenoData = Biobase::AnnotatedDataFrame(pd),
      featureData = Biobase::AnnotatedDataFrame(feature_data),
      lowerDetectionLimit = 0.5,
      expressionFamily = VGAM::negbinomial.size()
    )

    result$stage <- "estimating_size_factors"
    log_message("STAGE: estimating size factors")
    cds <- BiocGenerics::estimateSizeFactors(cds)

    result$stage <- "detecting_genes"
    log_message("STAGE: detecting expressed genes")
    cds <- monocle::detectGenes(cds, min_expr = 0.5)

    result$stage <- "compatibility_dispersion"
    log_message("STAGE: creating pooled Monocle-compatible dispersion model")
    cds <- estimate_dispersions_blind_compat(cds, min_cells_detected = 1L, remove_outliers = TRUE, verbose = FALSE)
    cds <- monocle::setOrderingFilter(cds, ordering_genes)

    result$stage <- "DDRTree_reduction"
    log_message("STAGE: running DDRTree with sample adjustment")
    reduction_note <- "Sample-adjusted DDRTree reduction completed."
    cds <- tryCatch(
      monocle::reduceDimension(
        cds, max_components = 2, reduction_method = "DDRTree", norm_method = "log",
        residualModelFormulaStr = "~sample", verbose = FALSE
      ),
      error = function(e) {
        result$reduction_fallback_used <<- TRUE
        reduction_note <<- paste0("Sample-adjusted DDRTree failed; fallback without residual model used: ", conditionMessage(e))
        log_message(reduction_note)
        monocle::reduceDimension(cds, max_components = 2, reduction_method = "DDRTree", norm_method = "log", verbose = FALSE)
      }
    )
    writeLines(reduction_note, REDUCTION_TXT)

    result$stage <- "DDRTRee_graph_pseudotime"
    log_message("STAGE: calculating direct principal-graph pseudotime (orderCells bypassed)")
    ordering_result <- run_ddrtree_graph_pseudotime(cds, root_score_column = root_column)
    cds <- ordering_result$cds
    result$root_state <- ordering_result$root_state
    result$root_node <- ordering_result$root_node

    state_medians <- tapply(
      as.numeric(Biobase::pData(cds)[[root_column]]),
      as.integer(Biobase::pData(cds)$State),
      stats::median, na.rm = TRUE
    )
    state_medians <- state_medians[is.finite(state_medians)]
    write.csv(data.frame(
      State = as.integer(names(state_medians)),
      Median_G2M_root_score = as.numeric(state_medians),
      Is_selected_root = as.integer(names(state_medians)) == result$root_state
    ), ROOT_CSV, row.names = FALSE)

    write.csv(data.frame(Gene = ordering_genes), ORDERING_GENES_TRAJECTORY_CSV, row.names = FALSE)
    writeLines(c(
      paste0("Root state: ", result$root_state),
      paste0("Root principal node: ", result$root_node),
      "Root rule: state with the lowest median G2M_root_score; then lowest-median-G2M principal node inside that state.",
      "Pseudotime: direct weighted shortest-path distance from the selected root node on the DDRTree principal graph plus cell-to-nearest-node distance.",
      "This is an inferred transcriptional continuum, not observed time.",
      reduction_note
    ), ROOT_TXT)

    saveRDS(cds, CDS_RDS)
    result$trajectory_success <- TRUE
    result$stage <- "trajectory_completed"
    log_message("TRAJECTORY SUCCESSFUL | Root state = ", result$root_state, " | Root node = ", result$root_node)
  }, error = function(e) {
    result$error <<- conditionMessage(e)
    result$stage <<- paste0("failed_at_", result$stage)
    writeLines(c("CORE ERROR", paste0("Stage: ", result$stage), paste0("Message: ", conditionMessage(e)), "", "Call stack:", capture.output(sys.calls())), ERROR_TXT)
    log_message("CORE ERROR: ", conditionMessage(e))
  })

  if (isTRUE(result$trajectory_success)) {
    tryCatch({
      result$stage <- "creating_figures"
      log_message("STAGE: creating Figure 3A, Figure 3B and Figure 8")
      create_all_figures(cds = cds, target_gene = target_gene, root_state = result$root_state)
      result$figure_success <- TRUE
      result$stage <- "figures_completed"
      log_message("FIGURES SUCCESSFUL")
    }, error = function(e) {
      result$figure_error <<- conditionMessage(e)
      log_message("FIGURE ERROR: ", conditionMessage(e))
    })
  }

  write_status()
  cat("\n============================================================\n",
      "MONOCLE 2 DDRTREE GRAPH-PSEUDOTIME RUN FINISHED\n",
      "Trajectory completed: ", result$trajectory_success, "\n",
      "Figures completed: ", result$figure_success, "\n",
      "Output folder:\n", PSEUDOTIME_MONOCLE_DIR,
      "\n============================================================\n", sep = "")

  if (!isTRUE(result$trajectory_success) || !isTRUE(result$figure_success)) {
    warning("The run did not fully complete. Read Monocle2_error.txt, Monocle2_status.txt and Monocle2_run.log in Pseudotime_Monocle.")
  }
  invisible(result)
}

# EXECUTE FINAL ROBUST PART 2

run_trajectory_and_figures()
