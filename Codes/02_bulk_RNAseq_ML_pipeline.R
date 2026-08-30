# STEP 01: Setup packages, paths, and helper functions.

if (!exists("project_path")) {
  project_path <- getwd()
}

setwd(project_path)
invisible(gc())
set.seed(123)
options(stringsAsFactors = FALSE)
options(scipen = 100)
options(timeout = 3600)

data_path <- file.path("Data", "Bulk")
results_path <- file.path("Results", "ML")
manuscript_figures_path <- results_path
FIGURE_DIR <- results_path

scrna_figdir <- file.path("Results", "scRNAseq", "figures")

dir.create(data_path, recursive = TRUE, showWarnings = FALSE)
dir.create(results_path, recursive = TRUE, showWarnings = FALSE)
dir.create(scrna_figdir, recursive = TRUE, showWarnings = FALSE)

TARGET_GENE <- "MYBL2"

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

# STEP 02: Prepare bulk expression matrices and split datasets.

step02_dir <- make_stage_dir(2, "Bulk_Preprocessing_Train_Validation_External")

train1_expr  <- file.path(data_path, "GSE127165_raw_counts_GRCh38.p13_NCBI.csv")
train2_expr  <- file.path(data_path, "GSE142083_raw_counts_GRCh38.p13_NCBI.csv")
test_expr    <- file.path(data_path, "GSE130605_raw_counts_GRCh38.p13_NCBI.csv")

train1_annot <- file.path(data_path, "GSE127165_annot.csv")
train2_annot <- file.path(data_path, "GSE142083_annot.csv")
test_annot   <- file.path(data_path, "GSE130605_annot.csv")

pheno_tr1    <- file.path(data_path, "Pheno_Data_GSE127165.csv")
pheno_tr2    <- file.path(data_path, "Pheno_Data_GSE142083.csv")
pheno_te     <- file.path(data_path, "Pheno_Data_GSE130605.csv")

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
  
  if (length(keep_genes) == 0) {
    stop("No genes left after zero filtering. zero_cutoff = ", zero_cutoff)
  }
  
  df_filtered <- df[, c(meta_cols, keep_genes), with = FALSE]
  data.table::fwrite(df_filtered, output_file)
  
  output_file
}

s2_tr1_cpm <- file.path(step02_dir, "tr1_cpm_geneid.csv")
s2_tr2_cpm <- file.path(step02_dir, "tr2_cpm_geneid.csv")
s2_te_cpm  <- file.path(step02_dir, "ext_cpm_geneid.csv")

calc_cpm_tmm_nonnegative(train1_expr, s2_tr1_cpm)
calc_cpm_tmm_nonnegative(train2_expr, s2_tr2_cpm)
calc_cpm_tmm_nonnegative(test_expr,   s2_te_cpm)

s2_tr1_sym <- file.path(step02_dir, "tr1_cpm_symbol.csv")
s2_tr2_sym <- file.path(step02_dir, "tr2_cpm_symbol.csv")
s2_te_sym  <- file.path(step02_dir, "ext_cpm_symbol.csv")

aggregate_geneid_to_symbol(train1_annot, s2_tr1_cpm, s2_tr1_sym)
aggregate_geneid_to_symbol(train2_annot, s2_tr2_cpm, s2_tr2_sym)
aggregate_geneid_to_symbol(test_annot,   s2_te_cpm,  s2_te_sym)

s2_tr1_trp <- file.path(step02_dir, "tr1_cpm_t.csv")
s2_tr2_trp <- file.path(step02_dir, "tr2_cpm_t.csv")
s2_te_trp  <- file.path(step02_dir, "ext_cpm_t.csv")

transpose_symbol_matrix(s2_tr1_sym, s2_tr1_trp)
transpose_symbol_matrix(s2_tr2_sym, s2_tr2_trp)
transpose_symbol_matrix(s2_te_sym,  s2_te_trp)

s2_tr1_meta <- file.path(step02_dir, "tr1_cpm_meta.csv")
s2_tr2_meta <- file.path(step02_dir, "tr2_cpm_meta.csv")
s2_te_meta  <- file.path(step02_dir, "ext_cpm_meta.csv")

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

s2_merged_csv <- file.path(step02_dir, "disc_cpm.csv")
data.table::fwrite(S_tr_merged, s2_merged_csv)

s2_merged_nozero_csv <- file.path(step02_dir, "disc_cpm_nozero.csv")

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

training_file_cpm <- file.path(step02_dir, "train_cpm.csv")
validation_file_cpm <- file.path(step02_dir, "valid_cpm.csv")

write.csv(training_set_cpm, training_file_cpm, row.names = FALSE)
write.csv(validation_set_cpm, validation_file_cpm, row.names = FALSE)

external_meta <- data.table::fread(s2_te_meta)
kept_discovery_genes <- setdiff(colnames(training_set_cpm), c("Sample", "group", "batch"))
external_common_genes <- intersect(kept_discovery_genes, setdiff(names(external_meta), req_cols))

if (length(external_common_genes) < 2) {
  stop("Too few common genes between discovery train set and external GSE130605.")
}

external_set_cpm <- external_meta[, c(req_cols, external_common_genes), with = FALSE]
external_file_cpm <- file.path(step02_dir, "ext_cpm.csv")
data.table::fwrite(external_set_cpm, external_file_cpm)

# STEP 03: Run differential expression analysis and generate Figure 2.

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

train_log2_file <- file.path(step03_dir, "train_log2.csv")
valid_log2_file <- file.path(step03_dir, "valid_log2.csv")
external_log2_file <- file.path(step03_dir, "ext_log2.csv")

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

write.csv(TT_full, file.path(step03_dir, "deg_all.csv"), row.names = FALSE)
write.csv(TT_deg, file.path(step03_dir, "deg_sig.csv"), row.names = FALSE)
writeLines(Bulk_UP, file.path(step03_dir, "deg_up.txt"))
writeLines(Bulk_DOWN, file.path(step03_dir, "deg_down.txt"))

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
  filename_stem = "fig02a_volcano",
  dir_path = FIGURE_DIR,
  width = FIG_SINGLE_W,
  height = FIG_SINGLE_H
)

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
  filename_stem = "fig02b_heatmap",
  dir_path = FIGURE_DIR,
  width = 10.80,
  height = 8.20
)

p_fig2b <- patchwork::wrap_elements(
  full = heatmap_grob_with_tag
)

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
  filename_stem = "fig02_deg",
  dir_path = FIGURE_DIR,
  width = FIGURE2_COMBINED_W,
  height = FIGURE2_COMBINED_H
)

# STEP 04: Run WGCNA and generate Figure 3.

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

write.csv(sft_table, file.path(step04_dir, "wgcna_power.csv"), row.names = FALSE)
writeLines(paste0("Chosen softPower = ", softPower), file.path(step04_dir, "wgcna_power.txt"))

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

write.csv(modTraitCor, file.path(step04_dir, "wgcna_cor.csv"))
write.csv(modTraitPval, file.path(step04_dir, "wgcna_p.csv"))

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

write.csv(mod_table, file.path(step04_dir, paste0(targetModule, "_gs_mm.csv")), row.names = FALSE)
write.csv(data.frame(Gene = WGCNA_Strong_Tumor_Module_Genes), file.path(step04_dir, "wgcna_genes.csv"), row.names = FALSE)

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
  filename_stem = "fig03a_samples",
  dir_path = step04_dir,
  width = 13.50,
  height = 5.80
)

fig3a_file <- file.path(
  step04_dir,
  "fig03a_samples.png"
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
  filename_stem = "fig03b_power",
  dir_path = step04_dir,
  width = 9.60,
  height = 4.70
)

fig3b_file <- file.path(
  step04_dir,
  "fig03b_power.png"
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
  filename_stem = "fig03c_modules",
  dir_path = step04_dir,
  width = 8.20,
  height = 5.80
)

fig3c_file <- file.path(
  step04_dir,
  "fig03c_modules.png"
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
  filename_stem = "fig03d_trait",
  dir_path = step04_dir,
  width = 8.50,
  height = 10.20
)

fig3d_file <- file.path(
  step04_dir,
  "fig03d_trait.png"
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
  filename_stem = "fig03e_mm_gs",
  dir_path = step04_dir,
  width = 6.20,
  height = 5.80
)

fig3e_file <- file.path(
  step04_dir,
  "fig03e_mm_gs.png"
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
  "fig03_wgcna.png"
)

fig3_tiff_file <- file.path(
  step04_dir,
  "fig03_wgcna.tiff"
)

fig3_pdf_file <- file.path(
  step04_dir,
  "fig03_wgcna.pdf"
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
  "fig03_wgcna.png"
)

# STEP 05: Identify overlapping candidate genes.

if (!exists("results_path", inherits = FALSE)) {
  results_path <- file.path("Results", "ML")
}

if (!exists("scrna_figdir", inherits = FALSE)) {
  scrna_figdir <- file.path("Results", "scRNAseq", "figures")
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

bulk_up_file <- file.path(step03_dir, "deg_up.txt")
wgcna_file   <- file.path(step04_dir, "wgcna_genes.csv")

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
  stop("wgcna_genes.csv must contain a column named Gene.")
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
  file.path(step05_dir, "overlap_genes.csv"),
  row.names = FALSE
)

writeLines(Triple_Overlap_Genes, file.path(step05_dir, "overlap_genes.txt"))

if (length(Triple_Overlap_Genes) < 2) {
  stop("Too few triple-overlap genes for LASSO. Check DEG/WGCNA/scRNA thresholds.")
}

# STEP 06: Perform LASSO feature selection.

if (!exists("results_path", inherits = FALSE)) {
  results_path <- file.path("Results", "ML")
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
  train_log2_file <- file.path(results_path, "train_log2.csv")
}

if (!exists("valid_log2_file", inherits = FALSE)) {
  valid_log2_file <- file.path(results_path, "valid_log2.csv")
}

if (!exists("external_log2_file", inherits = FALSE)) {
  external_log2_file <- file.path(
    results_path,
    "ext_log2.csv"
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

triple_file <- file.path(step05_dir, "overlap_genes.txt")
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
  file.path(step06_dir, "lasso_input_genes.csv"),
  row.names = FALSE
)
writeLines(common_model_genes, file.path(step06_dir, "lasso_input_genes.txt"))

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
  file.path(step06_dir, "lasso_coef.csv"),
  row.names = FALSE
)

writeLines(lasso_coef_table$Gene, file.path(step06_dir, "lasso_genes.txt"))

lasso_genes <- lasso_coef_table$Gene

# STEP 07: Perform internal ROC analysis.

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

lasso_genes <- clean_gene_vector(readLines(file.path(step06_dir, "lasso_genes.txt"), warn = FALSE))

common_model_genes <- clean_gene_vector(readLines(file.path(step06_dir, "lasso_input_genes.txt"), warn = FALSE))
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

write.csv(roc_table, file.path(step07_dir, "lasso_auc.csv"), row.names = FALSE)

final_biomarkers <- roc_table %>%
  dplyr::filter(Pass_ROC_Cutoff) %>%
  dplyr::pull(Gene)

if (length(final_biomarkers) < 1) {
  write.csv(
    roc_table,
    file.path(step07_dir, "lasso_auc_none.csv"),
    row.names = FALSE
  )
  stop(
    "No LASSO gene passed the strict ROC cutoff. Required: Train_AUC >= ",
    ROC_AUC_MIN,
    " AND Valid_AUC >= ",
    ROC_AUC_MIN
  )
}

writeLines(final_biomarkers, file.path(step07_dir, "biomarkers.txt"))

write.csv(
  roc_table %>% dplyr::filter(Gene %in% final_biomarkers),
  file.path(step07_dir, "biomarkers_auc.csv"),
  row.names = FALSE
)

# STEP 08: Validate biomarker expression and generate Figure 5.

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

final_biomarkers <- clean_gene_vector(readLines(file.path(step07_dir, "biomarkers.txt"), warn = FALSE))

writeLines(final_biomarkers, file.path(step08_dir, "expr_genes.txt"))

expr_long <- dplyr::bind_rows(
  build_long_df(train_df, "Training", final_biomarkers),
  build_long_df(valid_df, "Validation", final_biomarkers)
)

if (nrow(expr_long) < 1) {
  stop("Expression table is empty. Final biomarker(s) were not found in train/validation matrices.")
}

expr_long$Gene <- factor(expr_long$Gene, levels = final_biomarkers)

write.csv(expr_long, file.path(step08_dir, "expr_long.csv"), row.names = FALSE)

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

options(stringsAsFactors = FALSE, scipen = 100, width = 140)
set.seed(123)

results_path <- file.path("Results", "ML")
FIGURE_DIR <- results_path

scrna_figdir <- file.path("Results", "scRNAseq", "figures")

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
  "deg_up.txt"
)

wgcna_file <- file.path(
  results_path,
  "wgcna_genes.csv"
)

triple_file <- file.path(
  results_path,
  "overlap_genes.txt"
)

lasso_genes_file <- file.path(
  results_path,
  "lasso_genes.txt"
)

final_biomarkers_file <- file.path(
  results_path,
  "biomarkers.txt"
)

train_log2_file <- file.path(
  results_path,
  "train_log2.csv"
)

valid_log2_file <- file.path(
  results_path,
  "valid_log2.csv"
)

external_log2_file <- file.path(
  results_path,
  "ext_log2.csv"
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
  stop("wgcna_genes.csv must contain a column named 'Gene'.")
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
  filename_stem = "fig05a_venn",
  dir_path = FIGURE_DIR,
  width = FIG_DOUBLE_W,
  height = FIG_DOUBLE_H
)

fig5a_file <- file.path(
  FIGURE_DIR,
  "fig05a_venn.png"
)

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
  filename_stem = "fig05b_lasso",
  dir_path = FIGURE_DIR,
  width = FIG_DOUBLE_W,
  height = FIG_DOUBLE_H
)

fig5b_file <- file.path(
  FIGURE_DIR,
  "fig05b_lasso.png"
)

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
  filename_stem = "fig05c_roc",
  dir_path = FIGURE_DIR,
  width = max(FIG_DOUBLE_W, min(3, length(roc_plots)) * FIG_DOUBLE_W),
  height = ceiling(length(roc_plots) / 3) * FIG_DOUBLE_H
)

fig5c_file <- file.path(
  FIGURE_DIR,
  "fig05c_roc.png"
)

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
  filename_stem = "fig05d_expr",
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
  "fig05d_expr.png"
)

rm(list = intersect(
  c("train_df", "valid_df", "roc_train", "roc_valid", "p_fig5c",
    "fig5c_img", "fig5_combined_img"),
  ls()
))

options(stringsAsFactors = FALSE, scipen = 100)

results_path <- file.path("Results", "ML")
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

gene <- "MYBL2"

train_file <- file.path(results_path, "train_log2.csv")
valid_file <- file.path(results_path, "valid_log2.csv")

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
  filename_stem = "fig05c_roc",
  dir_path = FIGURE_DIR,
  width = 7.2,
  height = 6.5
)

fig5c_file <- file.path(
  FIGURE_DIR,
  "fig05c_roc.png"
)

fig5a_file <- file.path(
  FIGURE_DIR,
  "fig05a_venn.png"
)

fig5b_file <- file.path(
  FIGURE_DIR,
  "fig05b_lasso.png"
)

fig5d_file <- file.path(
  FIGURE_DIR,
  "fig05d_expr.png"
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
  filename_stem = "fig05_screening",
  dir_path = FIGURE_DIR
)

# STEP 09: Run external machine-learning validation and generate Figure 6A.

if (!exists("results_path", inherits = FALSE)) {
  results_path <- file.path("Results", "ML")
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
  train_log2_file <- file.path(results_path, "train_log2.csv")
}

if (!exists("external_log2_file", inherits = FALSE)) {
  external_log2_file <- file.path(
    results_path,
    "ext_log2.csv"
  )
}

final_biomarkers_file <- file.path(
  step07_dir,
  "biomarkers.txt"
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

# STEP 09: Run external machine-learning validation and generate Figure 6A.

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
  file.path(step09_dir, "ext_ml_genes.csv"),
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
  file.path(step09_dir, "ext_ml_metrics.csv"),
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
  "fig06a_ext_ml_roc.png"
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

