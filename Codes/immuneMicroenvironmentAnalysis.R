# STEP 01: Setup settings, packages, and paths.

rm(list = ls())
invisible(gc())

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

# STEP 02: Define helper functions.

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

# STEP 03: Prepare pooled tumor cohort.

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
    "tumor_meta.csv"
  ),
  row.names = FALSE
)

write.csv(
  estimate_cpm_sample_by_gene,
  file.path(step12_tabledir, "estimate_input.csv")
)

saveRDS(
  list(
    expression_log2_TMM_CPM_no_batch_correction = expr_sample_by_gene,
    ESTIMATE_input_nonlog_TMM_CPM_no_batch_correction = estimate_cpm_sample_by_gene,
    sample_metadata = sample_metadata,
    global_MYBL2_median = global_mybl2_median,
    input_TMM_CPM_file = discovery_file
  ),
  file.path(step12_rdsdir, "tumor_expr_obj.rds")
)

# STEP 04: Run ssGSEA immune scoring.

official_immune28 <- read_local_charoentong_immune28()
immune28_sets <- official_immune28$gene_sets
immune28_source <- CHAROENTONG_REFERENCE_CITATION
immune28_source_names <- official_immune28$source_names

official_gmt_file <- file.path(
  step12_tabledir,
  "immune28.gmt"
)
write_gmt(immune28_sets, official_gmt_file)

local_workbook_copy <- file.path(
  step12_tabledir,
  "immune28_source.xlsx"
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
    "immune28_source.csv"
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
  file.path(step12_tabledir, "immune28_info.csv"),
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
  file.path(step12_tabledir, "immune28_overlap.csv"),
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
    "\n\nInspect: ", file.path(step12_tabledir, "immune28_overlap.csv")
  )
}

expr_gene_by_sample <- t(expr_sample_by_gene)
ssgsea_scores <- run_ssgsea_safe(expr_gene_by_sample, immune28_sets_used)

ssgsea_scores <- ssgsea_scores[immune28_labels, sample_metadata$Sample, drop = FALSE]

write.csv(
  ssgsea_scores,
  file.path(step12_tabledir, "ssgsea_scores.csv")
)

saveRDS(
  list(
    immune28_sets = immune28_sets,
    immune28_sets_used = immune28_sets_used,
    ssGSEA_scores = ssgsea_scores,
    sample_metadata = sample_metadata,
    source = immune28_source
  ),
  file.path(step12_rdsdir, "ssgsea_obj.rds")
)

# STEP 05: Compare MYBL2 groups.

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
  file.path(step12_tabledir, "ssgsea_group_tests.csv"),
  row.names = FALSE
)

write.csv(
  immune_wilcox %>% dplyr::filter(Significant_FDR),
  file.path(step12_tabledir, "ssgsea_group_sig.csv"),
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
  "fig06b_ssgsea_box.png"
)

save_plot(
  p_immune_box,
  basename(fig_06b_file),
  width = 16.8,
  height = 10.9
)

# STEP 06: Test MYBL2 immune correlations.

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
  file.path(step12_tabledir, "ssgsea_cor.csv"),
  row.names = FALSE
)

write.csv(
  immune_cor %>% dplyr::filter(Significant_FDR_rule),
  file.path(step12_tabledir, "ssgsea_cor_sig.csv"),
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
  "fig06c_ssgsea_cor.png"
)

save_plot(
  p_immune_cor,
  basename(fig_06c_file),
  width = 9.4,
  height = 16.9
)

# STEP 07: Run ESTIMATE analysis.

estimate_gene_by_sample <- t(estimate_cpm_sample_by_gene)

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
  file.path(step12_tabledir, "estimate_scores.csv"),
  row.names = FALSE
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
  file.path(step12_tabledir, "estimate_tests.csv"),
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
  "fig06d_estimate_group.png"
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
  file.path(step12_tabledir, "estimate_cor.csv"),
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
  "fig06e_estimate_cor.png"
)

save_plot(
  p_estimate_cor,
  basename(fig_06e_file),
  width = 14.0,
  height = 7.2
)

# STEP 08: Build final Figure 6.

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
  "fig06a_ext_roc.png"
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
  file.path(step12_tabledir, "ext_roc_auc.csv"),
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
  "fig06_immune.png"
)

save_magick_all_formats(
  image_object = fig6_combined_img,
  filename_stem = "fig06_immune",
  dir_path = FIGURE_DIR
)
