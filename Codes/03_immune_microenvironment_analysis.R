# LSCC MYBL2 immune-microenvironment analysis
# Includes ssGSEA, ESTIMATE, tumor-purity adjustment, and dataset sensitivity.

# Step 1: clean environment, settings, packages, and paths

rm(list = ls())
gc()

set.seed(123)
options(stringsAsFactors = FALSE)
options(scipen = 100)
options(timeout = 7200)

TARGET_GENE <- "MYBL2"

# Figure settings

# Fonts and export quality
FONT_FAMILY <- "Arial"
FIG_DPI <- 600
FIG_BACKGROUND <- "white"

# Figure export policy
# FALSE = do not save any intermediate figure PNGs.
# The final composite figure is saved explicitly near the end of this script.
SAVE_INTERMEDIATE_FIGURES <- FALSE

# Standard manuscript figure dimensions (inches)
FIG_SINGLE_W <- 3.50
FIG_SINGLE_H <- 4.20
FIG_DOUBLE_W <- 7.20
FIG_DOUBLE_H <- 5.40
FIGURE2_COMBINED_W <- 14.00
FIGURE2_COMBINED_H <- 6.40

# Typography (pt)
BASE_TEXT_PT <- 9.5
AXIS_TEXT_PT <- 8.5
AXIS_TITLE_PT <- 10
LEGEND_TEXT_PT <- 8.5
LEGEND_TITLE_PT <- 9
PLOT_TITLE_PT <- 10
PANEL_TAG_PT <- 14

# Shared visual rules for every Figure 6 panel.
PANEL_BORDER_LWD <- 0.45
AXIS_LWD <- 0.40
TICK_LWD <- 0.40
GEOM_LWD <- 0.55
GRID_LWD <- 0.30
TICK_LENGTH_PT <- 2
PANEL_SPACING_PT <- 6

# Raster-composed panel tags and frames (panels a-e).
RASTER_PANEL_TAG_SIZE <- 145
RASTER_PANEL_TAG_OFFSET_X <- 100
RASTER_PANEL_TAG_OFFSET_Y <- 82
RASTER_PANEL_BORDER_PX <- 3
RASTER_PANEL_INSET_X <- 280
# Smaller vertical inset lets the actual plots occupy more of each final panel.
RASTER_PANEL_INSET_Y <- 80

# Heatmap typography
HEATMAP_FONT_PT <- 8
HEATMAP_ROW_FONT_PT <- 7
HEATMAP_COL_FONT_PT <- 7
HEATMAP_LEGEND_FONT_PT <- 8

# Repeated semantic colours
# Navy / rose palette for immune-microenvironment panels b-e.
# MYBL2-low and negative correlations use indigo-blue; MYBL2-high and
# positive correlations use rose. Panel a uses complementary navy / copper.
COL_NORMAL <- "#5B6FB2"   # indigo-blue: MYBL2-low / negative correlations
COL_TUMOR <- "#CE5C8A"    # rose: MYBL2-high / positive correlations
COL_MYBL2_LOW <- COL_NORMAL
COL_MYBL2_HIGH <- COL_TUMOR
COL_UP <- COL_TUMOR
COL_DOWN <- COL_NORMAL
COL_NS <- "#D1D3DA"

CELLCHAT_PALETTE <- grDevices::colorRampPalette(
  c("#F4F1F8", "#DDD5EA", "#BDA9D5", "#8B70B4", "#573D7E")
)(100)

# Equal internal whitespace around every ggplot panel.
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

# Statistical thresholds.
# BH-FDR < 0.05 defines statistical significance after multiple-testing correction.
# |rho| > 0.30 is NOT an additional inferential significance requirement; it is
# used only to describe/highlight correlations of at least moderate magnitude in
# the unadjusted descriptive analyses.
FDR_CUTOFF <- 0.05
CORRELATION_RHO_CUTOFF <- 0.30

# Figure refinement settings.
# Panel B is deliberately restricted to exactly 15 signatures for manuscript
# readability. The 15 FDR-significant signatures are shown when at least 15
# are available; otherwise, the panel is completed with the next FDR-ranked
# signatures. All 28 results remain available in the exported CSV tables.
PANEL_B_N_SIGNATURES <- 15L
# Complete signature names are retained in panel B. Wrapping avoids
# abbreviations while keeping the labels readable at manuscript scale.
PANEL_B_X_TEXT_SIZE <- 8.0
PANEL_B_X_TEXT_ANGLE <- 45
PANEL_B_LABEL_WRAP_WIDTH <- 19
PANEL_B_BOX_WIDTH <- 0.42
PANEL_B_DODGE_WIDTH <- 0.52
PANEL_B_POINT_SIZE <- 0.72

# Compact correlation display.
# The full 28-cell statistical analysis is preserved in the CSV outputs.
# The strongest 12 purity-adjusted associations are displayed in Panel C.
# The same Top-12 limit may also be used in the standalone sensitivity figure,
# but that sensitivity figure is NOT part of the final combined A/B/C figure.
PANEL_CD_N_CELLS <- 12L

# Official Charoentong et al. Cell Reports 2017 supplementary file already
# downloaded locally by the user. This script reads the local file directly and
# never attempts to download it from Cell Press / Elsevier.
CHAROENTONG_MMC3_FILE <- if (file.exists("E:/LSCC/Results_LSCC/mmc3.xlsx")) {
  "E:/LSCC/Results_LSCC/mmc3.xlsx"
} else {
  "D:/LSCC/Results_LSCC/mmc3.xlsx"
}
CHAROENTONG_SHEET <- "Sheet1"
CHAROENTONG_REFERENCE_CITATION <- paste0(
  "Charoentong et al. 2017, Cell Reports 18:248-262; ",
  "Table S6 (Pancancer immune metagenes), local mmc3.xlsx"
)
MIN_SIGNATURE_GENE_OVERLAP <- 3

# INPUT AND OUTPUT PATHS
# The bulk-RNA-seq/ML pipeline has already created the corrected matrices here.
ml_input_dir <- if (dir.exists("E:/LSCC/Results_LSCC/ML")) {
  "E:/LSCC/Results_LSCC/ML"
} else {
  "D:/LSCC/Results_LSCC/ML"
}

# Every output from this immune analysis is saved directly in this exact folder.
# No figures/, tables/, rds/, Step_*, manuscript-copy, or style subfolders
# are created.
immune_output_dir <- file.path(
  "E:/LSCC/Results_LSCC/ML",
  "Purity_Adjusted_Immune_Analysis"
)
dir.create(immune_output_dir, recursive = TRUE, showWarnings = FALSE)

# These aliases deliberately point to the same one output folder so all PNG,
# CSV, RDS, TXT, and session-information files are stored together.
step12_dir <- immune_output_dir
step12_figdir <- immune_output_dir
step12_tabledir <- immune_output_dir
step12_rdsdir <- immune_output_dir
manuscript_figures_path <- immune_output_dir

# PACKAGES
# All required packages must already be installed in the active R environment.

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

# Step 2: helper functions

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

# Complete publication-friendly signature names for panel B.
# Labels are NOT abbreviated; they are simply wrapped onto multiple lines.
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

# Partial Spearman correlation used for tumor-purity adjustment:
# - rank-transform x and y;
# - rank-transform numeric covariates (e.g. estimated tumor purity);
# - retain categorical covariates (e.g. Dataset) as factors;
# - regress ranked x and ranked y on the covariates;
# - correlate the residuals and calculate the partial-correlation P value.
partial_spearman <- function(x, y, covariates) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  y <- suppressWarnings(as.numeric(as.character(y)))

  covariates <- as.data.frame(
    covariates,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  dat <- data.frame(
    X = x,
    Y = y,
    covariates,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  numeric_cols <- vapply(dat, is.numeric, logical(1))
  finite_ok <- rep(TRUE, nrow(dat))

  for (nm in names(dat)[numeric_cols]) {
    finite_ok <- finite_ok & is.finite(dat[[nm]])
  }

  keep <- complete.cases(dat) & finite_ok
  dat <- dat[keep, , drop = FALSE]

  if (
    nrow(dat) < 6 ||
    length(unique(dat$X)) < 2 ||
    length(unique(dat$Y)) < 2
  ) {
    return(
      data.frame(
        n = nrow(dat),
        partial_rho = NA_real_,
        p_value = NA_real_,
        df = NA_real_,
        covariate_df = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  dat$X_rank <- rank(dat$X, ties.method = "average")
  dat$Y_rank <- rank(dat$Y, ties.method = "average")

  cov_names <- colnames(covariates)
  cov_ranked <- data.frame(row.names = seq_len(nrow(dat)))

  for (nm in cov_names) {
    z <- dat[[nm]]

    if (is.numeric(z)) {
      if (length(unique(z)) < 2) next
      cov_ranked[[nm]] <- rank(z, ties.method = "average")
    } else {
      z <- factor(z)
      if (nlevels(z) < 2) next
      cov_ranked[[nm]] <- z
    }
  }

  if (ncol(cov_ranked) == 0) {
    raw <- safe_spearman(dat$X, dat$Y)
    return(
      data.frame(
        n = nrow(dat),
        partial_rho = raw$rho,
        p_value = raw$p_value,
        df = nrow(dat) - 2,
        covariate_df = 0,
        stringsAsFactors = FALSE
      )
    )
  }

  model_df <- data.frame(
    X_rank = dat$X_rank,
    Y_rank = dat$Y_rank,
    cov_ranked,
    check.names = FALSE
  )

  cov_formula <- paste(names(cov_ranked), collapse = " + ")

  fit_x <- stats::lm(
    stats::as.formula(paste("X_rank ~", cov_formula)),
    data = model_df
  )

  fit_y <- stats::lm(
    stats::as.formula(paste("Y_rank ~", cov_formula)),
    data = model_df
  )

  rx <- stats::residuals(fit_x)
  ry <- stats::residuals(fit_y)

  partial_rho <- suppressWarnings(
    stats::cor(rx, ry, method = "pearson")
  )

  mm <- stats::model.matrix(
    stats::as.formula(paste("~", cov_formula)),
    data = model_df
  )

  k <- qr(mm)$rank - 1L
  df_t <- nrow(model_df) - k - 2L

  if (!is.finite(partial_rho) || df_t <= 0) {
    p_value <- NA_real_
  } else if (abs(partial_rho) >= 1) {
    p_value <- 0
  } else {
    t_stat <- partial_rho * sqrt(df_t / (1 - partial_rho^2))
    p_value <- 2 * stats::pt(-abs(t_stat), df = df_t)
  }

  data.frame(
    n = nrow(model_df),
    partial_rho = partial_rho,
    p_value = p_value,
    df = df_t,
    covariate_df = k,
    stringsAsFactors = FALSE
  )
}

# Approximate 95% CI for a partial correlation using the Fisher-z transform.
# partial_spearman() returns df = n - k - 2; therefore the Fisher-z standard
# error is approximated by 1/sqrt(n-k-3) = 1/sqrt(df-1). This is used for
# visualization only; P values/FDR remain those calculated above.
partial_rho_fisher_ci <- function(rho, df, conf_level = 0.95) {
  rho <- suppressWarnings(as.numeric(rho))
  df <- suppressWarnings(as.numeric(df))

  lower <- rep(NA_real_, length(rho))
  upper <- rep(NA_real_, length(rho))

  valid <- is.finite(rho) &
    is.finite(df) &
    abs(rho) < 1 &
    df > 1

  if (any(valid)) {
    z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)
    z_r <- atanh(rho[valid])
    se_z <- 1 / sqrt(df[valid] - 1)
    lower[valid] <- tanh(z_r - z_crit * se_z)
    upper[valid] <- tanh(z_r + z_crit * se_z)
  }

  boundary <- is.finite(rho) & abs(rho) >= 1
  if (any(boundary)) {
    bounded <- pmax(-1, pmin(1, rho[boundary]))
    lower[boundary] <- bounded
    upper[boundary] <- bounded
  }

  data.frame(
    CI_low = lower,
    CI_high = upper,
    stringsAsFactors = FALSE
  )
}

save_plot <- function(plot_obj, filename, width = 8, height = 6) {
  # In the final-image-only version, intermediate/legacy figures are calculated
  # as ggplot objects when needed by downstream composition, but are NOT saved.
  if (!isTRUE(SAVE_INTERMEDIATE_FIGURES)) {
    return(invisible(NULL))
  }

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

  # In this one-folder version, the source is normally already the final file.
  # Avoid trying to copy a file onto itself.
  if (identical(
    normalizePath(from_file, winslash = "/", mustWork = FALSE),
    normalizePath(to_file, winslash = "/", mustWork = FALSE)
  )) {
    return(invisible(to_file))
  }

  file.copy(from = from_file, to = to_file, overwrite = TRUE)
  invisible(to_file)
}

# RASTER FIGURE-COMPOSITION HELPERS
# These helpers implement the final Figure 6 image layout.

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

  # Identical usable canvas and inner whitespace for panels a-e.
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

  # No outer raster border is added. This intentionally removes the black
  # separator/frame lines between A/B/C/D in the combined figure.
  add_panel_label_img(
    img = img,
    label = label,
    size = label_size
  )
}

save_magick_all_formats <- function(image_object,
                                    filename_stem,
                                    dir_path) {
  if (!isTRUE(SAVE_INTERMEDIATE_FIGURES)) {
    return(invisible(list(PNG = NA_character_)))
  }

  dir.create(
    dir_path,
    recursive = TRUE,
    showWarnings = FALSE
  )

  png_file <- file.path(
    dir_path,
    paste0(
      filename_stem,
      ".png"
    )
  )

  magick::image_write(
    image_object,
    path = png_file,
    format = "png",
    density = paste0(
      FIG_DPI,
      "x",
      FIG_DPI
    )
  )

  invisible(
    list(
      PNG = png_file
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

      # Main titles: same font, bold state, size, and centered position.
      plot.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PLOT_TITLE_PT,
        hjust = 0.5,
        vjust = 0.5,
        color = "black"
      ),
      plot.title.position = "panel",

      # Axis labels, numbers, and ticks: one visual standard.
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

      # Facet titles use the same Arial bold treatment.
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

      # Legends use one consistent typography and spacing convention.
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

      # Every ggplot panel has the same white field and thin black border.
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

      # One internal margin and panel spacing standard.
      panel.spacing = grid::unit(PANEL_SPACING_PT, "pt"),
      plot.margin = MANUSCRIPT_MARGIN
    )
}

# MINIMAL THEME FOR THE MAIN PURITY-FIRST FIGURE 6

make_minimal_main_theme <- function(
    base_size = 10,
    legend_position = "bottom",
    show_y_grid = TRUE) {

  ggplot2::theme_classic(
    base_size = base_size,
    base_family = FONT_FAMILY
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = FONT_FAMILY,
        color = "black"
      ),
      axis.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = 10.2,
        color = "black"
      ),
      axis.text = ggplot2::element_text(
        family = FONT_FAMILY,
        size = 8.7,
        color = "black"
      ),
      axis.line = ggplot2::element_line(
        color = "black",
        linewidth = 0.42
      ),
      axis.ticks = ggplot2::element_line(
        color = "black",
        linewidth = 0.35
      ),
      axis.ticks.length = grid::unit(1.8, "pt"),
      panel.grid.major.y = if (show_y_grid) {
        ggplot2::element_line(
          color = "grey93",
          linewidth = 0.28
        )
      } else {
        ggplot2::element_blank()
      },
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = 10.2,
        color = "black",
        hjust = 0
      ),
      legend.position = legend_position,
      legend.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = 8.4
      ),
      legend.text = ggplot2::element_text(
        family = FONT_FAMILY,
        size = 8.1
      ),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = 15,
        color = "black"
      ),
      plot.tag.position = c(0.01, 0.99),
      plot.margin = ggplot2::margin(
        6, 7, 6, 7,
        unit = "pt"
      )
    )
}

# OFFICIAL LOCAL CHAROENTONG / TCIA 28-IMMUNE-CELL METAGENE SOURCE
# Local file structure confirmed by the user:
#   Sheet1
#   Row 1: "Table S6. List of Pancancer immune metagenes."
#   Row 3: Metagene | Cell type | Immunity
#   Data begin below that header row.
# Sheet2 and Sheet3 are empty and are intentionally ignored.

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

# These aliases allow harmless differences in spacing, capitalization, plural
# form, and the original workbook typo "memeory" while retaining the standard
# 28 labels in the final figures and tables.
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

# Sanity check for the known spelling used in the official local workbook.
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

  # Read Sheet1 with no assumed header, then locate the actual row that contains
  # the official column names "Metagene" and "Cell type".
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

# GSVA compatibility wrapper for ssGSEA (not ordinary GSVA).
# Recent GSVA versions require ssgseaParam(); older versions use the legacy API.
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

# Input expected by IOBR: genes in rows and samples in columns.
# ESTIMATE is supplied directly with original non-log TMM-CPM values.
# No batch-effect correction is applied.
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

# Step 3: pooled tumor cohort without batch-effect correction

# This input is created by Step 02.5 of the second bulk RNA-seq / ML pipeline.
# It contains TMM-normalized, non-negative CPM values for the merged discovery
# datasets. No batch-effect correction is performed in this immune analysis.
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

# Pool all tumor samples from both discovery datasets.
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

# Original non-log TMM-CPM values: used for ESTIMATE.
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

# ssGSEA input: log2(TMM-CPM + 1), with no batch-effect correction.
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

# Global MYBL2 median defines MYBL2-high and MYBL2-low groups.
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

# Step 4: ssgsea of 28 immune cell types

# Read the locally stored official Charoentong Table S6 source.
# The original mmc3.xlsx file is never downloaded from the Internet.
official_immune28 <- read_local_charoentong_immune28()
immune28_sets <- official_immune28$gene_sets
immune28_source <- CHAROENTONG_REFERENCE_CITATION
immune28_source_names <- official_immune28$source_names

official_gmt_file <- file.path(
  step12_tabledir,
  "Charoentong_2017_28_immune_cell_metagenes_official_source.gmt"
)
write_gmt(immune28_sets, official_gmt_file)

# Keep an exact local copy of the original workbook with the analysis outputs.
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

# Keep only official Charoentong genes represented in the LSCC bulk expression matrix.
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

# Require a minimum overlap, but use three genes rather than five because some
# rare immune signatures are sparsely represented in bulk microarray data.
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

# Enforce intended order of the standard 28 immune cells.
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

# Step 5: mybl2-high vs mybl2-low immune comparison, bh-fdr, and boxplots

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

# Apply BH correction across all 28 immune-signature Wilcoxon tests.
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

# The descriptive ssGSEA heatmap is intentionally omitted from Figure 6 and is
# not exported as a supplementary figure. The complete 28-cell ssGSEA matrix
# remains available in ssGSEA_28_immune_cell_scores_tumor_samples.csv.

# Supplementary raw ssGSEA panel: exactly 15 signatures for a compact, readable
# manuscript panel. When at least 15 signatures meet BH-FDR < 0.05, the 15
# most significant are plotted. If fewer than 15 pass BH-FDR, the remaining
# slots are completed with the next FDR-ranked signatures and are shown without
# significance stars.
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
    "Supplementary raw ssGSEA panel displays the ",
    PANEL_B_N_SIGNATURES,
    " most significant immune-cell signatures (BH-FDR < 0.05)."
  )
} else {
  boxplot_cells <- immune_wilcox_ranked %>%
    dplyr::slice_head(n = PANEL_B_N_SIGNATURES) %>%
    dplyr::pull(Immune_cell) %>%
    as.character()

  boxplot_panel_note <- paste0(
    "Supplementary raw ssGSEA panel displays the top ",
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
  "Supplementary_Figure_Raw_A_Top15_BH_FDR_ssGSEA_immune_cell_boxplot.png"
)

# Aspect ratio matches panel b: 7400 x 4800. Larger device height makes
# the boxplot itself occupy the final panel rather than leaving blank space.
save_plot(
  p_immune_box,
  basename(fig_06b_file),
  width = 16.8,
  height = 10.9
)

# Step 6: mybl2–immune-cell spearman correlation with bh-fdr

# All 28 correlations are calculated. BH-FDR is applied across all 28
# correlation P-values. A correlation is highlighted only when:
# |rho| > 0.30 AND BH-FDR < 0.05.
immune_cor <- immune_long %>%
  dplyr::group_by(Immune_cell) %>%
  dplyr::summarise(
    tmp = list(safe_spearman(MYBL2_expression, ssGSEA_score)),
    .groups = "drop"
  ) %>%
  tidyr::unnest(tmp) %>%
  dplyr::mutate(
    BH_FDR = p.adjust(p_value, method = "BH"),
    Significant_BH_FDR =
      is.finite(BH_FDR) &
      BH_FDR < FDR_CUTOFF,
    Moderate_Magnitude =
      is.finite(rho) &
      abs(rho) > CORRELATION_RHO_CUTOFF,
    Highlighted_Moderate_FDR =
      Significant_BH_FDR &
      Moderate_Magnitude,
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
  immune_cor %>% dplyr::filter(Significant_BH_FDR),
  file.path(
    step12_tabledir,
    "MYBL2_ssGSEA28_Spearman_BH_FDR_significant_correlations.csv"
  ),
  row.names = FALSE
)

write.csv(
  immune_cor %>% dplyr::filter(Highlighted_Moderate_FDR),
  file.path(
    step12_tabledir,
    "MYBL2_ssGSEA28_Spearman_ModerateMagnitude_and_BH_FDR_highlighted.csv"
  ),
  row.names = FALSE
)

# Supplementary raw correlation panel retains all 28 signatures.
# Statistical significance is defined by BH-FDR < 0.05. For visualization only,
# correlations are coloured/labeled when they also have |rho| > 0.30, indicating
# at least moderate magnitude. The rho cutoff is therefore descriptive, not an
# additional inferential significance threshold.
immune_cor_plot <- immune_cor %>%
  dplyr::filter(is.finite(rho)) %>%
  dplyr::mutate(
    Immune_cell_label = factor(
      as.character(Immune_cell),
      levels = rev(immune28_labels)
    ),
    Correlation_class = dplyr::case_when(
      Highlighted_Moderate_FDR & rho > 0 ~ "Moderate + BH-FDR positive",
      Highlighted_Moderate_FDR & rho < 0 ~ "Moderate + BH-FDR negative",
      TRUE ~ "Not highlighted"
    ),
    rho_label = ifelse(
      Highlighted_Moderate_FDR,
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
    data = dplyr::filter(immune_cor_plot, Highlighted_Moderate_FDR),
    aes(label = rho_label),
    size = 3.10,
    color = "black",
    fontface = "plain"
  ) +
  scale_color_manual(
    values = c(
      "Moderate + BH-FDR positive" = COL_TUMOR,
      "Moderate + BH-FDR negative" = COL_NORMAL,
      "Not highlighted" = "grey78"
    ),
    breaks = c(
      "Moderate + BH-FDR positive",
      "Moderate + BH-FDR negative",
      "Not highlighted"
    ),
    name = "Raw descriptive highlight"
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
  "Supplementary_Figure_Raw_C_MYBL2_ssGSEA28_BH_FDR_correlation_bubbleplot.png"
)

# Aspect ratio matches panel c: 4000 x 7200.
save_plot(
  p_immune_cor,
  basename(fig_06c_file),
  width = 9.4,
  height = 16.9
)

# Step 7: estimate scores, group comparisons, and correlations

# ESTIMATE is run directly on original non-log TMM-CPM values.
# ssGSEA uses log2(TMM-CPM + 1); no batch-effect correction is applied.
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

# BH-FDR is also applied to the three ESTIMATE group-comparison tests.
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
  "Supplementary_Figure_Raw_B_ESTIMATE_score_MYBL2_high_low_BH_FDR_comparison.png"
)

# Aspect ratio matches panel d: 7000 x 3600.
save_plot(
  p_estimate_group,
  basename(fig_06d_file),
  width = 13.4,
  height = 6.9
)

# BH-FDR is also applied across the three ESTIMATE correlation tests.
estimate_cor <- estimate_long %>%
  dplyr::group_by(TME_score) %>%
  dplyr::summarise(
    tmp = list(safe_spearman(MYBL2_expression, Score)),
    .groups = "drop"
  ) %>%
  tidyr::unnest(tmp) %>%
  dplyr::mutate(
    BH_FDR = p.adjust(p_value, method = "BH"),
    Significant_BH_FDR =
      is.finite(BH_FDR) &
      BH_FDR < FDR_CUTOFF,
    Moderate_Magnitude =
      is.finite(rho) &
      abs(rho) > CORRELATION_RHO_CUTOFF,
    Highlighted_Moderate_FDR =
      Significant_BH_FDR &
      Moderate_Magnitude,
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
  "Supplementary_Figure_Raw_D_MYBL2_ESTIMATE_score_BH_FDR_correlations.png"
)

# Aspect ratio matches panel e: 7000 x 3600.
save_plot(
  p_estimate_cor,
  basename(fig_06e_file),
  width = 14.0,
  height = 7.2
)

# Step 7B: tumor-purity-adjusted analysis

# Published ESTIMATE-derived tumor-purity transformation:
# Tumor purity = cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)
ESTIMATE_PURITY_INTERCEPT <- 0.6049872018
ESTIMATE_PURITY_SLOPE <- 0.0001467884

estimate_scores$Tumor_Purity_ESTIMATE <- cos(
  ESTIMATE_PURITY_INTERCEPT +
    ESTIMATE_PURITY_SLOPE *
    estimate_scores$ESTIMATE_score
)

if (
  any(
    !is.finite(estimate_scores$Tumor_Purity_ESTIMATE) |
    estimate_scores$Tumor_Purity_ESTIMATE < 0 |
    estimate_scores$Tumor_Purity_ESTIMATE > 1
  )
) {
  stop("Non-finite or out-of-range ESTIMATE-derived tumor-purity values were detected.")
}

write.csv(
  estimate_scores,
  file.path(
    step12_tabledir,
    "01_ESTIMATE_Derived_Tumor_Purity_Per_Sample.csv"
  ),
  row.names = FALSE
)

# 12.7B.1. MYBL2 vs estimated tumor purity

mybl2_purity_cor <- safe_spearman(
  estimate_scores$MYBL2_expression,
  estimate_scores$Tumor_Purity_ESTIMATE
)

mybl2_purity_table <- data.frame(
  N = mybl2_purity_cor$n,
  Spearman_Rho = mybl2_purity_cor$rho,
  P_Value = mybl2_purity_cor$p_value,
  stringsAsFactors = FALSE
)

write.csv(
  mybl2_purity_table,
  file.path(
    step12_tabledir,
    "02_MYBL2_vs_ESTIMATE_Derived_Tumor_Purity_Spearman.csv"
  ),
  row.names = FALSE
)

# 12.7B.2. Raw vs purity-adjusted MYBL2–ssGSEA associations

immune_long_purity <- immune_long %>%
  dplyr::left_join(
    estimate_scores %>%
      dplyr::select(
        Sample,
        Tumor_Purity_ESTIMATE
      ),
    by = "Sample"
  )

immune_cells_for_adjustment <- unique(
  as.character(
    immune_long_purity$Immune_cell
  )
)

immune_adjusted_rows <- lapply(
  immune_cells_for_adjustment,
  function(cell_name) {

    dat <- immune_long_purity %>%
      dplyr::filter(
        as.character(Immune_cell) == cell_name
      )

    raw_res <- safe_spearman(
      dat$MYBL2_expression,
      dat$ssGSEA_score
    )

    purity_res <- partial_spearman(
      x = dat$MYBL2_expression,
      y = dat$ssGSEA_score,
      covariates = data.frame(
        Tumor_Purity_ESTIMATE = dat$Tumor_Purity_ESTIMATE
      )
    )

    purity_dataset_res <- partial_spearman(
      x = dat$MYBL2_expression,
      y = dat$ssGSEA_score,
      covariates = data.frame(
        Tumor_Purity_ESTIMATE = dat$Tumor_Purity_ESTIMATE,
        Dataset = dat$Dataset
      )
    )

    data.frame(
      Immune_cell = cell_name,
      N = raw_res$n,
      Raw_Spearman_Rho = raw_res$rho,
      Raw_P_Value = raw_res$p_value,
      Purity_Adjusted_N = purity_res$n,
      Purity_Adjusted_df = purity_res$df,
      Purity_Adjusted_Partial_Spearman_Rho = purity_res$partial_rho,
      Purity_Adjusted_P_Value = purity_res$p_value,
      Purity_Plus_Dataset_Adjusted_N = purity_dataset_res$n,
      Purity_Plus_Dataset_Adjusted_df = purity_dataset_res$df,
      Purity_Plus_Dataset_Adjusted_Partial_Spearman_Rho =
        purity_dataset_res$partial_rho,
      Purity_Plus_Dataset_Adjusted_P_Value =
        purity_dataset_res$p_value,
      stringsAsFactors = FALSE
    )
  }
)

immune_results <- dplyr::bind_rows(
  immune_adjusted_rows
) %>%
  dplyr::mutate(
    Raw_BH_FDR = p.adjust(
      Raw_P_Value,
      method = "BH"
    ),
    Purity_Adjusted_BH_FDR = p.adjust(
      Purity_Adjusted_P_Value,
      method = "BH"
    ),
    Purity_Plus_Dataset_Adjusted_BH_FDR = p.adjust(
      Purity_Plus_Dataset_Adjusted_P_Value,
      method = "BH"
    ),
    Raw_Significant_FDR =
      is.finite(Raw_BH_FDR) &
      Raw_BH_FDR < FDR_CUTOFF,
    Raw_Moderate_Magnitude =
      is.finite(Raw_Spearman_Rho) &
      abs(Raw_Spearman_Rho) > CORRELATION_RHO_CUTOFF,
    Raw_Highlighted_Moderate_FDR =
      Raw_Significant_FDR &
      Raw_Moderate_Magnitude,
    Purity_Adjusted_Significant_FDR =
      is.finite(Purity_Adjusted_BH_FDR) &
      Purity_Adjusted_BH_FDR < FDR_CUTOFF,
    Purity_Plus_Dataset_Adjusted_Significant_FDR =
      is.finite(Purity_Plus_Dataset_Adjusted_BH_FDR) &
      Purity_Plus_Dataset_Adjusted_BH_FDR < FDR_CUTOFF
  )

write.csv(
  immune_results,
  file.path(
    step12_tabledir,
    "03_ssGSEA28_Raw_vs_PurityAdjusted_PartialSpearman.csv"
  ),
  row.names = FALSE
)

write.csv(
  immune_results %>%
    dplyr::select(
      Immune_cell,
      N,
      Purity_Adjusted_N,
      Purity_Adjusted_df,
      Purity_Adjusted_Partial_Spearman_Rho,
      Purity_Adjusted_P_Value,
      Purity_Adjusted_BH_FDR,
      Purity_Adjusted_Significant_FDR
    ) %>%
    dplyr::arrange(
      Purity_Adjusted_BH_FDR
    ),
  file.path(
    step12_tabledir,
    "04_ssGSEA28_Primary_PurityAdjusted_Results.csv"
  ),
  row.names = FALSE
)

write.csv(
  immune_results %>%
    dplyr::select(
      Immune_cell,
      N,
      Purity_Plus_Dataset_Adjusted_N,
      Purity_Plus_Dataset_Adjusted_df,
      Purity_Plus_Dataset_Adjusted_Partial_Spearman_Rho,
      Purity_Plus_Dataset_Adjusted_P_Value,
      Purity_Plus_Dataset_Adjusted_BH_FDR,
      Purity_Plus_Dataset_Adjusted_Significant_FDR
    ) %>%
    dplyr::arrange(
      Purity_Plus_Dataset_Adjusted_BH_FDR
    ),
  file.path(
    step12_tabledir,
    "05_ssGSEA28_Purity_Plus_Dataset_Adjusted_Sensitivity.csv"
  ),
  row.names = FALSE
)

# 12.7B.3. ESTIMATE components: raw vs purity-adjusted

estimate_component_rows <- list()
idx <- 1L

for (
  component in c(
    "Immune_score",
    "Stromal_score"
  )
) {

  raw_res <- safe_spearman(
    estimate_scores$MYBL2_expression,
    estimate_scores[[component]]
  )

  adj_res <- partial_spearman(
    x = estimate_scores$MYBL2_expression,
    y = estimate_scores[[component]],
    covariates = data.frame(
      Tumor_Purity_ESTIMATE =
        estimate_scores$Tumor_Purity_ESTIMATE
    )
  )

  estimate_component_rows[[idx]] <- data.frame(
    Score = component,
    N = raw_res$n,
    Raw_Spearman_Rho = raw_res$rho,
    Raw_P_Value = raw_res$p_value,
    Purity_Adjusted_Partial_Spearman_Rho =
      adj_res$partial_rho,
    Purity_Adjusted_P_Value =
      adj_res$p_value,
    Note =
      "Partial Spearman adjusted for ESTIMATE-derived tumor purity",
    stringsAsFactors = FALSE
  )

  idx <- idx + 1L
}

raw_estimate_score <- safe_spearman(
  estimate_scores$MYBL2_expression,
  estimate_scores$ESTIMATE_score
)

estimate_component_rows[[idx]] <- data.frame(
  Score = "ESTIMATE_score",
  N = raw_estimate_score$n,
  Raw_Spearman_Rho = raw_estimate_score$rho,
  Raw_P_Value = raw_estimate_score$p_value,
  Purity_Adjusted_Partial_Spearman_Rho = NA_real_,
  Purity_Adjusted_P_Value = NA_real_,
  Note = paste0(
    "Not adjusted for purity because Tumor_Purity_ESTIMATE is a deterministic ",
    "transformation of ESTIMATE_score."
  ),
  stringsAsFactors = FALSE
)

estimate_component_results <- dplyr::bind_rows(
  estimate_component_rows
) %>%
  dplyr::mutate(
    Raw_BH_FDR =
      p.adjust(
        Raw_P_Value,
        method = "BH"
      ),
    Purity_Adjusted_BH_FDR =
      p.adjust(
        Purity_Adjusted_P_Value,
        method = "BH"
      )
  )

write.csv(
  estimate_component_results,
  file.path(
    step12_tabledir,
    "06_ESTIMATE_Components_Raw_vs_PurityAdjusted_Correlations.csv"
  ),
  row.names = FALSE
)

# 12.7B.4. MAIN FIGURE PANEL b — MYBL2 vs tumor purity

rho_purity <- mybl2_purity_table$Spearman_Rho[1]
p_purity <- mybl2_purity_table$P_Value[1]

p_purity_main <- ggplot(
  estimate_scores,
  aes(
    x = MYBL2_expression,
    y = Tumor_Purity_ESTIMATE,
    color = MYBL2_group,
    shape = Dataset
  )
) +
  geom_smooth(
    data = estimate_scores,
    aes(
      x = MYBL2_expression,
      y = Tumor_Purity_ESTIMATE
    ),
    inherit.aes = FALSE,
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black",
    fill = "grey82",
    alpha = 0.35,
    linewidth = 0.62,
    show.legend = FALSE
  ) +
  geom_point(
    size = 2.10,
    alpha = 0.78,
    stroke = 0.15
  ) +
  annotate(
    "label",
    x = -Inf,
    y = Inf,
    hjust = -0.04,
    vjust = 1.08,
    label = paste0(
      "Spearman rho = ",
      sprintf("%.3f", rho_purity),
      "\n",
      ifelse(
        p_purity < 0.001,
        "P < 0.001",
        paste0("P = ", sprintf("%.4f", p_purity))
      )
    ),
    size = 3.1,
    linewidth = 0,
    fill = "white",
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "MYBL2-low" = COL_MYBL2_LOW,
      "MYBL2-high" = COL_MYBL2_HIGH
    )
  ) +
  scale_shape_manual(
    values = c(
      "GSE127165" = 16,
      "GSE142083" = 17
    )
  ) +
  guides(
    color = guide_legend(
      title.position = "top",
      nrow = 1,
      order = 1
    ),
    shape = guide_legend(
      title.position = "top",
      nrow = 1,
      order = 2
    )
  ) +
  labs(
    x = "MYBL2 expression",
    y = "ESTIMATE-derived tumor purity",
    color = "MYBL2 group",
    shape = "Dataset"
  ) +
  make_minimal_main_theme(
    base_size = 10,
    legend_position = "bottom",
    show_y_grid = TRUE
  ) +
  theme(
    legend.box = "horizontal",
    legend.spacing.x = grid::unit(5, "pt")
  )

fig_06b_purity_file <- file.path(
  step12_figdir,
  "Figure_06B_MYBL2_vs_ESTIMATE_Derived_Tumor_Purity_MINIMAL.png"
)

save_plot(
  p_purity_main,
  basename(fig_06b_purity_file),
  width = 7.4,
  height = 4.8
)

# 12.7B.5. MAIN FIGURE PANEL c — COMPACT purity-adjusted forest plot
# Display-only selection: strongest 12 of the 28 tested immune-cell associations.
# All 28 tests remain in the exported result tables and BH-FDR is still calculated
# across the full family of 28 tests.

purity_ci <- partial_rho_fisher_ci(
  rho = immune_results$Purity_Adjusted_Partial_Spearman_Rho,
  df = immune_results$Purity_Adjusted_df,
  conf_level = 0.95
)

adjusted_plot_df <- immune_results %>%
  dplyr::mutate(
    Immune_cell = as.character(Immune_cell),
    Partial_rho = Purity_Adjusted_Partial_Spearman_Rho,
    CI_low = purity_ci$CI_low,
    CI_high = purity_ci$CI_high,
    BH_FDR = Purity_Adjusted_BH_FDR,
    Significant = is.finite(BH_FDR) & BH_FDR < FDR_CUTOFF,
    Correlation_class = dplyr::case_when(
      Significant & Partial_rho > 0 ~ "Significant positive",
      Significant & Partial_rho < 0 ~ "Significant negative",
      TRUE ~ "Not significant"
    )
  ) %>%
  dplyr::filter(is.finite(Partial_rho)) %>%
  # Compact display: strongest absolute adjusted associations, regardless of sign.
  dplyr::arrange(dplyr::desc(abs(Partial_rho)), BH_FDR) %>%
  dplyr::slice_head(n = PANEL_CD_N_CELLS) %>%
  dplyr::arrange(Partial_rho) %>%
  dplyr::mutate(
    Immune_cell = factor(
      Immune_cell,
      levels = Immune_cell
    )
  )

write.csv(
  adjusted_plot_df %>%
    dplyr::mutate(Immune_cell = as.character(Immune_cell)),
  file.path(
    step12_tabledir,
    paste0(
      "Figure_06C_Top",
      PANEL_CD_N_CELLS,
      "_PurityAdjusted_PartialSpearman_95CI.csv"
    )
  ),
  row.names = FALSE
)

p_immune_adjusted <- ggplot(
  adjusted_plot_df,
  aes(
    x = Partial_rho,
    y = Immune_cell
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey68",
    linewidth = 0.48
  ) +
  # Thicker 95% CI lines improve readability at manuscript scale.
  geom_segment(
    aes(
      x = CI_low,
      xend = CI_high,
      y = Immune_cell,
      yend = Immune_cell,
      color = Correlation_class
    ),
    linewidth = 1.35,
    alpha = 0.92,
    lineend = "round",
    na.rm = TRUE
  ) +
  geom_point(
    aes(color = Correlation_class),
    size = 3.7,
    alpha = 0.98
  ) +
  scale_color_manual(
    values = c(
      "Significant positive" = COL_TUMOR,
      "Significant negative" = COL_NORMAL,
      "Not significant" = "grey62"
    ),
    breaks = c(
      "Significant positive",
      "Significant negative",
      "Not significant"
    ),
    labels = c(
      "Sig. positive",
      "Sig. negative",
      "NS"
    ),
    name = NULL
  ) +
  # Explicit labels prevent long floating-point strings on the x axis.
  scale_x_continuous(
    breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
    labels = function(x) sprintf("%.1f", x),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        size = 3.0,
        linewidth = 1.0
      )
    )
  ) +
  labs(
    x = "Purity-adjusted partial rho (95% CI)",
    y = NULL
  ) +
  make_minimal_main_theme(
    base_size = 10.5,
    legend_position = "bottom",
    show_y_grid = FALSE
  ) +
  theme(
    axis.text.x = element_text(
      size = 9.3,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 10.3,
      color = "black",
      lineheight = 0.92
    ),
    axis.title.x = element_text(
      size = 10.0,
      face = "bold"
    ),
    legend.text = element_text(size = 6.9),
    legend.key.width = grid::unit(6, "pt"),
    legend.key.height = grid::unit(6, "pt"),
    legend.spacing.x = grid::unit(1.5, "pt"),
    legend.box.spacing = grid::unit(0, "pt"),
    legend.margin = ggplot2::margin(-2, 0, 0, 0, unit = "pt"),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
    legend.justification = "center",
    plot.margin = ggplot2::margin(1, 3, 1, 2, unit = "pt")
  ) +
  coord_cartesian(
    xlim = c(-0.65, 0.65),
    clip = "off"
  )

fig_06c_adjusted_file <- file.path(
  step12_figdir,
  paste0(
    "Figure_06C_MYBL2_ssGSEA_Top",
    PANEL_CD_N_CELLS,
    "_PurityAdjusted_FOREST_95CI_MINIMAL.png"
  )
)

save_plot(
  p_immune_adjusted,
  basename(fig_06c_adjusted_file),
  width = 9.2,
  height = 5.3
)

# 12.7B.6. MAIN FIGURE PANEL d — purity-adjusted Immune/Stromal scores

make_residualized_rank_df <- function(
    data,
    score_col,
    display_name) {

  dat <- data.frame(
    Sample = data$Sample,
    Dataset = data$Dataset,
    MYBL2_group = data$MYBL2_group,
    MYBL2_expression =
      suppressWarnings(as.numeric(data$MYBL2_expression)),
    Score =
      suppressWarnings(as.numeric(data[[score_col]])),
    Tumor_Purity_ESTIMATE =
      suppressWarnings(as.numeric(data$Tumor_Purity_ESTIMATE)),
    stringsAsFactors = FALSE
  )

  keep <- complete.cases(dat) &
    is.finite(dat$MYBL2_expression) &
    is.finite(dat$Score) &
    is.finite(dat$Tumor_Purity_ESTIMATE)

  dat <- dat[keep, , drop = FALSE]

  dat$MYBL2_rank <- rank(
    dat$MYBL2_expression,
    ties.method = "average"
  )

  dat$Score_rank <- rank(
    dat$Score,
    ties.method = "average"
  )

  dat$Purity_rank <- rank(
    dat$Tumor_Purity_ESTIMATE,
    ties.method = "average"
  )

  fit_x <- stats::lm(
    MYBL2_rank ~ Purity_rank,
    data = dat
  )

  fit_y <- stats::lm(
    Score_rank ~ Purity_rank,
    data = dat
  )

  dat$MYBL2_rank_residual <- stats::residuals(fit_x)
  dat$Score_rank_residual <- stats::residuals(fit_y)
  dat$TME_score <- display_name

  dat
}

adjusted_estimate_plot_df <- dplyr::bind_rows(
  make_residualized_rank_df(
    estimate_scores,
    "Stromal_score",
    "Stromal score"
  ),
  make_residualized_rank_df(
    estimate_scores,
    "Immune_score",
    "Immune score"
  )
)

adjusted_estimate_plot_df$TME_score <- factor(
  adjusted_estimate_plot_df$TME_score,
  levels = c(
    "Stromal score",
    "Immune score"
  )
)

adjusted_estimate_labels <- estimate_component_results %>%
  dplyr::filter(
    Score %in% c(
      "Stromal_score",
      "Immune_score"
    )
  ) %>%
  dplyr::mutate(
    TME_score = dplyr::recode(
      Score,
      Stromal_score = "Stromal score",
      Immune_score = "Immune score"
    ),
    Label = paste0(
      "partial rho = ",
      sprintf(
        "%.2f",
        Purity_Adjusted_Partial_Spearman_Rho
      ),
      "\n",
      vapply(
        Purity_Adjusted_BH_FDR,
        format_fdr,
        character(1)
      )
    )
  )

p_estimate_adjusted <- ggplot(
  adjusted_estimate_plot_df,
  aes(
    x = MYBL2_rank_residual,
    y = Score_rank_residual,
    color = MYBL2_group,
    shape = Dataset
  )
) +
  geom_smooth(
    mapping = aes(
      x = MYBL2_rank_residual,
      y = Score_rank_residual,
      group = TME_score
    ),
    inherit.aes = FALSE,
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black",
    fill = "grey82",
    alpha = 0.35,
    linewidth = 0.60,
    show.legend = FALSE
  ) +
  geom_point(
    size = 2.0,
    alpha = 0.76,
    stroke = 0.15
  ) +
  geom_label(
    data = adjusted_estimate_labels,
    aes(
      x = -Inf,
      y = Inf,
      label = Label
    ),
    inherit.aes = FALSE,
    hjust = -0.04,
    vjust = 1.07,
    size = 3.0,
    linewidth = 0,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(
    ~ TME_score,
    scales = "free_y",
    ncol = 1
  ) +
  scale_color_manual(
    values = c(
      "MYBL2-low" = COL_MYBL2_LOW,
      "MYBL2-high" = COL_MYBL2_HIGH
    )
  ) +
  scale_shape_manual(
    values = c(
      "GSE127165" = 16,
      "GSE142083" = 17
    )
  ) +
  guides(
    color = "none",
    shape = "none"
  ) +
  labs(
    x = "Residualized rank of MYBL2",
    y = "Residualized rank of ESTIMATE component"
  ) +
  make_minimal_main_theme(
    base_size = 10,
    legend_position = "none",
    show_y_grid = TRUE
  ) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 10.4,
      hjust = 0
    ),
    panel.spacing.y = grid::unit(12, "pt")
  )

fig_06d_adjusted_file <- file.path(
  step12_figdir,
  "Figure_06D_MYBL2_Immune_Stromal_PurityAdjusted_MINIMAL.png"
)

save_plot(
  p_estimate_adjusted,
  basename(fig_06d_adjusted_file),
  width = 7.2,
  height = 9.2
)

# 12.7B.7. SUPPLEMENTARY — COMPACT panels C and D
# C = Top-12 purity-adjusted partial Spearman forest plot with approximate 95% CI.
# D = Top-12 purity+dataset-adjusted sensitivity bubble plot.
# Selection affects visualization only; all 28-cell result tables are unchanged.

# Supplementary panel C: compact forest plot, using the same Top-12 selection
# and 95% CI display as the primary purity-adjusted panel C.

supp_c_forest_df <- adjusted_plot_df

p_supp_c_forest <- ggplot(
  supp_c_forest_df,
  aes(
    x = Partial_rho,
    y = Immune_cell
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.48,
    color = "grey68"
  ) +
  geom_segment(
    aes(
      x = CI_low,
      xend = CI_high,
      y = Immune_cell,
      yend = Immune_cell,
      color = Correlation_class
    ),
    linewidth = 1.35,
    alpha = 0.92,
    lineend = "round",
    na.rm = TRUE
  ) +
  geom_point(
    aes(color = Correlation_class),
    size = 3.7,
    alpha = 0.98
  ) +
  scale_color_manual(
    values = c(
      "Significant positive" = COL_TUMOR,
      "Significant negative" = COL_NORMAL,
      "Not significant" = "grey62"
    ),
    breaks = c(
      "Significant positive",
      "Significant negative",
      "Not significant"
    ),
    labels = c(
      "Sig. positive",
      "Sig. negative",
      "NS"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
    labels = function(x) sprintf("%.1f", x),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  # Tight y expansion reduces unused space above and below the 12 displayed cells.
  scale_y_discrete(
    # Minimal expansion makes the 12 cells fill the panel vertically and removes
    # the excess blank space above/below the first and last immune-cell rows.
    expand = expansion(add = c(0.03, 0.03))
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        size = 3.0,
        linewidth = 1.0
      )
    )
  ) +
  labs(
    x = "Purity-adjusted partial rho (95% CI)",
    y = NULL,
    title = NULL
  ) +
  make_minimal_main_theme(
    base_size = 10.5,
    legend_position = "bottom",
    show_y_grid = FALSE
  ) +
  theme(
    axis.text.x = element_text(size = 9.3, color = "black"),
    axis.text.y = element_text(
      size = 10.8,
      color = "black",
      lineheight = 0.94
    ),
    axis.title.x = element_text(size = 10.2, face = "bold"),
    legend.text = element_text(size = 6.6),
    legend.key.width = grid::unit(5.0, "pt"),
    legend.key.height = grid::unit(5.0, "pt"),
    legend.spacing.x = grid::unit(0.8, "pt"),
    legend.box.spacing = grid::unit(0, "pt"),
    legend.margin = ggplot2::margin(-5, 0, 0, 0, unit = "pt"),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
    legend.box = "horizontal",
    legend.justification = "center",
    plot.margin = ggplot2::margin(0.5, 2.5, 0.5, 1.5, unit = "pt")
  ) +
  coord_cartesian(
    xlim = c(-0.65, 0.65),
    clip = "off"
  )

fig_supp_c_compact_file <- file.path(
  step12_figdir,
  paste0(
    "Supplementary_Figure_C_Top",
    PANEL_CD_N_CELLS,
    "_PurityAdjusted_Forest_95CI.png"
  )
)

save_plot(
  p_supp_c_forest,
  basename(fig_supp_c_compact_file),
  width = 8.8,
  height = 5.4
)

# Supplementary panel D: compact Top-12 sensitivity analysis. The ranking is
# based on absolute purity+dataset-adjusted partial rho; BH-FDR remains based on
# all 28 tests and is represented by point colour.

sensitivity_ci <- partial_rho_fisher_ci(
  rho = immune_results$Purity_Plus_Dataset_Adjusted_Partial_Spearman_Rho,
  df = immune_results$Purity_Plus_Dataset_Adjusted_df,
  conf_level = 0.95
)

sensitivity_plot_df <- immune_results %>%
  dplyr::mutate(
    Immune_cell = as.character(Immune_cell),
    Partial_rho = Purity_Plus_Dataset_Adjusted_Partial_Spearman_Rho,
    CI_low = sensitivity_ci$CI_low,
    CI_high = sensitivity_ci$CI_high,
    BH_FDR = Purity_Plus_Dataset_Adjusted_BH_FDR,
    Significant = is.finite(BH_FDR) & BH_FDR < FDR_CUTOFF,
    Correlation_class = dplyr::case_when(
      Significant & Partial_rho > 0 ~ "Significant positive",
      Significant & Partial_rho < 0 ~ "Significant negative",
      TRUE ~ "Not significant"
    )
  ) %>%
  dplyr::filter(is.finite(Partial_rho)) %>%
  dplyr::arrange(dplyr::desc(abs(Partial_rho)), BH_FDR) %>%
  dplyr::slice_head(n = PANEL_CD_N_CELLS) %>%
  dplyr::arrange(Partial_rho) %>%
  dplyr::mutate(
    Immune_cell = factor(
      Immune_cell,
      levels = Immune_cell
    )
  )

write.csv(
  sensitivity_plot_df %>%
    dplyr::mutate(Immune_cell = as.character(Immune_cell)),
  file.path(
    step12_tabledir,
    paste0(
      "Supplementary_Panel_D_Top",
      PANEL_CD_N_CELLS,
      "_PurityPlusDataset_Adjusted.csv"
    )
  ),
  row.names = FALSE
)

p_sensitivity_adjusted <- ggplot(
  sensitivity_plot_df,
  aes(
    x = Partial_rho,
    y = Immune_cell
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.48,
    color = "grey68"
  ) +
  # Slightly thicker CI lines than before for clearer print reproduction.
  geom_segment(
    aes(
      x = CI_low,
      xend = CI_high,
      y = Immune_cell,
      yend = Immune_cell,
      color = Correlation_class
    ),
    linewidth = 1.15,
    alpha = 0.86,
    lineend = "round",
    na.rm = TRUE
  ) +
  geom_point(
    aes(
      size = abs(Partial_rho),
      color = Correlation_class
    ),
    alpha = 0.97
  ) +
  scale_color_manual(
    values = c(
      "Significant positive" = COL_TUMOR,
      "Significant negative" = COL_NORMAL,
      "Not significant" = "grey68"
    ),
    breaks = c(
      "Significant positive",
      "Significant negative",
      "Not significant"
    ),
    labels = c(
      "Sig. positive",
      "Sig. negative",
      "NS"
    ),
    name = NULL
  ) +
  scale_size_continuous(
    range = c(3.1, 6.2),
    breaks = c(0.15, 0.25, 0.35),
    labels = function(x) sprintf("%.2f", x),
    name = "|rho|"
  ) +
  # Five explicit rounded ticks prevent the long decimal labels seen previously.
  scale_x_continuous(
    breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
    labels = function(x) sprintf("%.1f", x),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  guides(
    size = guide_legend(
      title.position = "top",
      nrow = 1,
      order = 1
    ),
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 2,
      override.aes = list(size = 3.1)
    )
  ) +
  labs(
    x = "Purity + dataset-adjusted partial rho",
    y = NULL,
    title = NULL
  ) +
  make_minimal_main_theme(
    base_size = 10.5,
    legend_position = "bottom",
    show_y_grid = FALSE
  ) +
  theme(
    axis.text.x = element_text(size = 9.3, color = "black"),
    axis.text.y = element_text(
      size = 10.4,
      color = "black",
      lineheight = 0.92
    ),
    axis.title.x = element_text(size = 10.0, face = "bold"),
    legend.title = element_text(size = 6.8, face = "bold"),
    legend.text = element_text(size = 6.7),
    legend.key.width = grid::unit(6, "pt"),
    legend.key.height = grid::unit(6, "pt"),
    legend.spacing.x = grid::unit(1.2, "pt"),
    legend.box.spacing = grid::unit(0, "pt"),
    legend.margin = ggplot2::margin(-2, 0, 0, 0, unit = "pt"),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
    legend.box = "horizontal",
    legend.justification = "center",
    plot.margin = ggplot2::margin(1, 3, 1, 2, unit = "pt")
  ) +
  coord_cartesian(
    xlim = c(-0.65, 0.65),
    clip = "off"
  )

fig_supp_sensitivity_file <- file.path(
  step12_figdir,
  paste0(
    "Supplementary_Figure_D_Top",
    PANEL_CD_N_CELLS,
    "_PurityPlusDataset_Adjusted_Compact.png"
  )
)

save_plot(
  p_sensitivity_adjusted,
  basename(fig_supp_sensitivity_file),
  width = 9.4,
  height = 5.3
)

# Step 8: final combined manuscript figure 6 — purity first
# Main Figure 6:
#   a) External validation ROC (RF + SVM)
#   b) MYBL2 vs ESTIMATE-derived tumor purity
#   c) Compact Top-12 (of 28 tested) purity-adjusted partial Spearman forest plot
#   d) Purity-adjusted MYBL2 vs Immune/Stromal scores
# The original raw high/low ssGSEA and ESTIMATE plots remain exported as
# supplementary/descriptive results, but are no longer the central main figure.

# PANEL A. REBUILD EXTERNAL ROC PANEL (NO TITLE BANNER / UPDATED COLOURS)
# This panel is generated directly from the training and independent GSE130605
# matrices so no highlighted title from a previously exported image is retained.
# The model specifications match the manuscript Methods: RF (500 trees) and
# linear-kernel SVM, trained on MYBL2 in the training cohort and evaluated in
# GSE130605. Standardisation parameters are estimated only from training data.

COL_RF_ROC  <- "#274C77"  # deep navy
COL_SVM_ROC <- "#D06B43"  # warm copper

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

# Supports both numeric group coding (1 = normal, 2 = tumor) and text labels.
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
  aes(
    x = FPR,
    y = Sensitivity,
    colour = Model
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey70",
    linewidth = 0.45
  ) +
  geom_step(
    linewidth = 0.95
  ) +
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  scale_colour_manual(
    values = setNames(
      c(
        COL_RF_ROC,
        COL_SVM_ROC
      ),
      levels(roc_df$Model)
    ),
    name = NULL
  ) +
  guides(
    colour = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  ) +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  make_minimal_main_theme(
    base_size = 10,
    legend_position = c(0.68, 0.18),
    show_y_grid = TRUE
  ) +
  theme(
    legend.background = element_rect(
      fill = "white",
      color = NA
    ),
    legend.key.width = grid::unit(8, "pt")
  )

fig6a_external_file <- file.path(
  step12_figdir,
  "Figure_06A_External_validation_RF_SVM_ROC_MINIMAL.png"
)

save_plot(
  p_external_roc,
  basename(fig6a_external_file),
  width = 4.7,
  height = 4.8
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

# FINAL MAIN FIGURE 6 — NEW MINIMAL DIRECT PATCHWORK LAYOUT
#            a                          b
#     compact ROC              wide purity scatter
#            c                          d
#   Top-12 forest plot         Immune / Stromal stacked
# No magick::image_extent() or image_append() is used for the MAIN figure.
# Therefore the large blank raster spaces of the old combined image are removed.

p_a_main <- p_external_roc +
  labs(tag = "a") +
  theme(
    plot.margin = ggplot2::margin(6, 8, 6, 8, unit = "pt")
  )

p_b_main <- p_purity_main +
  labs(tag = "b") +
  theme(
    plot.margin = ggplot2::margin(6, 8, 6, 8, unit = "pt")
  )

p_c_main <- p_immune_adjusted +
  labs(tag = "c") +
  theme(
    plot.margin = ggplot2::margin(6, 10, 6, 8, unit = "pt")
  )

p_d_main <- p_estimate_adjusted +
  labs(tag = "d") +
  theme(
    plot.margin = ggplot2::margin(6, 8, 6, 10, unit = "pt")
  )

top_row <- (
  p_a_main +
    p_b_main
) +
  patchwork::plot_layout(
    widths = c(
      0.38,
      0.62
    )
  )

bottom_row <- (
  p_c_main +
    p_d_main
) +
  patchwork::plot_layout(
    widths = c(
      0.47,
      0.53
    )
  )

figure6_minimal <- (
  top_row /
    bottom_row
) +
  patchwork::plot_layout(
    heights = c(
      0.42,
      0.58
    )
  ) &
  theme(
    plot.background = element_rect(
      fill = "white",
      color = NA
    )
  )

figure_06_file <- file.path(
  FIGURE_DIR,
  "Figure_06_MAIN_PURITY_FIRST_MINIMAL_NEW_LAYOUT.png"
)

# Legacy/main-layout PNG exports are intentionally disabled.
# figure6_minimal is still constructed above for reproducibility, but no old PNG
# is written when SAVE_INTERMEDIATE_FIGURES = FALSE.
if (isTRUE(SAVE_INTERMEDIATE_FIGURES)) {
  ggsave(
    filename = figure_06_file,
    plot = figure6_minimal,
    width = 13.2,
    height = 11.8,
    units = "in",
    dpi = FIG_DPI,
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    filename = file.path(
      FIGURE_DIR,
      "Figure_06_MAIN_PURITY_FIRST_MINIMAL_NEW_LAYOUT_WIDE.png"
    ),
    plot = figure6_minimal,
    width = 14.5,
    height = 10.8,
    units = "in",
    dpi = FIG_DPI,
    bg = "white",
    limitsize = FALSE
  )
}

# FINAL COMBINED FIGURE — REFINED BALANCED VERSION
# a = SHORTER WIDE TOP / b = 42% LEFT / c = 58% RIGHT
# FINAL LAYOUT:
#   ┌──────────────────────────────────────────────────────────────┐
#   │                         PANEL A                              │
#   │                    shorter + wide                           │
#   └──────────────────────────────────────────────────────────────┘
#   ┌──────────────────────────┬───────────────────────────────────┐
#   │ B1  Stromal score        │                                   │
#   │                          │                                   │
#   │ B2  Immune score         │            PANEL C                │
#   │                          │       compact forest plot          │
#   │ B3  ESTIMATE score       │                                   │
#   └──────────────────────────┴───────────────────────────────────┘
#          42% width                         58% width
# Only the final composite PNG is saved.

# PANEL A — shorter, wide, compact legend, smaller/more wrapped x labels

supp_a_main <- p_immune_box +
  # Stronger wrapping reduces the horizontal crowding of long immune-cell names.
  scale_x_discrete(
    labels = function(x) stringr::str_wrap(
      as.character(x),
      width = 13
    )
  ) +
  labs(tag = "a") +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      title.position = "left"
    )
  ) +
  theme(
    # FIX: keep the Panel-a legend fully INSIDE the plot canvas.
    # y must stay below 1.00; anchoring at the legend top makes the whole
    # legend extend downward instead of being clipped above the figure.
    legend.position = c(0.50, 0.985),
    legend.direction = "horizontal",
    legend.justification = c(0.5, 1),
    legend.box.just = "center",
    # Compact, non-bold legend for Panel a.
    legend.title = element_text(
      family = FONT_FAMILY,
      face = "plain",
      size = 6.8
    ),
    legend.text = element_text(
      family = FONT_FAMILY,
      size = 6.5
    ),
    legend.key.width = grid::unit(5.5, "pt"),
    legend.key.height = grid::unit(5.5, "pt"),
    legend.spacing.x = grid::unit(1.2, "pt"),
    legend.box.spacing = grid::unit(0, "pt"),
    # No negative margins: negative legend margins can push grobs outside
    # the allocated patchwork cell on some ggplot2/grid versions.
    legend.margin = ggplot2::margin(
      1.5, 4, 1.5, 4,
      unit = "pt"
    ),
    legend.box.margin = ggplot2::margin(
      0, 0, 0, 0,
      unit = "pt"
    ),
    legend.background = ggplot2::element_rect(
      fill = "white",
      color = NA
    ),
    plot.tag = element_text(
      family = FONT_FAMILY,
      face = "plain",
      size = 13.5,
      color = "black"
    ),
    plot.tag.position = c(0.003, 0.995),
    axis.text.x = element_text(
      family = FONT_FAMILY,
      size = 6.7,
      angle = 40,
      hjust = 1,
      vjust = 1,
      lineheight = 0.88,
      color = "black"
    ),
    axis.text.y = element_text(
      family = FONT_FAMILY,
      size = 8.2,
      color = "black"
    ),
    axis.title.y = element_text(
      family = FONT_FAMILY,
      size = 9.5,
      face = "bold"
    ),
    # Give Panel A a small top safety margin so the tag and legend remain
    # completely inside the final exported image.
    plot.margin = ggplot2::margin(
      5, 4, 2, 4,
      unit = "pt"
    )
  )

# PANEL B — three ESTIMATE plots stacked vertically with small spacing

make_vertical_estimate_panel <- function(
    score_label,
    show_x = FALSE,
    show_y_title = FALSE) {

  panel_data <- estimate_long %>%
    dplyr::filter(
      as.character(TME_score) == score_label
    )

  panel_stat <- estimate_panel_stats %>%
    dplyr::filter(
      as.character(TME_score) == score_label
    )

  if (nrow(panel_data) == 0 || nrow(panel_stat) == 0) {
    stop("Could not build vertical Panel B component for: ", score_label)
  }

  local_range <- pmax(
    panel_stat$y_max - panel_stat$y_min,
    1
  )

  p <- ggplot(
    panel_data,
    aes(
      x = MYBL2_group,
      y = Score,
      fill = MYBL2_group
    )
  ) +
    geom_violin(
      trim = FALSE,
      alpha = 0.82,
      color = NA,
      width = 0.88
    ) +
    geom_boxplot(
      width = 0.17,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.50
    ) +
    geom_jitter(
      width = 0.095,
      size = 0.76,
      alpha = 0.67,
      color = "black"
    ) +
    geom_segment(
      data = panel_stat,
      aes(
        x = 1,
        xend = 2,
        y = y_line,
        yend = y_line
      ),
      inherit.aes = FALSE,
      linewidth = 0.36,
      color = "black"
    ) +
    geom_text(
      data = panel_stat,
      aes(
        x = 1.5,
        y = y_text,
        label = fdr_label
      ),
      inherit.aes = FALSE,
      family = FONT_FAMILY,
      size = 2.55,
      color = "black"
    ) +
    facet_wrap(
      ~ TME_score,
      scales = "free_y",
      nrow = 1
    ) +
    scale_fill_manual(
      values = c(
        "MYBL2-low" = COL_MYBL2_LOW,
        "MYBL2-high" = COL_MYBL2_HIGH
      )
    ) +
    labs(
      x = NULL,
      y = if (show_y_title) "ESTIMATE-derived score" else NULL,
      title = NULL
    ) +
    make_clean_theme(
      base_size = 8.8,
      legend_position = "none"
    ) +
    theme(
      strip.text = element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = 8.7,
        color = "black"
      ),
      strip.background = element_rect(
        fill = "white",
        color = "black",
        linewidth = 0.38
      ),
      axis.text.y = element_text(
        family = FONT_FAMILY,
        size = 7.2,
        color = "black"
      ),
      axis.title.y = element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = 8.3,
        color = "black"
      ),
      axis.text.x = element_text(
        family = FONT_FAMILY,
        size = 7.2,
        color = "black"
      ),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.40
      ),
      # Slightly smaller internal margins; spacing is handled by explicit
      # spacer rows between B1/B2/B3 below.
      plot.margin = ggplot2::margin(
        0.5, 2, 0.5, 2,
        unit = "pt"
      )
    ) +
    coord_cartesian(
      ylim = c(
        panel_stat$y_min - 0.06 * local_range,
        panel_stat$y_text + 0.08 * local_range
      ),
      clip = "off"
    )

  if (!show_x) {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  p
}

b_stromal <- make_vertical_estimate_panel(
  score_label = "Stromal score",
  show_x = FALSE,
  show_y_title = FALSE
) +
  labs(tag = "b") +
  theme(
    plot.tag = element_text(
      family = FONT_FAMILY,
      face = "plain",
      size = 13.5,
      color = "black"
    ),
    plot.tag.position = c(0.002, 0.995)
  )

b_immune <- make_vertical_estimate_panel(
  score_label = "Immune score",
  show_x = FALSE,
  show_y_title = TRUE
)

b_estimate <- make_vertical_estimate_panel(
  score_label = "ESTIMATE score",
  show_x = TRUE,
  show_y_title = FALSE
)

# Small explicit gaps make the three B plots visually separate without wasting space.
supp_b_vertical <- (
  b_stromal /
    patchwork::plot_spacer() /
    b_immune /
    patchwork::plot_spacer() /
    b_estimate
) +
  patchwork::plot_layout(
    heights = c(
      1,
      0.035,
      1,
      0.035,
      1
    )
  )

# PANEL C — tighter vertical fill, cell names offset from axis, compact legend

supp_c_main <- p_supp_c_forest +
  labs(
    tag = "c",
    x = "purity-adjusted partial rho (95% ci)"
  ) +
  theme(
    axis.text.x = element_text(
      family = FONT_FAMILY,
      size = 8.8,
      color = "black"
    ),
    axis.text.y = element_text(
      family = FONT_FAMILY,
      size = 9.5,
      color = "black",
      lineheight = 0.93,
      # Move immune-cell names slightly away from the y-axis line.
      margin = ggplot2::margin(
        0, 6, 0, 0,
        unit = "pt"
      )
    ),
    axis.title.x = element_text(
      family = FONT_FAMILY,
      size = 9.3,
      face = "plain",
      color = "black",
      margin = ggplot2::margin(
        6, 0, 3, 0,
        unit = "pt"
      )
    ),
    legend.position = "bottom",
    legend.text = element_text(
      family = FONT_FAMILY,
      size = 6.5
    ),
    legend.key.width = grid::unit(4.8, "pt"),
    legend.key.height = grid::unit(4.8, "pt"),
    legend.spacing.x = grid::unit(0.7, "pt"),
    # Keep a small clean gap between the x-axis title and legend.
    legend.box.spacing = grid::unit(2, "pt"),
    legend.margin = ggplot2::margin(
      0, 0, 0, 0,
      unit = "pt"
    ),
    legend.box.margin = ggplot2::margin(
      0, 0, 0, 0,
      unit = "pt"
    ),
    plot.tag = element_text(
      family = FONT_FAMILY,
      face = "plain",
      size = 13.5,
      color = "black"
    ),
    plot.tag.position = c(0.003, 0.995),
    plot.margin = ggplot2::margin(
      0, 2, 5, 4,
      unit = "pt"
    )
  )

# LOWER ROW — exact requested balance: B = 42%, C = 58%

supp_bottom_row <- (
  supp_b_vertical |
    supp_c_main
) +
  patchwork::plot_layout(
    widths = c(
      0.42,
      0.58
    )
  )

# FINAL COMPOSITION
# Panel a is shorter than before, with a small added gap above the lower b/c row.
# The overall canvas is also shorter, reducing the absolute height of Panel C.

# Add a small explicit vertical gap between Panel a and the lower b/c row.
supp_combined_final <- (
  supp_a_main /
    patchwork::plot_spacer() /
    supp_bottom_row
) +
  patchwork::plot_layout(
    heights = c(
      0.29,
      0.025,
      0.685
    )
  ) &
  theme(
    plot.background = element_rect(
      fill = "white",
      color = NA
    )
  )

# =============================================================================
# THE ONLY FIGURE PNG SAVED BY THIS SCRIPT
# =============================================================================

FINAL_FIGURE_FILE <- file.path(
  FIGURE_DIR,
  "Figure_06_FINAL_REFINED_A_SHORT_B42_C58.png"
)

ggsave(
  filename = FINAL_FIGURE_FILE,
  plot = supp_combined_final,
  width = 11.5,
  height = 11.3,
  units = "in",
  dpi = FIG_DPI,
  bg = "white",
  limitsize = FALSE
)

# Step 9: analysis summary and session information

n_raw_sig <- sum(
  immune_results$Raw_Significant_FDR,
  na.rm = TRUE
)

n_raw_highlight <- sum(
  immune_results$Raw_Highlighted_Moderate_FDR,
  na.rm = TRUE
)

n_purity_sig <- sum(
  immune_results$Purity_Adjusted_Significant_FDR,
  na.rm = TRUE
)

n_sensitivity_sig <- sum(
  immune_results$Purity_Plus_Dataset_Adjusted_Significant_FDR,
  na.rm = TRUE
)

n_persist <- sum(
  immune_results$Raw_Significant_FDR &
    immune_results$Purity_Adjusted_Significant_FDR,
  na.rm = TRUE
)

n_attenuated <- sum(
  immune_results$Raw_Significant_FDR &
    !immune_results$Purity_Adjusted_Significant_FDR,
  na.rm = TRUE
)

summary_lines <- c(
  "MYBL2-associated immune microenvironment analysis completed.",
  "",
  "MAIN FIGURE 6:",
  "a = External validation ROC (RF + SVM).",
  "b = MYBL2 vs ESTIMATE-derived tumor purity.",
  "c = Compact Top-12 display of tumor-purity-adjusted MYBL2–immune-cell partial Spearman correlations (all 28 tested/exported).",
  "d = Main analysis: tumor-purity-adjusted MYBL2 vs Immune/Stromal scores. The purity+dataset sensitivity analysis is retained only as a standalone supplementary output and is not included in the final combined A/B/C figure.",
  "",
  "FINAL COMBINED FIGURE:",
  "Only one final composite PNG is saved: Figure_06_FINAL_REFINED_A_SHORT_B42_C58.png.",
  "Final layout: shorter wide A on top; lower row uses B = 42% and C = 58%; B contains three vertically stacked plots with small gaps; no Panel D is included.",
  "Purity + dataset sensitivity statistics remain in the analysis tables/results; no separate sensitivity PNG is saved in the final-image-only version.",
  "",
  "RAW RESULTS:",
  "Original raw MYBL2-high vs MYBL2-low ssGSEA and ESTIMATE calculations are retained for descriptive analysis, but intermediate PNG files are not saved.",
  "",
  paste0(
    "Pooled tumor samples analysed: ",
    nrow(sample_metadata)
  ),
  paste0(
    "GSE127165 tumor samples: ",
    sum(sample_metadata$Dataset == "GSE127165")
  ),
  paste0(
    "GSE142083 tumor samples: ",
    sum(sample_metadata$Dataset == "GSE142083")
  ),
  "Batch correction: none; pooled TMM-CPM values were analysed without correction.",
  "ESTIMATE input: original non-log TMM-CPM values from the pooled discovery matrix.",
  paste0(
    "Global MYBL2 median on log2(TMM-CPM + 1): ",
    format(
      global_mybl2_median,
      digits = 8
    )
  ),
  paste0(
    "MYBL2-low samples: ",
    sum(sample_metadata$MYBL2_group == "MYBL2-low")
  ),
  paste0(
    "MYBL2-high samples: ",
    sum(sample_metadata$MYBL2_group == "MYBL2-high")
  ),
  paste0(
    "28 official Charoentong immune-signature source: ",
    immune28_source
  ),
  "",
  paste0(
    "MYBL2 vs ESTIMATE-derived tumor purity: Spearman rho = ",
    sprintf(
      "%.4f",
      rho_purity
    ),
    "; p = ",
    format(
      p_purity,
      scientific = TRUE,
      digits = 5
    ),
    "."
  ),
  paste0(
    "Raw MYBL2–ssGSEA associations significant at BH-FDR < 0.05: ",
    n_raw_sig,
    " of ",
    nrow(immune_results),
    "."
  ),
  paste0(
    "Raw correlations additionally showing |rho| > 0.30 (descriptive moderate-magnitude highlight): ",
    n_raw_highlight,
    " of ",
    nrow(immune_results),
    "."
  ),
  paste0(
    "Purity-adjusted MYBL2–ssGSEA associations significant at BH-FDR < 0.05: ",
    n_purity_sig,
    " of ",
    nrow(immune_results),
    "."
  ),
  paste0(
    "Associations significant both before and after purity adjustment: ",
    n_persist,
    "."
  ),
  paste0(
    "Previously significant associations attenuated after purity adjustment: ",
    n_attenuated,
    "."
  ),
  paste0(
    "Sensitivity analysis controlling for purity + dataset: ",
    n_sensitivity_sig,
    " significant associations at BH-FDR < 0.05."
  ),
  "",
  "Interpretation:",
  "The main figure intentionally establishes tumor-purity confounding before presenting immune associations.",
  "The unadjusted negative immune/stromal signals are descriptive and should not be interpreted as generalized immune depletion.",
  "Tumor purity may explain a substantial component of the initially observed negative associations.",
  "ssGSEA and ESTIMATE remain computational estimates and do not establish causal immune suppression by MYBL2."
)

writeLines(
  summary_lines,
  file.path(
    step12_dir,
    "Figure_06_PURITY_FIRST_analysis_summary.txt"
  )
)

write.csv(
  data.frame(
    Panel = c(
      "a",
      "b",
      "c",
      "d"
    ),
    Content = c(
      "External validation ROC (RF + SVM)",
      "MYBL2 vs ESTIMATE-derived tumor purity",
      paste0("Top-", PANEL_CD_N_CELLS, " purity-adjusted partial Spearman forest plot (95% CI; 28 tested)"),
      "Purity-adjusted MYBL2 vs Immune/Stromal scores"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(
    step12_tabledir,
    "Figure_06_MAIN_PURITY_FIRST_Panel_Order.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    step12_dir,
    "Figure_06_PURITY_FIRST_sessionInfo.txt"
  )
)

cat("\n============================================================\n")
cat("MYBL2 IMMUNE MICROENVIRONMENT ANALYSIS FINISHED SUCCESSFULLY\n")
cat("PURITY-ADJUSTED IMMUNE ANALYSIS\n")
cat("============================================================\n")
cat("Main Figure 6 panel order:\n")
cat("a = External validation ROC (RF + SVM)\n")
cat("b = MYBL2 vs ESTIMATE-derived tumor purity\n")
cat("c = Purity-adjusted MYBL2–28 immune-cell correlations\n")
cat("d = Purity-adjusted MYBL2 vs Immune/Stromal scores\n\n")
cat("Raw high/low immune results were retained as supplementary outputs.\n\n")
cat(
  "MYBL2 vs purity: rho = ",
  sprintf("%.4f", rho_purity),
  "; p = ",
  format(
    p_purity,
    scientific = TRUE,
    digits = 5
  ),
  "\n",
  sep = ""
)
cat(
  "Purity-adjusted significant immune associations: ",
  n_purity_sig,
  " / ",
  nrow(immune_results),
  "\n",
  sep = ""
)
cat(
  "Purity + dataset sensitivity significant associations: ",
  n_sensitivity_sig,
  " / ",
  nrow(immune_results),
  "\n\n",
  sep = ""
)
cat("Pooled input CPM file:\n", discovery_file, "\n", sep = "")
cat("Output directory (all files):\n", immune_output_dir, "\n", sep = "")
cat("Final combined Figure 6:\n", figure_06_file, "\n", sep = "")
cat("============================================================\n")
