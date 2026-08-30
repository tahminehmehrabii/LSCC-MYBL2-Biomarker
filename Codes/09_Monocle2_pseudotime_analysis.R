# LSCC Monocle 2 pseudotime analysis
# Performs DDRTree reduction, graph-distance pseudotime, and MYBL2 trajectory analysis.

rm(list = ls(all.names = TRUE))
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 100,
  timeout = 7200,
  width = 160
)
set.seed(123)

# Step 1: version check and paths

if (as.character(getRversion()) != "4.4.3") {
  stop(
    "Run this script only in R 4.4.3. Current R version: ",
    as.character(getRversion())
  )
}

# Analysis output directory
PSEUDOTIME_MONOCLE_DIR <- "D:/LSCC/Results_LSCC/Pseudotime_Monocle"

# Dedicated Monocle 2.34.0 library used only in the fresh trajectory/figure RGui
# session. This is NOT used in the preparation session.
MONOCLE_LIB <- "D:/LSCC/Monocle2_R44_library"

# Prepared input created by Script 08
INPUT_RDS <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "Monocle2_full_analysis_input.rds"
)

# Monocle trajectory-stage output files.
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

# Figure filename stems. Every one is saved as PNG, TIFF and PDF.
FIGURE_3A_STEM <- "Figure_3A_HighCNV_DDRTree_pseudotime"
FIGURE_3B_STEM <- "Figure_3B_MYBL2_across_pseudotime"
FIGURE_8_STEM  <- "Figure_8_Combined_Pseudotime_MYBL2"

dir.create(
  PSEUDOTIME_MONOCLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Step 2: analysis settings

TARGET_GENE <- "MYBL2"

# Balances cell representation across each sample–subcluster combination.
MAX_CELLS_PER_SAMPLE_SUBCLUSTER <- 350L

# Ordering genes come from precomputed strict markers.
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

# The Monocle state with the LOWEST median score for these G2/M genes becomes
# the root state.
G2M_ROOT_GENE_SET <- c(
  "CDK1", "CCNB1", "CCNB2", "CDC20", "CDC25C", "CENPA", "CENPE",
  "CENPF", "CENPM", "CENPU", "CKS1B", "CKS2", "MKI67", "NUSAP1",
  "PLK1", "PRC1", "TOP2A", "TPX2", "UBE2C", "BIRC5"
)

LOESS_SPAN <- c(
  "Cluster 0" = 0.90,
  "Cluster 8" = 1.00
)

# Approximate pointwise confidence interval shown around each LOESS fit.
LOESS_CI_LEVEL <- 0.95

# Confidence ribbons are shown only within the central, data-supported portion
# of each cluster's pseudotime distribution. This prevents visually unstable
# endpoint ribbons driven by sparse observations.
LOESS_RIBBON_Q_LOWER <- 0.05
LOESS_RIBBON_Q_UPPER <- 0.95
LOESS_RIBBON_MIN_LOCAL_CELLS <- 20L
LOESS_RIBBON_WINDOW_FRACTION <- 0.12

# Step 3: standard manuscript figure settings

# Font and export quality.
FONT_FAMILY <- "Arial"
FIG_DPI <- 600L
FIG_BACKGROUND <- "white"

# Standard manuscript dimensions in inches.
FIG_SINGLE_W <- 3.50
FIG_SINGLE_H <- 4.20
FIG_DOUBLE_W <- 7.20
FIG_DOUBLE_H <- 5.40

# Three-panel pseudotime figures. Extra width is reserved for side legends beside each trajectory panel.
FIGURE3_COMBINED_W <- 14.00
FIGURE3_COMBINED_H <- 7.20

# Five-panel combined pseudotime/MYBL2 figure.
# Extra width reserves clean side legends beside panels a-c.
FIGURE8_COMBINED_W <- 14.00
FIGURE8_COMBINED_H <- 12.40

# Typography in points.
BASE_TEXT_PT <- 9.5
AXIS_TEXT_PT <- 8.5
AXIS_TITLE_PT <- 10
LEGEND_TEXT_PT <- 8.5
LEGEND_TITLE_PT <- 9
PLOT_TITLE_PT <- 10
PANEL_TAG_PT <- 14

# Line widths.
PANEL_BORDER_LWD <- 0.45
AXIS_LWD <- 0.40
GEOM_LWD <- 0.55
GRID_LWD <- 0.30

# Heatmap typography retained for consistency with other manuscript figures.
HEATMAP_FONT_PT <- 8
HEATMAP_ROW_FONT_PT <- 7
HEATMAP_COL_FONT_PT <- 7
HEATMAP_LEGEND_FONT_PT <- 8

# Repeated manuscript semantic colours.
COL_NORMAL <- "#7470B2"
COL_TUMOR <- "#D81B60"
COL_MYBL2_LOW <- COL_NORMAL
COL_MYBL2_HIGH <- COL_TUMOR
COL_UP <- COL_TUMOR
COL_DOWN <- COL_NORMAL
COL_NS <- "#B8B8B8"

# Retained for future CellChat panels.
CELLCHAT_PALETTE <- grDevices::colorRampPalette(
  c("#E8F6F4", "#B8E1DA", "#73C6BE", "#2D9D8C", "#006E63")
)(100)

# Pseudotime and state palettes.
PSEUDOTIME_COLOURS <- c(
  "#132B43", "#2166AC", "#4EB3D3",
  "#A1DAB4", "#FEE08B", "#D73027"
)

# Okabe-Ito colourblind-friendly palette. State 1 and State 6 are
# deliberately blue and vermilion, avoiding the former orange-on-orange pairing.
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

# Step 4: general helper functions

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

# Step 5: Monocle and figure helper functions

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

  # Monocle usually stores coordinates as components x cells.
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

  # A confidence ribbon is only informative where the local pseudotime density
  # is adequate. We therefore trim the plotted ribbon to the central 5–95%
  # pseudotime interval and retain only grid locations with sufficient nearby cells.
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

  # Separate non-contiguous valid intervals so ggplot never bridges a sparse gap.
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

  # Convert a graph and a coordinate table into drawable edges. The helper first
  # matches vertices by name and only uses vertex order when the graph and
  # coordinate table have the exact same number of nodes.
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

    # In some Monocle objects, the principal-node graph has generic vertex
    # names while the principal-coordinate matrix has no row names. In that
    # case, a one-to-one vertex-order mapping is valid only when sizes match.
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

  # -------------------------------------------------------------------------
  # First preference: Monocle's projected DDRTree graph. This is the graph
  # used by Monocle for cell projection and typically aligns directly to cell
  # barcodes and reducedDimS coordinates.
  # -------------------------------------------------------------------------

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

  # -------------------------------------------------------------------------
  # Second preference: principal-node coordinates plus the minimum-spanning
  # tree. These are the compact DDRTree backbone nodes rather than cell points.
  # -------------------------------------------------------------------------

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

  # -------------------------------------------------------------------------
  # Third preference: minSpanningTree mapped directly to reduced cell
  # coordinates. This preserves the learned graph when only a cell-level MST
  # was stored in the Monocle object.
  # -------------------------------------------------------------------------

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
# Step 6: igraph compatibility for Monocle 2
# Monocle 2 was written against older igraph releases that used the argument
# `neimode` in dfs(). igraph >= 1.3 renamed this argument to `mode`; in igraph
# 2.x, `neimode` is defunct. This compatibility wrapper translates the legacy
# argument before Monocle calls the current igraph dfs() implementation.

install_igraph_legacy_monocle_compatibility_patch <- function() {
  # Monocle 2.34.0 contains a small number of calls written for old igraph
  # versions: bare dfs(..., neimode = ...) and nei(...).  A namespace is a
  # locked environment in R, therefore it is unsafe to add new bindings such
  # as "dfs" or ".nei" to the Monocle namespace.  Instead, this patch rewrites
  # ONLY existing Monocle functions that contain those legacy calls.
  #
  # This is session-only: installed package files are never changed.

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
    # Never try to create a new name in a locked namespace.  That was the
    # cause of: "cannot add bindings to a locked environment".
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

    # igraph >= 2.1 removed bare nei().  Use the package-internal helper
    # explicitly, avoiding the need to create .nei in Monocle's namespace.
    updated_body_text <- gsub(
      "(?<![[:alnum:]_.:])nei\\s*\\(",
      "igraph:::.nei(",
      updated_body_text,
      perl = TRUE
    )

    # igraph >= 1.3 renamed the dfs argument.
    updated_body_text <- gsub(
      "\\bneimode\\s*=",
      "mode =",
      updated_body_text,
      perl = TRUE
    )

    # Resolve legacy bare dfs() calls directly through igraph rather than
    # attempting to create an imported dfs binding inside Monocle.
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

# Step 7: load Monocle packages

load_monocle_environment <- function() {
  # In auto-rebuild mode, Seurat/SeuratObject may already be loaded only to
  # reconstruct the missing prepared input. The Monocle analysis below uses
  # explicit namespaces and the compatibility dispersion model, so this is safe.
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

  # igraph 2.0.3 was the originally tested version. Newer igraph 2.x
  # releases retain the functions used in this workflow, so do not block a
  # valid current installation solely because its patch/minor version differs.
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

# Step 8: figure settings and DDRTree backbone extraction

# Extra width gives every top-row panel enough room for its own side legend.
FIGURE3_COMBINED_W <- 17.40
# Keep the three-panel DDRTree figure tall enough for readable side legends.
FIGURE3A_COMBINED_H <- 6.60
# MYBL2 has only two wide panels; a lower height removes unused vertical space
# without excluding any cells or changing the fitted curves.
FIGURE3B_COMBINED_H <- 4.45
FIGURE8_COMBINED_W <- 17.40
# The combined figure is exported taller so it occupies roughly half a Word page more comfortably.
FIGURE8_COMBINED_H <- 11.20

# DDRTree backbone extraction
# This function intentionally NEVER draws a projected-cell graph or a cell-level
# MST. Those graphs contain one edge per cell and produce the unwanted zig-zag
# appearance. Only compact DDRTree principal nodes are used. When an exact stored
# principal-node tree is unavailable, a compact MST is reconstructed solely from
# the principal-node coordinates; no individual cells are connected by lines.
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

    # DDRTree principal coordinates are commonly 2 x K; transpose to K x 2.
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

    # Use explicit node labels if available; otherwise use matched vertex order
    # only when graph and principal-node coordinate counts are identical.
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

    # Remove duplicate undirected edges.
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

  # Guard against accidentally treating a cell-level coordinate matrix as a principal-node matrix.
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

  # First use an exact stored tree only when its vertices match the compact principal nodes.
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

  # When Monocle does not retain the principal-node tree separately, reconstruct
  # a compact minimum-spanning tree only among principal nodes. This keeps the
  # backbone smooth and branch-level; it never links individual cells.
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

  # ---------------------------------------------------------------------------
  # Extract a visible DDRTree / MST backbone for all trajectory panels.
  # This does not alter trajectory inference; it is only a graphical overlay.
  # ---------------------------------------------------------------------------

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

  # ---------------------------------------------------------------------------
  # Identify and mark the pseudotime start within the selected root state.
  # The selected point is the lowest-pseudotime cell in the G2/M-defined root state.
  # ---------------------------------------------------------------------------

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

  # Compact root annotation: smaller, closer to the root cell, and kept
  # inside the panel so it informs the reader without dominating panel a.
  x_span <- diff(trajectory_x_limits)
  y_span <- diff(trajectory_y_limits)
  x_mid <- mean(trajectory_x_limits)
  y_mid <- mean(trajectory_y_limits)

  x_direction <- if (root_point$DDRTree_1 >= x_mid) -1 else 1
  y_direction <- if (root_point$DDRTree_2 >= y_mid) -1 else 1

  root_label_x <- root_point$DDRTree_1 + x_direction * max(0.08 * x_span, 0.32)
  root_label_y <- root_point$DDRTree_2 + y_direction * max(0.07 * y_span, 0.28)

  # Keep the compact label safely within the panel interior.
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

  # Each top-row panel carries its own compact legend at the right side.
  # Legends are deliberately not collected below the panels.
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

  # ---------------------------------------------------------------------------
  # Figure 3A: Pseudotime, Monocle state, selected High-CNV subclusters
  # ---------------------------------------------------------------------------

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

  # Independent side legends remain beside their corresponding trajectory panels.
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

  # ---------------------------------------------------------------------------
  # MYBL2 data and Figure 3B
  # ---------------------------------------------------------------------------

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

  # Each MYBL2 panel receives its own data-supported pseudotime range.
  # This removes blank space while retaining the original numeric pseudotime scale.

  # Scale the y-axis using actual observations and only confidence intervals
  # that passed the data-support filter. Sparse, trimmed endpoints cannot inflate
  # the plotting range.
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

    # Use the observed pseudotime span of this specific cluster, with a small
    # proportional padding. This avoids unused horizontal space, especially
    # for Cluster 11, while preserving true pseudotime values and tick labels.
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

  # ---------------------------------------------------------------------------
  # Figure 8: five-panel combined final manuscript figure
  # ---------------------------------------------------------------------------

  # A compact lower row prevents excess blank vertical space in panels d and e.
  # All observed cells and the full LOESS fits remain plotted; only panel height changes.
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

# Step 9: DDRTree graph-distance pseudotime
# Rationale:
# Monocle 2.34.0's orderCells()/extract_ddrtree_ordering() was written against
# old igraph vertex-sequence semantics. In this R 4.4.3 / igraph 2.0.3 session,
# it still stops with “incorrect number of dimensions” even after project2MST()
# has been repaired. DDRTree reduction itself completed successfully.
# This final runner keeps Monocle's DDRTree embedding and principal MST, then
# calculates graph pseudotime directly from the learned principal tree. It does
# NOT call monocle::orderCells(), project2MST(), or extract_ddrtree_ordering().
# The graph-distance pseudotime is therefore stable across current igraph APIs.

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

  # A non-branching tree is one trajectory segment.
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

  # Assign each branch node to its closest already-labelled neighbour in graph
  # distance from the selected root; this makes every principal node usable.
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

  # Defensive final fill, normally unnecessary for a connected tree.
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

  # Assign every cell to its nearest principal DDRTree node.
  z_norm <- colSums(Z * Z)
  y_norm <- colSums(Y * Y)
  squared_distance <- outer(y_norm, z_norm, "+") - 2 * crossprod(Y, Z)
  squared_distance[squared_distance < 0 & squared_distance > -1e-10] <- 0
  closest_node_index <- apply(squared_distance, 2L, which.min)
  closest_node_name <- principal_names[closest_node_index]
  cell_to_node_distance <- sqrt(pmax(
    squared_distance[cbind(closest_node_index, seq_len(n_cells))], 0
  ))

  # Initial structural states are independent of root selection. Use temporary
  # root distances from the first node only for deterministic branch-node labels.
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

  # Select the graph node inside the chosen state with the lowest median root
  # score among its assigned cells. Nodes without cells are not candidates.
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

  # Recalculate structural states with the final root for deterministic labels,
  # then retain the state selected by the specified lowest-median-G2/M rule.
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

  # Keep the learned principal graph in the CDS and add transparent metadata.
  monocle::minSpanningTree(cds) <- tree

  list(
    cds = cds,
    root_state = as.integer(root_state),
    root_node = root_node,
    state_medians = state_medians,
    node_table = node_table
  )
}

# Final replacement for the full trajectory runner. The auto-discovery wrapper
# defined above calls this object when the prepared input exists.
# Step 10: run trajectory analysis and generate figures

run_trajectory_and_figures <- function() {
  cat(
    "\n============================================================\n",
    "MONOCLE 2 DDRTREE AND GRAPH-DISTANCE PSEUDOTIME\n",
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

# Execute analysis
run_trajectory_and_figures()
