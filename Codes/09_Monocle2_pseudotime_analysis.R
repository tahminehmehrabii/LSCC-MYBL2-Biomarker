# STEP 01: Set paths and parameters.
# STEP 02: Define helper functions.
# STEP 03: Load Monocle environment.
# STEP 04: Build trajectory and figures.
# STEP 05: Run analysis.

rm(list = ls(all.names = TRUE))
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 100,
  timeout = 7200,
  width = 160
)
set.seed(123)


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
  "monocle_input.rds"
)

CELL_BALANCE_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "cell_balance.csv"
)
ORDERING_GENES_INPUT_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "ordering_genes.csv"
)
PREP_PROVENANCE_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "input_provenance.txt"
)
PREP_STATUS_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "input_status.txt"
)
PREP_SESSION_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "input_session.txt"
)

STATUS_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "status.txt"
)
ERROR_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "error.txt"
)
RUN_LOG <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "run.log"
)
SESSION_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "session.txt"
)
COMPATIBILITY_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "compat.txt"
)
REDUCTION_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "reduction.txt"
)
CDS_RDS <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "monocle_cds.rds"
)
TRAJECTORY_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "pseudotime_meta.csv"
)
ORDERING_GENES_TRAJECTORY_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "ordering_genes_used.csv"
)
ROOT_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "root_state.csv"
)
ROOT_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "root_rule.txt"
)
MYBL2_CELL_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "mybl2_cells.csv"
)
LOESS_CURVE_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "mybl2_loess.csv"
)
BACKBONE_STATUS_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "backbone_status.txt"
)
BACKBONE_SEGMENTS_CSV <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "backbone_segments.csv"
)

FIGURE_3A_STEM <- "fig03a_pseudotime"
FIGURE_3B_STEM <- "fig03b_mybl2"
FIGURE_8_STEM  <- "fig08_pseudotime_mybl2"

dir.create(
  PSEUDOTIME_MONOCLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


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

log_message <- function(...) invisible(NULL)


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
  } else {
  }
  
  invisible(length(patched_functions) > 0L)
}


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
  
  
  install_igraph_legacy_monocle_compatibility_patch()
  
  invisible(TRUE)
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


run_trajectory_and_figures <- function() {
  
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
  write_status <- function() invisible(NULL)
  
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
      
      
      log_message("TRAJECTORY SUCCESSFUL | Root state = ", root_state)
    },
    error = function(e) {
      result$error <<- conditionMessage(e)
      result$stage <<- paste0("failed_at_", result$stage)
      
      
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
  
  
  if (!isTRUE(result$trajectory_success) || !isTRUE(result$figure_success)) {
    warning(
      "The run did not fully complete. Read these files inside Pseudotime_Monocle:\n",
      "error.txt\n",
      "status.txt\n",
      "run.log"
    )
  }
}


FIGURE_REFRESH_STATUS_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "figure_refresh.txt"
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
      "monocle_input.RDS"
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
  
  if (!file.exists(CDS_RDS)) {
    stop("Completed CellDataSet is missing:\n", CDS_RDS)
  }
  
  load_monocle_environment()
  
  cds <- readRDS(CDS_RDS)
  
  if (!inherits(cds, "CellDataSet")) {
    stop("Saved monocle_cds.rds is not a valid Monocle CellDataSet.")
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
  
  
  
  invisible(TRUE)
}


AUTO_REBUILD_STATUS_TXT <- file.path(
  PSEUDOTIME_MONOCLE_DIR,
  "auto_rebuild.txt"
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
  
  
  rm(
    highcnv_obj, highcnv_de_obj, counts_all, counts, metadata,
    phenotype_data, markers, marker_data
  )
  gc()
  
  
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
    
    
    return(original_run_trajectory_and_figures())
  }
  
  
  rebuild_prepared_input_from_real_sources()
  
  if (!file.exists(INPUT_RDS)) {
    stop(
      "Auto-rebuild finished without creating the required input RDS:\n",
      INPUT_RDS
    )
  }
  
  original_run_trajectory_and_figures()
}


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
  
  
  result <- list(
    trajectory_success = FALSE, figure_success = FALSE,
    reduction_fallback_used = FALSE, error = NA_character_,
    figure_error = NA_character_, stage = "started", root_state = NA_integer_,
    root_node = NA_character_, n_cells = 0L, n_ordering_genes = 0L,
    monocle_version = as.character(utils::packageVersion("monocle")),
    igraph_version = as.character(utils::packageVersion("igraph"))
  )
  write_status <- function() invisible(NULL)
  
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
    
    saveRDS(cds, CDS_RDS)
    result$trajectory_success <- TRUE
    result$stage <- "trajectory_completed"
    log_message("TRAJECTORY SUCCESSFUL | Root state = ", result$root_state, " | Root node = ", result$root_node)
  }, error = function(e) {
    result$error <<- conditionMessage(e)
    result$stage <<- paste0("failed_at_", result$stage)
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
  
  if (!isTRUE(result$trajectory_success) || !isTRUE(result$figure_success)) {
    warning("The run did not fully complete. Read error.txt, status.txt and run.log in Pseudotime_Monocle.")
  }
  invisible(result)
}

run_trajectory_and_figures()
