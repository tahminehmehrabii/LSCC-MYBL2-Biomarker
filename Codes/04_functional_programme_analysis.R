# LSCC high-CNV malignant-cell functional-programme analysis
# Scores functional programmes and generates Figure 7.

# Step 1: setup, paths, and figure settings

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

# Standard manuscript figure settings

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

# Larger final Figure 7 canvas: actual panel content is rendered larger.
FIGURE7_COMBINED_W <- 16.50
FIGURE7_COMBINED_H <- 12.20

# Real panel-row allocations. These control how much of the enlarged canvas
# is assigned to the top-row panels (a/b) and bottom-row panel (c).
FIGURE7_TOP_ROW_HEIGHT <- 1.15
# Shorter bottom row removes excess empty framed space around panel c.
FIGURE7_BOTTOM_ROW_HEIGHT <- 0.82

BASE_TEXT_PT <- 9.5
AXIS_TEXT_PT <- 8.5
AXIS_TITLE_PT <- 10
LEGEND_TEXT_PT <- 8.5
LEGEND_TITLE_PT <- 9
PLOT_TITLE_PT <- 10
PANEL_TAG_PT <- 14

# Shared visual settings used by every Figure 7 panel.
# PANEL_BORDER_LWD frames the data region; PANEL_SECTION_BORDER_LWD frames
# each complete a/b/c section, including its tag and legend.
PANEL_BORDER_LWD <- 0.45
PANEL_SECTION_BORDER_LWD <- 0.45
AXIS_LWD <- 0.40
TICK_LWD <- 0.40
GEOM_LWD <- 0.55
GRID_LWD <- 0.30
TICK_LENGTH_PT <- 2
# Smaller inter-component gap leaves more area for actual panel content.
PANEL_SPACING_PT <- 3
PANEL_TAG_POSITION <- c(0.045, 0.955)

# Vertical DotPlot settings: one MYBL2 column across all ten subtypes.
# This replaces the empty single horizontal dot row in panel b.
VERTICAL_DOTPLOT_X_LIMITS <- c(0.55, 1.45)
VERTICAL_DOTPLOT_POINT_RANGE <- c(3.0, 14.0)

# Panel-b-only clearance: reserves a clean top-left area for tag b, away from
# the first y-axis subtype label (Cluster 0).
PANEL_B_TAG_POSITION <- c(0.018, 0.988)
PANEL_B_TAG_MARGIN <- ggplot2::margin(
  t = 13,
  r = 4,
  b = 3,
  l = 19,
  unit = "pt"
)

# Panel-a-only clearance: reserves a clean top-left area for tag a, away from
# the t-SNE embedding and its nearby point labels.
PANEL_A_TAG_POSITION <- c(0.018, 0.988)
PANEL_A_TAG_MARGIN <- ggplot2::margin(
  t = 13,
  r = 4,
  b = 3,
  l = 19,
  unit = "pt"
)

HEATMAP_FONT_PT <- 8
# Heatmap axes and legends follow the same text hierarchy as the other panels.
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

# Reduced shared margin: Figure 7 equivalent of decreasing a raster inset.
# It gives more of every framed section to its actual plot, axes and legend.
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

      # All panel titles use the same font, bold state, size, and placement.
      plot.title = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PLOT_TITLE_PT,
        hjust = 0.5,
        vjust = 0.5,
        color = "black"
      ),
      plot.title.position = "panel",

      # All panel tags use the same Arial-bold style and top-left inset position.
      plot.tag = ggplot2::element_text(
        family = FONT_FAMILY,
        face = "bold",
        size = PANEL_TAG_PT,
        hjust = 0,
        vjust = 1,
        color = "black"
      ),
      plot.tag.position = PANEL_TAG_POSITION,

      # All axis lines, ticks, and number labels follow the same standard.
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

      # All legends use identical typography, key size, and spacing.
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

      # Each complete panel section receives an outer black frame, while the
      # data area keeps its inner thin black frame. The outer frame encloses
      # the panel tag, axes, plot, titles, and legends as one clear section.
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

      # Equal internal panel margin and distance between panels.
      panel.spacing = grid::unit(PANEL_SPACING_PT, "pt"),
      plot.margin = MANUSCRIPT_MARGIN
    )
}

# Figure 7 layout settings.
# Panel b receives more width because the MYBL2 DotPlot is now vertical;
# panel a remains large enough for a readable t-SNE embedding.
TSNE_PANEL_WIDTH <- 1.20
MYBL2_PANEL_WIDTH <- 0.80

# Panel c is intentionally narrower than the full figure and centered between
# two empty spacers, so the heatmap no longer looks excessively stretched.
PANEL_C_LEFT_SPACER_WIDTH <- 0.35
PANEL_C_HEATMAP_WIDTH <- 1.30
PANEL_C_RIGHT_SPACER_WIDTH <- 0.35

# Panel-c-only clearance: reserves a clean top-left area for tag c, away from
# the first y-axis programme label ("Stress response").
PANEL_C_TAG_POSITION <- c(0.018, 0.988)
PANEL_C_TAG_MARGIN <- ggplot2::margin(
  t = 13,
  r = 4,
  b = 3,
  l = 19,
  unit = "pt"
)

# Existing High-CNV clustering is preserved exactly as previously defined.
RESOLUTION_MALIGNANT_FROM_PREVIOUS_PIPELINE <- 0.5

FIXED_CLUSTER_ORDER <- c(
  "Cluster 0", "Cluster 1", "Cluster 2", "Cluster 3",
  "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 8",
  "Cluster 10", "Cluster 11"
)

# Panel c is intentionally restricted to these two High-CNV malignant subclusters.
# Panels a and b continue to display all ten existing subclusters.
PANEL_C_CLUSTER_ORDER <- c("Cluster 0", "Cluster 8")

# Publication-ready labels for the functional-programme heatmap.
# The order is defined from bottom to top on the y-axis.
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

# Delete only prior versions of the functional Figure 7.
CLEAN_OLD_FUNCTIONAL_FIGURES <- TRUE

if (CLEAN_OLD_FUNCTIONAL_FIGURES) {
  old_functional_figure_files <- list.files(
    OUTPUT_DIR,
    pattern = "^Figure_07_HighCNV_MYBL2_functional_programmes[.](png|tiff|pdf)$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(old_functional_figure_files) > 0) {
    unlink(old_functional_figure_files, force = TRUE)
  }
}

# Optional direct paths. Keep NA unless automatic detection fails.
HIGHCNV_OBJECT_RDS_OVERRIDE <- NA_character_
HIGHCNV_DE_OBJECT_RDS_OVERRIDE <- NA_character_

mark_step_done("01.1")

# Step 2: detect completed prior-pipeline rds files

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

# Step 3: package check and loading

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

# Step 4: helper functions

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

  # PNG: 600 dpi
  ggplot2::ggsave(
    filename = png_file,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = FIG_DPI,
    bg = FIG_BACKGROUND,
    limitsize = FALSE
  )

  # TIFF: 600 dpi with LZW compression
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

  # PDF: vector version
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

# Step 5: load and validate existing high-cnv objects

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

# Use the original High-CNV object exclusively for the pre-existing embedding,
# while using the DE object for expression-based downstream analyses.
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

# Step 6: functional programmes and figure 7

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

# Figure 7a: Existing High-CNV subtype embedding.
# The automatic technical title generated by DimPlot is explicitly removed.
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

# Figure 7b: MYBL2 DotPlot. Seurat::DotPlot calculates average expression and
# percent expressed; it is redrawn here for the manuscript layout.
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
  # One MYBL2 column is displayed vertically across all ten subtypes, using
  # the full height of panel b rather than one short horizontal row of dots.
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
    # Compact legends at the right keep panel b tall and prevent a wide,
    # visually heavy legend strip below the DotPlot.
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

# Figure 7c: descriptive functional-programme comparison for Cluster 0 and
# Cluster 8 only. The programme Z-score scale remains based on all ten
# subclusters so colours remain comparable with the full analysis.
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
    "Panel c could not retain the expected subclusters: ",
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

# Figure 7 has no in-image title. Panel letters are placed directly inside their
# own panels, so a, b, and c share the same font, size, and top-left inset.
# plot.background adds an outer thin black frame to each complete panel section,
# creating visible divisions between panels a, b, and c.
# Panel a alone gets a dedicated top-left clearance zone so tag a remains
# separate from the t-SNE embedding and the nearby cluster labels.
p_subclusters_tagged <- p_subclusters +
  ggplot2::labs(tag = "a") +
  ggplot2::theme(
    plot.tag.position = PANEL_A_TAG_POSITION,
    plot.margin = PANEL_A_TAG_MARGIN
  )

# Panel b alone gets a dedicated top-left clearance zone so tag b does not
# overlap the Cluster 0 y-axis label.
p_mybl2_only_tagged <- p_mybl2_only +
  ggplot2::labs(tag = "b") +
  ggplot2::theme(
    plot.tag.position = PANEL_B_TAG_POSITION,
    plot.margin = PANEL_B_TAG_MARGIN
  )

# Panel c alone gets a larger top/left margin and a dedicated tag position.
# This prevents the c tag from overlapping the top programme label.
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

# Panel c is centered and deliberately narrower than the full figure.
# The spacers have no border; only the heatmap itself remains framed as panel c.
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
  filename_stem = "Figure_07_HighCNV_MYBL2_functional_programmes",
  dir_path = OUTPUT_DIR,
  width = FIGURE7_COMBINED_W,
  height = FIGURE7_COMBINED_H
)

mark_step_done("01.6")

# Step 7: export objects and final status

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
  "Figure 7a: High-CNV malignant-subcluster embedding.",
  "Figure 7b: MYBL2 DotPlot across all High-CNV subclusters.",
  "Figure 7c: descriptive functional-programme heatmap for Cluster 0 and Cluster 8 only.",
  "Functional programmes were gene-wise Z-scored before averaging and are descriptive, not inferential statistical comparisons.",
  "Custom programme gene sets: Stress_response, Senescence_like, and Epithelial_differentiation; exact genes are exported separately."
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
cat("Figure 7: Figure_07_HighCNV_MYBL2_functional_programmes.{png,tiff,pdf}\n")
cat("============================================================\n")

mark_step_done("01.7")

message(
  "All numbered steps finished. Outputs are in: ",
  OUTPUT_DIR,
  "\nFigure 7 includes all ten High-CNV malignant subtypes in panels a-b; panel c includes Cluster 0 and Cluster 8 only."
)
