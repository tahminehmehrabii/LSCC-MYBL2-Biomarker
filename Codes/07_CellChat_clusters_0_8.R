rm(list = ls(all.names = TRUE))
invisible(gc())

if (!exists("project_path")) {
  project_path <- getwd()
}

setwd(project_path)

set.seed(123)
options(
  stringsAsFactors = FALSE,
  scipen = 100,
  timeout = 7200,
  width = 160
)

# STEP 01: Set paths and analysis parameters.

BASE_RESULTS_DIR <- file.path("Results")
SC_PROJECT_DIR <- file.path("Results", "scRNAseq", "GSE206332")
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

# STEP 02: Define common figure settings.

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

# STEP 03: Load required R packages.

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

# STEP 04: Define helper functions.

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

# STEP 05: Load or build the CellChat object.

cellchat_rds_candidates <- unique(c(
  file.path(BASE_RESULTS_DIR, "CellChat", "cellchat_c0_c8.rds"),
  file.path(BASE_RESULTS_DIR, "CellChat", "CellChat_object.rds"),
  file.path(BASE_RESULTS_DIR, "CellChat", "LSCC_Stage02_selected_HighCNV_CellChat.rds"),
  file.path(BASE_RESULTS_DIR, "LSCC_02_HighCNV_CellChat", "LSCC_Stage02_selected_HighCNV_CellChat.rds"),
  file.path(BASE_RESULTS_DIR, "CellChat", "GSE206332_all_annotated_cells_existing_resolution0.5_HighCNV_CellChat.rds"),
  file.path(BASE_RESULTS_DIR, "CellChat", "GSE206332_all_annotated_cells_existing_resolution0.5_HighCNV_CellChat_core.rds"),
  list.files(BASE_RESULTS_DIR, pattern = "CellChat.*\\.rds$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
))

cellchat_rds_candidates <- cellchat_rds_candidates[file.exists(cellchat_rds_candidates)]

cellchat <- NULL
if (!FORCE_RECOMPUTE_CELLCHAT && length(cellchat_rds_candidates) > 0L) {
  for (candidate in cellchat_rds_candidates) {
    candidate_obj <- safe_read_rds(candidate)
    
    if (is_compatible_cellchat(candidate_obj)) {
      cellchat <- candidate_obj
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
    file.path(OUTPUT_DIR, "cc_counts_build.csv"),
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
}

STABLE_CELLCHAT_RDS <- file.path(OUTPUT_DIR, "cellchat_c0_c8.rds")
saveRDS(cellchat, STABLE_CELLCHAT_RDS)

# STEP 06: Validate CellChat results and export communication tables.

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
  file.path(OUTPUT_DIR, "cc_counts.csv"),
  row.names = FALSE
)

all_communications <- tryCatch(
  CellChat::subsetCommunication(cellchat),
  error = function(e) data.frame()
)

if (nrow(all_communications) > 0L) {
  write.csv(
    all_communications,
    file.path(OUTPUT_DIR, "cc_all_comm.csv"),
    row.names = FALSE
  )
}

# STEP 07: Prepare network weights and ligand-receptor tables.

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
  file.path(OUTPUT_DIR, "c0_c8_weights.csv"),
  row.names = FALSE
)

write.csv(
  bubble_table,
  file.path(OUTPUT_DIR, "c0_c8_lr_top5.csv"),
  row.names = FALSE
)

probability_limits <- safe_limits(bubble_table$prob)
significance_limits <- safe_limits(bubble_table$Minus_log10_P)
edge_limits <- safe_limits(network_table$Weight)

# STEP 08: Define plotting functions.

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

# STEP 09: Build and save Figure 13.

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

save_plot_all_formats(
  plot_obj = figure_13_combined,
  filename_stem = "fig13_cellchat_c0_c8",
  dir_path = OUTPUT_DIR,
  width = FIGURE13_COMBINED_W,
  height = FIGURE13_COMBINED_H
)
