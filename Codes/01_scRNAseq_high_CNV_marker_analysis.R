# LSCC single-cell RNA-seq, inferCNV, and high-CNV marker analysis
# Dataset: GSE206332

# Step 1: clean environment

rm(list = ls())
gc()

set.seed(123)
options(stringsAsFactors = FALSE)
options(scipen = 100)

# Step 2: libraries

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

# Step 3: paths

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

# Step 4: parameters

# QC
QC_MIN_FEATURES <- 200
QC_MAX_MT <- 20

# Dimensionality reduction / clustering
N_HVG <- 2000
NPCS <- 40
HARMONY_DIMS <- 1:30
RESOLUTION_GLOBAL <- 0.8
RESOLUTION_MALIGNANT <- 0.5

# Marker detection
FINDALL_LOGFC <- 0.10
FINDALL_MINPCT <- 0.10
GLOBAL_MARKER_SIG_LOGFC <- 0.25
GLOBAL_MARKER_SIG_PADJ <- 0.05

# High-CNV malignant marker detection after low-CNV cluster removal
HIGH_CNV_CLUSTER_MARKER_LOGFC <- 1
HIGH_CNV_CLUSTER_MARKER_MINPCT <- 0.25
HIGH_CNV_CLUSTER_MARKER_PADJ <- 0.05

# Automatic low-CNV malignant cluster selection
LOW_CNV_AUTO_QUANTILE <- 0.25
LOW_CNV_MIN_CLUSTERS <- 1
LOW_CNV_MAX_FRACTION <- 0.35

# Reference selection mode for refined inferCNV
USE_MANUAL_LOW_CNV_CLUSTERS <- FALSE
MANUAL_LOW_CNV_CLUSTERS <- c(
  # Example: "Cluster 0", "Cluster 2", "Cluster 4"
)

# inferCNV
INITIAL_INFERCNV_CUTOFF <- 0.1
REFINED_INFERCNV_CUTOFF <- 0.1
MIN_LOW_CNV_REFERENCE_CELLS <- 10

# Figure settings
FIG_DPI <- 600
FONT_FAMILY <- "Arial"

# Step 5: helper functions

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

# Step 6: sample info

sample_info <- data.frame(
  gsm = c("GSM6251294", "GSM6251297", "GSM6251300"),
  sample = c("LSCC1", "LSCC2", "LSCC3"),
  stringsAsFactors = FALSE
)

# Step 7: load raw files

if (length(list.files(scrna_rawdir, recursive = TRUE, full.names = TRUE)) == 0) {
  tar_file <- file.choose()
  untar(tar_file, exdir = scrna_rawdir)
}

raw_files <- list.files(scrna_rawdir, recursive = TRUE, full.names = TRUE)

# Step 8: create seurat objects and apply basic qc filters

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

# Step 9: merge, normalize, hvg, pca, harmony, clustering, t-sne

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

# Step 10: findallmarkers for global clusters

# Seurat v5 can keep several RNA layers after merging samples.
# To avoid empty marker tables caused by unjoined layers, this step builds a
# dedicated joined assay for differential expression only.

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

# Keep the old filename so older downstream notes/paths do not break.
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

# Remove temporary object to reduce memory usage.
rm(seurat_de_obj)
gc()

# Step 11: simple cell type annotation using canonical markers

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

# Step 12: initial infercnv and cnv score for malignant cell detection

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

# Step 13: final cell type label with malignant cells

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

# Step 14: dotplot of canonical cell-type markers

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

# Step 15: re-clustering of malignant epithelial cells

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

# Step 16: initial cnv violin + weak low-cnv reference selection + refined infercnv

# Correct logic:
# 1) Use the INITIAL CNV score (CNV_score) after malignant re-clustering.
# 2) Draw the initial CNV violin plot to inspect weak/low-CNV clusters.
# 3) Select weak/low-CNV clusters from the INITIAL CNV summary.
# 4) Use these weak/low-CNV clusters as reference groups for refined inferCNV.
# 5) Only after this step, calculate Refined_CNV_score.

# Step 17: initial cnv score violin plot for reference selection

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

# Manuscript panel E filename for the combined A-G figure
save_plot(
  p_initial_cnv_violin,
  "Figure_SC_01E_Initial_CNV_score_violin_for_reference_selection.png",
  width = 8,
  height = 6
)

# Step 18: select weak / low-cnv clusters from initial cnv score

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

# Keep selected clusters in numeric cluster order.
low_cnv_clusters <- cluster_levels_num[cluster_levels_num %in% low_cnv_clusters]

low_cnv_reference_cells <- colnames(malignant_obj)[
  malignant_obj$Cluster_label %in% low_cnv_clusters
]

# Ensure enough reference cells by adding the next lowest-initial-CNV clusters if needed.
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

# Keep older filename for compatibility.
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

# Step 19: run refined infercnv using initial weak low-cnv clusters as reference

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

# Step 20: violin plot of refined cnv scores

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

# Keep the older filename for compatibility with previous notes.
save_plot(
  p_refined_cnv_violin,
  "Figure_SC_01E_Refined_CNV_score_violin_malignant_subclusters.png",
  width = 8,
  height = 6
)

# Manuscript panel G filename for the combined A-G figure
save_plot(
  p_refined_cnv_violin,
  "Figure_SC_01G_Refined_CNV_score_violin_malignant_subclusters.png",
  width = 8,
  height = 6
)

# Step 21: use original infercnv output image as panel f

# The custom pheatmap heatmap is intentionally removed.
# Instead, the original inferCNV output image generated by inferCNV is used.

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

# Step 22: remove low-cnv malignant clusters

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

# Step 23: high-cnv cluster-specific marker detection after low-cnv removal

# Article-like marker strategy:
# 1) Low-CNV malignant clusters have already been selected from the INITIAL CNV score.
# 2) These Low-CNV clusters have already been used as reference for refined inferCNV.
# 3) Low-CNV clusters are removed from the malignant object.
# 4) FindAllMarkers is run among the remaining High-CNV malignant subclusters only.
# This means marker genes are detected across the strong/high-CNV malignant
# clusters, not by comparing High-CNV malignant cells against Low-CNV cells.
# No maximum gene cap is applied; all genes passing logFC >= 1 and adjusted p-value filtering are saved.

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

# Step 24: prepare high-cnv malignant object for de analysis

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

# Step 25: run findallmarkers among remaining high-cnv malignant clusters

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

# Step 26: final filtering and top-gene selection

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

# If FDR filtering is too strict, use the top ranked markers from the full table.
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

# Keep all unique genes passing the logFC and adjusted p-value filters.
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

# Step 27: save final marker files

# Main full High-CNV cluster marker table
write.csv(
  high_cnv_cluster_markers_all,
  file.path(
    scrna_figdir,
    "High_CNV_Malignant_Cluster_FindAllMarkers_after_lowCNV_removed_logFC1_all.csv"
  ),
  row.names = FALSE
)

# Main final detailed CSV for downstream analysis
write.csv(
  final_high_cnv_markers,
  file.path(scrna_figdir, "Final_High_CNV_Malignant_Genes.csv"),
  row.names = FALSE
)

# Final gene names only
write.csv(
  final_high_cnv_gene_names,
  file.path(scrna_figdir, "Final_High_CNV_Malignant_Gene_Names.csv"),
  row.names = FALSE
)

writeLines(
  final_high_cnv_gene_names$Gene,
  file.path(scrna_figdir, "Final_High_CNV_Malignant_Gene_Names.txt")
)

# Compatibility files for downstream bulk/overlap code
write.csv(
  final_high_cnv_markers,
  file.path(scrna_figdir, "High_CNV_Malignant_FindAllMarkers_after_lowCNV_removed_logFC1.csv"),
  row.names = FALSE
)

writeLines(
  final_high_cnv_gene_names$Gene,
  file.path(scrna_figdir, "High_CNV_Malignant_marker_genes_after_lowCNV_removed_logFC1.txt")
)

# Older names also kept for safety
write.csv(
  final_high_cnv_markers,
  file.path(scrna_figdir, "High_CNV_Malignant_specific_markers_FINAL.csv"),
  row.names = FALSE
)

writeLines(
  final_high_cnv_gene_names$Gene,
  file.path(scrna_figdir, "High_CNV_Malignant_marker_genes_FINAL.txt")
)

# Step 28: save corrected objects

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

# Step 29: save objects

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

# Step 30: combined manuscript figure

# Figure order for the main manuscript text:
# A = global clusters
# B = final annotated cell types
# C = DotPlot of canonical cell-type markers
# D = malignant epithelial subclusters
# E = initial CNV score violin plot for weak/low-CNV reference selection
# F = original inferCNV heatmap image from refined inferCNV
# G = refined CNV score violin plot for final Low-CNV/High-CNV separation
# Improved layout:
# Row 1: A | B
# Row 2: C (full width)
# Row 3: D | E
# Row 4: F (full width, large)
# Row 5: G (full width)
# To make the combined panel cleaner, redundant legends are removed from
# panels A, B, D, E, and G. Panel C keeps its DotPlot legends, and panel F keeps
# the original inferCNV legend embedded inside the inferCNV image.

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
