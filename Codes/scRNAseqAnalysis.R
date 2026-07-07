# Step 1: CLEAN ENVIRONMENT

rm(list = ls())
invisible(gc())

set.seed(123)
options(stringsAsFactors = FALSE)
options(scipen = 100)

# Step 2: LIBRARIES

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

# Step 3: PATHS

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

# Step 4: PARAMETERS

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

# Step 5: HELPER FUNCTIONS

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

# Step 6: SAMPLE INFO

sample_info <- data.frame(
  gsm = c("GSM6251294", "GSM6251297", "GSM6251300"),
  sample = c("LSCC1", "LSCC2", "LSCC3"),
  stringsAsFactors = FALSE
)

# Step 7: LOAD RAW FILES

if (length(list.files(scrna_rawdir, recursive = TRUE, full.names = TRUE)) == 0) {
  tar_file <- file.choose()
  untar(tar_file, exdir = scrna_rawdir)
}

raw_files <- list.files(scrna_rawdir, recursive = TRUE, full.names = TRUE)

# Step 8: CREATE SEURAT OBJECTS AND APPLY BASIC QC FILTERS

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
  file.path(scrna_figdir, "qc_report.csv"),
  row.names = FALSE
)

# Step 9: MERGE, NORMALIZE, HVG, PCA, HARMONY, CLUSTERING, t-SNE

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
  "fig01a_global.png",
  width = 8,
  height = 6
)

# Step 10: FINDALLMARKERS FOR GLOBAL CLUSTERS

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
  file.path(scrna_figdir, "global_markers_all.csv"),
  row.names = FALSE
)

write.csv(
  global_cluster_markers_all,
  file.path(scrna_figdir, "global_markers_compat.csv"),
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
  file.path(scrna_figdir, "global_markers_sig.csv"),
  row.names = FALSE
)

rm(seurat_de_obj)
invisible(gc())

# Step 11: SIMPLE CELL TYPE ANNOTATION USING CANONICAL MARKERS

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
  file.path(scrna_figdir, "celltypes.csv"),
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
  "fig01b_celltypes_initial.png",
  width = 8,
  height = 6
)

# Step 12: INITIAL inferCNV AND CNV SCORE FOR MALIGNANT CELL DETECTION

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
  file.path(scrna_figdir, "cnv_initial_report.csv"),
  row.names = FALSE
)

# Step 13: FINAL CELL TYPE LABEL WITH MALIGNANT CELLS

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
  "fig01b_celltypes.png",
  width = 8,
  height = 6
)

# Step 14: DOTPLOT OF CANONICAL CELL-TYPE MARKERS

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
    "fig01c_markers.png",
    width = 12,
    height = 6
  )
  
} else {
  
  p_dotplot <- ggplot() +
    annotate("text", x = 0, y = 0, label = "Not enough marker genes for DotPlot") +
    theme_void()
}

# Step 15: RE-CLUSTERING OF MALIGNANT EPITHELIAL CELLS

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
  "fig01d_malignant.png",
  width = 8,
  height = 6
)

cnv_initial_table <- malignant_obj@meta.data %>%
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
  cnv_initial_table,
  file.path(scrna_figdir, "cnv_initial_clusters.csv"),
  row.names = FALSE
)

# Step 16: INITIAL CNV VIOLIN + WEAK LOW-CNV REFERENCE SELECTION + REFINED inferCNV

# Step 17: INITIAL CNV SCORE VIOLIN PLOT FOR REFERENCE SELECTION

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
  "fig01e_cnv_initial_check.png",
  width = 8,
  height = 6
)

save_plot(
  p_initial_cnv_violin,
  "fig01e_cnv_initial.png",
  width = 8,
  height = 6
)

# Step 18: SELECT WEAK / LOW-CNV CLUSTERS FROM INITIAL CNV SCORE

cnv_initial_ranked <- cnv_initial_table %>%
  dplyr::arrange(Median_initial_CNV_score, Cluster_num)

n_malignant_clusters <- nrow(cnv_initial_ranked)

if (n_malignant_clusters < 2) {
  stop("At least two malignant subclusters are required for low-CNV reference selection.")
}

low_cnv_auto_cutoff <- as.numeric(
  quantile(
    cnv_initial_ranked$Median_initial_CNV_score,
    probs = LOW_CNV_AUTO_QUANTILE,
    na.rm = TRUE
  )
)

max_low_cnv_clusters <- floor(n_malignant_clusters * LOW_CNV_MAX_FRACTION)
max_low_cnv_clusters <- max(LOW_CNV_MIN_CLUSTERS, max_low_cnv_clusters)
max_low_cnv_clusters <- min(max_low_cnv_clusters, n_malignant_clusters - 1)

low_cnv_selection_table <- cnv_initial_ranked %>%
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
    as.character(cnv_initial_ranked$Cluster_label)
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
    low_cnv_clusters <- cnv_initial_ranked %>%
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
  
  ordered_clusters <- as.character(cnv_initial_ranked$Cluster_label)
  
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
  file.path(scrna_figdir, "low_cnv_selection.csv"),
  row.names = FALSE
)

write.csv(
  weak_low_cnv_reference_clusters_df,
  file.path(scrna_figdir, "low_cnv_clusters.csv"),
  row.names = FALSE
)

write.csv(
  low_cnv_selection_table,
  file.path(scrna_figdir, "low_cnv_auto.csv"),
  row.names = FALSE
)

if (length(low_cnv_reference_cells) < MIN_LOW_CNV_REFERENCE_CELLS) {
  stop(
    "Too few weak/low-CNV reference cells after selection. Found: ",
    length(low_cnv_reference_cells)
  )
}

# Step 19: RUN REFINED inferCNV USING INITIAL WEAK LOW-CNV CLUSTERS AS REFERENCE

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

refined_infercnv_dir <- file.path(infercnv_out, "refined_infercnv")
dir.create(refined_infercnv_dir, recursive = TRUE, showWarnings = FALSE)

anno_refined_file <- file.path(refined_infercnv_dir, "anno_refined.txt")

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

cnv_refined_table <- malignant_obj@meta.data %>%
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
  cnv_refined_table,
  file.path(scrna_figdir, "cnv_refined_clusters.csv"),
  row.names = FALSE
)

write.csv(
  malignant_obj@meta.data,
  file.path(scrna_rdsdir, "malignant_cnv_meta.csv"),
  row.names = TRUE
)

# Step 20: VIOLIN PLOT OF REFINED CNV SCORES

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
  "fig01g_cnv_refined_compat.png",
  width = 8,
  height = 6
)

save_plot(
  p_refined_cnv_violin,
  "fig01g_cnv_refined.png",
  width = 8,
  height = 6
)

# Step 21: USE ORIGINAL inferCNV OUTPUT IMAGE AS PANEL F

old_custom_heatmap_file <- file.path(
  scrna_figdir,
  "fig01f_infercnv.png"
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
  "fig01f_infercnv.png"
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

# Step 22: REMOVE LOW-CNV MALIGNANT CLUSTERS

high_cnv_malignant_obj <- subset(
  malignant_obj,
  subset = CNV_class == "High-CNV malignant"
)

write.csv(
  malignant_obj@meta.data,
  file.path(scrna_rdsdir, "malignant_all_meta.csv"),
  row.names = TRUE
)

write.csv(
  high_cnv_malignant_obj@meta.data,
  file.path(scrna_rdsdir, "highcnv_meta.csv"),
  row.names = TRUE
)

# Step 23: HIGH-CNV CLUSTER-SPECIFIC MARKER DETECTION AFTER LOW-CNV REMOVAL

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

# Step 24: Prepare High-CNV malignant object for DE analysis

high_cnv_malignant_obj$Cluster_label <- droplevels(
  high_cnv_malignant_obj$Cluster_label
)

high_cnv_de_obj <- make_DE_assay(
  high_cnv_malignant_obj,
  assay_name = "RNA_DE_HIGH_CNV"
)

Idents(high_cnv_de_obj) <- "Cluster_label"

if (length(levels(Idents(high_cnv_de_obj))) < 2) {
  stop(
    "FindAllMarkers needs at least two High-CNV malignant clusters after low-CNV removal. ",
    "Current number of High-CNV clusters: ",
    length(levels(Idents(high_cnv_de_obj)))
  )
}

# Step 25: Run FindAllMarkers among remaining High-CNV malignant clusters

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
    "highcnv_markers_all.csv"
  ),
  row.names = FALSE
)

# Step 26: Final filtering and top-gene selection

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

# Step 27: Save final marker files

write.csv(
  high_cnv_cluster_markers_all,
  file.path(
    scrna_figdir,
    "highcnv_markers_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  final_high_cnv_markers,
  file.path(scrna_figdir, "highcnv_markers.csv"),
  row.names = FALSE
)

write.csv(
  final_high_cnv_gene_names,
  file.path(scrna_figdir, "highcnv_genes.csv"),
  row.names = FALSE
)

writeLines(
  final_high_cnv_gene_names$Gene,
  file.path(scrna_figdir, "highcnv_genes.txt")
)

write.csv(
  final_high_cnv_markers,
  file.path(scrna_figdir, "highcnv_markers_compat.csv"),
  row.names = FALSE
)

writeLines(
  final_high_cnv_gene_names$Gene,
  file.path(scrna_figdir, "highcnv_genes_compat.txt")
)

write.csv(
  final_high_cnv_markers,
  file.path(scrna_figdir, "highcnv_markers_final.csv"),
  row.names = FALSE
)

writeLines(
  final_high_cnv_gene_names$Gene,
  file.path(scrna_figdir, "highcnv_genes_final.txt")
)

# Step 28: Save corrected objects

saveRDS(
  high_cnv_de_obj,
  file.path(scrna_rdsdir, "highcnv_de_obj.rds")
)

saveRDS(
  high_cnv_malignant_obj,
  file.path(scrna_rdsdir, "highcnv_obj.rds")
)

# Step 29: SAVE OBJECTS

saveRDS(
  seurat_obj,
  file.path(scrna_rdsdir, "seurat_final.rds")
)

saveRDS(
  malignant_obj,
  file.path(scrna_rdsdir, "malignant_cnv_obj.rds")
)

saveRDS(
  high_cnv_malignant_obj,
  file.path(scrna_rdsdir, "highcnv_obj.rds")
)

# Step 30: COMBINED MANUSCRIPT FIGURE

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
  filename = file.path(scrna_figdir, "fig01_combined.png"),
  plot = fig_sc_combined,
  width = 13.5,
  height = 22.5,
  dpi = FIG_DPI,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(scrna_figdir, "fig01_combined_compat.png"),
  plot = fig_sc_combined,
  width = 13.5,
  height = 22.5,
  dpi = FIG_DPI,
  bg = "white",
  limitsize = FALSE
)

# Step 32: Rebuild manuscript figures

rm(list = ls())
invisible(gc())

options(
  stringsAsFactors = FALSE,
  scipen = 100,
  width = 150
)

# Step 1: REQUIRED PACKAGES

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

# Step 2: PATHS

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
  "fig_rebuild"
)

dir.create(
  figure_out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Step 3: FIGURE SETTINGS

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

# Step 4: HELPER FUNCTIONS

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

# Step 5: LOAD EXISTING OBJECTS ONLY

seurat_rds_file <- first_existing(
  c(
    file.path(
      scrna_rdsdir,
      "seurat_final.rds"
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
      "malignant_cnv_obj.rds"
    ),
    find_first_file(
      directory = scrna_rdsdir,
      pattern = "MALIGNANT_all_object.*\\.rds$"
    )
  )
)

if (is.na(seurat_rds_file)) {
  stop(
    "Could not find seurat_final.rds in:\n",
    scrna_rdsdir
  )
}

if (is.na(malignant_rds_file)) {
  stop(
    "Could not find malignant_cnv_obj.rds in:\n",
    scrna_rdsdir
  )
}

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

# Step 6: PANEL A — GLOBAL t-SNE CLUSTERS

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
  filename_stem = "fig01a_global",
  dir_path = figure_out_dir,
  width = 8.2,
  height = 6.3
)

# Step 7: PANEL B — FINAL CELL-TYPE ANNOTATION

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
  filename_stem = "fig01b_celltypes",
  dir_path = figure_out_dir,
  width = 8.2,
  height = 6.3
)

# Step 8: PANEL C — MALIGNANT EPITHELIAL SUBCLUSTERS

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
  filename_stem = "fig01c_malignant",
  dir_path = figure_out_dir,
  width = 8.2,
  height = 6.3
)

# Step 9: PANEL D — INITIAL CNV SCORE VIOLIN PLOT

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
  filename_stem = "fig01d_cnv_initial",
  dir_path = figure_out_dir,
  width = 8.2,
  height = 6.3
)

# Step 10: PANEL E — LARGE REFINED inferCNV IMAGE

infercnv_png_source <- first_existing(
  c(
    file.path(
      scrna_figdir,
      "fig01f_infercnv.png"
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
  filename_stem = "fig01e_infercnv",
  dir_path = figure_out_dir,
  width = 14.5,
  height = 11.5
)

# Step 11: COMBINED MANUSCRIPT FIGURE

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
  filename_stem = "fig01_rebuild",
  dir_path = figure_out_dir,
  width = 16.5,
  height = 20.8
)
