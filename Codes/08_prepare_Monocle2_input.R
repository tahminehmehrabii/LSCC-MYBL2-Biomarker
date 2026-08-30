rm(list = ls(all.names = TRUE))
invisible(gc())

options(
  stringsAsFactors = FALSE,
  scipen = 100,
  timeout = 7200,
  width = 160
)

set.seed(123)

# STEP 01: Set paths and parameters.

if (as.character(getRversion()) != "4.4.3") {
  stop("Use R 4.4.3. Current version: ", as.character(getRversion()))
}

SC_PROJECT_DIR <- "D:/LSCC/ScRNAseq_Results/GSE206332"
SC_RESULTS_DIR <- file.path(SC_PROJECT_DIR, "Results")
SC_RDS_DIR <- file.path(SC_RESULTS_DIR, "rds")
SC_FIG_DIR <- file.path(SC_RESULTS_DIR, "figures")

OUTPUT_DIR <- "D:/LSCC/Results_LSCC/Pseudotime_Monocle"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

INPUT_RDS <- file.path(OUTPUT_DIR, "monocle_input.rds")
CELL_BALANCE_CSV <- file.path(OUTPUT_DIR, "cell_balance.csv")
ORDERING_GENES_CSV <- file.path(OUTPUT_DIR, "ordering_genes.csv")

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

# STEP 02: Load packages.

required_packages <- c("Seurat", "SeuratObject", "Matrix")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0L) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

# STEP 03: Define helper functions.

find_first_existing <- function(candidates, label) {
  existing <- candidates[file.exists(candidates)]
  
  if (length(existing) == 0L) {
    stop(label, " was not found. Checked:\n", paste(candidates, collapse = "\n"))
  }
  
  normalizePath(existing[1], winslash = "/", mustWork = TRUE)
}

cluster_number <- function(x) {
  x <- trimws(as.character(x))
  output <- suppressWarnings(as.integer(sub("^.*?([0-9]+)\\s*$", "\\1", x)))
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
    stop("The count matrix is not a sparse Matrix object.")
  }
  
  if (length(x@x) == 0L) {
    stop("The count matrix has no non-zero values.")
  }
  
  if (any(!is.finite(x@x)) || any(x@x < 0)) {
    stop("The count matrix contains invalid values.")
  }
  
  if (any(abs(x@x - round(x@x)) > 1e-8)) {
    stop("The count matrix is not raw integer counts.")
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
    single_layer <- NULL
    
    if (length(layer_names) == 1L) {
      single_layer <- tryCatch(
        SeuratObject::LayerData(
          object,
          assay = assay,
          layer = layer_names[1],
          fast = FALSE
        ),
        error = function(e) NULL
      )
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
      
      return(as(single_layer[, retained_cells, drop = FALSE], "dgCMatrix"))
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
    stop("Fewer than five G2/M root-score genes are present.")
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
    stop("The marker CSV lacks a required column.")
  }
  
  required
}

# STEP 04: Locate input files.

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

# STEP 05: Load and validate objects.

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
  stop("Fewer than 100 shared High-CNV cells were found.")
}

analysis_obj <- subset(highcnv_de_obj, cells = shared_cells)

source_meta <- as.data.frame(
  highcnv_obj@meta.data[shared_cells, , drop = FALSE],
  stringsAsFactors = FALSE
)

required_metadata <- c("sample", "Cluster_label", "CNV_class")
missing_metadata <- setdiff(required_metadata, colnames(source_meta))

if (length(missing_metadata) > 0L) {
  stop("Missing metadata: ", paste(missing_metadata, collapse = ", "))
}

analysis_obj$sample <- as.character(source_meta$sample)
analysis_obj$Cluster_label <- as.character(source_meta$Cluster_label)
analysis_obj$CNV_class <- as.character(source_meta$CNV_class)

if (any(is.na(analysis_obj$CNV_class) | analysis_obj$CNV_class != "High-CNV malignant")) {
  stop("Only High-CNV malignant cells are allowed.")
}

analysis_obj$Malignant_subcluster <- factor(
  normalise_cluster(analysis_obj$Cluster_label),
  levels = FIXED_CLUSTER_ORDER
)

if (any(is.na(analysis_obj$Malignant_subcluster))) {
  stop("Unrecognised Cluster_label values were found.")
}

observed_subclusters <- unique(as.character(analysis_obj$Malignant_subcluster))

if (!setequal(observed_subclusters, FIXED_CLUSTER_ORDER)) {
  stop(
    "The expected ten High-CNV subclusters are not all present. Observed: ",
    paste(sort(observed_subclusters), collapse = ", ")
  )
}

ASSAY_NAME <- if ("RNA_DE_HIGH_CNV" %in% names(analysis_obj@assays)) {
  "RNA_DE_HIGH_CNV"
} else if ("RNA" %in% names(analysis_obj@assays)) {
  "RNA"
} else {
  stop("Neither RNA_DE_HIGH_CNV nor RNA assay is available.")
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
  stop("Raw counts and metadata are not aligned.")
}

if (!(TARGET_GENE %in% rownames(count_matrix))) {
  stop(TARGET_GENE, " is absent from the raw count matrix.")
}

# STEP 06: Balance cells.

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
  stop("Counts and metadata are not aligned after balancing.")
}

if (!setequal(unique(as.character(metadata$Malignant_subcluster)), FIXED_CLUSTER_ORDER)) {
  stop("At least one expected High-CNV subcluster disappeared after balancing.")
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

# STEP 07: Select ordering genes.

markers <- read.csv(
  MARKER_CSV,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(markers) == 0L) {
  stop("The marker CSV has no rows.")
}

marker_columns <- get_marker_columns(markers)

markers_clean <- data.frame(
  Gene = as.character(markers[[marker_columns["Gene"]]]),
  cluster = as.character(markers[[marker_columns["cluster"]]]),
  p_val_adj = suppressWarnings(as.numeric(markers[[marker_columns["p_val_adj"]]])),
  avg_logFC = suppressWarnings(as.numeric(markers[[marker_columns["avg_logFC"]]])),
  pct.1 = suppressWarnings(as.numeric(markers[[marker_columns["pct.1"]]])),
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
  stop("No marker genes passed the required filters.")
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
    "Too few ordering genes remained: ",
    length(ordering_genes),
    ". Required minimum: ",
    MIN_ORDERING_GENES
  )
}

write.csv(
  data.frame(Gene = ordering_genes),
  ORDERING_GENES_CSV,
  row.names = FALSE
)

# STEP 08: Save Monocle input.

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
  stop("Final phenotype metadata and count matrix are not aligned.")
}

saveRDS(
  list(
    schema_version = "LSCC_Monocle2_input_v2",
    count_mat = as(count_matrix, "dgCMatrix"),
    pd_df = phenotype_data,
    ordering_genes = ordering_genes,
    target_gene = TARGET_GENE,
    root_column = "G2M_root_score",
    expected_subclusters = FIXED_CLUSTER_ORDER,
    max_cells_per_sample_subcluster = MAX_CELLS_PER_SAMPLE_SUBCLUSTER,
    root_rule = "Root = state with the lowest median G2M_root_score.",
    source_highcnv_object_rds = HIGHCNV_OBJECT_RDS,
    source_highcnv_de_object_rds = HIGHCNV_DE_OBJECT_RDS,
    source_marker_csv = MARKER_CSV,
    raw_count_assay = ASSAY_NAME
  ),
  INPUT_RDS
)

stopifnot(file.exists(INPUT_RDS))
