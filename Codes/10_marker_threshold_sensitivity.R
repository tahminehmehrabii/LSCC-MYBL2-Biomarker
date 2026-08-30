# LSCC single-cell marker-threshold sensitivity analysis
# Compares permissive, primary, and stringent marker definitions.

rm(list = ls())
gc()

set.seed(123)

options(
  stringsAsFactors = FALSE,
  scipen = 100,
  width = 200
)

# Step 1: packages

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "data.table",
  "dplyr",
  "ggplot2",
  "patchwork",
  "scales",
  "tidyr"
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
    "Install these packages first:\n",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({

  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)

})

# Step 2: paths

BASE_RESULTS_DIR <- "E:/LSCC/Results_LSCC/ML"

SCRNA_RESULTS_DIR <-
  "E:/LSCC/ScRNAseq_Results/GSE206332/Results"

SCRNA_RDS_FILE <- file.path(
  SCRNA_RESULTS_DIR,
  "rds",
  "HIGH_CNV_MALIGNANT_object_DE_fixed.rds"
)

# Original 707-gene list.
# This is used to verify that Primary exactly reproduces the original analysis.

ORIGINAL_707_FILE <- file.path(
  SCRNA_RESULTS_DIR,
  "figures",
  "Final_High_CNV_Malignant_Gene_Names.csv"
)

BULK_UP_FILE <- file.path(
  BASE_RESULTS_DIR,
  "DEG_up_logFC2_adjP0.01.txt"
)

WGCNA_FILE <- file.path(
  BASE_RESULTS_DIR,
  "WGCNA_feature_genes.csv"
)

OUTPUT_DIR <- file.path(
  BASE_RESULTS_DIR,
  "Marker_Threshold_Sensitivity"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Step 3: check input files

required_files <- c(
  SCRNA_RDS_FILE,
  ORIGINAL_707_FILE,
  BULK_UP_FILE,
  WGCNA_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  stop(
    "Required file(s) missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

# Step 4: sensitivity definitions

# IMPORTANT:
# Primary is the exact definition that reproduced the saved 707 genes.
# Only log2FC is changed in the sensitivity analysis.
# min.pct and adjusted P remain fixed.
# Marker IDs are extracted from FindAllMarkers row names to match the original pipeline.

marker_definitions <- data.frame(

  Definition = factor(
    c(
      "Permissive",
      "Primary",
      "Stringent"
    ),
    levels = c(
      "Permissive",
      "Primary",
      "Stringent"
    )
  ),

  logfc_threshold = c(
    1.50,
    2.00,
    2.50
  ),

  min_pct = c(
    0.50,
    0.50,
    0.50
  ),

  adjusted_p_threshold = c(
    0.001,
    0.001,
    0.001
  ),

  stringsAsFactors = FALSE
)

TARGET_GENE <- "MYBL2"

FIG_DPI <- 600

FONT_FAMILY <- "Arial"

if (.Platform$OS.type == "windows") {

  try(
    grDevices::windowsFonts(
      Arial = grDevices::windowsFont("Arial")
    ),
    silent = TRUE
  )
}

# Step 5: helper functions

clean_gene_vector <- function(x) {

  x <- toupper(
    trimws(
      as.character(x)
    )
  )

  sort(
    unique(
      x[
        !is.na(x) &
          nzchar(x)
      ]
    )
  )
}

read_gene_list <- function(path) {

  if (grepl("\\.csv$", path, ignore.case = TRUE)) {

    dat <- as.data.frame(
      data.table::fread(path),
      check.names = FALSE
    )

    gene_col <- intersect(
      c(
        "Gene",
        "Genes",
        "gene",
        "Gene.symbol",
        "Symbol",
        "symbol"
      ),
      colnames(dat)
    )

    if (length(gene_col) == 0) {

      stop(
        "No recognized gene column in: ",
        path
      )
    }

    return(
      clean_gene_vector(
        dat[[gene_col[[1]]]]
      )
    )
  }

  clean_gene_vector(
    readLines(
      path,
      warn = FALSE
    )
  )
}

extract_marker_gene <- function(marker_table) {

  # IMPORTANT:
  # The ORIGINAL single-cell pipeline created the final "Gene" field from
  # Marker identifiers are retained from FindAllMarkers row names.
  #
  # To reproduce the original analysis exactly, use the row names here.

  ids <- rownames(marker_table)

  if (is.null(ids) || length(ids) != nrow(marker_table)) {
    stop("FindAllMarkers row names are unavailable or inconsistent.")
  }

  as.character(ids)
}

get_fc_column <- function(marker_table) {

  if ("avg_log2FC" %in% colnames(marker_table)) {

    return("avg_log2FC")
  }

  if ("avg_logFC" %in% colnames(marker_table)) {

    return("avg_logFC")
  }

  stop(
    "Neither avg_log2FC nor avg_logFC exists in FindAllMarkers output."
  )
}

safe_name <- function(x) {

  gsub(
    "[^A-Za-z0-9]+",
    "_",
    as.character(x)
  )
}

save_figure <- function(
    plot_object,
    stem,
    width = 11.5,
    height = 6) {

  ggplot2::ggsave(
    filename = file.path(
      OUTPUT_DIR,
      paste0(stem, ".png")
    ),
    plot = plot_object,
    width = width,
    height = height,
    dpi = FIG_DPI,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = file.path(
      OUTPUT_DIR,
      paste0(stem, ".tiff")
    ),
    plot = plot_object,
    width = width,
    height = height,
    dpi = FIG_DPI,
    compression = "lzw",
    bg = "white"
  )

  tryCatch(
    {

      ggplot2::ggsave(
        filename = file.path(
          OUTPUT_DIR,
          paste0(stem, ".pdf")
        ),
        plot = plot_object,
        width = width,
        height = height,
        device = grDevices::cairo_pdf,
        bg = "white"
      )

    },
    error = function(e) {

      warning(
        "PDF export failed: ",
        conditionMessage(e)
      )
    }
  )
}

theme_manuscript <- function() {

  theme_bw(
    base_size = 11,
    base_family = FONT_FAMILY
  ) +

    theme(

      text = element_text(
        color = "black",
        family = FONT_FAMILY
      ),

      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0.5
      ),

      plot.tag = element_text(
        face = "bold",
        size = 14
      ),

      axis.title = element_text(
        face = "bold"
      ),

      axis.text = element_text(
        color = "black"
      ),

      panel.grid.minor =
        element_blank(),

      panel.grid.major =
        element_line(
          color = "grey92",
          linewidth = 0.3
        ),

      panel.border =
        element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.5
        ),

      legend.title =
        element_text(
          face = "bold"
        ),

      plot.background =
        element_rect(
          fill = "white",
          color = NA
        )
    )
}

# Step 6: load original 707 genes

original_707_genes <- read_gene_list(
  ORIGINAL_707_FILE
)

cat("\n============================================================\n")
cat("ORIGINAL HIGH-CNV MARKER LIST\n")
cat("============================================================\n")

cat(
  "Original unique genes:",
  length(original_707_genes),
  "\n"
)

cat(
  "MYBL2 present:",
  TARGET_GENE %in% original_707_genes,
  "\n"
)

if (length(original_707_genes) != 707) {

  stop(
    "Original marker file does not contain exactly 707 unique genes."
  )
}

# Step 7: load bulk deg and wgcna gene sets

Bulk_UP <- read_gene_list(
  BULK_UP_FILE
)

WGCNA_Genes <- read_gene_list(
  WGCNA_FILE
)

cat("\nBulk upregulated genes:", length(Bulk_UP), "\n")
cat("WGCNA genes:", length(WGCNA_Genes), "\n")

# Step 8: load high-cnv seurat object

high_cnv_de_obj <- readRDS(
  SCRNA_RDS_FILE
)

if (!inherits(high_cnv_de_obj, "Seurat")) {

  stop(
    "The High-CNV RDS file is not a Seurat object."
  )
}

if (!"Cluster_label" %in%
    colnames(high_cnv_de_obj@meta.data)) {

  stop(
    "Cluster_label is missing from metadata."
  )
}

if ("RNA_DE_HIGH_CNV" %in%
    names(high_cnv_de_obj@assays)) {

  DefaultAssay(
    high_cnv_de_obj
  ) <- "RNA_DE_HIGH_CNV"

} else {

  stop(
    "RNA_DE_HIGH_CNV assay is missing."
  )
}

high_cnv_de_obj <- tryCatch(

  JoinLayers(
    high_cnv_de_obj,
    assay = "RNA_DE_HIGH_CNV"
  ),

  error = function(e) {

    high_cnv_de_obj
  }
)

high_cnv_de_obj$Cluster_label <-
  droplevels(
    factor(
      high_cnv_de_obj$Cluster_label
    )
  )

Idents(
  high_cnv_de_obj
) <- "Cluster_label"

cat("\n============================================================\n")
cat("HIGH-CNV OBJECT\n")
cat("============================================================\n")

cat(
  "Cells:",
  ncol(high_cnv_de_obj),
  "\n"
)

cat(
  "Genes:",
  nrow(high_cnv_de_obj),
  "\n"
)

cat(
  "Clusters:",
  nlevels(Idents(high_cnv_de_obj)),
  "\n\n"
)

print(
  table(
    Idents(high_cnv_de_obj)
  )
)

# Step 9: run sensitivity analysis

summary_rows <- list()

membership_rows <- list()

mybl2_rows <- list()

all_results <- list()

for (
  i in seq_len(
    nrow(marker_definitions)
  )
) {

  definition <-
    as.character(
      marker_definitions$Definition[i]
    )

  logfc_i <-
    marker_definitions$logfc_threshold[i]

  minpct_i <-
    marker_definitions$min_pct[i]

  padj_i <-
    marker_definitions$adjusted_p_threshold[i]

  cat("\n\n############################################################\n")
  cat("RUNNING:", definition, "\n")
  cat("log2FC >=", logfc_i, "\n")
  cat("min.pct =", minpct_i, "\n")
  cat("adjusted P <", padj_i, "\n")
  cat("############################################################\n")

  markers_all <- FindAllMarkers(

    object = high_cnv_de_obj,

    assay = "RNA_DE_HIGH_CNV",

    only.pos = TRUE,

    test.use = "wilcox",

    logfc.threshold = logfc_i,

    min.pct = minpct_i,

    return.thresh = 1,

    verbose = FALSE
  )

  if (nrow(markers_all) == 0) {

    markers_all$Gene <- character(0)
    markers_all$Seurat_gene_symbol <- character(0)

    markers_sig <-
      markers_all

  } else {

    # Preserve Seurat's biological gene-symbol column for auditing, if present.
    if ("gene" %in% colnames(markers_all)) {
      markers_all$Seurat_gene_symbol <-
        toupper(trimws(as.character(markers_all$gene)))
    } else if ("Gene" %in% colnames(markers_all)) {
      markers_all$Seurat_gene_symbol <-
        toupper(trimws(as.character(markers_all$Gene)))
    } else {
      markers_all$Seurat_gene_symbol <- NA_character_
    }

    # CRITICAL FIX:
    # Use rownames exactly as the original High-CNV marker pipeline did.
    markers_all$Gene <-
      toupper(
        trimws(
          extract_marker_gene(
            markers_all
          )
        )
      )

    fc_col <-
      get_fc_column(
        markers_all
      )

    # IMPORTANT:
    # Explicitly reproduce the final filtering rule.

    markers_sig <-
      markers_all %>%

      dplyr::filter(

        !is.na(.data[[fc_col]]),

        !is.na(p_val_adj),

        .data[[fc_col]] >= logfc_i,

        p_val_adj < padj_i
      )
  }

  sc_genes <-
    clean_gene_vector(
      markers_sig$Gene
    )

  sc_bulk_genes <-
    sort(
      intersect(
        sc_genes,
        Bulk_UP
      )
    )

  final_genes <-
    sort(
      intersect(
        sc_bulk_genes,
        WGCNA_Genes
      )
    )

  # Direct intersection check

  direct_final <- sort(
    Reduce(
      intersect,
      list(
        sc_genes,
        Bulk_UP,
        WGCNA_Genes
      )
    )
  )

  if (!identical(
    final_genes,
    direct_final
  )) {

    stop(
      "Sequential and direct intersections disagree for ",
      definition
    )
  }

  # ---------------------------------------------------------------------------
  # PRIMARY MUST EXACTLY REPRODUCE ORIGINAL 707 GENES
  # ---------------------------------------------------------------------------

  if (definition == "Primary") {

    primary_exact_match <-
      setequal(
        sc_genes,
        original_707_genes
      )

    cat("\nPRIMARY VALIDATION\n")

    cat(
      "Primary genes:",
      length(sc_genes),
      "\n"
    )

    cat(
      "Original genes:",
      length(original_707_genes),
      "\n"
    )

    cat(
      "Exact same gene set:",
      primary_exact_match,
      "\n"
    )

    if (!primary_exact_match) {

      stop(
        paste0(
          "\nPRIMARY DID NOT reproduce the original 707-gene list.\n",
          "Primary genes = ",
          length(sc_genes),
          "\nOriginal genes = ",
          length(original_707_genes),
          "\nThe script is already using the original rowname-based marker IDs. ",
          "If this still fails, verify that the RDS and saved 707-marker CSV ",
          "come from the same analysis run."
        )
      )
    }
  }

  # ---------------------------------------------------------------------------
  # SAVE DEFINITION-SPECIFIC FILES
  # ---------------------------------------------------------------------------

  definition_file <-
    safe_name(
      definition
    )

  write.csv(
    markers_all,
    file.path(
      OUTPUT_DIR,
      paste0(
        definition_file,
        "_FindAllMarkers_all.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    markers_sig,
    file.path(
      OUTPUT_DIR,
      paste0(
        definition_file,
        "_significant_scRNA_markers.csv"
      )
    ),
    row.names = FALSE
  )

  writeLines(
    sc_genes,
    file.path(
      OUTPUT_DIR,
      paste0(
        definition_file,
        "_scRNA_marker_genes.txt"
      )
    )
  )

  writeLines(
    sc_bulk_genes,
    file.path(
      OUTPUT_DIR,
      paste0(
        definition_file,
        "_after_BulkUP.txt"
      )
    )
  )

  writeLines(
    final_genes,
    file.path(
      OUTPUT_DIR,
      paste0(
        definition_file,
        "_final_after_WGCNA.txt"
      )
    )
  )

  # ---------------------------------------------------------------------------
  # SUMMARY
  # ---------------------------------------------------------------------------

  summary_rows[[definition]] <-
    data.frame(

      Definition = definition,

      logfc_threshold = logfc_i,

      min_pct = minpct_i,

      adjusted_p_threshold = padj_i,

      N_scRNA_markers =
        length(sc_genes),

      N_after_scRNA_BulkUP =
        length(sc_bulk_genes),

      N_after_scRNA_BulkUP_WGCNA =
        length(final_genes),

      MYBL2_in_scRNA =
        TARGET_GENE %in% sc_genes,

      MYBL2_after_BulkUP =
        TARGET_GENE %in% sc_bulk_genes,

      MYBL2_final =
        TARGET_GENE %in% final_genes,

      stringsAsFactors = FALSE
    )

  # ---------------------------------------------------------------------------
  # MEMBERSHIP TABLE
  # ---------------------------------------------------------------------------

  all_union_genes <-
    sort(
      unique(
        c(
          sc_genes,
          sc_bulk_genes,
          final_genes
        )
      )
    )

  membership_rows[[definition]] <-
    data.frame(

      Definition = definition,

      Gene = all_union_genes,

      In_scRNA_markers =
        all_union_genes %in%
        sc_genes,

      After_BulkUP =
        all_union_genes %in%
        sc_bulk_genes,

      Final_after_WGCNA =
        all_union_genes %in%
        final_genes,

      stringsAsFactors = FALSE
    )

  # ---------------------------------------------------------------------------
  # MYBL2 DETAILS
  # ---------------------------------------------------------------------------

  mybl2_detail <-
    markers_sig %>%
    dplyr::filter(
      Gene == TARGET_GENE
    )

  if (nrow(mybl2_detail) > 0) {

    mybl2_detail$Definition <-
      definition

    mybl2_detail$MYBL2_final <-
      TARGET_GENE %in%
      final_genes

    mybl2_rows[[definition]] <-
      mybl2_detail

  } else {

    mybl2_rows[[definition]] <-
      data.frame(

        Definition =
          definition,

        Gene =
          TARGET_GENE,

        cluster =
          NA_character_,

        p_val =
          NA_real_,

        avg_log2FC =
          NA_real_,

        pct.1 =
          NA_real_,

        pct.2 =
          NA_real_,

        p_val_adj =
          NA_real_,

        MYBL2_final =
          FALSE,

        stringsAsFactors =
          FALSE
      )
  }

  all_results[[definition]] <-
    list(

      markers_all =
        markers_all,

      markers_significant =
        markers_sig,

      sc_genes =
        sc_genes,

      sc_bulk_genes =
        sc_bulk_genes,

      final_genes =
        final_genes
    )

  cat(
    "\nSignificant unique scRNA markers (original rowname convention):",
    length(sc_genes),
    "\n"
  )

  if ("Seurat_gene_symbol" %in% colnames(markers_sig)) {
    biological_symbols <- clean_gene_vector(markers_sig$Seurat_gene_symbol)
    cat(
      "Unique biological gene symbols in the same marker rows:",
      length(biological_symbols),
      "\n"
    )
  }

  cat(
    "After bulk DEG intersection:",
    length(sc_bulk_genes),
    "\n"
  )

  cat(
    "After WGCNA intersection:",
    length(final_genes),
    "\n"
  )

  cat(
    "MYBL2 final:",
    TARGET_GENE %in% final_genes,
    "\n"
  )
}

# Step 10: combine results

summary_table <-
  bind_rows(
    summary_rows
  ) %>%

  mutate(
    Definition = factor(
      Definition,
      levels = levels(
        marker_definitions$Definition
      )
    )
  ) %>%

  arrange(
    Definition
  )

membership_table <-
  bind_rows(
    membership_rows
  ) %>%

  mutate(
    Definition = factor(
      Definition,
      levels = levels(
        marker_definitions$Definition
      )
    )
  ) %>%

  arrange(
    Definition,
    Gene
  )

mybl2_table <-
  bind_rows(
    mybl2_rows
  ) %>%

  mutate(
    Definition = factor(
      Definition,
      levels = levels(
        marker_definitions$Definition
      )
    )
  ) %>%

  arrange(
    Definition
  )

# Step 11: final candidate membership across thresholds

final_union <-
  sort(
    unique(
      unlist(
        lapply(
          all_results,
          `[[`,
          "final_genes"
        )
      )
    )
  )

if (length(final_union) > 0) {

  final_membership <-
    expand.grid(

      Gene =
        final_union,

      Definition =
        levels(
          marker_definitions$Definition
        ),

      stringsAsFactors =
        FALSE
    )

  # Robust list indexing without rowwise factor/list ambiguity.
  final_membership$Retained <- mapply(
    FUN = function(g, d) {
      g %in% all_results[[as.character(d)]]$final_genes
    },
    g = final_membership$Gene,
    d = final_membership$Definition,
    USE.NAMES = FALSE
  )

} else {

  final_membership <-
    data.frame(

      Gene =
        character(0),

      Definition =
        character(0),

      Retained =
        logical(0)
    )
}

# Step 12: save tables

write.csv(

  summary_table,

  file.path(
    OUTPUT_DIR,
    "Table_S1_marker_threshold_sensitivity_summary.csv"
  ),

  row.names = FALSE
)

write.csv(

  membership_table,

  file.path(
    OUTPUT_DIR,
    "All_gene_stage_membership.csv"
  ),

  row.names = FALSE
)

write.csv(

  mybl2_table,

  file.path(
    OUTPUT_DIR,
    "Table_S2_MYBL2_marker_statistics.csv"
  ),

  row.names = FALSE
)

write.csv(

  final_membership,

  file.path(
    OUTPUT_DIR,
    "Final_candidate_membership_across_thresholds.csv"
  ),

  row.names = FALSE
)

# Step 13: figure — sequential filtering

count_long <-
  summary_table %>%

  select(

    Definition,

    `scRNA markers` =
      N_scRNA_markers,

    `After bulk DEG intersection` =
      N_after_scRNA_BulkUP,

    `After WGCNA intersection` =
      N_after_scRNA_BulkUP_WGCNA
  ) %>%

  pivot_longer(

    cols =
      -Definition,

    names_to =
      "Filtering_step",

    values_to =
      "Gene_count"
  )

count_long$Filtering_step <-
  factor(

    count_long$Filtering_step,

    levels = c(

      "scRNA markers",

      "After bulk DEG intersection",

      "After WGCNA intersection"
    )
  )

p_count <-

  ggplot(

    count_long,

    aes(
      x = Definition,
      y = Gene_count,
      fill = Filtering_step
    )
  ) +

  geom_col(

    position =
      position_dodge(
        width = 0.78
      ),

    width = 0.70,

    color = "black",

    linewidth = 0.25
  ) +

  geom_text(

    aes(
      label =
        scales::comma(
          Gene_count
        )
    ),

    position =
      position_dodge(
        width = 0.78
      ),

    vjust = -0.35,

    size = 3.2
  ) +

  scale_y_continuous(

    expand =
      expansion(
        mult = c(
          0,
          0.15
        )
      )
  ) +

  labs(

    x =
      "Single-cell marker definition",

    y =
      "Number of genes",

    fill =
      "Filtering step",

    title =
      "Sensitivity of sequential candidate selection to the single-cell marker threshold"
  ) +

  theme_manuscript() +

  theme(
    legend.position =
      "bottom"
  )

save_figure(

  plot_object =
    p_count,

  stem =
    "Figure_S1_MYBL2_marker_threshold_sensitivity",

  width =
    9,

  height =
    6
)

# Step 14: robustness result

all_robust <-
  all(
    summary_table$MYBL2_final
  )

primary_n <-
  summary_table$N_scRNA_markers[
    summary_table$Definition ==
      "Primary"
  ]

primary_reproduces_707 <-
  length(primary_n) == 1 &&
  primary_n == 707 &&
  setequal(
    all_results$Primary$sc_genes,
    original_707_genes
  )

# Step 15: analysis summary

threshold_text <-
  paste0(

    "log2FC thresholds of 1.5, 2.0, and 2.5, ",

    "with min.pct fixed at 0.50 and adjusted P < 0.001"
  )

if (all_robust) {

  results_sentence <-
    paste0(

      "Sensitivity analysis using ",
      threshold_text,
      " yielded ",

      paste(
        summary_table$N_scRNA_markers,
        collapse = ", "
      ),

      " high-CNV malignant epithelial marker genes under the permissive, ",
      "primary, and stringent definitions, respectively. ",

      "The primary definition reproduced the original 707-gene marker set exactly. ",

      "After sequential intersection with bulk upregulated DEGs and ",
      "WGCNA-derived tumor-associated genes, MYBL2 was retained under all ",
      "three definitions, supporting the robustness of its selection to ",
      "variation in the single-cell marker effect-size threshold."
    )

} else {

  failed <-
    paste(

      as.character(
        summary_table$Definition[
          !summary_table$MYBL2_final
        ]
      ),

      collapse = ", "
    )

  results_sentence <-
    paste0(

      "MYBL2 was not retained after the complete sequential intersection ",
      "under the following marker definition(s): ",
      failed,
      ". The threshold dependence should therefore be reported transparently."
    )
}

writeLines(
  c(
    "MARKER-THRESHOLD SENSITIVITY SUMMARY",
    results_sentence
  ),
  file.path(
    OUTPUT_DIR,
    "Marker_threshold_sensitivity_summary.txt"
  )
)

# Step 16: save complete results

saveRDS(

  list(

    marker_definitions =
      marker_definitions,

    original_707_genes =
      original_707_genes,

    primary_reproduces_707 =
      primary_reproduces_707,

    summary_table =
      summary_table,

    mybl2_table =
      mybl2_table,

    final_membership =
      final_membership,

    results =
      all_results,

    input_files =
      list(

        scRNA_RDS =
          SCRNA_RDS_FILE,

        original_707 =
          ORIGINAL_707_FILE,

        Bulk_UP =
          BULK_UP_FILE,

        WGCNA =
          WGCNA_FILE
      )
  ),

  file.path(
    OUTPUT_DIR,
    "Marker_threshold_sensitivity_results.rds"
  )
)

writeLines(

  capture.output(
    sessionInfo()
  ),

  file.path(
    OUTPUT_DIR,
    "sessionInfo.txt"
  )
)

# Step 17: final console output

cat("\n\n")
cat("======================================================================\n")
cat("MYBL2 MARKER-THRESHOLD SENSITIVITY ANALYSIS COMPLETED\n")
cat("======================================================================\n\n")

print(
  summary_table,
  row.names = FALSE
)

cat("\n")
cat(
  "Primary exactly reproduces original 707 genes: ",
  primary_reproduces_707,
  "\n",
  sep = ""
)

cat(
  "MYBL2 retained after final intersection under ALL definitions: ",
  all_robust,
  "\n",
  sep = ""
)

cat("\nOutput folder:\n")
cat(OUTPUT_DIR, "\n")

cat("\n======================================================================\n")
