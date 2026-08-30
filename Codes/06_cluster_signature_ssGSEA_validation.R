rm(list = ls(all.names = TRUE))
gc()
set.seed(123)

# Step 1: Setup
options(stringsAsFactors = FALSE, scipen = 100, timeout = 7200, width = 160)

if (as.character(getRversion()) != "4.4.3") {
  stop("Run this script only in R 4.4.3. Current R version: ", as.character(getRversion()))
}

BASE_RESULTS_DIR <- "D:/LSCC/Results_LSCC"
BULK_ML_DIR <- file.path(BASE_RESULTS_DIR, "ML")
SC_PROJECT_DIR <- "D:/LSCC/ScRNAseq_Results/GSE206332"
SC_RESULTS_DIR <- file.path(SC_PROJECT_DIR, "Results")
SC_FIG_DIR <- file.path(SC_RESULTS_DIR, "figures")
SC_RDS_DIR <- file.path(SC_RESULTS_DIR, "rds")
OUTPUT_DIR <- file.path(BASE_RESULTS_DIR, "LSCC_ssGSEA_C0_C8")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

TARGET_GENE <- "MYBL2"
SELECTED_SUBCLUSTERS <- c("Cluster 0", "Cluster 8")
MARKER_PADJ_CUTOFF <- 0.05
MARKER_LOGFC_CUTOFF <- 1.00
MARKER_PCT_CUTOFF <- 0.25
TOP_SIGNATURE_GENES <- 50L
MIN_SIGNATURE_GENES <- 10L
REMOVE_MITOCHONDRIAL_GENES <- TRUE
REMOVE_RIBOSOMAL_GENES <- TRUE

DEVELOPMENT_CPM_CANDIDATES <- c(
  file.path(BULK_ML_DIR, "train_discovery_merged_CPM_nonnegative.csv"),
  file.path(BULK_ML_DIR, "train_discovery_merged_CPM_nonnegative_noZero.csv")
)
EXTERNAL_CPM_CANDIDATES <- file.path(BULK_ML_DIR, "external_GSE130605_CPM.csv")
MARKER_CSV_CANDIDATES <- c(
  file.path(SC_FIG_DIR, "High_CNV_Malignant_Cluster_FindAllMarkers_after_lowCNV_removed_logFC1_all.csv"),
  file.path(SC_FIG_DIR, "High_CNV_Malignant_FindAllMarkers_after_lowCNV_removed_logFC1.csv"),
  file.path(SC_FIG_DIR, "Final_High_CNV_Malignant_Genes.csv"),
  file.path(SC_RDS_DIR, "High_CNV_Malignant_Cluster_FindAllMarkers_after_lowCNV_removed_logFC1_all.csv")
)

SIG_GENES_CSV <- file.path(OUTPUT_DIR, "sig_genes_c0_c8.csv")
SIG_GMT <- file.path(OUTPUT_DIR, "sig_c0_c8.gmt")
DISC_SCORES_CSV <- file.path(OUTPUT_DIR, "disc_scores.csv")
EXT_SCORES_CSV <- file.path(OUTPUT_DIR, "ext_scores.csv")
DISC_TESTS_CSV <- file.path(OUTPUT_DIR, "disc_tests.csv")
EXT_TESTS_CSV <- file.path(OUTPUT_DIR, "ext_tests.csv")
DISC_CORR_CSV <- file.path(OUTPUT_DIR, "disc_mybl2_corr.csv")
FINAL_FIGURE_STEM <- "fig08_ssgsea_c0_c8"

FONT_FAMILY <- "Arial"
FIG_DPI <- 600L
FIG_BACKGROUND <- "white"
FIGURE_W <- 12.8
FIGURE_H <- 11.2
BASE_TEXT_PT <- 9.5
AXIS_TEXT_PT <- 8.5
AXIS_TITLE_PT <- 10
LEGEND_TEXT_PT <- 8.5
LEGEND_TITLE_PT <- 9
PLOT_TITLE_PT <- 10
PANEL_TAG_PT <- 14
PANEL_BORDER_LWD <- 0.45
AXIS_LWD <- 0.40
GRID_LWD <- 0.30

COL_NORMAL <- "#2A9D8F"
COL_TUMOR <- "#E76F51"
COL_CLUSTER0 <- "#287A78"
COL_CLUSTER8 <- "#825B9B"
COL_FIT <- "#303030"
COL_BAND <- "#D9D9D9"
COL_EXTERNAL_PANEL_BG <- "#F3F7FB"
COL_EXTERNAL_PANEL_BORDER <- "#4C78A8"
DISCOVERY_COHORT_HEADER <- "Discovery cohort (GSE127165 + GSE142083)"
EXTERNAL_COHORT_HEADER <- "Independent validation cohort (GSE130605)"
MANUSCRIPT_MARGIN <- ggplot2::margin(t = 5.5, r = 6, b = 5.5, l = 6, unit = "pt")

# Step 2: Load packages
required_packages <- c("data.table", "dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "GSVA")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages) > 0L) {
  stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

if (.Platform$OS.type == "windows") {
  try(grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial")), silent = TRUE)
}

# Step 3: Define helpers
theme_manuscript <- function(show_grid = FALSE, legend_position = "right", panel_fill = FIG_BACKGROUND, panel_border_colour = "black") {
  ggplot2::theme_bw(base_size = BASE_TEXT_PT, base_family = FONT_FAMILY) +
    ggplot2::theme(
      text = ggplot2::element_text(family = FONT_FAMILY, color = "black"),
      plot.title = ggplot2::element_text(family = FONT_FAMILY, face = "bold", size = PLOT_TITLE_PT, hjust = 0.5, color = "black"),
      plot.tag = ggplot2::element_text(family = FONT_FAMILY, face = "bold", size = PANEL_TAG_PT, color = "black"),
      plot.tag.position = c(0.01, 0.99),
      axis.title = ggplot2::element_text(family = FONT_FAMILY, face = "bold", size = AXIS_TITLE_PT, color = "black"),
      axis.text = ggplot2::element_text(family = FONT_FAMILY, size = AXIS_TEXT_PT, color = "black"),
      axis.line = ggplot2::element_line(color = "black", linewidth = AXIS_LWD),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = AXIS_LWD),
      axis.ticks.length = grid::unit(2, "pt"),
      legend.position = legend_position,
      legend.title = ggplot2::element_text(family = FONT_FAMILY, face = "bold", size = LEGEND_TITLE_PT),
      legend.text = ggplot2::element_text(family = FONT_FAMILY, size = LEGEND_TEXT_PT),
      legend.background = ggplot2::element_rect(fill = FIG_BACKGROUND, color = NA),
      legend.key = ggplot2::element_rect(fill = FIG_BACKGROUND, color = NA),
      panel.background = ggplot2::element_rect(fill = panel_fill, color = NA),
      plot.background = ggplot2::element_rect(fill = FIG_BACKGROUND, color = NA),
      panel.border = ggplot2::element_rect(color = panel_border_colour, fill = NA, linewidth = PANEL_BORDER_LWD),
      panel.grid.major = if (show_grid) ggplot2::element_line(color = "grey92", linewidth = GRID_LWD) else ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = MANUSCRIPT_MARGIN
    )
}

save_plot_all_formats <- function(plot_obj, filename_stem, width = FIGURE_W, height = FIGURE_H) {
  filename_stem <- sub("\\.(png|tiff|tif|pdf)$", "", filename_stem, ignore.case = TRUE)
  png_file <- file.path(OUTPUT_DIR, paste0(filename_stem, ".png"))
  tiff_file <- file.path(OUTPUT_DIR, paste0(filename_stem, ".tiff"))
  pdf_file <- file.path(OUTPUT_DIR, paste0(filename_stem, ".pdf"))
  ggplot2::ggsave(png_file, plot_obj, width = width, height = height, dpi = FIG_DPI, bg = FIG_BACKGROUND, limitsize = FALSE)
  ggplot2::ggsave(tiff_file, plot_obj, width = width, height = height, dpi = FIG_DPI, compression = "lzw", bg = FIG_BACKGROUND, limitsize = FALSE)
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(pdf_file, plot_obj, width = width, height = height, device = pdf_device, bg = FIG_BACKGROUND, limitsize = FALSE)
  invisible(list(png = png_file, tiff = tiff_file, pdf = pdf_file))
}

find_first_existing <- function(candidates, label) {
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) stop(label, " was not found. Checked:\n", paste(candidates, collapse = "\n"))
  normalizePath(existing[1], winslash = "/", mustWork = TRUE)
}

cluster_number <- function(x) {
  x <- trimws(as.character(x))
  value <- suppressWarnings(as.integer(sub("^.*?([0-9]+)\\s*$", "\\1", x)))
  value[!grepl("[0-9]+\\s*$", x)] <- NA_integer_
  value
}

normalise_cluster_label <- function(x) {
  x <- trimws(as.character(x))
  id <- cluster_number(x)
  x[!is.na(id)] <- paste0("Cluster ", id[!is.na(id)])
  x
}

make_group_factor <- function(x) {
  raw <- as.character(x)
  raw_lower <- tolower(trimws(raw))
  if (all(raw %in% c("1", "2"))) return(factor(raw, levels = c("1", "2"), labels = c("Normal", "Tumor")))
  known <- c("normal", "tumor", "non", "lscc", "cancer", "margin")
  if (all(raw_lower %in% known)) {
    mapped <- ifelse(raw_lower %in% c("tumor", "lscc", "cancer"), "Tumor", "Normal")
    return(factor(mapped, levels = c("Normal", "Tumor")))
  }
  values <- sort(unique(raw))
  if (length(values) != 2L) stop("The group column must contain exactly two classes. Found: ", paste(values, collapse = ", "))
  factor(ifelse(raw == values[1], "Normal", "Tumor"), levels = c("Normal", "Tumor"))
}

read_bulk_cpm <- function(path, cohort_label) {
  dat <- as.data.frame(data.table::fread(path, check.names = FALSE), check.names = FALSE, stringsAsFactors = FALSE)
  missing_columns <- setdiff(c("Sample", "group"), colnames(dat))
  if (length(missing_columns) > 0L) stop(cohort_label, " is missing: ", paste(missing_columns, collapse = ", "))
  dat$Sample <- trimws(as.character(dat$Sample))
  dat$Group <- make_group_factor(dat$group)
  dat$Batch <- if ("batch" %in% colnames(dat)) as.character(dat$batch) else cohort_label
  gene_columns <- setdiff(colnames(dat), c("Sample", "group", "batch", "Group", "Batch"))
  if (length(gene_columns) < 100L) stop(cohort_label, " has too few gene columns.")
  expr <- as.data.frame(lapply(dat[, gene_columns, drop = FALSE], function(x) suppressWarnings(as.numeric(as.character(x)))), check.names = FALSE)
  colnames(expr) <- toupper(colnames(expr))
  expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]
  expr_matrix <- as.matrix(expr)
  storage.mode(expr_matrix) <- "numeric"
  expr_matrix[!is.finite(expr_matrix)] <- 0
  if (any(expr_matrix < 0, na.rm = TRUE)) stop(cohort_label, " contains negative values.")
  expr_gene_by_sample <- t(log2(expr_matrix + 1))
  rownames(expr_gene_by_sample) <- colnames(expr)
  colnames(expr_gene_by_sample) <- dat$Sample
  list(
    metadata = data.frame(Sample = dat$Sample, Group = dat$Group, Batch = dat$Batch, stringsAsFactors = FALSE),
    expression = expr_gene_by_sample,
    source = path
  )
}

run_ssgsea_safe <- function(expression_gene_by_sample, gene_sets, min_size = MIN_SIGNATURE_GENES) {
  expression_gene_by_sample <- as.matrix(expression_gene_by_sample)
  storage.mode(expression_gene_by_sample) <- "numeric"
  expression_gene_by_sample[!is.finite(expression_gene_by_sample)] <- 0
  new_result <- tryCatch({
    param <- GSVA::ssgseaParam(exprData = expression_gene_by_sample, geneSets = gene_sets, minSize = min_size, maxSize = 500, normalize = TRUE)
    GSVA::gsva(param, verbose = FALSE)
  }, error = function(e) e)
  if (!inherits(new_result, "error")) return(as.matrix(new_result))
  old_result <- tryCatch({
    GSVA::gsva(expr = expression_gene_by_sample, gset.idx.list = gene_sets, method = "ssgsea", kcdf = "Gaussian", min.sz = min_size, max.sz = 500, ssgsea.norm = TRUE, verbose = FALSE)
  }, error = function(e) e)
  if (inherits(old_result, "error")) stop("ssGSEA failed. New API: ", new_result$message, " Legacy API: ", old_result$message)
  as.matrix(old_result)
}

safe_wilcox <- function(x, group) {
  keep <- is.finite(x) & !is.na(group)
  x <- x[keep]
  group <- droplevels(factor(group[keep]))
  if (length(levels(group)) < 2L || min(table(group)) < 2L) return(NA_real_)
  tryCatch(stats::wilcox.test(x ~ group, exact = FALSE)$p.value, error = function(e) NA_real_)
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5L) return(data.frame(rho = NA_real_, p_value = NA_real_, n = sum(keep)))
  out <- suppressWarnings(stats::cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  data.frame(rho = as.numeric(out$estimate), p_value = as.numeric(out$p.value), n = sum(keep))
}

format_fdr <- function(x) {
  if (!is.finite(x)) return("FDR = NA")
  if (x < 0.001) return("FDR < 0.001")
  paste0("FDR = ", sprintf("%.3f", x))
}

safe_annotation_y <- function(y) {
  y <- y[is.finite(y)]
  if (length(y) == 0L) return(1)
  span <- diff(range(y))
  if (!is.finite(span) || span <= 0) span <- max(abs(y), 1) * 0.15
  max(y) + 0.10 * span
}

make_header <- function(label, text_colour = "black") {
  patchwork::wrap_elements(
    full = grid::textGrob(
      label = label, x = 0.5, y = 0.5, just = "centre",
      gp = grid::gpar(fontfamily = FONT_FAMILY, fontface = "bold", fontsize = 11, col = text_colour)
    )
  )
}

write_gmt <- function(gene_sets, file_path) {
  lines <- vapply(names(gene_sets), function(nm) paste(c(nm, "Marker_signature_without_MYBL2", gene_sets[[nm]]), collapse = "\t"), character(1))
  writeLines(lines, file_path)
  invisible(file_path)
}

# Step 4: Build signatures
build_signatures <- function(marker_csv, expression_gene_names) {
  marker_table <- read.csv(marker_csv, stringsAsFactors = FALSE, check.names = FALSE)
  required_columns <- c("Gene", "cluster", "p_val_adj", "pct.1")
  missing_columns <- setdiff(required_columns, colnames(marker_table))
  if (length(missing_columns) > 0L) stop("Marker CSV is missing: ", paste(missing_columns, collapse = ", "))
  fc_col <- if ("avg_log2FC" %in% colnames(marker_table)) "avg_log2FC" else if ("avg_logFC" %in% colnames(marker_table)) "avg_logFC" else NA_character_
  if (is.na(fc_col)) stop("Marker CSV needs avg_log2FC or avg_logFC.")
  markers <- data.frame(
    Gene = toupper(trimws(as.character(marker_table$Gene))),
    Subcluster = normalise_cluster_label(marker_table$cluster),
    p_val_adj = suppressWarnings(as.numeric(marker_table$p_val_adj)),
    avg_logFC = suppressWarnings(as.numeric(marker_table[[fc_col]])),
    pct.1 = suppressWarnings(as.numeric(marker_table$pct.1)),
    stringsAsFactors = FALSE
  )
  markers <- markers %>%
    dplyr::filter(
      Subcluster %in% SELECTED_SUBCLUSTERS,
      !is.na(Gene), nzchar(Gene), Gene != TARGET_GENE,
      is.finite(p_val_adj), p_val_adj < MARKER_PADJ_CUTOFF,
      is.finite(avg_logFC), avg_logFC >= MARKER_LOGFC_CUTOFF,
      is.finite(pct.1), pct.1 >= MARKER_PCT_CUTOFF
    )
  if (REMOVE_MITOCHONDRIAL_GENES) markers <- markers %>% dplyr::filter(!grepl("^MT-", Gene))
  if (REMOVE_RIBOSOMAL_GENES) markers <- markers %>% dplyr::filter(!grepl("^RPL|^RPS", Gene))
  markers <- markers %>%
    dplyr::filter(Gene %in% expression_gene_names) %>%
    dplyr::group_by(Subcluster, Gene) %>%
    dplyr::slice_min(order_by = p_val_adj, n = 1L, with_ties = FALSE) %>%
    dplyr::ungroup()
  gene_cluster_counts <- markers %>% dplyr::distinct(Subcluster, Gene) %>% dplyr::count(Gene, name = "n_clusters")
  signature_table <- markers %>%
    dplyr::left_join(gene_cluster_counts, by = "Gene") %>%
    dplyr::filter(n_clusters == 1L) %>%
    dplyr::arrange(Subcluster, p_val_adj, dplyr::desc(avg_logFC), dplyr::desc(pct.1), Gene) %>%
    dplyr::group_by(Subcluster) %>%
    dplyr::slice_head(n = TOP_SIGNATURE_GENES) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Subcluster) %>%
    dplyr::mutate(Signature_rank = dplyr::row_number(), Signature = paste0(Subcluster, " signature")) %>%
    dplyr::ungroup() %>%
    dplyr::select(Signature, Subcluster, Signature_rank, Gene, p_val_adj, avg_logFC, pct.1)
  counts <- signature_table %>% dplyr::count(Signature, Subcluster, name = "n_genes")
  missing_subclusters <- setdiff(SELECTED_SUBCLUSTERS, counts$Subcluster)
  small_subclusters <- counts$Subcluster[counts$n_genes < MIN_SIGNATURE_GENES]
  if (length(missing_subclusters) > 0L || length(small_subclusters) > 0L) {
    stop("Too few eligible genes. Missing: ", paste(missing_subclusters, collapse = ", "), "; small: ", paste(small_subclusters, collapse = ", "))
  }
  gene_sets <- split(signature_table$Gene, signature_table$Signature)
  gene_sets <- lapply(gene_sets, unique)
  gene_sets <- gene_sets[paste0(SELECTED_SUBCLUSTERS, " signature")]
  list(signature_table = signature_table, gene_sets = gene_sets)
}

# Step 5: Test scores
make_group_tests <- function(scores_long, cohort_label) {
  scores_long %>%
    dplyr::group_by(Signature) %>%
    dplyr::group_modify(~ {
      normal_values <- .x$Score[.x$Group == "Normal"]
      tumor_values <- .x$Score[.x$Group == "Tumor"]
      data.frame(
        Cohort = cohort_label,
        n_Normal = sum(.x$Group == "Normal"),
        n_Tumor = sum(.x$Group == "Tumor"),
        Median_Normal = stats::median(normal_values, na.rm = TRUE),
        Median_Tumor = stats::median(tumor_values, na.rm = TRUE),
        Median_difference_Tumor_minus_Normal = stats::median(tumor_values, na.rm = TRUE) - stats::median(normal_values, na.rm = TRUE),
        Wilcoxon_P_value = safe_wilcox(.x$Score, .x$Group),
        stringsAsFactors = FALSE
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(BH_FDR = stats::p.adjust(Wilcoxon_P_value, method = "BH"))
}

make_correlation_tests <- function(scores_long) {
  scores_long %>%
    dplyr::filter(Group == "Tumor") %>%
    dplyr::group_by(Signature) %>%
    dplyr::group_modify(~ {
      corr <- safe_spearman(.x$Score, .x$MYBL2_expression)
      data.frame(n_Tumor = corr$n, Spearman_rho = corr$rho, P_value = corr$p_value, stringsAsFactors = FALSE)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(BH_FDR = stats::p.adjust(P_value, method = "BH"))
}

# Step 6: Plot figures
make_boxplot <- function(scores_long, test_table, signature_name, panel_tag = NULL, y_axis_title = "ssGSEA score", panel_fill = FIG_BACKGROUND, panel_border_colour = "black") {
  data_one <- scores_long %>% dplyr::filter(Signature == signature_name)
  test_one <- test_table %>% dplyr::filter(Signature == signature_name)
  if (nrow(data_one) == 0L || nrow(test_one) != 1L) stop("Cannot create box plot for ", signature_name, ".")
  annotation <- paste0("\u0394median = ", sprintf("%.2f", test_one$Median_difference_Tumor_minus_Normal), "\n", format_fdr(test_one$BH_FDR))
  ggplot2::ggplot(data_one, ggplot2::aes(x = Group, y = Score, colour = Group)) +
    ggplot2::geom_boxplot(width = 0.48, fill = "white", outlier.shape = NA, linewidth = 0.85) +
    ggplot2::geom_jitter(width = 0.10, height = 0, size = 1.55, alpha = 0.68) +
    ggplot2::annotate("text", x = 1.5, y = safe_annotation_y(data_one$Score), label = annotation, family = FONT_FAMILY, size = 2.9, lineheight = 0.95) +
    ggplot2::scale_colour_manual(values = c("Normal" = COL_NORMAL, "Tumor" = COL_TUMOR)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.18))) +
    ggplot2::labs(title = signature_name, tag = panel_tag, x = NULL, y = y_axis_title, colour = NULL) +
    theme_manuscript(show_grid = FALSE, legend_position = "none", panel_fill = panel_fill, panel_border_colour = panel_border_colour)
}

make_correlation_plot <- function(scores_long, test_table, signature_name, panel_tag = NULL, point_colour = COL_CLUSTER0, y_axis_title = "ssGSEA score") {
  data_one <- scores_long %>% dplyr::filter(Group == "Tumor", Signature == signature_name)
  test_one <- test_table %>% dplyr::filter(Signature == signature_name)
  if (nrow(data_one) == 0L || nrow(test_one) != 1L) stop("Cannot create correlation plot for ", signature_name, ".")
  x_range <- range(data_one$MYBL2_expression, na.rm = TRUE)
  y_range <- range(data_one$Score, na.rm = TRUE)
  x_span <- diff(x_range)
  y_span <- diff(y_range)
  if (!is.finite(x_span) || x_span <= 0) x_span <- 1
  if (!is.finite(y_span) || y_span <= 0) y_span <- 1
  annotation <- paste0("\u03c1 = ", sprintf("%.2f", test_one$Spearman_rho), "\n", format_fdr(test_one$BH_FDR))
  ggplot2::ggplot(data_one, ggplot2::aes(x = MYBL2_expression, y = Score)) +
    ggplot2::geom_point(color = point_colour, size = 1.70, alpha = 0.78) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.78, color = COL_FIT, fill = COL_BAND, alpha = 0.65) +
    ggplot2::annotate("text", x = x_range[1] + 0.05 * x_span, y = y_range[2] - 0.04 * y_span, label = annotation, hjust = 0, vjust = 1, family = FONT_FAMILY, size = 2.95, lineheight = 0.95) +
    ggplot2::labs(title = signature_name, tag = panel_tag, x = NULL, y = y_axis_title) +
    theme_manuscript(show_grid = TRUE, legend_position = "none")
}

# Step 7: Run analysis
run_analysis <- function() {
  discovery_cpm_file <- find_first_existing(DEVELOPMENT_CPM_CANDIDATES, "Discovery CPM matrix")
  external_cpm_file <- find_first_existing(EXTERNAL_CPM_CANDIDATES, "External CPM matrix")
  marker_csv <- find_first_existing(MARKER_CSV_CANDIDATES, "High-CNV marker CSV")
  
  discovery <- read_bulk_cpm(discovery_cpm_file, "Discovery")
  external <- read_bulk_cpm(external_cpm_file, "External_GSE130605")
  
  if (!(TARGET_GENE %in% rownames(discovery$expression))) stop(TARGET_GENE, " is absent from discovery matrix.")
  if (!(TARGET_GENE %in% rownames(external$expression))) stop(TARGET_GENE, " is absent from external matrix.")
  
  shared_genes <- intersect(rownames(discovery$expression), rownames(external$expression))
  signature_object <- build_signatures(marker_csv, shared_genes)
  signature_table <- signature_object$signature_table
  gene_sets <- signature_object$gene_sets
  signature_columns <- names(gene_sets)
  
  write.csv(signature_table, SIG_GENES_CSV, row.names = FALSE)
  write_gmt(gene_sets, SIG_GMT)
  
  discovery_gene_sets <- lapply(gene_sets, intersect, y = rownames(discovery$expression))
  external_gene_sets <- lapply(gene_sets, intersect, y = rownames(external$expression))
  
  discovery_score_matrix <- run_ssgsea_safe(discovery$expression, discovery_gene_sets)
  external_score_matrix <- run_ssgsea_safe(external$expression, external_gene_sets)
  
  discovery_scores <- as.data.frame(t(discovery_score_matrix), check.names = FALSE)
  external_scores <- as.data.frame(t(external_score_matrix), check.names = FALSE)
  discovery_scores$Sample <- rownames(discovery_scores)
  external_scores$Sample <- rownames(external_scores)
  
  discovery_scores <- dplyr::left_join(discovery$metadata, discovery_scores, by = "Sample")
  external_scores <- dplyr::left_join(external$metadata, external_scores, by = "Sample")
  discovery_scores$MYBL2_expression <- as.numeric(discovery$expression[TARGET_GENE, discovery_scores$Sample])
  external_scores$MYBL2_expression <- as.numeric(external$expression[TARGET_GENE, external_scores$Sample])
  
  discovery_scores <- discovery_scores[, c("Sample", "Group", "Batch", "MYBL2_expression", signature_columns), drop = FALSE]
  external_scores <- external_scores[, c("Sample", "Group", "Batch", "MYBL2_expression", signature_columns), drop = FALSE]
  write.csv(discovery_scores, DISC_SCORES_CSV, row.names = FALSE)
  write.csv(external_scores, EXT_SCORES_CSV, row.names = FALSE)
  
  discovery_long <- discovery_scores %>%
    tidyr::pivot_longer(cols = dplyr::all_of(signature_columns), names_to = "Signature", values_to = "Score") %>%
    dplyr::mutate(Cohort = "Discovery", Group = factor(Group, levels = c("Normal", "Tumor")))
  external_long <- external_scores %>%
    tidyr::pivot_longer(cols = dplyr::all_of(signature_columns), names_to = "Signature", values_to = "Score") %>%
    dplyr::mutate(Cohort = "External GSE130605", Group = factor(Group, levels = c("Normal", "Tumor")))
  
  discovery_tests <- make_group_tests(discovery_long, "Discovery")
  external_tests <- make_group_tests(external_long, "External_GSE130605")
  discovery_corr <- make_correlation_tests(discovery_long)
  write.csv(discovery_tests, DISC_TESTS_CSV, row.names = FALSE)
  write.csv(external_tests, EXT_TESTS_CSV, row.names = FALSE)
  write.csv(discovery_corr, DISC_CORR_CSV, row.names = FALSE)
  
  cluster0 <- "Cluster 0 signature"
  cluster8 <- "Cluster 8 signature"
  final_figure <- (
    make_header(DISCOVERY_COHORT_HEADER) +
      make_boxplot(discovery_long, discovery_tests, cluster0, panel_tag = "a") +
      make_boxplot(discovery_long, discovery_tests, cluster8) +
      make_header(EXTERNAL_COHORT_HEADER, text_colour = COL_EXTERNAL_PANEL_BORDER) +
      make_boxplot(external_long, external_tests, cluster0, panel_tag = "b", panel_fill = COL_EXTERNAL_PANEL_BG, panel_border_colour = COL_EXTERNAL_PANEL_BORDER) +
      make_boxplot(external_long, external_tests, cluster8, panel_fill = COL_EXTERNAL_PANEL_BG, panel_border_colour = COL_EXTERNAL_PANEL_BORDER) +
      make_correlation_plot(discovery_long, discovery_corr, cluster0, panel_tag = "c", point_colour = COL_CLUSTER0) +
      make_correlation_plot(discovery_long, discovery_corr, cluster8, point_colour = COL_CLUSTER8)
  ) + patchwork::plot_layout(design = "AA\nBC\nDD\nEF\nGH", heights = c(0.10, 1, 0.10, 1, 1.08))
  
  save_plot_all_formats(final_figure, FINAL_FIGURE_STEM)
  invisible(TRUE)
}

run_analysis()
