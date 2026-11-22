required_packages <- c(
  "data.table", "ggplot2", "ggrepel", "clusterProfiler", "org.Hs.eg.db",
  "patchwork", "ggwordcloud", "GSEABase", "TissueEnrich", "SummarizedExperiment"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(patchwork)
  library(ggwordcloud)
  library(GSEABase)
  library(TissueEnrich)
  library(SummarizedExperiment)
})

options(timeout = 10000)
theme_set(theme_minimal(base_family = "Arial"))

OMICS_TYPE <- "clinical"
TARGET_ALLERGY <- "any_allergy"
COVARS <- c(
  "age", "sex", "ethnicity", "smoke", "alcohol", "IPAQ", "BMI", "TDI",
  "hypertension_status", "DM_status"
)
CATEGORICAL_COVARS <- c("sex", "ethnicity", "smoke", "alcohol", "IPAQ", "hypertension_status", "DM_status")
FDR_CUTOFF <- 0.05
BETA_CUTOFF <- 0
TOP_N_LABEL <- 20L
DO_MET_NORMALIZE <- FALSE
OUTPUT_DIR <- Sys.getenv("ANALYSIS_OUTPUT_DIR", unset = "aggregate_outputs")
RUN_PATHWAY_ENRICHMENT <- TRUE
RUN_TISSUE_ENRICHMENT <- TRUE
TISSUE_ENRICH_DATASET <- 1L
TISSUE_GENE_TYPE <- 1L
PATHWAY_TOP_KEEP <- 4L
MAX_LOGP <- 250
GENE_ID_TYPE <- "SYMBOL"

if (!OMICS_TYPE %in% c("metabolic", "proteomic", "clinical", "inflam")) {
  stop("OMICS_TYPE must be one of: metabolic, proteomic, clinical, inflam.", call. = FALSE)
}
if (!exists("analysis_inputs", envir = .GlobalEnv, inherits = FALSE)) {
  stop("Provide analysis_inputs in memory in the authorized analysis environment.", call. = FALSE)
}
if (!is.list(analysis_inputs) || is.null(analysis_inputs$baseline) || is.null(analysis_inputs$omics)) {
  stop("analysis_inputs must contain baseline and omics data objects.", call. = FALSE)
}

assert_columns <- function(x, required, object_name) {
  missing_columns <- setdiff(required, names(x))
  if (length(missing_columns) > 0L) {
    stop(object_name, " is missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

as_numeric_checked <- function(x, column_name) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    x[trimws(x) == ""] <- NA_character_
    converted <- suppressWarnings(as.numeric(x))
    invalid <- !is.na(x) & is.na(converted)
    if (any(invalid)) {
      stop("Non-numeric values found in ", column_name, call. = FALSE)
    }
    return(converted)
  }
  if (!is.numeric(x)) {
    stop(column_name, " must be numeric.", call. = FALSE)
  }
  as.numeric(x)
}

as_binary_checked <- function(x, column_name) {
  value <- as_numeric_checked(x, column_name)
  invalid <- !is.na(value) & !value %in% c(0, 1)
  if (any(invalid)) {
    stop(column_name, " must be coded 0/1/NA.", call. = FALSE)
  }
  as.integer(value)
}

validate_ids <- function(x, object_name) {
  if (anyNA(x$eid) || anyDuplicated(x$eid)) {
    stop(object_name, " must have unique, non-missing participant keys.", call. = FALSE)
  }
  invisible(TRUE)
}

prepare_inputs <- function(baseline, omics, exposure, covariates) {
  base <- as.data.table(copy(baseline))
  omic <- as.data.table(copy(omics))
  if (anyDuplicated(names(base)) || anyDuplicated(names(omic))) {
    stop("Input objects contain duplicated column names.", call. = FALSE)
  }
  assert_columns(base, c("eid", exposure, covariates), "analysis_inputs$baseline")
  assert_columns(omic, "eid", "analysis_inputs$omics")
  base[, eid := as.character(eid)]
  omic[, eid := as.character(eid)]
  validate_ids(base, "Baseline input")
  validate_ids(omic, "Omics input")
  base[, (exposure) := as_binary_checked(get(exposure), exposure)]
  for (v in CATEGORICAL_COVARS) {
    if (v %in% covariates) base[[v]] <- as.factor(base[[v]])
  }
  numeric_covars <- setdiff(covariates, CATEGORICAL_COVARS)
  for (v in numeric_covars) base[, (v) := as_numeric_checked(get(v), v)]
  feature_names <- setdiff(names(omic), c("eid", exposure, covariates))
  if (!length(feature_names)) stop("No omics features were found.", call. = FALSE)
  overlap <- intersect(feature_names, c(exposure, covariates))
  if (length(overlap)) stop("Omics feature names conflict with model variables.", call. = FALSE)
  for (v in feature_names) omic[, (v) := as_numeric_checked(get(v), v)]
  omic <- omic[, c("eid", feature_names), with = FALSE]
  pheno <- base[, c("eid", exposure, covariates), with = FALSE]
  pheno <- pheno[complete.cases(pheno)]
  if (length(unique(pheno[[exposure]])) != 2L) {
    stop("Both exposure groups must be present after baseline complete-case filtering.", call. = FALSE)
  }
  merged <- merge(pheno, omic, by = "eid", all = FALSE, sort = FALSE)
  if (!nrow(merged)) stop("No overlapping participants between baseline and omics inputs.", call. = FALSE)
  list(data = merged, features = feature_names)
}

normalize_metabolomics <- function(data, features) {
  matrix <- as.matrix(data[, ..features])
  storage.mode(matrix) <- "double"
  reference <- apply(matrix, 2L, stats::median, na.rm = TRUE)
  reference[!is.finite(reference) | reference == 0] <- NA_real_
  quotient <- sweep(matrix, 2L, reference, "/")
  factors <- apply(quotient, 1L, stats::median, na.rm = TRUE)
  factors[!is.finite(factors) | factors <= 0] <- NA_real_
  normalized <- sweep(matrix, 1L, factors, "/")
  if (any(normalized < 0, na.rm = TRUE)) {
    stop("PQN-normalized metabolite values contain negatives; log10 transformation is invalid.", call. = FALSE)
  }
  logged <- log10(normalized + 1e-6)
  scaled <- scale(logged)
  scaled[, apply(scaled, 2L, function(z) all(!is.finite(z)))] <- NA_real_
  data[, (features) := as.data.table(scaled)]
  data
}

fit_feature_model <- function(data, feature, exposure, covariates) {
  required <- c(exposure, covariates, feature)
  d <- copy(data[, ..required])
  setnames(d, feature, "feature_value")
  d <- d[complete.cases(d)]
  d <- d[is.finite(feature_value)]
  if (nrow(d) <= length(covariates) + 5L) {
    return(data.table(feature = feature, beta = NA_real_, se = NA_real_, p_value = NA_real_, n_sample = nrow(d), model_status = "insufficient_sample"))
  }
  if (length(unique(d$feature_value)) < 2L || stats::sd(d$feature_value) == 0) {
    return(data.table(feature = feature, beta = NA_real_, se = NA_real_, p_value = NA_real_, n_sample = nrow(d), model_status = "zero_variance"))
  }
  formula <- stats::reformulate(c("feature_value", covariates), response = exposure)
  fit <- tryCatch(
    suppressWarnings(stats::glm(formula, data = d, family = stats::binomial(link = "logit"))),
    error = function(e) NULL
  )
  if (is.null(fit) || !isTRUE(fit$converged)) {
    return(data.table(feature = feature, beta = NA_real_, se = NA_real_, p_value = NA_real_, n_sample = nrow(d), model_status = "model_failed_or_not_converged"))
  }
  coefficient_table <- summary(fit)$coefficients
  if (!"feature_value" %in% rownames(coefficient_table)) {
    return(data.table(feature = feature, beta = NA_real_, se = NA_real_, p_value = NA_real_, n_sample = nrow(d), model_status = "feature_coefficient_unavailable"))
  }
  beta <- unname(coefficient_table["feature_value", "Estimate"])
  standard_error <- unname(coefficient_table["feature_value", "Std. Error"])
  p_value <- unname(coefficient_table["feature_value", "Pr(>|z|)"])
  if (!is.finite(beta) || !is.finite(standard_error) || !is.finite(p_value)) {
    return(data.table(feature = feature, beta = NA_real_, se = NA_real_, p_value = NA_real_, n_sample = nrow(d), model_status = "non_finite_estimate"))
  }
  data.table(feature = feature, beta = beta, se = standard_error, p_value = p_value, n_sample = nrow(d), model_status = "ok")
}

make_volcano_plot <- function(results, omics_type, exposure, fdr_cutoff, top_n, max_logp, output_dir) {
  plot_data <- copy(results)
  plot_data[, log_fdr_plot := pmin(log_fdr, max_logp)]
  plot_data[!is.finite(log_fdr_plot), log_fdr_plot := max_logp]
  plot_data <- plot_data[is.finite(beta) & is.finite(log_fdr_plot)]
  if (!nrow(plot_data)) return(invisible(NULL))
  significant_labels <- plot_data[group != "Not significant"][order(fdr, -abs(beta))]
  labels <- significant_labels[seq_len(min(top_n, nrow(significant_labels)))]
  beta_range <- range(plot_data$beta, finite = TRUE)
  x_span <- diff(beta_range)
  if (!is.finite(x_span) || x_span == 0) x_span <- max(abs(beta_range), 1)
  y_max <- max(plot_data$log_fdr_plot, na.rm = TRUE)
  y_max <- ifelse(is.finite(y_max) && y_max > 0, y_max, 1)
  plot <- ggplot(plot_data, aes(x = beta, y = log_fdr_plot, color = group)) +
    geom_point(size = 2.4, alpha = 0.75) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray35", linewidth = 0.45) +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "gray35", linewidth = 0.45) +
    geom_point(data = labels, aes(fill = group), shape = 21, size = 4, colour = "black", stroke = 0.45, show.legend = FALSE) +
    geom_label_repel(data = labels, aes(label = feature), color = "black", size = 3, box.padding = 0.4,
                     point.padding = 0.25, segment.color = "gray35", max.overlaps = 50,
                     fill = grDevices::adjustcolor("white", alpha.f = 0.85), show.legend = FALSE) +
    annotate("text", x = beta_range[1] + 0.02 * x_span, y = 0.82 * y_max,
             label = paste(sum(plot_data$group == "Negative association"), "negative"),
             color = "#8FA8C8", hjust = 0, size = 4) +
    annotate("text", x = beta_range[2] - 0.02 * x_span, y = 0.82 * y_max,
             label = paste(sum(plot_data$group == "Positive association"), "positive"),
             color = "#C97C7C", hjust = 1, size = 4) +
    scale_color_manual(values = c("Negative association" = "#8FA8C8", "Not significant" = "#CCCCCC", "Positive association" = "#C97C7C")) +
    scale_fill_manual(values = c("Negative association" = "#8FA8C8", "Positive association" = "#C97C7C")) +
    labs(x = "Beta (log odds ratio)", y = "-log10(FDR)",
         title = paste(tools::toTitleCase(omics_type), "association with", exposure),
         color = "Association") +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = c(0.86, 0.9),
          plot.title = element_text(hjust = 0.5), panel.border = element_rect(colour = "black", fill = NA))
  ggplot2::ggsave(file.path(output_dir, paste0(exposure, "_", omics_type, "_volcano.pdf")), plot,
                  width = 8, height = 8)
  invisible(plot)
}

parse_ratio <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)[[1L]]
  if (length(parts) != 2L) return(NA_real_)
  values <- suppressWarnings(as.numeric(parts))
  if (any(!is.finite(values)) || values[1L] <= 0 || values[2L] <= 0) return(NA_real_)
  values
}

standardize_enrichment <- function(enrichment, database) {
  if (is.null(enrichment)) return(data.table())
  result <- as.data.table(as.data.frame(enrichment))
  if (!nrow(result)) return(data.table())
  if ("p.adjust" %in% names(result)) setnames(result, "p.adjust", "padj")
  if (!"padj" %in% names(result) || !"Description" %in% names(result)) return(data.table())
  result[, padj := as.numeric(padj)]
  if (!"Count" %in% names(result)) result[, Count := NA_real_]
  if (!"BgRatio" %in% names(result)) result[, BgRatio := NA_character_]
  result[, bg_count := vapply(BgRatio, function(z) parse_ratio(z)[1L], numeric(1))]
  result[, RichFactor := ifelse(is.finite(bg_count) & bg_count > 0, Count / bg_count, NA_real_)]
  result[, pval_log := -log10(pmax(padj, .Machine$double.xmin))]
  result[, Database := database]
  result
}

map_gene_symbols <- function(gene_symbols, universe_symbols) {
  gene_symbols <- unique(toupper(trimws(as.character(gene_symbols))))
  universe_symbols <- unique(toupper(trimws(as.character(universe_symbols))))
  gene_symbols <- gene_symbols[!is.na(gene_symbols) & nzchar(gene_symbols)]
  universe_symbols <- universe_symbols[!is.na(universe_symbols) & nzchar(universe_symbols)]
  if (!length(gene_symbols) || !length(universe_symbols)) return(NULL)
  mapping <- suppressMessages(clusterProfiler::bitr(
    unique(c(gene_symbols, universe_symbols)), fromType = GENE_ID_TYPE,
    toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db, drop = TRUE
  ))
  if (is.null(mapping) || !nrow(mapping)) return(NULL)
  mapping <- unique(as.data.table(mapping))
  gene_ids <- unique(mapping[toupper(get(GENE_ID_TYPE)) %in% gene_symbols, ENTREZID])
  universe_ids <- unique(mapping[toupper(get(GENE_ID_TYPE)) %in% universe_symbols, ENTREZID])
  gene_ids <- gene_ids[!is.na(gene_ids)]
  universe_ids <- universe_ids[!is.na(universe_ids)]
  gene_ids <- intersect(gene_ids, universe_ids)
  if (!length(gene_ids) || !length(universe_ids)) return(NULL)
  list(gene = gene_ids, universe = universe_ids)
}

run_pathway_enrichment <- function(gene_symbols, universe_symbols, regulation, output_dir, exposure) {
  mapped <- map_gene_symbols(gene_symbols, universe_symbols)
  if (is.null(mapped)) {
    message("No mapped genes for ", regulation, " pathway enrichment.")
    return(list())
  }
  enrichment_objects <- list(
    GOBP = tryCatch(clusterProfiler::enrichGO(gene = mapped$gene, universe = mapped$universe,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
      pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", readable = TRUE), error = function(e) NULL),
    GOMF = tryCatch(clusterProfiler::enrichGO(gene = mapped$gene, universe = mapped$universe,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID", ont = "MF",
      pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", readable = TRUE), error = function(e) NULL),
    KEGG = tryCatch(clusterProfiler::enrichKEGG(gene = mapped$gene, universe = mapped$universe,
      organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH"), error = function(e) NULL),
    WikiPathways = tryCatch(clusterProfiler::enrichWP(gene = mapped$gene, universe = mapped$universe,
      organism = "Homo sapiens", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH"), error = function(e) NULL)
  )
  if (requireNamespace("ReactomePA", quietly = TRUE)) {
    enrichment_objects$Reactome <- tryCatch(ReactomePA::enrichPathway(
      gene = mapped$gene, universe = mapped$universe, organism = "human",
      pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", readable = TRUE
    ), error = function(e) NULL)
  }
  standardized <- list()
  for (database in names(enrichment_objects)) {
    object <- enrichment_objects[[database]]
    result <- standardize_enrichment(object, database)
    standardized[[database]] <- result
    output_file <- file.path(output_dir, paste0(exposure, "_proteomic_", regulation, "_", database, "_result.csv"))
    fwrite(result, output_file)
    if (nrow(result)) {
      plot_data <- result[is.finite(padj) & padj < FDR_CUTOFF][order(padj)][seq_len(min(15L, .N))]
      if (nrow(plot_data)) {
        dot <- tryCatch(clusterProfiler::dotplot(object, showCategory = min(15L, nrow(plot_data)),
                                                title = paste(database, regulation)), error = function(e) NULL)
        if (!is.null(dot)) {
          ggsave(file.path(output_dir, paste0(exposure, "_proteomic_", regulation, "_", database, "_dot.pdf")),
                 dot, width = 9, height = 7)
        }
      }
    }
  }
  standardized
}

read_top_pathways <- function(enrichment_list, top_k, fdr_cutoff) {
  selected <- lapply(names(enrichment_list), function(database) {
    dt <- copy(enrichment_list[[database]])
    if (!nrow(dt)) return(NULL)
    dt <- dt[is.finite(padj) & padj < fdr_cutoff][order(padj)]
    if (!nrow(dt)) return(NULL)
    dt <- dt[seq_len(min(top_k, nrow(dt)))]
    dt[, PathLabel := paste0(Database, ": ", Description)]
    dt
  })
  selected <- Filter(Negate(is.null), selected)
  if (!length(selected)) return(data.table())
  combined <- rbindlist(selected, fill = TRUE)
  combined[, Database := factor(Database, levels = c("GOBP", "GOMF", "KEGG", "Reactome", "WikiPathways"))]
  setorder(combined, Database, padj)
  combined
}

make_wordcloud_plot <- function(pathway_data, seed = 123L) {
  if (is.null(pathway_data) || !nrow(pathway_data)) return(ggplot() + theme_void())
  words <- tolower(unlist(strsplit(paste(pathway_data$Description, collapse = " "), "[^A-Za-z0-9]+")))
  stop_words <- c("the", "of", "and", "in", "to", "by", "via", "for", "with", "from", "its", "their",
                  "a", "is", "are", "was", "were", "upon", "process", "regulation", "positive", "negative",
                  "cell", "cellular", "protein", "mediated", "signal", "signaling", "pathway", "system", "cytokine")
  words <- words[nchar(words) > 2L & !words %in% stop_words]
  frequency <- as.data.table(as.data.frame(table(words)))
  setnames(frequency, c("words", "Freq"), c("word", "freq"))
  frequency[, freq := as.numeric(freq)]
  frequency <- frequency[freq >= 2]
  if (!nrow(frequency)) return(ggplot() + theme_void())
  set.seed(seed)
  ggplot(frequency, aes(label = word, size = freq, color = word)) +
    ggwordcloud::geom_text_wordcloud_area(rm_outside = TRUE, shape = "circle") +
    scale_size_area(max_size = 12) +
    scale_color_brewer(palette = "Set2", guide = "none") +
    theme_void()
}

plot_pathway_summary <- function(pathway_data, regulation, output_dir, exposure) {
  if (is.null(pathway_data) || !nrow(pathway_data)) return(invisible(NULL))
  pathway_data <- copy(pathway_data)
  pathway_data[, PathLabel := factor(PathLabel, levels = rev(unique(PathLabel)))]
  pathway_data[, text_x := 0.03]
  colors_up <- c(GOBP = "#883838", GOMF = "#D05050", KEGG = "#F28040", Reactome = "#F8B080", WikiPathways = "#F9E0C8")
  colors_down <- c(GOBP = "#4070A8", GOMF = "#7098C8", KEGG = "#B0D0E8", Reactome = "#40A8B0", WikiPathways = "#8060A8")
  color_map <- if (regulation == "up") colors_up else colors_down
  bar <- ggplot(pathway_data, aes(x = pval_log, y = PathLabel, fill = Database)) +
    geom_col(width = 0.8, alpha = 0.85) +
    geom_text(aes(x = text_x, label = PathLabel), hjust = 0, size = 3.2, color = "black") +
    scale_fill_manual(values = color_map, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    coord_cartesian(clip = "off") +
    labs(x = "-log10(FDR)", y = NULL, title = paste(toupper(regulation), "pathways")) +
    theme_classic() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(10, 20, 10, 10), legend.position = "right")
  cloud <- make_wordcloud_plot(pathway_data)
  combined <- bar + cloud + plot_layout(widths = c(3, 1))
  ggsave(file.path(output_dir, paste0(exposure, "_proteomic_", regulation, "_pathway_summary.pdf")),
         combined, width = 10, height = 9)
  fwrite(pathway_data, file.path(output_dir, paste0(exposure, "_proteomic_", regulation, "_top_pathways.csv")))
  invisible(combined)
}

run_tissue_enrichment <- function(significant_symbols, universe_symbols, output_dir, exposure,
                                  dataset, gene_type) {
  significant_symbols <- unique(toupper(na.omit(as.character(significant_symbols))))
  universe_symbols <- unique(toupper(na.omit(as.character(universe_symbols))))
  significant_symbols <- intersect(significant_symbols, universe_symbols)
  if (!length(significant_symbols)) {
    message("No significant protein symbols available for tissue enrichment.")
    return(invisible(NULL))
  }
  input_set <- GSEABase::GeneSet(geneIds = significant_symbols, organism = "Homo Sapiens",
                                 geneIdType = GSEABase::SymbolIdentifier())
  background_set <- GSEABase::GeneSet(geneIds = universe_symbols, organism = "Homo Sapiens",
                                     geneIdType = GSEABase::SymbolIdentifier())
  tissue_result <- TissueEnrich::teEnrichment(
    inputGenes = input_set, rnaSeqDataset = dataset,
    tissueSpecificGeneType = gene_type, multiHypoCorrection = TRUE,
    backgroundGenes = background_set
  )
  if (is.null(tissue_result) || !length(tissue_result) || is.null(tissue_result[[1L]])) {
    message("Tissue enrichment returned no result.")
    return(invisible(NULL))
  }
  se <- tissue_result[[1L]]
  enrichment <- as.data.frame(SummarizedExperiment::assay(se))
  rownames(enrichment) <- as.character(SummarizedExperiment::rowData(se)[, 1L])
  colnames(enrichment) <- as.character(SummarizedExperiment::colData(se)[, 1L])
  enrichment$Tissue <- rownames(enrichment)
  result <- as.data.table(enrichment)
  if (!"Log10PValue" %in% names(result)) stop("TissueEnrich output does not contain Log10PValue.", call. = FALSE)
  result[, Log10PValue := as.numeric(Log10PValue)]
  result <- result[is.finite(Log10PValue)]
  significant <- result[Log10PValue > -log10(FDR_CUTOFF)]
  nonsignificant <- result[Log10PValue <= -log10(FDR_CUTOFF)][order(-Log10PValue)]
  if (nrow(nonsignificant) > 5L) nonsignificant <- nonsignificant[seq_len(5L)]
  result <- rbindlist(list(significant, nonsignificant), fill = TRUE)
  setorder(result, -Log10PValue)
  result[, id := .I]
  fwrite(result, file.path(output_dir, paste0(exposure, "_proteomic_tissue_enrichment.csv")))
  if (!nrow(result)) return(invisible(NULL))
  threshold <- -log10(FDR_CUTOFF)
  circle <- ggplot(result, aes(x = id, y = Log10PValue, fill = Log10PValue)) +
    geom_col(width = 0.72) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "#E03E3E") +
    scale_fill_gradientn(colours = c("#B7A1C9", "#8FA8C8", "#A8B89A", "#D7C9A4", "#D8B384", "#D6A5A4")) +
    coord_polar(start = 0) +
    labs(title = "Tissue-specific gene enrichment", y = "-log10(adjusted P-value)", x = NULL) +
    theme_void() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  ggsave(file.path(output_dir, paste0(exposure, "_proteomic_tissue_enrichment_circle.pdf")),
         circle, width = 10, height = 10, device = grDevices::cairo_pdf)
  bar <- ggplot(result, aes(x = reorder(Tissue, Log10PValue), y = Log10PValue, fill = Log10PValue)) +
    geom_col(width = 0.72) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
    scale_fill_gradientn(colours = c("#B7A1C9", "#8FA8C8", "#A8B89A", "#D7C9A4", "#D8B384", "#D6A5A4")) +
    coord_flip() +
    labs(title = "Tissue-specific gene enrichment", x = NULL, y = "-log10(adjusted P-value)") +
    theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
  ggsave(file.path(output_dir, paste0(exposure, "_proteomic_tissue_enrichment_bar.pdf")),
         bar, width = 9, height = 7)
  if ("fold.change" %in% names(result)) {
    fold_plot <- ggplot(result, aes(x = reorder(Tissue, fold.change), y = fold.change)) +
      geom_col(width = 0.72, fill = "#B9B0A3") + coord_flip() +
      labs(title = "Tissue enrichment fold change", x = NULL, y = "Fold change") +
      theme_bw() + theme(panel.grid = element_blank())
    ggsave(file.path(output_dir, paste0(exposure, "_proteomic_tissue_enrichment_fold_change.pdf")),
           fold_plot, width = 9, height = 7)
  }
  invisible(result)
}

omics_type <- match.arg(OMICS_TYPE, c("metabolic", "proteomic", "clinical", "inflam"))
omics_key <- switch(
  omics_type,
  proteomic = "proteomics",
  metabolic = "metabolomics",
  clinical = "clinical",
  inflam = "inflammation"
)
if (is.null(analysis_inputs$omics[[omics_key]])) {
  stop("Missing omics input for OMICS_TYPE: ", OMICS_TYPE, call. = FALSE)
}
prepared <- prepare_inputs(analysis_inputs$baseline, analysis_inputs$omics[[omics_key]], TARGET_ALLERGY, COVARS)
omics_data <- prepared$data
omics_features <- prepared$features

if (omics_type == "metabolic" && DO_MET_NORMALIZE) {
  omics_data <- normalize_metabolomics(omics_data, omics_features)
}

message("Analyzing ", length(omics_features), " features across ", nrow(omics_data), " overlapping participants.")
start_time <- Sys.time()
association_results <- rbindlist(lapply(seq_along(omics_features), function(i) {
  if (i %% 50L == 0L || i == length(omics_features)) {
    message("Processed ", i, "/", length(omics_features), " features.")
  }
  fit_feature_model(omics_data, omics_features[i], TARGET_ALLERGY, COVARS)
}), fill = TRUE)
message("Association analysis elapsed time: ", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1), " seconds.")

if (!any(!is.na(association_results$p_value))) stop("No valid feature association models were fitted.", call. = FALSE)
association_results[, fdr := NA_real_]
association_results[!is.na(p_value), fdr := stats::p.adjust(p_value, method = "BH")]
association_results[, log_fdr := -log10(fdr)]
association_results[, log_fdr_capped := fifelse(is.infinite(log_fdr), MAX_LOGP, log_fdr)]
association_results[, group := fcase(
  !is.na(fdr) & fdr < FDR_CUTOFF & beta > BETA_CUTOFF, "Positive association",
  !is.na(fdr) & fdr < FDR_CUTOFF & beta < -BETA_CUTOFF, "Negative association",
  default = "Not significant"
)]

if (any(!is.finite(association_results$beta)) || any(!is.finite(association_results$p_value))) {
  stop("Non-finite association estimates remain after model filtering.", call. = FALSE)
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
association_file <- file.path(OUTPUT_DIR, paste0(TARGET_ALLERGY, "_", omics_type, "_assoc_result.csv"))
fwrite(association_results, association_file)
significant_features <- unique(association_results[group != "Not significant", feature])
fwrite(data.table(feature = significant_features),
       file.path(OUTPUT_DIR, paste0(TARGET_ALLERGY, "_", omics_type, "_significant_features.csv")))
make_volcano_plot(association_results, omics_type, TARGET_ALLERGY, FDR_CUTOFF,
                  TOP_N_LABEL, MAX_LOGP, OUTPUT_DIR)

if (omics_type == "proteomic" && RUN_PATHWAY_ENRICHMENT) {
  if (GENE_ID_TYPE != "SYMBOL") stop("Pathway enrichment currently expects gene symbols.", call. = FALSE)
  tested_symbols <- unique(association_results[!is.na(p_value), feature])
  up_symbols <- unique(association_results[group == "Positive association", feature])
  down_symbols <- unique(association_results[group == "Negative association", feature])
  fwrite(data.table(feature = up_symbols), file.path(OUTPUT_DIR, paste0(TARGET_ALLERGY, "_proteomic_up_features.csv")))
  fwrite(data.table(feature = down_symbols), file.path(OUTPUT_DIR, paste0(TARGET_ALLERGY, "_proteomic_down_features.csv")))
  enrichment_up <- run_pathway_enrichment(up_symbols, tested_symbols, "up", OUTPUT_DIR, TARGET_ALLERGY)
  enrichment_down <- run_pathway_enrichment(down_symbols, tested_symbols, "down", OUTPUT_DIR, TARGET_ALLERGY)
  top_up <- read_top_pathways(enrichment_up, PATHWAY_TOP_KEEP, FDR_CUTOFF)
  top_down <- read_top_pathways(enrichment_down, PATHWAY_TOP_KEEP, FDR_CUTOFF)
  plot_pathway_summary(top_up, "up", OUTPUT_DIR, TARGET_ALLERGY)
  plot_pathway_summary(top_down, "down", OUTPUT_DIR, TARGET_ALLERGY)
  combined_summary <- rbindlist(list(
    if (nrow(top_up)) copy(top_up)[, Regulation := "up"],
    if (nrow(top_down)) copy(top_down)[, Regulation := "down"]
  ), fill = TRUE)
  if (nrow(combined_summary)) {
    fwrite(combined_summary, file.path(OUTPUT_DIR, paste0(TARGET_ALLERGY, "_proteomic_pathway_top", PATHWAY_TOP_KEEP, "_summary.csv")))
  }
}

if (omics_type == "proteomic" && RUN_TISSUE_ENRICHMENT) {
  run_tissue_enrichment(
    significant_symbols = significant_features,
    universe_symbols = unique(association_results[!is.na(p_value), feature]),
    output_dir = OUTPUT_DIR,
    exposure = TARGET_ALLERGY,
    dataset = TISSUE_ENRICH_DATASET,
    gene_type = TISSUE_GENE_TYPE
  )
}

message("Analysis completed. Association results: ", nrow(association_results),
        "; significant features: ", length(significant_features), ".")
