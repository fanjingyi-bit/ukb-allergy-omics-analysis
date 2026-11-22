required_packages <- c("data.table", "glmnet", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Install required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(ggplot2)
})

OMICS_TYPE <- "proteomics"
EXPOSURE <- "any_allergy"
COVARIATES <- c("age", "sex", "BMI", "smoke", "alcohol", "ethnicity", "TDI", "IPAQ", "hypertension_status", "DM_status")
NUMERIC_COVARIATES <- c("age", "BMI", "TDI")
FACTOR_COVARIATES <- c("sex", "smoke", "alcohol", "ethnicity", "IPAQ", "hypertension_status", "DM_status")
TRAIN_FRACTION <- 0.70
RANDOM_SEED <- 123L
N_FOLDS <- 10L
MAX_PROTEIN_CANDIDATES <- 200L
PROTEIN_BONFERRONI_ALPHA <- 0.05
PROTEIN_MIN_ABS_COEF <- 0.01
METABOLOMICS_FDR_THRESHOLD <- 0.05
OUTPUT_DIR <- Sys.getenv("ANALYSIS_OUTPUT_DIR", unset = "score_outputs")

assert_columns <- function(x, required, object_name) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(object_name, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

as_numeric_strict <- function(x, object_name) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    x[trimws(x) == ""] <- NA_character_
    y <- suppressWarnings(as.numeric(x))
    if (any(!is.na(x) & is.na(y))) stop(object_name, " contains non-numeric values", call. = FALSE)
    return(y)
  }
  if (!is.numeric(x)) stop(object_name, " must be numeric", call. = FALSE)
  as.numeric(x)
}

validate_ids <- function(x, object_name) {
  if (anyNA(x$eid) || anyDuplicated(x$eid)) stop(object_name, " must have unique non-missing participant keys", call. = FALSE)
  invisible(TRUE)
}

make_feature_formula <- function(response, feature, covariates) {
  terms <- c(feature, covariates)
  if (any(grepl("`", c(response, terms), fixed = TRUE))) stop("Column names containing backticks are not supported", call. = FALSE)
  quoted_terms <- paste0("`", terms, "`")
  stats::as.formula(paste0("`", response, "` ~ ", paste(quoted_terms, collapse = " + ")))
}

screen_feature <- function(feature, training_data, exposure, covariates) {
  fit <- tryCatch(
    suppressWarnings(stats::glm(make_feature_formula(exposure, feature, covariates), data = training_data, family = stats::binomial())),
    error = function(e) NULL
  )
  if (is.null(fit)) return(data.table(Feature = feature, Estimate = NA_real_, P_value = NA_real_))
  sm <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  feature_row <- if (is.null(sm)) integer(0) else which(rownames(sm) %in% c(feature, paste0("`", feature, "`")))
  if (length(feature_row) != 1L) return(data.table(Feature = feature, Estimate = NA_real_, P_value = NA_real_))
  data.table(Feature = feature, Estimate = unname(sm[feature_row, "Estimate"]), P_value = unname(sm[feature_row, "Pr(>|z|)"]))
}

draw_coefficient_plot <- function(coef_data, output_file, max_features = Inf) {
  plot_data <- copy(coef_data)
  if (is.finite(max_features) && nrow(plot_data) > max_features) {
    plot_data <- plot_data[order(-abs(Coefficient))][seq_len(max_features)]
  }
  if (!nrow(plot_data)) return(invisible(NULL))
  plot_data <- plot_data[order(Coefficient)]
  plot_data[, Feature := factor(Feature, levels = Feature)]
  max_abs <- max(abs(plot_data$Coefficient))
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
  p <- ggplot(plot_data, aes(x = Feature, y = Coefficient, fill = Coefficient)) +
    geom_col(width = 0.7) +
    scale_fill_gradientn(colors = c("#2C7BB6", "#6BAED6", "#9ECAE1", "#FFFFFF", "#FDD49E", "#E66101", "#D7191C"), limits = c(-max_abs, max_abs), name = "Coefficient") +
    labs(x = NULL, y = "Coefficient") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
  grDevices::pdf(output_file, width = max(8, nrow(plot_data) * 0.25), height = 6)
  print(p)
  grDevices::dev.off()
  invisible(p)
}

draw_category_donut <- function(features, annotation_map, feature_column, category_column, palette, title, output_file) {
  if (is.null(annotation_map)) return(invisible(NULL))
  map <- as.data.table(copy(annotation_map))
  assert_columns(map, c(feature_column, category_column), "annotation_map")
  map <- unique(map[, .(Feature = as.character(get(feature_column)), Category = as.character(get(category_column)))])
  category <- map$Category[match(features, map$Feature)]
  category[is.na(category) | category == ""] <- "Other"
  category_counts <- table(category)
  counts <- data.table(Category = names(category_counts), Count = as.integer(category_counts))
  counts <- counts[Count > 0]
  if (!nrow(counts)) return(invisible(NULL))
  counts[, Percentage := round(Count / sum(Count) * 100, 1)]
  counts[, Category := factor(Category, levels = as.character(Category[order(-Percentage)]))]
  counts[, LegendLabel := paste0(as.character(Category), " (", Percentage, "%)")]
  cats <- levels(counts$Category)
  cols <- palette[cats]
  cols[is.na(cols)] <- "#B0B0B0"
  legend_labels <- stats::setNames(counts$LegendLabel, as.character(counts$Category))
  p <- ggplot(counts, aes(x = 2, y = Count, fill = Category)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y", start = 0) +
    scale_fill_manual(values = cols, labels = legend_labels, drop = FALSE) +
    labs(title = title, fill = "Category") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"), legend.position = "right") +
    xlim(0.5, 2.5)
  grDevices::pdf(output_file, width = 6, height = 6)
  print(p)
  grDevices::dev.off()
  invisible(p)
}

run_signature_construction <- function(omics_data, omics_type = c("proteomics", "metabolomics"), exposure = EXPOSURE, covariates = COVARIATES, annotation_map = NULL, annotation_feature_column = NULL, annotation_category_column = NULL, output_dir = OUTPUT_DIR, seed = RANDOM_SEED) {
  omics_type <- match.arg(omics_type)
  dt <- as.data.table(copy(omics_data))
  if (anyDuplicated(names(dt))) stop("Input contains duplicated column names", call. = FALSE)
  assert_columns(dt, c("eid", exposure, covariates), "omics_data")
  dt[, eid := as.character(eid)]
  validate_ids(dt, "omics_data")
  dt[, (exposure) := as_numeric_strict(get(exposure), exposure)]
  if (anyNA(dt[[exposure]]) || any(!dt[[exposure]] %in% c(0, 1))) stop("Exposure must be complete and coded 0/1", call. = FALSE)
  for (v in intersect(covariates, NUMERIC_COVARIATES)) dt[, (v) := as_numeric_strict(get(v), v)]
  for (v in intersect(covariates, FACTOR_COVARIATES)) dt[[v]] <- as.factor(dt[[v]])
  features <- setdiff(names(dt), c("eid", exposure, covariates))
  if (length(features) < 2L) stop("At least two omics features are required", call. = FALSE)
  dt[, (features) := lapply(.SD, function(x) as_numeric_strict(x, "omics feature")), .SDcols = features]
  model_columns <- c(exposure, covariates, features)
  if (anyNA(dt[, ..model_columns])) stop("Exposure, covariates, and omics features must be complete in the supplied imputed dataset", call. = FALSE)
  if (any(!vapply(dt[, ..features], function(x) all(is.finite(x)), logical(1)))) stop("Omics features must contain only finite values", call. = FALSE)
  if (length(unique(dt[[exposure]])) != 2L) stop("Both exposure groups are required", call. = FALSE)
  if (nrow(dt) < 30L) stop("Insufficient sample size for a 70/30 split and 10-fold cross-validation", call. = FALSE)
  for (v in intersect(covariates, FACTOR_COVARIATES)) {
    dt[[v]] <- droplevels(dt[[v]])
    if (nlevels(dt[[v]]) < 2L) stop("Categorical covariate has fewer than two levels: ", v, call. = FALSE)
  }
  for (v in intersect(covariates, NUMERIC_COVARIATES)) {
    if (any(!is.finite(dt[[v]]))) stop("Numeric covariate contains non-finite values: ", v, call. = FALSE)
  }
  set.seed(seed)
  train_index <- sample.int(nrow(dt), size = floor(0.70 * nrow(dt)))
  train_data <- copy(dt[train_index])
  valid_data <- copy(dt[-train_index])
  for (v in intersect(covariates, FACTOR_COVARIATES)) {
    train_data[[v]] <- droplevels(train_data[[v]])
    if (nlevels(train_data[[v]]) < 2L) stop("Training set has a single-level categorical covariate: ", v, call. = FALSE)
  }
  class_counts <- table(train_data[[exposure]])
  if (length(class_counts) != 2L || min(class_counts) < N_FOLDS) stop("The training set must contain at least 10 participants in each exposure class", call. = FALSE)
  if (length(unique(valid_data[[exposure]])) != 2L) stop("The validation split must contain both exposure classes", call. = FALSE)
  screen_list <- lapply(features, screen_feature, training_data = as.data.frame(train_data), exposure = exposure, covariates = covariates)
  screen_table <- rbindlist(screen_list, fill = TRUE)
  screen_table[, Abs_Estimate := abs(Estimate)]
  valid_screen <- screen_table[is.finite(P_value) & is.finite(Estimate)]
  if (!nrow(valid_screen)) stop("No valid adjusted univariate logistic models were fitted", call. = FALSE)
  if (omics_type == "proteomics") {
    threshold <- PROTEIN_BONFERRONI_ALPHA / length(features)
    candidates <- valid_screen[P_value < threshold & Abs_Estimate > PROTEIN_MIN_ABS_COEF]
    if (nrow(candidates) > MAX_PROTEIN_CANDIDATES) candidates <- candidates[order(-Abs_Estimate)][seq_len(MAX_PROTEIN_CANDIDATES)]
    candidates[, FDR := NA_real_]
  } else {
    valid_screen[, FDR := stats::p.adjust(P_value, method = "BH")]
    candidates <- valid_screen[FDR < METABOLOMICS_FDR_THRESHOLD]
    screen_table <- merge(screen_table, valid_screen[, .(Feature, FDR)], by = "Feature", all.x = TRUE, sort = FALSE)
  }
  if (nrow(candidates) < 2L) stop("Fewer than two candidate features passed the modality-specific screening threshold", call. = FALSE)
  candidate_features <- candidates$Feature
  x_train <- as.matrix(train_data[, ..candidate_features])
  storage.mode(x_train) <- "double"
  y_train <- as.integer(train_data[[exposure]])
  cv_fit <- glmnet::cv.glmnet(x = x_train, y = y_train, family = "binomial", alpha = 1, nfolds = N_FOLDS, type.measure = "deviance")
  selected_coef <- as.matrix(stats::coef(cv_fit, s = "lambda.min"))
  selected_features <- setdiff(rownames(selected_coef)[selected_coef[, 1L] != 0], "(Intercept)")
  if (!length(selected_features)) stop("LASSO selected no non-zero features at lambda.min", call. = FALSE)
  coefficients <- as.numeric(selected_coef[selected_features, 1L])
  names(coefficients) <- selected_features
  signature_coef <- data.table(Feature = selected_features, Coefficient = coefficients)
  train_matrix <- as.matrix(train_data[, ..selected_features])
  valid_matrix <- as.matrix(valid_data[, ..selected_features])
  all_matrix <- as.matrix(dt[, ..selected_features])
  storage.mode(train_matrix) <- "double"
  storage.mode(valid_matrix) <- "double"
  storage.mode(all_matrix) <- "double"
  raw_train <- as.vector(train_matrix %*% coefficients)
  raw_valid <- as.vector(valid_matrix %*% coefficients)
  raw_all <- as.vector(all_matrix %*% coefficients)
  score_mean <- mean(raw_train)
  score_sd <- stats::sd(raw_train)
  if (!is.finite(score_mean) || !is.finite(score_sd) || score_sd <= 0) stop("Training signature score has invalid standard deviation", call. = FALSE)
  score_name <- if (omics_type == "proteomics") "proteomic_signature" else "metabolic_signature"
  train_out <- data.table(eid = train_data$eid)
  train_out[, (score_name) := (raw_train - score_mean) / score_sd]
  valid_out <- data.table(eid = valid_data$eid)
  valid_out[, (score_name) := (raw_valid - score_mean) / score_sd]
  all_out <- data.table(eid = dt$eid)
  all_out[, (score_name) := (raw_all - score_mean) / score_sd]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- if (omics_type == "proteomics") "proteomic" else "metabolic"
  fwrite(signature_coef, file.path(output_dir, paste0(prefix, "_signature_coef.csv")))
  fwrite(screen_table, file.path(output_dir, paste0(prefix, "_feature_screen.csv")))
  fwrite(train_out, file.path(output_dir, paste0(prefix, "_signature_scores_train.csv")))
  fwrite(valid_out, file.path(output_dir, paste0(prefix, "_signature_scores_valid.csv")))
  fwrite(all_out, file.path(output_dir, paste0(prefix, "_signature_scores_all.csv")))
  saveRDS(list(omics_type = omics_type, exposure = exposure, covariates = covariates, selected_features = selected_features, coefficients = coefficients, score_mean = score_mean, score_sd = score_sd, seed = seed, train_fraction = TRAIN_FRACTION, screening = if (omics_type == "proteomics") "Bonferroni and absolute coefficient threshold, maximum 200 candidates" else "BH-FDR below 0.05", lambda = cv_fit$lambda.min), file.path(output_dir, paste0(prefix, "_signature_model.rds")))
  draw_coefficient_plot(signature_coef, file.path(output_dir, paste0(prefix, "_signature_coefficients.pdf")), if (omics_type == "proteomics") 50 else Inf)
  if (!is.null(annotation_map)) {
    if (is.null(annotation_feature_column) || is.null(annotation_category_column)) stop("Annotation column names must be supplied when annotation_map is provided", call. = FALSE)
    if (omics_type == "proteomics") {
      category_palette <- c("Cardiometabolic" = "#E5908F", "Oncology" = "#6C9DCE", "Neurology" = "#7BA37B", "Inflammation" = "#F09F60", "Other" = "#B0B0B0")
    } else {
      category_palette <- c("Cholesterol" = "#F7A3A3", "Glycolysis-related" = "#8EBDE0", "Inflammation" = "#F09F60", "Amino acids" = "#A0D29D", "Fluid balance" = "#8DD4D0", "Ketone bodies" = "#DCB0D6", "Fatty acids" = "#F8D276", "Ratio" = "#C0AAC4", "Lipid-related" = "#EDCEAD", "Other" = "#C4C4C4")
    }
    draw_category_donut(candidates$Feature, annotation_map, annotation_feature_column, annotation_category_column, category_palette, "Screened feature categories", file.path(output_dir, paste0(prefix, "_screened_categories.pdf")))
    draw_category_donut(selected_features, annotation_map, annotation_feature_column, annotation_category_column, category_palette, "Signature feature categories", file.path(output_dir, paste0(prefix, "_signature_categories.pdf")))
  }
  list(signature_coefficients = signature_coef, feature_screen = screen_table, candidates = candidates, selected_features = selected_features, train_scores = train_out, validation_scores = valid_out, all_scores = all_out, model = list(lambda = cv_fit$lambda.min, score_mean = score_mean, score_sd = score_sd))
}

if (!exists("analysis_inputs", envir = .GlobalEnv, inherits = FALSE)) stop("Provide analysis_inputs in memory before running this script", call. = FALSE)
if (is.null(analysis_inputs$baseline) || is.null(analysis_inputs$omics[[OMICS_TYPE]])) stop("analysis_inputs must contain baseline and the selected omics table", call. = FALSE)

signature_baseline <- data.table::as.data.table(data.table::copy(analysis_inputs$baseline))
signature_omics <- data.table::as.data.table(data.table::copy(analysis_inputs$omics[[OMICS_TYPE]]))
assert_columns(signature_baseline, c("eid", EXPOSURE, COVARIATES), "analysis_inputs$baseline")
assert_columns(signature_omics, "eid", paste0("analysis_inputs$omics$", OMICS_TYPE))
signature_baseline[, eid := as.character(eid)]
signature_omics[, eid := as.character(eid)]
omics_input <- merge(
  signature_baseline[, c("eid", EXPOSURE, COVARIATES), with = FALSE],
  signature_omics,
  by = "eid", all = FALSE, sort = FALSE
)

omics_annotation_map <- NULL
if (!is.null(analysis_inputs$annotation_maps)) {
  omics_annotation_map <- analysis_inputs$annotation_maps[[OMICS_TYPE]]
}
if (!is.null(omics_annotation_map)) {
  if (OMICS_TYPE == "proteomics") {
    signature_result <- run_signature_construction(omics_input, OMICS_TYPE, annotation_map = omics_annotation_map, annotation_feature_column = "Protein_abbreviation", annotation_category_column = "Protein_panel")
  } else {
    signature_result <- run_signature_construction(omics_input, OMICS_TYPE, annotation_map = omics_annotation_map, annotation_feature_column = "Metabolite_abbreviation", annotation_category_column = "Metabolite_category")
  }
} else {
  signature_result <- run_signature_construction(omics_input, OMICS_TYPE)
}
