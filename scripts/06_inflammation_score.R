required_packages <- c("data.table", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Install required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

EXPOSURE <- "any_allergy"
COVARIATES <- c("age", "sex", "BMI", "smoke", "alcohol", "ethnicity", "TDI", "IPAQ", "hypertension_status", "DM_status")
INFLAMMATORY_MARKERS <- c("WBC", "PLT", "NLR", "CRP")
RANDOM_SEED <- 123L
TRAIN_FRACTION <- 0.70
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

validate_ids <- function(x) {
  if (anyNA(x$eid) || anyDuplicated(x$eid)) stop("Input must have unique non-missing participant keys", call. = FALSE)
  invisible(TRUE)
}

decile_to_score <- function(x) {
  x <- as_numeric_strict(x, "inflammatory marker")
  if (any(!is.na(x) & !is.finite(x))) stop("Inflammatory markers must not contain infinite values", call. = FALSE)
  if (length(unique(x[!is.na(x)])) < 2L) stop("Each inflammatory marker must have at least two distinct observed values", call. = FALSE)
  decile_group <- dplyr::ntile(x, 10L)
  score_map <- c(-4, -3, -2, -1, 0, 0, 1, 2, 3, 4)
  as.numeric(score_map[decile_group])
}

build_inflammation_score <- function(inflammation_data, exposure = EXPOSURE, covariates = COVARIATES, markers = INFLAMMATORY_MARKERS, output_dir = OUTPUT_DIR, seed = RANDOM_SEED) {
  dt <- as.data.table(copy(inflammation_data))
  if (anyDuplicated(names(dt))) stop("Input contains duplicated column names", call. = FALSE)
  assert_columns(dt, c("eid", exposure, covariates, markers), "inflammation_data")
  dt[, eid := as.character(eid)]
  validate_ids(dt)
  dt[, (exposure) := as_numeric_strict(get(exposure), exposure)]
  if (any(!is.na(dt[[exposure]]) & !dt[[exposure]] %in% c(0, 1))) stop("Exposure must be coded 0/1/NA", call. = FALSE)
  for (v in intersect(c("age", "BMI", "TDI"), covariates)) dt[, (v) := as_numeric_strict(get(v), v)]
  for (v in intersect(c("sex", "smoke", "alcohol", "ethnicity", "IPAQ", "hypertension_status", "DM_status"), covariates)) dt[[v]] <- as.factor(dt[[v]])
  for (marker in markers) dt[, (marker) := as_numeric_strict(get(marker), marker)]
  for (marker in markers) {
    if (any(!is.na(dt[[marker]]) & !is.finite(dt[[marker]]))) stop(marker, " contains infinite values", call. = FALSE)
  }
  if (nrow(dt) < 2L) stop("At least two participants are required", call. = FALSE)
  for (marker in markers) dt[, (paste0("Score_", marker)) := decile_to_score(get(marker))]
  component_columns <- paste0("Score_", markers)
  dt[, INFLA_score := rowSums(.SD), .SDcols = component_columns]
  set.seed(seed)
  train_index <- sample.int(nrow(dt), size = floor(TRAIN_FRACTION * nrow(dt)))
  if (!length(train_index) || length(train_index) >= nrow(dt)) stop("The 70/30 split produced an empty training or validation set", call. = FALSE)
  training_scores <- dt$INFLA_score[train_index]
  score_mean <- mean(training_scores, na.rm = TRUE)
  score_sd <- stats::sd(training_scores, na.rm = TRUE)
  if (!is.finite(score_mean) || !is.finite(score_sd) || score_sd <= 0) stop("Training INFLA_score has invalid mean or standard deviation", call. = FALSE)
  dt[, INFLA_score_z := (INFLA_score - score_mean) / score_sd]
  keep_columns <- c("eid", exposure, covariates, "INFLA_score", "INFLA_score_z")
  final_scores <- dt[, ..keep_columns]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(final_scores, file.path(output_dir, "inflammation_final_clean_scores_for_survival.csv"))
  saveRDS(list(markers = markers, score_map = c(-4, -3, -2, -1, 0, 0, 1, 2, 3, 4), score_mean_train = score_mean, score_sd_train = score_sd, seed = seed, train_fraction = TRAIN_FRACTION), file.path(output_dir, "inflammation_score_parameters.rds"))
  list(scores = final_scores, score_mean_train = score_mean, score_sd_train = score_sd, train_index = train_index)
}

if (!exists("analysis_inputs", envir = .GlobalEnv, inherits = FALSE)) stop("Provide analysis_inputs in memory before running this script", call. = FALSE)
if (is.null(analysis_inputs$baseline) || is.null(analysis_inputs$omics$inflammation)) stop("analysis_inputs must contain baseline and omics$inflammation", call. = FALSE)

inflammation_baseline <- data.table::as.data.table(data.table::copy(analysis_inputs$baseline))
inflammation_omics <- data.table::as.data.table(data.table::copy(analysis_inputs$omics$inflammation))
assert_columns(inflammation_baseline, c("eid", EXPOSURE, COVARIATES), "analysis_inputs$baseline")
assert_columns(inflammation_omics, c("eid", INFLAMMATORY_MARKERS), "analysis_inputs$omics$inflammation")
inflammation_baseline[, eid := as.character(eid)]
inflammation_omics[, eid := as.character(eid)]
inflammation_input <- merge(
  inflammation_baseline[, c("eid", EXPOSURE, COVARIATES), with = FALSE],
  inflammation_omics,
  by = "eid", all = FALSE, sort = FALSE
)
inflammation_result <- build_inflammation_score(inflammation_input)
