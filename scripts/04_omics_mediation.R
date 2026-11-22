












MEDIATION_COVARIATES <- c(
  "age", "sex", "BMI", "smoke", "alcohol", "ethnicity", "TDI", "IPAQ",
  "hypertension_status", "DM_status"
)

CATEGORICAL_COVARIATES <- c(
  "sex", "smoke", "alcohol", "ethnicity", "IPAQ",
  "hypertension_status", "DM_status"
)

NUMERIC_COVARIATES <- c("age", "BMI", "TDI")


ANALYSIS_DISEASES <- c(
  "Diabetes", "Hypertension", "Arrhythmia", "Pulmonale", "COPD", "CKD",
  "Liver disease", "FCID", "Anxiety", "Anaemia", "Autoimmune", "Death"
)

BOOTSTRAP_REPS <- 1000L
BOOTSTRAP_CI_TYPE <- "perc"
BOOTSTRAP_CONFIDENCE_LEVEL <- 0.95

MIN_ANALYSIS_N <- 200L
MIN_EVENT_N <- 20L
DEFAULT_CUTOFF_DATE <- as.Date("2023-03-31")


assert_columns <- function(x, required, object_name) {
  missing_columns <- setdiff(required, names(x))
  if (length(missing_columns) > 0L) {
    stop(
      object_name, " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

parse_date_strict <- function(x, column_name) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))

  x_chr <- as.character(x)
  x_chr[is.na(x_chr) | trimws(x_chr) == ""] <- NA_character_

  parsed <- as.Date(x_chr, format = "%Y-%m-%d")
  needs_slash_format <- !is.na(x_chr) & is.na(parsed)
  if (any(needs_slash_format)) {
    parsed[needs_slash_format] <- as.Date(
      x_chr[needs_slash_format], format = "%Y/%m/%d"
    )
  }

  invalid <- !is.na(x_chr) & is.na(parsed)
  if (any(invalid)) {
    stop("Unparseable date values in ", column_name, call. = FALSE)
  }
  parsed
}

as_numeric_strict <- function(x, column_name) {
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
    stop(column_name, " must be numeric", call. = FALSE)
  }
  as.numeric(x)
}

as_binary_strict <- function(x, column_name) {
  x_num <- as_numeric_strict(x, column_name)
  invalid <- !is.na(x_num) & !(x_num %in% c(0, 1))
  if (any(invalid)) {
    stop(column_name, " must be coded 0/1/NA", call. = FALSE)
  }
  as.integer(x_num)
}

validate_unique_ids <- function(x, object_name) {
  if (anyNA(x$eid) || anyDuplicated(x$eid)) {
    stop(object_name, " must contain unique, non-missing IDs", call. = FALSE)
  }
  invisible(TRUE)
}


prepare_cohort <- function(cohort, exposure) {
  dt <- data.table::as.data.table(data.table::copy(cohort))

  if (anyDuplicated(names(dt))) {
    stop("cohort contains duplicated column names", call. = FALSE)
  }

  assert_columns(
    dt,
    c("eid", exposure, MEDIATION_COVARIATES, "base_date", "death_date"),
    "cohort"
  )

  dt[, eid := as.character(eid)]
  validate_unique_ids(dt, "cohort")

  dt[, (exposure) := as_binary_strict(get(exposure), exposure)]
  dt[, base_date := parse_date_strict(base_date, "base_date")]
  dt[, death_date := parse_date_strict(death_date, "death_date")]

  if (anyNA(dt$base_date)) {
    stop("cohort contains missing baseline dates", call. = FALSE)
  }

  for (v in NUMERIC_COVARIATES) {
    dt[, (v) := as_numeric_strict(get(v), v)]
  }
  for (v in CATEGORICAL_COVARIATES) {
    dt[[v]] <- as.factor(dt[[v]])
  }

  keep_columns <- unique(c("eid", exposure, MEDIATION_COVARIATES,
                           "base_date", "death_date"))
  dt[, ..keep_columns]
}

prepare_omics <- function(omics, mediator_names, omics_type, exposure) {
  dt <- data.table::as.data.table(data.table::copy(omics))

  if (anyDuplicated(names(dt))) {
    stop("omics contains duplicated column names", call. = FALSE)
  }
  assert_columns(dt, "eid", "omics")



  if (identical(omics_type, "inflammation")) {
    duplicated_suffix_columns <- grep("\\.\\.1$", names(dt), value = TRUE)
    if (length(duplicated_suffix_columns) > 0L) {
      dt[, (duplicated_suffix_columns) := NULL]
    }
  }

  dt[, eid := as.character(eid)]
  validate_unique_ids(dt, "omics")

  if (is.null(mediator_names)) {
    mediator_names <- setdiff(
      names(dt),
      c("eid", exposure, "Exposure", "Mediator", MEDIATION_COVARIATES,
        "base_date", "death_date", "lost_date")
    )
  }

  if (!length(mediator_names) || anyDuplicated(mediator_names)) {
    stop("Provide a non-empty vector of unique mediator column names", call. = FALSE)
  }
  assert_columns(dt, mediator_names, "omics")

  forbidden_names <- c(
    "eid", "Exposure", "Mediator", MEDIATION_COVARIATES,
    "base_date", "death_date", "surv_time", "event"
  )
  if (any(mediator_names %in% forbidden_names)) {
    stop("Mediator names overlap with model variables", call. = FALSE)
  }

  for (v in mediator_names) {
    dt[, (v) := as_numeric_strict(get(v), v)]
  }

  dt[, c("eid", mediator_names), with = FALSE]
}








prepare_survival_data <- function(data, disease, outcomes,
                                 cutoff_date, outcome_time_unit) {
  dt <- data.table::copy(data)
  dt[, followup_limit_days := as.numeric(cutoff_date - base_date)]

  if (disease == "Death") {
    death_days <- as.numeric(dt$death_date - dt$base_date)
    if (any(!is.na(death_days) & death_days < 0)) {
      stop("death_date precedes baseline date", call. = FALSE)
    }
    death_before_cutoff <- !is.na(death_days) &
      death_days >= 0 & death_days <= dt$followup_limit_days

    dt[, event := as.integer(death_before_cutoff)]
    dt[, surv_time_days := ifelse(
      death_before_cutoff, death_days, followup_limit_days
    )]
  } else {
    if (!disease %in% names(outcomes)) {
      stop("No outcome table supplied for: ", disease, call. = FALSE)
    }

    outcome_dt <- data.table::as.data.table(data.table::copy(outcomes[[disease]]))
    assert_columns(outcome_dt, c("eid", "time", "status"), paste0("outcomes[['", disease, "']]") )
    outcome_dt[, eid := as.character(eid)]
    validate_unique_ids(outcome_dt, paste0("outcome table: ", disease))

    outcome_dt[, time := as_numeric_strict(time, paste0(disease, " time"))]
    outcome_dt[, status := as_binary_strict(status, paste0(disease, " status"))]

    outcome_dt[, time_days := if (outcome_time_unit == "years") {
      time * 365.25
    } else {
      time
    }]

    outcome_dt <- outcome_dt[, .(eid, time_days, status)]
    dt <- merge(dt, outcome_dt, by = "eid", all.x = TRUE, sort = FALSE)
    dt[, followup_limit_days := as.numeric(cutoff_date - base_date)]

    valid_outcome <- !is.na(dt$status) & !is.na(dt$time_days) &
      is.finite(dt$time_days) & dt$time_days >= 0

    dt <- dt[valid_outcome]
    dt[, event := as.integer(status == 1L & time_days <= followup_limit_days)]
    dt[, surv_time_days := pmin(time_days, followup_limit_days)]
  }

  dt <- dt[is.finite(surv_time_days) & surv_time_days > 0 &
             is.finite(followup_limit_days) & followup_limit_days > 0]

  dt[, surv_time := surv_time_days / 365.25]
  dt[, event := as.integer(event)]
  columns_to_drop <- intersect(
    c("surv_time_days", "time_days", "status", "followup_limit_days"),
    names(dt)
  )
  if (length(columns_to_drop) > 0L) {
    dt[, (columns_to_drop) := NULL]
  }
  dt
}


fit_one_mediator <- function(data, disease, feature, exposure,
                             min_n, min_events, random_seed) {
  required <- c(exposure, MEDIATION_COVARIATES, feature, "surv_time", "event")
  x <- data.table::copy(data[, ..required])
  data.table::setnames(x, old = c(exposure, feature), new = c("Exposure", "Mediator"))

  x <- x[stats::complete.cases(x)]
  x <- x[is.finite(surv_time) & surv_time > 0 &
           is.finite(Mediator) & event %in% c(0L, 1L)]

  x[, Exposure := as_binary_strict(Exposure, "Exposure")]
  if (nrow(x) < min_n) stop("fewer than the required complete-case participants")
  if (!identical(sort(unique(x$Exposure)), c(0L, 1L))) {
    stop("both exposure groups (0 and 1) are required")
  }
  if (sum(x$event == 1L) < min_events) {
    stop("fewer than the required outcome events")
  }
  if (length(unique(x$Mediator)) < 2L ||
      !is.finite(stats::sd(x$Mediator)) || stats::sd(x$Mediator) <= 0) {
    stop("mediator has zero or invalid variance")
  }

  for (v in CATEGORICAL_COVARIATES) {
    x[[v]] <- droplevels(as.factor(x[[v]]))
    if (nlevels(x[[v]]) < 2L) {
      stop("categorical covariate has fewer than two levels: ", v)
    }
  }
  for (v in NUMERIC_COVARIATES) {
    if (!is.numeric(x[[v]]) || any(!is.finite(x[[v]]))) {
      stop("numeric covariate contains invalid values: ", v)
    }
  }


  x[, Mediator := as.numeric(scale(Mediator))]

  covariate_terms <- paste(MEDIATION_COVARIATES, collapse = " + ")
  mediator_formula <- stats::as.formula(
    paste("Mediator ~ Exposure +", covariate_terms),
    env = environment()
  )



  Surv <- survival::Surv
  outcome_formula <- stats::as.formula(
    paste0("Surv(surv_time, event) ~ Exposure + Mediator + ", covariate_terms),
    env = environment()
  )

  mediator_fit <- stats::lm(
    mediator_formula, data = x, na.action = stats::na.fail, model = TRUE
  )
  outcome_fit <- survival::survreg(
    outcome_formula, data = x, dist = "weibull",
    na.action = stats::na.fail, model = TRUE, x = TRUE, y = TRUE
  )

  if (anyNA(stats::coef(mediator_fit)) || any(!is.finite(stats::coef(mediator_fit)))) {
    stop("mediator model contains aliased or non-finite coefficients")
  }
  if (anyNA(stats::coef(outcome_fit)) || any(!is.finite(stats::coef(outcome_fit)))) {
    stop("Weibull outcome model contains aliased or non-finite coefficients")
  }

  set.seed(random_seed)
  mediation_fit <- mediation::mediate(
    model.m = mediator_fit,
    model.y = outcome_fit,
    treat = "Exposure",
    mediator = "Mediator",
    outcome = "surv_time",
    control.value = 0,
    treat.value = 1,
    boot = TRUE,
    sims = BOOTSTRAP_REPS,
    boot.ci.type = BOOTSTRAP_CI_TYPE,
    conf.level = BOOTSTRAP_CONFIDENCE_LEVEL,
    long = TRUE,
    dropobs = FALSE
  )


  if (!isTRUE(mediation_fit$boot) ||
      !identical(as.integer(mediation_fit$sims), BOOTSTRAP_REPS)) {
    stop("mediation package did not report the required 1,000 bootstrap draws")
  }

  required_draws <- c("d0.sims", "z0.sims", "tau.sims", "n0.sims")
  if (!all(required_draws %in% names(mediation_fit))) {
    stop("bootstrap simulation draws are missing from the mediation result")
  }

  draw_counts <- vapply(
    mediation_fit[required_draws],
    function(z) sum(is.finite(as.numeric(z))),
    integer(1)
  )
  draw_lengths <- vapply(mediation_fit[required_draws], length, integer(1))

  if (any(draw_lengths[1:3] != BOOTSTRAP_REPS) ||
      any(draw_counts[1:3] != BOOTSTRAP_REPS)) {
    stop(
      "ACME/ADE/total-effect bootstrap did not yield 1,000 finite draws; counts: ",
      paste(names(draw_counts)[1:3], draw_counts[1:3], sep = "=", collapse = "; ")
    )
  }

  mediator_path <- summary(mediator_fit)$coefficients["Exposure", ]
  mediator_coef <- unname(stats::coef(outcome_fit)["Mediator"])
  mediator_se <- sqrt(stats::vcov(outcome_fit)["Mediator", "Mediator"])

  if (!is.finite(mediator_coef) || !is.finite(mediator_se) || mediator_se <= 0) {
    stop("invalid mediator coefficient in the Weibull outcome model")
  }

  prop_valid <- unname(draw_counts["n0.sims"])
  prop_estimate <- if (is.finite(mediation_fit$n0)) {
    unname(mediation_fit$n0)
  } else {
    NA_real_
  }
  prop_p <- if (prop_valid == BOOTSTRAP_REPS) {
    unname(mediation_fit$n0.p)
  } else {
    NA_real_
  }

  data.table::data.table(
    Disease = disease,
    Feature = feature,
    Sample_N = nrow(x),
    Event_N = sum(x$event == 1L),
    ACME_Est = unname(mediation_fit$d0),
    ACME_2.5_ci = unname(mediation_fit$d0.ci[1]),
    ACME_97.5_ci = unname(mediation_fit$d0.ci[2]),
    ACME_P = unname(mediation_fit$d0.p),
    ADE_Est = unname(mediation_fit$z0),
    ADE_2.5_ci = unname(mediation_fit$z0.ci[1]),
    ADE_97.5_ci = unname(mediation_fit$z0.ci[2]),
    ADE_P = unname(mediation_fit$z0.p),
    Total_Est = unname(mediation_fit$tau.coef),
    Total_P = unname(mediation_fit$tau.p),
    Prop_Mediated = prop_estimate,
    Prop_Mediated_P = prop_p,
    Path_a_Est = unname(mediator_path["Estimate"]),
    Path_a_P = unname(mediator_path["Pr(>|t|)"]),
    Path_b_Est = mediator_coef,
    Path_b_P = 2 * stats::pnorm(-abs(mediator_coef / mediator_se)),
    Bootstrap_Method = "Nonparametric percentile bootstrap",
    Bootstrap_Requested = BOOTSTRAP_REPS,
    Bootstrap_Valid_ACME = unname(draw_counts["d0.sims"]),
    Bootstrap_Valid_ADE = unname(draw_counts["z0.sims"]),
    Bootstrap_Valid_Total = unname(draw_counts["tau.sims"]),
    Bootstrap_Valid_Prop = prop_valid
  )
}










run_omics_mediation <- function(
    cohort,
    omics,
    outcomes,
    omics_type = c("proteomics", "inflammation", "metabolomics"),
    exposure = "any_allergy",
    diseases = ANALYSIS_DISEASES,
    mediator_names = NULL,
    cutoff_date = DEFAULT_CUTOFF_DATE,
    outcome_time_unit = "days",
    min_n = MIN_ANALYSIS_N,
    min_events = MIN_EVENT_N,
    random_seed = 2026L,
    previous_results = NULL,
    verbose = TRUE) {

  omics_type <- match.arg(omics_type)

  if (is.null(outcome_time_unit) || length(outcome_time_unit) != 1L ||
      !outcome_time_unit %in% c("days", "years")) {
    stop("Set outcome_time_unit explicitly to 'days' or 'years'.", call. = FALSE)
  }
  if (!is.list(outcomes) || is.null(names(outcomes))) {
    stop("outcomes must be a named list keyed by the analysis disease labels.",
         call. = FALSE)
  }
  if (anyDuplicated(diseases) || !length(diseases)) {
    stop("diseases must be a non-empty vector of unique labels.", call. = FALSE)
  }
  if (!all(diseases %in% ANALYSIS_DISEASES)) {
    stop("Unknown disease labels: ",
         paste(setdiff(diseases, ANALYSIS_DISEASES), collapse = ", "),
         call. = FALSE)
  }
  if (!is.numeric(min_n) || length(min_n) != 1L || min_n < 2L ||
      !is.numeric(min_events) || length(min_events) != 1L || min_events < 1L) {
    stop("Invalid minimum sample/event thresholds.", call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE) ||
      !requireNamespace("survival", quietly = TRUE) ||
      !requireNamespace("mediation", quietly = TRUE) ||
      !requireNamespace("boot", quietly = TRUE)) {
    stop("Required packages: data.table, survival, mediation, boot.", call. = FALSE)
  }

  cutoff_date <- parse_date_strict(cutoff_date, "cutoff_date")
  if (length(cutoff_date) != 1L || is.na(cutoff_date)) {
    stop("cutoff_date must be one valid date.", call. = FALSE)
  }

  cohort_dt <- prepare_cohort(cohort, exposure)
  omics_dt <- prepare_omics(omics, mediator_names, omics_type, exposure)
  features <- setdiff(names(omics_dt), "eid")

  if (any(features %in% c(exposure, MEDIATION_COVARIATES,
                          "base_date", "death_date", "surv_time", "event"))) {
    stop("Mediator names conflict with exposure, covariates, or outcome variables.",
         call. = FALSE)
  }

  analysis_data <- merge(cohort_dt, omics_dt, by = "eid", all = FALSE, sort = FALSE)
  if (!nrow(analysis_data)) {
    stop("No overlapping participant IDs between cohort and omics inputs.",
         call. = FALSE)
  }

  results <- if (is.null(previous_results)) {
    data.table::data.table()
  } else {
    data.table::as.data.table(data.table::copy(previous_results))
  }

  if (nrow(results) > 0L) {
    assert_columns(results, c("Disease", "Feature", "Omics_Type"), "previous_results")
    if (anyNA(results$Omics_Type) || any(results$Omics_Type != omics_type)) {
      stop("previous_results contains a different or missing omics modality.", call. = FALSE)
    }
    if (anyNA(results$Disease) || anyNA(results$Feature)) {
      stop("previous_results contains missing Disease/Feature keys.", call. = FALSE)
    }
    if (anyDuplicated(results[, .(Disease, Feature)])) {
      stop("previous_results contains duplicated Disease/Feature keys.", call. = FALSE)
    }
  }

  failures <- list()
  failure_index <- 0L

  for (disease in diseases) {
    if (verbose) message("Outcome: ", disease, " | modality: ", omics_type)

    disease_data <- tryCatch(
      prepare_survival_data(
        analysis_data, disease, outcomes, cutoff_date, outcome_time_unit
      ),
      error = function(e) e
    )

    if (inherits(disease_data, "error")) {
      failure_index <- failure_index + 1L
      failures[[failure_index]] <- data.table::data.table(
        Disease = disease, Feature = NA_character_,
        Reason = conditionMessage(disease_data)
      )
      next
    }

    if (sum(disease_data$event == 1L) < min_events) {
      failure_index <- failure_index + 1L
      failures[[failure_index]] <- data.table::data.table(
        Disease = disease, Feature = NA_character_,
        Reason = "Insufficient outcome events after outcome preparation"
      )
      next
    }

    for (feature in features) {
      already_completed <- nrow(results) > 0L &&
        any(results$Disease == disease & results$Feature == feature)
      if (already_completed) next

      seed_key <- utf8ToInt(paste(omics_type, disease, feature, sep = "::"))
      feature_seed <- as.integer(
        (as.double(random_seed) +
           sum(as.double(seed_key) * seq_along(seed_key))) %% 2147483646 + 1
      )

      fit_result <- tryCatch(
        fit_one_mediator(
          disease_data, disease, feature, exposure,
          min_n, min_events, feature_seed
        ),
        error = function(e) e
      )

      if (inherits(fit_result, "error")) {
        failure_index <- failure_index + 1L
        failures[[failure_index]] <- data.table::data.table(
          Disease = disease, Feature = feature,
          Reason = conditionMessage(fit_result)
        )
      } else {
        fit_result[, Omics_Type := omics_type]
        results <- data.table::rbindlist(list(results, fit_result), fill = TRUE)
      }
    }
  }

  if (nrow(results) > 0L) {
    results[, ACME_P_FDR := stats::p.adjust(ACME_P, method = "BH"),
            by = .(Omics_Type, Disease)]
    data.table::setcolorder(
      results,
      c("Omics_Type", "Disease", "Feature",
        setdiff(names(results), c("Omics_Type", "Disease", "Feature")))
    )
    data.table::setorder(results, Disease, Feature)
  }

  failure_table <- if (length(failures) > 0L) {
    data.table::rbindlist(failures, fill = TRUE)
  } else {
    data.table::data.table(
      Disease = character(), Feature = character(), Reason = character()
    )
  }

  list(
    results = results,
    failures = failure_table,
    settings = list(
      omics_type = omics_type,
      exposure = exposure,
      diseases = diseases,
      cutoff_date = cutoff_date,
      outcome_time_unit = outcome_time_unit,
      bootstrap_method = "Nonparametric percentile bootstrap",
      bootstrap_replicates = BOOTSTRAP_REPS,
      bootstrap_confidence_level = BOOTSTRAP_CONFIDENCE_LEVEL,
      min_n = min_n,
      min_events = min_events,
      random_seed = random_seed
    )
  )
}

if (!exists("analysis_inputs", envir = .GlobalEnv, inherits = FALSE)) {
  stop("Provide analysis_inputs in memory before running this script.", call. = FALSE)
}
if (!is.list(analysis_inputs) || is.null(analysis_inputs$baseline) ||
    is.null(analysis_inputs$disease_by_system) || is.null(analysis_inputs$omics)) {
  stop("analysis_inputs must contain baseline, disease_by_system, and omics.", call. = FALSE)
}

mediation_baseline <- data.table::as.data.table(data.table::copy(analysis_inputs$baseline))
if (!"base_date" %in% names(mediation_baseline)) {
  if (!"baseline_date" %in% names(mediation_baseline)) {
    stop("baseline must contain baseline_date for mediation analysis.", call. = FALSE)
  }
  mediation_baseline[, base_date := baseline_date]
}
if (!"death_date" %in% names(mediation_baseline)) {
  mediation_baseline[, death_date := as.Date(NA)]
}

mediation_outcomes <- list()
for (system_name in names(analysis_inputs$disease_by_system)) {
  system_outcomes <- analysis_inputs$disease_by_system[[system_name]]
  for (disease_name in names(system_outcomes)) {
    if (disease_name %in% names(mediation_outcomes)) {
      stop("Disease labels must be unique across systems for mediation analysis: ", disease_name, call. = FALSE)
    }
    mediation_outcomes[[disease_name]] <- data.table::as.data.table(
      data.table::copy(system_outcomes[[disease_name]])
    )
  }
}

mediation_config <- analysis_inputs$config
if (is.null(mediation_config)) mediation_config <- list()
mediation_omics_type <- mediation_config$omics_type
if (is.null(mediation_omics_type)) mediation_omics_type <- "proteomics"
if (!mediation_omics_type %in% c("proteomics", "inflammation", "metabolomics")) {
  stop("config$omics_type must be proteomics, inflammation, or metabolomics.", call. = FALSE)
}
if (is.null(analysis_inputs$omics[[mediation_omics_type]])) {
  stop("Missing omics input for config$omics_type: ", mediation_omics_type, call. = FALSE)
}
mediation_diseases <- mediation_config$diseases
if (is.null(mediation_diseases)) {
  mediation_diseases <- intersect(ANALYSIS_DISEASES, c(names(mediation_outcomes), "Death"))
}
mediation_outcome_time_unit <- mediation_config$outcome_time_unit
if (is.null(mediation_outcome_time_unit)) mediation_outcome_time_unit <- "days"
mediation_cutoff_date <- mediation_config$cutoff_date
if (is.null(mediation_cutoff_date)) mediation_cutoff_date <- DEFAULT_CUTOFF_DATE

mediation_result <- run_omics_mediation(
  cohort = mediation_baseline,
  omics = analysis_inputs$omics[[mediation_omics_type]],
  outcomes = mediation_outcomes,
  omics_type = mediation_omics_type,
  diseases = mediation_diseases,
  cutoff_date = mediation_cutoff_date,
  outcome_time_unit = mediation_outcome_time_unit
)
