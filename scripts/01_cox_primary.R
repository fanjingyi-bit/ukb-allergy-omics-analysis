

























required_packages <- c("data.table", "survival", "matrixStats", "writexl")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(matrixStats)
  library(writexl)
})


TARGET_ALLERGY <- "any_allergy"
ADMIN_CUTOFF_DATE <- as.Date("2023-03-31")
PREVALENCE_THRESHOLD <- 0.005
MIN_STRATUM_SIZE <- 100L
OUTPUT_DIR <- Sys.getenv("ANALYSIS_OUTPUT_DIR", unset = "aggregate_outputs")

if (!exists("analysis_inputs", envir = .GlobalEnv, inherits = FALSE)) {
  stop("The secure analysis environment must provide `analysis_inputs` in memory.")
}
if (!is.list(analysis_inputs) || is.null(analysis_inputs$baseline) ||
    is.null(analysis_inputs$disease_by_system)) {
  stop("`analysis_inputs` must contain `baseline` and `disease_by_system`.")
}

as_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  as.Date(x, tryFormats = c("%Y-%m-%d", "%Y/%m/%d", "%d/%m/%Y"))
}

assert_columns <- function(x, required, object_name) {
  absent <- setdiff(required, names(x))
  if (length(absent) > 0L) {
    stop(object_name, " is missing required columns: ",
         paste(absent, collapse = ", "))
  }
}

safe_key <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  out <- gsub("^_+|_+$", "", out)
  ifelse(nchar(out) == 0L, "unnamed", out)
}

normalise_baseline <- function(x) {
  dt <- as.data.table(copy(x))
  required <- c("eid", TARGET_ALLERGY, "age", "sex", "ethnicity",
                "smoke", "alcohol", "IPAQ", "BMI", "TDI",
                "hypertension_status", "DM_status", "baseline_date")
  assert_columns(dt, required, "analysis_inputs$baseline")
  if (anyDuplicated(dt$eid)) {
    stop("Baseline input must contain exactly one row per participant.")
  }
  if (!"death_date" %in% names(dt)) dt[, death_date := as.Date(NA)]
  if (!"lost_date" %in% names(dt)) dt[, lost_date := as.Date(NA)]
  dt[, baseline_date := as_date_safe(baseline_date)]
  dt[, death_date := as_date_safe(death_date)]
  dt[, lost_date := as_date_safe(lost_date)]
  dt[, sex := suppressWarnings(as.numeric(as.character(sex)))]
  numeric_covariates <- c("age", "BMI", "TDI")
  dt[, (numeric_covariates) := lapply(.SD, function(v) suppressWarnings(as.numeric(as.character(v)))),
     .SDcols = numeric_covariates]
  dt[, (TARGET_ALLERGY) := suppressWarnings(as.numeric(as.character(get(TARGET_ALLERGY))))]
  if (any(!is.na(dt$sex) & !dt$sex %in% c(0, 1))) {
    stop("`sex` must be coded 0 (female) or 1 (male).")
  }
  if (any(!is.na(dt[[TARGET_ALLERGY]]) & !dt[[TARGET_ALLERGY]] %in% c(0, 1))) {
    stop("The exposure must be coded 0/1.")
  }
  dt[, censor_date := pmin(death_date, lost_date, ADMIN_CUTOFF_DATE, na.rm = TRUE)]
  dt[, max_follow_days := as.numeric(censor_date - baseline_date)]
  dt[is.na(baseline_date) | is.na(max_follow_days) | max_follow_days <= 0,
     max_follow_days := NA_real_]
  dt
}

validate_disease_list <- function(x) {
  if (!is.list(x) || length(x) == 0L || is.null(names(x)) ||
      any(names(x) == "")) {
    stop("`disease_by_system` must be a non-empty named list of named disease lists.")
  }
  for (sys in names(x)) {
    if (!is.list(x[[sys]]) || length(x[[sys]]) == 0L ||
        is.null(names(x[[sys]]) ) || any(names(x[[sys]]) == "")) {
      stop("Each system must contain a non-empty named list of disease outcomes: ", sys)
    }
    disease_keys <- vapply(names(x[[sys]]), safe_key, character(1))
    if (anyDuplicated(disease_keys)) {
      stop("Disease names collide after safe-name conversion in system: ", sys)
    }
    for (dis in names(x[[sys]])) {
      d <- as.data.table(x[[sys]][[dis]])
      assert_columns(d, c("eid", "time", "status"),
                     paste0("disease outcome ", sys, "/", dis))
      if (anyDuplicated(d$eid)) {
        stop("Duplicate participant keys in disease outcome: ", sys, "/", dis)
      }
    }
  }
  invisible(TRUE)
}

prepare_outcome <- function(outcome, baseline) {
  d <- as.data.table(copy(outcome))
  d[, time_days := suppressWarnings(as.numeric(time))]
  d[, event := suppressWarnings(as.integer(as.character(status)))]
  d <- d[, .(eid, time_days, event)]
  d <- merge(baseline[, .(eid, max_follow_days)], d,
             by = "eid", all.x = TRUE, sort = FALSE)
  invalid <- is.na(d$time_days) | is.na(d$event) |
    !d$event %in% c(0L, 1L) | is.na(d$max_follow_days)
  d[invalid, `:=`(time_days = NA_real_, event = NA_integer_)]
  beyond_censor <- !is.na(d$time_days) & d$time_days > d$max_follow_days
  d[beyond_censor, `:=`(time_days = max_follow_days, event = 0L)]
  d[time_days <= 0, `:=`(time_days = NA_real_, event = NA_integer_)]
  d[, .(eid, time_days, event)]
}

fit_cox_model <- function(data, time_col, event_col, exposure, covariates,
                          model_name, max_iter = 100L) {
  needed <- unique(c(time_col, event_col, exposure, covariates))
  d <- copy(as.data.table(data))
  d <- d[complete.cases(d[, ..needed])]
  d <- d[is.finite(get(time_col)) & get(time_col) > 0 &
           get(event_col) %in% c(0, 1)]
  if (nrow(d) < 20L || length(unique(d[[exposure]])) < 2L ||
      sum(d[[event_col]] == 1L) == 0L) return(NULL)
  rhs <- paste(c(exposure, covariates), collapse = " + ")
  f <- as.formula(paste0("survival::Surv(", time_col, ", ", event_col,") ~ ", rhs))
  fit <- tryCatch(
    survival::coxph(f, data = d, ties = "efron", x = TRUE, model = TRUE,
                    control = survival::coxph.control(iter.max = max_iter)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  b <- stats::coef(fit)
  ci <- tryCatch(stats::confint(fit), error = function(e) NULL)
  sm <- summary(fit)
  if (is.null(ci) || !exposure %in% names(b) ||
      !exposure %in% rownames(sm$coefficients) ||
      any(!is.finite(b)) || any(!is.finite(ci))) return(NULL)
  list(
    fit = fit,
    data = d,
    model = model_name,
    HR = unname(exp(b[[exposure]])),
    CI_lower = unname(exp(ci[exposure, 1L])),
    CI_upper = unname(exp(ci[exposure, 2L])),
    P_raw = unname(sm$coefficients[exposure, "Pr(>|z|)"]),
    N = nrow(d),
    Events = sum(d[[event_col]] == 1L)
  )
}

ph_diagnostics <- function(model_object, time_col, event_col, exposure,
                           covariates, time_scale_per_year = TRUE) {
  empty <- list(Global_PH_P = NA_real_, PH_satisfied = NA,
                time_interaction_HR_per_year = NA_real_,
                time_interaction_P = NA_real_, split_time_years = NA_real_,
                period1_HR = NA_real_, period2_HR = NA_real_)
  fit <- model_object$fit
  d <- model_object$data
  zph <- tryCatch(survival::cox.zph(fit, transform = "rank"),
                  error = function(e) NULL)
  if (is.null(zph) || !"GLOBAL" %in% rownames(zph$table)) return(empty)
  p_global <- unname(zph$table["GLOBAL", "p"])
  out <- empty
  out$Global_PH_P <- p_global
  out$PH_satisfied <- is.finite(p_global) && p_global > 0.05
  if (!is.finite(p_global) || p_global > 0.05) return(out)



  tt_formula <- stats::update.formula(
    stats::formula(fit), paste(". ~ . + tt(", exposure, ")")
  )
  tt_fit <- tryCatch(
    survival::coxph(
      tt_formula, data = d,
      tt = function(x, t, ...) x * if (time_scale_per_year) t / 365.25 else t,
      ties = "efron", control = survival::coxph.control(iter.max = 100)
    ), error = function(e) NULL
  )
  if (!is.null(tt_fit)) {
    tt_name <- paste0("tt(", exposure, ")")
    tt_coef <- stats::coef(tt_fit)
    tt_sum <- summary(tt_fit)$coefficients
    if (tt_name %in% names(tt_coef) && tt_name %in% rownames(tt_sum)) {
      out$time_interaction_HR_per_year <- unname(exp(tt_coef[[tt_name]]))
      out$time_interaction_P <- unname(tt_sum[tt_name, "Pr(>|z|)"])
    }
  }

  cut_time <- stats::median(d[[time_col]], na.rm = TRUE)
  if (is.finite(cut_time) && cut_time > min(d[[time_col]]) &&
      cut_time < max(d[[time_col]])) {
    split_data <- tryCatch(
      survival::survSplit(
        stats::as.formula(paste0("survival::Surv(", time_col, ", ",
                                 event_col, ") ~ .")),
        data = d, cut = cut_time, episode = "time_period"
      ), error = function(e) NULL
    )
    if (!is.null(split_data) && "time_period" %in% names(split_data)) {
      split_data[, time_period := factor(time_period)]
      split_rhs <- paste(c(paste0(exposure, " * time_period"), covariates),
                         collapse = " + ")
      split_formula <- stats::as.formula(
        paste0("survival::Surv(tstart, ", time_col, ", ", event_col,
               ") ~ ", split_rhs)
      )
      split_fit <- tryCatch(
        survival::coxph(split_formula, data = split_data, ties = "efron",
                        control = survival::coxph.control(iter.max = 100)),
        error = function(e) NULL
      )
      if (!is.null(split_fit)) {
        co <- stats::coef(split_fit)
        inter_name <- paste0(exposure, ":time_period2")
        if (exposure %in% names(co)) {
          out$period1_HR <- unname(exp(co[[exposure]]))
          if (inter_name %in% names(co)) {
            out$period2_HR <- unname(exp(co[[exposure]] + co[[inter_name]]))
          }
        }
        out$split_time_years <- cut_time / 365.25
      }
    }
  }
  out
}

run_triplet <- function(data, time_col, event_col, exposure, demo_covariates,
                        full_covariates) {
  definitions <- list(
    Crude = character(0),
    Demographic_adjusted = demo_covariates,
    Fully_adjusted = full_covariates
  )
  fitted <- lapply(names(definitions), function(nm) {
    fit_cox_model(data, time_col, event_col, exposure,
                  definitions[[nm]], nm)
  })
  names(fitted) <- names(definitions)
  if (any(vapply(fitted, is.null, logical(1)))) return(NULL)
  rows <- lapply(fitted, function(m) {
    ph <- ph_diagnostics(m, time_col, event_col, exposure,
                         covariates = definitions[[m$model]])
    data.table(
      Model = m$model, HR = m$HR, CI_lower = m$CI_lower,
      CI_upper = m$CI_upper, P_raw = m$P_raw,
      N = m$N, Events = m$Events,
      Global_PH_P = ph$Global_PH_P,
      PH_satisfied = ph$PH_satisfied,
      time_interaction_HR_per_year = ph$time_interaction_HR_per_year,
      time_interaction_P = ph$time_interaction_P,
      split_time_years = ph$split_time_years,
      period1_HR = ph$period1_HR,
      period2_HR = ph$period2_HR
    )
  })
  list(results = rbindlist(rows), fits = fitted)
}

make_phewas_wide <- function(long_data) {
  if (nrow(long_data) == 0L) return(data.table())
  long_data[, FDR_P_value := p.adjust(P_raw, method = "BH"), by = Sex_stratum]
  long_data[, HR_CI := sprintf("%.2f (%.2f-%.2f)", HR, CI_lower, CI_upper)]
  wide <- data.table::dcast(
    long_data,
    System + Disease + Code ~ Sex_stratum,
    value.var = c("Total_cases", "No", "HR_CI", "P_raw", "FDR_P_value")
  )
  rename_map <- c(
    Total_cases_Overall = "Overall_Total_cases", No_Overall = "Overall_No.",
    HR_CI_Overall = "Overall_HR (95% CI)", P_raw_Overall = "Overall_P_value",
    FDR_P_value_Overall = "Overall_FDR_P_value",
    Total_cases_Male = "Male_Total_cases", No_Male = "Male_No.",
    HR_CI_Male = "Male_HR (95% CI)", P_raw_Male = "Male_P_value",
    FDR_P_value_Male = "Male_FDR_P_value",
    Total_cases_Female = "Female_Total_cases", No_Female = "Female_No.",
    HR_CI_Female = "Female_HR (95% CI)", P_raw_Female = "Female_P_value",
    FDR_P_value_Female = "Female_FDR_P_value"
  )
  for (old in names(rename_map)) {
    if (old %in% names(wide)) data.table::setnames(wide, old, rename_map[[old]])
  }
  desired <- c(
    "Disease", "Code", "Overall_Total_cases", "Overall_No.",
    "Overall_HR (95% CI)", "Overall_P_value", "Overall_FDR_P_value",
    "Male_Total_cases", "Male_No.", "Male_HR (95% CI)",
    "Male_P_value", "Male_FDR_P_value", "Female_Total_cases", "Female_No.",
    "Female_HR (95% CI)", "Female_P_value", "Female_FDR_P_value", "System"
  )
  data.table::setcolorder(wide, intersect(desired, names(wide)))
  data.table::setorder(wide, System, Disease)
  wide
}

validate_disease_list(analysis_inputs$disease_by_system)
baseline <- normalise_baseline(analysis_inputs$baseline)
sex_strata <- list(Overall = NA_real_, Male = 1, Female = 0)
all_model_results <- list()
all_phewas_rows <- list()
result_index <- 0L

for (system_name in names(analysis_inputs$disease_by_system)) {
  disease_list <- analysis_inputs$disease_by_system[[system_name]]
  system_key <- safe_key(system_name)
  message("Processing system: ", system_name)

  system_data <- copy(baseline)
  disease_map <- list()
  for (disease_name in names(disease_list)) {
    disease_key <- safe_key(disease_name)
    code <- paste(system_key, disease_key, sep = "_")
    time_col <- paste0("outcome_", length(disease_map) + 1L, "_time")
    event_col <- paste0("outcome_", length(disease_map) + 1L, "_event")
    outcome <- prepare_outcome(disease_list[[disease_name]], baseline)
    data.table::setnames(outcome, c("time_days", "event"), c(time_col, event_col))
    system_data <- merge(system_data, outcome, by = "eid",
                         all.x = TRUE, sort = FALSE)
    disease_map[[length(disease_map) + 1L]] <- list(
      system = system_name, disease = gsub("_", " ", disease_name),
      code = code, time_col = time_col, event_col = event_col
    )
  }

  for (sex_name in names(sex_strata)) {
    sex_value <- sex_strata[[sex_name]]
    dt_stratum <- if (is.na(sex_value)) copy(system_data) else
      system_data[sex == sex_value]
    if (nrow(dt_stratum) < MIN_STRATUM_SIZE) next

    if (sex_name == "Overall") {
      demo_covariates <- c("age", "sex", "ethnicity")
      full_covariates <- c("age", "sex", "ethnicity", "smoke", "alcohol",
                           "IPAQ", "BMI", "TDI", "hypertension_status", "DM_status")
    } else {
      demo_covariates <- c("age", "ethnicity")
      full_covariates <- c("age", "ethnicity", "smoke", "alcohol", "IPAQ",
                           "BMI", "TDI", "hypertension_status", "DM_status")
    }

    keep_diseases <- list()
    for (disease in disease_map) {
      time_col <- disease$time_col
      event_col <- disease$event_col
      valid <- dt_stratum[!is.na(get(time_col)) & get(time_col) > 0 &
                            !is.na(get(event_col)) & get(event_col) %in% c(0, 1)]
      n_events <- sum(valid[[event_col]] == 1L)
      prevalence <- n_events / nrow(dt_stratum)
      if (prevalence < PREVALENCE_THRESHOLD) next
      keep_diseases[[length(keep_diseases) + 1L]] <- disease

      exposure_valid <- valid[!is.na(get(TARGET_ALLERGY))]
      total_cases <- sum(exposure_valid[[event_col]] == 1L)
      exposed_cases <- sum(exposure_valid[[event_col]] == 1L &
                             exposure_valid[[TARGET_ALLERGY]] == 1L)
      fitted <- run_triplet(
        valid, time_col, event_col, TARGET_ALLERGY,
        demo_covariates, full_covariates
      )
      if (is.null(fitted)) next

      result_index <- result_index + 1L
      model_rows <- copy(fitted$results)
      model_rows[, `:=`(System = system_name, Disease = disease$disease,
                        Code = disease$code, Sex_stratum = sex_name,
                        Exposure = TARGET_ALLERGY)]
      all_model_results[[result_index]] <- model_rows

      full_row <- model_rows[Model == "Fully_adjusted"]
      if (nrow(full_row) == 1L) {
        all_phewas_rows[[length(all_phewas_rows) + 1L]] <- data.table(
          System = system_name, Disease = disease$disease, Code = disease$code,
          Sex_stratum = sex_name, Total_cases = total_cases,
          No = exposed_cases, HR = full_row$HR,
          CI_lower = full_row$CI_lower, CI_upper = full_row$CI_upper,
          P_raw = full_row$P_raw
        )
      }
    }



    if (length(keep_diseases) > 0L) {
      event_cols <- vapply(keep_diseases, `[[`, character(1), "event_col")
      time_cols <- vapply(keep_diseases, `[[`, character(1), "time_col")
      required_outcomes <- c(event_cols, time_cols)
      comp <- copy(dt_stratum)
      comp <- comp[complete.cases(comp[, ..required_outcomes])]
      for (ec in event_cols) comp <- comp[get(ec) %in% c(0, 1)]
      for (tc in time_cols) comp <- comp[is.finite(get(tc)) & get(tc) > 0]
      if (nrow(comp) > 0L) {
        event_matrix <- as.matrix(comp[, ..event_cols])
        time_matrix <- as.matrix(comp[, ..time_cols])
        comp_event <- as.integer(rowSums(event_matrix) >= 1L)
        event_time_matrix <- time_matrix
        event_time_matrix[event_matrix != 1L] <- Inf
        first_event_time <- matrixStats::rowMins(event_time_matrix)
        last_followup_time <- matrixStats::rowMaxs(time_matrix)
        comp[, composite_event := comp_event]
        comp[, composite_time := ifelse(comp_event == 1L,
                                         first_event_time, last_followup_time)]
        comp <- comp[is.finite(composite_time) & composite_time > 0]
        event_rate <- if (nrow(comp) > 0L) mean(comp$composite_event) else NA_real_
        if (nrow(comp) >= 20L && is.finite(event_rate) &&
            event_rate >= PREVALENCE_THRESHOLD && event_rate <= 0.995) {
          composite_code <- paste0(system_key, "_Any")
          composite_disease <- paste0(system_name, " disease composite")
          comp_total_cases <- sum(comp$composite_event == 1L)
          comp_exposed_cases <- sum(comp$composite_event == 1L &
                                      comp[[TARGET_ALLERGY]] == 1L, na.rm = TRUE)
          comp_fit <- run_triplet(
            comp, "composite_time", "composite_event", TARGET_ALLERGY,
            demo_covariates, full_covariates
          )
          if (!is.null(comp_fit)) {
            result_index <- result_index + 1L
            comp_rows <- copy(comp_fit$results)
            comp_rows[, `:=`(System = system_name, Disease = composite_disease,
                             Code = composite_code, Sex_stratum = sex_name,
                             Exposure = TARGET_ALLERGY)]
            all_model_results[[result_index]] <- comp_rows
            comp_full <- comp_rows[Model == "Fully_adjusted"]
            if (nrow(comp_full) == 1L) {
              all_phewas_rows[[length(all_phewas_rows) + 1L]] <- data.table(
                System = system_name, Disease = composite_disease,
                Code = composite_code, Sex_stratum = sex_name,
                Total_cases = comp_total_cases, No = comp_exposed_cases,
                HR = comp_full$HR, CI_lower = comp_full$CI_lower,
                CI_upper = comp_full$CI_upper, P_raw = comp_full$P_raw
              )
            }
          }
        }
      }
    }
  }
}

if (length(all_model_results) == 0L || length(all_phewas_rows) == 0L) {
  stop("No valid Cox models were produced. Check input completeness and event coding.")
}

model_results <- rbindlist(all_model_results, fill = TRUE)
phewas_long <- rbindlist(all_phewas_rows, fill = TRUE)
phewas_wide <- make_phewas_wide(phewas_long)
significant_overall <- if ("Overall_FDR_P_value" %in% names(phewas_wide)) {
  phewas_wide[!is.na(Overall_FDR_P_value) & Overall_FDR_P_value < 0.05]
} else data.table()



dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_xlsx(as.data.frame(model_results),
           file.path(OUTPUT_DIR, paste0(TARGET_ALLERGY, "_Cox_model_details.xlsx")))
write_xlsx(as.data.frame(phewas_wide),
           file.path(OUTPUT_DIR, paste0(TARGET_ALLERGY, "_Cox_PheWAS_all.xlsx")))
write_xlsx(as.data.frame(significant_overall),
           file.path(OUTPUT_DIR, paste0(TARGET_ALLERGY, "_Cox_PheWAS_overall_FDR_lt_0.05.xlsx")))

message("Conventional Cox analysis completed.")
message("Models: ", nrow(model_results), "; PheWAS rows: ", nrow(phewas_wide),
        "; Overall FDR-significant rows: ", nrow(significant_overall))
