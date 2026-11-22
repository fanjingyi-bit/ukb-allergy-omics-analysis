






















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
ALLERGY_DATE_COLUMN <- "any_allergy_first_date"
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
  required <- c("eid", TARGET_ALLERGY, ALLERGY_DATE_COLUMN,
                "age", "sex", "ethnicity", "smoke", "alcohol", "IPAQ",
                "BMI", "TDI", "hypertension_status", "DM_status", "baseline_date")
  assert_columns(dt, required, "analysis_inputs$baseline")
  if (anyDuplicated(dt$eid)) {
    stop("Baseline input must contain exactly one row per participant.")
  }
  if (!"death_date" %in% names(dt)) dt[, death_date := as.Date(NA)]
  if (!"lost_date" %in% names(dt)) dt[, lost_date := as.Date(NA)]
  dt[, baseline_date := as_date_safe(baseline_date)]
  dt[, death_date := as_date_safe(death_date)]
  dt[, lost_date := as_date_safe(lost_date)]
  dt[, (ALLERGY_DATE_COLUMN) := as_date_safe(get(ALLERGY_DATE_COLUMN))]
  dt[, sex := suppressWarnings(as.numeric(as.character(sex)))]
  numeric_covariates <- c("age", "BMI", "TDI")
  dt[, (numeric_covariates) := lapply(.SD, function(v) suppressWarnings(as.numeric(as.character(v)))),
     .SDcols = numeric_covariates]
  dt[, (TARGET_ALLERGY) := suppressWarnings(as.numeric(as.character(get(TARGET_ALLERGY))))]
  if (any(!is.na(dt$sex) & !dt$sex %in% c(0, 1))) {
    stop("`sex` must be coded 0 (female) or 1 (male).")
  }
  if (any(!is.na(dt[[TARGET_ALLERGY]]) & !dt[[TARGET_ALLERGY]] %in% c(0, 1))) {
    stop("The exposure indicator must be coded 0/1.")
  }
  dt[, censor_date := pmin(death_date, lost_date, ADMIN_CUTOFF_DATE, na.rm = TRUE)]
  dt[, max_follow_days := as.numeric(censor_date - baseline_date)]
  dt[is.na(baseline_date) | is.na(max_follow_days) | max_follow_days <= 0,
     max_follow_days := NA_real_]
  dt
}

validate_disease_list <- function(x) {
  if (!is.list(x) || length(x) == 0L || is.null(names(x)) || any(names(x) == "")) {
    stop("`disease_by_system` must be a non-empty named list of named disease lists.")
  }
  for (sys in names(x)) {
    if (!is.list(x[[sys]]) || length(x[[sys]]) == 0L ||
        is.null(names(x[[sys]])) || any(names(x[[sys]]) == "")) {
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

make_tv_intervals <- function(baseline, covariates) {
  d <- copy(as.data.table(baseline))
  d <- d[!is.na(max_follow_days) & max_follow_days > 0]
  missing_onset_date <- sum(d[[TARGET_ALLERGY]] == 1 &
                              is.na(d[[ALLERGY_DATE_COLUMN]]), na.rm = TRUE)
  if (missing_onset_date > 0L) {
    warning(missing_onset_date,
            " participants have exposure=1 but no first-exposure date; the source date-based rule classifies missing dates as never exposed. Verify this before interpreting the time-varying analysis.")
  }
  d[, exposure_group := data.table::fifelse(
    is.na(get(ALLERGY_DATE_COLUMN)), "never",
    data.table::fifelse(get(ALLERGY_DATE_COLUMN) <= baseline_date,
                        "baseline_exposed", "incident_exposed")
  )]

  static <- d[exposure_group %in% c("baseline_exposed", "never")]
  static_out <- static[, c("eid", covariates), with = FALSE]
  static_out[, `:=`(tstart = rep(0, nrow(static)),
                    tstop = static$max_follow_days,
                    allergy_tv = as.integer(static$exposure_group == "baseline_exposed"))]

  incident <- d[exposure_group == "incident_exposed"]
  if (nrow(incident) > 0L) {
    incident[, allergy_days := as.numeric(
      get(ALLERGY_DATE_COLUMN) - baseline_date
    )]
    seg1 <- incident[, c("eid", covariates), with = FALSE]
    seg1[, `:=`(tstart = 0,
                tstop = pmin(incident$allergy_days, incident$max_follow_days),
                allergy_tv = 0L)]
    incident_exposed <- incident[allergy_days < max_follow_days]
    seg2 <- incident_exposed[, c("eid", covariates), with = FALSE]
    seg2[, `:=`(tstart = incident_exposed$allergy_days,
                tstop = incident_exposed$max_follow_days,
                allergy_tv = rep(1L, nrow(incident_exposed)))]
    interval_data <- rbindlist(list(static_out, seg1, seg2), fill = TRUE)
  } else {
    interval_data <- static_out
  }

  data.table::setcolorder(interval_data,
                          c("eid", "tstart", "tstop", "allergy_tv", covariates))
  interval_data <- interval_data[is.finite(tstart) & is.finite(tstop) &
                                   tstart < tstop & allergy_tv %in% c(0L, 1L)]
  interval_data
}

attach_outcome_to_intervals <- function(intervals, outcome_data) {
  d <- merge(intervals, outcome_data, by = "eid", all = FALSE,
             sort = FALSE)
  d <- d[!is.na(time_days) & is.finite(time_days) & time_days > 0 &
           !is.na(event) & event %in% c(0L, 1L)]
  original_stop <- d$tstop
  d[, tstop := pmin(tstop, time_days)]
  d[, event_interval := as.integer(event == 1L & original_stop >= time_days &
                                     tstart < time_days)]
  d <- d[tstart < tstop]
  d[, event := event_interval]
  d[, c("time_days", "event_interval") := NULL]
  d
}

fit_tv_model <- function(data, event_col = "event", covariates = character(),
                         model_name = "Model") {
  needed <- unique(c("tstart", "tstop", event_col, "allergy_tv", covariates))
  d <- copy(as.data.table(data))
  d <- d[complete.cases(d[, ..needed])]
  d <- d[is.finite(tstart) & is.finite(tstop) & tstart < tstop &
           get(event_col) %in% c(0, 1)]
  if (nrow(d) < 20L || length(unique(d$allergy_tv)) < 2L ||
      sum(d[[event_col]] == 1L) == 0L) return(NULL)
  rhs <- paste(c("allergy_tv", covariates), collapse = " + ")
  f <- stats::as.formula(paste0("survival::Surv(tstart, tstop, ",
                                event_col, ") ~ ", rhs))
  fit <- tryCatch(
    survival::coxph(f, data = d, ties = "efron", x = TRUE, model = TRUE,
                    control = survival::coxph.control(iter.max = 200)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  b <- stats::coef(fit)
  ci <- tryCatch(stats::confint(fit), error = function(e) NULL)
  sm <- summary(fit)
  if (is.null(ci) || !"allergy_tv" %in% names(b) ||
      !"allergy_tv" %in% rownames(sm$coefficients) ||
      any(!is.finite(b)) || any(!is.finite(ci))) return(NULL)
  list(
    fit = fit, data = d, model = model_name,
    HR = unname(exp(b[["allergy_tv"]])),
    CI_lower = unname(exp(ci["allergy_tv", 1L])),
    CI_upper = unname(exp(ci["allergy_tv", 2L])),
    P_raw = unname(sm$coefficients["allergy_tv", "Pr(>|z|)"]),
    N_intervals = nrow(d), N_participants = data.table::uniqueN(d$eid),
    Events = sum(d[[event_col]] == 1L)
  )
}

ph_diagnostics_tv <- function(model_object, event_col = "event", covariates = character()) {
  out <- list(Global_PH_P = NA_real_, PH_satisfied = NA,
              time_interaction_HR_per_year = NA_real_,
              time_interaction_P = NA_real_, split_time_years = NA_real_,
              period1_HR = NA_real_, period2_HR = NA_real_)
  if (is.null(model_object)) return(out)
  fit <- model_object$fit
  d <- model_object$data
  zph <- tryCatch(survival::cox.zph(fit, transform = "rank"),
                  error = function(e) NULL)
  if (is.null(zph) || !"GLOBAL" %in% rownames(zph$table)) return(out)
  p_global <- unname(zph$table["GLOBAL", "p"])
  out$Global_PH_P <- p_global
  out$PH_satisfied <- is.finite(p_global) && p_global > 0.05
  if (!is.finite(p_global) || p_global > 0.05) return(out)



  tt_formula <- stats::as.formula(paste(
    paste(deparse(stats::formula(fit)), collapse = ""),
    "+ tt(allergy_tv)"
  ))
  tt_fit <- tryCatch(
    survival::coxph(
      tt_formula, data = d,
      tt = function(x, t, ...) x * (t / 365.25),
      ties = "efron", control = survival::coxph.control(iter.max = 200)
    ), error = function(e) NULL
  )
  if (!is.null(tt_fit)) {
    tt_coef <- stats::coef(tt_fit)
    tt_sum <- summary(tt_fit)$coefficients
    tt_name <- "tt(allergy_tv)"
    if (tt_name %in% names(tt_coef) && tt_name %in% rownames(tt_sum)) {
      out$time_interaction_HR_per_year <- unname(exp(tt_coef[[tt_name]]))
      out$time_interaction_P <- unname(tt_sum[tt_name, "Pr(>|z|)"])
    }
  }

  cut_time <- stats::median(d$tstop, na.rm = TRUE)
  if (is.finite(cut_time) && cut_time > min(d$tstart) && cut_time < max(d$tstop)) {
    d_for_split <- copy(d)
    d_for_split[, split_event := get(event_col)]
    split_data <- tryCatch(
      survival::survSplit(
        survival::Surv(tstart, tstop, split_event) ~ .,
        data = d_for_split, cut = cut_time, episode = "time_period"
      ), error = function(e) NULL
    )
    if (!is.null(split_data) && "time_period" %in% names(split_data)) {
      split_data[, time_period := factor(time_period)]
      split_rhs <- paste(c("allergy_tv * time_period", covariates), collapse = " + ")
      split_formula <- stats::as.formula(
        paste0("survival::Surv(tstart, tstop, split_event) ~ ", split_rhs)
      )
      split_fit <- tryCatch(
        survival::coxph(split_formula, data = split_data, ties = "efron",
                        control = survival::coxph.control(iter.max = 200)),
        error = function(e) NULL
      )
      if (!is.null(split_fit)) {
        co <- stats::coef(split_fit)
        interaction_name <- "allergy_tv:time_period2"
        if ("allergy_tv" %in% names(co)) {
          out$period1_HR <- unname(exp(co[["allergy_tv"]]))
          if (interaction_name %in% names(co)) {
            out$period2_HR <- unname(exp(co[["allergy_tv"]] + co[[interaction_name]]))
          }
        }
        out$split_time_years <- cut_time / 365.25
      }
    }
  }
  out
}

run_tv_triplet <- function(data, demo_covariates, full_covariates) {
  definitions <- list(
    Crude = character(0),
    Demographic_adjusted = demo_covariates,
    Fully_adjusted = full_covariates
  )
  fitted <- lapply(names(definitions), function(nm) {
    fit_tv_model(data, covariates = definitions[[nm]], model_name = nm)
  })
  names(fitted) <- names(definitions)
  if (any(vapply(fitted, is.null, logical(1)))) return(NULL)
  rows <- lapply(fitted, function(m) {
    ph <- ph_diagnostics_tv(m, covariates = definitions[[m$model]])
    data.table(
      Model = m$model, HR = m$HR, CI_lower = m$CI_lower,
      CI_upper = m$CI_upper, P_raw = m$P_raw,
      N_intervals = m$N_intervals, N_participants = m$N_participants,
      Events = m$Events, Global_PH_P = ph$Global_PH_P,
      PH_satisfied = ph$PH_satisfied,
      time_interaction_HR_per_year = ph$time_interaction_HR_per_year,
      time_interaction_P = ph$time_interaction_P,
      split_time_years = ph$split_time_years,
      period1_HR = ph$period1_HR, period2_HR = ph$period2_HR
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



interval_covariates <- c("age", "sex", "ethnicity", "smoke", "alcohol", "IPAQ",
                         "BMI", "TDI", "hypertension_status", "DM_status")
sex_strata <- list(Overall = NA_real_, Male = 1, Female = 0)
all_model_results <- list()
all_phewas_rows <- list()
model_index <- 0L

for (system_name in names(analysis_inputs$disease_by_system)) {
  system_key <- safe_key(system_name)
  disease_list <- analysis_inputs$disease_by_system[[system_name]]
  message("Processing system: ", system_name)


  prepared_outcomes <- list()
  disease_map <- list()
  for (disease_name in names(disease_list)) {
    disease_key <- safe_key(disease_name)
    disease_code <- paste(system_key, disease_key, sep = "_")
    prepared_outcomes[[disease_name]] <- prepare_outcome(
      disease_list[[disease_name]], baseline
    )
    disease_map[[disease_name]] <- list(
      system = system_name, disease = gsub("_", " ", disease_name),
      code = disease_code
    )
  }

  for (sex_name in names(sex_strata)) {
    sex_value <- sex_strata[[sex_name]]
    dt_stratum <- if (is.na(sex_value)) copy(baseline) else
      baseline[sex == sex_value]
    if (nrow(dt_stratum) < MIN_STRATUM_SIZE) next

    tv_intervals <- make_tv_intervals(dt_stratum, interval_covariates)
    if (nrow(tv_intervals) == 0L) next
    if (sex_name == "Overall") {
      demo_covariates <- c("age", "sex", "ethnicity")
      full_covariates <- interval_covariates
    } else {
      demo_covariates <- c("age", "ethnicity")
      full_covariates <- c("age", "ethnicity", "smoke", "alcohol", "IPAQ",
                           "BMI", "TDI", "hypertension_status", "DM_status")
    }

    keep_diseases <- character(0)
    system_outcomes <- data.table(eid = dt_stratum$eid)

    for (disease_name in names(prepared_outcomes)) {
      outcome <- copy(prepared_outcomes[[disease_name]])
      outcome <- outcome[eid %in% dt_stratum$eid]
      valid_outcome <- outcome[!is.na(time_days) & is.finite(time_days) &
                                 time_days > 0 & !is.na(event) & event %in% c(0, 1)]
      prevalence <- sum(valid_outcome$event == 1L) / nrow(dt_stratum)
      if (!is.finite(prevalence) || prevalence < PREVALENCE_THRESHOLD) next
      keep_diseases <- c(keep_diseases, disease_name)

      interval_outcome <- attach_outcome_to_intervals(tv_intervals, valid_outcome)
      if (nrow(interval_outcome) == 0L) next
      total_cases <- sum(valid_outcome$event == 1L)


      exposed_cases <- sum(interval_outcome$event == 1L &
                             interval_outcome$allergy_tv == 1L)
      fitted <- run_tv_triplet(interval_outcome, demo_covariates, full_covariates)
      if (is.null(fitted)) next

      disease <- disease_map[[disease_name]]
      model_index <- model_index + 1L
      model_rows <- copy(fitted$results)
      model_rows[, `:=`(System = system_name, Disease = disease$disease,
                        Code = disease$code, Sex_stratum = sex_name,
                        Exposure = TARGET_ALLERGY)]
      all_model_results[[model_index]] <- model_rows

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

      disease_key <- safe_key(disease_name)
      time_col <- paste0(disease_key, "_time")
      event_col <- paste0(disease_key, "_event")
      data.table::setnames(outcome, c("time_days", "event"), c(time_col, event_col))
      system_outcomes <- merge(
        system_outcomes,
        outcome[, c("eid", time_col, event_col), with = FALSE],
        by = "eid", all.x = TRUE, sort = FALSE
      )
    }




    if (length(keep_diseases) > 0L) {
      time_cols <- paste0(vapply(keep_diseases, safe_key, character(1)), "_time")
      event_cols <- paste0(vapply(keep_diseases, safe_key, character(1)), "_event")
      required_cols <- c(time_cols, event_cols)
      if (all(required_cols %in% names(system_outcomes))) {
        comp <- copy(system_outcomes)
        comp <- comp[complete.cases(comp[, ..required_cols])]
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
          if (nrow(comp) > 0L) {
            comp_intervals <- attach_outcome_to_intervals(
              tv_intervals,
              comp[, .(eid, time_days = composite_time,
                       event = composite_event)]
            )
            if (nrow(comp_intervals) > 0L) {
              comp_total_cases <- sum(comp$composite_event == 1L)
              comp_exposed_cases <- sum(comp_intervals$event == 1L &
                                          comp_intervals$allergy_tv == 1L)
              comp_fit <- run_tv_triplet(comp_intervals, demo_covariates, full_covariates)
              if (!is.null(comp_fit)) {
                composite_code <- paste0(system_key, "_Any")
                composite_disease <- paste0(system_name, " disease composite")
                model_index <- model_index + 1L
                comp_rows <- copy(comp_fit$results)
                comp_rows[, `:=`(System = system_name,
                                 Disease = composite_disease,
                                 Code = composite_code,
                                 Sex_stratum = sex_name,
                                 Exposure = TARGET_ALLERGY)]
                all_model_results[[model_index]] <- comp_rows
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
  }
}

if (length(all_model_results) == 0L || length(all_phewas_rows) == 0L) {
  stop("No valid time-varying Cox models were produced. Check dates and outcome inputs.")
}

model_results <- rbindlist(all_model_results, fill = TRUE)
phewas_long <- rbindlist(all_phewas_rows, fill = TRUE)
phewas_wide <- make_phewas_wide(phewas_long)
significant_overall <- if ("Overall_FDR_P_value" %in% names(phewas_wide)) {
  phewas_wide[!is.na(Overall_FDR_P_value) & Overall_FDR_P_value < 0.05]
} else data.table()



dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_xlsx(as.data.frame(model_results),
           file.path(OUTPUT_DIR, paste0(TARGET_ALLERGY, "_TV_Cox_model_details.xlsx")))
write_xlsx(as.data.frame(phewas_wide),
           file.path(OUTPUT_DIR, paste0(TARGET_ALLERGY, "_TV_Cox_PheWAS_all.xlsx")))
write_xlsx(as.data.frame(significant_overall),
           file.path(OUTPUT_DIR, paste0(TARGET_ALLERGY,
                                        "_TV_Cox_PheWAS_overall_FDR_lt_0.05.xlsx")))

message("Time-varying Cox sensitivity analysis completed.")
message("Models: ", nrow(model_results), "; PheWAS rows: ", nrow(phewas_wide),
        "; Overall FDR-significant rows: ", nrow(significant_overall))
