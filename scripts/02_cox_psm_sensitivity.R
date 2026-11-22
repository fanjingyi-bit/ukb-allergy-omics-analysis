
















required_packages <- c("data.table", "survival", "MatchIt", "cobalt", "writexl")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(MatchIt)
  library(cobalt)
  library(writexl)
})


TARGET_ALLERGY <- "any_allergy"
ADMIN_CUTOFF_DATE <- as.Date("2023-03-31")
PREVALENCE_THRESHOLD <- 0.005
MATCH_RATIO <- 3L
PS_CALIPER <- 0.2
MIN_EXPOSED_FOR_MATCH <- 10L
MIN_UNEXPOSED_FOR_MATCH <- 30L
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
  required <- c("eid", TARGET_ALLERGY, "age", "BMI", "sex",
                "ethnicity", "smoke", "alcohol", "IPAQ", "TDI",
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
  dt[, sex_code := suppressWarnings(as.numeric(as.character(sex)))]
  numeric_covariates <- c("age", "BMI", "TDI")
  dt[, (numeric_covariates) := lapply(.SD, function(v) suppressWarnings(as.numeric(as.character(v)))),
     .SDcols = numeric_covariates]
  if (any(!is.na(dt$sex_code) & !dt$sex_code %in% c(0, 1))) {
    stop("`sex` must be coded 0 (female) or 1 (male).")
  }
  dt[, sex_factor := factor(sex_code, levels = c(0, 1),
                            labels = c("Female", "Male"))]
  dt[, (TARGET_ALLERGY) := suppressWarnings(as.numeric(as.character(get(TARGET_ALLERGY))))]
  if (any(!is.na(dt[[TARGET_ALLERGY]]) & !dt[[TARGET_ALLERGY]] %in% c(0, 1))) {
    stop("The exposure must be coded 0/1.")
  }
  factor_columns <- c("sex_factor", "ethnicity", "smoke", "alcohol", "IPAQ",
                      "hypertension_status", "DM_status")
  dt[, (factor_columns) := lapply(.SD, as.factor), .SDcols = factor_columns]
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
  d[, time_years := time_days / 365.25]
  d[, .(eid, time_years, event)]
}

get_balance_table <- function(match_fit) {
  balance_object <- tryCatch(
    cobalt::bal.tab(match_fit, un = TRUE, stats = "mean.diffs",
                    thresholds = c(m = 0.1)),
    error = function(e) NULL
  )
  if (is.null(balance_object) || is.null(balance_object$Balance)) {
    return(data.table(Covariate = character(), SMD_before = numeric(),
                      SMD_after = numeric()))
  }
  balance <- as.data.frame(balance_object$Balance)
  balance$Covariate <- rownames(balance)
  before_col <- intersect(c("Diff.Un", "Diff"), names(balance))
  after_col <- intersect(c("Diff.Adj", "Diff"), names(balance))
  before <- if (length(before_col)) balance[[before_col[1L]]] else rep(NA_real_, nrow(balance))
  after <- if (length(after_col)) balance[[after_col[1L]]] else rep(NA_real_, nrow(balance))
  out <- data.table(
    Covariate = as.character(balance$Covariate),
    SMD_before = as.numeric(before),
    SMD_after = as.numeric(after)
  )
  out[, `:=`(SMD_before = round(SMD_before, 3),
             SMD_after = round(SMD_after, 3))]
  out
}

fit_conditional_cox <- function(data, time_col, event_col, covariates = character()) {
  required <- unique(c(time_col, event_col, "case", "subclass", covariates))
  d <- copy(as.data.table(data))
  d <- d[complete.cases(d[, ..required])]
  d <- d[is.finite(get(time_col)) & get(time_col) > 0 &
           get(event_col) %in% c(0, 1)]
  if (nrow(d) < 20L || length(unique(d$case)) < 2L ||
      sum(d[[event_col]] == 1L) == 0L) return(NULL)
  rhs <- paste(c("case", covariates, "strata(subclass)"), collapse = " + ")
  f <- as.formula(paste0("survival::Surv(", time_col, ", ", event_col,") ~ ", rhs))
  fit <- tryCatch(
    survival::coxph(f, data = d, ties = "efron", x = TRUE, model = TRUE,
                    control = survival::coxph.control(iter.max = 100)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  b <- stats::coef(fit)
  ci <- tryCatch(stats::confint(fit), error = function(e) NULL)
  sm <- summary(fit)
  if (is.null(ci) || !"case" %in% names(b) || !"case" %in% rownames(sm$coefficients) ||
      any(!is.finite(b)) || any(!is.finite(ci))) return(NULL)
  list(fit = fit, data = d, HR = unname(exp(b[["case"]])),
       CI_lower = unname(exp(ci["case", 1L])),
       CI_upper = unname(exp(ci["case", 2L])),
       P_raw = unname(sm$coefficients["case", "Pr(>|z|)"]),
       N = nrow(d), Events = sum(d[[event_col]] == 1L))
}

psm_ph_diagnostics <- function(model_object, time_col, event_col,
                               additional_covariates = character()) {
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
    "+ tt(case)"
  ))
  tt_fit <- tryCatch(
    survival::coxph(tt_formula, data = d, tt = function(x, t, ...) x * t,
                    ties = "efron", control = survival::coxph.control(iter.max = 100)),
    error = function(e) NULL
  )
  if (!is.null(tt_fit)) {
    tt_coef <- stats::coef(tt_fit)
    tt_sum <- summary(tt_fit)$coefficients
    tt_name <- "tt(case)"
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
      split_rhs <- paste(c("case * time_period", additional_covariates,
                           "strata(subclass)"), collapse = " + ")
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
        interaction_name <- "case:time_period2"
        if ("case" %in% names(co)) {
          out$period1_HR <- unname(exp(co[["case"]]))
          if (interaction_name %in% names(co)) {
            out$period2_HR <- unname(exp(co[["case"]] + co[[interaction_name]]))
          }
        }
        out$split_time_years <- cut_time
      }
    }
  }
  out
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
strata_config <- list(
  Overall = list(sex_code = NA_real_,
                 match_covariates = c("age", "BMI", "sex_factor", "ethnicity", "smoke",
                                      "alcohol", "IPAQ", "TDI", "hypertension_status", "DM_status"),
                 exact_formula = ~ sex_factor + ethnicity + smoke),
  Male = list(sex_code = 1,
              match_covariates = c("age", "BMI", "ethnicity", "smoke", "alcohol",
                                   "IPAQ", "TDI", "hypertension_status", "DM_status"),
              exact_formula = ~ ethnicity + smoke),
  Female = list(sex_code = 0,
                match_covariates = c("age", "BMI", "ethnicity", "smoke", "alcohol",
                                     "IPAQ", "TDI", "hypertension_status", "DM_status"),
                exact_formula = ~ ethnicity + smoke)
)



matched_sets <- list()
balance_tables <- list()
for (stratum_name in names(strata_config)) {
  cfg <- strata_config[[stratum_name]]
  d <- copy(baseline)
  if (!is.na(cfg$sex_code)) d <- d[sex_code == cfg$sex_code]
  d[, case := get(TARGET_ALLERGY)]
  d <- d[!is.na(case)]
  match_covs <- cfg$match_covariates
  d <- d[complete.cases(d[, ..match_covs]) & is.finite(age) &
           is.finite(BMI) & is.finite(TDI)]
  exposed_n <- sum(d$case == 1L)
  unexposed_n <- sum(d$case == 0L)
  if (exposed_n < MIN_EXPOSED_FOR_MATCH || unexposed_n < MIN_UNEXPOSED_FOR_MATCH) {
    warning("Insufficient sample for matching in stratum: ", stratum_name)
    next
  }
  if (any(vapply(match_covs, function(v) length(unique(d[[v]])) < 2L,
                 logical(1)))) {
    warning("A matching covariate has only one level in stratum: ", stratum_name)
    next
  }
  ps_formula <- stats::as.formula(paste("case ~", paste(match_covs, collapse = " + ")))
  match_fit <- tryCatch(
    MatchIt::matchit(
      formula = ps_formula, data = d, method = "nearest", ratio = MATCH_RATIO,
      estimand = "ATT", distance = "glm", link = "linear.logit",
      exact = cfg$exact_formula, caliper = PS_CALIPER, std.caliper = TRUE,
      discard = "both", replace = FALSE
    ), error = function(e) {
      warning("Matching failed in ", stratum_name, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(match_fit)) next
  matched_data <- as.data.table(MatchIt::match.data(match_fit))
  if (nrow(matched_data) == 0L) next
  matched_sets[[stratum_name]] <- matched_data[, .(eid, subclass)]

  smd <- get_balance_table(match_fit)
  smd[, `:=`(Sex_stratum = stratum_name)]
  balance_values <- smd$SMD_after[is.finite(smd$SMD_after)]
  max_smd <- if (length(balance_values)) max(abs(balance_values)) else NA_real_
  balance_flag <- if (length(balance_values) && all(abs(balance_values) < 0.1)) {
    "Balanced"
  } else {
    "Not balanced or unavailable"
  }
  smd[, `:=`(Balance_status = balance_flag, Max_abs_SMD_after = max_smd)]
  balance_tables[[stratum_name]] <- smd
  message(stratum_name, " matching: ", nrow(matched_data),
          " participants; max |SMD| after matching = ",
          ifelse(is.finite(max_smd), round(max_smd, 3), NA))
}

model_results <- list()
phewas_rows <- list()
model_index <- 0L

for (system_name in names(analysis_inputs$disease_by_system)) {
  system_key <- safe_key(system_name)
  disease_list <- analysis_inputs$disease_by_system[[system_name]]
  message("Processing system: ", system_name)

  for (disease_name in names(disease_list)) {
    disease_key <- safe_key(disease_name)
    disease_code <- paste(system_key, disease_key, sep = "_")
    outcome <- prepare_outcome(disease_list[[disease_name]], baseline)
    full_data <- merge(baseline, outcome, by = "eid", all.x = TRUE,
                       sort = FALSE)
    full_data[, case := get(TARGET_ALLERGY)]
    full_data <- full_data[!is.na(case) & !is.na(event) &
                             event %in% c(0, 1) & is.finite(time_years) & time_years > 0]

    for (stratum_name in names(strata_config)) {
      if (is.null(matched_sets[[stratum_name]])) next
      cfg <- strata_config[[stratum_name]]
      stratum_data <- if (is.na(cfg$sex_code)) copy(full_data) else
        full_data[sex_code == cfg$sex_code]
      if (nrow(stratum_data) == 0L) next
      prevalence <- sum(stratum_data$event == 1L) / nrow(stratum_data)
      if (!is.finite(prevalence) || prevalence < PREVALENCE_THRESHOLD) next

      match_covs <- cfg$match_covariates
      stratum_data <- stratum_data[complete.cases(stratum_data[, ..match_covs])]
      matched_data <- merge(stratum_data, matched_sets[[stratum_name]],
                            by = "eid", all = FALSE, sort = FALSE)
      if (nrow(matched_data) < 40L || length(unique(matched_data$case)) < 2L) next



      candidate_adjusters <- c("alcohol", "IPAQ", "TDI", "hypertension_status")
      adjusters <- candidate_adjusters[vapply(candidate_adjusters, function(v) {
        length(unique(stats::na.omit(matched_data[[v]]))) >= 2L
      }, logical(1))]

      crude_fit <- fit_conditional_cox(matched_data, "time_years", "event")
      adjusted_fit <- fit_conditional_cox(matched_data, "time_years", "event", adjusters)
      if (is.null(crude_fit) || is.null(adjusted_fit)) next

      ph <- psm_ph_diagnostics(crude_fit, "time_years", "event", character(0))
      exposed_events <- sum(matched_data$event == 1L & matched_data$case == 1L)
      total_events <- sum(matched_data$event == 1L)
      case_counts <- matched_data[, .(
        N = .N, Events = sum(event == 1L)
      ), by = case]
      case_event_rate <- if (any(case_counts$case == 1L)) {
        100 * case_counts[case == 1L, Events / N]
      } else NA_real_
      control_event_rate <- if (any(case_counts$case == 0L)) {
        100 * case_counts[case == 0L, Events / N]
      } else NA_real_

      smd <- copy(balance_tables[[stratum_name]])
      smd_max <- if (nrow(smd) && any(is.finite(smd$SMD_after))) {
        max(abs(smd$SMD_after[is.finite(smd$SMD_after)]))
      } else NA_real_
      smd_status <- if (is.finite(smd_max) && smd_max < 0.1) "Balanced" else
        "Not balanced or unavailable"

      model_index <- model_index + 1L
      model_results[[model_index]] <- data.table(
        System = system_name, Disease = disease_name, Code = disease_code,
        Sex_stratum = stratum_name, Exposure = TARGET_ALLERGY,
        N_matched = adjusted_fit$N, Events = total_events,
        Exposed_events = exposed_events,
        Crude_HR = crude_fit$HR, Crude_CI_lower = crude_fit$CI_lower,
        Crude_CI_upper = crude_fit$CI_upper, Crude_P = crude_fit$P_raw,
        Adjusted_HR = adjusted_fit$HR,
        Adjusted_CI_lower = adjusted_fit$CI_lower,
        Adjusted_CI_upper = adjusted_fit$CI_upper,
        Adjusted_P = adjusted_fit$P_raw,
        Adjusted_covariates = paste(adjusters, collapse = ";"),
        Case_event_rate_percent = case_event_rate,
        Control_event_rate_percent = control_event_rate,
        Global_PH_P = ph$Global_PH_P, PH_satisfied = ph$PH_satisfied,
        time_interaction_HR_per_year = ph$time_interaction_HR_per_year,
        time_interaction_P = ph$time_interaction_P,
        split_time_years = ph$split_time_years,
        period1_HR = ph$period1_HR, period2_HR = ph$period2_HR,
        Balance_status = smd_status, Max_abs_SMD_after = smd_max
      )

      phewas_rows[[length(phewas_rows) + 1L]] <- data.table(
        System = system_name, Disease = disease_name, Code = disease_code,
        Sex_stratum = stratum_name, Total_cases = total_events,
        No = exposed_events, HR = adjusted_fit$HR,
        CI_lower = adjusted_fit$CI_lower, CI_upper = adjusted_fit$CI_upper,
        P_raw = adjusted_fit$P_raw
      )
    }
  }
}

if (length(model_results) == 0L || length(phewas_rows) == 0L) {
  stop("No valid matched Cox models were produced. Check matching and outcome inputs.")
}

model_results_dt <- rbindlist(model_results, fill = TRUE)
phewas_long <- rbindlist(phewas_rows, fill = TRUE)
phewas_wide <- make_phewas_wide(phewas_long)
significant_overall <- if ("Overall_FDR_P_value" %in% names(phewas_wide)) {
  phewas_wide[!is.na(Overall_FDR_P_value) & Overall_FDR_P_value < 0.05]
} else data.table()


balance_output <- rbindlist(balance_tables, fill = TRUE)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_xlsx(as.data.frame(model_results_dt),
           file.path(OUTPUT_DIR, paste0(TARGET_ALLERGY, "_PSM_Cox_model_details.xlsx")))
write_xlsx(as.data.frame(balance_output),
           file.path(OUTPUT_DIR, paste0(TARGET_ALLERGY, "_PSM_balance_SMD.xlsx")))
write_xlsx(as.data.frame(phewas_wide),
           file.path(OUTPUT_DIR, paste0(TARGET_ALLERGY, "_PSM_Cox_PheWAS_all.xlsx")))
write_xlsx(as.data.frame(significant_overall),
           file.path(OUTPUT_DIR, paste0(TARGET_ALLERGY,
                                        "_PSM_Cox_PheWAS_overall_FDR_lt_0.05.xlsx")))

message("PSM-Cox sensitivity analysis completed.")
message("Matched Cox models: ", nrow(model_results_dt),
        "; PheWAS rows: ", nrow(phewas_wide),
        "; Overall FDR-significant rows: ", nrow(significant_overall))
