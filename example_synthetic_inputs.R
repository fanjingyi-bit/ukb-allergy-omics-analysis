required_packages <- c("data.table")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install data.table before running this example.")
suppressPackageStartupMessages(library(data.table))

set.seed(2026)
n <- 120L
eid <- sprintf("SYN%04d", seq_len(n))
exposure <- rep(c(0L, 1L), each = n / 2L)
baseline_date <- as.Date("2010-01-01") + sample(0:30, n, replace = TRUE)

baseline <- data.table(
  eid = eid,
  any_allergy = exposure,
  age = round(rnorm(n, 55, 8), 1),
  sex = sample(c(0L, 1L), n, replace = TRUE),
  BMI = round(rnorm(n, 27, 4), 1),
  smoke = factor(sample(c("Never", "Former", "Current"), n, replace = TRUE)),
  alcohol = factor(sample(c("No", "Yes"), n, replace = TRUE)),
  ethnicity = factor(sample(c("Group_A", "Group_B"), n, replace = TRUE)),
  TDI = round(rnorm(n, 0, 2), 2),
  IPAQ = factor(sample(c("Low", "Moderate", "High"), n, replace = TRUE)),
  hypertension_status = factor(sample(c(0, 1), n, replace = TRUE)),
  DM_status = factor(sample(c(0, 1), n, replace = TRUE)),
  baseline_date = baseline_date,
  death_date = as.Date(NA),
  lost_date = as.Date(NA),
  any_allergy_first_date = as.Date(NA)
)

exposed_ids <- which(exposure == 1L)
baseline$any_allergy_first_date[exposed_ids] <- ifelse(
  seq_along(exposed_ids) %% 2L == 0L,
  as.character(baseline_date[exposed_ids] - sample(30:900, length(exposed_ids), replace = TRUE)),
  as.character(baseline_date[exposed_ids] + sample(30:900, length(exposed_ids), replace = TRUE))
)
baseline$any_allergy_first_date <- as.Date(baseline$any_allergy_first_date, origin = "1970-01-01")

anxiety_outcome <- data.table(
  eid = eid,
  time = sample(180:3000, n, replace = TRUE),
  status = rbinom(n, 1, 0.25)
)

make_omics <- function(prefix, feature_count) {
  values <- replicate(feature_count, rnorm(n))
  values <- as.data.table(values)
  setnames(values, paste0(prefix, "_", seq_len(feature_count)))
  values[, eid := eid]
  setcolorder(values, c("eid", setdiff(names(values), "eid")))
  values
}

proteomics <- make_omics("protein", 8L)
metabolomics <- make_omics("metabolite", 8L)
clinical_omics <- make_omics("clinical_feature", 5L)
inflammation <- data.table(
  eid = eid,
  WBC = rgamma(n, 8, 1),
  PLT = rgamma(n, 250, 1),
  NLR = rgamma(n, 2, 1),
  CRP = rgamma(n, 3, 1)
)

analysis_inputs <- list(
  baseline = baseline,
  disease_by_system = list(
    Synthetic_System = list(Anxiety = anxiety_outcome)
  ),
  omics = list(
    proteomics = proteomics,
    metabolomics = metabolomics,
    inflammation = inflammation,
    clinical = clinical_omics
  ),
  config = list(
    omics_type = "proteomics",
    outcome_time_unit = "days",
    cutoff_date = as.Date("2023-03-31")
  ),
  annotation_maps = list(
    proteomics = NULL,
    metabolomics = NULL
  )
)

stopifnot(
  is.list(analysis_inputs),
  all(c("baseline", "disease_by_system", "omics", "config") %in% names(analysis_inputs)),
  "eid" %in% names(analysis_inputs$baseline),
  all(c("eid", "time", "status") %in% names(analysis_inputs$disease_by_system$Synthetic_System$Anxiety)),
  all(c("WBC", "PLT", "NLR", "CRP") %in% names(analysis_inputs$omics$inflammation))
)

message("Synthetic analysis_inputs object created. No real participant data are used.")
