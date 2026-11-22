# Input schema

## 1. Common R object

The common in-memory interface is `analysis_inputs`.

```r
analysis_inputs <- list(
  baseline = baseline_table,
  disease_by_system = disease_outcome_groups,
  omics = list(
    proteomics = proteomics_table,
    metabolomics = metabolomics_table,
    inflammation = inflammation_table,
    clinical = clinical_omics_table
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
```

Only include the omics tables needed for a given analysis. The example uses synthetic identifiers and values only.

## 2. Baseline table

One row per participant. Required fields:

| Field | Type | Meaning |
|---|---|---|
| `eid` | character or integer-like | Participant key; converted to character internally |
| `any_allergy` | 0/1 | Binary exposure |
| `age` | numeric | Age covariate |
| `sex` | 0/1 | Sex coding expected by the Cox and PSM scripts |
| `BMI` | numeric | BMI covariate |
| `smoke` | categorical | Smoking covariate |
| `alcohol` | categorical | Alcohol covariate |
| `ethnicity` | categorical | Ethnicity covariate |
| `TDI` | numeric | Townsend deprivation index |
| `IPAQ` | categorical | Physical activity covariate |
| `hypertension_status` | categorical/binary | Baseline hypertension status |
| `DM_status` | categorical/binary | Baseline diabetes status |
| `baseline_date` | Date or ISO date string | Baseline assessment date |

Optional fields:

| Field | Type | Used by |
|---|---|---|
| `death_date` | Date or missing | Cox censoring and mediation Death outcome |
| `lost_date` | Date or missing | Cox censoring |
| `any_allergy_first_date` | Date or missing | Time-varying Cox |

The mediation script creates `base_date` from `baseline_date` if `base_date` is not supplied.

## 3. Disease outcome groups

`disease_by_system` is a named list. Each system contains named disease outcome tables. The repository contains no real disease group names or outcome data.

```r
disease_by_system <- list(
  System_A = list(
    Outcome_A = data.frame(
      eid = c("SYN001", "SYN002"),
      time = c(365, 730),
      status = c(1L, 0L)
    )
  )
)
```

Each outcome table must contain:

- `eid`: unique participant key.
- `time`: follow-up time. Cox scripts expect days.
- `status`: event indicator coded 0/1.

For mediation, disease tables are flattened by disease label at runtime. Disease labels must be unique across systems and must match the labels supported by the mediation script. The Death outcome is derived from `death_date` in the baseline table.

## 4. Omics tables

Each omics table is a wide table with one row per participant and `eid` as the key. Feature columns must be numeric.

| `analysis_inputs$omics` key | Minimum fields | Used by |
|---|---|---|
| `proteomics` | `eid` plus protein feature columns | 04, 05, 09 |
| `metabolomics` | `eid` plus metabolite feature columns | 04, 05, 09 |
| `inflammation` | `eid`, `WBC`, `PLT`, `NLR`, `CRP`, plus optional mediator features | 04, 06, 09 |
| `clinical` | `eid` plus clinical feature columns | 09 |

The baseline table supplies exposure and covariates. Scripts 05 and 06 merge these baseline fields with the selected omics table in memory using `eid`.

For 04, `config$omics_type` selects `proteomics`, `metabolomics`, or `inflammation`. For 09, the script's `OMICS_TYPE` maps to the corresponding key above.

## 5. Configuration

- `config$omics_type`: `proteomics`, `metabolomics`, or `inflammation` for mediation.
- `config$outcome_time_unit`: `days` or `years` for mediation. Cox input `time` remains in days.
- `config$cutoff_date`: optional administrative cutoff date.
- `config$diseases`: optional vector of mediation outcome labels. If omitted, labels are selected from the supplied outcomes and the supported mediation outcome list.

## 6. File-based interfaces

The PLINK2 and MR scripts remain file-based and use environment variables. Paths are supplied only at runtime in the secure analysis environment. Never commit those files or the actual values of the environment variables.
