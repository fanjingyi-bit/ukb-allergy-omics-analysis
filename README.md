# UK Biobank Allergy and Multi-omics Analysis Code

This repository contains analysis scripts only. It does not contain UK Biobank participant-level data, phenotype files, disease outcome tables, genotype files, pQTL files, local data paths, credentials, or project-specific disease mappings.

## Repository structure

- `scripts/01_cox_primary.R`: primary Cox analysis.
- `scripts/02_cox_psm_sensitivity.R`: propensity-score-matched Cox sensitivity analysis.
- `scripts/03_cox_time_varying_sensitivity.R`: time-varying exposure Cox sensitivity analysis.
- `scripts/04_omics_mediation.R`: proteomic, inflammatory, and metabolomic mediation analysis.
- `scripts/05_signature_proteomics_metabolomics.R`: proteomic and metabolomic signature construction.
- `scripts/06_inflammation_score.R`: inflammatory score construction.
- `scripts/07_plink2_gwas_pipeline.sh`: PLINK2 genotype QC, relatedness filtering, PCA, GWAS, and clumping.
- `scripts/08_gwas_mr_druggable_targets.R`: GWAS summary-statistic preparation, pQTL MR, and target annotation.
- `scripts/09_omics_association_enrichment.R`: omics association, pathway enrichment, and tissue enrichment.
- `INPUT_SCHEMA.md`: common input object and file-interface specification.
- `example_synthetic_inputs.R`: synthetic, non-UKB example showing the common R input structure.
- `.gitignore`: excludes data and analysis outputs from version control.

## Common R input object

Scripts 01-06 and 09 use a common in-memory object named `analysis_inputs`. Prepare this object inside the authorized analysis environment before running a script. The object has these top-level elements:

- `baseline`: one row per participant, keyed by `eid`.
- `disease_by_system`: a named nested list of disease outcome tables. Actual disease groupings and outcomes are supplied by the authorized user at runtime and are not included in this repository.
- `omics`: a named list containing `proteomics`, `metabolomics`, `inflammation`, and/or `clinical` feature tables. Each table is keyed by `eid`.
- `config`: optional settings such as mediation `omics_type`, `outcome_time_unit`, `cutoff_date`, and `diseases`.
- `annotation_maps`: optional public feature annotation tables for signature plots.

The baseline table uses `eid` as its canonical participant key. The Cox scripts use this same key internally. `baseline_date` is the canonical baseline date; the mediation script creates its internal `base_date` alias from `baseline_date` when needed.

Disease outcome tables use `eid`, `time`, and `status`. For the Cox scripts, `time` is in days. The mediation script reads `config$outcome_time_unit`, which must be explicitly set to `days` or `years` when it differs from its default.

The synthetic example is for understanding the interface only. All identifiers and values in `example_synthetic_inputs.R` are randomly generated; no UK Biobank participant-level data or real participant records were used to create it. Its purpose is to demonstrate the input-object structure without distributing study data. It is not a substitute for approved UK Biobank data and is not intended to reproduce study estimates.

## Example R setup

```r
source("example_synthetic_inputs.R")

# Select one modality for the script being run.
analysis_inputs$config$omics_type <- "proteomics"
analysis_inputs$config$outcome_time_unit <- "days"
analysis_inputs$config$cutoff_date <- as.Date("2023-03-31")

# Run the relevant script in the same R session.
source("scripts/01_cox_primary.R")
```

For the signature script, set `OMICS_TYPE` in `scripts/05_signature_proteomics_metabolomics.R` to `proteomics` or `metabolomics`. For the association/enrichment script, set `OMICS_TYPE` to `proteomic`, `metabolic`, `clinical`, or `inflam`. The corresponding table is selected from `analysis_inputs$omics`.

Scripts 05 and 06 write participant-level score files to `ANALYSIS_OUTPUT_DIR`. These outputs must remain inside the approved secure analysis environment and must not be uploaded to GitHub, Zenodo, or included in a manuscript supplement unless UK Biobank has approved their release.

## PLINK2 and MR input interfaces

The genotype pipeline and summary-statistic MR pipeline use file-based environment-variable interfaces because their inputs are large files or directory-based reference resources. Set paths only at runtime in the approved secure environment. Do not place actual paths, filenames identifying participants, genotype files, phenotype files, covariate files, pQTL files, or generated outputs in this repository.

The PLINK2 script requires `OUTPUT_DIR`, `PLINK1_PATH`, `PLINK2_PATH`, `INPUT_BED_PREFIX`, `PHENO_FILE`, and `COVAR_FILE`. Optional settings include `THREADS`, `MEMORY_MB`, `KING_CUTOFF`, `PHENO_NAME`, `GWAS_PREFIX`, and `GWAS_SEED`.

The MR script requires `GWAS_FILE`, `DBSNP_FILE`, `PQTL_DIR`, `TARGET_REFERENCE_FILE`, `UP_GENES_FILE`, `DOWN_GENES_FILE`, `OUTPUT_DIR`, `GWAS_BUILD`, `DBSNP_BUILD`, and `PQTL_BUILD`. Additional optional environment variables are documented by the script's configuration block.

## Data protection

- Do not commit participant-level data or derived participant-level outputs.
- Do not commit local paths, credentials, access tokens, project identifiers, or private disease mapping files.
- Keep all UK Biobank data and intermediate files within the approved UK Biobank analysis environment.
- Review the complete Git history, not only the latest file versions, before publication.
- Keep the UK Biobank Git Audit Tool output folder (`ukb_audit_reports/`) out of version control.
- Follow current UK Biobank requirements for code repository registration and Git Audit Tool review before making code public.
- Only release aggregate outputs that have passed the applicable UK Biobank disclosure/output checks.

## Reproducibility

Record the R version, package versions, PLINK version, genome build, input schema version, and analysis configuration used for each run. The code repository does not include UK Biobank data, so exact numerical reproduction requires the same approved input data and access permissions.

## Citation and availability

A final Code Availability statement should identify the public repository URL and the archived Zenodo DOI after the repository has been reviewed and released. Do not cite a placeholder DOI or publish the repository before the required UK Biobank review is complete.
