# Repository audit summary

## Scope

This audit covers the nine source scripts uploaded for this code package and the minimal interface/documentation changes required to standardize their inputs. The analysis methods, model families, thresholds, and parameter settings were not intentionally redesigned.

## Repository contents

- 01-03: Cox primary, PSM sensitivity, and time-varying Cox scripts. Their participant key was standardized to `eid` to match the common input schema.
- 04: mediation script. A runtime adapter now reads the common `analysis_inputs` object, flattens the supplied outcome groups, and selects the requested omics modality from `analysis_inputs$omics`.
- 05: proteomic/metabolomic signature script. A runtime adapter merges the common baseline and selected omics table in memory.
- 06: inflammation score script. A runtime adapter merges the common baseline and inflammation table in memory.
- 07: PLINK2 GWAS shell pipeline. Its environment-variable/file interface is retained.
- 08: GWAS/MR and target annotation script. Its environment-variable/file interface is retained.
- 09: omics association and enrichment script. It now selects the relevant omics table from the common input list.

## Static checks completed

- No personal Windows drive paths or previously supplied server absolute paths were found in the nine scripts.
- No actual participant-level input files, disease outcome tables, disease-by-system mappings, genotype files, or pQTL files are included.
- No `participant_id` key remains in the scripts; `eid` is the documented canonical in-memory key.
- R source files passed delimiter/string-balance checks.
- The PLINK2 shell script passed `bash -n`.
- `.gitignore` excludes common participant-level, genomic, result, and audit-output file types.
- The synthetic input example uses generated identifiers and values only.

## Not verified

- R and Rscript were not available in the build environment, so no R parser or runtime execution test was performed.
- No analysis was run against UK Biobank data. Numerical reproducibility has not been established.
- The UK Biobank Git Audit Tool has not been run on the final Git repository and commit history.
- The authors have not yet selected a software license.
- GitHub publication, Zenodo archival, and DOI minting have not been performed.

## Important operational note

Scripts 05 and 06 generate participant-level score outputs when run. Those outputs must remain inside the authorized secure analysis environment and must never be added to this repository or uploaded to Zenodo. The same restriction applies to any PLINK2 intermediate files, phenotype/covariate files, disease outcomes, pQTL data, or other participant-level derivatives.
