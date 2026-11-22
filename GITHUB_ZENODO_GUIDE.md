# GitHub and Zenodo release guide

## Before publishing

1. Confirm that the scripts and all Git history contain no participant-level data, UK Biobank identifiers, local paths, credentials, private project names, or restricted disease mapping files.
2. Run the UK Biobank Git Audit Tool and complete any required AMS repository registration before making the repository public.
3. Run the scripts in the authorized analysis environment and record software/package versions. The synthetic example checks only the input object structure; it does not validate study estimates.
4. Confirm that the corresponding study outputs are releasable under UK Biobank disclosure controls.
5. Choose and add a license only after the authors agree on the license. No license is included by default.

## Publish to GitHub

1. Sign in to GitHub and create a new repository named `ukb-allergy-omics-analysis`.
2. Keep it private until the UK Biobank code review and author approval are complete.
3. Upload the contents of this repository, preserving the `scripts/` folder and the Markdown documentation.
4. Check the repository file list and `.gitignore`. Do not upload any real input data, generated score tables, GWAS results, logs, or intermediate files.
5. Review the commit history and repository settings. After the required checks, change visibility to public if approved.
6. Copy the canonical repository URL, for example `https://github.com/OWNER/ukb-allergy-omics-analysis`.

## Archive a release on Zenodo and obtain a DOI

1. Sign in to Zenodo using the same GitHub account.
2. Open Zenodo's GitHub integration settings and grant access to the approved repository.
3. In the GitHub repository, create a version tag such as `v1.0.0` and publish a GitHub Release for that tag.
4. Wait for Zenodo to archive the release. Open the resulting Zenodo record and review its metadata.
5. Add the final title, author list, publication date, version, description, keywords, related article information, and license consistent with the manuscript and author approvals.
6. Publish the Zenodo record. Zenodo will assign a version-specific DOI and a concept DOI for the software record.
7. Copy the DOI from the published Zenodo record. Use the version-specific DOI when citing the exact code version used for the manuscript; the concept DOI can point to the latest version.
8. Add the GitHub URL and DOI to the manuscript's Code Availability statement. Do not use a draft DOI or an unpublished placeholder.

## Suggested Code Availability wording

"The analysis code supporting the findings of this study is available at [GitHub repository URL] and archived at Zenodo ([version-specific DOI]). UK Biobank individual-level data are not publicly shared and are available to approved researchers subject to UK Biobank access procedures."

Replace the bracketed fields only after the repository and Zenodo record have been published. Confirm the final wording against the journal's current instructions and the applicable UK Biobank terms.
