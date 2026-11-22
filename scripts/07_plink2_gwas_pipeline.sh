#!/usr/bin/env bash
set -Eeuo pipefail

: "${OUTPUT_DIR:?}"
: "${PLINK1_PATH:?}"
: "${PLINK2_PATH:?}"
: "${INPUT_BED_PREFIX:?}"
: "${PHENO_FILE:?}"
: "${COVAR_FILE:?}"

THREADS="${THREADS:-4}"
MEMORY_MB="${MEMORY_MB:-300000}"
KING_CUTOFF="${KING_CUTOFF:-0.0884}"
PHENO_NAME="${PHENO_NAME:-any_allergy}"
GWAS_PREFIX="${GWAS_PREFIX:-ukb_gwas_results}"
GWAS_SEED="${GWAS_SEED:-9912}"
OUTPUT_DIR="$(mkdir -p "${OUTPUT_DIR}" && cd "${OUTPUT_DIR}" && pwd)"
GWAS_OUT="${OUTPUT_DIR}/${GWAS_PREFIX}"
GWAS_RESULT="${GWAS_OUT}.${PHENO_NAME}.glm.logistic.hybrid"
LOG_FILE="${OUTPUT_DIR}/gwas_run.log"

[[ -x "${PLINK1_PATH}" ]] || { printf '%s\n' "PLINK1_PATH is not executable." >&2; exit 1; }
[[ -x "${PLINK2_PATH}" ]] || { printf '%s\n' "PLINK2_PATH is not executable." >&2; exit 1; }
[[ -s "${PHENO_FILE}" ]] || { printf '%s\n' "PHENO_FILE is missing or empty." >&2; exit 1; }
[[ -s "${COVAR_FILE}" ]] || { printf '%s\n' "COVAR_FILE is missing or empty." >&2; exit 1; }

for chr in $(seq 1 22); do
    prefix="${INPUT_BED_PREFIX}${chr}_v3"
    for ext in bed bim fam; do
        [[ -s "${prefix}.${ext}" ]] || { printf '%s\n' "Missing genotype input for chromosome ${chr}: ${prefix}.${ext}" >&2; exit 1; }
    done
done

for chr in $(seq 1 22); do
    prefix="${OUTPUT_DIR}/cleaned_chr${chr}"
    if [[ ! -s "${prefix}.bed" || ! -s "${prefix}.bim" || ! -s "${prefix}.fam" ]]; then
        rm -f "${prefix}.bed" "${prefix}.bim" "${prefix}.fam" "${prefix}.log" "${prefix}.nosex"
        "${PLINK2_PATH}" \
            --bfile "${INPUT_BED_PREFIX}${chr}_v3" \
            --set-all-var-ids '@:#$r,$a' \
            --new-id-max-allele-len 178 missing \
            --rm-dup force-first \
            --make-bed \
            --memory "${MEMORY_MB}" \
            --threads "${THREADS}" \
            --out "${prefix}"
    fi
done

for chr in $(seq 1 22); do
    prefix="${OUTPUT_DIR}/clean_prune_chr${chr}"
    if [[ ! -s "${prefix}.prune.in" ]]; then
        "${PLINK1_PATH}" \
            --bfile "${OUTPUT_DIR}/cleaned_chr${chr}" \
            --indep-pairwise 50 5 0.2 \
            --memory "${MEMORY_MB}" \
            --threads "${THREADS}" \
            --out "${prefix}"
    fi
    [[ -s "${prefix}.prune.in" ]] || { printf '%s\n' "LD-pruning output missing for chromosome ${chr}." >&2; exit 1; }
done

PRUNED_LIST="${OUTPUT_DIR}/all_pruned_snps.txt"
if [[ ! -s "${PRUNED_LIST}" ]]; then
    tmp_list="${PRUNED_LIST}.tmp"
    : > "${tmp_list}"
    for chr in $(seq 1 22); do
        cat "${OUTPUT_DIR}/clean_prune_chr${chr}.prune.in" >> "${tmp_list}"
    done
    LC_ALL=C sort -u "${tmp_list}" > "${PRUNED_LIST}"
    rm -f "${tmp_list}"
fi
[[ -s "${PRUNED_LIST}" ]] || { printf '%s\n' "The combined LD-pruned variant list is empty." >&2; exit 1; }

MERGE_LIST="${OUTPUT_DIR}/clean_chr_list.txt"
if [[ ! -s "${OUTPUT_DIR}/ukb_all_cleaned.pgen" || ! -s "${OUTPUT_DIR}/ukb_all_cleaned.pvar" || ! -s "${OUTPUT_DIR}/ukb_all_cleaned.psam" ]]; then
    : > "${MERGE_LIST}"
    for chr in $(seq 1 22); do
        printf '%s\n' "${OUTPUT_DIR}/cleaned_chr${chr}" >> "${MERGE_LIST}"
    done
    rm -f "${OUTPUT_DIR}/ukb_all_cleaned.pgen" "${OUTPUT_DIR}/ukb_all_cleaned.pvar" "${OUTPUT_DIR}/ukb_all_cleaned.psam"
    "${PLINK2_PATH}" \
        --pmerge-list "${MERGE_LIST}" bfile \
        --make-pgen \
        --memory "${MEMORY_MB}" \
        --threads "${THREADS}" \
        --out "${OUTPUT_DIR}/ukb_all_cleaned"
fi

if [[ ! -s "${OUTPUT_DIR}/ukb_all_pruned.pgen" || ! -s "${OUTPUT_DIR}/ukb_all_pruned.pvar" || ! -s "${OUTPUT_DIR}/ukb_all_pruned.psam" ]]; then
    "${PLINK2_PATH}" \
        --pfile "${OUTPUT_DIR}/ukb_all_cleaned" \
        --extract "${PRUNED_LIST}" \
        --make-pgen \
        --memory "${MEMORY_MB}" \
        --threads "${THREADS}" \
        --out "${OUTPUT_DIR}/ukb_all_pruned"
fi

if [[ ! -s "${OUTPUT_DIR}/ukb_unrelated.pgen" || ! -s "${OUTPUT_DIR}/ukb_unrelated.pvar" || ! -s "${OUTPUT_DIR}/ukb_unrelated.psam" ]]; then
    "${PLINK2_PATH}" \
        --pfile "${OUTPUT_DIR}/ukb_all_pruned" \
        --king-cutoff "${KING_CUTOFF}" \
        --make-pgen \
        --memory "${MEMORY_MB}" \
        --threads "${THREADS}" \
        --out "${OUTPUT_DIR}/ukb_unrelated"
fi

[[ -s "${OUTPUT_DIR}/ukb_unrelated.king.cutoff.in.id" ]] || { printf '%s\n' "KING retained-sample list was not produced." >&2; exit 1; }

if [[ ! -s "${OUTPUT_DIR}/unrelated_pca.eigenvec" ]]; then
    "${PLINK2_PATH}" \
        --pfile "${OUTPUT_DIR}/ukb_unrelated" \
        --pca approx 10 \
        --memory "${MEMORY_MB}" \
        --threads "${THREADS}" \
        --out "${OUTPUT_DIR}/unrelated_pca"
fi

[[ -s "${OUTPUT_DIR}/unrelated_pca.eigenvec" ]] || { printf '%s\n' "PCA eigenvector file is missing." >&2; exit 1; }

COVAR_FINAL="${OUTPUT_DIR}/covariates_final.txt"
if [[ ! -s "${COVAR_FINAL}" ]]; then
    COVAR_TMP="${COVAR_FINAL}.tmp"
    awk '
        BEGIN { FS = "[ \t]+"; OFS = "\t" }
        FNR == NR {
            if (FNR == 1) next
            key = $1 SUBSEP $2
            keep[key] = 1
            for (i = 1; i <= 10; i++) pc[key, i] = $(i + 2)
            next
        }
        FNR == 1 {
            print "FID", "IID", "sex", "age", "BMI", "TDI", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"
            next
        }
        NF >= 6 {
            key = $1 SUBSEP $2
            if (key in keep) {
                printf "%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5, $6
                for (i = 1; i <= 10; i++) printf "\t%s", pc[key, i]
                printf "\n"
                retained++
            } else {
                missing++
            }
        }
        END {
            if (retained == 0) exit 2
            if (missing > 0) printf "Covariate rows excluded because no unrelated-sample PCA was available: %d\n", missing > "/dev/stderr"
        }
    ' "${OUTPUT_DIR}/unrelated_pca.eigenvec" "${COVAR_FILE}" > "${COVAR_TMP}"
    mv "${COVAR_TMP}" "${COVAR_FINAL}"
fi

[[ -s "${COVAR_FINAL}" ]] || { printf '%s\n' "Final covariate file is missing or empty." >&2; exit 1; }

if [[ ! -s "${GWAS_RESULT}" ]]; then
    "${PLINK2_PATH}" \
        --pfile "${OUTPUT_DIR}/ukb_unrelated" \
        --pheno "${PHENO_FILE}" \
        --pheno-name "${PHENO_NAME}" \
        --covar "${COVAR_FINAL}" \
        --covar-name sex,age,BMI,TDI,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
        --maf 0.01 \
        --geno 0.05 \
        --hwe 1e-6 keep-fewhet \
        --vif 1000 \
        --glm firth-fallback hide-covar cols=chrom,pos,ref,alt,omitted,a1freq,nobs,beta,se,tz,p \
        --seed "${GWAS_SEED}" \
        --memory "${MEMORY_MB}" \
        --threads "${THREADS}" \
        --out "${GWAS_OUT}" > "${LOG_FILE}" 2>&1
fi

[[ -s "${GWAS_RESULT}" ]] || { printf '%s\n' "Expected PLINK2 GWAS result is missing: ${GWAS_RESULT}" >&2; exit 1; }

if [[ ! -s "${OUTPUT_DIR}/gwas_clumped.clumps" ]]; then
    "${PLINK2_PATH}" \
        --pfile "${OUTPUT_DIR}/ukb_unrelated" \
        --clump "${GWAS_RESULT}" \
        --clump-p1 5e-8 \
        --clump-p2 1e-5 \
        --clump-r2 0.1 \
        --clump-kb 250 \
        --clump-id-field ID \
        --clump-p-field P \
        --memory "${MEMORY_MB}" \
        --threads "${THREADS}" \
        --out "${OUTPUT_DIR}/gwas_clumped"
fi

[[ -s "${OUTPUT_DIR}/gwas_clumped.clumps" ]] || { printf '%s\n' "LD-clumping output is missing." >&2; exit 1; }
printf '%s\n' "GWAS and LD clumping completed."
