required_packages <- c("data.table", "R.utils", "TwoSampleMR", "MRPRESSO", "dplyr", "readxl", "ggplot2", "ggrepel", "scales", "psych")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(MRPRESSO)
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
})

get_env <- function(name, default = "") {
  value <- Sys.getenv(name, unset = default)
  if (!nzchar(value)) stop("Set environment variable: ", name, call. = FALSE)
  value
}

get_num_env <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (length(value) != 1L || !is.finite(value)) stop("Invalid numeric environment variable: ", name, call. = FALSE)
  value
}

normalize_chr <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y <- sub("^CHR", "", y)
  y[y %in% c("M", "MTDNA")] <- "MT"
  y
}

require_columns <- function(x, cols, label) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

compress_file <- function(path) {
  if (!file.exists(path)) stop("Cannot compress missing file: ", basename(path), call. = FALSE)
  R.utils::gzip(path, destname = paste0(path, ".gz"), overwrite = TRUE, remove = TRUE)
  paste0(path, ".gz")
}

GWAS_FILE <- get_env("GWAS_FILE")
DBSNP_FILE <- get_env("DBSNP_FILE")
PQTL_DIR <- get_env("PQTL_DIR")
TARGET_REFERENCE_FILE <- get_env("TARGET_REFERENCE_FILE")
UP_GENES_FILE <- get_env("UP_GENES_FILE")
DOWN_GENES_FILE <- get_env("DOWN_GENES_FILE")
OUTPUT_DIR <- get_env("OUTPUT_DIR")
GWAS_BUILD <- toupper(get_env("GWAS_BUILD"))
DBSNP_BUILD <- toupper(get_env("DBSNP_BUILD"))
PQTL_BUILD <- toupper(get_env("PQTL_BUILD"))

if (!identical(GWAS_BUILD, DBSNP_BUILD) || !identical(GWAS_BUILD, PQTL_BUILD)) {
  stop("GWAS, dbSNP, and pQTL genome builds must be identical.", call. = FALSE)
}
for (path in c(GWAS_FILE, DBSNP_FILE, TARGET_REFERENCE_FILE, UP_GENES_FILE, DOWN_GENES_FILE)) {
  if (!file.exists(path)) stop("Input file does not exist.", call. = FALSE)
}
if (!dir.exists(PQTL_DIR)) stop("PQTL_DIR does not exist.", call. = FALSE)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

DBSNP_CHUNK_SIZE <- as.integer(get_num_env("DBSNP_CHUNK_SIZE", 300000))
DBSNP_POSITION_OFFSET <- as.integer(get_num_env("DBSNP_POSITION_OFFSET", 1))
DBSNP_HAS_HEADER <- tolower(Sys.getenv("DBSNP_HAS_HEADER", unset = "false")) %in% c("true", "1", "yes")
MHC_START <- get_num_env("MHC_START", 25500000)
MHC_END <- get_num_env("MHC_END", 34000000)
MHC_CHROMOSOME <- normalize_chr(Sys.getenv("MHC_CHROMOSOME", unset = "6"))
MRPRESSO_DRAWS <- as.integer(get_num_env("MRPRESSO_DRAWS", 1000))
OUTCOME_NCASE <- suppressWarnings(as.numeric(Sys.getenv("OUTCOME_NCASE", unset = "")))
OUTCOME_NCONTROL <- suppressWarnings(as.numeric(Sys.getenv("OUTCOME_NCONTROL", unset = "")))
OUTCOME_PREVALENCE <- suppressWarnings(as.numeric(Sys.getenv("OUTCOME_PREVALENCE", unset = "")))
PQTL_N_COLUMN <- Sys.getenv("PQTL_N_COLUMN", unset = "N")

if (DBSNP_CHUNK_SIZE < 1L || MHC_START < 0 || MHC_END <= MHC_START || MRPRESSO_DRAWS < 100L) {
  stop("Invalid chunk, MHC, or MR-PRESSO settings.", call. = FALSE)
}

GWAS <- fread(GWAS_FILE, showProgress = TRUE)
if ("#CHROM" %in% names(GWAS)) setnames(GWAS, "#CHROM", "CHROM")
require_columns(GWAS, c("CHROM", "POS", "ID", "REF", "A1", "OMITTED", "A1_FREQ", "BETA", "SE", "P", "OBS_CT"), "GWAS summary statistics")
GWAS[, CHR_KEY := normalize_chr(CHROM)]
GWAS[, POS := suppressWarnings(as.numeric(POS))]
GWAS[, KEY := paste(CHR_KEY, POS, sep = ":")]
GWAS <- GWAS[!is.na(POS) & is.finite(POS) & POS > 0 & !is.na(CHR_KEY)]
TARGET_KEYS <- unique(GWAS$KEY)

DBSNP_CONNECTION <- gzfile(DBSNP_FILE, open = "rt")
if (DBSNP_HAS_HEADER) readLines(DBSNP_CONNECTION, n = 1L, warn = FALSE)
DBSNP_MATCHES <- data.table(key = character(), rsid = character())
repeat {
  LINES <- readLines(DBSNP_CONNECTION, n = DBSNP_CHUNK_SIZE, warn = FALSE)
  if (!length(LINES)) break
  CHUNK <- tryCatch(
    fread(text = paste(LINES, collapse = "\n"), sep = "\t", header = FALSE,
          select = c(2L, 3L, 5L), showProgress = FALSE),
    error = function(e) stop("Unable to parse dbSNP chunk: ", conditionMessage(e), call. = FALSE)
  )
  if (ncol(CHUNK) != 3L) stop("dbSNP columns 2, 3, and 5 could not be parsed.", call. = FALSE)
  setnames(CHUNK, c("chrom", "pos", "rsid"))
  CHUNK[, chrom := normalize_chr(chrom)]
  CHUNK[, pos := suppressWarnings(as.numeric(pos)) + DBSNP_POSITION_OFFSET]
  CHUNK[, key := paste(chrom, pos, sep = ":")]
  CHUNK <- CHUNK[!is.na(rsid) & nzchar(rsid) & key %chin% TARGET_KEYS, .(key, rsid = as.character(rsid))]
  if (nrow(CHUNK)) {
    DBSNP_MATCHES <- unique(rbindlist(list(DBSNP_MATCHES, CHUNK), use.names = TRUE))
  }
  rm(CHUNK, LINES)
  gc(FALSE)
}
close(DBSNP_CONNECTION)

if (nrow(DBSNP_MATCHES)) {
  DBSNP_MAP <- DBSNP_MATCHES[, .(rsid_count = uniqueN(rsid), rsid = if (uniqueN(rsid) == 1L) rsid[1L] else NA_character_), by = key]
  GWAS[DBSNP_MAP, rsid := i.rsid, on = c("KEY" = "key")]
} else {
  GWAS[, rsid := NA_character_]
}
GWAS[, final_SNP := fifelse(!is.na(rsid) & nzchar(rsid), rsid, as.character(ID))]
GWAS[, `:=`(
  A1 = toupper(as.character(A1)),
  A2 = toupper(as.character(OMITTED)),
  freq = suppressWarnings(as.numeric(A1_FREQ)),
  b = suppressWarnings(as.numeric(BETA)),
  se = suppressWarnings(as.numeric(SE)),
  p = suppressWarnings(as.numeric(P)),
  n = suppressWarnings(as.numeric(OBS_CT))
)]
GWAS_VALID <- GWAS[
  !is.na(final_SNP) & nzchar(final_SNP) &
    !is.na(A1) & nzchar(A1) & !is.na(A2) & nzchar(A2) & A1 != A2 &
    grepl("^[ACGT]+$", A1) & grepl("^[ACGT]+$", A2) &
    is.finite(freq) & freq > 0 & freq < 1 &
    is.finite(b) & is.finite(se) & se > 0 &
    is.finite(p) & p >= 0 & p <= 1 & is.finite(n) & n > 0
]
if (!nrow(GWAS_VALID)) stop("No valid GWAS variants remain after allele and statistic checks.", call. = FALSE)

SMR_MIXED <- GWAS_VALID[, .(SNP = final_SNP, A1, A2, freq, b, se, p, n)]
SMR_RSID <- GWAS_VALID[!is.na(rsid) & nzchar(rsid), .(SNP = rsid, A1, A2, freq, b, se, p, n)]
SMR_RSID <- SMR_RSID[order(p)][!duplicated(SNP)]
SMR_RSID_FULL <- GWAS_VALID[!is.na(rsid) & nzchar(rsid), .(SNP = rsid, chrom = CHROM, POS, A1, A2, freq, b, se, p, n)]
SMR_RSID_FULL <- SMR_RSID_FULL[order(p)][!duplicated(SNP)]
SMR_FULL <- GWAS_VALID[, .(SNP = final_SNP, chrom = CHROM, POS, A1, A2, freq, b, se, p, n)]

MIXED_FILE <- file.path(OUTPUT_DIR, "GWAS_mixed_SMR_format.ma")
RSID_FILE <- file.path(OUTPUT_DIR, "GWAS_rsid_only.txt")
FULL_FILE <- file.path(OUTPUT_DIR, "GWAS_full_variant_format.txt")
RSID_FULL_FILE <- file.path(OUTPUT_DIR, "GWAS_rsid_with_coordinates.txt")
fwrite(SMR_MIXED, MIXED_FILE, sep = "\t", quote = FALSE, na = "NA")
fwrite(SMR_RSID, RSID_FILE, sep = "\t", quote = FALSE, na = "NA")
fwrite(SMR_FULL, FULL_FILE, sep = "\t", quote = FALSE, na = "NA")
fwrite(SMR_RSID_FULL, RSID_FULL_FILE, sep = "\t", quote = FALSE, na = "NA")
OLD_WORKING_DIRECTORY <- getwd()
setwd(OUTPUT_DIR)
utils::tar(
  tarfile = "GWAS_mixed_SMR_format.tar.gz",
  files = basename(MIXED_FILE), compression = "gzip", tar = "internal"
)
setwd(OLD_WORKING_DIRECTORY)
RSID_GZ <- compress_file(RSID_FILE)
FULL_GZ <- compress_file(FULL_FILE)
RSID_FULL_GZ <- compress_file(RSID_FULL_FILE)

OUTCOME_RAW <- copy(SMR_RSID_FULL)
OUTCOME_RAW[, units := "log odds"]
if (is.finite(OUTCOME_NCASE) && OUTCOME_NCASE > 0 && is.finite(OUTCOME_NCONTROL) && OUTCOME_NCONTROL > 0) {
  OUTCOME_RAW[, `:=`(ncase = OUTCOME_NCASE, ncontrol = OUTCOME_NCONTROL)]
}
OUTCOME_ARGS <- list(
  dat = as.data.frame(OUTCOME_RAW), type = "outcome", snp_col = "SNP",
  beta_col = "b", se_col = "se", effect_allele_col = "A1",
  other_allele_col = "A2", pval_col = "p", eaf_col = "freq",
  samplesize_col = "n", units_col = "units", chr_col = "chrom",
  pos_col = "POS", min_pval = 1e-200
)
if ("ncase" %in% names(OUTCOME_RAW)) OUTCOME_ARGS$ncase_col <- "ncase"
if ("ncontrol" %in% names(OUTCOME_RAW)) OUTCOME_ARGS$ncontrol_col <- "ncontrol"
OUTCOME_DATA <- do.call(TwoSampleMR::format_data, OUTCOME_ARGS)
if (is.finite(OUTCOME_PREVALENCE) && OUTCOME_PREVALENCE > 0 && OUTCOME_PREVALENCE < 1) {
  OUTCOME_DATA$prevalence.outcome <- OUTCOME_PREVALENCE
}
OUTCOME_DATA <- as.data.table(OUTCOME_DATA)
setorder(OUTCOME_DATA, pval.outcome)
OUTCOME_DATA <- OUTCOME_DATA[!duplicated(SNP)]

PQTL_FILES <- list.files(PQTL_DIR, pattern = "-cis_unclump\\.csv$", full.names = TRUE, recursive = FALSE)
if (!length(PQTL_FILES)) stop("No pQTL input files matched the required filename pattern.", call. = FALSE)
MR_RESULT_LIST <- list()
FAILURE_LIST <- list()
MR_INDEX <- 0L
FAILURE_INDEX <- 0L

for (PQTL_FILE in PQTL_FILES) {
  PROTEIN_NAME <- sub("-cis_unclump\\.csv$", "", basename(PQTL_FILE))
  tryCatch({
    PQTL <- fread(PQTL_FILE, showProgress = FALSE)
    require_columns(PQTL, c("rsids", "Beta", "SE", "effectAllele", "otherAllele", "Pval", "Chrom", "Pos", "ImpMAF"), PROTEIN_NAME)
    PQTL[, `:=`(
      Chrom_key = normalize_chr(Chrom),
      Pos_num = suppressWarnings(as.numeric(Pos)),
      Pval_num = suppressWarnings(as.numeric(Pval)),
      Beta_num = suppressWarnings(as.numeric(Beta)),
      SE_num = suppressWarnings(as.numeric(SE)),
      ImpMAF_num = suppressWarnings(as.numeric(ImpMAF))
    )]
    PQTL <- PQTL[!(Chrom_key == MHC_CHROMOSOME & Pos_num >= MHC_START & Pos_num <= MHC_END)]
    PQTL <- PQTL[
      !is.na(rsids) & nzchar(rsids) & is.finite(Pos_num) &
        is.finite(Pval_num) & Pval_num >= 0 & Pval_num <= 1 &
        is.finite(Beta_num) & is.finite(SE_num) & SE_num > 0 &
        is.finite(ImpMAF_num) & ImpMAF_num > 0 & ImpMAF_num < 1
    ]
    if (!nrow(PQTL)) stop("No pQTL variants remain after QC and MHC exclusion.")
    setorder(PQTL, Pval_num)
    PQTL <- PQTL[!duplicated(rsids)]
    setnames(PQTL, c("Beta_num", "SE_num", "Pval_num", "ImpMAF_num", "Chrom_key", "Pos_num"),
             c("Beta_use", "SE_use", "Pval_use", "MAF_use", "Chrom_use", "Pos_use"))
    EXPOSURE_ARGS <- list(
      dat = as.data.frame(PQTL), type = "exposure", snp_col = "rsids",
      beta_col = "Beta_use", se_col = "SE_use", effect_allele_col = "effectAllele",
      other_allele_col = "otherAllele", pval_col = "Pval_use", eaf_col = "MAF_use",
      chr_col = "Chrom_use", pos_col = "Pos_use", min_pval = 1e-200
    )
    HAS_PQTL_N <- PQTL_N_COLUMN %in% names(PQTL)
    if (HAS_PQTL_N) EXPOSURE_ARGS$samplesize_col <- PQTL_N_COLUMN
    EXPOSURE_DATA <- do.call(TwoSampleMR::format_data, EXPOSURE_ARGS)
    EXPOSURE_DATA <- as.data.table(EXPOSURE_DATA)
    setorder(EXPOSURE_DATA, pval.exposure)
    EXPOSURE_DATA <- EXPOSURE_DATA[!duplicated(SNP)]
    if (!nrow(EXPOSURE_DATA)) stop("Exposure formatting returned no variants.")

    HARMONIZED <- tryCatch(
      TwoSampleMR::harmonise_data(as.data.frame(EXPOSURE_DATA), as.data.frame(OUTCOME_DATA), action = 2),
      error = function(e) stop("Harmonisation failed: ", conditionMessage(e))
    )
    HARMONIZED <- as.data.table(HARMONIZED)
    if (!nrow(HARMONIZED)) stop("No variants remained after harmonisation.")
    if ("mr_keep" %in% names(HARMONIZED)) HARMONIZED <- HARMONIZED[mr_keep %in% TRUE]
    HARMONIZED[, F_statistic := (beta.exposure / se.exposure)^2]
    HARMONIZED <- HARMONIZED[is.finite(F_statistic) & F_statistic > 10]
    if (!nrow(HARMONIZED)) stop("No instruments passed F-statistic > 10.")
    F_MEAN <- mean(HARMONIZED$F_statistic)
    F_MEDIAN <- median(HARMONIZED$F_statistic)

    MR_RESULT <- TwoSampleMR::mr(as.data.frame(HARMONIZED))
    if (!nrow(MR_RESULT)) stop("MR returned no estimates.")

    HET_Q <- HET_DF <- HET_P <- EGGER_INTERCEPT <- EGGER_SE <- EGGER_P <- PRESSO_GLOBAL_P <- NA_real_
    STEIGER_DIRECTION <- STEIGER_P <- NA_real_
    if (nrow(HARMONIZED) >= 3L) {
      HET <- tryCatch(TwoSampleMR::mr_heterogeneity(as.data.frame(HARMONIZED)), error = function(e) NULL)
      if (!is.null(HET) && nrow(HET)) {
        HET_ROW <- HET[HET$method == "Inverse variance weighted", , drop = FALSE]
        if (!nrow(HET_ROW)) HET_ROW <- HET[1L, , drop = FALSE]
        HET_Q <- HET_ROW$Q[1L]
        HET_DF <- HET_ROW$Q_df[1L]
        HET_P <- HET_ROW$Q_pval[1L]
      }
      EGGER <- tryCatch(TwoSampleMR::mr_pleiotropy_test(as.data.frame(HARMONIZED)), error = function(e) NULL)
      if (!is.null(EGGER) && nrow(EGGER)) {
        EGGER_INTERCEPT <- EGGER$egger_intercept[1L]
        EGGER_SE <- EGGER$se[1L]
        EGGER_P <- EGGER$pval[1L]
      }
    }
    if (nrow(HARMONIZED) >= 4L) {
      PRESSO <- tryCatch(
        MRPRESSO::mr_presso(
          BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome", SdExposure = "se.exposure",
          data = as.data.frame(HARMONIZED), OUTLIERtest = TRUE,
          DISTORTIONtest = TRUE, NbDistribution = MRPRESSO_DRAWS,
          SignifThreshold = 0.05
        ), error = function(e) NULL
      )
      if (!is.null(PRESSO)) {
        PRESSO_GLOBAL_P <- tryCatch(
          as.numeric(PRESSO[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]),
          error = function(e) NA_real_
        )
      }
    }
    STEIGER_AVAILABLE <- HAS_PQTL_N &&
      is.finite(OUTCOME_NCASE) && OUTCOME_NCASE > 0 &&
      is.finite(OUTCOME_NCONTROL) && OUTCOME_NCONTROL > 0 &&
      is.finite(OUTCOME_PREVALENCE) && OUTCOME_PREVALENCE > 0 && OUTCOME_PREVALENCE < 1
    if (STEIGER_AVAILABLE) {
      HARMONIZED$prevalence.outcome <- OUTCOME_PREVALENCE
      STEIGER <- tryCatch(TwoSampleMR::steiger_filtering(as.data.frame(HARMONIZED)), error = function(e) NULL)
      if (!is.null(STEIGER) && nrow(STEIGER)) {
        STEIGER_DIRECTION <- mean(as.logical(STEIGER$steiger_dir), na.rm = TRUE)
        STEIGER_P <- median(STEIGER$steiger_pval, na.rm = TRUE)
      }
    }

    MR_RESULT <- as.data.table(MR_RESULT)
    MR_RESULT[, `:=`(
      Protein = PROTEIN_NAME,
      Instrument_N = nrow(HARMONIZED),
      Mean_F_statistic = F_MEAN,
      Median_F_statistic = F_MEDIAN,
      Heterogeneity_Q = HET_Q,
      Heterogeneity_df = HET_DF,
      Heterogeneity_P = HET_P,
      Egger_intercept = EGGER_INTERCEPT,
      Egger_intercept_SE = EGGER_SE,
      Egger_intercept_P = EGGER_P,
      MRPRESSO_global_P = PRESSO_GLOBAL_P,
      Steiger_direction_proportion = STEIGER_DIRECTION,
      Steiger_median_P = STEIGER_P
    )]
    MR_INDEX <- MR_INDEX + 1L
    MR_RESULT_LIST[[MR_INDEX]] <- MR_RESULT
  }, error = function(e) {
    FAILURE_INDEX <<- FAILURE_INDEX + 1L
    FAILURE_LIST[[FAILURE_INDEX]] <<- data.table(
      Protein = PROTEIN_NAME, Stage = "MR", Reason = conditionMessage(e)
    )
  })
}

if (!length(MR_RESULT_LIST)) stop("No pQTL protein produced a valid MR result.", call. = FALSE)
MR_RESULTS <- rbindlist(MR_RESULT_LIST, fill = TRUE)
setorder(MR_RESULTS, Protein, method, pval)
fwrite(MR_RESULTS, file.path(OUTPUT_DIR, "Batch_MR_Results.csv"))
if (length(FAILURE_LIST)) {
  MR_FAILURES <- rbindlist(FAILURE_LIST, fill = TRUE)
} else {
  MR_FAILURES <- data.table(Protein = character(), Stage = character(), Reason = character())
}
fwrite(MR_FAILURES, file.path(OUTPUT_DIR, "MR_Failures.csv"))

MR_RESULTS[, method_priority := fcase(
  method == "Inverse variance weighted", 1L,
  method == "Wald ratio", 2L,
  method == "MR Egger", 3L,
  default = 4L
)]
PRIMARY_RESULTS <- MR_RESULTS[order(Protein, method_priority, pval)][!duplicated(Protein)]
PRIMARY_RESULTS[, FDR := p.adjust(pval, method = "BH")]
CANDIDATES <- PRIMARY_RESULTS[pval < 0.05]
CANDIDATES[, `:=`(
  OR = exp(b),
  OR_95_CI_Lower = exp(b - 1.96 * se),
  OR_95_CI_Upper = exp(b + 1.96 * se)
)]
CANDIDATE_OUTPUT <- CANDIDATES[, .(
  Protein, Method = method, SNP_Count = nsnp, Beta = b, SE = se,
  P_value = pval, FDR, OR, OR_95_CI_Lower, OR_95_CI_Upper
)]
fwrite(CANDIDATE_OUTPUT, file.path(OUTPUT_DIR, "Candidate_Proteins.csv"))

STANDARD <- as.data.table(readxl::read_excel(TARGET_REFERENCE_FILE, sheet = 1))
require_columns(STANDARD, c("hgnc_names", "ensembl_gene_id", "druggability_tier", "description", "chr_b37"), "Target reference table")
STANDARD <- STANDARD[!is.na(hgnc_names) & nzchar(trimws(as.character(hgnc_names)))]
STANDARD[, hgnc_key := toupper(trimws(as.character(hgnc_names)))]
if (anyDuplicated(STANDARD$hgnc_key)) stop("Target reference contains duplicated HGNC symbols.", call. = FALSE)
REF_SYMBOLS <- STANDARD$hgnc_key

extract_gene_token <- function(protein, index) {
  parts <- strsplit(as.character(protein), "_", fixed = TRUE)[[1L]]
  if (length(parts) >= index && nzchar(parts[index])) toupper(parts[index]) else NA_character_
}

lookup_gene <- function(gene) {
  if (is.na(gene) || !nzchar(gene)) return(NA_character_)
  idx <- match(toupper(gene), REF_SYMBOLS)
  if (is.na(idx)) NA_character_ else as.character(STANDARD$hgnc_names[idx])
}

if (nrow(CANDIDATES)) {
  CANDIDATE_ANNOTATION <- copy(CANDIDATES)
  CANDIDATE_ANNOTATION[, gene1 := vapply(Protein, extract_gene_token, character(1), index = 3L)]
  CANDIDATE_ANNOTATION[, gene2 := vapply(Protein, extract_gene_token, character(1), index = 4L)]
  CANDIDATE_ANNOTATION[, matched_gene1 := vapply(gene1, lookup_gene, character(1))]
  CANDIDATE_ANNOTATION[, matched_gene2 := vapply(gene2, lookup_gene, character(1))]
  CANDIDATE_ANNOTATION[, matched_hgnc := fifelse(!is.na(matched_gene1), matched_gene1, matched_gene2)]
  ANNOTATED <- merge(CANDIDATE_ANNOTATION, STANDARD, by.x = "matched_hgnc", by.y = "hgnc_names", all.x = TRUE, sort = FALSE)
  ANNOTATED[, tier_order := match(druggability_tier, c("Tier 1", "Tier 2", "Tier 3A", "Tier 3B"))]
  setorder(ANNOTATED, tier_order, Protein, na.last = TRUE)
  TARGET_OUTPUT <- ANNOTATED[, .(
    Protein, Method = method, SNP_Count = nsnp, Beta = b, SE = se,
    P_value = pval, FDR, OR = exp(b), OR_95_CI_Lower = exp(b - 1.96 * se),
    OR_95_CI_Upper = exp(b + 1.96 * se), Extracted_Gene1 = gene1,
    Extracted_Gene2 = gene2, Matched_HGNC = matched_hgnc,
    Ensembl_ID = ensembl_gene_id, Druggability_Tier = druggability_tier,
    Description = description, Chromosome = chr_b37
  )]
} else {
  TARGET_OUTPUT <- data.table()
}
fwrite(TARGET_OUTPUT, file.path(OUTPUT_DIR, "Potential_Druggable_Targets.csv"))

PLOT_DATA <- copy(MR_RESULTS)
PLOT_DATA[, method_priority := fcase(
  method == "Inverse variance weighted", 1L,
  method == "Wald ratio", 2L,
  method == "MR Egger", 3L,
  default = 4L
)]
PLOT_DATA <- PLOT_DATA[order(Protein, method_priority, pval)][!duplicated(Protein)]
if (nrow(PLOT_DATA)) {
  GENE_MAP <- if (nrow(TARGET_OUTPUT)) unique(TARGET_OUTPUT[, .(Protein, Extracted_Gene1)]) else data.table(Protein = character(), Extracted_Gene1 = character())
  PLOT_DATA <- merge(PLOT_DATA, GENE_MAP, by = "Protein", all.x = TRUE, sort = FALSE)
  PLOT_DATA[, Label := fifelse(!is.na(Extracted_Gene1) & nzchar(Extracted_Gene1), Extracted_Gene1, Protein)]
  setorder(PLOT_DATA, method_priority, pval)
  PLOT_DATA <- PLOT_DATA[!duplicated(Label)]
  PLOT_DATA[, P_value := pmax(as.numeric(pval), .Machine$double.xmin)]
  PLOT_DATA[, Direction := fcase(
    b > 0 & P_value < 0.05, "Positive (Risk)",
    b < 0 & P_value < 0.05, "Negative (Protective)",
    default = "Not significant"
  )]
  TOP_GENES <- PLOT_DATA[Direction != "Not significant"][order(P_value), .SD[seq_len(min(.N, 10L))], by = Direction]
  MORANDI_COLORS <- c("Not significant" = "#C8C8C8", "Positive (Risk)" = "#C15B5B", "Negative (Protective)" = "#76A08A")
  VOLCANO <- ggplot(PLOT_DATA, aes(x = b, y = -log10(P_value))) +
    geom_point(data = PLOT_DATA[Direction == "Not significant"], aes(color = Direction), alpha = 0.4, size = 1.5) +
    geom_point(data = PLOT_DATA[Direction != "Not significant"], aes(color = Direction), alpha = 0.9, size = 3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_label_repel(data = TOP_GENES, aes(label = Label), size = 3.5,
                     fill = scales::alpha("white", 0.85), color = "black",
                     label.size = 0.3, box.padding = 0.5, point.padding = 0.3,
                     max.overlaps = 50, segment.color = "grey50") +
    scale_color_manual(values = MORANDI_COLORS) +
    theme_bw(base_size = 14) +
    labs(title = "Volcano Plot of MR Results", x = "Beta (log odds)",
         y = "-log10(P-value)", color = "Causal Effect") +
    theme(legend.position = "top", plot.title = element_text(hjust = 0.5))
  ggsave(file.path(OUTPUT_DIR, "Volcano_Plot_MR_Top20.pdf"), plot = VOLCANO, width = 10, height = 7)
}

UP_GENES <- fread(UP_GENES_FILE)
DOWN_GENES <- fread(DOWN_GENES_FILE)
require_columns(UP_GENES, "gene", "upregulated gene list")
require_columns(DOWN_GENES, "gene", "downregulated gene list")
CANDIDATE_GENES <- if (nrow(TARGET_OUTPUT)) unique(na.omit(toupper(TARGET_OUTPUT$Extracted_Gene1))) else character()
UP_SET <- unique(na.omit(toupper(as.character(UP_GENES$gene))))
DOWN_SET <- unique(na.omit(toupper(as.character(DOWN_GENES$gene))))
INTERSECTIONS <- data.table(
  Direction = c(rep("Upregulated", length(intersect(CANDIDATE_GENES, UP_SET))),
                rep("Downregulated", length(intersect(CANDIDATE_GENES, DOWN_SET)))),
  Gene = c(intersect(CANDIDATE_GENES, UP_SET), intersect(CANDIDATE_GENES, DOWN_SET))
)
fwrite(INTERSECTIONS, file.path(OUTPUT_DIR, "Candidate_Gene_Expression_Intersections.csv"))

SUMMARY <- data.table(
  Metric = c("GWAS variants retained", "GWAS variants with unambiguous rsID", "pQTL files found",
             "Proteins with MR results", "Nominally significant primary MR candidates",
             "Druggable target annotations matched", "Upregulated gene intersections",
             "Downregulated gene intersections"),
  Value = c(nrow(GWAS_VALID), sum(!is.na(GWAS_VALID$rsid)), length(PQTL_FILES),
            uniqueN(MR_RESULTS$Protein), nrow(CANDIDATES),
            if (nrow(TARGET_OUTPUT)) sum(!is.na(TARGET_OUTPUT$Matched_HGNC)) else 0L,
            length(intersect(CANDIDATE_GENES, UP_SET)),
            length(intersect(CANDIDATE_GENES, DOWN_SET)) )
)
fwrite(SUMMARY, file.path(OUTPUT_DIR, "Analysis_Summary.csv"))
