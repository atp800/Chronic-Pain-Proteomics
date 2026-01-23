# ======================================================
# DIA Replicate QC — relaxed, inclusion-first (Stage-aware)
# Assess → (optional) Normalize → (optional) Auto-Trim → Export
# Stage folders: "Raw QC" (always), "After norm" (if used), "After trim" (if used), "Final Outputs" (always)
# PCA logic unchanged; final-only PDF; per-protein status across stages.
# ======================================================

# -------------------- User Parameters --------------------
file_path       <- "Chronic Code/Data_To_Clean.xlsx"
sheet_name      <- "Soluble - Time"
id_cols         <- c("Sample_ID", "Replicate")      # MUST exist in sheet
ignore_cols     <- c("Condition")                   # optional metadata; safe if missing

# If your sheet already contains log2 intensities (e.g., Spectronaut PG.Log2Quantity), set TRUE
already_log2    <- FALSE

# Toggle median normalization (per-sample, per-replicate) in log2 space
# TRUE = force; FALSE = skip; "auto" = apply only if QC improves (see heuristic below)
DO_NORMALIZE    <- "auto"

# Trimming strategy:
# FALSE = no trimming; TRUE = force trimming; "auto" = apply if heuristics indicate benefit
DO_EXTREME_TRIM <- "auto"

# Extreme-value rules (used for trimming when applied)
extreme_iqr_k             <- 3.0    # extreme if outside [Q1-3*IQR, Q3+3*IQR]
protein_extreme_frac_drop <- 0.10   # drop protein if ≥10% of samples are extreme

# --- RELAXED thresholds (flagging/QC bands) ---
coverage_min     <- 0.125   # keep if protein quantified in ≥50% of samples (matrix/PCA only)
icc_reliable_min <- 0.50
icc_inspect_min  <- 0.30
cv_reliable_max  <- 35
cv_inspect_max   <- 60
ba_abs_bias_max  <- 0.60   # |mean(Rep1-Rep2)| in log2 (~1.52x)
alpha_mcnemar    <- 0.001  # missingness asymmetry

# Output folders & naming
base_dir   <- "Chronic Code/Clean Output/x"
date_str   <- format(Sys.Date(), "%Y-%m-%d")
root_dir   <- file.path(base_dir, date_str, sheet_name)
if (!dir.exists(root_dir)) dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)

# Stage subfolders (created only when relevant)
raw_qc_dir     <- file.path(root_dir, "Raw QC");        dir.create(raw_qc_dir, showWarnings = FALSE)
after_norm_dir <- file.path(root_dir, "After norm")     # created only if DO_NORMALIZE is applied
after_trim_dir <- file.path(root_dir, "After trim")     # created only if trim actually applied
final_dir      <- file.path(root_dir, "Final Outputs"); dir.create(final_dir, showWarnings = FALSE)

DATA_STUB  <- tools::file_path_sans_ext(basename(file_path))
NAME_BASE  <- paste0(DATA_STUB, "_", sheet_name)  # e.g., "20251009_CleanUp_Acute-solubles"

# Writing options
WRITE_PDF          <- TRUE     # final-only combined PDF
ALLOW_INSTALL_PKGS <- FALSE

# Reproducibility
set.seed(20250827L)

# -------------------- Libraries -------------------------
need_pkgs <- c("readxl","dplyr","tidyr","readr","ggplot2","irr",
               "patchwork","openxlsx","rlang","RColorBrewer","pheatmap")
for (p in need_pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  if (!ok && ALLOW_INSTALL_PKGS) install.packages(p)
}
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(irr);    library(patchwork); library(openxlsx); library(rlang)
  library(RColorBrewer); library(pheatmap)
})

# -------------------- Resolve function conflicts --------------------
conflicted::conflicts_prefer(dplyr::intersect)
conflicted::conflicts_prefer(dplyr::setdiff)
conflicted::conflicts_prefer(dplyr::union)
conflicted::conflicts_prefer(dplyr::filter)


# -------------------- Helpers ---------------------------
# zeros/non-positive → NA (only used when data are linear), then safe log2
safe_log2 <- function(x) { x[x <= 0] <- NA_real_; log2(x) }
median_normalize <- function(v) v - stats::median(v, na.rm = TRUE)
NA_STRING <- "NA"
stage_stub <- function(stage_dir, stem) file.path(stage_dir, paste0(NAME_BASE, "_", stem))

# keep numeric columns with ≥3 finite values and non-zero variance
drop_zero_var <- function(X) {
  X <- as.data.frame(X)
  num <- vapply(X, is.numeric, logical(1))
  X <- X[, num, drop = FALSE]
  keep <- vapply(X, function(v) sum(is.finite(v)) >= 3 && stats::sd(v, na.rm = TRUE) > 0, logical(1))
  X[, keep, drop = FALSE]
}

# QC per protein (relaxed bands, flag only; robust ICC)
qc_per_protein_relaxed <- function(w_df) {
  w_df %>%
    dplyr::group_by(Protein) %>%
    dplyr::summarise(Rep1_vec = list(Rep1), Rep2_vec = list(Rep2), .groups = "drop") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      v1 = list(as.numeric(unlist(Rep1_vec))),
      v2 = list(as.numeric(unlist(Rep2_vec))),
      mask_pair = list(stats::complete.cases(v1, v2)),
      v1p = list(v1[unlist(mask_pair)]),
      v2p = list(v2[unlist(mask_pair)]),
      n_pairs   = length(v1p),
      ICC_value = {
        if (n_pairs >= 2) {
          obj <- try(suppressWarnings(
            irr::icc(cbind(unlist(v1p), unlist(v2p)),
                     model="twoway", type="agreement", unit="single")
          ), silent = TRUE)
          if (!inherits(obj, "try-error") && !is.null(obj$value)) obj$value else NA_real_
        } else NA_real_
      },
      CV_med_pct = if (n_pairs >= 1) {
        a_lin <- 2^unlist(v1p); b_lin <- 2^unlist(v2p)
        m <- (a_lin + b_lin)/2
        s <- sqrt(((a_lin - m)^2 + (b_lin - m)^2))
        stats::median(100 * (s / m), na.rm = TRUE)
      } else NA_real_,
      Diff_log2_mean = if (n_pairs >= 1) mean(unlist(v1p) - unlist(v2p), na.rm = TRUE) else NA_real_,
      Diff_log2_sd   = if (n_pairs >= 2) stats::sd(unlist(v1p) - unlist(v2p), na.rm = TRUE) else NA_real_,
      BA_LOA_lower   = if (!is.na(Diff_log2_sd)) Diff_log2_mean - 1.96 * Diff_log2_sd else NA_real_,
      BA_LOA_upper   = if (!is.na(Diff_log2_sd)) Diff_log2_mean + 1.96 * Diff_log2_sd else NA_real_,
      Miss_b = sum(!is.na(v1) &  is.na(v2)),
      Miss_c = sum( is.na(v1) & !is.na(v2)),
      Miss_p = { disc <- Miss_b + Miss_c; if (disc > 0) min(1, 2 * stats::pbinom(min(Miss_b, Miss_c), disc, 0.5)) else NA_real_ }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      BAD_BIAS = !is.na(Diff_log2_mean) & abs(Diff_log2_mean) > ba_abs_bias_max,
      BAD_MISS = !is.na(Miss_p) & (Miss_p < alpha_mcnemar),
      ICC_band = dplyr::case_when(
        is.na(ICC_value)              ~ "drop",
        ICC_value >= icc_reliable_min ~ "reliable",
        ICC_value >= icc_inspect_min  ~ "inspect",
        TRUE                          ~ "drop"
      ),
      CV_band  = dplyr::case_when(
        is.na(CV_med_pct)             ~ "drop",
        CV_med_pct <= cv_reliable_max ~ "reliable",
        CV_med_pct <= cv_inspect_max  ~ "inspect",
        TRUE                          ~ "drop"
      ),
      QC_Category_relaxed = dplyr::case_when(
        BAD_BIAS | BAD_MISS                             ~ "Drop",
        ICC_band == "drop" | CV_band == "drop"         ~ "Drop",
        ICC_band == "reliable" & CV_band == "reliable" ~ "Reliable",
        TRUE                                            ~ "Inspect"
      )
    )
}

# Core QC + plots (returns ggplots and summaries)
run_qc <- function(w_df, label) {
  pair_tbl <- w_df %>%
    dplyr::mutate(pair_ok = stats::complete.cases(Rep1, Rep2)) %>%
    dplyr::group_by(Sample_ID) %>%
    dplyr::summarise(n_pairs = sum(pair_ok), .groups = "drop")
  
  per_sample_rho <- w_df %>%
    dplyr::group_by(Sample_ID) %>%
    dplyr::summarise(
      Spearman_rho = {
        ok <- stats::complete.cases(Rep1, Rep2)
        if (sum(ok) >= 2) suppressWarnings(stats::cor(Rep1[ok], Rep2[ok], method = "spearman")) else NA_real_
      },
      .groups = "drop"
    )
  
  all_pairs <- w_df %>%
    dplyr::transmute(Mean_log2 = rowMeans(cbind(Rep1, Rep2), na.rm = TRUE),
                     Diff_log2 = Rep1 - Rep2) %>%
    dplyr::filter(is.finite(Mean_log2), is.finite(Diff_log2))
  ba_bias <- mean(all_pairs$Diff_log2, na.rm = TRUE)
  ba_sd   <- stats::sd(all_pairs$Diff_log2, na.rm = TRUE)
  ba_lwr  <- ba_bias - 1.96 * ba_sd
  ba_upr  <- ba_bias + 1.96 * ba_sd
  ba_text <- sprintf("Bias = %.3f | LoA [%.3f, %.3f] | SD = %.3f", ba_bias, ba_lwr, ba_upr, ba_sd)
  
  p_spearman <- per_sample_rho %>%
    dplyr::filter(is.finite(Spearman_rho)) %>%
    ggplot2::ggplot(ggplot2::aes(Spearman_rho)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::geom_vline(xintercept = c(0.8, 0.9), linetype = "dashed") +
    ggplot2::labs(title = paste0("Per-sample Spearman rho (", label, ")"),
                  x = "Spearman rho", y = "count")
  
  p_ba <- ggplot2::ggplot(all_pairs, ggplot2::aes(Mean_log2, Diff_log2)) +
    ggplot2::geom_point(alpha = 0.3, size = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::geom_hline(yintercept = ba_bias) +
    ggplot2::geom_hline(yintercept = c(ba_lwr, ba_upr), linetype = "dashed") +
    ggplot2::labs(title = paste0("Global Bland–Altman (", label, ")"),
                  subtitle = ba_text, x = "Mean (log2)", y = "Rep1 - Rep2 (log2)")
  
  pool_long <- w_df %>%
    dplyr::select(Rep1, Rep2) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "Rep", values_to = "Log2") %>%
    dplyr::filter(is.finite(Log2))
  p_dens <- ggplot2::ggplot(pool_long, ggplot2::aes(Log2, colour = Rep)) +
    ggplot2::geom_density() +
    ggplot2::labs(title = paste0("Pooled intensity distributions (", label, ")"), x = "log2 intensity")
  
  list(
    summary = tibble::tibble(
      Label = label,
      Global_Spearman_rho = suppressWarnings(stats::cor(w_df$Rep1, w_df$Rep2,
                                                        use = "pairwise.complete.obs",
                                                        method = "spearman")),
      Global_BA_bias = ba_bias,
      Global_BA_LOA_lower = ba_lwr,
      Global_BA_LOA_upper = ba_upr
    ),
    plots = list(spearman = p_spearman, bland_altman = p_ba, density = p_dens),
    per_sample_rho = per_sample_rho,
    per_protein = qc_per_protein_relaxed(w_df),
    pair_check = pair_tbl,
    all_pairs = all_pairs
  )
}

# Wrapper to run QC and save 3 core plots to a given folder
qc_stage <- function(w_df, label, out_dir) {
  qc <- run_qc(w_df, label)
  ggplot2::ggsave(paste0(stage_stub(out_dir, label), "_Spearman_Hist.png"), qc$plots$spearman,     width = 6.5, height = 4.5, dpi = 300)
  ggplot2::ggsave(paste0(stage_stub(out_dir, label), "_BlandAltman.png"),   qc$plots$bland_altman, width = 7.5, height = 4.5, dpi = 300)
  ggplot2::ggsave(paste0(stage_stub(out_dir, label), "_Density.png"),       qc$plots$density,      width = 6.5, height = 4.5, dpi = 300)
  qc
}

# ---------- PCA helpers (UNCHANGED LOGIC) ----------
make_pca_matrix <- function(long_df, id_cols_no_rep = c("Sample_ID"), cov_min = coverage_min) {
  wide <- long_df %>%
    dplyr::select(dplyr::all_of(c(id_cols_no_rep, "Protein", "Combined"))) %>%
    tidyr::pivot_wider(names_from = Protein, values_from = Combined)
  meta <- wide[, id_cols_no_rep, drop = FALSE]
  mat  <- as.data.frame(wide[, setdiff(names(wide), id_cols_no_rep), drop = FALSE])
  
  cov <- colMeans(is.finite(as.matrix(mat)), na.rm = TRUE)
  keep <- names(cov)[cov >= cov_min]
  if (length(keep) < 2L) {
    ord <- order(cov, decreasing = TRUE, na.last = NA)
    keep <- names(cov)[ord[seq_len(min(200, length(ord)))]]
  }
  mat <- as.data.frame(mat[, keep, drop = FALSE])
  
  for (j in seq_along(mat)) {
    cj <- mat[[j]]
    med <- stats::median(cj, na.rm = TRUE)
    cj[!is.finite(cj)] <- med
    mat[[j]] <- cj
  }
  mat <- drop_zero_var(mat)
  list(X = as.matrix(mat), meta = meta, kept_proteins = keep, coverage = cov[keep])
}

run_pca <- function(X) {
  if (!is.matrix(X) || ncol(X) < 2L) stop("PCA: need at least 2 variables.")
  prcomp(X, center = TRUE, scale. = TRUE)
}
var_explained <- function(pca) (pca$sdev^2) / sum(pca$sdev^2)

# -------------------- Load & Prepare ---------------------

# 1. Read the sheet
raw0 <- readxl::read_excel(path = file_path, sheet = sheet_name)

# 2. Identify id/meta/protein columns
ignore_cols <- intersect(ignore_cols, names(raw0))
missing_ids <- setdiff(id_cols, names(raw0))
if (length(missing_ids))
  stop("Missing id_cols: ", paste(missing_ids, collapse = ", "))

meta_cols    <- unique(c(id_cols, ignore_cols))
protein_cols <- setdiff(names(raw0), meta_cols)
if (!length(protein_cols))
  stop("No protein columns detected.")

# 3. Coerce protein columns to numeric
raw <- raw0 %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(protein_cols),
      ~ suppressWarnings(as.numeric(as.character(.)))
    )
  )

# 4. Handle explicit 0 and NaN as missing values (count for audit)
zeros_to_NA_count <- NA_integer_
nan_to_NA_count   <- NA_integer_

before_zero <- sum(unlist(lapply(raw[protein_cols],
                                 function(v) sum(v == 0, na.rm = TRUE))), na.rm = TRUE)
before_nan  <- sum(unlist(lapply(raw[protein_cols],
                                 function(v) sum(is.nan(v), na.rm = TRUE))), na.rm = TRUE)

raw[protein_cols] <- lapply(raw[protein_cols], function(v) {
  v[is.nan(v)] <- NA_real_   # NaN → NA
  v[v == 0]    <- NA_real_   # 0   → NA
  v
})

zeros_to_NA_count <- before_zero
nan_to_NA_count   <- before_nan

cat("Converted", zeros_to_NA_count, "zeros and", nan_to_NA_count, "NaNs to NA.\n")

# 5. Optional: further clean only if not already log2
if (!already_log2) {
  before_nonpos <- sum(unlist(lapply(raw[protein_cols],
                                     function(v) sum(is.finite(v) & v <= 0, na.rm = TRUE))), na.rm = TRUE)
  raw[protein_cols] <- lapply(raw[protein_cols], function(v) {
    v[v <= 0] <- NA_real_   # remove negatives
    v
  })
  after_nonpos <- sum(unlist(lapply(raw[protein_cols],
                                    function(v) sum(is.finite(v) & v <= 0, na.rm = TRUE))), na.rm = TRUE)
  zeros_to_NA_count <- zeros_to_NA_count + (before_nonpos - after_nonpos)
  cat("Removed additional", before_nonpos - after_nonpos, "negative or zero values before log2.\n")
}



# long format
long <- raw %>%
  tidyr::pivot_longer(dplyr::all_of(protein_cols), names_to = "Protein", values_to = "Intensity")

# standardize replicate labels and deduplicate tuple
long_clean <- long %>%
  dplyr::mutate(Replicate = trimws(as.character(.data[["Replicate"]]))) %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(c("Sample_ID","Replicate","Protein")))) %>%
  dplyr::summarise(Intensity = dplyr::first(stats::na.omit(Intensity)), .groups = "drop")

# wide Rep1/Rep2
w <- long_clean %>%
  tidyr::pivot_wider(id_cols = c("Sample_ID","Protein"),
                     names_from = Replicate, values_from = Intensity,
                     names_prefix = "Rep")

# harmonize column names to Rep1/Rep2 if needed
rep_cols_now <- setdiff(names(w), c("Sample_ID","Protein"))
col_1 <- rep_cols_now[grepl("^Rep\\s*1\\s*$", rep_cols_now, ignore.case = TRUE)]
col_2 <- rep_cols_now[grepl("^Rep\\s*2\\s*$", rep_cols_now, ignore.case = TRUE)]
if (length(col_1) == 1 && col_1 != "Rep1") names(w)[names(w) == col_1] <- "Rep1"
if (length(col_2) == 1 && col_2 != "Rep2") names(w)[names(w) == col_2] <- "Rep2"
if (!all(c("Rep1","Rep2") %in% names(w))) stop("Rep1/Rep2 columns not found. Check 'Replicate' contains only 1 and 2.")

# ---------- Log2 transform if needed (zeros already NA) ----------
log2_applied <- FALSE
if (!already_log2) {
  w <- w %>% dplyr::mutate(Rep1 = safe_log2(Rep1), Rep2 = safe_log2(Rep2))
  log2_applied <- TRUE
}

# -------------------- STAGE 1: RAW QC (always; pre-norm, pre-trim) ---
w_raw <- w
qc_before <- qc_stage(w_raw, "before_norm", raw_qc_dir)

# -------------------- STAGE 2: Decide normalization ------------------
DO_NORMALIZE_APPLIED <- FALSE
DECISION_NORM_MODE   <- if (identical(DO_NORMALIZE, "auto")) "auto" else "forced"

if (identical(DO_NORMALIZE, "auto")) {
  # Candidate normalization
  w_norm_cand <- w_raw %>% dplyr::mutate(Rep1 = median_normalize(Rep1), Rep2 = median_normalize(Rep2))
  qc_norm_cand <- qc_stage(w_norm_cand, "after_norm", raw_qc_dir)  # keep candidate plots under Raw QC for audit
  # Heuristic: adopt normalization if it IMPROVES overall QC
  med_rho_before <- suppressWarnings(median(qc_before$per_sample_rho$Spearman_rho, na.rm = TRUE))
  med_rho_after  <- suppressWarnings(median(qc_norm_cand$per_sample_rho$Spearman_rho, na.rm = TRUE))
  bias_before    <- abs(qc_before$summary$Global_BA_bias[1])
  bias_after     <- abs(qc_norm_cand$summary$Global_BA_bias[1])
  improve_rho    <- isTRUE(med_rho_after >= (med_rho_before + 0.02))
  reduce_bias    <- isTRUE(bias_after <= (bias_before - 0.10))
  # apply if either improvement criterion is met
  if (isTRUE(improve_rho) || isTRUE(reduce_bias)) {
    DO_NORMALIZE_APPLIED <- TRUE
    w_norm  <- w_norm_cand
    qc_after <- qc_norm_cand
    # also mirror the candidate plots into the stage folder
    dir.create(after_norm_dir, showWarnings = FALSE, recursive = TRUE)
    file.copy(list.files(raw_qc_dir, pattern = "after_norm_.*\\.png$", full.names = TRUE),
              after_norm_dir, overwrite = TRUE)
  } else {
    w_norm  <- w_raw
    qc_after <- qc_before
  }
} else if (isTRUE(DO_NORMALIZE)) {
  dir.create(after_norm_dir, showWarnings = FALSE, recursive = TRUE)
  w_norm <- w_raw %>% dplyr::mutate(Rep1 = median_normalize(Rep1), Rep2 = median_normalize(Rep2))
  qc_after <- qc_stage(w_norm, "after_norm", after_norm_dir)
  DO_NORMALIZE_APPLIED <- TRUE
} else {
  w_norm <- w_raw
  qc_after <- qc_before
}

# -------------------- STAGE 3: TRIM decision/apply (as needed) ------
# Compute Combined AFTER_NORM (or RAW if no norm) for extreme rule
long_src_after <- w_norm %>%
  dplyr::mutate(Combined = rowMeans(cbind(Rep1, Rep2), na.rm = TRUE),
                Combined = dplyr::if_else(is.nan(Combined), NA_real_, Combined)) %>%
  dplyr::select(Sample_ID, Protein, Rep1, Rep2, Combined)

# Extreme detection per protein
prot_iqr_after <- long_src_after %>%
  dplyr::group_by(Protein) %>%
  dplyr::summarise(Q1 = stats::quantile(Combined, 0.25, na.rm = TRUE),
                   Q3 = stats::quantile(Combined, 0.75, na.rm = TRUE),
                   IQR = Q3 - Q1, .groups = "drop") %>%
  dplyr::mutate(Low = Q1 - extreme_iqr_k * IQR, High = Q3 + extreme_iqr_k * IQR)

long_flagged_after <- long_src_after %>%
  dplyr::left_join(prot_iqr_after, by = "Protein") %>%
  dplyr::mutate(IsExtreme = is.finite(Combined) & (Combined < Low | Combined > High))

extreme_by_protein_after <- long_flagged_after %>%
  dplyr::group_by(Protein) %>%
  dplyr::summarise(
    n_extreme = sum(IsExtreme, na.rm = TRUE),
    n_nonNA   = sum(is.finite(Combined)),
    frac_extreme = ifelse(n_nonNA > 0, n_extreme / n_nonNA, 0),
    .groups = "drop"
  )

# Auto decision heuristic (only consulted if DO_EXTREME_TRIM == "auto")
median_rho_after <- suppressWarnings(median(qc_after$per_sample_rho$Spearman_rho, na.rm = TRUE))
ba_bias_after    <- qc_after$summary$Global_BA_bias[1]
frac_over_rule   <- mean(extreme_by_protein_after$frac_extreme >= protein_extreme_frac_drop, na.rm = TRUE)

if (identical(DO_EXTREME_TRIM, "auto")) {
  DO_EXTREME_TRIM_EFFECTIVE <- isTRUE(median_rho_after < 0.85) ||
    isTRUE(abs(ba_bias_after) > ba_abs_bias_max) ||
    isTRUE(frac_over_rule > 0.05)
  DECISION_SOURCE <- "auto"
} else {
  DO_EXTREME_TRIM_EFFECTIVE <- isTRUE(DO_EXTREME_TRIM)
  DECISION_SOURCE <- "forced"
}

TRIM_APPLIED <- isTRUE(DO_EXTREME_TRIM_EFFECTIVE)

if (TRIM_APPLIED) {
  dir.create(after_trim_dir, showWarnings = FALSE, recursive = TRUE)
  prot_ext_summary <- long_flagged_after %>%
    dplyr::group_by(Protein) %>%
    dplyr::summarise(
      n_extreme = sum(IsExtreme, na.rm = TRUE),
      n_nonNA   = sum(is.finite(Combined)),
      frac_extreme = ifelse(n_nonNA > 0, n_extreme / n_nonNA, 0),
      Remove_Protein = frac_extreme >= protein_extreme_frac_drop,
      .groups = "drop"
    )
  kept_after_extreme   <- prot_ext_summary %>% dplyr::filter(!Remove_Protein) %>% dplyr::pull(Protein)
  removed_extremes_tbl <- prot_ext_summary %>% dplyr::filter(Remove_Protein) %>% dplyr::arrange(dplyr::desc(frac_extreme))
  
  w_stage <- w_norm %>% dplyr::filter(Protein %in% kept_after_extreme)
  qc_trim <- qc_stage(w_stage, "after_trim", after_trim_dir)
  qc_stage_final <- qc_trim
  stage_name <- "after_trim"
  stage_dir  <- after_trim_dir
} else {
  removed_extremes_tbl <- tibble::tibble(Protein = character(), n_extreme = integer(),
                                         n_nonNA = integer(), frac_extreme = numeric(),
                                         Remove_Protein = logical())
  w_stage <- w_norm
  qc_stage_final <- qc_after
  stage_name <- if (DO_NORMALIZE_APPLIED) "after_norm" else "before_norm"
  stage_dir  <- if (DO_NORMALIZE_APPLIED) after_norm_dir else raw_qc_dir
}

# -------------------- Per-protein status across stages --------------
pp_before <- qc_before$per_protein %>%
  dplyr::select(Protein, QC_Category_relaxed) %>%
  dplyr::rename(QC_Category_relaxed_before = QC_Category_relaxed)

pp_after <- if (DO_NORMALIZE_APPLIED) {
  qc_after$per_protein %>%
    dplyr::select(Protein, QC_Category_relaxed) %>%
    dplyr::rename(QC_Category_relaxed_after_norm = QC_Category_relaxed)
} else NULL

pp_trim <- if (TRIM_APPLIED) {
  qc_trim$per_protein %>%
    dplyr::select(Protein, QC_Category_relaxed) %>%
    dplyr::rename(QC_Category_relaxed_after_trim = QC_Category_relaxed)
} else NULL

pp_master <- pp_before
if (!is.null(pp_after)) pp_master <- dplyr::full_join(pp_master, pp_after, by = "Protein")
if (!is.null(pp_trim))  pp_master <- dplyr::full_join(pp_master,  pp_trim,  by = "Protein")

pp_master <- pp_master %>%
  dplyr::left_join(
    if (TRIM_APPLIED) removed_extremes_tbl %>% dplyr::select(Protein, Remove_Protein) else
      tibble::tibble(Protein = character(), Remove_Protein = logical()),
    by = "Protein"
  ) %>%
  dplyr::mutate(Removed_by_trim = dplyr::coalesce(Remove_Protein, FALSE),
                Removed_stage   = dplyr::if_else(Removed_by_trim, "after_trim", NA_character_)) %>%
  dplyr::select(-Remove_Protein)

# -------------------- Build FINAL stage matrices/plots --------------
long_src_for_export <- w_stage %>%
  dplyr::mutate(Combined = rowMeans(cbind(Rep1, Rep2), na.rm = TRUE),
                Combined = dplyr::if_else(is.nan(Combined), NA_real_, Combined)) %>%
  dplyr::select(Sample_ID, Protein, Rep1, Rep2, Combined)

coverage_tbl <- long_src_for_export %>%
  dplyr::group_by(Protein) %>%
  dplyr::summarise(Coverage = mean(is.finite(Combined)), .groups = "drop")
keep_by_cov <- coverage_tbl %>% dplyr::mutate(Keep_by_coverage = Coverage >= coverage_min)
kept_proteins <- keep_by_cov %>% dplyr::filter(Keep_by_coverage) %>% dplyr::pull(Protein)

# PCA on final
pca_in <- make_pca_matrix(
  long_src_for_export %>% dplyr::filter(Protein %in% kept_proteins) %>%
    dplyr::select(Sample_ID, Protein, Combined),
  id_cols_no_rep = "Sample_ID", cov_min = coverage_min
)
pca <- run_pca(pca_in$X)
ve <- round(100 * var_explained(pca)[1:2], 1)
scores <- as.data.frame(pca$x[, 1:2, drop = FALSE]); scores$Sample_ID <- pca_in$meta$Sample_ID
p_pca_scores <- ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2)) + ggplot2::geom_point(alpha = 0.85) +
  ggplot2::labs(title = paste0(NAME_BASE, " — PCA scores (", stage_name, ")"),
                x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"))
ve_all <- round(100 * var_explained(pca), 1)
p_pca_var <- ggplot2::ggplot(data.frame(PC = seq_along(ve_all), Var = ve_all),
                             ggplot2::aes(PC, Var)) + ggplot2::geom_col() +
  ggplot2::labs(title = paste0(NAME_BASE, " — PCA variance (%) (", stage_name, ")"),
                x = "PC", y = "% variance")

# Sample correlation heatmap (final stage) + capture for PDF
wide_source <- long_src_for_export %>%
  dplyr::filter(Protein %in% kept_proteins) %>%
  dplyr::group_by(Sample_ID, Protein) %>%
  dplyr::summarise(Combined = mean(Combined, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Protein, values_from = Combined)
mat_corr <- as.matrix(wide_source[, setdiff(names(wide_source), "Sample_ID")])
for (j in seq_len(ncol(mat_corr))) {
  med <- stats::median(mat_corr[, j], na.rm = TRUE)
  idx <- !is.finite(mat_corr[, j])
  if (any(idx)) mat_corr[idx, j] <- med
}
cormat <- suppressWarnings(stats::cor(t(mat_corr), method = "spearman"))
rownames(cormat) <- colnames(cormat) <- wide_source$Sample_ID

# --- Quantify inter-sample variability from heatmap correlations ---
sample_mean_corr <- rowMeans(cormat, na.rm = TRUE)
mean_corr <- mean(sample_mean_corr, na.rm = TRUE)
sd_corr   <- sd(sample_mean_corr, na.rm = TRUE)
cv_corr   <- sd_corr / mean_corr * 100

cat("\n=== Inter-sample correlation metrics ===\n")
cat("Mean inter-sample correlation:", round(mean_corr, 3),
    "| SD:", round(sd_corr, 3),
    "| CV:", round(cv_corr, 1), "%\n")

# Optional thresholds for QC interpretation
if (mean_corr >= 0.9 && cv_corr < 5) {
  cat("→ Excellent consistency across samples — no detectable batch effect.\n")
} else if (mean_corr >= 0.85 && cv_corr < 10) {
  cat("→ Acceptable variability; check if biological grouping explains structure.\n")
} else {
  cat("→ Elevated variability — potential batch or normalization issue.\n")
}


# Save heatmap PNG for Excel
pheatmap::pheatmap(cormat,
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "correlation",
                   color = colorRampPalette(brewer.pal(9, "RdBu"))(255),
                   filename = file.path(final_dir, paste0(NAME_BASE, "_SampleCorr_Heatmap_", stage_name, ".png")),
                   width = 9, height = 9)

# Convert pheatmap to a grob so it can be embedded in the PDF
library(grid)
heatmap_obj  <- pheatmap::pheatmap(cormat,
                                   clustering_distance_rows = "correlation",
                                   clustering_distance_cols = "correlation",
                                   color = colorRampPalette(brewer.pal(9, "RdBu"))(255),
                                   silent = TRUE)
heatmap_grob <- grid::grid.grabExpr(grid::grid.draw(heatmap_obj$gtable))

# -------------------- Save FINAL figures (final-only PDF) -----------
p_spearman_final <- qc_stage_final$plots$spearman
p_ba_final       <- qc_stage_final$plots$bland_altman
p_dens_final     <- qc_stage_final$plots$density

if (WRITE_PDF) {
  combo <- (p_spearman_final | p_ba_final) /
    (p_dens_final     | p_pca_var)  /
    (p_pca_scores     | patchwork::wrap_elements(heatmap_grob))
  ggplot2::ggsave(file.path(final_dir, paste0(NAME_BASE, "_QC_Summary_", stage_name, ".pdf")),
                  combo, width = 13, height = 13)
}

# -------------------- Tables & Exports -------------------
# Meta table (one row per Sample_ID, keep first non-NA group if present)
meta_only <- raw %>%
  dplyr::select(dplyr::all_of(id_cols), dplyr::any_of(ignore_cols)) %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(setdiff(id_cols, "Replicate")))) %>%
  dplyr::summarise(dplyr::across(dplyr::any_of(ignore_cols), ~ { v <- stats::na.omit(.x); if (length(v)==0) NA else v[1] }),
                   .groups = "drop") %>%
  dplyr::distinct(dplyr::across(dplyr::all_of(setdiff(id_cols, "Replicate"))), .keep_all = TRUE)

# Final correlation inputs (kept proteins only)
cor_long <- long_src_for_export %>%
  dplyr::left_join(meta_only, by = setdiff(id_cols, "Replicate")) %>%
  dplyr::left_join(coverage_tbl, by = "Protein") %>%
  dplyr::left_join(qc_stage_final$per_protein %>% dplyr::select(Protein, QC_Category_relaxed), by = "Protein") %>%
  dplyr::left_join(keep_by_cov %>% dplyr::select(Protein, Keep_by_coverage), by = "Protein")

cor_matrix <- wide_source %>%
  dplyr::left_join(meta_only, by = "Sample_ID") %>%
  dplyr::select(Sample_ID, dplyr::any_of(ignore_cols), dplyr::everything())

# -------------------- Consolidated Final Results + Audit workbook ----
fix_write <- function(wb, sheet, x) openxlsx::writeData(wb, sheet, x = as.data.frame(x), na.string = NA_STRING)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "README")
fix_write(wb, "README", tibble::tibble(
  Section = c("Purpose","Final Stage","Settings","Normalization decision","Trim decision","Provenance"),
  Details = c(
    "Stage-aware replicate QC with optional normalization and auto-trim; PCA unchanged; final-only PDF; this workbook consolidates audit + key data sheets.",
    stage_name,
    paste0("coverage_min=", coverage_min,
           "; already_log2=", already_log2,
           "; DO_NORMALIZE=", DO_NORMALIZE, " (applied? ", DO_NORMALIZE_APPLIED, ")",
           "; DO_EXTREME_TRIM=", as.character(DO_EXTREME_TRIM)),
    if (identical(DO_NORMALIZE, "auto"))
      paste0("median rho before/after = ",
             sprintf("%.3f", suppressWarnings(median(qc_before$per_sample_rho$Spearman_rho, na.rm = TRUE))), " / ",
             sprintf("%.3f", suppressWarnings(median(qc_after$per_sample_rho$Spearman_rho,  na.rm = TRUE))),
             "; |BA bias| before/after = ",
             sprintf("%.3f", abs(qc_before$summary$Global_BA_bias[1])), " / ",
             sprintf("%.3f", abs(qc_after$summary$Global_BA_bias[1])),
             "; normalization applied? ", DO_NORMALIZE_APPLIED)
    else
      paste0("Normalization forced to ", isTRUE(DO_NORMALIZE), "."),
    if (identical(DO_EXTREME_TRIM, "auto"))
      paste0("Auto metrics — median Spearman(after)=", sprintf("%.3f", median_rho_after),
             "; |BA bias(after)|=", sprintf("%.3f", abs(ba_bias_after)),
             "; %proteins over extreme rule=", sprintf("%.1f%%", 100*frac_over_rule),
             "; TRIM_APPLIED=", TRIM_APPLIED)
    else
      paste0("Trimming forced: ", isTRUE(DO_EXTREME_TRIM), "."),
    paste("file_path:", file_path, "| sheet:", sheet_name, "| generated:", as.character(Sys.time()))
  )
))

# --- Crucial audit sheets & data you actually use
openxlsx::addWorksheet(wb, "QC_Summary_Final");         fix_write(wb, "QC_Summary_Final", qc_stage_final$summary)
openxlsx::addWorksheet(wb, "PerSample_Spearman_Final"); fix_write(wb, "PerSample_Spearman_Final", qc_stage_final$per_sample_rho)
openxlsx::addWorksheet(wb, "QC_PerProtein_Final");      fix_write(wb, "QC_PerProtein_Final", qc_stage_final$per_protein)
openxlsx::addWorksheet(wb, "Coverage_Final");           fix_write(wb, "Coverage_Final", coverage_tbl)
openxlsx::addWorksheet(wb, "Kept_by_Coverage_Final");   fix_write(wb, "Kept_by_Coverage_Final", keep_by_cov)
openxlsx::addWorksheet(wb, "PerProtein_Status_All");    fix_write(wb, "PerProtein_Status_All", pp_master)
openxlsx::addWorksheet(wb, "Correlation_Matrix");       fix_write(wb, "Correlation_Matrix", cor_matrix)
openxlsx::addWorksheet(wb, "Correlation_Long");         fix_write(wb, "Correlation_Long", cor_long)

# Embed heatmap PNG if present
png_heat <- file.path(final_dir, paste0(NAME_BASE, "_SampleCorr_Heatmap_", stage_name, ".png"))
if (file.exists(png_heat)) {
  openxlsx::addWorksheet(wb, "Heatmap")
  openxlsx::writeData(wb, "Heatmap", data.frame(Note = "Sample×Sample Spearman correlation (final stage)"))
  openxlsx::insertImage(wb, "Heatmap", png_heat, startRow = 3, startCol = 1, width = 18, height = 18, dpi = 300)
}

for (sh in openxlsx::sheets(wb)) {
  openxlsx::setColWidths(wb, sh, cols = 1:200, widths = "auto")
  openxlsx::addStyle(wb, sh, style = openxlsx::createStyle(textDecoration = "bold"),
                     rows = 1, cols = 1:200, gridExpand = TRUE)
}
final_xlsx <- file.path(final_dir, paste0(NAME_BASE, "_FINAL_Results_Audit_", stage_name, ".xlsx"))
openxlsx::saveWorkbook(wb, final_xlsx, overwrite = TRUE)

# -------------------- Console Summary -------------------
n_samples   <- dplyr::n_distinct(cor_matrix$Sample_ID)
n_kept_cov  <- length(kept_proteins)
final_pdf   <- file.path(final_dir, paste0(NAME_BASE, "_QC_Summary_", stage_name, ".pdf"))
preferred   <- paste0(basename(final_xlsx), " → sheet: 'Correlation_Matrix'")

cat("\n=== SUMMARY (stage-aware) ===\n")
cat("Normalization mode:", DECISION_NORM_MODE, "| Applied?:", DO_NORMALIZE_APPLIED, "\n")
cat("Trim decision mode:", DECISION_SOURCE,     "| Applied?:", TRIM_APPLIED, "\n")
cat("Median Spearman (", stage_name, "): ",
    sprintf("%.3f", median(qc_stage_final$per_sample_rho$Spearman_rho, na.rm = TRUE)),
    "\n", sep = "")
cat("Bland-Altman | bias:", sprintf("%.3f", qc_stage_final$summary$Global_BA_bias),
    "LoA:[", sprintf("%.3f", qc_stage_final$summary$Global_BA_LOA_lower), ",",
    sprintf("%.3f", qc_stage_final$summary$Global_BA_LOA_upper), "]\n")
cat("coverage_min:", coverage_min, "| Kept proteins (final by coverage):", n_kept_cov, "| Samples:", n_samples, "\n")
cat("Final stage:", stage_name, "\n")
cat("Final outputs dir:", final_dir, "\n")
cat("PDF saved:", final_pdf, "\n")
if (!already_log2) {
  cat("Pre-log2 cleanup: zeros->NA count (approx):", zeros_to_NA_count, "| log2 applied?:", log2_applied, "\n")
} else {
  cat("Input declared already_log2=TRUE; no log2 transformation applied.\n")
}
cat("\n>>> Recommended file for downstream analysis:", preferred, "\n")


