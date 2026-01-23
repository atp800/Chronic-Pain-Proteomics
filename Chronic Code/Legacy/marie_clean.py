import os
import sys
import warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Try to import pingouin for ICC, handle gracefully if missing
try:
    import pingouin as pg
    HAS_PINGOUIN = True
except ImportError:
    HAS_PINGOUIN = False
    warnings.warn("Package 'pingouin' not found. ICC values will be NA. Install via `pip install pingouin`.")

# ======================================================
# DIA Replicate QC — Python Port
# Assess → (optional) Normalize → (optional) Auto-Trim → Export
# ======================================================

# -------------------- User Parameters --------------------
# Note: Python uses forward slashes or raw strings (r"...") for Windows paths
FILE_PATH       = "Chronic Code/Data_To_Clean.xlsx"
SHEET_NAME      = "Soluble - Group"
ID_COLS         = ["Sample_ID", "Replicate"]      # MUST exist in sheet
IGNORE_COLS     = ["Condition"]                   # optional metadata

# If your sheet already contains log2 intensities, set TRUE
ALREADY_LOG2    = True

# Toggle median normalization: True, False, or "auto"
DO_NORMALIZE    = "auto"

# Trimming strategy: True, False, or "auto"
DO_EXTREME_TRIM = "auto"

# Extreme-value rules
EXTREME_IQR_K             = 3.0
PROTEIN_EXTREME_FRAC_DROP = 0.10

# --- RELAXED thresholds ---
COVERAGE_MIN     = 0.125
ICC_RELIABLE_MIN = 0.50
ICC_INSPECT_MIN  = 0.30
CV_RELIABLE_MAX  = 35
CV_INSPECT_MAX   = 60
BA_ABS_BIAS_MAX  = 0.60
ALPHA_MCNEMAR    = 0.001

# Output folders
BASE_DIR = "Chronic Code/Clean Output"
DATE_STR = pd.Timestamp.now().strftime("%Y-%m-%d")
ROOT_DIR = os.path.join(BASE_DIR, DATE_STR, SHEET_NAME)

os.makedirs(ROOT_DIR, exist_ok=True)

# Stage subfolders
RAW_QC_DIR     = os.path.join(ROOT_DIR, "Raw QC")
AFTER_NORM_DIR = os.path.join(ROOT_DIR, "After norm")
AFTER_TRIM_DIR = os.path.join(ROOT_DIR, "After trim")
FINAL_DIR      = os.path.join(ROOT_DIR, "Final Outputs")

os.makedirs(RAW_QC_DIR, exist_ok=True)
os.makedirs(FINAL_DIR, exist_ok=True)

DATA_STUB = os.path.splitext(os.path.basename(FILE_PATH))[0]
NAME_BASE = f"{DATA_STUB}_{SHEET_NAME}"

WRITE_PDF = True
SEED = 20250827
np.random.seed(SEED)

# -------------------- Helpers ---------------------------

def safe_log2(x):
    """Log2 that handles zeros/negative by masking."""
    # Create a copy to avoid SettingWithCopy warnings on slices
    x = x.copy()
    mask = (x > 0)
    out = np.full(x.shape, np.nan)
    out[mask] = np.log2(x[mask])
    return out

def median_normalize(v):
    """Subtract median (ignoring NaNs)."""
    return v - np.nanmedian(v)

def drop_zero_var(df):
    """Drop numeric columns with zero variance or insufficient data."""
    numerics = df.select_dtypes(include=[np.number])
    keep = []
    for col in numerics.columns:
        v = numerics[col]
        if np.sum(np.isfinite(v)) >= 3 and np.nanstd(v) > 0:
            keep.append(col)
    return numerics[keep]

def calculate_icc_for_group(g):
    """Helper to calculate QC metrics for a single protein group."""
    # g is a DataFrame with Rep1, Rep2 for one Protein
    v1 = g['Rep1'].values.astype(float)
    v2 = g['Rep2'].values.astype(float)
    
    # Mask for pairs
    mask_pair = np.isfinite(v1) & np.isfinite(v2)
    v1p = v1[mask_pair]
    v2p = v2[mask_pair]
    n_pairs = len(v1p)
    
    # ICC
    icc_val = np.nan
    if n_pairs >= 2 and HAS_PINGOUIN:
        try:
            # Create long format for pingouin
            data_tmp = pd.DataFrame({
                'Targets': np.tile(np.arange(n_pairs), 2),
                'Raters': np.repeat(['A', 'B'], n_pairs),
                'Ratings': np.concatenate([v1p, v2p])
            })
            # Equivalent to R irr::icc(model="twoway", type="agreement", unit="single") -> ICC2
            icc_res = pg.intraclass_correlation(data=data_tmp, targets='Targets', raters='Raters', ratings='Ratings')
            # Extract ICC2 (Single random raters)
            icc_val = icc_res.set_index('Type').loc['ICC2', 'ICC']
        except:
            icc_val = np.nan
            
    # CV (Median of percentages)
    cv_med_pct = np.nan
    if n_pairs >= 1:
        a_lin = 2**v1p
        b_lin = 2**v2p
        m = (a_lin + b_lin) / 2.0
        # Avoid division by zero
        m[m == 0] = np.nan
        s = np.sqrt(((a_lin - m)**2 + (b_lin - m)**2)) # SD of 2 values
        cvs = 100 * (s / m)
        cv_med_pct = np.nanmedian(cvs)
        
    # Bland Altman Stats
    diffs = v1p - v2p
    diff_mean = np.mean(diffs) if n_pairs >= 1 else np.nan
    diff_sd   = np.std(diffs, ddof=1) if n_pairs >= 2 else np.nan
    ba_low    = diff_mean - 1.96 * diff_sd if not np.isnan(diff_sd) else np.nan
    ba_high   = diff_mean + 1.96 * diff_sd if not np.isnan(diff_sd) else np.nan
    
    # Missingness
    # Rep1 present, Rep2 missing
    miss_b = np.sum((~np.isnan(v1)) & (np.isnan(v2)))
    # Rep1 missing, Rep2 present
    miss_c = np.sum((np.isnan(v1)) & (~np.isnan(v2)))
    disc = miss_b + miss_c
    if disc > 0:
        # McNemar exact equivalent: 2 * binomial_cdf(min(b,c), n=disc, p=0.5)
        miss_p = min(1.0, 2 * stats.binom.cdf(min(miss_b, miss_c), disc, 0.5))
    else:
        miss_p = np.nan

    return pd.Series({
        'n_pairs': n_pairs,
        'ICC_value': icc_val,
        'CV_med_pct': cv_med_pct,
        'Diff_log2_mean': diff_mean,
        'Diff_log2_sd': diff_sd,
        'BA_LOA_lower': ba_low,
        'BA_LOA_upper': ba_high,
        'Miss_p': miss_p
    })

def qc_per_protein_relaxed(w_df):
    """Compute per-protein metrics."""
    # Group by Protein and apply calculation
    # Using a list comprehension or apply can be slow, but is most direct translation.
    # We optimize by grouping and applying the function defined above.
    
    metrics = w_df.groupby("Protein").apply(calculate_icc_for_group)
    metrics = metrics.reset_index()
    
    # Classify
    conditions = [
        (abs(metrics['Diff_log2_mean']) > BA_ABS_BIAS_MAX),
        (metrics['Miss_p'] < ALPHA_MCNEMAR)
    ]
    # Note: NaNs in conditions propagate as False usually in numpy comparison, 
    # but let's handle explicitly.
    bad_bias = (metrics['Diff_log2_mean'].notna() & (metrics['Diff_log2_mean'].abs() > BA_ABS_BIAS_MAX))
    bad_miss = (metrics['Miss_p'].notna() & (metrics['Miss_p'] < ALPHA_MCNEMAR))
    
    # ICC Band
    icc_band = []
    for x in metrics['ICC_value']:
        if pd.isna(x): icc_band.append("drop")
        elif x >= ICC_RELIABLE_MIN: icc_band.append("reliable")
        elif x >= ICC_INSPECT_MIN: icc_band.append("inspect")
        else: icc_band.append("drop")
    metrics['ICC_band'] = icc_band

    # CV Band
    cv_band = []
    for x in metrics['CV_med_pct']:
        if pd.isna(x): cv_band.append("drop")
        elif x <= CV_RELIABLE_MAX: cv_band.append("reliable")
        elif x <= CV_INSPECT_MAX: cv_band.append("inspect")
        else: cv_band.append("drop")
    metrics['CV_band'] = cv_band
    
    # Final Category
    cats = []
    for i in range(len(metrics)):
        if bad_bias[i] or bad_miss[i]:
            cats.append("Drop")
        elif icc_band[i] == "drop" or cv_band[i] == "drop":
            cats.append("Drop")
        elif icc_band[i] == "reliable" and cv_band[i] == "reliable":
            cats.append("Reliable")
        else:
            cats.append("Inspect")
    metrics['QC_Category_relaxed'] = cats
    
    return metrics

def run_qc(w_df, label):
    # Pair check
    w_df['pair_ok'] = w_df['Rep1'].notna() & w_df['Rep2'].notna()
    pair_tbl = w_df.groupby("Sample_ID")['pair_ok'].sum().reset_index(name='n_pairs')
    
    # Spearman Rho per sample
    rhos = []
    ids = []
    for sid, grp in w_df.groupby("Sample_ID"):
        sub = grp.dropna(subset=['Rep1', 'Rep2'])
        r = np.nan
        if len(sub) >= 2:
            r, _ = stats.spearmanr(sub['Rep1'], sub['Rep2'])
        rhos.append(r)
        ids.append(sid)
    per_sample_rho = pd.DataFrame({'Sample_ID': ids, 'Spearman_rho': rhos})
    
    # Global pairs for Bland-Altman
    # Filter finite values
    valid = w_df.dropna(subset=['Rep1', 'Rep2'])
    mean_log2 = (valid['Rep1'] + valid['Rep2']) / 2
    diff_log2 = valid['Rep1'] - valid['Rep2']
    
    all_pairs = pd.DataFrame({'Mean_log2': mean_log2, 'Diff_log2': diff_log2})
    
    ba_bias = all_pairs['Diff_log2'].mean()
    ba_sd   = all_pairs['Diff_log2'].std(ddof=1)
    ba_lwr  = ba_bias - 1.96 * ba_sd
    ba_upr  = ba_bias + 1.96 * ba_sd
    ba_text = f"Bias = {ba_bias:.3f} | LoA [{ba_lwr:.3f}, {ba_upr:.3f}] | SD = {ba_sd:.3f}"
    
    # --- Plots ---
    # 1. Spearman Histogram
    fig_sp, ax_sp = plt.subplots(figsize=(6.5, 4.5))
    rho_clean = per_sample_rho['Spearman_rho'].dropna()
    if len(rho_clean) > 0:
        sns.histplot(rho_clean, bins=30, ax=ax_sp)
        ax_sp.axvline(0.8, linestyle="--", color='k')
        ax_sp.axvline(0.9, linestyle="--", color='k')
    ax_sp.set_title(f"Per-sample Spearman rho ({label})")
    ax_sp.set_xlabel("Spearman rho")
    plt.close(fig_sp) # Close to prevent display, we just want object if possible, but matplotlib relies on state. 
                      # We will keep the figure object for saving.

    # 2. Bland-Altman
    fig_ba, ax_ba = plt.subplots(figsize=(7.5, 4.5))
    ax_ba.scatter(all_pairs['Mean_log2'], all_pairs['Diff_log2'], alpha=0.3, s=3)
    ax_ba.axhline(0, linestyle=":", color='gray')
    ax_ba.axhline(ba_bias, color='red')
    ax_ba.axhline(ba_lwr, linestyle="--", color='red')
    ax_ba.axhline(ba_upr, linestyle="--", color='red')
    ax_ba.set_title(f"Global Bland–Altman ({label})\n{ba_text}")
    ax_ba.set_xlabel("Mean (log2)")
    ax_ba.set_ylabel("Rep1 - Rep2 (log2)")
    plt.close(fig_ba)
    
    # 3. Density
    fig_dens, ax_dens = plt.subplots(figsize=(6.5, 4.5))
    # Melt for density
    long_pool = w_df[['Rep1', 'Rep2']].melt(var_name='Rep', value_name='Log2').dropna()
    sns.kdeplot(data=long_pool, x='Log2', hue='Rep', ax=ax_dens)
    ax_dens.set_title(f"Pooled intensity distributions ({label})")
    plt.close(fig_dens)
    
    # Global Spearman
    valid_global = w_df.dropna(subset=['Rep1', 'Rep2'])
    glob_rho, _ = stats.spearmanr(valid_global['Rep1'], valid_global['Rep2'])

    summary = pd.DataFrame([{
        'Label': label,
        'Global_Spearman_rho': glob_rho,
        'Global_BA_bias': ba_bias,
        'Global_BA_LOA_lower': ba_lwr,
        'Global_BA_LOA_upper': ba_upr
    }])
    
    return {
        'summary': summary,
        'plots': {'spearman': fig_sp, 'bland_altman': fig_ba, 'density': fig_dens},
        'per_sample_rho': per_sample_rho,
        'per_protein': qc_per_protein_relaxed(w_df), # Calculation
        'pair_check': pair_tbl,
        'all_pairs': all_pairs
    }

def qc_stage(w_df, label, out_dir):
    qc = run_qc(w_df, label)
    # Save plots
    stub = os.path.join(out_dir, f"{NAME_BASE}_{label}")
    qc['plots']['spearman'].savefig(f"{stub}_Spearman_Hist.png", dpi=300)
    qc['plots']['bland_altman'].savefig(f"{stub}_BlandAltman.png", dpi=300)
    qc['plots']['density'].savefig(f"{stub}_Density.png", dpi=300)
    return qc

def make_pca_matrix(df, id_cols_no_rep, cov_min):
    # Pivot
    wide = df.pivot(index=id_cols_no_rep, columns='Protein', values='Combined')
    
    # Filter by coverage
    cov = wide.notna().mean(axis=0)
    keep = cov[cov >= cov_min].index.tolist()
    
    if len(keep) < 2:
        # Keep top 200 if too few
        keep = cov.sort_values(ascending=False).head(200).index.tolist()
        
    mat = wide[keep]
    
    # Impute missing with median of column
    for col in mat.columns:
        med = mat[col].median()
        mat[col] = mat[col].fillna(med)
        
    # Drop zero variance
    mat_clean = drop_zero_var(mat)
    
    return mat_clean, keep

# -------------------- Load & Prepare ---------------------

print("Loading data...")
raw0 = pd.read_excel(FILE_PATH, sheet_name=SHEET_NAME)

# Check columns
missing_ids = [c for c in ID_COLS if c not in raw0.columns]
if missing_ids:
    sys.exit(f"Missing id_cols: {missing_ids}")

meta_cols = list(set(ID_COLS + IGNORE_COLS).intersection(raw0.columns))
protein_cols = [c for c in raw0.columns if c not in meta_cols]

if not protein_cols:
    sys.exit("No protein columns detected.")

# Numeric conversion
raw = raw0.copy()
for col in protein_cols:
    raw[col] = pd.to_numeric(raw[col], errors='coerce')

# Pre-process: 0 and NaN -> NaN
# In pandas, 0 is often distinct from NaN.
zeros_to_NA_count = 0
before_zero = (raw[protein_cols] == 0).sum().sum()
before_nan  = raw[protein_cols].isna().sum().sum()

raw[protein_cols] = raw[protein_cols].replace(0, np.nan)

# Update count logic to match R script
# R logic: count zeros, count nan, convert.
# Here we just counted.
zeros_to_NA_count = before_zero 
print(f"Converted {zeros_to_NA_count} zeros and {before_nan} NaNs to NA.")

if not ALREADY_LOG2:
    # Check negatives
    neg_mask = (raw[protein_cols] <= 0)
    before_nonpos = neg_mask.sum().sum()
    raw[protein_cols] = raw[protein_cols].mask(neg_mask, np.nan)
    
    after_nonpos = (raw[protein_cols] <= 0).sum().sum() # Should be 0
    zeros_to_NA_count += (before_nonpos - after_nonpos)
    print(f"Removed additional {before_nonpos} negative or zero values before log2.")

# Long format
long_df = raw.melt(id_vars=meta_cols, value_vars=protein_cols, 
                   var_name="Protein", value_name="Intensity")

# Standardize Replicate and deduplicate
long_df['Replicate'] = long_df['Replicate'].astype(str).str.strip()
long_clean = long_df.dropna(subset=['Intensity']).groupby(
    ["Sample_ID", "Replicate", "Protein"], as_index=False
)['Intensity'].first()

# Wide Rep1/Rep2
w = long_clean.pivot(index=["Sample_ID", "Protein"], columns="Replicate", values="Intensity").reset_index()

# Harmonize column names (look for "1" and "2" approx)
found_cols = w.columns.tolist()
rep1_col = next((c for c in found_cols if "1" in str(c) and c not in ["Sample_ID", "Protein"]), None)
rep2_col = next((c for c in found_cols if "2" in str(c) and c not in ["Sample_ID", "Protein"]), None)

if rep1_col: w.rename(columns={rep1_col: 'Rep1'}, inplace=True)
if rep2_col: w.rename(columns={rep2_col: 'Rep2'}, inplace=True)

if 'Rep1' not in w.columns or 'Rep2' not in w.columns:
    sys.exit("Rep1/Rep2 columns not found. Check Replicate column content.")

log2_applied = False
if not ALREADY_LOG2:
    w['Rep1'] = safe_log2(w['Rep1'])
    w['Rep2'] = safe_log2(w['Rep2'])
    log2_applied = True

# -------------------- STAGE 1: RAW QC --------------------
print("Running Raw QC...")
w_raw = w.copy()
qc_before = qc_stage(w_raw, "before_norm", RAW_QC_DIR)

# -------------------- STAGE 2: Normalization -------------
print("Evaluating Normalization...")
DO_NORMALIZE_APPLIED = False
DECISION_NORM_MODE = "auto" if DO_NORMALIZE == "auto" else "forced"

w_norm = w_raw.copy()
qc_after = qc_before

if DO_NORMALIZE == "auto":
    # Candidate
    w_cand = w_raw.copy()
    w_cand['Rep1'] = median_normalize(w_cand['Rep1'])
    w_cand['Rep2'] = median_normalize(w_cand['Rep2'])
    qc_cand = qc_stage(w_cand, "after_norm", RAW_QC_DIR) # Save to Raw QC for comparison
    
    med_rho_before = qc_before['per_sample_rho']['Spearman_rho'].median()
    med_rho_after  = qc_cand['per_sample_rho']['Spearman_rho'].median()
    bias_before    = abs(qc_before['summary']['Global_BA_bias'].iloc[0])
    bias_after     = abs(qc_cand['summary']['Global_BA_bias'].iloc[0])
    
    improve_rho = (med_rho_after >= (med_rho_before + 0.02))
    reduce_bias = (bias_after <= (bias_before - 0.10))
    
    if improve_rho or reduce_bias:
        DO_NORMALIZE_APPLIED = True
        w_norm = w_cand
        qc_after = qc_cand
        # Copy plots
        if not os.path.exists(AFTER_NORM_DIR): os.makedirs(AFTER_NORM_DIR)
        import shutil
        for f in os.listdir(RAW_QC_DIR):
            if "after_norm" in f and f.endswith(".png"):
                shutil.copy2(os.path.join(RAW_QC_DIR, f), AFTER_NORM_DIR)
    else:
        # Revert
        qc_after = qc_before

elif DO_NORMALIZE is True:
    if not os.path.exists(AFTER_NORM_DIR): os.makedirs(AFTER_NORM_DIR)
    w_norm['Rep1'] = median_normalize(w_norm['Rep1'])
    w_norm['Rep2'] = median_normalize(w_norm['Rep2'])
    qc_after = qc_stage(w_norm, "after_norm", AFTER_NORM_DIR)
    DO_NORMALIZE_APPLIED = True

# -------------------- STAGE 3: Trimming ------------------
print("Evaluating Trimming...")
# Calc combined for extreme detection
w_norm['Combined'] = w_norm[['Rep1', 'Rep2']].mean(axis=1)

# Extreme detection per protein
prot_iqr = w_norm.groupby('Protein')['Combined'].quantile([0.25, 0.75]).unstack()
prot_iqr['IQR'] = prot_iqr[0.75] - prot_iqr[0.25]
prot_iqr['Low'] = prot_iqr[0.25] - EXTREME_IQR_K * prot_iqr['IQR']
prot_iqr['High'] = prot_iqr[0.75] + EXTREME_IQR_K * prot_iqr['IQR']

w_flagged = w_norm.merge(prot_iqr[['Low', 'High']], on='Protein', how='left')
w_flagged['IsExtreme'] = (w_flagged['Combined'] < w_flagged['Low']) | (w_flagged['Combined'] > w_flagged['High'])

# Summary per protein
def summarize_extreme(x):
    n_ext = x['IsExtreme'].sum()
    n_fin = x['Combined'].notna().sum()
    return pd.Series({'n_extreme': n_ext, 'n_nonNA': n_fin})

ext_summary = w_flagged.groupby('Protein').apply(summarize_extreme).reset_index()
ext_summary['frac_extreme'] = np.where(ext_summary['n_nonNA'] > 0, ext_summary['n_extreme'] / ext_summary['n_nonNA'], 0)

# Heuristics
med_rho_after = qc_after['per_sample_rho']['Spearman_rho'].median()
ba_bias_after = qc_after['summary']['Global_BA_bias'].iloc[0]
frac_over_rule = (ext_summary['frac_extreme'] >= PROTEIN_EXTREME_FRAC_DROP).mean()

DECISION_SOURCE = "forced"
DO_EXTREME_TRIM_EFFECTIVE = False

if DO_EXTREME_TRIM == "auto":
    DECISION_SOURCE = "auto"
    if (med_rho_after < 0.85) or (abs(ba_bias_after) > BA_ABS_BIAS_MAX) or (frac_over_rule > 0.05):
        DO_EXTREME_TRIM_EFFECTIVE = True
else:
    DO_EXTREME_TRIM_EFFECTIVE = bool(DO_EXTREME_TRIM)

TRIM_APPLIED = DO_EXTREME_TRIM_EFFECTIVE

removed_extremes_tbl = pd.DataFrame()
w_stage = w_norm.copy()
qc_stage_final = qc_after
stage_name = "after_norm" if DO_NORMALIZE_APPLIED else "before_norm"
stage_dir = AFTER_NORM_DIR if DO_NORMALIZE_APPLIED else RAW_QC_DIR

if TRIM_APPLIED:
    if not os.path.exists(AFTER_TRIM_DIR): os.makedirs(AFTER_TRIM_DIR)
    
    ext_summary['Remove_Protein'] = ext_summary['frac_extreme'] >= PROTEIN_EXTREME_FRAC_DROP
    kept_prots = ext_summary.loc[~ext_summary['Remove_Protein'], 'Protein'].tolist()
    removed_extremes_tbl = ext_summary.loc[ext_summary['Remove_Protein']].sort_values('frac_extreme', ascending=False)
    
    w_stage = w_norm[w_norm['Protein'].isin(kept_prots)].copy()
    qc_trim = qc_stage(w_stage, "after_trim", AFTER_TRIM_DIR)
    
    qc_stage_final = qc_trim
    stage_name = "after_trim"
    stage_dir = AFTER_TRIM_DIR

# -------------------- Per-protein status tracking --------
pp_before = qc_before['per_protein'][['Protein', 'QC_Category_relaxed']].rename(
    columns={'QC_Category_relaxed': 'QC_Category_relaxed_before'}
)

pp_after = None
if DO_NORMALIZE_APPLIED:
    pp_after = qc_after['per_protein'][['Protein', 'QC_Category_relaxed']].rename(
        columns={'QC_Category_relaxed': 'QC_Category_relaxed_after_norm'}
    )

pp_trim = None
if TRIM_APPLIED:
    pp_trim = qc_trim['per_protein'][['Protein', 'QC_Category_relaxed']].rename(
        columns={'QC_Category_relaxed': 'QC_Category_relaxed_after_trim'}
    )

pp_master = pp_before
if pp_after is not None:
    pp_master = pp_master.merge(pp_after, on='Protein', how='outer')
if pp_trim is not None:
    pp_master = pp_master.merge(pp_trim, on='Protein', how='outer')

# Add trim flag
if TRIM_APPLIED and not removed_extremes_tbl.empty:
    rem_map = removed_extremes_tbl[['Protein', 'Remove_Protein']]
    pp_master = pp_master.merge(rem_map, on='Protein', how='left')
    pp_master['Removed_by_trim'] = pp_master['Remove_Protein'].fillna(False)
    pp_master['Removed_stage'] = np.where(pp_master['Removed_by_trim'], "after_trim", None)
    pp_master.drop(columns=['Remove_Protein'], inplace=True)
else:
    pp_master['Removed_by_trim'] = False
    pp_master['Removed_stage'] = np.nan

# -------------------- Final PCA & Heatmap ----------------
print(f"Generating Final outputs (Stage: {stage_name})...")
w_stage['Combined'] = w_stage[['Rep1', 'Rep2']].mean(axis=1)

# Coverage logic for PCA/Correlation
coverage_series = w_stage.groupby('Protein')['Combined'].apply(lambda x: x.notna().mean())
coverage_tbl = coverage_series.reset_index(name='Coverage')
kept_by_cov = coverage_tbl[coverage_tbl['Coverage'] >= COVERAGE_MIN]['Protein'].tolist()

# Prepare PCA input
pca_df = w_stage[w_stage['Protein'].isin(kept_by_cov)].copy()
X_clean, kept_pca_prots = make_pca_matrix(pca_df, 'Sample_ID', COVERAGE_MIN)

# Run PCA
pca = PCA(n_components=min(5, len(X_clean.columns), len(X_clean)-1))
# Standardize (scale=True in R)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_clean)
pca_coords = pca.fit_transform(X_scaled)
var_exp = pca.explained_variance_ratio_ * 100

pca_scores = pd.DataFrame(pca_coords, columns=[f"PC{i+1}" for i in range(pca_coords.shape[1])])
pca_scores['Sample_ID'] = X_clean.index.values

# Heatmap Logic
wide_final = w_stage[w_stage['Protein'].isin(kept_by_cov)].pivot(
    index='Sample_ID', columns='Protein', values='Combined'
)
# Median impute for correlation calc (per column)
wide_imp = wide_final.copy()
for c in wide_imp.columns:
    wide_imp[c] = wide_imp[c].fillna(wide_imp[c].median())

# Transpose for R-like cor(t(mat)) -> Cor between samples
# Spearman
cormat = wide_imp.T.corr(method='spearman')

# Inter-sample variability
tri_idx = np.tril_indices_from(cormat, k=-1)
sample_corrs = cormat.values[tri_idx]
mean_corr = np.nanmean(sample_corrs)
sd_corr = np.nanstd(sample_corrs)
cv_corr = (sd_corr / mean_corr * 100) if mean_corr != 0 else np.nan

print(f"\n=== Inter-sample correlation metrics ===")
print(f"Mean: {mean_corr:.3f} | SD: {sd_corr:.3f} | CV: {cv_corr:.1f}%")

# Save Heatmap PNG
heatmap_path = os.path.join(FINAL_DIR, f"{NAME_BASE}_SampleCorr_Heatmap_{stage_name}.png")
plt.figure(figsize=(9, 9))
sns.clustermap(cormat, cmap="RdBu_r", center=0, 
               row_cluster=True, col_cluster=True,
               dendrogram_ratio=(.1, .1), cbar_pos=(0, .02, .03, .2))
plt.savefig(heatmap_path, dpi=300)
plt.close() # clustermap creates its own figure

# -------------------- Final PDF Composition --------------
if WRITE_PDF:
    pdf_path = os.path.join(FINAL_DIR, f"{NAME_BASE}_QC_Summary_{stage_name}.pdf")
    with PdfPages(pdf_path) as pdf:
        fig = plt.figure(figsize=(13, 13))
        gs = gridspec.GridSpec(3, 2, figure=fig)
        
        # 1. Spearman (Re-plot using final stats)
        ax1 = fig.add_subplot(gs[0, 0])
        rho_vals = qc_stage_final['per_sample_rho']['Spearman_rho'].dropna()
        sns.histplot(rho_vals, bins=30, ax=ax1)
        ax1.axvline(0.8, ls="--", c='k'); ax1.axvline(0.9, ls="--", c='k')
        ax1.set_title(f"Per-sample Spearman ({stage_name})")
        
        # 2. Bland-Altman
        ax2 = fig.add_subplot(gs[0, 1])
        # Re-calc strings
        ba_b = qc_stage_final['summary']['Global_BA_bias'].iloc[0]
        ba_l = qc_stage_final['summary']['Global_BA_LOA_lower'].iloc[0]
        ba_u = qc_stage_final['summary']['Global_BA_LOA_upper'].iloc[0]
        ap = qc_stage_final['all_pairs']
        ax2.scatter(ap['Mean_log2'], ap['Diff_log2'], alpha=0.3, s=2)
        ax2.axhline(ba_b, c='r'); ax2.axhline(ba_l, ls='--', c='r'); ax2.axhline(ba_u, ls='--', c='r')
        ax2.set_title(f"BA Bias={ba_b:.3f} LoA=[{ba_l:.3f}, {ba_u:.3f}]")
        
        # 3. Density
        ax3 = fig.add_subplot(gs[1, 0])
        lp = w_stage[['Rep1', 'Rep2']].melt(var_name='variable', value_name='value').dropna()
        sns.kdeplot(data=lp, x='value', hue='variable', ax=ax3)
        ax3.set_title("Intensity Density")
        
        # 4. PCA Variance
        ax4 = fig.add_subplot(gs[1, 1])
        ax4.bar(range(1, len(var_exp)+1), var_exp)
        ax4.set_xlabel("PC"); ax4.set_ylabel("% Variance")
        ax4.set_title(f"PCA Variance ({stage_name})")
        
        # 5. PCA Scores
        ax5 = fig.add_subplot(gs[2, 0])
        ax5.scatter(pca_scores['PC1'], pca_scores['PC2'], alpha=0.8)
        ax5.set_xlabel(f"PC1 ({var_exp[0]:.1f}%)"); ax5.set_ylabel(f"PC2 ({var_exp[1]:.1f}%)")
        ax5.set_title("PCA Scores")
        
        # 6. Placeholder for Heatmap (Inserting the image is complex in matplotlib grid, 
        #    we usually just show the correlation matrix directly or skip)
        ax6 = fig.add_subplot(gs[2, 1])
        sns.heatmap(cormat, cmap="RdBu_r", center=0, cbar=False, ax=ax6)
        ax6.set_title("Sample Correlation Matrix")
        
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()
    
    print(f"PDF saved: {pdf_path}")

# -------------------- Excel Export -----------------------
# Prepare Meta table
meta_only = raw0[ID_COLS + [c for c in IGNORE_COLS if c in raw0.columns]].drop_duplicates(subset=ID_COLS[0])

# Final outputs for Excel
cor_long = w_stage.merge(coverage_tbl, on="Protein", how="left")
cor_long = cor_long.merge(qc_stage_final['per_protein'][['Protein', 'QC_Category_relaxed']], on="Protein", how="left")
cor_long['Keep_by_coverage'] = cor_long['Coverage'] >= COVERAGE_MIN

cor_matrix_df = wide_final.reset_index().merge(meta_only, on="Sample_ID", how="left")

final_xlsx = os.path.join(FINAL_DIR, f"{NAME_BASE}_FINAL_Results_Audit_{stage_name}.xlsx")

print(f"Writing Excel to {final_xlsx}...")
with pd.ExcelWriter(final_xlsx, engine='openpyxl') as writer:
    # README
    readme_data = pd.DataFrame({
        "Section": ["Purpose", "Final Stage", "Normalization", "Trimming", "Metrics"],
        "Details": [
            "Python port of DIA Replicate QC.",
            stage_name,
            f"Mode: {DECISION_NORM_MODE}, Applied: {DO_NORMALIZE_APPLIED}",
            f"Mode: {DECISION_SOURCE}, Applied: {TRIM_APPLIED}",
            f"Median Rho: {qc_stage_final['per_sample_rho']['Spearman_rho'].median():.3f}"
        ]
    })
    readme_data.to_excel(writer, sheet_name="README", index=False)
    
    # Audit Sheets
    qc_stage_final['summary'].to_excel(writer, sheet_name="QC_Summary_Final", index=False)
    qc_stage_final['per_sample_rho'].to_excel(writer, sheet_name="PerSample_Spearman_Final", index=False)
    qc_stage_final['per_protein'].to_excel(writer, sheet_name="QC_PerProtein_Final", index=False)
    coverage_tbl.to_excel(writer, sheet_name="Coverage_Final", index=False)
    pp_master.to_excel(writer, sheet_name="PerProtein_Status_All", index=False)
    
    # Data
    cor_matrix_df.to_excel(writer, sheet_name="Correlation_Matrix", index=False)
    cor_long.to_excel(writer, sheet_name="Correlation_Long", index=False)

    # Insert Image (Heatmap)
    try:
        wb = writer.book
        ws = wb.create_sheet("Heatmap")
        ws['A1'] = "Sample Correlation Heatmap (Final Stage)"
        from openpyxl.drawing.image import Image
        img = Image(heatmap_path)
        ws.add_image(img, 'A3')
    except Exception as e:
        print(f"Could not insert heatmap image into Excel: {e}")

print("\n=== SUMMARY (stage-aware) ===")
print(f"Final stage: {stage_name}")
print(f"Final outputs dir: {FINAL_DIR}")
print(f"Preferred file: {os.path.basename(final_xlsx)}")