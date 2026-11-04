import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import glob
from matplotlib.backends.backend_pdf import PdfPages
# from sklearn.feature_selection import VarianceThreshold
# from helpers import impute_missing_vals, derive_subject_ids
from helpers import derive_subject_ids

##################################################################
# SET LOCATION OF R INSTALLATION
##################################################################

# OPTION A: Manually to your R installation path here
r_home_path = None                                                      # Example: r_home_path = r"C:\Program Files\R\R-4.3.1"

# OPTION B: Auto-detect standard Windows R location if not set above
if not r_home_path:
    try:
        possible_paths = glob.glob(r"C:\Program Files\R\R-*")
        if possible_paths:
            # Pick the highest version number
            possible_paths.sort(reverse=True)
            r_home_path = possible_paths[0]
            print(f"Auto-detected R installation: {r_home_path}")
    except Exception:
        pass

if r_home_path and os.path.exists(r_home_path):
    # Set R_HOME environment variable
    os.environ['R_HOME'] = r_home_path
    
    # Set R_USER to current user to avoid permission errors
    os.environ['R_USER'] = os.path.expanduser("~")
    
    # Add R architecture-specific bin to PATH (Required for DLL loading)
    # Checks for x64 (64-bit) or i386 (32-bit)
    r_bin_x64 = os.path.join(r_home_path, 'bin', 'x64')
    r_bin_root = os.path.join(r_home_path, 'bin')
    
    if os.path.exists(r_bin_x64):
        os.environ['PATH'] = r_bin_x64 + ";" + os.environ['PATH']
    elif os.path.exists(r_bin_root):
        os.environ['PATH'] = r_bin_root + ";" + os.environ['PATH']
    else:
        print("Warning: Could not find R/bin folder. rpy2 might fail.")
else:
    print("Warning: Could not determine R_HOME. If rpy2 fails, set 'r_home_path' manually in the script.")


# R Imports
import rpy2
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.packages import importr

##########################################################################
# Version details for reproducibility/citations
print("\n--- R and Limma Version Information ---")

r_version = ro.r('R.version.string')[0]
print(f"R version: {r_version}")

limma_version = ro.r('as.character(packageVersion("limma"))')[0]
print(f"limma version: {limma_version}")

session_info = "\n".join(ro.r('capture.output(sessionInfo())'))
print(f"Session Info:\n{session_info}")

rpy2_version = rpy2.__version__
print(f"rpy2 version: {rpy2_version}")

##########################################################################

def run_limma(df_subset, group_col, protein_cols, baseline_name, compare_name, output_dir, is_paired, is_delta, config):

    print("\n\n------------------------------------------------------------------")
    print(f"Running Limma: {baseline_name} vs {compare_name}")
    print("Mode: PAIRED Analysis (grouping by subject ID)" if is_paired else "Mode: INDEPENDENT Analysis")
    

    # R-Python DataFrame Conversion    
    pandas2ri.activate()
    limma = importr('limma')
    base = importr('base')
    stats = importr('stats')


    ######### PREPARE DATA ###########
    subjects = None
    if is_paired:
        sample_ids = df_subset[config["ID_COLUMN"]]
        subjects = pd.Series(
            derive_subject_ids(sample_ids, delimiter=config["ID_DELIMITER"]),
            index=df_subset.index, dtype=object
        )
        known = sample_ids.notna() & subjects.str.strip().ne("")
        if not known.all():
            print(f"Excluding {int((~known).sum())} rows without a participant ID.")
            df_subset = df_subset.loc[known]
            subjects = subjects.loc[known]

    groups = df_subset[group_col].astype(str)
    baseline, comparison = str(baseline_name), str(compare_name)
    if baseline == comparison or set(groups) != {baseline, comparison}:
        print("Skipping limma: both distinct levels of the requested comparison are required.")
        return

    X_limma = df_subset[protein_cols].replace([np.inf, -np.inf], np.nan)
    if is_paired:
        if pd.MultiIndex.from_arrays([subjects, groups]).has_duplicates:
            print("Skipping paired limma: duplicate participant/visit. Check IDs or resolve replicate samples.")
            return
        # With one row per visit, two observations mean this person has both
        # measurements for this protein. Mask singletons before fitting.
        complete = X_limma.groupby(subjects).count().eq(2)
        replicated = complete.sum().ge(2)
        X_limma = X_limma.where(complete.loc[subjects].set_axis(X_limma.index))
        requirement = "two complete participant pairs"
    else:
        replicated = X_limma.groupby(groups).count().ge(2).all()
        requirement = "two observations in each group"

    keep = replicated & X_limma.var(ddof=0).gt(0)
    X_limma = X_limma.loc[:, keep]
    print(
        f"Retained {X_limma.shape[1]}/{len(protein_cols)} proteins "
        f"with at least {requirement} and observed variation."
    )
    if X_limma.empty:
        print("Skipping limma: no eligible proteins. No new limma results were written.")
        return

    print(
        f"Missing protein values retained: "
        f"{int(X_limma.isna().sum().sum())}"
    )


    #########################################
    print("\n--- DATA SCALE DIAGNOSTICS ---")
    data_max = X_limma.max().max()
    data_min = X_limma.min().min()
    print(f"Data Max: {data_max:.2f}, Data Min: {data_min:.2f}")
    
    if data_max > 100:
        print("!!! WARNING: Data magnitude suggests RAW INTENSITY, not Log2.")
        print("Limma results will be invalid unless you Log2 transform.")
########################################

    # Rows in R are proteins. Design vectors follow the same sample order.
    ro.globalenv['expression_data'] = X_limma.T
    # Carry IDs as limma annotations: lmFit can drop row names for a single protein.
    ro.globalenv['protein_names'] = ro.StrVector(X_limma.columns.astype(str).tolist())
    ro.globalenv['comparison'] = ro.IntVector(groups.eq(comparison).astype(int).tolist())
    formula_terms = ["1"]
    if is_paired:
        ro.globalenv['subjects'] = ro.StrVector(subjects.tolist())
        ro.r('subjects <- as.factor(subjects)')
        formula_terms.append("subjects")
        if config.get("COVARIATE_COLUMNS"):
            print("Participant-constant covariates cannot be estimated separately from participant effects.")

    # Add covariates
    covariates = config.get("COVARIATE_COLUMNS", [])
    for cov in covariates:
        # Clean variable name for R
        clean_cov_name = "cov_" + str(cov).replace(" ", "_").replace("-", "_").replace("(", "").replace(")", "")
        
        # Check if the covariate is numeric (e.g. age) or categorical (e.g. sex)
        if pd.api.types.is_numeric_dtype(df_subset[cov]):
            # Fill missing numeric values with mean (Limma fails if design matrix has NAs)
            vals = df_subset[cov].fillna(df_subset[cov].mean()).values
            ro.globalenv[clean_cov_name] = ro.FloatVector(vals)
            print(f"Added numeric covariate: {cov}")
        else:
            # Fill missing categorical values with "Unknown" and convert to R Factor
            vals = df_subset[cov].fillna("Unknown").astype(str).values
            ro.globalenv[clean_cov_name] = ro.StrVector(vals)
            ro.r(f'{clean_cov_name} <- as.factor({clean_cov_name})')
            print(f"Added categorical covariate: {cov}")
            
        formula_terms.append(clean_cov_name)

    # Enter the target effect after all nuisance terms: lmFit then marks it
    # non-estimable if the observed participant/covariate columns explain it.
    formula_terms.append("comparison")
    formula_str = " + ".join(formula_terms)
    print(f"Creating Design Matrix (~ {formula_str})...")
    ro.r(f'design <- model.matrix(~ {formula_str})')
    compare_string = f"{comparison} - {baseline}"
    n_testable = int(ro.r("""
        fit <- lmFit(as.matrix(expression_data), design)
        fit$genes <- data.frame(Protein=protein_names)
        testable <- fit$df.residual > 0 & is.finite(fit$coefficients[, "comparison"])
        sum(testable)
    """)[0])
    print(
        f"Testable proteins for {compare_string}: {n_testable}/{X_limma.shape[1]} "
        "(estimable effect and positive residual degrees of freedom)."
    )
    if n_testable == 0:
        print("Skipping limma: no testable proteins. No new limma results were written.")
        return

    ro.r('fit <- eBayes(fit[testable, ])')
    r_results = ro.r('topTable(fit, coef="comparison", number=Inf, sort.by="P")')
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_results = ro.conversion.rpy2py(r_results)
    
    # Extract any covariate results from the same set of tested proteins.
    cov_results_dict = {}
    
    # Get any columns in the design matrix that start with "cov_"
    cov_cols = list(ro.r('grep("^cov_", colnames(design), value=TRUE)'))
    
    for c_col in cov_cols:
        r_cov_res = ro.r(f'topTable(fit, coef="{c_col}", number=Inf, sort.by="P")')
        with localconverter(ro.default_converter + pandas2ri.converter):
            cov_results_dict[c_col] = ro.conversion.rpy2py(r_cov_res)


    # Output generation
    file_prefix = f"Limma_{baseline_name}_vs_{compare_name}"
    
    # 1. Volcano Plot
    pdf_path = os.path.join(output_dir, f"{file_prefix}_Volcano.pdf")
    with PdfPages(pdf_path) as pdf:
        plt.figure(figsize=(10, 8))
        sns.scatterplot(data=df_results, x='logFC', y='adj.P.Val', alpha=0.5, edgecolor=None)
        plt.yscale('log')
        plt.gca().invert_yaxis()
        plt.axhline(0.05, color='r', linestyle='--', alpha=0.5, label='p=0.05')
        plt.axvline(config["LOG_FC_THRESHOLD"], color='b', linestyle='--', alpha=0.5)
        plt.axvline(-config["LOG_FC_THRESHOLD"], color='b', linestyle='--', alpha=0.5)
        plt.title(f"Volcano Plot: {compare_string}")
        plt.xlabel("Log2 Fold Change")
        plt.ylabel("Adjusted P-Value (Log Scale)")
        plt.legend()
        pdf.savefig()
        plt.close()

    # 2. Excel File
    excel_path = os.path.join(output_dir, f"{file_prefix}_Results.xlsx")
    final_cols = [c for c in['Protein', 'logFC', 'P.Value', 'adj.P.Val', 'AveExpr', 'B', 't'] if c in df_results.columns]
    
    sig_mask = (df_results['P.Value'] < config["P_THRESHOLD"]) & (df_results['logFC'].abs() > config["LOG_FC_THRESHOLD"])
    df_sig = df_results.loc[sig_mask, final_cols]

    sig_mask_adj = (df_results['adj.P.Val'] < config["P_THRESHOLD"]) & (df_results['logFC'].abs() > config["LOG_FC_THRESHOLD"])
    df_sig_adj = df_results.loc[sig_mask_adj, final_cols]

    print(f"Significant proteins (raw p-value < {config['P_THRESHOLD']} & abs(logFC) > {config['LOG_FC_THRESHOLD']}): {df_sig.shape[0]}")
    print(f"Significant proteins (adjusted p-value < {config['P_THRESHOLD']} & abs(logFC) > {config['LOG_FC_THRESHOLD']}): {df_sig_adj.shape[0]}")

    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        df_results[final_cols].sort_values(by='adj.P.Val').to_excel(writer, sheet_name='All_Proteins', index=False)
        df_sig.sort_values(by='P.Value').to_excel(writer, sheet_name='Significant_Raw_PVal', index=False)
        df_sig_adj.sort_values(by='adj.P.Val').to_excel(writer, sheet_name='Significant_Adj_PVal', index=False)

        # Add a new tab for each covariate
        for cov_name, df_cov in cov_results_dict.items():
            sheet_name = f"Cov_Effect_{cov_name.replace('cov_', '')}"[:31]  # truncate if over 31 char excel sheet name limit
            c_cols = [c for c in ['Protein', 'logFC', 'P.Value', 'adj.P.Val', 'AveExpr', 't'] if c in df_cov.columns]
            df_cov[c_cols].to_excel(writer, sheet_name=sheet_name, index=False)
