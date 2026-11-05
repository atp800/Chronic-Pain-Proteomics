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
from limma_reporting import RESULT_COLUMNS, attach_support, protein_support, write_results

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
    excluded_id_rows = 0
    if is_paired:
        sample_ids = df_subset[config["ID_COLUMN"]]
        subjects = pd.Series(
            derive_subject_ids(sample_ids, delimiter=config["ID_DELIMITER"]),
            index=df_subset.index, dtype=object
        )
        known = sample_ids.notna() & subjects.str.strip().ne("")
        if not known.all():
            excluded_id_rows = int((~known).sum())
            print(f"Excluding {int((~known).sum())} rows without a participant ID.")
            df_subset = df_subset.loc[known]
            subjects = subjects.loc[known]

    groups = df_subset[group_col].astype(str)
    baseline, comparison = str(baseline_name), str(compare_name)
    if baseline == comparison or set(groups) != {baseline, comparison}:
        print("Skipping limma: both distinct levels of the requested comparison are required.")
        return

    observed = df_subset[protein_cols].replace([np.inf, -np.inf], np.nan)
    X_limma = observed.copy()
    complete = None
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
    # IDs are used for reporting distinct participants in independent comparisons;
    # delta rows already contain participant IDs, so do not split them a second time.
    report_subjects = subjects
    if report_subjects is None:
        id_col = config.get('ID_COLUMN')
        if id_col in df_subset:
            report_subjects = df_subset[id_col].astype('string')
            if not is_delta:
                report_subjects = report_subjects.str.split(config.get('ID_DELIMITER', '-'), n=1, regex=False).str[0]
            report_subjects = report_subjects.mask(report_subjects.str.strip().eq(''))
        else:
            report_subjects = pd.Series(pd.NA, index=df_subset.index, dtype='string')
    support = protein_support(observed, X_limma, groups, report_subjects, baseline, comparison,
                              replicated, keep, paired=is_paired, delta=is_delta,
                              complete_pairs=None if complete is None else complete.sum())
    info = pd.DataFrame([
        {'Setting': 'Comparison', 'Value': f'{comparison} - {baseline}'},
        {'Setting': 'Mode', 'Value': 'Delta group comparison' if is_delta else ('Paired longitudinal' if is_paired else 'Independent group comparison')},
        {'Setting': 'Observed counts', 'Value': 'Finite protein values before eligibility filtering and before masking unpaired measurements; zeros count as observed. Delta counts describe observed changes, not individual endpoints.'},
        {'Setting': 'Participant counts', 'Value': 'Distinct IDs using the configured delimiter; delta IDs are already participant IDs. Blank when an observed value has an unknown ID. Counts describe the supplied input, not upstream imputation.'},
        {'Setting': 'Fit observation counts', 'Value': 'Observed rows supplied to lmFit for this protein, after masking unpaired measurements. Zero for proteins excluded before fitting.'},
        {'Setting': 'Complete pairs', 'Value': 'Longitudinal: participants observed at both compared visits. Delta: participants with a finite change (both endpoints observed in the main pipeline). Not applicable to single-visit group comparisons.'},
        {'Setting': 'Missing proportions', 'Value': 'Missing/non-finite protein values divided by input rows in each group, before masking unpaired measurements; proportions run from 0 to 1.'},
        {'Setting': 'Protein imputation', 'Value': 'None in limma; protein imputed proportions are zero. Missing covariates are handled separately by the existing mean/Unknown replacement.'},
        {'Setting': 'Excluded rows with unknown participant ID', 'Value': excluded_id_rows},
        {'Setting': 'Residual_DF', 'Value': 'Per-protein residual degrees of freedom from lmFit, including the fitted covariates/participant effects. Blank if not fitted.'},
        {'Setting': 'Moderated_DF', 'Value': 'Per-protein total degrees of freedom from eBayes used for moderated t-tests and confidence intervals. Blank if not tested.'},
        {'Setting': 'CI.L / CI.R', 'Value': 'Individual 95% moderated confidence limits from limma topTable for the reported logFC coefficient. Not simultaneous or selection-adjusted intervals across proteins; retain adj.P.Val for multiple-testing assessment.'},
        {'Setting': 'Effect scale', 'Value': 'Comparison minus baseline on the supplied scale. With log2 input: log2 fold change; in delta mode: between-group difference in log2 changes. Covariate sheets show that covariate coefficient.'},
        {'Setting': 'Covariates', 'Value': ', '.join(config.get('COVARIATE_COLUMNS', [])) or 'None'},
    ])
    file_prefix = f"Limma_{baseline_name}_vs_{compare_name}"
    excel_path = os.path.join(output_dir, f"{file_prefix}_Results.xlsx")
    empty_results = attach_support(pd.DataFrame(columns=RESULT_COLUMNS), support)
    X_limma = X_limma.loc[:, keep]
    print(
        f"Retained {X_limma.shape[1]}/{len(protein_cols)} proteins "
        f"with at least {requirement} and observed variation."
    )
    if X_limma.empty:
        write_results(excel_path, empty_results, {}, support, info, config)
        print("Skipping limma: no eligible proteins. Support and exclusion reasons saved; result sheets are empty.")
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
    with localconverter(ro.default_converter + pandas2ri.converter):
        fit_support = ro.conversion.rpy2py(ro.r('data.frame(Protein=fit$genes$Protein, Residual_DF=fit$df.residual, Effect_Estimable=is.finite(fit$coefficients[, "comparison"]), Testable=testable)'))
    support = support.set_index('Protein')
    fit_support = fit_support.set_index('Protein')
    support.loc[fit_support.index, 'Residual_DF'] = fit_support['Residual_DF']
    support.loc[fit_support.index, 'Status'] = np.where(
        ~fit_support['Effect_Estimable'], 'Comparison effect is not estimable',
        np.where(fit_support['Residual_DF'].le(0), 'No residual degrees of freedom', 'Tested'))
    support = support.reset_index()
    print(
        f"Testable proteins for {compare_string}: {n_testable}/{X_limma.shape[1]} "
        "(estimable effect and positive residual degrees of freedom)."
    )
    if n_testable == 0:
        write_results(excel_path, attach_support(pd.DataFrame(columns=RESULT_COLUMNS), support), {}, support, info, config)
        print("Skipping limma: no testable proteins. Support and exclusion reasons saved; result sheets are empty.")
        return

    ro.r('fit <- eBayes(fit[testable, ])')
    ro.r('fit$genes$Residual_DF <- fit$df.residual; fit$genes$Moderated_DF <- fit$df.total')
    r_results = ro.r('topTable(fit, coef="comparison", number=Inf, sort.by="P", confint=0.95)')
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_results = ro.conversion.rpy2py(r_results)
    moderated_df = df_results.set_index('Protein')['Moderated_DF']
    support['Moderated_DF'] = support['Protein'].map(moderated_df)
    df_results = attach_support(df_results, support)
    
    # Extract any covariate results from the same set of tested proteins.
    cov_results_dict = {}
    
    # Get any columns in the design matrix that start with "cov_"
    cov_cols = list(ro.r('grep("^cov_", colnames(design), value=TRUE)'))
    
    for c_col in cov_cols:
        r_cov_res = ro.r(f'topTable(fit, coef="{c_col}", number=Inf, sort.by="P", confint=0.95)')
        with localconverter(ro.default_converter + pandas2ri.converter):
            cov_results_dict[c_col] = attach_support(ro.conversion.rpy2py(r_cov_res), support)


    # Output generation
    
    # 1. Volcano Plot
    pdf_path = os.path.join(output_dir, f"{file_prefix}_Volcano.pdf")
    with PdfPages(pdf_path) as pdf:
        plt.figure(figsize=(10, 8))
        sns.scatterplot(data=df_results, x='logFC', y='adj.P.Val', alpha=0.5, edgecolor=None)
        plt.yscale('log')
        plt.gca().invert_yaxis()
        plt.axhline(0.05, color='r', linestyle='--', alpha=0.5, label='p=0.05')
        plt.axvline(config["LOG_FC_THRESHOLD"], color='b', linestyle='--', alpha=0.5, label=f'|logFC| = {config["LOG_FC_THRESHOLD"]}')
        plt.axvline(-config["LOG_FC_THRESHOLD"], color='b', linestyle='--', alpha=0.5)
        plt.title(f"Volcano Plot: {compare_string}")
        plt.xlabel("Log2 Fold Change")
        plt.ylabel("Adjusted P-Value (Log Scale)")
        plt.legend()
        pdf.savefig()
        plt.close()

    # 2. Excel File
    write_results(excel_path, df_results, cov_results_dict, support, info, config)
