import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import glob
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.feature_selection import VarianceThreshold
from helpers import impute_missing_vals, derive_subject_ids

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
    X_limma = df_subset[protein_cols]
    
    # Remove zero-variance features
    print("Removing zero-variance proteins")
    X_temp_for_var = X_limma.fillna(0)                                    # Temporarily fill NaNs for variance calculation
    selector = VarianceThreshold(threshold=0)                             # Increase threshold to 0.1 to remove low variance features
    try:
        selector.fit(X_temp_for_var)
        cols_kept = X_limma.columns[selector.get_support()]
        X_limma = X_limma[cols_kept]
        print(f"Proteins after dropping low variance: {X_limma.shape[1]}")
    except ValueError:
        print("Warning: Variance filter failed (variance too low globally?), skipping.")


    # Impute missing values (unless already done in Delta)
    if not is_delta and X_limma.isna().sum().sum() > 0:
        print("Imputing missing values...")
        X_limma = impute_missing_vals(X_limma, shift=1.8, width=0.3)


    #########################################
    print("\n--- DATA SCALE DIAGNOSTICS ---")
    data_max = X_limma.max().max()
    data_min = X_limma.min().min()
    print(f"Data Max: {data_max:.2f}, Data Min: {data_min:.2f}")
    
    if data_max > 100:
        print("!!! WARNING: Data magnitude suggests RAW INTENSITY, not Log2.")
        print("Limma results will be invalid unless you Log2 transform.")
########################################

    # Need to transpose: Rows = Proteins, Columns = Samples
    r_expression_matrix = X_limma.T

    # Prep groups
    clean_groups =[str(g).replace(" ", "_").replace("-", "_") for g in df_subset[group_col]]
    r_groups = ro.StrVector(clean_groups)
    
    # Prepare IDs
    raw_sample_ids = df_subset[config["ID_COLUMN"]].tolist()

    if is_paired:
        derived_subjects = derive_subject_ids(raw_sample_ids, delimiter=config["ID_DELIMITER"])
        r_subjects = ro.StrVector(derived_subjects)
    else:
        r_subjects = ro.StrVector(raw_sample_ids)


    # ######### CREATE DESIGN MATRIX XXXXXXXXXXXX
    ro.globalenv['expression_data'] = r_expression_matrix
    ro.globalenv['groups'] = r_groups
    ro.globalenv['subjects'] = r_subjects
    ro.r('groups <- as.factor(groups)')

    if is_paired:
        print("Creating Paired Design Matrix (~ 0 + groups + subjects)...")
        
        if len(set(derived_subjects)) == len(derived_subjects):                     # Validation: Check if pairing is actually possible
             print("CRITICAL WARNING: Derived Subject IDs are all unique.")
             print("Limma cannot pair samples if every ID is different.")
             ro.r('design <- model.matrix(~ 0 + groups)')
        else:
            ro.r('subjects <- as.factor(subjects)')                                 # force R to treat subjects as a factor
            ro.r('design <- model.matrix(~ 0 + groups + subjects)')
    else:
        print("Creating Independent Design Matrix...")
        ro.r('design <- model.matrix(~ 0 + groups)')
        



# 1. Clean the baseline and compare strings for the contrast equation
    clean_base = str(baseline_name).replace(" ", "_").replace("-", "_")
    clean_compare = str(compare_name).replace(" ", "_").replace("-", "_")
    
    # 2. Get column names directly from the R object 'design'
    col_names = list(ro.r('colnames(design)'))
    
    # Create unique clean groups (Assumes 'clean_groups' list is defined earlier in your script)
    unique_groups_clean = sorted(list(set(clean_groups)))
    new_col_names =[]

    # 3. Robust column renaming logic
    for c in col_names:
        found_match = False
        for g in unique_groups_clean:
            # Strictly check if the column is exactly "groups" + the group name
            # This prevents substring errors (e.g. "Resp" matching inside "groupsNon_Resp")
            expected_col_name = f"groups{g}"
            if c == expected_col_name:
                new_col_names.append(g)
                found_match = True
                break
        if not found_match:
            new_col_names.append(c)
            
    # Assign new names back to the R object 'design'
    ro.globalenv['new_cols'] = ro.StrVector(new_col_names)
    ro.r('colnames(design) <- new_cols')

    # 4. Fit Linear Model
    ro.r('expression_data <- as.matrix(expression_data)')
    ro.r('fit <- lmFit(expression_data, design)')
    
    # 5. Define contrasts explicitly using baseline and compare names
    compare_string = f"{clean_compare} - {clean_base}"
    print(f"Defining Contrast: {compare_string}")
    
    ro.r(f'contrast.matrix <- makeContrasts("{compare_string}", levels=design)')
    ro.r('fit2 <- contrasts.fit(fit, contrast.matrix)')
    ro.r('fit2 <- eBayes(fit2)')
    
    # 6. Extract results back to Python/Pandas
    r_results = ro.r('topTable(fit2, number=Inf, sort.by="P")')
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_results = ro.conversion.rpy2py(r_results)
    
    df_results.index.name = "Protein"
    df_results.reset_index(inplace=True)

    

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

    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        df_results[final_cols].sort_values(by='adj.P.Val').to_excel(writer, sheet_name='All_Proteins', index=False)
        df_sig.sort_values(by='P.Value').to_excel(writer, sheet_name='Significant_Raw_PVal', index=False)
        df_sig_adj.sort_values(by='adj.P.Val').to_excel(writer, sheet_name='Significant_Adj_PVal', index=False)