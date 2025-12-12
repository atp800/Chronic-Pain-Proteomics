import os
import sys
import glob
import datetime
import json
import hashlib
import warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
import openpyxl
from openpyxl.utils import get_column_letter

# Logsitic Regression (LASSO)
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import roc_curve, auc

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# R
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter



filepath = os.path.abspath(__file__)







##################################################################
###### Set lcoation of R installation for running subscript ######
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

##################################################################








# Activate pandas to R conversion
pandas2ri.activate()
limma = importr('limma')
base = importr('base')

# ==========================================
# 1. LOAD AND PREPARE DATA
# ==========================================

def load_and_clean_data(filepath):
    """
    Loads data, handles 0s, log-transforms, normalizes, and imputes.
    Assumes structure: Subject, Time, Group, PPT, Replicate, Protein1, Protein2...
    """
    df = pd.read_csv(filepath)

    # separate metadata and protein data
    metadata_cols = ['Subject_ID', 'Timepoint', 'Group', 'PPT', 'Replicate']
    protein_data = df.drop(columns=metadata_cols)
    metadata = df[metadata_cols]

    # A. Zero Handling: Convert 0 in proteins to NaN (LOD), but keep 0 in PPT
    protein_data = protein_data.replace(0, np.nan)

    # B. Log2 Transformation
    protein_log = np.log2(protein_data)

    # C. Normalization: Sample Median Centering
    # Subtract the median of each sample from all proteins in that sample
    sample_medians = protein_log.median(axis=1)
    protein_norm = protein_log.sub(sample_medians, axis=0)

    # D. Imputation (MinProb / 1st Percentile Strategy)
    # If a protein is missing > 50% overall, drop it
    valid_proteins = protein_norm.columns[protein_norm.isnull().mean() < 0.5]
    protein_norm = protein_norm[valid_proteins]
    
    # Fill remaining NaNs with the 1st percentile of that specific protein (simulating LOD)
    for col in protein_norm.columns:
        min_val = np.nanpercentile(protein_norm[col], 1)
        protein_norm[col] = protein_norm[col].fillna(min_val)

    # Recombine
    df_clean = pd.concat([metadata, protein_norm], axis=1)
    
    return df_clean

def handle_replicates(df):
    """
    Checks PCA for batch effects, then averages Replicates A and B.
    """
    # Quick PCA check
    features = df.drop(columns=['Subject_ID', 'Timepoint', 'Group', 'PPT', 'Replicate'])
    pca = PCA(n_components=2)
    components = pca.fit_transform(StandardScaler().fit_transform(features))
    
    # Plot Replicate Batch Effect
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=components[:,0], y=components[:,1], hue=df['Replicate'])
    plt.title("QC: PCA by Replicate (Check for Batch Effect)")
    plt.show()

    # Average Replicates (A & B) by Subject and Timepoint
    # Numeric columns (Proteins + PPT) get averaged
    df_averaged = df.groupby(['Subject_ID', 'Timepoint', 'Group']).mean(numeric_only=True).reset_index()
    
    return df_averaged

# ==========================================
# 2. STATISTICAL INFERENCE (LIMMA via RPY2)
# ==========================================

def run_limma_paired(df, target_group="Responder"):
    """
    Obj 3: Compare D0 vs D14 in Responders (Paired by Subject).
    """
    # Subset: Specific Group, D0 and D14 only
    subset = df[(df['Group'] == target_group) & (df['Timepoint'].isin(['D0', 'D14']))].copy()
    subset = subset.sort_values(['Subject_ID', 'Timepoint']) # Ensure sorting for blocking

    # Prepare R variables
    # Transpose: Limma expects Rows=Proteins, Cols=Samples
    protein_matrix = subset.drop(columns=['Subject_ID', 'Timepoint', 'Group', 'PPT']).T
    
    # Create Design Matrix in R
    # Formula: ~ 0 + Timepoint + Subject_ID
    subjects = ro.StrVector(subset['Subject_ID'].astype(str).tolist())
    times = ro.StrVector(subset['Timepoint'].astype(str).tolist())
    
    design = ro.r['model.matrix'](ro.Formula('~ 0 + times + subjects'))
    colnames = ro.r['colnames'](design)
    
    # Fit Linear Model
    fit = limma.lmFit(protein_matrix, design)
    
    # Define Contrast (D14 - D0)
    # Note: Check print(colnames) to ensure exact string match. Usually 'timesD14' and 'timesD0'
    contrast_matrix = limma.makeContrasts(Diff = "timesD14 - timesD0", levels=design)
    
    fit2 = limma.contrasts_fit(fit, contrast_matrix)
    fit2 = limma.eBayes(fit2)
    
    # Get Results
    results = limma.topTable(fit2, coef="Diff", number=float('inf'), sort_by="P")
    
    # Convert back to Pandas
    with localconverter(ro.default_converter + pandas2ri.converter):
        pd_results = ro.conversion.rpy2py(results)
    
    pd_results['Gene'] = results.rownames
    return pd_results

def run_limma_groups(df, timepoint="D14"):
    """
    Obj 4: Compare Responder vs Non-Responder at D14.
    """
    subset = df[df['Timepoint'] == timepoint].copy()
    
    protein_matrix = subset.drop(columns=['Subject_ID', 'Timepoint', 'Group', 'PPT']).T
    groups = ro.StrVector(subset['Group'].tolist())
    
    design = ro.r['model.matrix'](ro.Formula('~ 0 + groups'))
    
    fit = limma.lmFit(protein_matrix, design)
    
    # Contrast: Responder - NonResponder
    # Note: R typically names these 'groupsResponder' and 'groupsNonResponder'
    contrast_matrix = limma.makeContrasts(Diff = "groupsResponder - groupsNonResponder", levels=design)
    
    fit2 = limma.contrasts_fit(fit, contrast_matrix)
    fit2 = limma.eBayes(fit2)
    
    results = limma.topTable(fit2, coef="Diff", number=float('inf'), sort_by="P")
    
    with localconverter(ro.default_converter + pandas2ri.converter):
        pd_results = ro.conversion.rpy2py(results)
        
    pd_results['Gene'] = results.rownames
    return pd_results

# ==========================================
# 3. PREDICTIVE MODELING (LASSO)
# ==========================================

def run_lasso_classifier(df, timepoint="D14"):
    """
    Obj 4/Bonus: Find biomarker panel to classify Responder vs Non-Responder.
    """
    subset = df[df['Timepoint'] == timepoint].copy()
    
    X = subset.drop(columns=['Subject_ID', 'Timepoint', 'Group', 'PPT'])
    y = subset['Group'].apply(lambda x: 1 if x == 'Responder' else 0)
    
    # Standardize features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Logistic Regression with L1 Penalty (LASSO)
    # CV=5 automatically finds best regularization strength (C)
    clf = LogisticRegressionCV(cv=5, penalty='l1', solver='liblinear', random_state=42, max_iter=10000)
    clf.fit(X_scaled, y)
    
    # Extract coefficients
    coefs = pd.Series(clf.coef_[0], index=X.columns)
    important_features = coefs[coefs != 0].sort_values(ascending=False)
    
    # Plot ROC
    y_pred_prob = clf.predict_proba(X_scaled)[:, 1]
    fpr, tpr, _ = roc_curve(y, y_pred_prob)
    roc_auc = auc(fpr, tpr)
    
    plt.figure(figsize=(6, 6))
    plt.plot(fpr, tpr, label=f'AUC = {roc_auc:.2f}')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.title(f'LASSO ROC Curve ({timepoint})')
    plt.legend()
    plt.show()
    
    return important_features

# ==========================================
# 4. CLINICAL CORRELATION (SPEARMAN)
# ==========================================

def run_correlation(df, timepoint="D14"):
    """
    Obj 6: Correlate proteins with %PPT score.
    """
    subset = df[df['Timepoint'] == timepoint].copy()
    proteins = subset.drop(columns=['Subject_ID', 'Timepoint', 'Group', 'PPT'])
    scores = subset['PPT']
    
    results = []
    for protein in proteins.columns:
        rho, pval = spearmanr(proteins[protein], scores, nan_policy='omit')
        if abs(rho) > 0.5 and pval < 0.05: # Filter for interesting hits
            results.append({'Protein': protein, 'Rho': rho, 'P_Value': pval})
            
    return pd.DataFrame(results).sort_values('P_Value')

# ==========================================
# 5. TRAJECTORY PLOTTING
# ==========================================

def plot_trajectory(df, proteins_of_interest):
    """
    Obj 5: Visualize the 'Ramp Up' of top hits over time.
    """
    # Define time order
    time_order = ['D0', 'D3', 'D7', 'D9', 'D11', 'D14', 'D25']
    
    # Melt dataframe for seaborn
    melted = df.melt(id_vars=['Subject_ID', 'Timepoint', 'Group'], 
                     value_vars=proteins_of_interest, 
                     var_name='Protein', value_name='Abundance')
    
    plt.figure(figsize=(12, 6))
    sns.lineplot(data=melted, x='Timepoint', y='Abundance', hue='Group', 
                 style='Protein', markers=True, dashes=False, 
                 hue_order=['Non-Responder', 'Responder'], 
                 order=time_order) # Explicit time sort
    plt.title("Trajectory of Top Biomarkers")
    plt.grid(True, linestyle='--')
    plt.show()

# ==========================================
# MAIN EXECUTION
# ==========================================

if __name__ == "__main__":
    # 1. PREP
    # df_clean = load_and_clean_data('your_data.csv')
    # df_final = handle_replicates(df_clean)
    
    # --- MOCK DATA GENERATION FOR DEMO PURPOSES ---
    # (Generating dummy data so you can run this script immediately to see logic)
    print("Generating mock data...")
    cols = ['Subject_ID', 'Timepoint', 'Group', 'PPT', 'Replicate'] + [f'Prot_{i}' for i in range(1, 101)]
    rows = []
    for subj in range(1, 11): # 10 subjects
        grp = 'Responder' if subj <= 5 else 'Non-Responder'
        for time in ['D0', 'D14', 'D25']:
            for rep in ['A', 'B']:
                row = [subj, time, grp, np.random.randint(10, 90), rep]
                row += np.random.normal(10, 2, 100).tolist() # Random protein data
                rows.append(row)
    df_raw = pd.DataFrame(rows, columns=cols)
    # ----------------------------------------------

    # 1. Clean & Average Replicates
    # Skipping load_and_clean_data for mock, going straight to averaging
    df_final = handle_replicates(df_raw) 

    # 2. Run Limma (Time: D0 vs D14)
    print("\nRunning Limma (D0 vs D14 Responders)...")
    res_time = run_limma_paired(df_final, target_group="Responder")
    top_time_hits = res_time[res_time['adj.P.Val'] < 0.05]['Gene'].tolist()
    print(f"Found {len(top_time_hits)} proteins changing over time.")
    print(res_time.head(3))

    # 3. Run Limma (Group: Resp vs NonResp)
    print("\nRunning Limma (Responder vs Non-Responder at D14)...")
    res_group = run_limma_groups(df_final, timepoint="D14")
    top_group_hits = res_group[res_group['adj.P.Val'] < 0.05]['Gene'].tolist()
    print(f"Found {len(top_group_hits)} proteins differing between groups.")
    print(res_group.head(3))

    # 4. Run LASSO
    print("\nRunning LASSO Classifier...")
    lasso_hits = run_lasso_classifier(df_final, timepoint="D14")
    print("Top Predictive Proteins:")
    print(lasso_hits)

    # 5. Run Correlation
    print("\nRunning Spearman Correlation with %PPT...")
    corr_results = run_correlation(df_final, timepoint="D14")
    print(corr_results.head())

    # 6. Plot Trajectory of top hit (if any exist, otherwise random)
    targets = (top_time_hits + top_group_hits)[:3] 
    if not targets: targets = ['Prot_1', 'Prot_2'] # Fallback for mock data
    
    print("\nPlotting trajectories...")
    plot_trajectory(df_final, targets)

