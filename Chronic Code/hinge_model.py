import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.regression.mixed_linear_model import MixedLM
from itertools import product
import os


'''
Creates a mixed-effects hinge model for each protein, trying multiple candidate peak days and selecting the best fit based on AIC
Runs across every 

'''

# ---------------------------
# USER PARAMETERS
# ---------------------------
INPUT_PATH = "INPUT_PATH"   # path to input CSV
SHEET_NAME = "SHEET_NAME"   # sheet name in Excel file
OUTPUT_PATH = "OUTPUT_PATH" # path to save outputs
PEAK_CANDIDATES = [11, 14, 17]  # candidate peak days
TOP_N = 5  # number of best/worst proteins to plot

# ---------------------------
# HELPER FUNCTIONS
# ---------------------------

def create_hinge_terms(df, peak_day):
    """Create hinge terms for the mixed effect model"""
    df['t_minus'] = np.minimum(0, df['Time'] - peak_day)
    df['t_plus'] = np.maximum(0, df['Time'] - peak_day)
    return df

def fit_hinge_mixed_model(df, protein_col, peak_day):
    """Fit mixed-effects hinge model for one protein"""
    df = create_hinge_terms(df.copy(), peak_day)
    
    # formula: protein ~ t_minus + t_plus + (1|Sample)
    endog = df[protein_col].values
    exog = df[['t_minus', 't_plus']]
    exog = sm.add_constant(exog)
    
    # random intercept per sample
    model = MixedLM(endog, exog, groups=df['Sample'])
    try:
        result = model.fit(reml=False)
        aic = result.aic
        return result, aic
    except:
        return None, np.inf

def select_best_peak(df, protein_col, candidate_peaks):
    """Try multiple peaks, return best-fit model and peak"""
    best_model = None
    best_peak = None
    best_aic = np.inf
    for peak in candidate_peaks:
        model, aic = fit_hinge_mixed_model(df, protein_col, peak)
        if aic < best_aic:
            best_model = model
            best_peak = peak
            best_aic = aic
    return best_model, best_peak, best_aic

def calculate_r2(result, df, protein_col):
    """Compute marginal R^2 (approx)"""
    # variance of fitted values over total variance
    fitted = result.fittedvalues
    var_fitted = np.var(fitted)
    var_total = np.var(df[protein_col].values)
    return var_fitted / var_total

# ---------------------------
# MAIN SCRIPT
# ---------------------------

# load data
df = pd.read_excel(INPUT_PATH, sheet_name=SHEET_NAME)

# identify protein columns
exclude_cols = ['Sample', 'Time', 'Group']
protein_cols = [c for c in df.columns if c not in exclude_cols]

# separate groups
groups = df['Group'].unique()
results_all = {}

import statsmodels.api as sm
import warnings
warnings.filterwarnings("ignore")

for group in groups:
    df_group = df[df['Group'] == group]
    results_list = []
    
    for protein in protein_cols:
        # fit hinge model for each protein
        model, peak, aic = select_best_peak(df_group, protein, PEAK_CANDIDATES)
        if model is not None:
            # extract slopes and intercept
            coefs = model.params
            # get p-values for slopes
            pvals = model.pvalues
            r2 = calculate_r2(model, df_group, protein)
            results_list.append({
                'Protein': protein,
                'Peak': peak,
                'AIC': aic,
                'Beta0': coefs.get('const', np.nan),
                'Beta1_pre': coefs.get('t_minus', np.nan),
                'Beta2_post': coefs.get('t_plus', np.nan),
                'P_Beta1': pvals.get('t_minus', np.nan),
                'P_Beta2': pvals.get('t_plus', np.nan),
                'R2': r2
            })
        else:
            results_list.append({
                'Protein': protein,
                'Peak': np.nan,
                'AIC': np.nan,
                'Beta0': np.nan,
                'Beta1_pre': np.nan,
                'Beta2_post': np.nan,
                'P_Beta1': np.nan,
                'P_Beta2': np.nan,
                'R2': np.nan
            })
    
    results_df = pd.DataFrame(results_list)
    results_all[group] = results_df

    # ---------------------------
    # VISUALISATIONS
    # ---------------------------
    # sort by R2
    best_proteins = results_df.sort_values('R2', ascending=False).head(TOP_N)['Protein']
    worst_proteins = results_df.sort_values('R2', ascending=True).head(TOP_N)['Protein']
    
    # function to plot proteins
    def plot_protein(df_group, protein, peak, save_path):
        df_plot = df_group[['Time', 'Sample', protein]].copy()
        df_plot = create_hinge_terms(df_plot, peak)
        plt.figure(figsize=(6,4))
        for s in df_plot['Sample'].unique():
            df_s = df_plot[df_plot['Sample'] == s]
            plt.plot(df_s['Time'], df_s[protein], 'o-', alpha=0.5, label=s)
        # overlay fitted hinge
        x_vals = np.linspace(df_plot['Time'].min(), df_plot['Time'].max(), 100)
        y_vals = []
        b0 = results_df.loc[results_df['Protein']==protein, 'Beta0'].values[0]
        b1 = results_df.loc[results_df['Protein']==protein, 'Beta1_pre'].values[0]
        b2 = results_df.loc[results_df['Protein']==protein, 'Beta2_post'].values[0]
        for x in x_vals:
            t_minus = min(0, x-peak)
            t_plus = max(0, x-peak)
            y_vals.append(b0 + b1*t_minus + b2*t_plus)
        plt.plot(x_vals, y_vals, 'r--', lw=2, label='Hinge fit')
        plt.title(f"{protein} ({group})")
        plt.xlabel("Time")
        plt.ylabel("Abundance")
        plt.legend(fontsize=6)
        plt.tight_layout()
        plt.savefig(os.path.join(save_path, f"{protein}_{group}.png"))
        plt.close()
    
    # make folder for plots
    plot_dir = os.path.join(OUTPUT_PATH, f"Plots_{group}")
    os.makedirs(plot_dir, exist_ok=True)
    
    for protein in list(best_proteins) + list(worst_proteins):
        peak = results_df.loc[results_df['Protein']==protein, 'Peak'].values[0]
        if not np.isnan(peak):
            plot_protein(df_group, protein, peak, plot_dir)

# ---------------------------
# SAVE RESULTS TO EXCEL
# ---------------------------
excel_path = os.path.join(OUTPUT_PATH, "Hinge_Model_Results.xlsx")
with pd.ExcelWriter(excel_path, engine='xlsxwriter') as writer:
    for group, df_res in results_all.items():
        df_res.to_excel(writer, sheet_name=group, index=False)

print("Mixed-effects hinge modeling complete. Results saved to Excel and plots folder.")
