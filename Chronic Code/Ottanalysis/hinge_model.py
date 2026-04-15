from logging import config

import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import ConvergenceWarning
import os
import warnings
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from word2number import w2n
import re
from helpers import convert_timepoint_to_number, derive_subject_ids, impute_missing_vals
import numpy as np
import traceback


'''
Creates a mixed-effects hinge model for each protein (or fixed-effects as a fallback), trying multiple candidate peak days and selecting the best fit based on AIC:

- finds an average peak (or trough) amongst participants 
- creates an initial upward or downward slope to that that best fits the data from all participants 
- creates an opposite slope on the other side of the peak/troughin the same way 
- looks at how far the participant's actual datapoints are from the created slopes looks at the gradient of the slopes
- generates a p-value from how significantly the slopes differ from 0 


Output should include (for each protein):
- average peak/trough
- gradient of upward slope
- gradient of downward slope
- how well the model fits the data
- p-values for anything it makes sense to have them for
- visualisation of the best fitting models (separate file)
'''


# Suppress convergence warnings - handling failed models with try/except
warnings.simplefilter('ignore', ConvergenceWarning)




def plot_best_fit(protein_name, df_long, best_peak, fit_results, output_path, subject_id_col, time_col, time_numeric):
    """
    Generates and saves a plot visualizing the raw data and the fitted hinge model.
    """
    # Extract coefficients
    intercept = fit_results.params['Intercept']
    slope1 = fit_results.params['slope1']
    slope2_change = fit_results.params['slope2_change']
    slope2 = slope1 + slope2_change
    
# Generate predicted values for the fixed effect line
    numeric_timepoints = sorted(df_long[time_numeric].unique())
    pred_vals = []
    
    for t_num in numeric_timepoints:
        slope1_val = t_num - best_peak
        slope2_change_val = max(0, t_num - best_peak)
        
        pred = intercept + (slope1 * slope1_val) + (slope2_change * slope2_change_val)
        pred_vals.append(pred)

    # Plotting
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot individual participant data
    sns.lineplot(data=df_long, x=time_numeric, y='value', units=subject_id_col, estimator=None,
                 color='grey', alpha=0.3, ax=ax, label='_nolegend_')

    # Plot the fixed-effect (average) hinge model
    ax.plot(numeric_timepoints, pred_vals, 'r-', linewidth=3, label='Fitted Hinge Model (Average)')
    
    # Add a vertical line at the peak
    ax.axvline(x=best_peak, color='b', linestyle='--', label=f'Peak/Trough at {best_peak:.2f}')
    
    # Map numeric timepoints back to original timpoint labels for x-axis labels
    time_label_map = dict(zip(df_long[time_numeric], df_long[time_col]))
    tick_locations = sorted(time_label_map.keys())
    tick_labels = [time_label_map[t] for t in tick_locations]
    ax.set_xticks(tick_locations)
    ax.set_xticklabels(tick_labels)

    ax.set_title(f'Mixed-Effects Hinge Model for: {protein_name}\n(Best Fit AIC: {fit_results.aic:.2f})')
    ax.set_xlabel('Time')
    ax.set_ylabel('Log2 Abundance')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close(fig)


def run_hinge_model(df, protein_cols, time_col, id_col, candidate_peaks, num_plots, output_dir, subject_id_delimiter):
    """
    Fits a mixed-effects hinge model for each protein

    For each protein, iterates through candidate peaks, selects the
    best model based on AIC, and saves the results and visualisations.
    """
    print("Creating hinge models (attempts mixed effects, falls back to fixed effects on failure)")

    df_copy = df.copy()
    
    # Extract subject IDs from id_col
    try:
        df_copy['subject_id_col'] = derive_subject_ids(df_copy[id_col], subject_id_delimiter)
    except Exception as e:
        print(f"Error extracting subject IDs: {e}")
        return

    # Create numeric time column
    df_copy['time_numeric'] = df_copy[time_col].astype(str).apply(convert_timepoint_to_number)

    # Drop rows where conversion failed
    failed_conversions = df_copy['time_numeric'].isna().sum()
    if failed_conversions > 0:
        print(f"Dropped rows where time couldn't be extracted - try reformating to timepoints to integers: {failed_conversions}")
        df_copy.dropna(subset=['time_numeric'], inplace=True)

    # Impute missing values in protein columns (imputes within subject)
    df_copy[protein_cols] = (
        df_copy
        .groupby('subject_id_col')[protein_cols]
        .apply(impute_missing_vals)
        .reset_index(level=0, drop=True)
    )

    # Convert candidate peaks
    candidate_peaks_numeric = [convert_timepoint_to_number(str(peak)) for peak in candidate_peaks]
    candidate_peaks_numeric = [p for p in candidate_peaks_numeric if p is not None] # Remove candidate peaks that couldn't be covnerted to a numeric value
    if not candidate_peaks_numeric:
        print("Candidate peaks couldn't be converted to numeric timepoints - try reformatting timepoints to integers")
        return

    all_results = []
    failed_protein_count = 0        # number of proteins where code fails to fit mixed effects model
    no_variance_count = 0

    # Use tqdm for a progress bar
    for protein in tqdm(protein_cols, desc="Processing proteins"):
        df_long = df_copy[['subject_id_col', time_col, 'time_numeric', protein]].copy()            # reformat dataframe to long format
        df_long.rename(columns={protein: 'value'}, inplace=True)                    # rename the protein column to a generic name for modelling
        df_long.dropna(subset=['value', 'time_numeric'], inplace=True)                    # drop columns where time or protein emasurment is missing
        
        if df_long['time_numeric'].nunique() < 3: # Need at least 3 unique time points
            print(f"Skipping {protein} - not enough unique time points to fit a hinge model")
            continue

        if df_long['value'].std() == 0:
            print(f"Skipping {protein} - no variance")
            no_variance_count += 1
            continue

        best_model = {'aic': float('inf')}

        # Iterate through numeric candidate peaks to find the best fit
        for peak in candidate_peaks_numeric:
            df_long['slope1'] = df_long['time_numeric'] - peak                            # Initial slope (before the peak)    - previously defined as time_numeric, now defined as relative to the peak to make th epredictors less correlated and improve model convergence
            df_long['slope2_change'] = (df_long['time_numeric'] - peak).clip(lower=0)
            
            mdf = None
            model_type = "None"


            try:
                ########## First attempt to fit mixed-effects model ############
                md = smf.mixedlm("value ~ slope1 + slope2_change", df_long, groups=df_long['subject_id_col'])
                mdf = md.fit(method=["lbfgs"])
                
                # FIX: statsmodels returns NaN for AIC if the covariance matrix is singular
                if pd.isna(mdf.aic):
                    raise ValueError("Mixed model converged but produced NaN AIC (singular covariance).")
                    
                model_type = "Mixed-Effects" # Added hyphen to match your later pandas filter
                
            except (ValueError, ZeroDivisionError, Exception) as e:
                # Silenced the print statement here so it doesn't spam your console 2000 times
                # print(f"Could not fit mixed-effects model for {protein} at peak {peak}: {e}")
                
                try:
                    ########## Fallback - attempt to fit fixed-effects (OLS) model ############
                    md_ols = smf.ols("value ~ slope1 + slope2_change", data=df_long)
                    mdf = md_ols.fit()
                    
                    # Also ensure OLS didn't return NaN
                    if pd.isna(mdf.aic):
                        continue
                        
                    model_type = "Fixed-Effects"
                except Exception as e:
                    # print(f"Could not fit fixed-effects model for {protein} at peak {peak}: {e}")
                    continue    # skip to next candidate peak if both models fail


            # Check if current model is better than the previous best
            if mdf.aic < best_model['aic']:
                best_model = {
                    'aic': mdf.aic,
                    'results': mdf,
                    'peak': peak,
                    'model_type': model_type
                }


        
        # Store the results of best fitting model for the protein
        if 'results' in best_model:
            res = best_model['results']
            slope1_coef = res.params.get('slope1', 0)
            slope2_change_coef = res.params.get('slope2_change', 0)

            numeric_to_label_map = dict(zip(df_copy['time_numeric'], df_copy[time_col]))
            best_peak_label = numeric_to_label_map.get(best_model['peak'], str(best_model['peak']))
            
            all_results.append({
                'Protein': protein,
                'Model_Type': best_model['model_type'],
                'Best_Peak': best_peak_label,
                'Best_Peak_Numeric': best_model['peak'],
                'AIC': best_model['aic'],
                'Slope_1_Gradient': slope1_coef,
                'Slope_2_Gradient': slope1_coef + slope2_change_coef,
                'Slope_2_Change': slope2_change_coef,
                'Slope_1_PValue': res.pvalues.get('slope1', 1),
                'Slope_2_Change_PValue': res.pvalues.get('slope2_change', 1),
                'Model_Converged': getattr(res, 'converged', 'N/A')             # fixed-effects/OLS models don't have convergence
            })
        else:
            failed_protein_count += 1
            print(f"Failed to fit any model for {protein} - skipping. Total failed so far: {failed_protein_count}")


    # Results print statements
    if not all_results:
        print("Hinge model analysis finished, but no valid models could be fitted")
        ############ DIAGNOSTICS #####################
        print(f"{no_variance_count} proteins were skipped due to no variance in measurements")
        print(f"Number of unique subjects: {df_copy['subject_id_col'].nunique()}")

        print("NaNs:", df_long.isna().sum())
        print("Unique times:", df_long['time_numeric'].unique())
        print("Subjects:", df_long['subject_id_col'].nunique())
        print("Value std:", df_long['value'].std())
        ############################################### 
        return
    if failed_protein_count > 0:
        print(f"\nNote: {failed_protein_count} proteins were skipped as no model (mixed or fixed) could be fitted.")
        ############ DIAGNOSTICS #####################
        print(f"{no_variance_count} proteins were skipped due to no variance in measurements")
        print(f"Number of unique subjects: {df_copy['subject_id_col'].nunique()}")

        print("NaNs:", df_long.isna().sum())
        print("Unique times:", df_long['time_numeric'].unique())
        print("Subjects:", df_long['subject_id_col'].nunique())
        print("Value std:", df_long['value'].std())
        ############################################### 
    else:
        print("\nHinge model analysis finished with successful fits for all proteins")



    # Create and save the final results dataframe
    df_results = pd.DataFrame(all_results)
    df_results.sort_values(by='Slope_2_Change_PValue', ascending=True, inplace=True)    # Sort by the most significant hinge effect
    output_filepath = os.path.join(output_dir, "hinge_model_results.xlsx")

    df_mixed_only = df_results[df_results['Model_Type'] == 'Mixed-Effects'][['Protein', 'Slope_2_Change_PValue']]  # Get list of proteins with fitted mixed-effects models
    df_fixed_only = df_results[df_results['Model_Type'] == 'Fixed-Effects'][['Protein', 'Slope_2_Change_PValue']]  # Get list of proteins with fitted fixed-effects models

    with pd.ExcelWriter(output_filepath, engine='xlsxwriter') as writer:
        df_results.to_excel(writer, sheet_name='All_Model_Results', index=False)
        df_fixed_only.to_excel(writer, sheet_name='Fixed_Effects_Models', index=False)
        df_mixed_only.to_excel(writer, sheet_name='Mixed_Effects_Models', index=False)
    print(f"\nHinge model results saved to: {output_filepath}")


    #################################################
    # Generate and save plots for the top N proteins
    plots_dir = os.path.join(output_dir, "Hinge_Model_Plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    top_proteins_to_plot = df_results.head(num_plots)
    
    print(f"Creating plots for the top {len(top_proteins_to_plot)} models (based on gradient change p-value)")
    for _, row in tqdm(top_proteins_to_plot.iterrows(), total=len(top_proteins_to_plot), desc="Generating plots"):
        protein_name = row['Protein']
        best_peak = row['Best_Peak_Numeric']
        
        # Re-create the long dataframe for plotting
        df_plot_long = df_copy[['subject_id_col', time_col, 'time_numeric', protein_name]].copy()
        df_plot_long.rename(columns={protein_name: 'value'}, inplace=True)
        df_plot_long.dropna(inplace=True)
        
        df_plot_long['slope1'] = df_plot_long['time_numeric'] - best_peak
        df_plot_long['slope2_change'] = (df_plot_long['time_numeric'] - best_peak).clip(lower=0)
        
        # Re-fit the best models for plotting (better than storing whole models)
        try:
            model_type = row['Model_Type']
            if model_type == "Mixed-Effects":
                md = smf.mixedlm("value ~ slope1 + slope2_change", df_plot_long, groups=df_plot_long['subject_id_col'])
                mdf = md.fit(method=["lbfgs"])
            else:
                md_ols = smf.ols("value ~ slope1 + slope2_change", data=df_plot_long)
                mdf = md_ols.fit()

            plot_filename = f"{protein_name.replace('/', '_')}_hinge_plot.png"
            plot_filepath = os.path.join(plots_dir, plot_filename)
            
            plot_best_fit(protein_name, df_plot_long, best_peak, mdf, plot_filepath,
                          subject_id_col='subject_id_col',
                          time_col=time_col,
                          time_numeric='time_numeric')
        except Exception as e:
            print(f"Could not generate plot for {protein_name}: {e}")

    print(f"Plots saved in: {plots_dir}")