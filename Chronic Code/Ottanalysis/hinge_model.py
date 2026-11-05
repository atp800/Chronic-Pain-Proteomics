from logging import config

import pandas as pd
import statsmodels.formula.api as smf
import os
from tempfile import mkdtemp
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
Creates a mixed-effects hinge model for each protein, selecting its candidate peak by AIC.
If all mixed fits fail, retains an OLS fallback for description only:

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


def _fit_mixed_model(df_long):
    """Retry the same random-intercept model; return the first usable fit and its diagnostics."""
    attempts = []
    for optimizer in ("lbfgs", "bfgs", "cg"):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                # Start each attempt afresh, including after singular-matrix errors.
                model = smf.mixedlm(
                    "value ~ slope1 + slope2_change", df_long,
                    groups=df_long['subject_id_col']
                )
                result = model.fit(method=optimizer, reml=False, maxiter=200)
                if not result.converged:
                    raise ValueError("did not converge")
                if not np.isfinite(result.aic) or not np.isfinite(result.params).all():
                    raise ValueError("non-finite AIC or coefficients")
                if not np.isfinite(result.bse_fe).all() or (result.bse_fe <= 0).any():
                    raise ValueError("unusable fixed-effect standard errors")
                if not np.isfinite(result.pvalues[['slope1', 'slope2_change']]).all():
                    raise ValueError("non-finite slope p-values")
            except Exception as exc:
                result = None
                status = f"{type(exc).__name__}: {exc}"
            else:
                status = "accepted"
        warning_text = "; ".join(dict.fromkeys(str(item.message) for item in caught))
        attempts.append(f"{optimizer}: {status}" + (f" [{warning_text}]" if warning_text else ""))
        if result is not None:
            return result, optimizer, "; ".join(attempts)
    return None, None, "; ".join(attempts)




def plot_best_fit(protein_name, df_long, model_result, output_path, subject_id_col, time_col, time_numeric, observed_means):
    """
    Generates and saves a plot visualizing the raw data and the fitted hinge model.
    """
    # Use the exported fit exactly; a plotting-time refit could fail or change it.
    best_peak = model_result['Best_Peak_Numeric']
    intercept = model_result['Intercept']
    slope1 = model_result['Slope_1_Gradient']
    slope2_change = model_result['Slope_2_Change']
    
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

    # Plot observed sample means calculated before imputation.
    ax.plot(
        observed_means.index,
        observed_means.to_numpy(),
        color='darkgreen',
        linestyle='--',
        marker='o',
        linewidth=2,
        markersize=6,
        label='Observed Mean (Before Imputation)',
        zorder=4
    )
    
    # Add a vertical line at the peak
    ax.axvline(x=best_peak, color='b', linestyle='--', label=f'Hinge point at {best_peak:.2f}')
    
    # Map numeric timepoints back to original timpoint labels for x-axis labels
    time_label_map = dict(zip(df_long[time_numeric], df_long[time_col]))
    tick_locations = sorted(time_label_map.keys())
    tick_labels = [time_label_map[t] for t in tick_locations]
    ax.set_xticks(tick_locations)
    ax.set_xticklabels(tick_labels)

    ax.set_title(
        f'Mixed-Effects Hinge Model for: {protein_name}\n'
        f'(Best Fit AIC: {model_result["AIC"]:.2f})'
    )
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
    best mixed model based on AIC, and saves the results and visualisations.
    If no mixed fit succeeds, saves the best descriptive OLS fit without p-values.
    """
    print("Creating hinge models (mixed-effects optimisers: lbfgs, bfgs, cg; OLS fallback is descriptive only)")

    if not protein_cols:
        print("Hinge analysis skipped: no protein columns found.")
        return

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

    # Calculate observed means at each time point before imputation.
    # Missing measurements are excluded separately for each protein.
    observed_means = (
        df_copy
        .groupby('time_numeric')[protein_cols]
        .mean()
        .sort_index()
    )

    # Impute missing values in protein columns (imputes within subject)
    df_copy[protein_cols] = (
        df_copy
        .groupby('subject_id_col')[protein_cols]
        .apply(impute_missing_vals)
        .reset_index(level=0, drop=True)
    )

    # Convert candidate peaks
    candidate_peaks_numeric = [convert_timepoint_to_number(str(peak)) for peak in candidate_peaks]
    candidate_peaks_numeric = list(dict.fromkeys(p for p in candidate_peaks_numeric if p is not None)) # Remove candidate peaks that couldn't be covnerted to a numeric value, avoiding duplicates
    if not candidate_peaks_numeric:
        print("Candidate peaks couldn't be converted to numeric timepoints - try reformatting timepoints to integers")
        return

    all_results = []
    all_candidate_results = []

    # Map numeric times back to their original labels, such as "D7".
    time_label_map = dict(
        zip(df_copy['time_numeric'], df_copy[time_col])
    )
    
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

        best_by_type = {}

        # Iterate through numeric candidate peaks to find the best fit
        for peak in candidate_peaks_numeric:
            # Create the row before fitting so failed candidates remain visible - a failed fit leaves a visible row with a blank p-value
            candidate_record = {
                'Protein': protein,
                'Candidate_Time': time_label_map.get(peak, str(peak)),
                'Candidate_Time_Numeric': peak,
                'Model_Type': None,
                'Inference_Status': 'Unavailable',
                'Mixed_Optimizer': None,
                'Mixed_Fit_Diagnostics': '',
                'AIC': np.nan,
                'Slope_2_Change_PValue': np.nan,
                'Status': 'Fit failed'
            }
            all_candidate_results.append(candidate_record)

            # A hinge requires measurements on both sides.
            earliest_time = df_long['time_numeric'].min()
            latest_time = df_long['time_numeric'].max()

            if not earliest_time < peak < latest_time:
                candidate_record['Status'] = 'No data on both sides'
                continue


            df_long['slope1'] = df_long['time_numeric'] - peak                            # Initial slope (before the peak)    - previously defined as time_numeric, now defined as relative to the peak to make th epredictors less correlated and improve model convergence
            df_long['slope2_change'] = (df_long['time_numeric'] - peak).clip(lower=0)
            
            mdf, optimizer, diagnostics = _fit_mixed_model(df_long)
            candidate_record['Mixed_Optimizer'] = optimizer
            candidate_record['Mixed_Fit_Diagnostics'] = diagnostics
            model_type = "Mixed-Effects"
            if mdf is None:
                model_type = "Fixed-Effects"
                try:
                    mdf = smf.ols("value ~ slope1 + slope2_change", data=df_long).fit()
                except Exception as e:
                    candidate_record['Status'] = f'OLS also failed: {e}'
                    continue

            if not np.isfinite(mdf.aic) or not np.isfinite(mdf.params).all():
                candidate_record['Status'] = 'Non-finite AIC or coefficients'
                continue

            is_mixed = model_type == 'Mixed-Effects'
            inference_status = ('Mixed-effects: nominal p-values' if is_mixed
                                else 'Descriptive only: repeated visits ignored')

            # Store this candidate's results, whether or not it wins.
            candidate_record.update({
                'Model_Type': model_type,
                'Inference_Status': inference_status,
                'AIC': mdf.aic,
                'Slope_2_Change_PValue': mdf.pvalues.get(
                    'slope2_change', np.nan
                ) if is_mixed else np.nan,
                'Status': 'Fitted'
            })

            # Compare AIC only within model type. A fallback cannot displace a mixed fit.
            if mdf.aic < best_by_type.get(model_type, {}).get('aic', float('inf')):
                best_by_type[model_type] = {
                    'aic': mdf.aic,
                    'results': mdf,
                    'peak': peak,
                    'model_type': model_type,
                    'inference_status': inference_status,
                    'optimizer': optimizer,
                    'diagnostics': diagnostics,
                }


        
        # Store the results of best fitting model for the protein
        best_model = best_by_type.get('Mixed-Effects', best_by_type.get('Fixed-Effects', {}))
        if 'results' in best_model:
            res = best_model['results']
            is_mixed = best_model['model_type'] == 'Mixed-Effects'
            slope1_coef = res.params.get('slope1', 0)
            slope2_change_coef = res.params.get('slope2_change', 0)

            numeric_to_label_map = dict(zip(df_copy['time_numeric'], df_copy[time_col]))
            best_peak_label = numeric_to_label_map.get(best_model['peak'], str(best_model['peak']))
            
            all_results.append({
                'Protein': protein,
                'Model_Type': best_model['model_type'],
                'Inference_Status': best_model['inference_status'],
                'Mixed_Optimizer': best_model['optimizer'],
                'Mixed_Fit_Diagnostics': best_model['diagnostics'],
                'Best_Peak': best_peak_label,
                'Best_Peak_Numeric': best_model['peak'],
                'AIC': best_model['aic'],
                'Intercept': res.params['Intercept'],
                'Slope_1_Gradient': slope1_coef,
                'Slope_2_Gradient': slope1_coef + slope2_change_coef,
                'Slope_2_Change': slope2_change_coef,
                'Slope_1_PValue': res.pvalues.get('slope1', np.nan) if is_mixed else np.nan,
                'Slope_2_Change_PValue': res.pvalues.get('slope2_change', np.nan) if is_mixed else np.nan,
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
        # If no valid models were fitted, still save the candidate diagnostics for inspection:
        if all_candidate_results:
            diagnostic_path = os.path.join(
                output_dir,
                "hinge_candidate_diagnostics.xlsx"
            )

            df_diagnostics = pd.DataFrame(all_candidate_results)
            df_diagnostics['Is_Best_Hinge'] = False

            df_diagnostics.to_excel(
                diagnostic_path,
                sheet_name='All_Candidate_Timepoints',
                index=False,
                engine='xlsxwriter'
            )

            print(f"All fittings failed. Candidate diagnostics saved to: {diagnostic_path}")
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
    df_results.sort_values(by='Slope_2_Change_PValue', ascending=True, inplace=True)    # Descriptive fallbacks have blank p-values and sort last.
    output_filepath = os.path.join(output_dir, "hinge_model_results.xlsx")

    df_mixed_only = df_results[df_results['Model_Type'] == 'Mixed-Effects']
    df_fixed_only = df_results[df_results['Model_Type'] == 'Fixed-Effects']
    top_proteins_to_plot = df_mixed_only.head(num_plots)
    print(f"Mixed-effects results: {len(df_mixed_only)}; descriptive OLS fallbacks: {len(df_fixed_only)}.")


    # Convert the candidate records into a table.
    df_candidates = pd.DataFrame(all_candidate_results)

    # Look up the selected hinge time for each protein.
    best_peak_by_protein = df_results.set_index(
        'Protein'
    )['Best_Peak_Numeric']

    df_candidates['Is_Best_Hinge'] = (
        df_candidates['Status'].eq('Fitted')
        & df_candidates['Candidate_Time_Numeric'].eq(
            df_candidates['Protein'].map(best_peak_by_protein)
        )
    )

    df_candidates.sort_values(
        by=['Protein', 'Candidate_Time_Numeric'],
        inplace=True
    )

    with pd.ExcelWriter(output_filepath, engine='xlsxwriter') as writer:
        df_results.to_excel(
            writer, sheet_name='All_Model_Results', index=False)
        df_fixed_only.to_excel(
            writer, sheet_name='Fixed_Effects_Models', index=False)
        df_mixed_only.to_excel(
            writer, sheet_name='Mixed_Effects_Models', index=False)
        df_candidates.to_excel(
            writer, sheet_name='All_Candidate_Timepoints', index=False)
        top_proteins_to_plot.to_excel(
            writer, sheet_name='Top_Mixed_Models', index=False)
    print(f"\nHinge model results saved to: {output_filepath}")


    #################################################
    # Generate and save plots for the top N proteins
    plots_dir = os.path.join(output_dir, "Hinge_Model_Plots")
    os.makedirs(plots_dir, exist_ok=True)
    # Preserve old generated plots without leaving earlier fallback fits in the current selection.
    with os.scandir(plots_dir) as entries:
        previous_plots = [entry for entry in entries
                          if entry.is_file(follow_symlinks=False) and entry.name.endswith('_hinge_plot.png')]
    if previous_plots:
        archive_dir = mkdtemp(prefix='previous_', dir=plots_dir)
        for entry in previous_plots:
            os.replace(entry.path, os.path.join(archive_dir, entry.name))
        print(f"Previous hinge plots preserved in: {archive_dir}")
    
    print(f"Creating plots for the top {len(top_proteins_to_plot)} mixed-effects models (based on nominal gradient change p-value)")
    for _, row in tqdm(top_proteins_to_plot.iterrows(), total=len(top_proteins_to_plot), desc="Generating plots"):
        protein_name = row['Protein']
        
        # Re-create the long dataframe for plotting
        df_plot_long = df_copy[['subject_id_col', time_col, 'time_numeric', protein_name]].copy()
        df_plot_long.rename(columns={protein_name: 'value'}, inplace=True)
        df_plot_long.dropna(inplace=True)
        
        try:
            plot_filename = f"{protein_name.replace('/', '_')}_hinge_plot.png"
            plot_filepath = os.path.join(plots_dir, plot_filename)
            
            plot_best_fit(protein_name, 
                          df_plot_long, 
                          row,
                          plot_filepath,
                          subject_id_col='subject_id_col',
                          time_col=time_col,
                          time_numeric='time_numeric',
                          observed_means=observed_means[protein_name])
        except Exception as e:
            print(f"Could not generate plot for {protein_name}: {e}")

    print(f"Plots saved in: {plots_dir}")
    return output_filepath
