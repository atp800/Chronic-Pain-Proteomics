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
from helpers import convert_timepoint_to_number, derive_subject_ids
from hinge_inference import (DEFAULT_BOOTSTRAPS, DEFAULT_SEED, adjust_hinge_fdr,
                             bootstrap_hinge_test, unavailable_test, validate_test_settings)
import numpy as np
import traceback


'''
Fits a random-intercept hinge model to each protein's observed measurements.
The minimum-AIC candidate estimates the bend day; the slopes need not reverse.
Significance compares the full candidate search with a linear no-bend mixed
model using a parametric bootstrap, then adjusts across proteins using BH FDR.
Ordinary coefficient p-values are nominal descriptions of the selected fit.
If all mixed fits fail, the best OLS fallback is retained for description only.
'''


def _fit_mixed_model(df_long, *, hinge=True):
    """Retry the same random-intercept model; return the first usable fit and its diagnostics."""
    attempts = []
    for optimizer in ("lbfgs", "bfgs", "cg"):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                # Start each attempt afresh, including after singular-matrix errors.
                model = smf.mixedlm(
                    "value ~ slope1 + slope2_change" if hinge else "value ~ slope1", df_long,
                    groups=df_long['subject_id_col']
                )
                result = model.fit(method=optimizer, reml=False, maxiter=200)
                if not result.converged:
                    raise ValueError("did not converge")
                if not np.isfinite(result.aic) or not np.isfinite(result.params).all():
                    raise ValueError("non-finite AIC or coefficients")
                if not np.isfinite(result.bse_fe).all() or (result.bse_fe <= 0).any():
                    raise ValueError("unusable fixed-effect standard errors")
                slope_terms = ['slope1', 'slope2_change'] if hinge else ['slope1']
                if not np.isfinite(result.pvalues[slope_terms]).all():
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




def plot_best_fit(protein_name, df_long, model_result, output_path, subject_id_col, time_col, time_numeric, observed_means, *, rank_by_nominal=False):
    """
    Generates and saves a plot visualizing the raw data and the fitted hinge model.
    """
    # Use the exported fit exactly; a plotting-time refit could fail or change it.
    best_peak = model_result['Best_Peak_Numeric']
    intercept = model_result['Intercept']
    slope1 = model_result['Slope_1_Gradient']
    slope2_change = model_result['Slope_2_Change']
    
# Generate predicted values for the fixed effect line
    numeric_timepoints = sorted(set(df_long[time_numeric]) | {best_peak})
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
                 color='grey', alpha=0.3, marker='o', markersize=3, ax=ax, label='_nolegend_')

    # Plot the fixed-effect (average) hinge model
    ax.plot(numeric_timepoints, pred_vals, 'r-', linewidth=3, label='Fitted Hinge Model (Average)')

    # The participants contributing to each observed mean can differ across visits.
    ax.plot(
        observed_means.index,
        observed_means.to_numpy(),
        color='darkgreen',
        linestyle='--',
        marker='o',
        linewidth=2,
        markersize=6,
        label='Observed Mean (Available Measurements)',
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

    p_label = (
        'search p upper bound'
        if model_result['Search_PValue_Is_Upper_Bound']
        else 'search p'
    )

    ranking_label = (
        'nominal slope-change p-value (exploratory)'
        if rank_by_nominal
        else 'search-adjusted BH FDR'
    )

    ax.set_title(
        f'Mixed-Effects Hinge Model for: {protein_name}\n'
        f'Best Fit AIC: {model_result["AIC"]:.2f}; '
        f'nominal slope-change p: '
        f'{model_result["Slope_2_Change_PValue_Nominal"]:.3g}\n'
        f'{p_label}: {model_result["Hinge_Search_PValue"]:.3g}; '
        f'BH FDR: {model_result["Hinge_FDR"]:.3g}\n'
        f'Ranked by: {ranking_label}'
    )
    ax.set_xlabel('Time')
    ax.set_ylabel('Log2 Abundance')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close(fig)


def run_hinge_model(df, protein_cols, time_col, id_col, candidate_peaks, num_plots, output_dir, subject_id_delimiter,
                    *, n_bootstrap=DEFAULT_BOOTSTRAPS, random_seed=DEFAULT_SEED, fdr_threshold=.05, rank_by_nominal=False):
    """
    Fits a mixed-effects hinge model for each protein

    For each protein, iterates through candidate peaks, selects the
    best mixed model based on AIC, and saves the results and visualisations.
    Significance uses a full-search bootstrap and BH correction across proteins.
    If no mixed fit succeeds, saves the best descriptive OLS fit without p-values.
    Missing outcomes are omitted separately for each protein; incomplete participants remain.
    """
    validate_test_settings(n_bootstrap, random_seed, fdr_threshold)
    print("Creating hinge models (mixed-effects optimisers: lbfgs, bfgs, cg; OLS fallback is descriptive only)")
    print("Using observed protein measurements only; missing abundances are not imputed.")
    print(f"Whole-search bootstrap: up to {n_bootstrap} simulations per protein; minimum p = {1 / (n_bootstrap + 1):.6g}.")

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

    # Non-finite abundances are missing outcomes, not artificial low abundances.
    df_copy[protein_cols] = df_copy[protein_cols].replace([np.inf, -np.inf], np.nan)
    # Missing measurements are excluded separately for each protein.
    observed_means = (
        df_copy
        .groupby('time_numeric')[protein_cols]
        .mean()
        .sort_index()
    )

    # Convert candidate peaks
    candidate_peaks_numeric = [convert_timepoint_to_number(str(peak)) for peak in candidate_peaks]
    candidate_peaks_numeric = list(dict.fromkeys(p for p in candidate_peaks_numeric if p is not None)) # Remove candidate peaks that couldn't be covnerted to a numeric value, avoiding duplicates
    if not candidate_peaks_numeric:
        print("Candidate peaks couldn't be converted to numeric timepoints - try reformatting timepoints to integers")
        return

    all_results = []
    all_candidate_results = []
    all_coverage = []
    coverage_by_time = []
    analysis_info = pd.DataFrame([
        {'Setting': 'Missing outcomes', 'Value': 'Observed finite abundances only; no imputation. Incomplete participants retained.'},
        {'Setting': 'Mixed-model replication', 'Value': 'At least two participants with observations at two or more distinct times. This is a minimum safeguard, not a guarantee of precision.'},
        {'Setting': 'Candidate comparison', 'Value': 'All candidate hinges for a protein use the same observed rows; AIC compared within model type.'},
        {'Setting': 'Missingness assumption', 'Value': 'Likelihood inference assumes missing at random given the model. Abundance-dependent non-detection may still bias estimates.'},
        {'Setting': 'Search test', 'Value': 'Parametric bootstrap of maximum mixed-model likelihood gain over a linear no-bend mixed model; repeats all eligible candidate days.'},
        {'Setting': 'Bootstrap target', 'Value': n_bootstrap},
        {'Setting': 'Bootstrap seed', 'Value': random_seed},
        {'Setting': 'Minimum bootstrap p-value', 'Value': 1 / (n_bootstrap + 1)},
        {'Setting': 'Bootstrap failures/stopping', 'Value': 'Failed draws count as exceedances. Early stopping only after significance is impossible; uncomputed draws also count as exceedances, yielding a labelled upper bound.'},
        {'Setting': 'FDR method', 'Value': 'Benjamini-Hochberg (BH); assumes independent or suitably positively dependent protein tests.'},
        {'Setting': 'FDR family', 'Value': f'All {len(protein_cols)} requested proteins in this group/sheet. Skipped/unavailable tests count as p=1 for adjustment and remain blank in results.'},
        {'Setting': 'FDR threshold', 'Value': fdr_threshold},
        {'Setting': 'Reporting mode', 'Value': 'Nominal' if rank_by_nominal else 'FDR'},
        {'Setting': 'Ranking criterion', 'Value': 'AIC-selected hinge nominal slope-change p-value' if rank_by_nominal else 'Search-adjusted BH FDR, then whole-search p-value'},
        {'Setting': 'Significant_Hinges criterion', 'Value': f'Nominal slope-change p < {fdr_threshold}' if rank_by_nominal else f'Search-adjusted BH FDR <= {fdr_threshold}'},
        {'Setting': 'Candidate p-value columns', 'Value': 'One nominal slope-change p-value per configured candidate hinge; blanks indicate unavailable tests.'},
        {'Setting': 'Nominal inference warning', 'Value': 'Nominal coefficient p-values are unadjusted for hinge-day selection and testing multiple proteins. Best_Peak is an estimate, not a confirmed exact day.'},
    ])

    # Map numeric times back to their original labels, such as "D7".
    time_label_map = dict(
        zip(df_copy['time_numeric'], df_copy[time_col])
    )
    
    failed_protein_count = 0

    # Use tqdm for a progress bar
    for protein in tqdm(protein_cols, desc="Processing proteins"):
        df_long = df_copy[['subject_id_col', time_col, 'time_numeric', protein]].copy()            # reformat dataframe to long format
        df_long.rename(columns={protein: 'value'}, inplace=True)                    # rename the protein column to a generic name for modelling
        counts_by_time = df_long.groupby('time_numeric')['value'].agg(['count', 'size'])
        counts_by_time = counts_by_time.rename(columns={'count': 'N_Observed', 'size': 'N_Recorded_Visits'})
        counts_by_time['N_Missing'] = counts_by_time['N_Recorded_Visits'] - counts_by_time['N_Observed']
        counts_by_time['Time_Label'] = counts_by_time.index.map(time_label_map)
        counts_by_time.insert(0, 'Protein', protein)
        coverage_by_time.append(counts_by_time.reset_index())

        # Keep every available measurement, including a participant's lone observation.
        # The identical observed rows are used for every candidate, optimiser and fallback.
        df_long.dropna(subset=['value', 'time_numeric'], inplace=True)
        observations_per_person = df_long.groupby('subject_id_col').size()
        times_per_person = df_long.groupby('subject_id_col')['time_numeric'].nunique()
        coverage = {
            'N_Observed': len(df_long),
            'N_Missing': len(df_copy) - len(df_long),
            'N_Participants_Observed': len(observations_per_person),
            'N_Participants_Repeated': int(times_per_person.ge(2).sum()),
            'N_Participants_Single_Observation': int(observations_per_person.eq(1).sum()),
            'N_Participants_Unobserved': df_copy['subject_id_col'].nunique() - len(observations_per_person),
            'N_Observed_Timepoints': df_long['time_numeric'].nunique(),
        }
        coverage_record = {'Protein': protein, **coverage, 'Status': 'No usable candidate fit'}
        all_coverage.append(coverage_record)
        protein_issue = None
        if coverage['N_Observed_Timepoints'] < 3:
            protein_issue = 'Fewer than three observed time points'
        elif coverage['N_Observed'] <= 3:
            protein_issue = 'Insufficient observed residual degrees of freedom'
        elif df_long['value'].nunique() < 2:
            protein_issue = 'No observed abundance variation'

        best_by_type = {}
        candidate_llfs = {}

        # Iterate through numeric candidate peaks to find the best fit
        for peak in candidate_peaks_numeric:
            # Create the row before fitting so failed candidates remain visible - a failed fit leaves a visible row with a blank p-value
            candidate_record = {
                'Protein': protein,
                **coverage,
                'Candidate_Time': time_label_map.get(peak, str(peak)),
                'Candidate_Time_Numeric': peak,
                'Model_Type': None,
                'Inference_Status': 'Unavailable',
                'Mixed_Optimizer': None,
                'Mixed_Fit_Diagnostics': '',
                'AIC': np.nan,
                'Slope_2_Change_PValue_Nominal': np.nan,
                'Status': 'Fit failed'
            }
            all_candidate_results.append(candidate_record)

            if protein_issue:
                candidate_record['Status'] = protein_issue
                coverage_record['Status'] = protein_issue
                continue

            # A hinge requires measurements on both sides.
            earliest_time = df_long['time_numeric'].min()
            latest_time = df_long['time_numeric'].max()

            if not earliest_time < peak < latest_time:
                candidate_record['Status'] = 'No data on both sides'
                continue


            df_long['slope1'] = df_long['time_numeric'] - peak                            # Initial slope (before the peak)    - previously defined as time_numeric, now defined as relative to the peak to make th epredictors less correlated and improve model convergence
            df_long['slope2_change'] = (df_long['time_numeric'] - peak).clip(lower=0)

            design = np.column_stack([np.ones(len(df_long)), df_long['slope1'], df_long['slope2_change']])
            if np.linalg.matrix_rank(design) < 3:
                candidate_record['Status'] = 'Observed hinge design is not estimable'
                continue
            candidate_llfs[peak] = np.nan
            
            if coverage['N_Participants_Repeated'] >= 2:
                mdf, optimizer, diagnostics = _fit_mixed_model(df_long)
            else:
                mdf, optimizer = None, None
                diagnostics = 'Mixed fit not attempted: fewer than two participants with repeated observed measurements'
            candidate_record['Mixed_Optimizer'] = optimizer
            candidate_record['Mixed_Fit_Diagnostics'] = diagnostics
            model_type = "Mixed-Effects"
            if mdf is not None:
                candidate_llfs[peak] = mdf.llf
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
            inference_status = ('Mixed-effects: candidate coefficients are nominal only' if is_mixed
                                else 'Descriptive only: repeated visits ignored')

            # Store this candidate's results, whether or not it wins.
            candidate_record.update({
                'Model_Type': model_type,
                'Inference_Status': inference_status,
                'AIC': mdf.aic,
                'Slope_2_Change_PValue_Nominal': mdf.pvalues.get(
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
            test = unavailable_test(n_bootstrap, random_seed, 'Not tested: descriptive OLS fallback')
            if is_mixed:
                test = bootstrap_hinge_test(
                    df_long, candidate_llfs, _fit_mixed_model, protein=protein,
                    n_bootstrap=n_bootstrap, seed=random_seed, alpha=fdr_threshold,
                )
            coverage_record['Status'] = 'Mixed-effects fitted' if is_mixed else 'Descriptive OLS fallback'
            slope1_coef = res.params.get('slope1', 0)
            slope2_change_coef = res.params.get('slope2_change', 0)

            numeric_to_label_map = dict(zip(df_copy['time_numeric'], df_copy[time_col]))
            best_peak_label = numeric_to_label_map.get(best_model['peak'], str(best_model['peak']))
            
            all_results.append({
                'Protein': protein,
                **coverage,
                **test,
                'Model_Type': best_model['model_type'],
                'Inference_Status': ('Mixed-effects: search-adjusted inference' if np.isfinite(test['Hinge_Search_PValue'])
                                     else 'Mixed-effects: search inference unavailable') if is_mixed else best_model['inference_status'],
                'Mixed_Optimizer': best_model['optimizer'],
                'Mixed_Fit_Diagnostics': best_model['diagnostics'],
                'Best_Peak': best_peak_label,
                'Best_Peak_Numeric': best_model['peak'],
                'AIC': best_model['aic'],
                'Intercept': res.params['Intercept'],
                'Slope_1_Gradient': slope1_coef,
                'Slope_2_Gradient': slope1_coef + slope2_change_coef,
                'Slope_2_Change': slope2_change_coef,
                'Slope_1_PValue_Nominal': res.pvalues.get('slope1', np.nan) if is_mixed else np.nan,
                'Slope_2_Change_PValue_Nominal': res.pvalues.get('slope2_change', np.nan) if is_mixed else np.nan,
                'Model_Converged': getattr(res, 'converged', 'N/A')             # fixed-effects/OLS models don't have convergence
            })
        else:
            failed_protein_count += 1
            print(f"Failed to fit any model for {protein} - skipping. Total failed so far: {failed_protein_count}")


    df_coverage = pd.DataFrame(all_coverage)
    df_time_coverage = pd.concat(coverage_by_time, ignore_index=True)

    if not all_results:
        print("Hinge model analysis finished, but no usable models could be fitted to the observed data.")
        diagnostic_path = os.path.join(output_dir, "hinge_candidate_diagnostics.xlsx")
        df_diagnostics = pd.DataFrame(all_candidate_results)
        df_diagnostics['Is_Best_Hinge'] = False
        with pd.ExcelWriter(diagnostic_path, engine='xlsxwriter') as writer:
            df_diagnostics.to_excel(writer, sheet_name='All_Candidate_Timepoints', index=False)
            df_coverage.to_excel(writer, sheet_name='Data_Coverage', index=False)
            df_time_coverage.to_excel(writer, sheet_name='Observed_By_Timepoint', index=False)
            analysis_info.to_excel(writer, sheet_name='Analysis_Info', index=False)
        print(f"Coverage and candidate diagnostics saved to: {diagnostic_path}")
        return
    if failed_protein_count > 0:
        print(f"\nNote: {failed_protein_count} proteins had insufficient observed data or no usable fit; see Data_Coverage and All_Candidate_Timepoints.")
    else:
        print("\nHinge model analysis finished with successful fits for all proteins")



    # The statistical fits and FDR calculation are identical in both modes.
    df_results = adjust_hinge_fdr(
        pd.DataFrame(all_results), len(protein_cols), fdr_threshold
    )
    df_candidates = pd.DataFrame(all_candidate_results)

    # Add nominal slope-change p-values for every configured candidate hinge.
    pivot_pvals = df_candidates.pivot(
        index='Protein',
        columns='Candidate_Time_Numeric',
        values='Slope_2_Change_PValue_Nominal'
    )

    # Put candidate columns in chronological order.
    pivot_pvals = pivot_pvals.reindex(
        columns=sorted(candidate_peaks_numeric)
    )

    # Label the new columns using the original timepoint labels.
    pivot_pvals.columns = [
        f"Peak_{time_label_map.get(t, str(t))}_Nominal_PValue"
        for t in pivot_pvals.columns
    ]

    # Merge all candidate p-values into the main results.
    df_results = df_results.merge(
        pivot_pvals,
        left_on='Protein',
        right_index=True,
        how='left',
        validate='one_to_one',
        sort=False
    )

    # Keep both significance definitions in the output workbook.
    df_results['Significant_Nominal'] = (
        df_results['Model_Type'].eq('Mixed-Effects')
        & df_results['Slope_2_Change_PValue_Nominal'].lt(fdr_threshold)
    )

    # Select the ranking and significance criteria.
    if rank_by_nominal:
        sort_columns = [
            'Slope_2_Change_PValue_Nominal',
            'Protein'
        ]
        ranking_column = 'Slope_2_Change_PValue_Nominal'
        significance_column = 'Significant_Nominal'
        ranking_description = 'nominal slope-change p-value (exploratory)'
    else:
        sort_columns = [
            'Hinge_FDR',
            'Hinge_Search_PValue'
        ]
        ranking_column = 'Hinge_FDR'
        significance_column = 'Significant_FDR'
        ranking_description = 'search-adjusted BH FDR'

    # Rank proteins according to the selected mode.
    df_results.sort_values(
        by=sort_columns,
        kind='stable',
        na_position='last',
        inplace=True
    )

    output_filepath = os.path.join(
        output_dir, 'hinge_model_results.xlsx'
    )

    # Keep mixed-effects and descriptive fallback models separate.
    df_mixed_only = df_results[
        df_results['Model_Type'] == 'Mixed-Effects'
    ]
    df_fixed_only = df_results[
        df_results['Model_Type'] == 'Fixed-Effects'
    ]

    # Select top N proteins using the configured ranking criterion.
    top_proteins_to_plot = df_mixed_only[
        df_mixed_only[ranking_column].notna()
    ].head(num_plots)

    # This determines which proteins go into Significant_Hinges.
    significant_mask = df_results[significance_column]

    print(
        f"Mixed-effects results: {len(df_mixed_only)}; "
        f"descriptive OLS fallbacks: {len(df_fixed_only)}."
    )
    print(
        f"FDR-significant proteins: "
        f"{int(df_results['Significant_FDR'].sum())}; "
        f"nominally significant proteins: "
        f"{int(df_results['Significant_Nominal'].sum())}."
    )
    print(f"Ranking mode: {ranking_description}.")


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
        df_results[significant_mask].to_excel(
            writer, sheet_name='Significant_Hinges', index=False)
        df_coverage.to_excel(writer, sheet_name='Data_Coverage', index=False)
        df_time_coverage.to_excel(writer, sheet_name='Observed_By_Timepoint', index=False)
        analysis_info.to_excel(writer, sheet_name='Analysis_Info', index=False)
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
    
    print(
        f"Creating plots for the top {len(top_proteins_to_plot)} "
        f"mixed-effects models ranked by {ranking_description} "
        f"(top-ranked does not imply significant)."
    )
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
                          observed_means=observed_means[protein_name],
                          rank_by_nominal=rank_by_nominal)
        except Exception as e:
            print(f"Could not generate plot for {protein_name}: {e}")

    print(f"Plots saved in: {plots_dir}")
    return output_filepath
