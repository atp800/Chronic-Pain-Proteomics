import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import ConvergenceWarning
import os
import warnings
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm



'''
Creates a mixed-effects hinge model for each protein, trying multiple candidate peak days and selecting the best fit based on AIC:

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

Gui options to add:
- hinge model model options (candidate peak days, number of proteins to include in visualisation)
- hinge model selection (checkbox to tick to run)
'''


# Suppress convergence warnings - handling failed models with try/except
warnings.simplefilter('ignore', ConvergenceWarning)


def plot_best_fit(protein_name, df_long, best_peak, fit_results, output_path):
    """
    Generates and saves a plot visualizing the raw data and the fitted hinge model.
    """
    time_col = df_long.columns[1] # Assumes time is the second column
    id_col = df_long.columns[0]   # Assumes id is the first column
    
    # Extract coefficients
    intercept = fit_results.params['Intercept']
    slope1 = fit_results.params['slope1']
    slope2_change = fit_results.params['slope2_change']
    slope2 = slope1 + slope2_change
    
    # Generate predicted values for the fixed effect line
    time_points = sorted(df_long[time_col].unique())
    pred_vals = []
    for t in time_points:
        if t <= best_peak:
            pred_vals.append(intercept + slope1 * t)
        else:
            pred_vals.append(intercept + slope1 * t + slope2_change * (t - best_peak))

    # Plotting
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot individual participant data
    sns.lineplot(data=df_long, x=time_col, y='value', units=id_col, estimator=None,
                 color='grey', alpha=0.3, ax=ax, label='_nolegend_')

    # Plot the fixed-effect (average) hinge model
    ax.plot(time_points, pred_vals, 'r-', linewidth=3, label='Fitted Hinge Model (Average)')
    
    # Add a vertical line at the peak
    ax.axvline(x=best_peak, color='b', linestyle='--', label=f'Peak/Trough at {best_peak:.2f}')
    
    ax.set_title(f'Mixed-Effects Hinge Model for: {protein_name}\n(Best Fit AIC: {fit_results.aic:.2f})')
    ax.set_xlabel('Time')
    ax.set_ylabel('Log2 Abundance')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close(fig)


def run_hinge_model(df, protein_cols, time_col, patient_id_col, candidate_peaks, num_plots, output_dir):
    """
    Fits a mixed-effects hinge model for each protein.

    For each protein, it iterates through candidate peaks (knots), selects the
    best model based on AIC, and saves the results and visualizations.
    """
    print("Starting Hinge Model analysis...")
    
    # Ensure time and peak values are numeric
    try:
        df[time_col] = pd.to_numeric(df[time_col], errors='coerce')
        candidate_peaks = [float(p) for p in candidate_peaks]
    except (ValueError, TypeError) as e:
        print(f"Error: Could not convert Time column or Candidate Peaks to numbers. Please check your data. Details: {e}")
        return

    all_results = []
    
    # Use tqdm for a progress bar
    for protein in tqdm(protein_cols, desc="Processing proteins"):
        
        # 1. Prepare data in "long" format for the current protein
        df_long = df[[patient_id_col, time_col, protein]].copy()
        df_long.rename(columns={protein: 'value'}, inplace=True)
        df_long.dropna(subset=['value', time_col], inplace=True)
        
        if len(df_long) < 3: # Not enough data to fit a model
            continue

        best_model = {'aic': float('inf')}

        # 2. Iterate through candidate peaks to find the best fit
        for peak in candidate_peaks:
            
            # Create the predictor variables for the hinge function
            df_long['slope1'] = df_long[time_col]
            # The change in slope occurs *after* the peak
            df_long['slope2_change'] = (df_long[time_col] - peak).clip(lower=0)
            
            try:
                # 3. Fit the mixed-effects model
                # Formula: value depends on slope1 and the change in slope after the peak
                # Random effect: each patient has their own intercept
                md = smf.mixedlm(
                    "value ~ slope1 + slope2_change", 
                    df_long, 
                    groups=df_long[patient_id_col]
                )
                mdf = md.fit(method=["lbfgs"]) # lbfgs is a robust optimizer

                # 4. Check if this model is better than the previous best
                if mdf.aic < best_model['aic']:
                    best_model = {
                        'aic': mdf.aic,
                        'results': mdf,
                        'peak': peak,
                        'data': df_long
                    }

            except (ValueError, ZeroDivisionError, Exception) as e:
                # This can happen if data for a given peak is singular or fails to converge
                # print(f"Could not fit model for {protein} at peak {peak}: {e}")
                pass
        
        # 5. Store the results of the best fitting model for the protein
        if 'results' in best_model:
            res = best_model['results']
            slope1_coef = res.params.get('slope1', 0)
            slope2_change_coef = res.params.get('slope2_change', 0)
            
            all_results.append({
                'Protein': protein,
                'Best_Peak': best_model['peak'],
                'AIC': best_model['aic'],
                'Slope_1_Gradient': slope1_coef,
                'Slope_2_Gradient': slope1_coef + slope2_change_coef,
                'Slope_1_PValue': res.pvalues.get('slope1', 1),
                'Slope_2_Change_PValue': res.pvalues.get('slope2_change', 1),
                'Model_Converged': res.converged
            })

    if not all_results:
        print("Hinge model analysis finished, but no valid models could be fitted.")
        return

    # 6. Create and save the final results dataframe
    df_results = pd.DataFrame(all_results)
    # Sort by the most significant "hinge" effect
    df_results.sort_values(by='Slope_2_Change_PValue', ascending=True, inplace=True)
    
    output_filepath = os.path.join(output_dir, "hinge_model_results.csv")
    df_results.to_csv(output_filepath, index=False)
    print(f"\nHinge model results saved to: {output_filepath}")

    # 7. Generate and save plots for the top N proteins
    plots_dir = os.path.join(output_dir, "Hinge_Model_Plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    top_proteins_to_plot = df_results.head(num_plots)
    
    print(f"Generating {len(top_proteins_to_plot)} plots...")
    for _, row in tqdm(top_proteins_to_plot.iterrows(), total=len(top_proteins_to_plot), desc="Generating plots"):
        protein_name = row['Protein']
        
        # Re-create the long dataframe for plotting
        df_plot_long = df[[patient_id_col, time_col, protein_name]].copy()
        df_plot_long.rename(columns={protein_name: 'value'}, inplace=True)
        df_plot_long.dropna(inplace=True)
        
        # We need to re-fit the best model to get the results object back for plotting
        # This is simpler than storing large model objects in memory
        best_peak = row['Best_Peak']
        df_plot_long['slope1'] = df_plot_long[time_col]
        df_plot_long['slope2_change'] = (df_plot_long[time_col] - best_peak).clip(lower=0)
        
        md = smf.mixedlm("value ~ slope1 + slope2_change", df_plot_long, groups=df_plot_long[patient_id_col])
        mdf = md.fit(method=["lbfgs"])
        
        plot_filename = f"{protein_name.replace('/', '_')}_hinge_plot.png"
        plot_filepath = os.path.join(plots_dir, plot_filename)
        
        plot_best_fit(protein_name, df_plot_long, best_peak, mdf, plot_filepath)
        
    print(f"Plots saved in: {plots_dir}")