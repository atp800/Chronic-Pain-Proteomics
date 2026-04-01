import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import LeaveOneOut
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import roc_curve, auc, confusion_matrix, accuracy_score, f1_score
from sklearn.base import clone
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import permutation_test_score
from sklearn.utils import resample

from helpers import groupwise_missing_filter, impute_missing_vals

def run_logistic_regression(df_subset, group_col, protein_cols, baseline_name, compare_name, output_dir):
    print("\n------------------------------------------------------------------")
    print(f"Running Logistic Regression avec: {compare_name} vs {baseline_name}")
    print("------------------------------------------------------------------")
    
    ########## PREPPARE DATA ############
    X = df_subset[protein_cols]
    y = df_subset[group_col]

    # More aggressive group-based missing protein filtering     -   REMOVE OR REMPLACE WITH EARLIER IMPUTATION???
    print(f"Original dimensionality: {X.shape[1]} proteins")
    X = groupwise_missing_filter(X, y, threshold=0.7)               # Keep proteins with <30% missing in at least one group
    print(f"Proteins after dropping >30% missing: {X.shape[1]}")
    
    # Remove low-variance features
    X_temp_for_var = X.fillna(0)                                    # Temporarily fill NaNs for variance calculation
    selector = VarianceThreshold(threshold=0.1) 
    try:
        selector.fit(X_temp_for_var)
        cols_kept = X.columns[selector.get_support()]
        X = X[cols_kept]
        print(f"Proteins after dropping low variance: {X.shape[1]}")
    except ValueError:
        print("Warning: Variance filter failed (variance too low globally?), skipping.")

    # Check for NaNs and impute
    if X.isna().sum().sum() > 0:
        print("Missing values detected: imputing with down-shifted normal distribution")
        X = impute_missing_vals(X)

    # Encode target
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    classes = le.classes_
    print(f"Classes: {classes}")

    if len(classes) != 2:
        print(f"Warning: {len(classes)} groups detected. This script is designed for 2 groups.")

    # Scaling
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    ########### TRAIN MODEL!!! ############
    print("Training Logistic Regression (sometimes it takes its sweet sweet time, but usually <5 mins)...")
    
    # Define the model structure (Nested CV: The model itself does CV to find C)
    clf = LogisticRegressionCV(
        Cs=10, 
        cv=5, 
        penalty='elasticnet', 
        solver='saga', 
        l1_ratios=[0.1, 0.5, 0.7, 0.9, 0.95, 1],
        random_state=42, 
        max_iter=10000, 
        scoring='accuracy',
        verbose=0,                  # increase for more print statements while running
        n_jobs=-1
    )

    print("Performing leave-one-out cross-validation...")
    clf.fit(X_scaled, y_encoded)
    cv_outer = LeaveOneOut()            # Simulates a test set for every single sample
    
    # Lists to store results
    y_pred_unbiased = []
    y_probs_unbiased =[]
    total_samples = len(X_scaled)
    
    # Loop through every sample
    for i, (train_index, test_index) in enumerate(cv_outer.split(X_scaled)):
        # Progress bar
        percent = (i + 1) / total_samples
        bar_length = 40
        filled_length = int(bar_length * percent)
        bar = '█' * filled_length + '-' * (bar_length - filled_length)
        sys.stdout.write(f'\rProgress: |{bar}| {int(percent * 100)}% ({i+1}/{total_samples})')      # \r returns cursor to start of line, allowing overwrite
        sys.stdout.flush()

        # Split Data
        X_train, X_test = X_scaled[train_index], X_scaled[test_index]
        y_train, y_test = y_encoded[train_index], y_encoded[test_index]
        
        # Clone the model (start fresh every time)
        model_fold = clone(clf) 
        model_fold.fit(X_train, y_train)
        
        # Predict
        prediction = model_fold.predict(X_test)[0]
        y_pred_unbiased.append(prediction)
        
        if len(classes) == 2:
            prob = model_fold.predict_proba(X_test)[0][1]
            y_probs_unbiased.append(prob)

    print("\nValidation Complete.\n")
    
    # Convert lists to numpy arrays for the metrics calculations
    y_pred_unbiased = np.array(y_pred_unbiased)
    if len(classes) == 2:
        y_probs_unbiased = np.array(y_probs_unbiased)

    print("\n\n--- Model Performance (Nested CV Estimate) ---")
    print(f"Accuracy: {accuracy_score(y_encoded, y_pred_unbiased):.4f}")
    print(f"F1 Score: {f1_score(y_encoded, y_pred_unbiased):.4f}")
    print(f"Confusion Matrix:\n{confusion_matrix(y_encoded, y_pred_unbiased)}\n")



    ######### EXTRACT FEATURES (and bootstrapping for feature stability ######### 
    best_C = clf.C_[0]
    if isinstance(best_C, np.ndarray): 
        best_C = best_C[0]
    best_l1 = clf.l1_ratio_[0]

    print(f"Best Regularistion Strength (C): {best_C}")
    print(f"Best L1 Ratio: {best_l1}")
    

    print("\nRunning bootstrapping for feature stability (calculating confidence intervals)")
    n_iterations = 500      # add to config??? icnrease to 1000 for more robust estimates???




    # Create a fixed model using the optimal hyperparameters found earlier to save time
    boot_model = LogisticRegression(
        C=best_C, 
        l1_ratio=best_l1, 
        penalty='elasticnet', 
        solver='saga', 
        max_iter=5000,
        random_state=42
    )

    bootstrapped_coefs = []
    
    for i in range(n_iterations):
        # Sample with replacement
        X_boot, y_boot = resample(X_scaled, y_encoded, random_state=i)
        boot_model.fit(X_boot, y_boot)
        bootstrapped_coefs.append(boot_model.coef_.flatten())

    bootstrapped_coefs = np.array(bootstrapped_coefs)
    
    # Calculate statistics
    mean_coefs = np.mean(bootstrapped_coefs, axis=0)
    lower_bounds = np.percentile(bootstrapped_coefs, 2.5, axis=0)  # 95% CI Lower
    upper_bounds = np.percentile(bootstrapped_coefs, 97.5, axis=0) # 95% CI Upper
    
    # Selection frequency: fraction of times coefficient is exactly not zero
    selection_freq = np.mean(bootstrapped_coefs != 0, axis=0) * 100

    # Determine "Significance" (CI does not cross zero AND selection > 0 in original model)
    original_coefs = clf.coef_.flatten()
    mask = original_coefs != 0
    current_protein_names = X.columns
    
    selected_proteins = np.array(current_protein_names)[mask]
    selected_coefs = original_coefs[mask]
    sel_mean = mean_coefs[mask]
    sel_lower = lower_bounds[mask]
    sel_upper = upper_bounds[mask]
    sel_freq = selection_freq[mask]
    
    # A feature/protein is "significant" if both limits of the CI are >0 or both are <0
    significant_ci = (sel_lower > 0) & (sel_upper > 0) | (sel_lower < 0) & (sel_upper < 0)

    print(f"Selected {len(selected_proteins)} proteins with non-zero coefficients.")

    logistic_df = pd.DataFrame({
        'Protein': selected_proteins,
        'Coefficient (Original)': selected_coefs,
        'Bootstrapped Mean': sel_mean,
        '95% CI Lower': sel_lower,
        '95% CI Upper': sel_upper,
        'Selection Frequency (%)': sel_freq,
        'Significant CI': significant_ci,
        'Abs_Coefficient': np.abs(selected_coefs)
    }).sort_values(by='Abs_Coefficient', ascending=False)
    
    file_prefix = f"LogReg_{compare_name}_vs_{baseline_name}"
    output_file_path = os.path.join(output_dir, f"{file_prefix}_Features.csv")
    logistic_df.drop(columns=['Abs_Coefficient']).to_csv(output_file_path, index=False)
    print(f"Saved features and stability metrics to: {output_file_path}")


    ######## OVERALL MODEL SIGNIFICANCE (permutation test)##########
    print("\nRunning permutation test to calculate P-Value for the model as a whole")       # Shuffles the Y labels to prove the model is significantly better than random guessing
    score, permutation_scores, pvalue = permutation_test_score(
        boot_model, X_scaled, y_encoded, scoring="roc_auc", cv=5, n_permutations=100, n_jobs=-1, random_state=42
    )
    print(f"Model Original AUC: {score:.4f}")
    print(f"Model Empirical P-value: {pvalue:.4f}")


    ######### VISUALISATION ###########
    pdf_path = os.path.join(output_dir, f"{file_prefix}_Plots.pdf")
    with PdfPages(pdf_path) as pdf:

        # Plot 1: ROC curve
        if len(classes) == 2:
            fpr, tpr, _ = roc_curve(y_encoded, y_probs_unbiased)
            roc_auc = auc(fpr, tpr)
            
            plt.figure(figsize=(8, 6))
            plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (CV AUC = {roc_auc:.2f})')
            plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title(f'ROC Curve (Nested CV): {classes[1]} vs {classes[0]}')
            plt.legend(loc="lower right")
            plt.grid(True, alpha=0.3)
            pdf.savefig()
            plt.close()

        # Plot II: Permutation distribution plot
        if len(classes) == 2:
            plt.figure(figsize=(8, 6))
            plt.hist(permutation_scores, bins=20, density=True, alpha=0.7, color='gray', label='Null Dist (Permuted)')
            plt.axvline(score, color='red', linestyle='dashed', linewidth=2, label=f'Real Model AUC = {score:.2f}\np-value = {pvalue:.4f}')
            plt.title('Overall Model Significance (Target Permutation Test)')
            plt.xlabel('Cross-Validated AUC')
            plt.ylabel('Frequency')
            plt.legend(loc='best')
            plt.grid(True, alpha=0.3)
            pdf.savefig()
            plt.close()

        # Plot C: Top coefficients with error bars
        if not logistic_df.empty:
            top_n = min(20, len(logistic_df))
            plot_df = logistic_df.head(top_n)
            
            plt.figure(figsize=(10, 8))
            
            # Use error bars to show 95% confidence intervals
            y_pos = np.arange(len(plot_df))
            x_err = [
                plot_df['Coefficient (Original)'] - plot_df['95% CI Lower'],
                plot_df['95% CI Upper'] - plot_df['Coefficient (Original)']
            ]
            
            plt.barh(y_pos, plot_df['Coefficient (Original)'], color=sns.color_palette('vlag', n_colors=len(plot_df)))
            plt.errorbar(plot_df['Coefficient (Original)'], y_pos, xerr=x_err, fmt='none', ecolor='black', capsize=4)
            
            plt.yticks(y_pos, plot_df['Protein'])
            plt.gca().invert_yaxis()  # sort with highest absolute value at the top
            
            plt.title(f'Top {top_n} Features with 95% Bootstrap CIs')
            plt.xlabel('Coefficient Value (Log Odds)')
            plt.axvline(0, color='k', linewidth=1)
            plt.grid(axis='x', alpha=0.3)
            plt.tight_layout()
            pdf.savefig()
            plt.close()
        
    print(f"Saved plots to: {pdf_path}")