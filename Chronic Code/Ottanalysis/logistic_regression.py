import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import warnings
from sklearn.exceptions import ConvergenceWarning

from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import GridSearchCV, LeaveOneOut, permutation_test_score
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_curve, auc, confusion_matrix, accuracy_score, f1_score
from sklearn.base import clone
from sklearn.linear_model import LogisticRegression
from sklearn.utils import resample

from logistic_preprocessing import NoEligibleProteins, ProteinFeatures

##########################################################################
# Version details for reproducibility/citations
print("\n--- Pandas and sklearn Version Information ---")

pandas_version = pd.__version__
print(f"pandas version: {pandas_version}")

import sklearn
sklearn_version = sklearn.__version__
print(f"scikit-learn version: {sklearn_version}")
########################################################################



# Suppress the max_iter ConvergenceWarnings from scikit-learn
warnings.filterwarnings("ignore", category=ConvergenceWarning)

# Preserve the original search values and tie order (ascending C within each L1 ratio).
_PARAMETER_GRID = [
    {"classifier__C": np.logspace(-4, 4, 10), "classifier__l1_ratio": [ratio]}
    for ratio in (0.1, 0.5, 0.7, 0.9, 0.95, 1.0)
]


def _fit_tuned_model(X, y, delta):
    pipeline = Pipeline([
        ("proteins", ProteinFeatures(delta=delta)),
        ("scale", StandardScaler()),
        ("classifier", LogisticRegression(
            penalty="elasticnet", solver="saga", max_iter=10000, random_state=42
        )),
    ])
    # GridSearchCV must wrap the whole pipeline so even inner validation rows
    # are excluded from selection, imputation and scaling.
    folds = min(5, int(np.bincount(y).min()))
    search = GridSearchCV(pipeline, _PARAMETER_GRID, cv=folds, scoring="accuracy",
                          n_jobs=-1, error_score="raise")
    return search.fit(X, y).best_estimator_


def _leave_one_out_probabilities(X, y, delta):
    probabilities = np.empty((len(X), 2))
    for i, (train, test) in enumerate(LeaveOneOut().split(X), start=1):
        filled = 40 * i // len(X)
        bar = '█' * filled + '-' * (40 - filled)
        print(f'\rProgress: |{bar}| {100 * i // len(X)}% ({i}/{len(X)})', end='', flush=True)
        model = _fit_tuned_model(X.iloc[train], y[train], delta)
        probabilities[test] = model.predict_proba(X.iloc[test])
    print("\nValidation Complete.\n")
    return probabilities


def run_logistic_regression(df_subset, group_col, protein_cols, baseline_name, compare_name, output_dir, *, delta_visits=None):
    print("\n------------------------------------------------------------------")
    print(f"Running Logistic Regression avec: {compare_name} vs {baseline_name}")
    print("------------------------------------------------------------------")
    
    y = df_subset[group_col]
    delta = delta_visits is not None
    if delta:
        if not delta_visits.index.equals(y.index) or not y.index.is_unique:
            raise ValueError("Delta endpoints and labels must have the same unique participant index.")
        X = delta_visits.loc[:, pd.IndexSlice[:, protein_cols]].replace([np.inf, -np.inf], np.nan)
        has_visits = X["baseline"].notna().any(axis=1) & X["followup"].notna().any(axis=1)
        if not has_visits.all():
            print(f"Excluding {int((~has_visits).sum())} participants without measurements at both visits.")
        X = X.loc[has_visits]
        y = y.loc[has_visits]
    else:
        X = df_subset[protein_cols].replace([np.inf, -np.inf], np.nan)

    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    classes = le.classes_
    print(f"Classes: {classes}")
    if len(classes) != 2 or np.bincount(y_encoded).min() < 3:
        print("Skipping logistic regression: nested leave-one-out validation requires at least three rows in each of two classes.")
        return

    print("Performing leave-one-out validation with preprocessing and tuning inside each training set...")
    try:
        probabilities = _leave_one_out_probabilities(X, y_encoded, delta)
        # The all-data refit supplies the reported coefficients and bootstrap feature space.
        full_model = _fit_tuned_model(X, y_encoded, delta)
    except NoEligibleProteins as exc:
        print(f"\nSkipping this logistic comparison: {exc}")
        return

    y_pred_cv = probabilities.argmax(axis=1)
    y_probs_cv = probabilities[:, 1]
    print("\n\n--- Model Performance (Nested CV Estimate) ---")
    print(f"Accuracy: {accuracy_score(y_encoded, y_pred_cv):.4f}")
    print(f"F1 Score: {f1_score(y_encoded, y_pred_cv):.4f}")
    print(f"Confusion Matrix:\n{confusion_matrix(y_encoded, y_pred_cv)}\n")
    clf = full_model.named_steps["classifier"]
    current_protein_names = full_model.named_steps["proteins"].proteins_
    X_full_scaled = full_model[:-1].transform(X)
    print(f"Proteins retained in the final fit: {len(current_protein_names)}")
    print(f"Best Regularisation Strength (C): {clf.C}")
    print(f"Best L1 Ratio: {clf.l1_ratio}")
    

    print("\nRunning bootstrapping for feature stability (calculating confidence intervals)")
    n_iterations = 500
    boot_model = clone(clf).set_params(max_iter=5000)

    bootstrapped_coefs = []
    
    for i in range(n_iterations):
        X_boot, y_boot = resample(X_full_scaled, y_encoded, random_state=i, stratify=y_encoded)
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
    original_coefs = clf.coef_.ravel()
    mask = original_coefs != 0
    
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


    # Hyperparameters remain fixed for this diagnostic; each fold and each
    # permutation refits the preprocessing using its own training labels.
    print("\nPermutation check at the final hyperparameters (tuning is not repeated)")
    valid_permutation = False
    try:
        score, permutation_scores, pvalue = permutation_test_score(
            full_model, X, y_encoded, scoring="roc_auc", cv=min(5, int(np.bincount(y_encoded).min())),
            n_permutations=100, n_jobs=-1, random_state=42
        )
        valid_permutation = np.isfinite(score) and np.isfinite(permutation_scores).all()
    except NoEligibleProteins as exc:
        print(f"Permutation result unavailable: {exc}")
    else:
        if not valid_permutation:
            print("Permutation result unavailable: at least one AUC was undefined.")
    if valid_permutation:
        print(f"Fixed-parameter CV AUC: {score:.4f}")
        print(f"Fixed-parameter permutation p-value: {pvalue:.4f}")


    ######### VISUALISATION ###########
    pdf_path = os.path.join(output_dir, f"{file_prefix}_Plots.pdf")
    with PdfPages(pdf_path) as pdf:

        # Plot 1: ROC curve
        if len(classes) == 2:
            fpr, tpr, _ = roc_curve(y_encoded, y_probs_cv)
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
        if valid_permutation:
            plt.figure(figsize=(8, 6))
            plt.hist(permutation_scores, bins=20, density=True, alpha=0.7, color='gray', label='Null Dist (Permuted)')
            plt.axvline(score, color='red', linestyle='dashed', linewidth=2, label=f'Real Model AUC = {score:.2f}\np-value = {pvalue:.4f}')
            plt.title('Permutation Test at Fixed Hyperparameters')
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
            
            y_pos = np.arange(len(plot_df))
            
            # Center the error bars on the Bootstrapped Mean (using np.abs to guarantee positive values)
            x_err = [
                np.abs(plot_df['Bootstrapped Mean'] - plot_df['95% CI Lower']),
                np.abs(plot_df['95% CI Upper'] - plot_df['Bootstrapped Mean'])
            ]
            
            # Draw the bars representing original coefficient
            plt.barh(y_pos, plot_df['Coefficient (Original)'], 
                     color=sns.color_palette('vlag', n_colors=len(plot_df)), 
                     label='Original Coefficient')
            
            # Draw the error bars centered on the Bootstrapped Mean (with a dot for the mean)
            plt.errorbar(plot_df['Bootstrapped Mean'], y_pos, xerr=x_err, 
                         fmt='o', color='black', ecolor='black', capsize=4, markersize=5, 
                         label='Bootstrap Mean & 95% CI')
            
            plt.yticks(y_pos, plot_df['Protein'])
            plt.gca().invert_yaxis()  # sort with highest absolute value at the top
            
            plt.title(f'Top {top_n} Features (Bars = Original, Points/Lines = Bootstrap 95% CI)')
            plt.xlabel('Coefficient Value (Log Odds)')
            plt.axvline(0, color='k', linewidth=1)
            plt.legend(loc='lower right')
            plt.grid(axis='x', alpha=0.3)
            plt.tight_layout()
            pdf.savefig()
            plt.close()
        
    print(f"Saved plots to: {pdf_path}")
