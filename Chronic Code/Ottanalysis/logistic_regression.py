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


    ######### EXTRACT FEATURES ######### 
    print(f"Best Regularization Strength (C): {clf.C_[0]}")
    print(f"Best L1 Ratio: {clf.l1_ratio_[0]}")
    
    if len(classes) == 2:
        coefs = clf.coef_.flatten()
    else:
        coefs = np.max(np.abs(clf.coef_), axis=0)

    mask = coefs != 0
    current_protein_names = X.columns
    selected_proteins = np.array(current_protein_names)[mask]
    selected_coefs = coefs[mask]
    
    print(f"Selected {len(selected_proteins)} proteins.")

    logistic_df = pd.DataFrame({
        'Protein': selected_proteins,
        'Coefficient': selected_coefs,
        'Abs_Coefficient': np.abs(selected_coefs)
    }).sort_values(by='Abs_Coefficient', ascending=False)
    
    file_prefix = f"LogReg_{compare_name}_vs_{baseline_name}"
    output_file_path = os.path.join(output_dir, f"{file_prefix}_Features.csv")
    logistic_df.to_csv(output_file_path, index=False)
    print(f"Saved features to: {output_file_path}")


    ######### VISUALISATION ###########
    pdf_path = os.path.join(output_dir, f"{file_prefix}_Plots.pdf")
    with PdfPages(pdf_path) as pdf:

        # Plot A: ROC Curve (Using Unbiased Probabilities)
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

        # Plot B: Top Coefficients
        if not logistic_df.empty:
            top_n = 20
            plot_df = logistic_df.head(top_n)
            
            plt.figure(figsize=(10, 8))
            # FIX: Added hue and legend=False
            sns.barplot(data=plot_df, x='Coefficient', y='Protein', hue='Protein', palette='vlag', legend=False)
            plt.title(f'Top {top_n} Logistic Coefficients (Model fitted on all data)')
            plt.xlabel('Coefficient Value (Log Odds)')
            plt.axvline(0, color='k', linewidth=1)
            plt.grid(axis='x', alpha=0.3)
            plt.tight_layout()
            pdf.savefig()
            plt.close()
        
    print(f"Saved plots to: {pdf_path}")