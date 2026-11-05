import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import warnings
from sklearn.exceptions import ConvergenceWarning

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, GroupKFold, LeaveOneGroupOut, LeaveOneOut
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_curve, roc_auc_score, confusion_matrix, accuracy_score, f1_score
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


def _encode_comparison(labels, baseline_name, compare_name):
    """Match the caller's comparison order, including labels selected as strings by the GUI."""
    if pd.isna(baseline_name) or pd.isna(compare_name):
        raise ValueError("Both baseline and comparison labels are required.")
    classes = np.asarray([str(baseline_name), str(compare_name)])
    if classes[0] == classes[1]:
        raise ValueError("Baseline and comparison labels must be distinct.")
    if labels.isna().any() or not labels.astype(str).isin(classes).all():
        raise ValueError("Outcome labels must contain only the requested baseline and comparison, without missing labels.")
    return labels.astype(str).eq(classes[1]).to_numpy(dtype=int), classes


def _complete_pair_mask(X, visits, participant_ids):
    if not participant_ids.index.equals(X.index):
        raise ValueError("Participant IDs and measurements must have the same row index.")
    known = participant_ids.notna() & participant_ids.astype(str).str.strip().ne("")
    if not known.all():
        print(f"Excluding {int((~known).sum())} rows without a participant ID.")
    if pd.MultiIndex.from_arrays([participant_ids[known], visits[known]]).has_duplicates:
        raise ValueError("Duplicate participant/visit. Check IDs or resolve replicate samples.")

    # The caller selects two distinct visits. Require one measured sample at each.
    observed = known & X.notna().any(axis=1)
    counts = participant_ids[observed].value_counts()
    complete = counts.index[counts.eq(2)]
    incomplete = participant_ids[known].nunique() - len(complete)
    if incomplete:
        print(f"Excluding {incomplete} participants without measured samples at both requested visits.")
    return observed & participant_ids.isin(complete)


def _validation_cv(y, participant_ids):
    if participant_ids is None:
        return min(5, int(np.bincount(y).min()))
    return GroupKFold(n_splits=min(5, len(np.unique(participant_ids))))


def _fit_tuned_model(X, y, delta, participant_ids=None):
    pipeline = Pipeline([
        ("proteins", ProteinFeatures(delta=delta)),
        ("scale", StandardScaler()),
        ("classifier", LogisticRegression(
            penalty="elasticnet", solver="saga", max_iter=10000, random_state=42
        )),
    ])
    # GridSearchCV must wrap the whole pipeline so even inner validation rows
    # are excluded from selection, imputation and scaling.
    search = GridSearchCV(pipeline, _PARAMETER_GRID, cv=_validation_cv(y, participant_ids), scoring="accuracy",
                          n_jobs=-1, error_score="raise")
    return search.fit(X, y, groups=participant_ids).best_estimator_


def _leave_one_out_probabilities(X, y, delta, participant_ids=None, *, show_progress=True):
    outer = LeaveOneOut() if participant_ids is None else LeaveOneGroupOut()
    n_folds = outer.get_n_splits(X, y, participant_ids)
    probabilities = np.empty((len(X), 2))
    for i, (train, test) in enumerate(outer.split(X, y, participant_ids), start=1):
        if show_progress:
            filled = 40 * i // n_folds
            bar = '█' * filled + '-' * (40 - filled)
            print(f'\rProgress: |{bar}| {100 * i // n_folds}% ({i}/{n_folds})', end='', flush=True)
        train_ids = None if participant_ids is None else participant_ids[train]
        model = _fit_tuned_model(X.iloc[train], y[train], delta, train_ids)
        probabilities[test] = model.predict_proba(X.iloc[test])
    if show_progress:
        print("\nValidation Complete.\n")
    return probabilities


def _participant_pairs(participant_ids):
    """Two input row positions per participant, after checking for complete pairs."""
    return np.argsort(participant_ids, kind="stable").reshape(-1, 2)


def _bootstrap_rows(y, participant_ids, n_iterations):
    if participant_ids is None:
        units, strata = np.arange(len(y)), y
    else:
        units = _participant_pairs(participant_ids)
        strata = None
    for seed in range(n_iterations):
        yield resample(units, replace=True, stratify=strata, random_state=seed).ravel()


def _coefficient_stability_table(protein_names, original_coefs, bootstrap_coefs):
    """Summarise all bootstrap draws for proteins with a nonzero full-data coefficient."""
    selected = original_coefs != 0
    lower, upper = np.percentile(bootstrap_coefs, [2.5, 97.5], axis=0)
    return pd.DataFrame({
        'Protein': np.asarray(protein_names)[selected],
        'Coefficient (Original)': original_coefs[selected],
        'Bootstrapped Mean': bootstrap_coefs.mean(axis=0)[selected],
        'Bootstrap 2.5th Percentile': lower[selected],
        'Bootstrap 97.5th Percentile': upper[selected],
        'Selection Frequency (%)': (bootstrap_coefs != 0).mean(axis=0)[selected] * 100,
    }).sort_values('Coefficient (Original)', key=np.abs, ascending=False)


def _coefficient_stability_figure(table, *, reference_class=None, positive_class=None):
    top = table.head(20)
    y = np.arange(len(top))
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.barh(y, top['Coefficient (Original)'], color=sns.color_palette('vlag', n_colors=len(top)),
            label='Full-data coefficient')
    # A sparse coefficient's mean can lie outside its central percentile range.
    # Draw the observed endpoints directly, including a tick for a zero-width range.
    lower = top['Bootstrap 2.5th Percentile']
    upper = top['Bootstrap 97.5th Percentile']
    ax.hlines(y, lower, upper, color='black', label='Bootstrap 2.5th–97.5th percentiles')
    ax.plot(lower, y, 'k|', upper, y, 'k|')
    ax.scatter(top['Bootstrapped Mean'], y, s=16, color='black', zorder=3, label='Bootstrap mean')
    ax.set_yticks(y, labels=top['Protein'])
    ax.invert_yaxis()
    ax.set_title(f'Coefficient stability (top {len(top)} by absolute full-data coefficient)\n'
                 'Preprocessing and tuning fixed across bootstrap refits')
    ax.set_xlabel('Coefficient (log odds)' if positive_class is None else
                  f'Coefficient (log odds of {positive_class} versus {reference_class})')
    ax.axvline(0, color='k', linewidth=1, zorder=0)
    ax.legend(loc='best')
    ax.grid(axis='x', alpha=0.3)
    fig.tight_layout()
    return fig


def _pooled_auc(y, probabilities):
    if not np.isfinite(probabilities).all():
        return np.nan
    return roc_auc_score(y, probabilities[:, 1])


def _nested_permutation_test(X, y, delta, participant_ids, observed_auc, *, n_permutations=100):
    """Compare the observed pooled held-out AUC with the same procedure under shuffled labels."""
    if not np.isfinite(observed_auc) or n_permutations < 1:
        raise ValueError("A finite observed AUC and at least one permutation are required.")
    rng = np.random.RandomState(42)
    pairs = None if participant_ids is None else _participant_pairs(participant_ids)
    scores = np.empty(n_permutations)
    # All fitting steps are deterministic. Count repeated labelings on every draw,
    # but run the nested search only once per distinct assignment for this comparison.
    scores_by_labels = {y.tobytes(): observed_auc}
    for i in range(n_permutations):
        print(f'\rPermutation {i + 1}/{n_permutations}', end='', flush=True)
        if pairs is None:
            labels = rng.permutation(y)
        else:
            labels = y.copy()
            for pair in pairs:
                labels[pair] = rng.permutation(y[pair])
        key = labels.tobytes()
        if key not in scores_by_labels:
            probabilities = _leave_one_out_probabilities(X, labels, delta, participant_ids, show_progress=False)
            score = _pooled_auc(labels, probabilities)
            if not np.isfinite(score):
                raise ValueError("A permuted held-out AUC was undefined.")
            scores_by_labels[key] = score
        scores[i] = scores_by_labels[key]
    print()

    # Count ties, including equal rank-based AUCs that differ only by roundoff.
    tied = np.isclose(scores, observed_auc, rtol=1e-14, atol=0)
    extreme = np.count_nonzero((scores > observed_auc) | tied)
    return scores, (extreme + 1) / (n_permutations + 1)


def run_logistic_regression(df_subset, group_col, protein_cols, baseline_name, compare_name, output_dir,
                            *, delta_visits=None, participant_ids=None):
    """Model comparison (1) versus baseline (0).

    For longitudinal visits, pass participant IDs indexed like df_subset.
    Each delta row already represents one person.
    """
    print("\n------------------------------------------------------------------")
    print(f"Running Logistic Regression avec: {compare_name} vs {baseline_name}")
    print("------------------------------------------------------------------")
    
    y = df_subset[group_col]
    delta = delta_visits is not None
    paired = participant_ids is not None
    if delta:
        if paired:
            raise ValueError("Each delta already represents one participant; paired visits apply only to longitudinal data.")
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

    if paired:
        requested = {str(baseline_name), str(compare_name)}
        if len(requested) != 2 or set(y.astype(str)) != requested:
            print("Skipping longitudinal logistic regression: both distinct requested visits are required.")
            return
        try:
            keep = _complete_pair_mask(X, y, participant_ids)
        except ValueError as exc:
            print(f"Skipping longitudinal logistic regression: {exc}")
            return
        X, y = X.loc[keep], y.loc[keep]
        participant_ids = participant_ids.loc[keep].to_numpy()
        n_participants = len(np.unique(participant_ids))
        print(f"Longitudinal logistic regression: {n_participants} complete participant pairs ({len(y)} visits).")

    try:
        y_encoded, classes = _encode_comparison(y, baseline_name, compare_name)
    except ValueError as exc:
        print(f"Skipping logistic regression: {exc}")
        return
    print(f"Class mapping: 0 = {classes[0]} (baseline), 1 = {classes[1]} (comparison/positive).")
    print(f"Positive coefficients increase predicted log odds of {classes[1]} versus {classes[0]}.")
    class_counts = np.bincount(y_encoded, minlength=2)
    if class_counts.min() < 3:
        minimum = "three complete pairs" if paired else "three rows in each of two classes"
        print(f"Skipping logistic regression: nested validation requires at least {minimum}.")
        return

    unit = "participant" if paired else "sample"
    cv_label = "Nested participant-held-out CV" if paired else "Nested CV"
    print(f"Leaving one {unit} out at a time, with preprocessing and tuning inside each training set...")
    try:
        probabilities = _leave_one_out_probabilities(X, y_encoded, delta, participant_ids)
        # The all-data refit supplies the reported coefficients and bootstrap feature space.
        full_model = _fit_tuned_model(X, y_encoded, delta, participant_ids)
    except NoEligibleProteins as exc:
        print(f"\nSkipping this logistic comparison: {exc}")
        return

    observed_auc = _pooled_auc(y_encoded, probabilities)
    if not np.isfinite(observed_auc):
        print("Skipping logistic comparison: the observed held-out ROC AUC was undefined.")
        return
    y_pred_cv = probabilities.argmax(axis=1)
    y_probs_cv = probabilities[:, 1]
    accuracy = accuracy_score(y_encoded, y_pred_cv)
    f1 = f1_score(y_encoded, y_pred_cv, pos_label=1, zero_division=0)
    confusion = pd.DataFrame(confusion_matrix(y_encoded, y_pred_cv, labels=[0, 1]),
                             index=classes, columns=classes).rename_axis(index='Actual', columns='Predicted')
    print(f"\n\n--- Model Performance ({cv_label} Estimate) ---")
    print(f"Pooled held-out ROC AUC: {observed_auc:.4f}")
    print(f"Accuracy: {accuracy:.4f}")
    print(f"F1 Score (positive class: {classes[1]}): {f1:.4f}")
    print(f"Confusion Matrix:\n{confusion}\n")
    clf = full_model.named_steps["classifier"]
    current_protein_names = full_model.named_steps["proteins"].proteins_
    X_full_scaled = full_model[:-1].transform(X)
    print(f"Proteins retained in the final fit: {len(current_protein_names)}")
    print(f"Best Regularisation Strength (C): {clf.C}")
    print(f"Best L1 Ratio: {clf.l1_ratio}")
    

    n_iterations = 500
    print(f"\nEstimating coefficient stability across {n_iterations} bootstrap refits with fixed preprocessing and tuning.")
    print("Selection frequency is the percentage of refits with a nonzero coefficient.")
    if paired:
        print("Resampling participants with both visits together.")
    boot_model = clone(clf).set_params(max_iter=5000)

    bootstrapped_coefs = []
    
    for rows in _bootstrap_rows(y_encoded, participant_ids, n_iterations):
        boot_model.fit(X_full_scaled[rows], y_encoded[rows])
        bootstrapped_coefs.append(boot_model.coef_.flatten())

    logistic_df = _coefficient_stability_table(current_protein_names, clf.coef_.ravel(), np.asarray(bootstrapped_coefs))
    logistic_df['Reference Class (0)'] = classes[0]
    logistic_df['Positive Class (1)'] = classes[1]
    print(f"Reporting {len(logistic_df)} proteins with nonzero coefficients in the full-data fit.")
    
    file_prefix = f"LogReg_{compare_name}_vs_{baseline_name}"
    output_file_path = os.path.join(output_dir, f"{file_prefix}_Features.csv")
    logistic_df.to_csv(output_file_path, index=False)
    print(f"Saved features and stability metrics to: {output_file_path}")
    # Persist the direction and F1 target even if the final model selects no proteins.
    performance_path = os.path.join(output_dir, f"{file_prefix}_Performance.csv")
    pd.DataFrame([{
        'Reference Class (0)': classes[0], 'Positive Class (1)': classes[1],
        'Reference Observations': class_counts[0], 'Positive Observations': class_counts[1],
        'Validation': cv_label, 'ROC AUC': observed_auc, 'Accuracy': accuracy,
        'F1 (Positive Class)': f1,
    }]).to_csv(performance_path, index=False)
    print(f"Saved class mapping and performance metrics to: {performance_path}")


    print("\nPermutation test of the pooled held-out ROC AUC, repeating protein selection and tuning in every training set.")
    if paired:
        print("Swapping visit labels only within participants, with participant-held-out validation.")
    permutation_scores = None
    try:
        permutation_scores, pvalue = _nested_permutation_test(X, y_encoded, delta, participant_ids, observed_auc)
    except ValueError as exc:
        print(f"\nPermutation result unavailable: {exc}")
    else:
        print(f"Nested CV permutation p-value ({len(permutation_scores)} draws): {pvalue:.4f}")


    ######### VISUALISATION ###########
    pdf_path = os.path.join(output_dir, f"{file_prefix}_Plots.pdf")
    with PdfPages(pdf_path) as pdf:

        # Plot 1: ROC curve
        if len(classes) == 2:
            fpr, tpr, _ = roc_curve(y_encoded, y_probs_cv)
            
            plt.figure(figsize=(8, 6))
            plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (CV AUC = {observed_auc:.2f})')
            plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title(f'ROC Curve ({cv_label}): {classes[1]} vs {classes[0]}')
            plt.legend(loc="lower right")
            plt.grid(True, alpha=0.3)
            pdf.savefig()
            plt.close()

        # Plot II: Permutation distribution plot
        if permutation_scores is not None:
            plt.figure(figsize=(8, 6))
            plt.hist(permutation_scores, bins=20, density=True, alpha=0.7, color='gray', label='Null Dist (Permuted)')
            plt.axvline(observed_auc, color='red', linestyle='dashed', linewidth=2, label=f'Observed AUC = {observed_auc:.2f}\np-value = {pvalue:.4f}')
            permutation_title = 'Within-Participant Permutation Test' if paired else 'Permutation Test'
            plt.title(f'{permutation_title} with Nested Model Selection')
            plt.xlabel('Pooled Held-Out ROC AUC')
            plt.ylabel('Density')
            plt.legend(loc='best')
            plt.grid(True, alpha=0.3)
            pdf.savefig()
            plt.close()

        if not logistic_df.empty:
            fig = _coefficient_stability_figure(logistic_df, reference_class=classes[0], positive_class=classes[1])
            pdf.savefig(fig)
            plt.close(fig)
        
    print(f"Saved plots to: {pdf_path}")
