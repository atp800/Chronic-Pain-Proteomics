"""Reporting only: describe the observed data and the fits actually used by limma."""

import numpy as np
import pandas as pd


RESULT_COLUMNS = ['Protein', 'logFC', 'CI.L', 'CI.R', 'P.Value', 'adj.P.Val', 'AveExpr', 'B', 't',
                  'Residual_DF', 'Moderated_DF']


def protein_support(observed, fitted_values, groups, subjects, baseline, comparison,
                    replicated, keep, *, paired, delta, complete_pairs=None):
    """One row per requested protein, including proteins excluded before fitting."""
    support = pd.DataFrame(index=observed.columns)
    support['Baseline_Group'] = baseline
    support['Comparison_Group'] = comparison
    for label, level in [('Baseline', baseline), ('Comparison', comparison)]:
        rows = groups.eq(level)
        values = observed.loc[rows]
        ids = subjects.loc[rows]
        support[f'N_Observed_{label}'] = values.notna().sum()
        support[f'N_Participants_Observed_{label}'] = [
            np.nan if ids[values[p].notna()].isna().any()
            else ids[values[p].notna()].nunique() for p in observed.columns
        ]
        support[f'N_Fit_Observations_{label}'] = fitted_values.loc[rows].notna().sum().where(keep, 0)
        support[f'Missing_Proportion_{label}'] = values.isna().mean()
        # This describes protein outcomes, not the separate handling of covariates.
        support[f'Protein_Imputed_Proportion_{label}'] = 0.
    support['N_Complete_Pairs'] = np.nan
    if paired:
        support['N_Complete_Pairs'] = complete_pairs
    elif delta:
        # Every finite input delta comes from both observed endpoints in the main pipeline.
        support['N_Complete_Pairs'] = [
            np.nan if subjects[observed[p].notna()].isna().any()
            else subjects[observed[p].notna()].nunique() for p in observed.columns
        ]
    support['Residual_DF'] = np.nan
    support['Moderated_DF'] = np.nan
    support['Status'] = 'Eligible for fitting'
    support.loc[~replicated, 'Status'] = ('Fewer than two complete pairs' if paired
                                         else 'Fewer than two observed values in at least one group')
    support.loc[replicated & ~keep, 'Status'] = 'No observed variation'
    support.index = support.index.astype(str)
    support.index.name = 'Protein'
    return support.reset_index()


def attach_support(results, support):
    """Join by protein ID, never by row position after topTable sorts the results."""
    columns = [c for c in support if c == 'Protein' or c not in results]
    return results.merge(support[columns], on='Protein', how='left', validate='one_to_one', sort=False)


def write_results(path, results, covariate_results, support, info, config):
    """Keep existing result sheets and add a complete support/eligibility audit."""
    final_cols = [c for c in RESULT_COLUMNS if c in results]
    final_cols += [c for c in support if c not in final_cols and c in results]
    raw = results['P.Value'].lt(config['P_THRESHOLD']) & results['logFC'].abs().gt(config['LOG_FC_THRESHOLD'])
    adjusted = results['adj.P.Val'].lt(config['P_THRESHOLD']) & results['logFC'].abs().gt(config['LOG_FC_THRESHOLD'])
    print(f"Significant proteins (raw p-value < {config['P_THRESHOLD']} & abs(logFC) > {config['LOG_FC_THRESHOLD']}): {int(raw.sum())}")
    print(f"Significant proteins (adjusted p-value < {config['P_THRESHOLD']} & abs(logFC) > {config['LOG_FC_THRESHOLD']}): {int(adjusted.sum())}")
    with pd.ExcelWriter(path, engine='openpyxl') as writer:
        results[final_cols].sort_values('adj.P.Val').to_excel(writer, sheet_name='All_Proteins', index=False)
        results.loc[raw, final_cols].sort_values('P.Value').to_excel(writer, sheet_name='Significant_Raw_PVal', index=False)
        results.loc[adjusted, final_cols].sort_values('adj.P.Val').to_excel(writer, sheet_name='Significant_Adj_PVal', index=False)
        for name, table in covariate_results.items():
            cols = [c for c in final_cols if c in table]
            table[cols].to_excel(writer, sheet_name=f"Cov_Effect_{name.replace('cov_', '')}"[:31], index=False)
        support.to_excel(writer, sheet_name='Protein_Support', index=False)
        info.to_excel(writer, sheet_name='Analysis_Info', index=False)
