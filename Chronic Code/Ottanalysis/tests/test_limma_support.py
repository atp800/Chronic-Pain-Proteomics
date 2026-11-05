import io
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
import pandas as pd
from scipy.stats import t

import limma


class LimmaSupportTests(unittest.TestCase):
    def setUp(self):
        self.enterContext(redirect_stdout(io.StringIO()))
        self.folder = Path(self.enterContext(TemporaryDirectory()))
        self.config = dict(ID_COLUMN='Sample_ID', ID_DELIMITER='-', COVARIATE_COLUMNS=[],
                           P_THRESHOLD=.05, LOG_FC_THRESHOLD=.58)
        self.data = pd.DataFrame({
            'Sample_ID': [f'p{i}-D14' for i in range(8)], 'level': ['ref'] * 4 + ['case'] * 4,
            'batch': [0, 0, 1, 1, 0, 0, 1, 1], 'order': [0, 1, 0, 1, 1, 0, 1, 0],
            'valid': [9., 10, 12, 11, 17, 15, 16, 18],
            'partial': [9., 10, 12, 11, 17, 15, np.nan, np.nan],
        })

    def run_comparison(self, data, proteins, *, paired=False, delta=False, covariates=()):
        original = data.copy(deep=True)
        self.config['COVARIATE_COLUMNS'] = list(covariates)
        limma.run_limma(data, 'level', proteins, 'ref', 'case', self.folder, paired, delta, self.config)
        pd.testing.assert_frame_equal(data, original)
        return pd.read_excel(self.folder / 'Limma_ref_vs_case_Results.xlsx', sheet_name=None)

    def test_paired_counts_separate_observed_values_from_complete_pairs_used(self):
        data = pd.DataFrame({
            'Sample_ID': ['a-ref', 'a-case', 'b-ref', 'b-case', 'c-ref', 'c-case', 'd-ref', 'd-case'],
            'level': ['ref', 'case'] * 4,
            'partial': [10., 11, 19, 22, 90, np.nan, np.nan, 80],
            'no_pairs': [10., np.nan, np.nan, 25, 10, np.nan, np.nan, 25],
        }).sample(frac=1, random_state=12)
        sheets = self.run_comparison(data, ['partial', 'no_pairs'], paired=True)
        row = sheets['All_Proteins'].set_index('Protein').loc['partial']
        for group in ['Baseline', 'Comparison']:
            self.assertEqual(row[f'N_Observed_{group}'], 3)
            self.assertEqual(row[f'N_Participants_Observed_{group}'], 3)
            self.assertEqual(row[f'N_Fit_Observations_{group}'], 2)
            self.assertEqual(row[f'Missing_Proportion_{group}'], .25)
            self.assertEqual(row[f'Protein_Imputed_Proportion_{group}'], 0)
        self.assertEqual(row.N_Complete_Pairs, 2)
        self.assertEqual(row.Residual_DF, 1)
        self.assertAlmostEqual(row.logFC, 2)
        excluded = sheets['Protein_Support'].set_index('Protein').loc['no_pairs']
        self.assertEqual(excluded.N_Complete_Pairs, 0)
        self.assertEqual(excluded.N_Fit_Observations_Baseline, 0)
        self.assertTrue(np.isnan(excluded.Residual_DF))
        self.assertIn('Fewer than two complete pairs', excluded.Status)

    def test_independent_counts_and_intervals_follow_protein_ids_after_sorting(self):
        self.data.loc[0, 'partial'] = 0  # An observed zero is not missing.
        self.data.loc[6, 'partial'] = np.inf  # Non-finite outcomes are missing.
        sheets = self.run_comparison(self.data.sample(frac=1, random_state=4), ['partial', 'valid'])
        results = sheets['All_Proteins'].set_index('Protein')
        self.assertEqual(results.loc['partial', 'N_Observed_Baseline'], 4)
        self.assertEqual(results.loc['partial', 'N_Observed_Comparison'], 2)
        self.assertEqual(results.loc['partial', 'N_Participants_Observed_Comparison'], 2)
        self.assertEqual(results.loc['partial', 'Residual_DF'], 4)
        self.assertEqual(results.loc['valid', 'Residual_DF'], 6)
        self.assertTrue(results.N_Complete_Pairs.isna().all())
        self.assertEqual(results.loc['partial', 'Missing_Proportion_Comparison'], .5)
        with limma.localconverter(limma.ro.default_converter + limma.pandas2ri.converter):
            expected = limma.ro.conversion.rpy2py(limma.ro.r('''data.frame(
                Protein=fit$genes$Protein, Effect=fit$coefficients[, "comparison"],
                SE=fit$stdev.unscaled[, "comparison"] * sqrt(fit$s2.post),
                DF=fit$df.total)''')).set_index('Protein')
            original_stats = limma.ro.conversion.rpy2py(limma.ro.r(
                'topTable(fit, coef="comparison", number=Inf, sort.by="P", confint=FALSE)')).set_index('Protein')
        expected = expected.loc[results.index]
        margin = t.ppf(.975, expected.DF) * expected.SE
        np.testing.assert_allclose(results['CI.L'], expected.Effect - margin)
        np.testing.assert_allclose(results['CI.R'], expected.Effect + margin)
        np.testing.assert_allclose(results.Moderated_DF, expected.DF)
        for col in ['logFC', 'P.Value', 'adj.P.Val', 'AveExpr', 't', 'B']:
            np.testing.assert_allclose(results[col], original_stats.loc[results.index, col])
        for name in ['Significant_Raw_PVal', 'Significant_Adj_PVal']:
            subset = sheets[name].set_index('Protein')
            pd.testing.assert_frame_equal(subset, results.loc[subset.index], check_dtype=False)
        self.assertIn('Not simultaneous', sheets['Analysis_Info'].set_index('Setting').loc['CI.L / CI.R', 'Value'])

    def test_delta_counts_use_complete_endpoint_pairs_and_unsplit_participant_ids(self):
        data = pd.DataFrame({'Sample_ID': [f'subject-{i}' for i in range(6)],
                             'level': ['ref'] * 3 + ['case'] * 3,
                             'change': [1., 2, np.nan, 7, 8, 10]})
        sheets = self.run_comparison(data, ['change'], delta=True)
        row = sheets['All_Proteins'].iloc[0]
        self.assertEqual(row.N_Participants_Observed_Baseline, 2)
        self.assertEqual(row.N_Participants_Observed_Comparison, 3)
        self.assertEqual(row.N_Complete_Pairs, 5)
        self.assertEqual(row.Residual_DF, 3)
        self.assertAlmostEqual(row.logFC, (7 + 8 + 10) / 3 - 1.5)
        self.assertEqual(row.Protein_Imputed_Proportion_Baseline, 0)

    def test_covariate_df_and_unestimable_or_saturated_proteins_remain_visible(self):
        self.data['confounded'] = [9., 10, np.nan, np.nan, np.nan, np.nan, 16, 18]
        self.data['saturated'] = [9., 10, np.nan, np.nan, 17, np.nan, 16, np.nan]
        sheets = self.run_comparison(self.data, ['valid', 'confounded', 'saturated'], covariates=['batch', 'order'])
        support = sheets['Protein_Support'].set_index('Protein')
        self.assertEqual(support.loc['valid', 'Residual_DF'], 4)
        self.assertEqual(support.loc['confounded', 'Status'], 'Comparison effect is not estimable')
        self.assertEqual(support.loc['saturated', 'Residual_DF'], 0)
        self.assertEqual(support.loc['saturated', 'Status'], 'No residual degrees of freedom')
        self.assertTrue(support.loc[['confounded', 'saturated'], 'Moderated_DF'].isna().all())
        for name in ['All_Proteins', 'Cov_Effect_batch', 'Cov_Effect_order']:
            row = sheets[name].iloc[0]
            self.assertEqual(row.Protein, 'valid')
            self.assertEqual(row.Residual_DF, 4)
            self.assertTrue(row['CI.L'] < row.logFC < row['CI.R'])

    def test_no_eligible_or_testable_proteins_still_export_support(self):
        for label in ['none_observed', 'constant', 'confounded']:
            with self.subTest(label=label):
                data = self.data.copy()
                data['none_observed'] = np.nan
                data['constant'] = 5.
                data['confounded'] = [9., 10, np.nan, np.nan, np.nan, np.nan, 16, 18]
                sheets = self.run_comparison(data, [label], covariates=['batch'])
                for name in ['All_Proteins', 'Significant_Raw_PVal', 'Significant_Adj_PVal']:
                    self.assertTrue(sheets[name].empty)
                    self.assertIn('CI.L', sheets[name])
                    self.assertIn('N_Participants_Observed_Baseline', sheets[name])
                self.assertEqual(sheets['Protein_Support'].Protein.tolist(), [label])
                self.assertNotIn('Eligible for fitting', sheets['Protein_Support'].Status.tolist())

    def test_distinct_participants_are_not_inferred_from_row_counts_or_missing_ids(self):
        self.data.loc[1, 'Sample_ID'] = self.data.loc[0, 'Sample_ID']
        self.data.loc[4, 'Sample_ID'] = None
        sheets = self.run_comparison(self.data, ['valid'])
        row = sheets['All_Proteins'].iloc[0]
        self.assertEqual(row.N_Observed_Baseline, 4)
        self.assertEqual(row.N_Participants_Observed_Baseline, 3)
        self.assertTrue(np.isnan(row.N_Participants_Observed_Comparison))
        self.assertEqual(row.N_Observed_Comparison, 4)
