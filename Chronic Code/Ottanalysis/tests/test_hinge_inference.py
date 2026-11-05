import io
import runpy
import sys
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd

import hinge_inference as inference
import hinge_model as hinge
import compare_hinge_results as comparison
from config import CONFIG
from test_hinge_fallback import mixed_result, repeated_data


def null_result(**changes):
    result = SimpleNamespace(llf=-10., cov_re=np.array([[9.]]), scale=4.,
                             predict=lambda frame: 3 + .2 * frame.time_numeric.to_numpy())
    for name, value in changes.items():
        setattr(result, name, value)
    return result


class HingeBootstrapTests(unittest.TestCase):
    def setUp(self):
        self.frame = pd.DataFrame({
            'subject_id_col': ['a', 'a', 'b', 'b', 'b', 'c'],
            'time_numeric': [0., 4., 0., 2., 4., 2.],
            'value': [7., 3., 2., 6., 8., 9.],
        }, index=[11, 8, 6, 5, 3, 42])
        self.likelihoods = {2.: -8., 3.: -7.}

    def test_simulation_keeps_observed_rows_and_refits_null_and_every_day(self):
        null = null_result()
        simulated, peaks = [], []
        null_calls = 0

        def fit(frame, *, hinge=True):
            nonlocal null_calls
            pd.testing.assert_series_equal(frame.time_numeric, self.frame.time_numeric)
            pd.testing.assert_series_equal(frame.subject_id_col, self.frame.subject_id_col)
            if not hinge:
                null_calls += 1
                if null_calls > 1:
                    simulated.append(frame.value.to_numpy().copy())
                np.testing.assert_allclose(frame.slope1, frame.time_numeric - frame.time_numeric.mean())
                return null, 'bfgs', 'accepted'
            peak = frame.time_numeric.iloc[0] - frame.slope1.iloc[0]
            peaks.append(peak)
            np.testing.assert_allclose(frame.slope2_change, np.maximum(frame.time_numeric - peak, 0))
            # The simulated winning day differs from the observed winning day.
            return SimpleNamespace(llf=-8. if peak == 2 else -9.), 'bfgs', 'accepted'

        rng = Mock()
        rng.normal.side_effect = [np.array([1., 2., 3.]), np.arange(6.)] * 5
        original = self.frame.copy(deep=True)
        with patch.object(inference.np.random, 'default_rng', return_value=rng):
            result = inference.bootstrap_hinge_test(self.frame, self.likelihoods, fit, protein='P',
                                                    n_bootstrap=5, alpha=.99, show_progress=False)
        self.assertEqual(null_calls, 6)
        self.assertEqual(peaks, [2., 3.] * 5)
        expected = 3 + .2 * self.frame.time_numeric.to_numpy() + np.array([1, 1, 2, 2, 2, 3]) + np.arange(6)
        for values in simulated:
            np.testing.assert_allclose(values, expected)
        self.assertEqual(rng.normal.call_args_list[0].args, (0., 3., 3))
        self.assertEqual(rng.normal.call_args_list[1].args, (0., 2., 6))
        self.assertEqual(result['Hinge_Search_LR'], 6.)
        self.assertEqual(result['Hinge_Search_PValue'], 1 / 6)
        self.assertEqual(result['Bootstrap_Draws'], 5)
        pd.testing.assert_frame_equal(self.frame, original)

    def test_ties_add_one_and_failed_draws_never_create_small_pvalues(self):
        for scores, expected_p, failures in [([6., 1., 6. + 1e-12], .75, 0),
                                             ([1., ValueError('fit failed'), 1.], .5, 1)]:
            with self.subTest(scores=scores), patch.object(inference, '_refit_search', side_effect=scores) as refit:
                result = inference.bootstrap_hinge_test(
                    self.frame, self.likelihoods, lambda *a, **k: (null_result(), 'bfgs', ''),
                    protein='P', n_bootstrap=3, alpha=.99, show_progress=False)
            self.assertEqual(refit.call_count, 3)
            self.assertEqual(result['Hinge_Search_PValue'], expected_p)
            self.assertEqual(result['Bootstrap_Failures'], failures)
            self.assertEqual(result['Search_PValue_Is_Upper_Bound'], bool(failures))

    def test_early_stop_reports_upper_bound_not_a_stopped_point_estimate(self):
        with patch.object(inference, '_refit_search', return_value=7.) as refit:
            result = inference.bootstrap_hinge_test(
                self.frame, self.likelihoods, lambda *a, **k: (null_result(), 'bfgs', ''),
                protein='P', n_bootstrap=99, alpha=.05, show_progress=False)
        self.assertEqual(refit.call_count, 5)
        self.assertEqual(result['Bootstrap_Draws'], 5)
        self.assertEqual(result['Hinge_Search_PValue'], 1.)
        self.assertTrue(result['Search_PValue_Is_Upper_Bound'])
        self.assertIn('stopped', result['Hinge_Test_Status'])
        # Even if all 94 remaining draws were below the observed score, p would be .06.
        self.assertGreater((result['Bootstrap_Exceedances'] + 1) / 100, .05)

    def test_unavailable_observed_or_simulated_fits_cannot_claim_significance(self):
        with patch.object(inference, '_refit_search') as refit:
            fitter = Mock()
            result = inference.bootstrap_hinge_test(self.frame, {2.: -8., 3.: np.nan}, fitter,
                                                    protein='P', n_bootstrap=3, show_progress=False)
            fitter.assert_not_called()
            refit.assert_not_called()
            self.assertTrue(np.isnan(result['Hinge_Search_PValue']))
        for null in [None, null_result(cov_re=np.array([[-1.]])), null_result(llf=np.nan),
                     null_result(llf=0.)]:
            with self.subTest(null=null), patch.object(inference, '_refit_search') as refit:
                result = inference.bootstrap_hinge_test(
                    self.frame, self.likelihoods, lambda *a, **k: (null, None, 'failed'),
                    protein='P', n_bootstrap=3, show_progress=False)
                self.assertTrue(np.isnan(result['Hinge_Search_PValue']))
                refit.assert_not_called()
        with patch.object(inference, '_refit_search', side_effect=ValueError('all failed')):
            result = inference.bootstrap_hinge_test(
                self.frame, self.likelihoods, lambda *a, **k: (null_result(), 'bfgs', ''),
                protein='P', n_bootstrap=3, alpha=.99, show_progress=False)
        self.assertEqual(result['Bootstrap_Failures'], 3)
        self.assertTrue(np.isnan(result['Hinge_Search_PValue']))

    def test_zero_gain_returns_one_without_simulations(self):
        with patch.object(inference, '_refit_search') as refit:
            result = inference.bootstrap_hinge_test(
                self.frame, {2.: -10., 3.: -10.}, lambda *a, **k: (null_result(), 'bfgs', ''),
                protein='P', n_bootstrap=19, show_progress=False)
        refit.assert_not_called()
        self.assertEqual(result['Hinge_Search_PValue'], 1.)

    def test_seed_is_reproducible_and_protein_specific(self):
        draws = []
        for protein in ['P', 'Q', 'P']:
            series = []
            def record(frame, peaks, fit):
                series.append(frame.value.to_numpy().copy())
                return 0.
            with patch.object(inference, '_refit_search', side_effect=record):
                inference.bootstrap_hinge_test(
                    self.frame, self.likelihoods, lambda *a, **k: (null_result(), 'bfgs', ''),
                    protein=protein, n_bootstrap=3, seed=17, alpha=.99, show_progress=False)
            draws.append(np.array(series))
        np.testing.assert_array_equal(draws[0], draws[2])
        self.assertFalse(np.array_equal(draws[0], draws[1]))

    def test_bh_includes_unavailable_and_skipped_proteins_in_family(self):
        table = pd.DataFrame({'Model_Type': ['Mixed-Effects'] * 4 + ['Fixed-Effects'],
                              'Hinge_Search_PValue': [.001, .02, .04, np.nan, 1e-12]})
        adjusted = inference.adjust_hinge_fdr(table, family_size=6, alpha=.05)
        np.testing.assert_allclose(adjusted.Hinge_FDR, [.006, .06, .08, np.nan, np.nan], equal_nan=True)
        self.assertEqual(adjusted.Significant_FDR.tolist(), [True, False, False, False, False])

    def test_invalid_settings_are_rejected(self):
        for count, seed, alpha in [(0, 42, .05), (1.5, 42, .05), (True, 42, .05),
                                   (10, -1, .05), (10, 42, 0), (10, 42, np.nan)]:
            with self.subTest(settings=(count, seed, alpha)), self.assertRaises(ValueError):
                inference.validate_test_settings(count, seed, alpha)


class HingeInferenceIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.enterContext(redirect_stdout(io.StringIO()))
        self.folder = Path(self.enterContext(TemporaryDirectory()))

    def test_rankings_and_significance_use_corrected_search_not_nominal_pvalues(self):
        data = repeated_data().rename(columns={'mixed': 'nominal_winner'})
        data['search_winner'] = data.nominal_winner + 100
        data['absent'] = np.nan

        def fit(frame):
            model = mixed_result()
            model.pvalues.slope2_change = 1e-12 if frame.value.mean() < 100 else .04
            return model, 'bfgs', 'accepted'

        def bootstrap(frame, likelihoods, fitter, **kwargs):
            self.assertEqual(set(likelihoods), {2., 4.})
            self.assertEqual(kwargs['n_bootstrap'], 9)
            result = inference.unavailable_test(9, 42, 'Completed')
            result['Hinge_Search_PValue'] = .4 if frame.value.mean() < 100 else .001
            return result

        with patch.object(hinge, '_fit_mixed_model', side_effect=fit), \
             patch.object(hinge, 'bootstrap_hinge_test', side_effect=bootstrap), \
             patch.object(hinge, 'plot_best_fit') as plot:
            path = hinge.run_hinge_model(data, ['nominal_winner', 'search_winner', 'absent'],
                                        'Time', 'Sample_ID', ['D2', 'D4'], 1, self.folder, '-', n_bootstrap=9)
        sheets = pd.read_excel(path, sheet_name=None)
        results = sheets['All_Model_Results'].set_index('Protein')
        self.assertNotIn('Slope_2_Change_PValue', results)
        self.assertEqual(results.loc['search_winner', 'Hinge_FDR'], .003)
        self.assertEqual(list(sheets['Significant_Hinges'].Protein), ['search_winner'])
        self.assertEqual(list(sheets['Top_Mixed_Models'].Protein), ['search_winner'])
        self.assertEqual(plot.call_args.args[0], 'search_winner')

    def test_real_search_pipeline_including_one_candidate(self):
        data = repeated_data()
        for peaks in [['D4'], ['D2', 'D4']]:
            with self.subTest(peaks=peaks):
                path = hinge.run_hinge_model(data, ['mixed'], 'Time', 'Sample_ID', peaks,
                                            0, self.folder, '-', n_bootstrap=9, fdr_threshold=.9)
                result = pd.read_excel(path, sheet_name='All_Model_Results').iloc[0]
                self.assertTrue(np.isfinite(result.Hinge_Search_PValue), result.Hinge_Test_Status)
                self.assertGreaterEqual(result.Hinge_Search_PValue, .1)
                self.assertEqual(result.Hinge_Search_PValue, result.Hinge_FDR)
                self.assertGreater(result.Bootstrap_Draws, 0)
                self.assertGreater(result.Hinge_Search_LR, 0)

    def test_failed_observed_day_is_not_removed_to_rescue_significance(self):
        def fit(frame):
            if frame.time_numeric.iloc[0] - frame.slope1.iloc[0] == 2:
                return None, None, 'all optimisers failed'
            return mixed_result(), 'bfgs', 'accepted'
        with patch.object(hinge, '_fit_mixed_model', side_effect=fit), \
             patch.object(hinge, 'plot_best_fit') as plot:
            path = hinge.run_hinge_model(repeated_data(), ['mixed'], 'Time', 'Sample_ID',
                                        ['D2', 'D4'], 1, self.folder, '-', n_bootstrap=9)
        row = pd.read_excel(path, sheet_name='All_Model_Results').iloc[0]
        self.assertEqual(row.Model_Type, 'Mixed-Effects')
        self.assertTrue(np.isnan(row.Hinge_Search_PValue))
        self.assertTrue(np.isnan(row.Hinge_FDR))
        self.assertFalse(row.Significant_FDR)
        self.assertIn('eligible observed hinge fit failed', row.Hinge_Test_Status)
        plot.assert_not_called()

    def test_overlap_uses_fdr_and_refuses_legacy_nominal_only_workbooks(self):
        config = {'HINGE_GROUPS_TO_RUN': ['Resp', 'NonResp'], 'P_THRESHOLD': .05,
                  'OUTPUT_FOLDER': str(self.folder)}
        paths = {}
        for group in config['HINGE_GROUPS_TO_RUN']:
            table = pd.DataFrame({'Protein': ['nominal_only', 'corrected'],
                                  'Model_Type': ['Mixed-Effects'] * 2,
                                  'Slope_2_Change_PValue_Nominal': [1e-20, .03],
                                  'Hinge_Search_PValue': [.01, .02], 'Hinge_FDR': [.1, .04]})
            paths[group] = self.folder / f'{group}.xlsx'
            table.to_excel(paths[group], sheet_name='All_Model_Results', index=False)
        path = comparison.compare_hinge_results(paths, config)
        self.assertEqual(pd.read_excel(path, sheet_name='Common_Significant').Protein.tolist(), ['corrected'])
        table.drop(columns=['Hinge_Search_PValue', 'Hinge_FDR']).to_excel(
            paths['Resp'], sheet_name='All_Model_Results', index=False)
        with self.assertRaisesRegex(ValueError, 'rerun hinge analysis'):
            comparison.compare_hinge_results(paths, config)

    def test_main_passes_bootstrap_settings_and_fdr_threshold(self):
        data = repeated_data().assign(Condition='Resp')
        settings = dict(LAUNCH_GUI=False, RUN_HINGE=True, RUN_LIMMA=False, RUN_LOGISTIC_REGRESSION=False,
                        ANALYSIS_MODE='GROUP_COMPARISON', BASELINE_VAL='NonResp', COMPARE_VALS=['Resp'],
                        LOOP_VALS=[None], TIME_COLUMN='Time', CONDITION_COLUMN='Condition',
                        HINGE_GROUPS_TO_RUN=['Resp'], HINGE_PEAK_CANDIDATES=['D2', 'D4'],
                        HINGE_BOOTSTRAP_SIMULATIONS=23, HINGE_BOOTSTRAP_SEED=321, P_THRESHOLD=.1)
        with patch.dict(CONFIG, settings), patch('pandas.read_excel', return_value=data), \
             patch.dict(sys.modules, {'limma': SimpleNamespace(run_limma=Mock())}), \
             patch('os.makedirs'), patch.object(hinge, 'run_hinge_model', return_value=None) as fit, \
             patch.object(comparison, 'compare_hinge_results'):
            runpy.run_path(str(Path(__file__).resolve().parents[1] / 'main.py'))
        self.assertEqual(fit.call_args.kwargs['n_bootstrap'], 23)
        self.assertEqual(fit.call_args.kwargs['random_seed'], 321)
        self.assertEqual(fit.call_args.kwargs['fdr_threshold'], .1)
