import io
import unittest
import warnings
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import Mock, patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.tools.sm_exceptions import ConvergenceWarning

import hinge_model as hinge
from compare_hinge_results import compare_hinge_results


def mixed_result(**changes):
    result = SimpleNamespace(
        converged=True, aic=10000.,
        params=pd.Series({'Intercept': 12., 'slope1': 2., 'slope2_change': -3., 'Group Var': 1.}),
        bse_fe=pd.Series([.2, .1, .15]),
        pvalues=pd.Series({'Intercept': .001, 'slope1': .01, 'slope2_change': .02}),
    )
    for name, value in changes.items():
        setattr(result, name, value)
    return result


def repeated_data():
    rng = np.random.default_rng(172)
    rows = []
    for person in range(8):
        offset = rng.normal(scale=2.)
        for time in [0, 2, 4, 6, 8]:
            value = 12 + offset + .7 * time - 1.3 * max(time - 4, 0) + rng.normal(scale=.3)
            rows.append([f'p{person}-D{time}', f'D{time}', value])
    return pd.DataFrame(rows, columns=['Sample_ID', 'Time', 'mixed'])


class MixedRetryTests(unittest.TestCase):
    def setUp(self):
        self.data = repeated_data().rename(columns={'mixed': 'value'})
        self.data['subject_id_col'] = self.data.Sample_ID.str.split('-').str[0]
        self.data['slope1'] = np.tile([0, 2, 4, 6, 8], 8) - 4
        self.data['slope2_change'] = self.data.slope1.clip(lower=0)

    def test_first_usable_fit_stops_retries(self):
        result = mixed_result()
        with patch.object(hinge.smf, 'mixedlm') as create:
            create.return_value.fit.return_value = result
            fitted, optimizer, log = hinge._fit_mixed_model(self.data)
        self.assertIs(fitted, result)
        self.assertEqual(optimizer, 'lbfgs')
        self.assertEqual(create.call_count, 1)
        self.assertIn('lbfgs: accepted', log)

    def test_failure_modes_retry_same_model_and_record_warnings(self):
        invalid_p = mixed_result()
        invalid_p.pvalues['slope2_change'] = np.nan
        invalid_params = mixed_result()
        invalid_params.params['slope1'] = np.inf
        cases = [
            np.linalg.LinAlgError('singular matrix'),
            mixed_result(converged=False), mixed_result(aic=np.nan),
            mixed_result(bse_fe=np.array([.1, np.nan, .2])),
            mixed_result(bse_fe=np.array([.1, 0., .2])), invalid_p, invalid_params,
        ]
        for first in cases:
            with self.subTest(first=first):
                models = []
                usable = mixed_result()

                def create(formula, data, groups):
                    self.assertEqual(formula, 'value ~ slope1 + slope2_change')
                    pd.testing.assert_frame_equal(data, self.data)
                    pd.testing.assert_series_equal(groups, self.data.subject_id_col)
                    model = Mock()
                    def fit(**kwargs):
                        if len(models) == 1:
                            warnings.warn('first attempt warning', ConvergenceWarning)
                            if isinstance(first, Exception):
                                raise first
                            return first
                        return usable
                    model.fit.side_effect = fit
                    models.append(model)
                    return model

                with patch.object(hinge.smf, 'mixedlm', side_effect=create):
                    result, optimizer, log = hinge._fit_mixed_model(self.data)
                self.assertIs(result, usable)
                self.assertEqual(optimizer, 'bfgs')
                self.assertEqual(len(models), 2)
                for model, method in zip(models, ['lbfgs', 'bfgs']):
                    model.fit.assert_called_once_with(method=method, reml=False, maxiter=200)
                self.assertIn('first attempt warning', log)
                self.assertIn('bfgs: accepted', log)

    def test_third_optimizer_is_tried_before_giving_up(self):
        for succeeds in [True, False]:
            with self.subTest(succeeds=succeeds), patch.object(hinge.smf, 'mixedlm') as create:
                usable = mixed_result()
                create.return_value.fit.side_effect = [
                    ValueError('first failed'), mixed_result(converged=False),
                    usable if succeeds else ValueError('third failed'),
                ]
                result, optimizer, log = hinge._fit_mixed_model(self.data)
                self.assertEqual([c.kwargs['method'] for c in create.return_value.fit.call_args_list],
                                 ['lbfgs', 'bfgs', 'cg'])
                self.assertIs(result, usable if succeeds else None)
                self.assertEqual(optimizer, 'cg' if succeeds else None)
                self.assertIn('cg:', log)

    def test_real_mixed_fit_can_be_recovered_by_retry(self):
        create = hinge.smf.mixedlm
        reference = create('value ~ slope1 + slope2_change', self.data,
                           groups=self.data.subject_id_col).fit(method='bfgs', reml=False, maxiter=200)
        calls = []

        def first_fails(*args, **kwargs):
            calls.append(kwargs)
            if len(calls) == 1:
                raise np.linalg.LinAlgError('forced first-attempt failure')
            return create(*args, **kwargs)

        with patch.object(hinge.smf, 'mixedlm', side_effect=first_fails):
            fitted, optimizer, log = hinge._fit_mixed_model(self.data)
        self.assertEqual(optimizer, 'bfgs')
        self.assertTrue(fitted.converged)
        pd.testing.assert_series_equal(fitted.params, reference.params)
        pd.testing.assert_series_equal(fitted.pvalues, reference.pvalues)
        self.assertIn('forced first-attempt failure', log)


class HingeReportingTests(unittest.TestCase):
    def setUp(self):
        self.output = self.enterContext(redirect_stdout(io.StringIO()))
        self.folder = Path(self.enterContext(TemporaryDirectory()))
        self.data = repeated_data()
        self.data['fallback'] = self.data['mixed'] + 100

    def run_hinge(self, proteins, peaks=('D2', 'D4'), num_plots=20):
        return hinge.run_hinge_model(self.data, proteins, 'Time', 'Sample_ID', peaks,
                                     num_plots, self.folder, '-')

    def test_fallback_cannot_win_selection_or_enter_top_plots(self):
        chosen = mixed_result()

        def fit(data):
            # OLS at D2 has much lower AIC, but the usable mixed fit at D4 must win.
            if data.value.mean() < 100 and data.slope1.iloc[0] == -4:
                return chosen, 'bfgs', 'lbfgs: failed; bfgs: accepted'
            return None, None, 'all mixed attempts failed'

        plots = self.folder / 'Hinge_Model_Plots'
        plots.mkdir()
        (plots / 'fallback_hinge_plot.png').write_bytes(b'previous plot')
        (plots / 'notes.txt').write_text('keep')
        with patch.object(hinge, '_fit_mixed_model', side_effect=fit) as fitter, \
             patch.object(hinge, 'plot_best_fit') as plot:
            result = self.run_hinge(['mixed', 'fallback'])
        self.assertEqual(fitter.call_count, 4)  # No additional fit during plotting.
        sheets = pd.read_excel(result, sheet_name=None)
        selected = sheets['All_Model_Results'].set_index('Protein')
        self.assertEqual(selected.at['mixed', 'Best_Peak_Numeric'], 4)
        self.assertEqual(selected.at['mixed', 'Mixed_Optimizer'], 'bfgs')
        self.assertEqual(selected.at['fallback', 'Model_Type'], 'Fixed-Effects')
        self.assertIn('Descriptive only', selected.at['fallback', 'Inference_Status'])
        self.assertTrue(selected.loc['fallback', ['Slope_1_PValue', 'Slope_2_Change_PValue']].isna().all())
        self.assertTrue(np.isfinite(selected.loc['fallback', ['Intercept', 'Slope_1_Gradient', 'Slope_2_Gradient']].astype(float)).all())
        candidates = sheets['All_Candidate_Timepoints']
        fixed = candidates[candidates.Model_Type.eq('Fixed-Effects')]
        self.assertTrue(fixed.Slope_2_Change_PValue.isna().all())
        self.assertLess(fixed.AIC.max(), selected.at['mixed', 'AIC'])
        self.assertEqual(list(sheets['Top_Mixed_Models'].Protein), ['mixed'])
        self.assertTrue(sheets['Fixed_Effects_Models'].Slope_2_Change_PValue.isna().all())
        self.assertEqual(plot.call_count, 1)
        self.assertEqual(plot.call_args.args[0], 'mixed')
        plotted_row = plot.call_args.args[2]
        self.assertEqual(plotted_row['Intercept'], chosen.params['Intercept'])
        self.assertEqual(plotted_row['Slope_2_Change'], chosen.params['slope2_change'])
        self.assertFalse((plots / 'fallback_hinge_plot.png').exists())
        archived = list(plots.glob('previous_*/fallback_hinge_plot.png'))
        self.assertEqual(len(archived), 1)
        self.assertEqual(archived[0].read_bytes(), b'previous plot')
        self.assertEqual((plots / 'notes.txt').read_text(), 'keep')

    def test_all_mixed_attempts_fail_but_descriptive_results_are_retained(self):
        with patch.object(hinge.smf, 'mixedlm', side_effect=np.linalg.LinAlgError('singular')) as fit, \
             patch.object(hinge, 'plot_best_fit') as plot:
            result = self.run_hinge(['fallback'], peaks=('D0', 'D2', 'D4', 'D8'))
        self.assertEqual(fit.call_count, 6)  # Three attempts at each eligible hinge.
        plot.assert_not_called()
        sheets = pd.read_excel(result, sheet_name=None)
        self.assertEqual(len(sheets['All_Model_Results']), 1)
        self.assertEqual(len(sheets['Fixed_Effects_Models']), 1)
        self.assertTrue(sheets['Mixed_Effects_Models'].empty)
        self.assertTrue(sheets['Top_Mixed_Models'].empty)
        candidates = sheets['All_Candidate_Timepoints']
        self.assertEqual(candidates.Is_Best_Hinge.sum(), 1)
        self.assertTrue(candidates.Slope_2_Change_PValue.isna().all())
        eligible = candidates[candidates.Status.eq('Fitted')]
        self.assertTrue(eligible.Mixed_Fit_Diagnostics.str.contains('cg: LinAlgError: singular').all())
        self.assertEqual(sheets['All_Model_Results'].AIC.iloc[0], eligible.AIC.min())

    def test_total_failure_still_exports_candidate_diagnostics(self):
        with patch.object(hinge.smf, 'mixedlm', side_effect=ValueError('mixed failed')), \
             patch.object(hinge.smf, 'ols', side_effect=ValueError('OLS failed')):
            self.assertIsNone(self.run_hinge(['mixed']))
        table = pd.read_excel(self.folder / 'hinge_candidate_diagnostics.xlsx')
        self.assertTrue(table.Slope_2_Change_PValue.isna().all())
        self.assertFalse(table.Is_Best_Hinge.any())
        self.assertTrue(table.Status.str.contains('OLS also failed').all())
        self.assertTrue(table.Mixed_Fit_Diagnostics.str.contains('cg:').all())

    def test_real_run_exports_and_plots_the_same_mixed_fit(self):
        original_plot = hinge.plot_best_fit
        captured = []

        def inspect_plot(*args, **kwargs):
            with patch.object(hinge.plt, 'close'):
                original_plot(*args, **kwargs)
                figure = plt.gcf()
            row = args[2]
            line = next(line for line in figure.axes[0].lines
                        if line.get_label() == 'Fitted Hinge Model (Average)')
            times = line.get_xdata()
            expected = (row.Intercept + row.Slope_1_Gradient * (times - row.Best_Peak_Numeric)
                        + row.Slope_2_Change * np.maximum(times - row.Best_Peak_Numeric, 0))
            np.testing.assert_allclose(line.get_ydata(), expected)
            captured.append(row)
            plt.close(figure)

        with patch.object(hinge, 'plot_best_fit', side_effect=inspect_plot):
            result = self.run_hinge(['mixed'])
        sheets = pd.read_excel(result, sheet_name=None)
        row = sheets['All_Model_Results'].iloc[0]
        self.assertEqual(row.Model_Type, 'Mixed-Effects')
        self.assertEqual(len(captured), 1)
        candidates = sheets['All_Candidate_Timepoints']
        mixed = candidates[candidates.Model_Type.eq('Mixed-Effects')]
        self.assertAlmostEqual(row.AIC, mixed.AIC.min())
        self.assertTrue((self.folder / 'Hinge_Model_Plots' / 'mixed_hinge_plot.png').is_file())

    def test_group_overlap_excludes_legacy_fallback_pvalues(self):
        config = {'HINGE_GROUPS_TO_RUN': ['Resp', 'NonResp'], 'P_THRESHOLD': .05,
                  'OUTPUT_FOLDER': str(self.folder)}
        for all_fixed in [False, True]:
            with self.subTest(all_fixed=all_fixed):
                paths = {}
                for group in config['HINGE_GROUPS_TO_RUN']:
                    table = pd.DataFrame({
                        'Protein': ['mixed', 'fallback'],
                        'Model_Type': ['Fixed-Effects' if all_fixed else 'Mixed-Effects', 'Fixed-Effects'],
                        'Slope_1_PValue': [.01, 1e-12],
                        'Slope_2_Change_PValue': [.02, 1e-12],
                    })
                    paths[group] = self.folder / f'{group}.xlsx'
                    table.to_excel(paths[group], sheet_name='All_Model_Results', index=False)
                result = compare_hinge_results(paths, config)
                sheets = pd.read_excel(result, sheet_name=None)
                self.assertEqual(list(sheets['Common_Significant'].Protein), [] if all_fixed else ['mixed'])
                self.assertTrue(sheets['Groups'].Proteins_Tested.eq(0 if all_fixed else 1).all())
                self.assertTrue(sheets['Groups'].Descriptive_Fallbacks.eq(2 if all_fixed else 1).all())
                reported = sheets['Group_1_Results']
                self.assertTrue(reported.loc[reported.Model_Type.eq('Fixed-Effects'), 'Slope_2_Change_PValue'].isna().all())
