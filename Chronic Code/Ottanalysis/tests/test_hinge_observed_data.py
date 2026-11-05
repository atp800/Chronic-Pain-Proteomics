import io
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import hinge_model as hinge
from hinge_inference import unavailable_test
from test_hinge_fallback import repeated_data


class HingeObservedDataTests(unittest.TestCase):
    def setUp(self):
        test = unavailable_test(9, 42, 'Completed')
        test['Hinge_Search_PValue'] = .02
        self.enterContext(patch.object(hinge, 'bootstrap_hinge_test', return_value=test))
        self.enterContext(redirect_stdout(io.StringIO()))
        self.folder = Path(self.enterContext(TemporaryDirectory()))
        self.data = repeated_data()

    def run_hinge(self, proteins=('mixed',), peaks=('D2', 'D4'), plots=0):
        return hinge.run_hinge_model(self.data, list(proteins), 'Time', 'Sample_ID',
                                     peaks, plots, self.folder, '-')

    def observed_rows(self):
        frame = self.data.rename(columns={'mixed': 'value'}).copy()
        frame['value'] = frame['value'].replace([np.inf, -np.inf], np.nan)
        frame['subject_id_col'] = frame.Sample_ID.str.split('-').str[0]
        frame['time_numeric'] = frame.Time.str[1:].astype(int)
        return frame.dropna(subset=['value'])

    def test_all_attempts_and_fallbacks_use_only_the_same_observed_rows(self):
        self.data.loc[[0, 1, 3, 4], 'mixed'] = np.nan  # p0: one observed value.
        self.data.loc[5:9, 'mixed'] = np.nan           # p1: no measurements.
        self.data.loc[10:14, 'mixed'] = [9., np.nan, np.nan, np.nan, 9.]
        self.data.loc[16:18, 'mixed'] = [0., np.inf, -np.inf]  # A measured zero stays zero.
        self.data.loc[23, 'mixed'] = np.nan
        original = self.data.copy(deep=True)
        observed = self.observed_rows()
        calls = []
        real_ols = hinge.smf.ols

        def check_rows(frame):
            pd.testing.assert_series_equal(frame.value, observed.value)
            pd.testing.assert_series_equal(frame.subject_id_col, observed.subject_id_col)
            self.assertEqual(int(frame.subject_id_col.eq('p0').sum()), 1)
            self.assertFalse(frame.subject_id_col.eq('p1').any())
            self.assertEqual(frame.loc[16, 'value'], 0.)
            self.assertEqual(set(frame.loc[frame.subject_id_col.eq('p2'), 'value']), {9.})

        def mixed(formula, frame, groups):
            check_rows(frame)
            pd.testing.assert_series_equal(groups, observed.subject_id_col)
            calls.append(('mixed', frame.slope1.iloc[0]))
            return Mock(fit=Mock(side_effect=np.linalg.LinAlgError('force retries and fallback')))

        def ols(formula, data):
            check_rows(data)
            calls.append(('ols', data.slope1.iloc[0]))
            return real_ols(formula, data=data)

        with patch.object(hinge.smf, 'mixedlm', side_effect=mixed), \
             patch.object(hinge.smf, 'ols', side_effect=ols):
            path = self.run_hinge(peaks=('D2', 'D4', 'D8'))
        self.assertEqual([name for name, _ in calls], ['mixed'] * 3 + ['ols'] + ['mixed'] * 3 + ['ols'])
        pd.testing.assert_frame_equal(self.data, original)
        sheets = pd.read_excel(path, sheet_name=None)
        coverage = sheets['Data_Coverage'].iloc[0]
        self.assertEqual(coverage.N_Observed, 25)
        self.assertEqual(coverage.N_Missing, 15)
        self.assertEqual(coverage.N_Participants_Observed, 7)
        self.assertEqual(coverage.N_Participants_Repeated, 6)
        self.assertEqual(coverage.N_Participants_Single_Observation, 1)
        self.assertEqual(coverage.N_Participants_Unobserved, 1)
        self.assertEqual(coverage.N_Observed_Timepoints, 5)
        counts = sheets['Observed_By_Timepoint'].set_index('time_numeric')
        expected_counts = observed.groupby('time_numeric').size()
        np.testing.assert_array_equal(counts.N_Observed, expected_counts)
        self.assertTrue(counts.N_Recorded_Visits.eq(8).all())
        self.assertEqual(counts.N_Missing.sum(), 15)
        candidates = sheets['All_Candidate_Timepoints']
        self.assertTrue(candidates.N_Observed.eq(25).all())
        self.assertEqual(candidates.loc[candidates.Candidate_Time_Numeric.eq(8), 'Status'].iloc[0],
                         'No data on both sides')
        self.assertTrue(candidates.Slope_2_Change_PValue_Nominal.isna().all())
        self.assertEqual(sheets['All_Model_Results'].N_Observed.iloc[0], 25)

    def test_incomplete_mixed_fit_and_plot_match_direct_observed_data_fit(self):
        self.data.loc[self.data.Time.eq('D4'), 'mixed'] = np.nan  # No observation at the hinge itself.
        self.data.loc[1:4, 'mixed'] = np.nan                     # Keep one p0 observation.
        self.data.loc[5:9, 'mixed'] = np.nan                     # Remove all p1 observations.
        expected = self.observed_rows()
        expected['slope1'] = expected.time_numeric - 4
        expected['slope2_change'] = expected.slope1.clip(lower=0)
        reference, _, _ = hinge._fit_mixed_model(expected)
        self.assertIsNotNone(reference)
        original_plot = hinge.plot_best_fit
        captured = []

        def inspect_plot(*args, **kwargs):
            data = args[1]
            pd.testing.assert_series_equal(data.value, expected.value)
            means = kwargs['observed_means']
            self.assertTrue(np.isnan(means.loc[4]))
            np.testing.assert_allclose(means.dropna(), expected.groupby('time_numeric').value.mean())
            with patch.object(hinge.plt, 'close'):
                original_plot(*args, **kwargs)
                figure = plt.gcf()
            self.addCleanup(plt.close, figure)
            lines = figure.axes[0].lines
            fitted = next(line for line in lines if line.get_label() == 'Fitted Hinge Model (Average)')
            self.assertIn(4, fitted.get_xdata())
            # The singleton is visibly plotted as a point at its observed abundance.
            singleton = [line for line in lines if len(line.get_xdata()) == 1]
            self.assertTrue(any(np.allclose(line.get_ydata(), [expected.loc[0, 'value']])
                                and line.get_marker() == 'o' for line in singleton))
            self.assertTrue(any(line.get_label() == 'Observed Mean (Available Measurements)' for line in lines))
            captured.append(args[2])

        with patch.object(hinge, 'plot_best_fit', side_effect=inspect_plot):
            path = self.run_hinge(peaks=('D4',), plots=1)
        self.assertEqual(len(captured), 1)
        row = pd.read_excel(path, sheet_name='All_Model_Results').iloc[0]
        self.assertEqual(row.Model_Type, 'Mixed-Effects')
        self.assertEqual(row.N_Observed, len(expected))
        np.testing.assert_allclose(
            [row.Intercept, row.Slope_1_Gradient, row.Slope_2_Change, row.AIC, row.Slope_2_Change_PValue_Nominal],
            [reference.params.Intercept, reference.params.slope1, reference.params.slope2_change,
             reference.aic, reference.pvalues.slope2_change], rtol=1e-12,
        )

    def test_insufficient_observed_information_is_reported_without_fabricated_fits(self):
        self.data['absent'] = np.nan
        self.data['two_times'] = self.data.mixed.where(self.data.Time.isin(['D0', 'D8']))
        self.data['saturated'] = self.data.mixed.where(self.data.index.isin([0, 2, 4]))
        self.data['flat'] = 9.
        self.data.loc[self.data.Time.eq('D2'), 'flat'] = np.nan
        proteins = ['absent', 'two_times', 'saturated', 'flat']
        with patch.object(hinge, '_fit_mixed_model') as mixed, patch.object(hinge.smf, 'ols') as ols:
            self.assertIsNone(self.run_hinge(proteins))
        mixed.assert_not_called()
        ols.assert_not_called()
        sheets = pd.read_excel(self.folder / 'hinge_candidate_diagnostics.xlsx', sheet_name=None)
        coverage = sheets['Data_Coverage'].set_index('Protein')
        self.assertEqual(coverage.at['absent', 'N_Observed'], 0)
        self.assertEqual(coverage.at['absent', 'N_Participants_Unobserved'], 8)
        self.assertEqual(coverage.at['two_times', 'Status'], 'Fewer than three observed time points')
        self.assertEqual(coverage.at['saturated', 'Status'], 'Insufficient observed residual degrees of freedom')
        self.assertEqual(coverage.at['flat', 'Status'], 'No observed abundance variation')
        candidates = sheets['All_Candidate_Timepoints']
        self.assertEqual(len(candidates), 8)
        self.assertTrue(candidates.Slope_2_Change_PValue_Nominal.isna().all())
        self.assertFalse(candidates.Is_Best_Hinge.any())

    def test_insufficient_observed_repeats_allow_only_descriptive_fallback(self):
        complete = self.data.mixed.copy()
        for add_one_repeat in [False, True]:
            with self.subTest(one_repeated_participant=add_one_repeat):
                self.data['mixed'] = np.nan
                rows = [person * 5 + person % 5 for person in range(8)]
                if add_one_repeat:
                    rows.append(1)
                self.data.loc[rows, 'mixed'] = complete.loc[rows]
                with patch.object(hinge, '_fit_mixed_model') as mixed, \
                     patch.object(hinge, 'plot_best_fit') as plot:
                    path = self.run_hinge(plots=20)
                mixed.assert_not_called()
                plot.assert_not_called()
                row = pd.read_excel(path, sheet_name='All_Model_Results').iloc[0]
                self.assertEqual(row.Model_Type, 'Fixed-Effects')
                self.assertEqual(row.N_Participants_Repeated, int(add_one_repeat))
                self.assertEqual(row.N_Observed, 8 + int(add_one_repeat))
                self.assertTrue(np.isnan(row.Slope_2_Change_PValue_Nominal))
                self.assertIn('fewer than two participants', row.Mixed_Fit_Diagnostics)
