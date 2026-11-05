import io
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pandas as pd
from sklearn.metrics import f1_score
from sklearn.model_selection import GridSearchCV

import logistic_regression as logistic


class LogisticDirectionTests(unittest.TestCase):
    def setUp(self):
        self.output = self.enterContext(redirect_stdout(io.StringIO()))
        self.enterContext(patch.object(logistic, '_PARAMETER_GRID', {
            'classifier__C': [1.0], 'classifier__l1_ratio': [.5],
        }))

        def serial_search(*args, **kwargs):
            kwargs['n_jobs'] = 1
            return GridSearchCV(*args, **kwargs)

        self.enterContext(patch.object(logistic, 'GridSearchCV', serial_search))
        self.X = pd.DataFrame({'marker': [1., 2, 3, 4, 9, 10, 11, 12]})
        self.encoded = np.repeat([0, 1], 4)

    def test_label_mapping_follows_requested_order_including_numeric_gui_labels(self):
        for baseline, target in [('T0', 'D14'), ('NonResp', 'Resp'), ('Z', 'A'), (14, 0)]:
            labels = pd.Series([target, baseline, target])
            encoded, classes = logistic._encode_comparison(labels, str(baseline), str(target))
            np.testing.assert_array_equal(encoded, [1, 0, 1])
            np.testing.assert_array_equal(classes, [str(baseline), str(target)])

    def test_invalid_labels_do_not_silently_become_baseline(self):
        for labels, baseline, target in [(['T0', None], 'T0', 'D14'),
                                         (['T0', 'D7'], 'T0', 'D14'),
                                         (['T0'], 'T0', 'T0'),
                                         (['T0'], None, 'T0')]:
            with self.subTest(labels=labels, baseline=baseline, target=target), self.assertRaises(ValueError):
                logistic._encode_comparison(pd.Series(labels), baseline, target)
        with patch.object(logistic, '_fit_tuned_model') as fit, \
             patch.object(logistic, '_leave_one_out_probabilities') as cv:
            logistic.run_logistic_regression(self.X.assign(group='unexpected'), 'group', ['marker'],
                                             'NonResp', 'Resp', 'unused')
        fit.assert_not_called()
        cv.assert_not_called()

    def test_reversing_comparison_reverses_real_coefficients_and_held_out_probabilities(self):
        labels = pd.Series(np.where(self.encoded, 'D14', 'T0'))
        ids = np.tile(np.arange(4), 2)
        coefficients, probabilities, aucs = [], [], []
        for baseline, target in [('T0', 'D14'), ('D14', 'T0')]:
            encoded, _ = logistic._encode_comparison(labels, baseline, target)
            model = logistic._fit_tuned_model(self.X, encoded, False, ids)
            probs = logistic._leave_one_out_probabilities(self.X, encoded, False, ids, show_progress=False)
            coefficients.append(model.named_steps['classifier'].coef_[0, 0])
            probabilities.append(probs)
            aucs.append(logistic._pooled_auc(encoded, probs))
        self.assertGreater(coefficients[0], 0)
        self.assertLess(coefficients[1], 0)
        np.testing.assert_allclose(coefficients[0], -coefficients[1], atol=.005)
        np.testing.assert_allclose(probabilities[0][:, 1], 1 - probabilities[1][:, 1], atol=.005)
        self.assertAlmostEqual(aucs[0], aucs[1])

    def test_all_four_analyses_use_comparison_class_for_fits_f1_roc_and_exports(self):
        # Unequal positive/negative F1 scores detect an accidentally reversed F1 target.
        positive_probs = np.array([.1, .7, .1, .1, .8, .3, .2, .9])
        probabilities = np.column_stack([1 - positive_probs, positive_probs])
        expected_f1 = f1_score(self.encoded, probabilities.argmax(axis=1))
        self.assertNotEqual(expected_f1, f1_score(self.encoded, probabilities.argmax(axis=1), pos_label=0))
        for mode in ['longitudinal', 'D14_group', 'baseline_group', 'delta_group']:
            baseline, target = ('T0', 'D14') if mode == 'longitudinal' else ('NonResp', 'Resp')
            data = self.X.assign(group=np.where(self.encoded, target, baseline))
            ids = pd.Series(np.tile(np.arange(4), 2)) if mode == 'longitudinal' else None
            endpoints = (pd.concat({'baseline': self.X * 0 + 10, 'followup': self.X + 10}, axis=1)
                         if mode == 'delta_group' else None)
            with self.subTest(mode=mode), TemporaryDirectory() as folder, \
                 patch.object(logistic, '_leave_one_out_probabilities', return_value=probabilities) as cv, \
                 patch.object(logistic, '_bootstrap_rows', return_value=[np.arange(8)] * 3) as bootstrap, \
                 patch.object(logistic, '_nested_permutation_test', return_value=(np.array([.2, .4, .6]), .25)) as permute, \
                 patch.object(logistic, 'roc_curve', wraps=logistic.roc_curve) as roc, \
                 patch.object(logistic, '_coefficient_stability_figure', wraps=logistic._coefficient_stability_figure) as plot:
                logistic.run_logistic_regression(data, 'group', ['marker'], baseline, target, folder,
                                                 participant_ids=ids, delta_visits=endpoints)
                np.testing.assert_array_equal(cv.call_args.args[1], self.encoded)
                np.testing.assert_array_equal(bootstrap.call_args.args[0], self.encoded)
                np.testing.assert_array_equal(permute.call_args.args[1], self.encoded)
                np.testing.assert_array_equal(roc.call_args.args[0], self.encoded)
                np.testing.assert_array_equal(roc.call_args.args[1], positive_probs)
                self.assertEqual(cv.call_args.args[2], mode == 'delta_group')
                if ids is not None:
                    np.testing.assert_array_equal(cv.call_args.args[3], ids)
                prefix = Path(folder) / f'LogReg_{target}_vs_{baseline}'
                features = pd.read_csv(str(prefix) + '_Features.csv')
                performance = pd.read_csv(str(prefix) + '_Performance.csv').iloc[0]
                self.assertTrue(features['Coefficient (Original)'].gt(0).all())
                self.assertTrue(features['Bootstrapped Mean'].gt(0).all())
                self.assertEqual(features['Positive Class (1)'].tolist(), [target])
                self.assertEqual(features['Reference Class (0)'].tolist(), [baseline])
                self.assertEqual(performance['Positive Class (1)'], target)
                self.assertEqual(performance['Reference Class (0)'], baseline)
                self.assertAlmostEqual(performance['F1 (Positive Class)'], expected_f1)
                self.assertEqual(plot.call_args.kwargs, {'reference_class': baseline, 'positive_class': target})
                self.assertTrue(Path(str(prefix) + '_Plots.pdf').is_file())
                self.assertIn(f'F1 Score (positive class: {target})', self.output.getvalue())

    def test_class_mapping_is_exported_even_when_no_protein_is_selected(self):
        data = self.X.assign(group=np.where(self.encoded, 'D14', 'T0'))
        model = logistic._fit_tuned_model(self.X, self.encoded, False)
        model.named_steps['classifier'].coef_[:] = 0
        probabilities = np.full((len(data), 2), .5)
        with TemporaryDirectory() as folder, \
             patch.object(logistic, '_fit_tuned_model', return_value=model), \
             patch.object(logistic, '_leave_one_out_probabilities', return_value=probabilities), \
             patch.object(logistic, '_bootstrap_rows', return_value=[np.arange(8)]), \
             patch.object(logistic, '_nested_permutation_test', side_effect=ValueError('unavailable')):
            logistic.run_logistic_regression(data, 'group', ['marker'], 'T0', 'D14', folder)
            features = pd.read_csv(Path(folder) / 'LogReg_D14_vs_T0_Features.csv')
            performance = pd.read_csv(Path(folder) / 'LogReg_D14_vs_T0_Performance.csv').iloc[0]
            self.assertTrue(features.empty)
            self.assertEqual(performance['Positive Class (1)'], 'D14')
            self.assertEqual(performance['Reference Class (0)'], 'T0')
            self.assertEqual(performance['F1 (Positive Class)'], 0)
