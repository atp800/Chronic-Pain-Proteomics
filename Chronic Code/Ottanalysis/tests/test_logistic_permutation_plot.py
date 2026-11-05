import io
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

import logistic_regression as logistic
from logistic_preprocessing import ProteinFeatures


class LogisticPermutationPlotTests(unittest.TestCase):
    def test_degenerate_and_ordinary_scores_finish_the_report_without_changing_results(self):
        X = pd.DataFrame({'marker': [1., 2, 3, 4, 9, 10, 11, 12]})
        y = np.repeat([0, 1], 4)
        data = X.assign(group=np.where(y, 'D14', 'T0'))
        model = Pipeline([
            ('proteins', ProteinFeatures()), ('scale', StandardScaler()),
            ('classifier', LogisticRegression(penalty='elasticnet', solver='saga',
                                               C=1, l1_ratio=.5, random_state=42)),
        ]).fit(X, y)
        # Only reporting is under test; statistical fitting/validation has separate tests.
        probabilities = np.column_stack([1 - y * .8 - .1, y * .8 + .1])
        cases = {
            'identical': np.full(100, .5),
            'one_ulp': np.tile([.5, np.nextafter(.5, 1.)], 50),
            'tiny_range': np.linspace(.5, .5 + 1e-15, 100),
            'zero_boundary': np.tile([0., np.nextafter(0., 1.)], 50),
            'one_boundary': np.tile([np.nextafter(1., 0.), 1.], 50),
            'ordinary': np.linspace(0, 1, 100),
        }
        real_hist = plt.hist
        for name, scores in cases.items():
            original = scores.copy()
            histograms = []

            def record_hist(values, *args, **kwargs):
                counts, edges, patches = real_hist(values, *args, **kwargs)
                histograms.append((np.array(values, copy=True), counts, edges))
                return counts, edges, patches

            with self.subTest(name=name), TemporaryDirectory() as folder, \
                 redirect_stdout(io.StringIO()) as output, \
                 patch.object(logistic, '_fit_tuned_model', return_value=model), \
                 patch.object(logistic, '_leave_one_out_probabilities', return_value=probabilities), \
                 patch.object(logistic, '_bootstrap_rows', return_value=[np.arange(8)]), \
                 patch.object(logistic, '_nested_permutation_test', return_value=(scores, .1234)), \
                 patch.object(plt, 'hist', side_effect=record_hist), \
                 patch.object(plt, 'axvline', wraps=plt.axvline) as observed_line, \
                 patch.object(logistic, '_coefficient_stability_figure', wraps=logistic._coefficient_stability_figure) as stability:
                logistic.run_logistic_regression(data, 'group', ['marker'], 'T0', 'D14', folder)
                np.testing.assert_array_equal(scores, original)
                self.assertEqual(len(histograms), 1)
                plotted, density, edges = histograms[0]
                np.testing.assert_array_equal(plotted, original)
                self.assertTrue(np.isfinite(density).all())
                self.assertTrue((np.diff(edges) > 0).all())
                self.assertAlmostEqual(np.sum(density * np.diff(edges)), 1.)
                # No draws, including boundary AUCs, are dropped from the histogram.
                counts, _ = np.histogram(original, bins=edges)
                self.assertEqual(counts.sum(), len(original))
                np.testing.assert_allclose(density * np.diff(edges) * len(original), counts)
                observed_line.assert_called_once()
                self.assertEqual(observed_line.call_args.args[0], 1.)
                self.assertIn('p-value = 0.1234', observed_line.call_args.kwargs['label'])
                self.assertIn('permutation p-value (100 draws): 0.1234', output.getvalue())
                stability.assert_called_once()  # The report continues after the histogram.
                report = Path(folder) / 'LogReg_D14_vs_T0_Plots.pdf'
                content = report.read_bytes()
                self.assertTrue(content.startswith(b'%PDF'))
                self.assertTrue(content.rstrip().endswith(b'%%EOF'))
                self.assertTrue((Path(folder) / 'LogReg_D14_vs_T0_Features.csv').is_file())
                self.assertTrue((Path(folder) / 'LogReg_D14_vs_T0_Performance.csv').is_file())
