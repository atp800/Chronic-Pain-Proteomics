import unittest

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from logistic_regression import _coefficient_stability_figure, _coefficient_stability_table


class LogisticStabilityTests(unittest.TestCase):
    def setUp(self):
        names = ["rare", "absent_from_full_fit", "negative", "positive", "mixed"]
        original = np.array([0.1, 0., -3., 4., 2.])
        draws = np.column_stack([
            np.r_[np.zeros(99), 100.], np.full(100, 4.), np.full(100, -2.),
            np.linspace(1., 3., 100), np.tile([-1., 0., 1., 0.], 25),
        ])
        self.table = _coefficient_stability_table(names, original, draws).reset_index(drop=True)

    def test_table_reports_selected_proteins_and_includes_zero_draws_in_every_summary(self):
        expected = pd.DataFrame({
            "Protein": ["positive", "negative", "mixed", "rare"],
            "Coefficient (Original)": [4., -3., 2., 0.1],
            "Bootstrapped Mean": [2., -2., 0., 1.],
            "Bootstrap 2.5th Percentile": [1.05, -2., -1., 0.],
            "Bootstrap 97.5th Percentile": [2.95, -2., 1., 0.],
            "Selection Frequency (%)": [100., 100., 50., 1.],
        })
        pd.testing.assert_frame_equal(self.table, expected)
        empty = _coefficient_stability_table(["protein"], np.zeros(1), np.ones((100, 1)))
        self.assertTrue(empty.empty)
        self.assertEqual(list(empty.columns), list(expected.columns))

    def test_plot_draws_actual_percentiles_even_when_mean_lies_outside_them(self):
        fig = _coefficient_stability_figure(self.table)
        self.addCleanup(plt.close, fig)
        ax = fig.axes[0]
        ranges, means = ax.collections
        endpoints = self.table[["Bootstrap 2.5th Percentile", "Bootstrap 97.5th Percentile"]]
        np.testing.assert_allclose(np.asarray(ranges.get_segments())[:, :, 0], endpoints)
        np.testing.assert_allclose(means.get_offsets()[:, 0], self.table["Bootstrapped Mean"])
        self.assertIn("percentiles", ranges.get_label())
        self.assertIn("Preprocessing and tuning fixed", ax.get_title())
