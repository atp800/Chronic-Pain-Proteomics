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
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.preprocessing import StandardScaler

import logistic_regression as logistic
from config import CONFIG
from helpers import prepare_delta_comparison
from logistic_preprocessing import ProteinFeatures


class LogisticValidationTests(unittest.TestCase):
    def setUp(self):
        # Use one candidate and one process so every real CV fit can be inspected.
        self.enterContext(patch.object(logistic, "_PARAMETER_GRID", {
            "classifier__C": [1.0], "classifier__l1_ratio": [0.5],
        }))

        def serial_search(*args, **kwargs):
            kwargs["n_jobs"] = 1
            return GridSearchCV(*args, **kwargs)

        self.enterContext(patch.object(logistic, "GridSearchCV", serial_search))
        self.output = self.enterContext(redirect_stdout(io.StringIO()))
        self.y = np.repeat([0, 1], 5)
        self.X = pd.DataFrame({
            "signal": [np.nan, 8, 10, 11, 9, 13, 14, 16, 14, 15],
            # 3/5 observed in one group. Some inner splits have 3/4, others 2/4.
            "borderline_coverage": [1, 2, 3] + [np.nan] * 7,
            "one_outlier": [0] * 9 + [100],
        }, index=list("abcdefghij"))

    def test_each_inner_fit_selects_imputes_and_scales_its_own_training_rows(self):
        fits, scaler_inputs = [], []
        fit_features, fit_scaler = ProteinFeatures.fit, StandardScaler.fit

        def record_features(transformer, X, y):
            fit_features(transformer, X, y)
            fits.append((X.index.tolist(), y.copy(), transformer))
            return transformer

        def record_scaler(scaler, X, *args, **kwargs):
            scaler_inputs.append(X.copy())
            return fit_scaler(scaler, X, *args, **kwargs)

        with patch.object(ProteinFeatures, "fit", record_features), patch.object(StandardScaler, "fit", record_scaler):
            logistic._fit_tuned_model(self.X, self.y, delta=False)

        inner = StratifiedKFold(5).split(self.X, self.y)
        expected_rows = [train for train, _ in inner] + [np.arange(len(self.X))]
        self.assertEqual(len(fits), len(expected_rows))
        self.assertEqual(len(scaler_inputs), len(expected_rows))
        for rows, (ids, labels, fitted), scaled_input in zip(expected_rows, fits, scaler_inputs):
            train = self.X.iloc[rows]
            self.assertEqual(ids, train.index.tolist())
            np.testing.assert_array_equal(labels, self.y[rows])
            expected = ["signal"]
            if train.loc[labels == 0, "borderline_coverage"].notna().mean() >= 0.7:
                expected.append("borderline_coverage")
            if train["one_outlier"].max() > 0:
                expected.append("one_outlier")
            self.assertEqual(fitted.proteins_.tolist(), expected)
            observed = train[expected]
            np.testing.assert_allclose(fitted.fill_mean_, observed.mean() - 1.8 * observed.std())
            np.testing.assert_allclose(fitted.fill_sd_, 0.3 * observed.std())
            np.testing.assert_allclose(scaled_input, fitted.transform(train))
            self.assertTrue(np.isfinite(scaled_input).all())

    def test_outer_prediction_ignores_its_label_and_training_constant_protein(self):
        X = self.X.assign(held_out_only=0.0)
        before = logistic._leave_one_out_probabilities(X, self.y, delta=False)
        X.loc["a", "held_out_only"] = 1e9
        changed_labels = self.y.copy()
        changed_labels[0] = 1
        after = logistic._leave_one_out_probabilities(X, changed_labels, delta=False)
        np.testing.assert_array_equal(before[0], after[0])

    def test_delta_fits_on_raw_training_endpoints_and_preserves_subject_pairs(self):
        rows = []
        for i, subject in enumerate(self.X.index):
            for day, value in [("T0", 10 + i), ("D14", 10 + 2 * i), ("D3", 1e9)]:
                rows.append({
                    "sample": f"{subject}-{day}", "day": day, "group": self.y[i],
                    "signal": np.nan if subject == "i" and day == "T0" else value,
                    "few_pairs": np.nan if day == "D14" and subject not in "af" else value,
                })
        raw = pd.DataFrame(rows).sample(frac=1, random_state=12)
        deltas, endpoints = prepare_delta_comparison(
            raw, "sample", "day", "group", ["signal", "few_pairs"], "-", "T0", "D14"
        )
        X = pd.concat(endpoints, axis=1)
        train, test = X.loc[list("abcdefgh")], X.loc[["i", "j"]].copy()
        self.assertTrue(np.isnan(deltas.set_index("sample").at["i", "signal"]))
        fitted = ProteinFeatures(delta=True).fit(train, self.y[:8])
        self.assertEqual(fitted.proteins_.tolist(), ["signal"])
        measured = np.r_[10 + np.arange(8), 10 + 2 * np.arange(8)]
        np.testing.assert_allclose(fitted.fill_mean_, measured.mean() - 1.8 * measured.std(ddof=1))
        np.testing.assert_allclose(fitted.fill_sd_, 0.3 * measured.std(ddof=1))
        np.testing.assert_allclose(fitted.transform(train).ravel(), np.arange(8))
        imputed, complete = fitted.transform(test).ravel()
        self.assertTrue(np.isfinite(imputed))
        self.assertEqual(complete, 9)
        test.loc["j", ("baseline", "signal")] = -1e9
        self.assertEqual(fitted.transform(test.iloc[::-1])[1, 0], imputed)
        test.loc["i", ("followup", "signal")] += 2
        self.assertAlmostEqual(fitted.transform(test.loc[["i"]])[0, 0], imputed + 2)

    def test_main_routes_observed_deltas_to_limma_and_raw_endpoints_to_logistic(self):
        data = pd.DataFrame({
            "Sample_ID": ["b-D14", "a-D3", "c-T0", "a-T0", "b-T0", "a-D14", "c-D14"],
            "Time_Point": ["D14", "D3", "T0", "T0", "T0", "D14", "D14"],
            "Condition": ["Resp", "NonResp", "Other", "NonResp", "Resp", "NonResp", "Other"],
            "Age": [25, 30, 26, 30, 25, 30, 26],
            "p": [15., 1e6, -1e6, 9., np.nan, 13., 1e6],
        })
        settings = dict(LAUNCH_GUI=False, ANALYSIS_MODE="DELTA", DELTA_BASELINE="T0", LOOP_VALS=["D14"],
                        BASELINE_VAL="NonResp", COMPARE_VALS=["Resp"], COVARIATE_COLUMNS=["Age"],
                        RUN_LIMMA=True, RUN_LOGISTIC_REGRESSION=True, RUN_HINGE=False)
        limma = Mock()
        with (
            patch.dict(CONFIG, settings),
            patch.dict(sys.modules, {"limma": SimpleNamespace(run_limma=limma)}),
            patch("pandas.read_excel", return_value=data), patch("os.makedirs"),
            patch.object(logistic, "run_logistic_regression") as classifier,
        ):
            runpy.run_path(str(Path(__file__).resolve().parents[1] / "main.py"))
        observed = limma.call_args.args[0].set_index("Sample_ID")
        self.assertEqual(observed.at["a", "p"], 4)
        self.assertTrue(np.isnan(observed.at["b", "p"]))
        labels = classifier.call_args.args[0]
        raw = classifier.call_args.kwargs["delta_visits"]
        self.assertEqual(classifier.call_args.args[2], ["p"])
        self.assertEqual(labels.index.tolist(), ["a", "b"])
        self.assertTrue(raw.index.equals(labels.index))
        np.testing.assert_allclose(raw["baseline"]["p"], [9., np.nan])
        np.testing.assert_allclose(raw["followup"]["p"], [13., 15.])

    def test_fold_without_eligible_proteins_does_not_produce_a_performance_score(self):
        data = pd.DataFrame({"constant": [3.] * 10, "group": self.y})
        with TemporaryDirectory() as output:
            logistic.run_logistic_regression(data, "group", ["constant"], 0, 1, output)
            self.assertEqual(list(Path(output).iterdir()), [])
        self.assertIn("No proteins passed", self.output.getvalue())
        self.assertNotIn("Accuracy:", self.output.getvalue())
