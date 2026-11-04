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
from sklearn.model_selection import GridSearchCV, StratifiedKFold, permutation_test_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils import resample

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
        # Each participant has one row of each class, in different orders at the two visits.
        self.participants = pd.Series(["p3", "p1", "p5", "p2", "p4", "p4", "p5", "p1", "p3", "p2"],
                                      index=self.X.index)

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

    def test_longitudinal_outer_predictions_and_inner_fits_hold_out_whole_participants(self):
        training_sets, predicted = [], []
        expected_probabilities = pd.DataFrame(index=self.X.index, columns=[0, 1], dtype=float)
        fit_features, predict = ProteinFeatures.fit, Pipeline.predict_proba

        def inspect_fit(transformer, X, y):
            ids = self.participants.loc[X.index]
            self.assertTrue(ids.value_counts().eq(2).all())
            training_sets.append(set(ids))
            return fit_features(transformer, X, y)

        def inspect_prediction(model, X):
            held_out = self.participants.loc[X.index]
            self.assertEqual(len(held_out), 2)
            self.assertEqual(held_out.nunique(), 1)
            remaining = set(self.participants) - set(held_out)
            self.assertTrue(all(inner <= remaining for inner in training_sets))
            self.assertEqual(training_sets[-1], remaining)  # refit on the outer training set
            self.assertCountEqual([remaining - inner for inner in training_sets[:-1]],
                                  [{participant} for participant in remaining])
            training_sets.clear()
            predicted.extend(X.index)
            result = predict(model, X)
            expected_probabilities.loc[X.index] = result
            return result

        with patch.object(ProteinFeatures, "fit", inspect_fit), patch.object(Pipeline, "predict_proba", inspect_prediction):
            probabilities = logistic._leave_one_out_probabilities(self.X, self.y, False, self.participants.to_numpy())
        self.assertCountEqual(predicted, self.X.index)
        np.testing.assert_allclose(probabilities, expected_probabilities)

    def test_bootstrap_draws_whole_pairs_and_preserves_independent_resampling(self):
        ids = self.participants.to_numpy()
        repeated = False
        for rows in logistic._bootstrap_rows(self.y, ids, 12):
            self.assertEqual(len(rows), len(ids))
            counts = np.bincount(rows, minlength=len(ids))
            for participant in np.unique(ids):
                first, second = counts[ids == participant]
                self.assertEqual(first, second)
            repeated |= counts.max() > 1
        self.assertTrue(repeated)

        for seed, rows in enumerate(logistic._bootstrap_rows(self.y, None, 3)):
            expected = resample(np.arange(len(self.y)), stratify=self.y, random_state=seed)
            np.testing.assert_array_equal(rows, expected)

    def test_paired_permutations_swap_only_within_people_and_refit_on_whole_training_pairs(self):
        swaps = []
        observed = pd.Series(self.y, index=self.X.index)
        # With only "signal" in this comparison, p3 has no measurement at its first visit.
        eligible_ids = self.participants[self.participants.ne("p3")]
        fit_features = ProteinFeatures.fit

        def inspect_fit(transformer, X, y):
            labels = pd.Series(y, index=X.index)
            ids = self.participants.loc[X.index]
            for _, pair in labels.groupby(ids):
                self.assertCountEqual(pair, [0, 1])
            changed = labels.ne(observed.loc[X.index])
            self.assertTrue(changed.groupby(ids).nunique().eq(1).all())
            swaps.append(changed.any() and not changed.all())
            return fit_features(transformer, X, y)

        def run_permutations(model, X, y, **options):
            ids = options["groups"]
            self.assertTrue(X.index.equals(eligible_ids.index))
            np.testing.assert_array_equal(ids, eligible_ids.to_numpy())
            for train, test in options["cv"].split(X, y, ids):
                self.assertFalse(set(ids[train]) & set(ids[test]))
            with patch.object(ProteinFeatures, "fit", inspect_fit):
                return permutation_test_score(model, X, y, **dict(options, n_jobs=1, n_permutations=5))

        data = self.X.assign(visit=self.y)
        with (TemporaryDirectory() as output, patch.object(logistic, "permutation_test_score", run_permutations),
              patch.object(logistic, "_bootstrap_rows", wraps=logistic._bootstrap_rows) as bootstrap):
            logistic.run_logistic_regression(data, "visit", ["signal"], 0, 1, output, participant_ids=self.participants)
            np.testing.assert_array_equal(bootstrap.call_args.args[1], eligible_ids.to_numpy())
            self.assertTrue((Path(output) / "LogReg_1_vs_0_Plots.pdf").is_file())
        self.assertTrue(any(swaps))

    def test_incomplete_or_ambiguous_pairs_do_not_fall_back_to_row_wise_validation(self):
        data = pd.DataFrame({"visit": ["T0", "D14"] * 4, "p": [7., 9, 6, 7, 3, np.nan, 5, 2]})
        ids = pd.Series(["a", "a", "b", "b", "c", "c", "d", None])
        keep = logistic._complete_pair_mask(data[["p"]], data.visit, ids)
        self.assertEqual(ids[keep].tolist(), ["a", "a", "b", "b"])
        with TemporaryDirectory() as output:
            logistic.run_logistic_regression(data, "visit", ["p"], "T0", "D14", output, participant_ids=ids)
            self.assertIn("at least three complete pairs", self.output.getvalue())
            ids.iloc[6] = "a"  # duplicate a's T0
            logistic.run_logistic_regression(data, "visit", ["p"], "T0", "D14", output, participant_ids=ids)
            self.assertIn("Duplicate participant/visit", self.output.getvalue())
            unmatched = pd.Series(list("abcdefgh"), index=data.index)
            logistic.run_logistic_regression(data, "visit", ["p"], "T0", "D14", output, participant_ids=unmatched)
            self.assertIn("0 complete participant pairs", self.output.getvalue())
            self.assertEqual(list(Path(output).iterdir()), [])

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
        self.assertIsNone(classifier.call_args.kwargs["participant_ids"])
        self.assertEqual(classifier.call_args.args[2], ["p"])
        self.assertEqual(labels.index.tolist(), ["a", "b"])
        self.assertTrue(raw.index.equals(labels.index))
        np.testing.assert_allclose(raw["baseline"]["p"], [9., np.nan])
        np.testing.assert_allclose(raw["followup"]["p"], [13., 15.])

    def test_main_only_passes_participant_ids_for_longitudinal_comparisons(self):
        data = pd.DataFrame({
            "Sample_ID": ["b|D14", "a|D3", "c|T0", "a|T0", "b|T0", "a|D14", "c|D14"],
            "Time_Point": ["D14", "D3", "T0", "T0", "T0", "D14", "D14"],
            "Condition": ["Resp", "Resp", "NonResp", "Resp", "Resp", "Resp", "NonResp"],
            "p": np.arange(7.),
        }, index=[18, 2, 35, 7, 41, 5, 11])
        common = dict(LAUNCH_GUI=False, ID_DELIMITER="|", ID_COLUMN="Sample_ID", CONDITION_COLUMN="Condition",
                      TIME_COLUMN="Time_Point", COVARIATE_COLUMNS=[], UNNEEDED_COLUMNS=[],
                      RUN_LIMMA=False, RUN_LOGISTIC_REGRESSION=True, RUN_HINGE=False)
        cases = [("LONGITUDINAL", ["Resp"], "T0", ["D14"]), ("GROUP_COMPARISON", ["T0", "D14"], "NonResp", ["Resp"])]
        for mode, loops, base, targets in cases:
            with (self.subTest(mode=mode), patch.dict(CONFIG, common, ANALYSIS_MODE=mode, BASELINE_VAL=base,
                                                       COMPARE_VALS=targets, LOOP_VALS=loops),
                  patch.dict(sys.modules, {"limma": SimpleNamespace(run_limma=Mock())}),
                  patch("pandas.read_excel", return_value=data), patch("os.makedirs"),
                  patch.object(logistic, "run_logistic_regression") as classifier):
                runpy.run_path(str(Path(__file__).resolve().parents[1] / "main.py"))
                self.assertEqual(classifier.call_count, len(loops))
                for call in classifier.call_args_list:
                    ids = call.kwargs["participant_ids"]
                    self.assertEqual(call.args[2], ["p"])
                    if mode == "LONGITUDINAL":
                        self.assertTrue(ids.index.equals(call.args[0].index))
                        self.assertEqual(ids.to_dict(), {18: "b", 7: "a", 41: "b", 5: "a"})
                    else:
                        self.assertIsNone(ids)

    def test_fold_without_eligible_proteins_does_not_produce_a_performance_score(self):
        data = pd.DataFrame({"constant": [3.] * 10, "group": self.y})
        with TemporaryDirectory() as output:
            logistic.run_logistic_regression(data, "group", ["constant"], 0, 1, output)
            self.assertEqual(list(Path(output).iterdir()), [])
        self.assertIn("No proteins passed", self.output.getvalue())
        self.assertNotIn("Accuracy:", self.output.getvalue())
