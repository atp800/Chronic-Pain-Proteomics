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
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import GridSearchCV, StratifiedKFold
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
        observed_aucs = []
        observed = pd.Series(self.y, index=self.X.index)
        # With only "signal" in this comparison, p3 has no measurement at its first visit.
        eligible_ids = self.participants[self.participants.ne("p3")]
        fit_features = ProteinFeatures.fit
        permutation_test = logistic._nested_permutation_test

        def inspect_fit(transformer, X, y):
            labels = pd.Series(y, index=X.index)
            ids = self.participants.loc[X.index]
            for _, pair in labels.groupby(ids):
                self.assertCountEqual(pair, [0, 1])
            changed = labels.ne(observed.loc[X.index])
            self.assertTrue(changed.groupby(ids).nunique().eq(1).all())
            swaps.append(changed.any() and not changed.all())
            return fit_features(transformer, X, y)

        def run_permutations(X, y, delta, ids, observed_auc, **kwargs):
            observed_aucs.append(observed_auc)
            self.assertTrue(X.index.equals(eligible_ids.index))
            np.testing.assert_array_equal(ids, eligible_ids.to_numpy())
            with patch.object(ProteinFeatures, "fit", inspect_fit):
                return permutation_test(X, y, delta, ids, observed_auc, n_permutations=5)

        data = self.X.assign(visit=self.y)
        with (TemporaryDirectory() as output, patch.object(logistic, "_nested_permutation_test", run_permutations),
              patch.object(logistic, "roc_curve", wraps=roc_curve) as plotted_roc,
              patch.object(logistic, "_bootstrap_rows", wraps=logistic._bootstrap_rows) as bootstrap):
            logistic.run_logistic_regression(data, "visit", ["signal"], 0, 1, output, participant_ids=self.participants)
            np.testing.assert_array_equal(bootstrap.call_args.args[1], eligible_ids.to_numpy())
            exported = pd.read_csv(Path(output) / "LogReg_1_vs_0_Features.csv")
            self.assertIn("Bootstrap 2.5th Percentile", exported)
            self.assertIn("Bootstrap 97.5th Percentile", exported)
            self.assertNotIn("Significant CI", exported)
            self.assertTrue((Path(output) / "LogReg_1_vs_0_Plots.pdf").is_file())
            self.assertEqual(observed_aucs, [roc_auc_score(*plotted_roc.call_args.args)])
            self.assertIn("Nested CV permutation p-value (5 draws)", self.output.getvalue())
        self.assertTrue(any(swaps))

    def test_group_and_delta_shuffles_retune_the_entire_grid_and_score_held_out_rows(self):
        raw = self.X.assign(on_off=[2., 4, 6, 8, 10] + [np.nan] * 5)
        baseline = pd.DataFrame(10., index=raw.index, columns=raw.columns)
        endpoints = pd.concat({"baseline": baseline, "followup": baseline + raw}, axis=1)
        grid = {"classifier__C": [0.05, 2.], "classifier__l1_ratio": [0.1, 1.]}
        rng = np.random.RandomState(42)
        assignments = [rng.permutation(self.y) for _ in range(2)]
        real_fit = logistic._fit_tuned_model

        for X, delta in [(raw, False), (endpoints, True)]:
            with self.subTest(delta=delta), patch.object(logistic, "_PARAMETER_GRID", grid):
                observed_probs = logistic._leave_one_out_probabilities(X, self.y, delta)
                observed_auc = roc_auc_score(self.y, observed_probs[:, 1])
                predictions, retained = [], []

                def inspect_fit(train, labels, is_delta, ids):
                    draw = assignments[len(predictions) // len(X)]
                    np.testing.assert_array_equal(labels, draw[X.index.get_indexer(train.index)])
                    pd.testing.assert_frame_equal(train, X.loc[train.index])
                    self.assertEqual(is_delta, delta)
                    self.assertIsNone(ids)
                    model = real_fit(train, labels, is_delta, ids)
                    marker = raw.loc[train.index, "on_off"].notna()
                    keep = marker.groupby(labels).mean().ge(0.7).any()
                    self.assertEqual("on_off" in model.named_steps["proteins"].proteins_, keep)
                    retained.append(keep)
                    test = X.loc[~X.index.isin(train.index)]
                    predictions.append((test.index.item(), model.predict_proba(test)[0, 1]))
                    return model

                with (patch.object(logistic, "_fit_tuned_model", inspect_fit),
                      patch.object(logistic, "GridSearchCV", wraps=logistic.GridSearchCV) as searches):
                    scores, pvalue = logistic._nested_permutation_test(X, self.y, delta, None, observed_auc,
                                                                       n_permutations=2)
                self.assertEqual(searches.call_count, len(assignments) * len(X))
                for call in searches.call_args_list:
                    self.assertEqual(call.args[1], grid)
                self.assertIn(False, retained)  # the full-data fit retains this marker under the real labels
                expected = []
                for i, labels in enumerate(assignments):
                    held_out = dict(predictions[i * len(X):(i + 1) * len(X)])
                    expected.append(roc_auc_score(labels, [held_out[row] for row in X.index]))
                np.testing.assert_allclose(scores, expected)
                self.assertAlmostEqual(pvalue, (1 + sum(s >= observed_auc - 1e-14 for s in expected)) / 3)

    def test_repeated_paired_assignments_count_every_draw_including_ties(self):
        y = np.array([0, 1] * 3)
        ids = np.repeat(["a", "b", "c"], 2)

        def predict(X, labels, delta, participant_ids, *, show_progress):
            # Perfect or reversed rankings, depending on the first pair's assignment.
            positive = labels if labels[0] == 0 else 1 - labels
            return np.column_stack([1 - positive, positive])

        with patch.object(logistic, "_leave_one_out_probabilities", side_effect=predict) as fit:
            scores, pvalue = logistic._nested_permutation_test(self.X.iloc[:6], y, False, ids, 1., n_permutations=100)
        self.assertEqual(len(scores), 100)
        self.assertEqual(set(scores), {0., 1.})
        self.assertEqual(fit.call_count, 7)  # eight possible assignments, observed score already known
        self.assertEqual(pvalue, (1 + np.count_nonzero(scores == 1.)) / 101)

    def test_no_permutation_pvalue_when_a_shuffled_fold_has_no_eligible_proteins(self):
        data = pd.DataFrame({"marker": [8., 9, 10, 11] + [np.nan] * 4, "group": [0] * 4 + [1] * 4})
        with TemporaryDirectory() as output, patch.object(logistic.plt, "hist") as histogram:
            logistic.run_logistic_regression(data, "group", ["marker"], 0, 1, output)
            self.assertTrue((Path(output) / "LogReg_1_vs_0_Features.csv").is_file())
            self.assertTrue((Path(output) / "LogReg_1_vs_0_Plots.pdf").is_file())
            histogram.assert_not_called()
        self.assertIn("Permutation result unavailable: No proteins passed", self.output.getvalue())
        self.assertNotIn("permutation p-value", self.output.getvalue())

    def test_eight_vs_four_uses_feasible_inner_folds_and_a_pooled_auc(self):
        y = np.repeat([0, 1], [8, 4])
        X = pd.DataFrame({"protein": [5., 8, 6, 9, 3, 7, 11, 4, 8, 5, 10, 6]})
        real_search = logistic.GridSearchCV
        folds = []

        def inspect_search(*args, **kwargs):
            folds.append(kwargs["cv"])
            return real_search(*args, **kwargs)

        with patch.object(logistic, "GridSearchCV", inspect_search):
            probs = logistic._leave_one_out_probabilities(X, y, False)[:, 1]
        self.assertEqual(folds, [4] * 8 + [3] * 4)
        self.assertTrue(np.isfinite(probs).all())
        differences = probs[y == 1, None] - probs[y == 0]
        auc = ((differences > 0).sum() + 0.5 * (differences == 0).sum()) / 32
        self.assertAlmostEqual(auc, roc_auc_score(y, probs))
        scores, pvalue = logistic._nested_permutation_test(X, y, False, None, auc, n_permutations=3)
        self.assertTrue(np.isfinite(scores).all())
        self.assertAlmostEqual(pvalue, (1 + (scores >= auc).sum()) / 4)

    def test_undefined_auc_on_a_later_draw_invalidates_the_whole_permutation_test(self):
        valid = np.full((len(self.y), 2), 0.5)
        invalid = valid.copy()
        invalid[2] = np.nan
        with (patch.object(logistic, "_leave_one_out_probabilities", side_effect=[valid, invalid]),
              self.assertRaisesRegex(ValueError, "permuted held-out AUC was undefined")):
            logistic._nested_permutation_test(self.X, self.y, False, None, 0.5, n_permutations=2)

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
                      RUN_LIMMA=False, RUN_LOGISTIC_REGRESSION=True, RUN_HINGE=False, LOGISTIC_PERMUTATIONS=1234)
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
                    self.assertEqual(call.kwargs['n_permutations'], 1234)
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
