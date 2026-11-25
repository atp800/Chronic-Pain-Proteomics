import unittest
import warnings
from unittest.mock import patch

import numpy as np
import pandas as pd
import statsmodels.api as sm
from threadpoolctl import threadpool_limits

import hinge_inference as inference
import hinge_model as hinge
from hinge_solver import fit_random_intercept


def observed_frame(seed=83, *, missing=False, intercept_sd=1.2):
    rng = np.random.default_rng(seed)
    times = np.tile([0., 2., 4., 6., 8.], 8)
    groups = np.repeat(np.arange(8), 5)
    values = (12 + .4 * times - .7 * np.maximum(times - 4, 0)
              + rng.normal(0, intercept_sd, 8)[groups] + rng.normal(0, .4, len(times)))
    frame = pd.DataFrame({'value': values, 'time_numeric': times, 'subject_id_col': groups})
    if missing:
        # One singleton, one absent participant, and unequal numbers of visits.
        frame = frame.drop([0, 1, 3, 4, 5, 6, 7, 8, 9, 12, 18, 26])
    return frame


def at_peak(frame, peak=4.):
    return frame.assign(slope1=frame.time_numeric - peak,
                        slope2_change=np.maximum(frame.time_numeric - peak, 0))


def reference_fit(frame, *, hinge=True):
    terms = ['slope1', 'slope2_change'] if hinge else ['slope1']
    model = sm.MixedLM(frame.value, sm.add_constant(frame[terms]), groups=frame.subject_id_col)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fit = model.fit(method='bfgs', reml=False, gtol=1e-10, maxiter=1000)
        if not fit.converged:
            fit = model.fit(method='powell', reml=False, maxiter=1000, ftol=1e-10)
    if not fit.converged or not np.isfinite(fit.llf) or not np.isfinite(fit.bse_fe).all():
        raise AssertionError('MixedLM reference fit is unusable')
    return fit


class RandomInterceptSolverTests(unittest.TestCase):
    def setUp(self):
        self.enterContext(threadpool_limits(limits=1))

    def test_likelihood_coefficients_and_wald_errors_match_mixedlm(self):
        for missing, intercept_sd in [(False, 1.2), (True, 1.2), (True, .15), (False, 8.)]:
            frame = at_peak(observed_frame(missing=missing, intercept_sd=intercept_sd))
            for use_hinge in [False, True]:
                with self.subTest(missing=missing, intercept_sd=intercept_sd, hinge=use_hinge):
                    fast = fit_random_intercept(frame, hinge=use_hinge)
                    reference = reference_fit(frame, hinge=use_hinge)
                    self.assertAlmostEqual(fast.llf, reference.llf, delta=2e-6)
                    self.assertAlmostEqual(fast.aic, reference.aic, delta=4e-6)
                    np.testing.assert_allclose(fast.params.iloc[:-1], reference.fe_params, rtol=1e-4, atol=2e-5)
                    np.testing.assert_allclose(fast.bse_fe, reference.bse_fe, rtol=2e-4, atol=2e-5)
                    np.testing.assert_allclose(fast.pvalues, reference.pvalues.iloc[:-1], rtol=2e-3, atol=2e-5)
                    np.testing.assert_allclose(fast.scale, reference.scale, rtol=2e-4, atol=1e-6)
                    np.testing.assert_allclose(fast.cov_re, reference.cov_re, rtol=2e-4, atol=1e-6)
                    design = sm.add_constant(frame[['slope1', 'slope2_change'] if use_hinge else ['slope1']])
                    np.testing.assert_allclose(fast.predict(frame), reference.predict(design), rtol=1e-5, atol=2e-5)

    def test_zero_variance_boundary_matches_direct_gaussian_likelihood(self):
        frame = at_peak(observed_frame())
        residuals = np.random.default_rng(4).normal(0, .4, (8, 5))
        residuals -= residuals.mean(axis=1, keepdims=True)
        frame['value'] = 12 + .4 * frame.slope1 - .7 * frame.slope2_change + residuals.ravel()
        fit = fit_random_intercept(frame)
        design = sm.add_constant(frame[['slope1', 'slope2_change']]).to_numpy()
        coefficients = np.linalg.lstsq(design, frame.value, rcond=None)[0]
        residual = frame.value.to_numpy() - design @ coefficients
        scale = residual @ residual / len(frame)
        llf = -.5 * len(frame) * (np.log(2 * np.pi * scale) + 1)
        self.assertLess(fit.cov_re.item(), 1e-9)
        self.assertAlmostEqual(fit.llf, llf, delta=1e-9)
        np.testing.assert_allclose(fit.params.iloc[:-1], coefficients, atol=1e-10)
        np.testing.assert_allclose(fit.bse_fe, np.sqrt(np.diag(scale * np.linalg.inv(design.T @ design))), rtol=1e-8)

    def test_row_order_and_outcome_units_preserve_the_fit(self):
        frame = at_peak(observed_frame(missing=True))
        original = fit_random_intercept(frame)
        changed = frame.sample(frac=1, random_state=32)
        changed['value'] = changed.value * 1000 + 1e6
        fit = fit_random_intercept(changed)
        expected = original.params.iloc[:-1].to_numpy() * 1000
        expected[0] += 1e6
        np.testing.assert_allclose(fit.params.iloc[:-1], expected, rtol=1e-8)
        np.testing.assert_allclose(fit.bse_fe, original.bse_fe * 1000, rtol=1e-6)
        self.assertAlmostEqual(fit.llf, original.llf - len(frame) * np.log(1000), delta=1e-7)

    def test_winning_hinge_and_bootstrap_likelihood_gains_match(self):
        frame = observed_frame(missing=True)
        peaks = [2., 4., 6.]
        fast = [fit_random_intercept(at_peak(frame, peak)) for peak in peaks]
        reference = [reference_fit(at_peak(frame, peak)) for peak in peaks]
        self.assertEqual(np.argmin([fit.aic for fit in fast]), np.argmin([fit.aic for fit in reference]))
        null = fit_random_intercept(inference._null_frame(frame), hinge=False)
        codes, names = pd.factorize(frame.subject_id_col, sort=True)
        rng = np.random.default_rng(46)
        for _ in range(6):
            simulated = frame.copy()
            simulated['value'] = (null.predict(inference._null_frame(frame))
                                  + rng.normal(0, np.sqrt(null.cov_re.item()), len(names))[codes]
                                  + rng.normal(0, np.sqrt(null.scale), len(frame)))
            fast_lr = inference._refit_search(simulated, peaks, hinge._fit_mixed_model)
            reference_lr = inference._refit_search(
                simulated, peaks, lambda data, **kw: (reference_fit(data, **kw), 'bfgs', 'reference'))
            self.assertAlmostEqual(fast_lr, reference_lr, delta=1e-5)

    def test_bootstrap_pvalue_and_failure_counts_match_on_the_same_draws(self):
        frame = observed_frame(missing=True)
        peaks = [2., 4.]
        likelihoods = {peak: fit_random_intercept(at_peak(frame, peak)).llf for peak in peaks}
        null = fit_random_intercept(inference._null_frame(frame), hinge=False)
        calls = []
        for simulated_fitter in [hinge._fit_mixed_model,
                                 lambda data, **kw: (reference_fit(data, **kw), 'bfgs', 'reference')]:
            # Hold the fitted simulation distribution fixed to isolate refitting differences.
            def fit(data, **kwargs):
                if not calls:
                    calls.append('observed null')
                    return null, 'profile_ml', 'accepted'
                return simulated_fitter(data, **kwargs)
            calls.clear()
            result = inference.bootstrap_hinge_test(frame, likelihoods, fit, protein='agreement',
                                                    n_bootstrap=9, alpha=.99, show_progress=False)
            if simulated_fitter is hinge._fit_mixed_model:
                expected = result
            else:
                for column in ['Hinge_Search_PValue', 'Bootstrap_Exceedances', 'Bootstrap_Draws', 'Bootstrap_Failures']:
                    self.assertEqual(result[column], expected[column])

    def test_normal_entry_uses_profile_and_failure_uses_existing_optimisers(self):
        frame = at_peak(observed_frame())
        with patch.object(hinge, '_fit_statsmodels_model', side_effect=AssertionError('unexpected fallback')):
            fit, optimizer, _ = hinge._fit_mixed_model(frame)
        self.assertTrue(fit.converged)
        self.assertEqual(optimizer, 'profile_ml')
        with patch.object(hinge, 'fit_random_intercept', side_effect=ValueError('forced failure')):
            fallback, optimizer, diagnostics = hinge._fit_mixed_model(frame)
        self.assertIsNotNone(fallback)
        self.assertIn(optimizer, ['lbfgs', 'bfgs', 'cg'])
        self.assertIn('forced failure', diagnostics)

    def test_invalid_designs_raise_for_the_caller_to_handle(self):
        frame = at_peak(observed_frame())
        cases = [frame.assign(value=np.nan), frame.assign(value=4.), frame.iloc[:3],
                 frame.assign(slope2_change=frame.slope1), frame.assign(subject_id_col=None),
                 frame.assign(subject_id_col=np.arange(len(frame)))]
        for data in cases:
            with self.subTest(rows=len(data)), self.assertRaises(ValueError):
                fit_random_intercept(data)
