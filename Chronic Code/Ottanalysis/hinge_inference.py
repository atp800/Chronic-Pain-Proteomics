"""Whole-search parametric bootstrap for the random-intercept hinge model.

All fits use ML, the same observed rows and the same candidate days. Simulations
come from a no-bend model with a linear time trend, fresh participant intercepts
and independent residuals. This is model-based calibration, conditional on the
observed measurement pattern; it does not correct non-random missingness.
"""

import hashlib
from numbers import Integral

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm


DEFAULT_BOOTSTRAPS = 9999
DEFAULT_SEED = 42


def validate_test_settings(n_bootstrap, seed, alpha):
    if isinstance(n_bootstrap, bool) or not isinstance(n_bootstrap, Integral) or n_bootstrap < 1:
        raise ValueError('Hinge bootstrap simulations must be a positive integer')
    if isinstance(seed, bool) or not isinstance(seed, Integral) or seed < 0:
        raise ValueError('Hinge bootstrap seed must be a nonnegative integer')
    if not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError('Hinge FDR threshold must be between 0 and 1, exclusive')


def unavailable_test(n_bootstrap, seed, reason):
    return {
        'Hinge_Search_LR': np.nan,
        'Hinge_Search_PValue': np.nan,
        'Search_PValue_Is_Upper_Bound': False,
        'Bootstrap_Target': n_bootstrap,
        'Bootstrap_Draws': 0,
        'Bootstrap_Exceedances': 0,
        'Bootstrap_Failures': 0,
        'Bootstrap_Seed': seed,
        'Null_Optimizer': None,
        'Null_Fit_Diagnostics': '',
        'Hinge_Test_Status': reason,
    }


def _likelihood_gain(null, alternative_llfs):
    values = np.asarray(alternative_llfs, dtype=float)
    if not len(values) or not np.isfinite(values).all() or not np.isfinite(null.llf):
        raise ValueError('Incomplete or non-finite mixed-model search')
    gain = 2 * (values.max() - null.llf)
    # Nested alternatives cannot have a worse maximum likelihood than the null.
    if gain < -1e-7 * (1 + abs(null.llf)):
        raise ValueError('Hinge likelihood below the nested null; optimisation unreliable')
    return max(0., gain)


def _null_frame(observed):
    frame = observed.copy()
    frame['slope1'] = frame['time_numeric'] - frame['time_numeric'].mean()
    return frame


def _refit_search(frame, peaks, fit_mixed):
    null, _, diagnostics = fit_mixed(_null_frame(frame), hinge=False)
    if null is None:
        raise ValueError('Simulated null fit failed: ' + diagnostics)
    likelihoods = []
    for peak in peaks:
        candidate = frame.copy()
        candidate['slope1'] = candidate['time_numeric'] - peak
        candidate['slope2_change'] = candidate['slope1'].clip(lower=0)
        fitted, _, diagnostics = fit_mixed(candidate)
        if fitted is None:
            raise ValueError(f'Simulated hinge fit failed at {peak}: {diagnostics}')
        likelihoods.append(fitted.llf)
    return _likelihood_gain(null, likelihoods)


def bootstrap_hinge_test(observed, candidate_llfs, fit_mixed, *, protein,
                         n_bootstrap=DEFAULT_BOOTSTRAPS, seed=DEFAULT_SEED,
                         alpha=.05, show_progress=True):
    """Calibrate the best-day LR; failed/uncomputed draws count as exceedances.

    candidate_llfs must include every design-eligible day, with NaN for failed
    observed fits. Do not select a smaller search after looking at fit outcomes.
    AIC selection and maximum likelihood selection coincide here because all
    candidate hinges have the same parameter count and use identical rows.
    """
    validate_test_settings(n_bootstrap, seed, alpha)
    result = unavailable_test(n_bootstrap, seed, 'Search test unavailable')
    if not candidate_llfs or not np.isfinite(list(candidate_llfs.values())).all():
        result['Hinge_Test_Status'] = 'Unavailable: at least one eligible observed hinge fit failed'
        return result
    frame = _null_frame(observed)
    null, optimizer, diagnostics = fit_mixed(frame, hinge=False)
    result.update(Null_Optimizer=optimizer, Null_Fit_Diagnostics=diagnostics)
    if null is None:
        result['Hinge_Test_Status'] = 'Unavailable: no usable no-bend mixed fit'
        return result
    try:
        observed_lr = _likelihood_gain(null, list(candidate_llfs.values()))
        variance = float(np.asarray(null.cov_re)[0, 0])
        residual_variance = float(null.scale)
        # predict supplies the population mean, NOT fitted participant BLUPs.
        mean = np.asarray(null.predict(frame), dtype=float)
        if (not np.isfinite([variance, residual_variance]).all() or variance < 0
                or residual_variance <= 0 or not np.isfinite(mean).all()):
            raise ValueError('Invalid null simulation parameters')
    except (ValueError, FloatingPointError) as exc:
        result['Hinge_Test_Status'] = f'Unavailable: {exc}'
        return result
    result['Hinge_Search_LR'] = observed_lr
    if observed_lr == 0:
        result.update(Hinge_Search_PValue=1., Hinge_Test_Status='No improvement over no-bend model')
        return result

    # Stable per-protein streams: reordering/filtering other proteins does not change a test.
    tag = int.from_bytes(hashlib.sha256(str(protein).encode('utf-8')).digest()[:8], 'little')
    rng = np.random.default_rng(np.random.SeedSequence([seed, tag]))
    groups, names = pd.factorize(frame['subject_id_col'], sort=True)
    exceedances = failures = 0
    first_failure = ''
    iterator = tqdm(range(n_bootstrap), desc=f'Bootstrap {protein}', leave=False, disable=not show_progress)
    for draw in iterator:
        simulated = frame.copy()
        intercepts = rng.normal(0., np.sqrt(variance), len(names))
        simulated['value'] = mean + intercepts[groups] + rng.normal(0., np.sqrt(residual_variance), len(frame))
        try:
            score = _refit_search(simulated, list(candidate_llfs), fit_mixed)
            if not np.isfinite(score):
                raise ValueError('Non-finite simulated statistic')
        except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
            failures += 1
            exceedances += 1  # Worst-case treatment; never drop or replace a failed draw.
            if not first_failure:
                first_failure = str(exc)
        else:
            exceedances += int(score >= observed_lr or np.isclose(score, observed_lr, rtol=1e-10, atol=1e-12))
        completed = draw + 1
        # Even zero further exceedances cannot yield p <= alpha (and BH q >= p).
        # Stop only then, reporting the upper bound obtained by counting EVERY
        # uncomputed draw as an exceedance. Never report a stopped point estimate.
        if (exceedances + 1) / (n_bootstrap + 1) > alpha:
            break
    uncomputed = n_bootstrap - completed
    result.update(
        Bootstrap_Draws=completed, Bootstrap_Exceedances=exceedances,
        Bootstrap_Failures=failures,
        Search_PValue_Is_Upper_Bound=bool(failures or uncomputed),
    )
    if failures == completed:
        result['Hinge_Test_Status'] = 'Unavailable: all attempted simulations failed; ' + first_failure
        return result
    result['Hinge_Search_PValue'] = (exceedances + uncomputed + 1) / (n_bootstrap + 1)
    result['Hinge_Test_Status'] = ('Conservative upper bound: stopped after ruling out significance'
                                   if uncomputed else 'Completed')
    if failures:
        result['Hinge_Test_Status'] += '; failed draws counted as exceedances; first failure: ' + first_failure
    return result


def adjust_hinge_fdr(results, family_size, alpha):
    """BH within one group/sheet; include skipped/unavailable proteins as p=1."""
    adjusted = results.copy()
    eligible = (adjusted['Model_Type'].eq('Mixed-Effects')
                & adjusted['Hinge_Search_PValue'].notna())
    pvalues = adjusted.loc[eligible, 'Hinge_Search_PValue'].to_numpy(dtype=float)
    if family_size < len(results) or not np.isfinite(pvalues).all() or ((pvalues < 0) | (pvalues > 1)).any():
        raise ValueError('Invalid hinge testing family or search p-values')
    adjusted['Hinge_FDR'] = np.nan
    if len(pvalues):
        family = np.r_[pvalues, np.ones(family_size - len(pvalues))]
        adjusted.loc[eligible, 'Hinge_FDR'] = multipletests(family, method='fdr_bh')[1][:len(pvalues)]
    adjusted['Significant_FDR'] = eligible & adjusted['Hinge_FDR'].le(alpha)
    return adjusted
