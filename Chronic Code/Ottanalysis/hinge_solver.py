"""Profiled maximum likelihood for a Gaussian random-intercept model.

For participant i, V_i = sigma² (I + lambda 11'), lambda >= 0. Within-person
contrasts have weight 1 and the participant mean has weight 1/(1+n_i*lambda).
GLS profiles out the fixed effects, and sigma² = residual quadratic form / n
profiles out the residual variance. Only log(1+lambda) needs optimisation.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import norm


@dataclass
class RandomInterceptFit:
    params: pd.Series
    bse_fe: pd.Series
    pvalues: pd.Series
    cov_re: np.ndarray
    scale: float
    llf: float
    aic: float
    converged: bool = True

    def predict(self, frame):
        """Population mean, used to simulate fresh participant effects."""
        terms = self.bse_fe.index[1:]
        return self.params['Intercept'] + frame[terms].to_numpy() @ self.params[terms].to_numpy()


def fit_random_intercept(frame, *, hinge=True):
    """Fit the same ML model and observed-information Wald errors as MixedLM.

    Inputs must contain only the protein's observed finite rows. Singular
    designs, degenerate residuals and unresolved optimisation raise ValueError
    or LinAlgError so the caller can try its existing MixedLM fallbacks.
    """
    terms = ['Intercept', 'slope1', 'slope2_change'] if hinge else ['Intercept', 'slope1']
    design = np.column_stack([np.ones(len(frame)), frame[terms[1:]].to_numpy(dtype=float)])
    values = frame['value'].to_numpy(dtype=float)
    groups, names = pd.factorize(frame['subject_id_col'], sort=True)
    n, p = design.shape
    if (n <= p or len(names) < 2 or (groups < 0).any()
            or not np.isfinite(design).all() or not np.isfinite(values).all()):
        raise ValueError('Insufficient or non-finite random-intercept data')
    if np.linalg.matrix_rank(design) != p:
        raise ValueError('Rank-deficient fixed-effect design')
    counts = np.bincount(groups).astype(float)
    if not (counts > 1).any():
        raise ValueError('Random-intercept variance requires repeated observations')

    # Centre/scale the outcome and orthogonalise the design for numerical stability.
    # Transform coefficients and their covariance back to the original units below.
    offset = values.mean()
    outcome_scale = values.std()
    if not np.isfinite(outcome_scale) or outcome_scale <= 0:
        raise ValueError('No observed abundance variation')
    y = (values - offset) / outcome_scale
    q, r = np.linalg.qr(design, mode='reduced')
    group_x = np.zeros((len(names), p))
    np.add.at(group_x, groups, q)
    group_x /= counts[:, None]
    group_y = np.bincount(groups, weights=y) / counts
    within_x = q - group_x[groups]
    within_y = y - group_y[groups]
    within_xx = within_x.T @ within_x
    within_xy = within_x.T @ within_y

    def evaluate(log_ratios):
        # Batch the initial profile scan into small matrix operations.
        ratios = np.expm1(log_ratios)
        weights = counts / (1 + counts * ratios[:, None])
        information = within_xx + np.einsum('gi,hg,gj->hij', group_x, weights, group_x)
        rhs = within_xy + (weights * group_y) @ group_x
        coefficients = np.linalg.solve(information, rhs[..., None])[..., 0]
        within_residual = within_y - coefficients @ within_x.T
        mean_residual = group_y - coefficients @ group_x.T
        residual_ss = (np.sum(within_residual ** 2, axis=1)
                       + np.sum(weights * mean_residual ** 2, axis=1))
        if not np.isfinite(residual_ss).all() or (residual_ss <= np.finfo(float).eps * n).any():
            raise ValueError('Degenerate residual variance in profiled likelihood')
        deviance = n * np.log(residual_ss / n) + np.log1p(counts * ratios[:, None]).sum(axis=1)
        return deviance, coefficients, residual_ss, information, mean_residual

    def objective(log_ratio):
        return evaluate(np.array([log_ratio]))[0].item()

    # A scalar profile can have several local minima for unequal group sizes.
    # Search every bracketed trough, explicitly compare the zero-variance boundary,
    # and refuse an optimum beyond the numerically supported range (lambda ~ 1e14).
    grid = np.linspace(0., 32., 65)
    deviances = evaluate(grid)[0]
    candidates = [(deviances[0], 0.)]
    for i in range(1, len(grid) - 1):
        if deviances[i] <= deviances[i - 1] and deviances[i] <= deviances[i + 1]:
            optimum = minimize_scalar(objective, bounds=(grid[i - 1], grid[i + 1]),
                                      method='bounded', options={'xatol': 1e-10, 'maxiter': 200})
            if not optimum.success or not np.isfinite(optimum.fun):
                raise ValueError('Profile likelihood optimisation did not converge')
            candidates.append((optimum.fun, optimum.x))
    # Include the first interval even when a shallow interior minimum lies close to zero.
    near_zero = minimize_scalar(objective, bounds=(grid[0], grid[1]), method='bounded',
                                options={'xatol': 1e-10, 'maxiter': 200})
    if not near_zero.success or not np.isfinite(near_zero.fun):
        raise ValueError('Boundary profile optimisation did not converge')
    candidates.append((near_zero.fun, near_zero.x))
    best_deviance, log_ratio = min(candidates)
    if deviances[0] - best_deviance <= 16 * np.finfo(float).eps * (1 + abs(deviances[0])):
        best_deviance, log_ratio = deviances[0], 0.
    if deviances[-1] <= best_deviance:
        raise ValueError('Variance-ratio optimum exceeds the profile search range')
    deviance, coefficients, residual_ss, information, mean_residual = (
        value[0] for value in evaluate(np.array([log_ratio]))
    )
    ratio = np.expm1(log_ratio)
    residual_scale = residual_ss / n

    # Observed information for (fixed effects, lambda), after profiling sigma².
    # The cross block accounts for uncertainty in lambda, matching MixedLM's
    # inverse joint Hessian rather than using a conditional GLS covariance.
    inverse_mean_variance = 1 / (1 + counts * ratio)
    weighted_counts = counts * inverse_mean_variance
    b = np.sum((weighted_counts * mean_residual) ** 2)
    gradient = (1 + ratio) * (weighted_counts.sum() - n * b / residual_ss)
    if ((ratio == 0 and gradient < -1e-6 * n)
            or (ratio > 0 and abs(gradient) > 1e-6 * n)):
        raise ValueError('Profile likelihood optimum failed the gradient check')
    d = 2 * np.sum(weighted_counts ** 3 * mean_residual ** 2)
    cross = group_x.T @ (weighted_counts ** 2 * mean_residual) / residual_scale
    curvature = -.5 * np.sum(weighted_counts ** 2) + .5 * n * (d / residual_ss - (b / residual_ss) ** 2)
    joint_information = np.empty((p + 1, p + 1))
    joint_information[:p, :p] = information / residual_scale
    joint_information[:p, p] = cross
    joint_information[p, :p] = cross
    joint_information[p, p] = curvature
    covariance = np.linalg.inv(joint_information)[:p, :p]
    transform = np.linalg.solve(r, np.eye(p)) * outcome_scale
    fixed_covariance = transform @ covariance @ transform.T
    variances = np.diag(fixed_covariance)
    if not np.isfinite(variances).all() or (variances <= 0).any():
        raise ValueError('Unusable fixed-effect standard errors')

    fixed_effects = transform @ coefficients
    fixed_effects[0] += offset
    scale = residual_scale * outcome_scale ** 2
    llf = -.5 * (deviance + n * (np.log(2 * np.pi) + 1 + 2 * np.log(outcome_scale)))
    params = pd.Series(np.r_[fixed_effects, ratio], index=terms + ['Group Var'])
    bse = pd.Series(np.sqrt(variances), index=terms)
    pvalues = pd.Series(2 * norm.sf(np.abs(fixed_effects / bse.to_numpy())), index=terms)
    if not np.isfinite(params).all() or not np.isfinite([scale, llf, scale * ratio]).all():
        raise ValueError('Non-finite fitted random-intercept parameters')
    return RandomInterceptFit(params=params, bse_fe=bse, pvalues=pvalues,
                              cov_re=np.array([[scale * ratio]]), scale=scale, llf=llf,
                              aic=-2 * llf + 2 * (p + 2))
