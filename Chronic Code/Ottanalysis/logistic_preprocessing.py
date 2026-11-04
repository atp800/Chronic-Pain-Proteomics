import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted

from helpers import IMPUTATION_SEED, IMPUTATION_SHIFT, IMPUTATION_WIDTH, groupwise_missing_filter


class NoEligibleProteins(ValueError):
    """A training fold cannot support a protein classifier."""


class ProteinFeatures(TransformerMixin, BaseEstimator):
    """Fit protein selection and abundance imputation on one training fold."""

    def __init__(self, delta=False):
        self.delta = delta

    def fit(self, X, y):
        observed = X["followup"] - X["baseline"] if self.delta else X
        eligible = groupwise_missing_filter(observed, np.asarray(y), threshold=0.7, verbose=False)
        variable = eligible.fillna(0).var(ddof=0).gt(0.1)
        self.proteins_ = eligible.columns[variable]
        if self.proteins_.empty:
            raise NoEligibleProteins("No proteins passed the training-fold coverage and variance filters.")

        if self.delta:
            abundances = pd.concat([X["baseline"][self.proteins_], X["followup"][self.proteins_]])
        else:
            abundances = X[self.proteins_]
        # Preserve the existing zero fallback when training spread is unestimable or zero.
        sd = abundances.std().fillna(0)
        mean = (abundances.mean() - IMPUTATION_SHIFT * sd).where(sd.gt(0), 0)
        self.fill_mean_ = mean.to_numpy(dtype=float)
        self.fill_sd_ = (IMPUTATION_WIDTH * sd).to_numpy(dtype=float)
        return self

    def _abundances(self, data, visit):
        values = data[self.proteins_].to_numpy(dtype=float, na_value=np.nan, copy=True)
        missing = np.isnan(values)
        if missing.any():
            # Stable draws per sample and visit, even when predicting one row at
            # a time. Other validation rows cannot change this sample's values.
            row_keys = pd.util.hash_pandas_object(data.index, index=False).to_numpy()
            for i in np.flatnonzero(missing.any(axis=1)):
                rng = np.random.default_rng([IMPUTATION_SEED, visit, int(row_keys[i])])
                draws = rng.normal(self.fill_mean_, self.fill_sd_)
                values[i, missing[i]] = draws[missing[i]]
        return values

    def transform(self, X):
        check_is_fitted(self, "proteins_")
        if self.delta:
            return self._abundances(X["followup"], visit=1) - self._abundances(X["baseline"], visit=0)
        return self._abundances(X, visit=0)
