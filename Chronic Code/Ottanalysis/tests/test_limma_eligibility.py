import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from limma import run_limma


def comparison_results(data, proteins, *, paired, covariates=()):
    config = {
        "ID_COLUMN": "Sample_ID",
        "ID_DELIMITER": "-",
        "COVARIATE_COLUMNS": covariates,
        "P_THRESHOLD": 0.05,
        "LOG_FC_THRESHOLD": 0.58,
    }
    with TemporaryDirectory() as output:
        run_limma(data, "level", proteins, "ref", "case", output, paired, False, config)
        workbook = Path(output) / "Limma_ref_vs_case_Results.xlsx"
        return pd.read_excel(workbook, sheet_name="All_Proteins", index_col="Protein")


class LimmaEligibilityTests(unittest.TestCase):
    def test_only_observed_pairs_contribute(self):
        missing = np.nan
        data = pd.DataFrame({
            "Sample_ID": ["a-ref", "a-case", "b-ref", "b-case", "c-ref", "c-case", "d-ref", "d-case"],
            "level": ["ref", "case"] * 4,
            "no_pairs": [10, missing, missing, 25, 10, missing, missing, 25],
            "one_pair": [10, 11, missing, 25, 10, missing, missing, 25],
            "two_pairs": [10, 11, 19, 22, 90, missing, missing, 80],
        }).sample(frac=1, random_state=12)
        results = comparison_results(data, ["no_pairs", "one_pair", "two_pairs"], paired=True)

        self.assertEqual(list(results.index), ["two_pairs"])
        self.assertAlmostEqual(results.at["two_pairs", "logFC"], (1 + 3) / 2)
        self.assertAlmostEqual(results.at["two_pairs", "AveExpr"], (10 + 11 + 19 + 22) / 4)
        self.assertTrue(np.isfinite(results.at["two_pairs", "adj.P.Val"]))

    def test_one_fully_observed_protein_keeps_its_identifier(self):
        data = pd.DataFrame({
            "Sample_ID": ["a-ref", "a-case", "b-ref", "b-case", "c-ref", "c-case"],
            "level": ["ref", "case"] * 3,
            "sex": ["F", "F", "M", "M", "F", "F"],
            "signal": [20., 22., 16., 17., 12., 15.],
            "reference_only": [11., np.nan, 12., np.nan, 13., np.nan],
        })
        results = comparison_results(data, ["signal", "reference_only"], paired=True, covariates=["sex"])

        self.assertEqual(list(results.index), ["signal"])
        self.assertAlmostEqual(results.at["signal", "logFC"], 2)

    def test_replication_estimability_and_residual_df_with_covariates(self):
        missing = np.nan
        data = pd.DataFrame({
            "Sample_ID": list("abcdefgh"), "level": ["ref"] * 4 + ["case"] * 4,
            "batch": [0, 0, 1, 1, 0, 0, 1, 1], "order": [0, 1, 0, 1, 1, 0, 1, 0],
            "valid": [9, 10, 12, 11, 17, 15, 16, 18],
            "one_group": [9, 10, 12, 11, missing, missing, missing, missing],
            "one_case": [9, 10, 12, 11, 17, missing, missing, missing],
            # In these observed rows, batch and group are identical.
            "confounded": [9, 10, missing, missing, missing, missing, 16, 18],
            # Four observations for four independent design columns.
            "saturated": [9, 10, missing, missing, 17, missing, 16, missing],
        })
        proteins = ["valid", "one_group", "one_case", "confounded", "saturated"]
        results = comparison_results(data, proteins, paired=False, covariates=["batch", "order"])

        self.assertEqual(list(results.index), ["valid"])
        self.assertAlmostEqual(results.at["valid", "logFC"], 6)
        self.assertTrue(np.isfinite(results.at["valid", "P.Value"]))
