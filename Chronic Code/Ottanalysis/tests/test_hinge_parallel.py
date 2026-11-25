import io
import os
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pandas as pd
from joblib import Parallel, delayed, effective_n_jobs, parallel_config
from threadpoolctl import threadpool_info, threadpool_limits

import hinge_model as hinge
from test_hinge_fallback import repeated_data


def inspect_protein_worker(frame, index):
    """Inspect libraries inside actual fits, including bootstrap refits."""
    counts = []
    original = hinge._fit_mixed_model

    def record_threads(*args, **kwargs):
        counts.extend(library['num_threads'] for library in threadpool_info())
        return original(*args, **kwargs)

    with patch.object(hinge, '_fit_mixed_model', side_effect=record_threads):
        result = hinge._fit_hinge_protein(index, f'p{index}', frame.copy(), [2., 4.],
                                         {2.: 'D2', 4.: 'D4'}, 8, 3, 42, .99)
    return os.getpid(), counts, result


class HingeParallelTests(unittest.TestCase):
    def setUp(self):
        self.enterContext(redirect_stdout(io.StringIO()))
        self.folder = Path(self.enterContext(TemporaryDirectory()))
        self.data = repeated_data()
        self.data['other'] = self.data.mixed + np.random.default_rng(5).normal(0, .1, len(self.data))
        self.data['missing'] = self.data.mixed
        self.data.loc[[0, 1, 3, 4, 10, 18], 'missing'] = np.nan
        self.data['descriptive'] = self.data.mixed.where(self.data.index.isin([i * 5 + i % 5 for i in range(8)]))
        self.data['absent'] = np.nan
        self.data['flat'] = 9.
        self.proteins = ['mixed', 'other', 'missing', 'descriptive', 'absent', 'flat']

    def test_workbooks_match_across_worker_counts_and_protein_order(self):
        tables = []
        original = self.data.copy(deep=True)
        for workers, proteins in [(1, self.proteins), (2, self.proteins), (3, self.proteins[::-1])]:
            directory = self.folder / str(workers)
            directory.mkdir()
            path = hinge.run_hinge_model(self.data, proteins, 'Time', 'Sample_ID', ['D2', 'D4'],
                                         0, directory, '-', n_bootstrap=9, fdr_threshold=.99, n_jobs=workers)
            workbook = pd.read_excel(path, sheet_name=None)
            settings = workbook['Analysis_Info'].set_index('Setting').Value
            self.assertEqual(int(settings['Protein workers']), min(workers, effective_n_jobs(-1)))
            self.assertEqual(int(settings['Numerical-library threads per worker']), 1)
            tables.append(workbook)

        for sheet, serial in tables[0].items():
            if sheet != 'Analysis_Info':
                pd.testing.assert_frame_equal(serial, tables[1][sheet], check_exact=True)
        for sheet in ['All_Model_Results', 'Data_Coverage', 'All_Candidate_Timepoints', 'Observed_By_Timepoint']:
            order = ['Protein']
            if sheet == 'All_Candidate_Timepoints':
                order.append('Candidate_Time_Numeric')
            if sheet == 'Observed_By_Timepoint':
                order.append('time_numeric')
            pd.testing.assert_frame_equal(tables[0][sheet].sort_values(order).reset_index(drop=True),
                                          tables[2][sheet].sort_values(order).reset_index(drop=True), check_exact=True)
        results = tables[0]['All_Model_Results'].set_index('Protein')
        self.assertTrue(results.loc[['mixed', 'other', 'missing'], 'Bootstrap_Draws'].gt(0).all())
        self.assertEqual(results.loc['descriptive', 'Model_Type'], 'Fixed-Effects')
        self.assertNotIn('absent', results.index)
        self.assertNotIn('flat', results.index)
        self.assertEqual(tables[0]['Data_Coverage'].Protein.tolist(), self.proteins)
        pd.testing.assert_frame_equal(self.data, original)

    def test_workers_and_serial_fits_limit_numeric_threads_to_one(self):
        frame = self.data.rename(columns={'mixed': 'value'}).copy()
        frame['subject_id_col'] = frame.Sample_ID.str.split('-').str[0]
        frame['time_numeric'] = frame.Time.str[1:].astype(float)
        with threadpool_limits(limits=2):
            pid, counts, _ = inspect_protein_worker(frame, 0)
            self.assertEqual(pid, os.getpid())
            self.assertTrue(counts)
            self.assertEqual(set(counts), {1})
        if effective_n_jobs(-1) >= 2:
            with parallel_config(backend='loky', inner_max_num_threads=1):
                results = Parallel(n_jobs=2, batch_size=1)(delayed(inspect_protein_worker)(frame, i) for i in range(4))
            self.assertEqual(len({pid for pid, _, _ in results}), 2)
            for pid, counts, result in results:
                self.assertNotEqual(pid, os.getpid())
                self.assertTrue(counts)
                self.assertEqual(set(counts), {1})
                self.assertIsNotNone(result[1])

    def test_worker_error_preserves_existing_output(self):
        existing = self.folder / 'hinge_model_results.xlsx'
        existing.write_bytes(b'previous completed results')

        def fail(*args, **kwargs):
            raise RuntimeError('forced worker failure')

        with patch.object(hinge, '_fit_hinge_protein', new=fail), self.assertRaisesRegex(RuntimeError, 'forced worker failure'):
            hinge.run_hinge_model(self.data, ['mixed', 'other'], 'Time', 'Sample_ID', ['D4'],
                                  0, self.folder, '-', n_bootstrap=3, n_jobs=2)
        self.assertEqual(existing.read_bytes(), b'previous completed results')
        self.assertEqual(list(self.folder.iterdir()), [existing])

    def test_invalid_worker_count_is_rejected_before_processing(self):
        for workers in [0, True, 1.5, None]:
            with self.subTest(workers=workers), self.assertRaisesRegex(ValueError, 'Hinge workers'):
                hinge.run_hinge_model(self.data, ['mixed'], 'Time', 'Sample_ID', ['D4'],
                                      0, self.folder, '-', n_bootstrap=3, n_jobs=workers)
