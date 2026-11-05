import io
from pathlib import Path
import runpy
import sys
from tempfile import TemporaryDirectory
import tkinter as tk
from types import SimpleNamespace
import unittest
from contextlib import redirect_stdout
from unittest.mock import Mock, patch

import pandas as pd

import gui
from config import CONFIG


class GUIColumnSelectionTests(unittest.TestCase):
    def setUp(self):
        self.folder = Path(self.enterContext(TemporaryDirectory()))
        self.data = pd.DataFrame([
            {'Sample_ID': f'p{person}-{visit}', 'Condition': 'Resp' if person > 2 else 'NonResp',
             'Time_Point': visit, 'Age': 30 + person, 'Gender': person % 2, 'Batch': person % 3,
             'protein_a': 8 + person + day / 10, 'protein_b': 12 - person + day / 20}
            for person in range(1, 5) for visit, day in [('T0', 0), ('D7', 7), ('D14', 14)]
        ])
        self.workbook = self.folder / 'data.xlsx'
        with pd.ExcelWriter(self.workbook) as writer:
            self.data.to_excel(writer, sheet_name='s', index=False)
            self.data[self.data.columns[::-1]].to_excel(writer, sheet_name='Reordered', index=False)
            self.data.drop(columns=['Age', 'Batch']).to_excel(writer, sheet_name='Missing', index=False)
        self.enterContext(patch.dict(CONFIG, {
            'INPUT_FILE_PATH': str(self.workbook), 'SHEET_NAME': 's',
            'OUTPUT_FOLDER': str(self.folder / 'output'), 'LAUNCH_GUI': True,
            'ID_COLUMN': 'Sample_ID', 'CONDITION_COLUMN': 'Condition', 'TIME_COLUMN': 'Time_Point',
            'COVARIATE_COLUMNS': ['Age', 'Gender'], 'UNNEEDED_COLUMNS': ['Batch'],
            'RUN_LIMMA': True, 'RUN_LOGISTIC_REGRESSION': True, 'RUN_HINGE': False,
        }))
        self.enterContext(redirect_stdout(io.StringIO()))
        self.errors = self.enterContext(patch.object(gui.messagebox, 'showerror'))
        self.app = self.new_app()

    def new_app(self):
        # Exercise real Tk selection behaviour without displaying a window.
        try:
            root = tk.Tk()
        except tk.TclError as exc:
            self.skipTest(f'Tk display unavailable: {exc}')
        root.withdraw()

        def close():
            try:
                root.destroy()
            except tk.TclError:
                pass  # on_submit already closed it.

        self.addCleanup(close)
        return gui.ProteomicsGUI(root)

    def select(self, listbox, *names):
        listbox.selection_clear(0, tk.END)
        for index, value in enumerate(listbox.get(0, tk.END)):
            if value in names:
                listbox.selection_set(index)
        self.assertCountEqual(self.app.get_listbox_vals(listbox), names)

    def assert_roles(self, covariates, ignored):
        self.assertCountEqual(self.app.get_listbox_vals(self.app.lb_covars), covariates)
        self.assertCountEqual(self.app.get_listbox_vals(self.app.lb_ignore), ignored)

    def switch_sheet(self, sheet):
        self.app.v_sheet.set(sheet)
        self.app.load_columns()

    def configure_comparison(self, mode, time='D14'):
        app = self.app
        if mode == 'GROUP_COMPARISON':
            app.notebook.select(0)
            app.v_group_base.set('NonResp')
            self.select(app.lb_grp_comp, 'Resp')
            self.select(app.lb_grp_times, time)
        elif mode == 'LONGITUDINAL':
            app.notebook.select(1)
            app.v_long_base_time.set('T0')
            self.select(app.lb_long_comp, 'D14')
            self.select(app.lb_long_grps, 'Resp')
        else:
            app.notebook.select(2)
            app.v_delta_base_time.set('T0')
            app.v_delta_base_grp.set('NonResp')
            self.select(app.lb_delta_comp_time, 'D14')
            self.select(app.lb_delta_comp_grp, 'Resp')

    def test_initial_load_restores_configured_roles(self):
        self.assert_roles(['Age', 'Gender'], ['Batch'])
        self.assertEqual(self.app.lb_covars.get(0, tk.END), tuple(self.data.columns))
        self.assertEqual(self.app.lb_ignore.get(0, tk.END), tuple(self.data.columns))

    def test_refresh_and_reordered_sheet_preserve_latest_choices_by_name(self):
        self.select(self.app.lb_covars, 'Gender')
        self.select(self.app.lb_ignore, 'Age', 'Batch')
        self.app.load_columns()
        self.assert_roles(['Gender'], ['Age', 'Batch'])
        self.switch_sheet('Reordered')
        self.assert_roles(['Gender'], ['Age', 'Batch'])
        self.assertEqual(self.app.lb_covars.get(0), 'protein_b')
        self.assertEqual(CONFIG['COVARIATE_COLUMNS'], ['Age', 'Gender'])  # Refresh does not submit.

    def test_explicitly_cleared_selections_stay_cleared(self):
        self.select(self.app.lb_covars)
        self.select(self.app.lb_ignore)
        for sheet in ['s', 'Reordered', 'Missing', 's']:
            self.switch_sheet(sheet)
            self.assert_roles([], [])
        self.configure_comparison('GROUP_COMPARISON')
        self.app.on_submit()
        self.errors.assert_not_called()
        self.assertEqual(CONFIG['COVARIATE_COLUMNS'], [])
        self.assertEqual(CONFIG['UNNEEDED_COLUMNS'], [])

    def test_temporarily_absent_columns_are_remembered_without_reselecting_cleared_columns(self):
        self.switch_sheet('Missing')
        self.assert_roles(['Gender'], [])
        self.select(self.app.lb_covars)  # Clear Gender; Age is currently absent.
        self.switch_sheet('Reordered')
        self.assert_roles(['Age'], ['Batch'])
        self.select(self.app.lb_covars)  # Now also clear Age explicitly.
        self.switch_sheet('s')
        self.assert_roles([], ['Batch'])

    def test_submission_uses_only_columns_present_on_current_sheet(self):
        self.switch_sheet('Missing')
        self.configure_comparison('GROUP_COMPARISON')
        self.app.on_submit()
        self.errors.assert_not_called()
        self.assertEqual(CONFIG['COVARIATE_COLUMNS'], ['Gender'])
        self.assertEqual(CONFIG['UNNEEDED_COLUMNS'], [])

    def test_changing_workbook_preserves_roles_for_matching_column_names(self):
        other = self.folder / 'other.xlsx'
        self.data[self.data.columns[::-1]].to_excel(other, sheet_name='Other', index=False)
        self.select(self.app.lb_covars, 'Gender')
        self.select(self.app.lb_ignore, 'Age', 'Batch')
        with patch.object(gui.filedialog, 'askopenfilename', return_value=str(other)):
            self.app.browse_in()
        self.assertEqual(self.app.v_sheet.get(), 'Other')
        self.assert_roles(['Gender'], ['Age', 'Batch'])

    def test_failed_header_read_leaves_selections_intact(self):
        self.select(self.app.lb_covars, 'Age')
        self.select(self.app.lb_ignore, 'Batch', 'Gender')
        with patch.object(gui.pd, 'read_excel', side_effect=PermissionError('Workbook unavailable')):
            self.app.load_columns()
        self.assert_roles(['Age'], ['Batch', 'Gender'])
        self.app.load_columns()
        self.assert_roles(['Age'], ['Batch', 'Gender'])

    def test_gui_to_analysis_keeps_metadata_out_of_proteins_in_all_four_comparisons(self):
        cases = [('GROUP_COMPARISON', 'D14'), ('GROUP_COMPARISON', 'T0'),
                 ('LONGITUDINAL', 'D14'), ('DELTA', 'D14')]
        for index, (mode, time) in enumerate(cases):
            with self.subTest(mode=mode, time=time):
                if index:
                    self.app = self.new_app()
                self.switch_sheet('Reordered')
                self.configure_comparison(mode, time)
                self.app.v_run_hinge.set(True)
                self.select(self.app.lb_hinge_groups, 'Resp')
                self.select(self.app.lb_hinge_peaks, 'D7')
                self.assert_roles(['Age', 'Gender'], ['Batch'])
                limma, logistic, hinge = Mock(), Mock(), Mock(return_value=None)
                modules = {
                    'limma': SimpleNamespace(run_limma=limma),
                    'logistic_regression': SimpleNamespace(run_logistic_regression=logistic),
                    'hinge_model': SimpleNamespace(run_hinge_model=hinge),
                    'compare_hinge_results': SimpleNamespace(compare_hinge_results=Mock()),
                }
                with patch.dict(sys.modules, modules), patch.object(gui, 'run_gui', side_effect=self.app.on_submit):
                    runpy.run_path(str(Path(__file__).resolve().parents[1] / 'main.py'))
                self.errors.assert_not_called()
                self.assertEqual(CONFIG['ANALYSIS_MODE'], mode)
                self.assertCountEqual(CONFIG['COVARIATE_COLUMNS'], ['Age', 'Gender'])
                self.assertEqual(CONFIG['UNNEEDED_COLUMNS'], ['Batch'])
                for call in [limma.call_args, logistic.call_args]:
                    self.assertCountEqual(call.args[2], ['protein_a', 'protein_b'])
                self.assertCountEqual(hinge.call_args.kwargs['protein_cols'], ['protein_a', 'protein_b'])
                self.assertCountEqual(limma.call_args.args[-1]['COVARIATE_COLUMNS'], ['Age', 'Gender'])
                self.assertTrue({'Age', 'Gender'} <= set(limma.call_args.args[0].columns))
                if mode == 'DELTA':
                    endpoints = logistic.call_args.kwargs['delta_visits']
                    self.assertCountEqual(endpoints['baseline'].columns, ['protein_a', 'protein_b'])
                    self.assertCountEqual(endpoints['followup'].columns, ['protein_a', 'protein_b'])
