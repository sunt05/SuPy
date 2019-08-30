
import io
import sys
import warnings
from time import time
from unittest import TestCase

import numpy as np
import pandas as pd

import supy as sp
import platform


class TestSuPy(TestCase):
    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)

    # test if supy_driver can be connected
    def test_is_driver_connected(self):
        s = sp._run.list_var_output
        self.assertTrue(isinstance(s[0], np.str_))

    # test if single-tstep mode can run
    def test_is_supy_running_single_step(self):
        df_state_init, df_forcing_tstep = sp.load_SampleData()
        df_forcing_part = df_forcing_tstep.iloc[:288 * 1]
        df_output, df_state = sp.run_supy(
            df_forcing_part, df_state_init,
            save_state=True)
        test_non_empty = np.all(
            [
                not df_output.empty,
                not df_state.empty,
            ]
        )
        self.assertTrue(test_non_empty)

    # test if multi-tstep mode can run
    def test_is_supy_running_multi_step(self):
        df_state_init, df_forcing_tstep = sp.load_SampleData()
        df_forcing_part = df_forcing_tstep.iloc[:]
        t_start = time()
        df_output, df_state = sp.run_supy(
            df_forcing_part, df_state_init)
        t_end = time()

        # only print to screen on macOS due incompatibility on Windows
        if platform.system() == 'Darwin':
            capturedOutput = io.StringIO()  # Create StringIO object
            sys.stdout = capturedOutput  # and redirect stdout.
            # Call function.
            print(f'Running time: {t_end-t_start:.2f} s')
            sys.stdout = sys.__stdout__                     # Reset redirect.
            # Now works as before.
            print('Captured:\n', capturedOutput.getvalue())

        test_non_empty = np.all(
            [
                not df_output.empty,
                not df_state.empty,
            ]
        )
        self.assertTrue(test_non_empty)

    # test if single-tstep and multi-tstep modes can produce the same SUEWS results
    def test_is_supy_euqal_mode(self):
        df_state_init, df_forcing_tstep = sp.load_SampleData()
        df_forcing_part = df_forcing_tstep.iloc[:288 * 1]
        # single-step results
        df_output_s, df_state_s = sp.run_supy(
            df_forcing_part, df_state_init,
            save_state=True)
        df_res_s = df_output_s\
            .loc[:, [
                'SUEWS',
                'DailyState',
                'snow',
            ]]\
            .fillna(-999.)\
            .sort_index(axis=1)\
            .round(6)\
            .applymap(lambda x: -999. if np.abs(x) > 3e4 else x)

        df_state_init, df_forcing_tstep = sp.load_SampleData()
        # multi-step results
        df_output_m, df_state_m = sp.run_supy(
            df_forcing_part, df_state_init,
            save_state=False)
        df_res_m = df_output_m\
            .loc[:, [
                'SUEWS',
                 'DailyState',
                 'snow',
            ]]\
            .fillna(-999.)\
            .sort_index(axis=1)\
            .round(6)\
            .applymap(lambda x: -999. if np.abs(x) > 3e4 else x)
        # print(df_res_m.iloc[:3, 86], df_res_s.iloc[:3, 86])
        pd.testing.assert_frame_equal(
            left=df_res_s,
            right=df_res_m,
        )
        # test_equal_mode = df_res_s.eq(df_res_m).all(None)
        # self.assertTrue(test_equal_mode)
