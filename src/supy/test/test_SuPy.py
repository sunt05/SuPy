import tempfile
from pathlib import Path
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
        warnings.simplefilter("ignore", category=ImportWarning)

    # test if supy_driver can be connected
    def test_is_driver_connected(self):
        s = sp._run.list_var_output
        self.assertTrue(isinstance(s[0], np.str_))

    # test if single-tstep mode can run
    def test_is_supy_running_single_step(self):
        df_state_init, df_forcing_tstep = sp.load_SampleData()
        df_forcing_part = df_forcing_tstep.iloc[: 288 * 1]
        df_output, df_state = sp.run_supy(
            df_forcing_part, df_state_init, save_state=True
        )
        test_non_empty = np.all([not df_output.empty, not df_state.empty,])
        self.assertTrue(test_non_empty)

    # test if multi-tstep mode can run
    def test_is_supy_running_multi_step(self):
        df_state_init, df_forcing_tstep = sp.load_SampleData()
        df_forcing_part = df_forcing_tstep.iloc[:]
        df_output, df_state = sp.run_supy(
            df_forcing_part, df_state_init, check_input=True
        )

        # # only print to screen on macOS due incompatibility on Windows
        # if platform.system() == "Darwin":
        #     # capturedOutput = io.StringIO()  # Create StringIO object
        #     # sys.stdout = capturedOutput  # and redirect stdout.
        #     # Call function.
        #     print(f"Running time: {t_end-t_start:.2f} s")
        #     # sys.stdout = sys.__stdout__  # Reset redirect.
        #     # Now works as before.
        #     # print("Captured:\n", capturedOutput.getvalue())

        test_non_empty = np.all([not df_output.empty, not df_state.empty,])
        self.assertTrue(test_non_empty)

    # test if multi-grid simulation can run in parallel
    def test_is_supy_sim_save_multi_grid_par(self):
        n_grid = 4
        df_state_init, df_forcing_tstep = sp.load_SampleData()
        df_state_init = pd.concat([df_state_init for x in range(n_grid)])
        df_state_init.index = pd.RangeIndex(n_grid, name="grid")
        df_forcing_part = df_forcing_tstep.iloc[:288*60]
        df_output, df_state = sp.run_supy(df_forcing_part, df_state_init)

        test_success_sim = np.all([not df_output.empty, not df_state.empty,])

        with tempfile.TemporaryDirectory() as dir_temp:
            list_outfile = sp.save_supy(
                df_output,
                df_state,
                path_dir_save=dir_temp,
                site="pytest",
                logging_level=10,
            )

        test_success_save = np.all([isinstance(fn, Path) for fn in list_outfile])
        self.assertTrue(test_success_sim and test_success_save)

        # # only print to screen on macOS due incompatibility on Windows
        # if platform.system() == "Darwin":
        #     # capturedOutput = io.StringIO()  # Create StringIO object
        #     # sys.stdout = capturedOutput  # and redirect stdout.
        #     # Call function.
        #     n_grid = df_state_init.index.size
        #     print(f"Running time: {t_end-t_start:.2f} s for {n_grid} grids")
        #     # sys.stdout = sys.__stdout__  # Reset redirect.
        #     # Now works as before.
        #     # print("Captured:\n", capturedOutput.getvalue())

        # test_non_empty = np.all([not df_output.empty, not df_state.empty,])
        # self.assertTrue(test_non_empty)

    # test if single-tstep and multi-tstep modes can produce the same SUEWS results
    def test_is_supy_euqal_mode(self):
        df_state_init, df_forcing_tstep = sp.load_SampleData()
        df_forcing_part = df_forcing_tstep.iloc[: 288 * 1]
        # single-step results
        df_output_s, df_state_s = sp.run_supy(
            df_forcing_part, df_state_init, save_state=True
        )
        df_res_s = (
            df_output_s.loc[:, ["SUEWS", "DailyState", "snow",]]
            .fillna(-999.0)
            .sort_index(axis=1)
            .round(6)
            .applymap(lambda x: -999.0 if np.abs(x) > 3e4 else x)
        )

        df_state_init, df_forcing_tstep = sp.load_SampleData()
        # multi-step results
        df_output_m, df_state_m = sp.run_supy(
            df_forcing_part, df_state_init, save_state=False
        )
        df_res_m = (
            df_output_m.loc[:, ["SUEWS", "DailyState", "snow",]]
            .fillna(-999.0)
            .sort_index(axis=1)
            .round(6)
            .applymap(lambda x: -999.0 if np.abs(x) > 3e4 else x)
        )
        # print(df_res_m.iloc[:3, 86], df_res_s.iloc[:3, 86])
        pd.testing.assert_frame_equal(
            left=df_res_s, right=df_res_m,
        )
        # test_equal_mode = df_res_s.eq(df_res_m).all(None)
        # self.assertTrue(test_equal_mode)

    # # test saving output files working
    # def test_is_supy_save_working(self):
    #     df_state_init, df_forcing_tstep = sp.load_SampleData()
    #     # df_state_init = pd.concat([df_state_init for x in range(6)])
    #     df_forcing_part = df_forcing_tstep.iloc[: 288 * 2]
    #     t_start = time()
    #     df_output, df_state = sp.run_supy(df_forcing_part, df_state_init)
    #     t_end = time()
    #     with tempfile.TemporaryDirectory() as dir_temp:
    #         list_outfile = sp.save_supy(df_output, df_state, path_dir_save=dir_temp)

    #     # only print to screen on macOS due incompatibility on Windows
    #     if platform.system() == "Darwin":
    #         capturedOutput = io.StringIO()  # Create StringIO object
    #         sys.stdout = capturedOutput  # and redirect stdout.
    #         # Call function.
    #         n_grid = df_state_init.index.size
    #         print(f"Running time: {t_end-t_start:.2f} s for {n_grid} grids")
    #         sys.stdout = sys.__stdout__  # Reset redirect.
    #         # Now works as before.
    #         print("Captured:\n", capturedOutput.getvalue())

    #     test_non_empty = np.all([isinstance(fn, Path) for fn in list_outfile])
    #     self.assertTrue(test_non_empty)

    # test saving output files working
    def test_is_checking_complete(self):
        df_state_init, df_forcing_tstep = sp.load_SampleData()
        dict_rules = sp._check.dict_rules

        # variables in loaded dataframe
        set_var_df_init = set(df_state_init.columns.get_level_values("var"))

        # variables in dict_rules
        set_var_dict_rules = set(list(dict_rules.keys()))

        # common variables
        set_var_common = set_var_df_init.intersection(set_var_dict_rules)

        # test if common variables are all those in `df_state_init`
        test_common_all = set_var_df_init == set_var_common
        self.assertTrue(test_common_all)

    # # test ERA5 forcing generation
    # def test_gen_forcing(self):
    #     import xarray as xr
    #     # mimic downloading
    #     dict_era5_file = sp.util.download_era5(
    #         57.7081, 11.9653, "20030101", "20031231", dir_save="./data_test"
    #     )
    #     list_fn_ml = [k for k in dict_era5_file.keys() if "ml" in k]
    #     list_fn_sfc = [k for k in dict_era5_file.keys() if "sfc" in k]
    #     # test forcing generation
    #     list_fn_fc = sp.util.gen_forcing_era5(
    #         57.7081, 11.9653, "20030101", "20031231", dir_save="./data_test"
    #     )
    #     df_forcing = sp.util.read_suews(list_fn_fc[0])
    #     ds_sfc = xr.open_mfdataset(list_fn_sfc)
    #     ser_t2 = ds_sfc.t2m.to_series()
    #     res_dif=((df_forcing.Tair + 273.15 - ser_t2.values) / 98).round(4)
    #     test_dif= -0.0066<res_dif.max()<-0.0063
    #     self.assertTrue(test_dif)
