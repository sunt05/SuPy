from unittest import TestCase
import warnings

import supy as sp
import numpy as np
# from pathlib import Path
# from pathlib import Path


class TestSuPy(TestCase):
    def setUp(self):
        warnings.simplefilter('ignore', category=ImportWarning)

    def test_is_driver_connected(self):
        s = sp.supy_run.list_var_output
        self.assertTrue(isinstance(s[0], np.str_))

    # def test_is_sampledata_init_loaded(self):
    #     path_SampleData = Path(sp.path_supy_module) / 'sample_run'
    #     ser_mod_cfg, df_state_init = sp.init_SUEWS_pd(path_SampleData)
    #     test_non_empty = np.all(
    #         [
    #             not ser_mod_cfg.empty,
    #             not df_state_init.empty,
    #         ]
    #     )
    #     self.assertTrue(test_non_empty)

    # def test_is_sampledata_forcing_loaded(self):
    #     ser_mod_cfg, df_state_init, df_forcing_tstep = sp.load_SampleData()
    #     test_non_empty = np.all(
    #         [
    #             not ser_mod_cfg.empty,
    #             not df_state_init.empty,
    #             not df_forcing_tstep.empty,
    #         ]
    #     )
    #     self.assertTrue(test_non_empty)

    def test_is_supy_running(self):
        df_state_init, df_forcing_tstep = sp.load_SampleData()
        df_forcing_part = df_forcing_tstep.iloc[:288 * 1]
        df_output, df_state = sp.run_suews_df(
            df_forcing_part, df_state_init,
            save_state=True)
        test_non_empty = np.all(
            [
                not df_output.empty,
                not df_state.empty,
            ]
        )
        self.assertTrue(test_non_empty)
