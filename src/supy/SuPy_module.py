# ###########################################################################
# SuPy: SUEWS for Python
#
# Authors:
# Ting Sun, ting.sun@reading.ac.uk
#
# History:
# 20 Jan 2018: first alpha release
# 01 Feb 2018: performance improvement
# 03 Feb 2018: improvement in output processing
# 08 Mar 2018: pypi packaging
# 04 Oct 2018: overhual of structure
# 05 Oct 2018: added sample run data
###########################################################################


from __future__ import division, print_function

import copy
import functools
from pathlib import Path

import numpy as np
import pandas as pd
from pandas import DataFrame as df

from .supy_env import path_supy_module
from .supy_load import (init_SUEWS_dict, load_SUEWS_Forcing_ESTM_df_raw,
                        load_SUEWS_Forcing_met_df_raw, resample_forcing_met,
                        resample_linear)
from .supy_post import pack_df_output, pack_df_output_array, pack_df_state
from .supy_run import suews_cal_tstep, suews_cal_tstep_multi

##############################################################################
# 1. compact wrapper for loading SUEWS settings


# convert dict_InitCond to pandas Series and DataFrame
# return pd.DataFrame
@functools.lru_cache(maxsize=16)
def init_SUEWS_pd(path_runcontrol):
    dict_mod_cfg, dict_state_init = init_SUEWS_dict(path_runcontrol)
    # ser_mod_cfg: all static model configuration info
    ser_mod_cfg = pd.Series(dict_mod_cfg)
    # df_state_init: initial conditions for SUEWS simulations
    df_state_init = df.from_dict(dict_state_init).T
    df_state_init.index.set_names('grid', inplace=True)

    return ser_mod_cfg, df_state_init


# load forcing datasets of `grid`
@functools.lru_cache(maxsize=16)
def load_SUEWS_Forcing_df_grid(path_runcontrol, grid):
    path_runcontrol = Path(path_runcontrol)
    ser_mod_cfg, df_state_init = init_SUEWS_pd(path_runcontrol)

    # load setting variables from ser_mod_cfg
    (
        filecode,
        kdownzen,
        tstep_met_in,
        tstep_ESTM_in,
        multiplemetfiles,
        multipleestmfiles,
        dir_input_cfg
    ) = ser_mod_cfg[
        [
            'filecode',
            'kdownzen',
            'resolutionfilesin',
            'resolutionfilesinestm',
            'multiplemetfiles',
            'multipleestmfiles',
            'fileinputpath'
        ]
    ]
    tstep_mod, lat, lon, alt, timezone = df_state_init.loc[
        grid,
        ['tstep', 'lat', 'lng', 'alt', 'timezone']]

    path_site = path_runcontrol.parent
    path_input = path_site / ser_mod_cfg['fileinputpath']

    # load raw data
    # met forcing
    df_forcing_met = load_SUEWS_Forcing_met_df_raw(
        path_input, filecode, grid, tstep_met_in, multiplemetfiles)

    # resample raw data from tstep_in to tstep_mod
    df_forcing_met_tstep = resample_forcing_met(
        df_forcing_met, tstep_met_in, tstep_mod,
        lat, lon, alt, timezone, kdownzen)

    # merge forcing datasets (met and ESTM)
    df_forcing_tstep = df_forcing_met_tstep.copy()

    # pack all records of `id` into `metforcingdata_grid` for AnOHM and others
    df_grp = df_forcing_tstep.groupby('id')
    dict_id_all = {xid: df_grp.get_group(xid)
                   for xid in df_forcing_tstep['id'].unique()}
    id_all = df_forcing_tstep['id'].apply(lambda xid: dict_id_all[xid])
    df_forcing_tstep = df_forcing_tstep.merge(
        id_all.to_frame(name='metforcingdata_grid'),
        left_index=True,
        right_index=True)

    # add Ts forcing for ESTM
    if df_state_init.iloc[0]['storageheatmethod'] == 4:
        # load ESTM forcing
        df_forcing_estm = load_SUEWS_Forcing_ESTM_df_raw(
            path_input, filecode, grid, tstep_ESTM_in, multipleestmfiles)
        # resample raw data from tstep_in to tstep_mod
        df_forcing_estm_tstep = resample_linear(
            df_forcing_estm, tstep_met_in, tstep_mod)
        df_forcing_tstep = df_forcing_tstep.merge(
            df_forcing_estm_tstep,
            left_on=['iy', 'id', 'it', 'imin'],
            right_on=['iy', 'id', 'it', 'imin'])
        # insert `ts5mindata_ir` into df_forcing_tstep
        ts_col = df_forcing_estm.columns[4:]
        df_forcing_tstep['ts5mindata_ir'] = (
            df_forcing_tstep.loc[:, ts_col].values.tolist())
        df_forcing_tstep['ts5mindata_ir'] = df_forcing_tstep[
            'ts5mindata_ir'].map(lambda x: np.array(x, order='F'))
    else:
        # insert some placeholder values
        df_forcing_tstep['ts5mindata_ir'] = df_forcing_tstep['temp_c']

    # new columns for later use in main calculation
    df_forcing_tstep[['iy', 'id', 'it', 'imin']] = df_forcing_tstep[[
        'iy', 'id', 'it', 'imin']].astype(np.int64)

    return df_forcing_tstep


# load sample data for quickly starting a demo run
def load_SampleData():
    path_SampleData = Path(path_supy_module) / 'sample_run'
    path_runcontrol = path_SampleData / 'RunControl.nml'
    ser_mod_cfg, df_state_init = init_SUEWS_pd(path_runcontrol)
    # path_input = path_runcontrol.parent / ser_mod_cfg['fileinputpath']
    df_forcing_tstep = load_SUEWS_Forcing_df_grid(
        path_runcontrol,
        df_state_init.index[0]
    )
    return ser_mod_cfg, df_state_init, df_forcing_tstep


# input processing code end here
##############################################################################


##############################################################################
# 2. compact wrapper for running a whole simulation
# # main calculation
# input as DataFrame
def run_suews_df(df_forcing, df_init_input, save_state=False):
    # save df_init without changing its original data
    # df.copy() in pandas does work as a standard python deepcopy
    df_init = pd.DataFrame(copy.deepcopy(df_init_input.to_dict()))
    # grid list determined by initial states
    grid_list = df_init.index

    # initialise dicts for holding results and model states
    dict_state = {}
    dict_output = {}

    if save_state:
        # use slower more functional single step wrapper
        # start tstep retrived from forcing data
        t_start = df_forcing.index[0]
        # convert df to dict with `itertuples` for better performance
        dict_forcing = {row.Index: row._asdict()
                        for row in df_forcing.itertuples()}

        # dict_state is used to save model states for later use
        dict_state = {(t_start, grid): series_state_init.to_dict()
                      for grid, series_state_init
                      in df_init.iterrows()}
        for tstep in df_forcing.index:
            # temporal loop
            # initialise output of tstep:
            # load met_forcing if the same across all grids:
            met_forcing_tstep = dict_forcing[tstep]
            # spatial loop
            for grid in grid_list:
                dict_state_start = dict_state[(tstep, grid)]
                # calculation at one step:
                # series_state_end, series_output_tstep = suews_cal_tstep_df(
                #     series_state_start, met_forcing_tstep)
                dict_state_end, dict_output_tstep = suews_cal_tstep(
                    dict_state_start, met_forcing_tstep,
                    save_state=save_state)
                # update output & model state at tstep for the current grid
                dict_output.update({(tstep, grid): dict_output_tstep})
                dict_state.update({(tstep + 1, grid): dict_state_end})

        # pack results as easier DataFrames
        df_output = pack_df_output(dict_output)
        df_state = pack_df_state(dict_state)

    else:
        # use higher level wrapper that calculate at a `block` level
        # for better performance
        dict_state = {grid: df_init.loc[grid]
                      for grid in grid_list}
        for grid in grid_list:
            df_state_start_grid = dict_state[grid]
            dict_state_end, dict_output_array = suews_cal_tstep_multi(
                df_state_start_grid, df_forcing)
            # update output & model state at tstep for the current grid
            dict_output.update({grid: dict_output_array})
            dict_state.update({grid: dict_state_end})

        # save results as time-aware DataFrame
        df_output0 = pack_df_output_array(dict_output, df_forcing)
        df_output = df_output0.replace(-999., np.nan)
        df_state = pd.DataFrame(dict_state).T

    return df_output, df_state
