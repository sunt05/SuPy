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
# 04 Oct 2018: overhaul of structure
# 05 Oct 2018: added sample run data
###########################################################################


from __future__ import division, print_function

import os
import sys
# import functools
from pathlib import Path
from typing import Tuple

import pandas
import pathlib


import numpy as np
import pandas as pd

from .supy_env import path_supy_module
from .supy_load import (load_InitialCond_grid_df, load_SUEWS_dict_ModConfig,
                        load_SUEWS_Forcing_ESTM_df_raw,
                        load_SUEWS_Forcing_met_df_raw,
                        load_df_state,
                        resample_forcing_met,
                        resample_linear)
from .supy_post import pack_df_output, pack_df_output_array, pack_df_state
from .supy_run import (pack_df_state_final, pack_grid_dict, suews_cal_tstep,
                       suews_cal_tstep_multi)
from .supy_save import get_save_info, save_df_output, save_df_state


##############################################################################
# 1. compact wrapper for loading SUEWS settings
# @functools.lru_cache(maxsize=16)
def init_supy(path_init: str)->pd.DataFrame:
    '''Initialise supy by loading initial model states.

    Parameters
    ----------
    path_init : str
        Path to a file that can initialise SuPy, which can be either of the follows:
            * SUEWS :ref:`RunControl.nml<suews:RunControl.nml>`: a namelist file for SUEWS configurations
            * SuPy `df_state.csv`: a CSV file including model states produced by a SuPy run via :py:func:`supy.save_supy`

    Returns
    -------
    df_state_init: pandas.DataFrame
        Initial model states.
        See `df_state_var` for details.

    Examples
    --------
    1. Use :ref:`RunControl.nml<suews:RunControl.nml>` to initialise SuPy

    >>> path_init = "~/SUEWS_sims/RunControl.nml"
    >>> df_state_init = supy.init_supy(path_init)

    2. Use ``df_state.csv`` to initialise SuPy

    >>> path_init = "~/SuPy_res/df_state_test.csv"
    >>> df_state_init = supy.init_supy(path_init)

    '''

    try:
        path_init_x = Path(path_init).expanduser().resolve()
    except FileNotFoundError:
        print('{path} does not exists!'.format(path=path_init_x))
    else:
        if path_init_x.suffix == '.nml':
            # SUEWS `RunControl.nml`:
            df_state_init = load_InitialCond_grid_df(path_init_x)
        elif path_init_x.suffix == '.csv':
            # SuPy `df_state.csv`:
            df_state_init = load_df_state(path_init_x)
        else:
            print('{path} is NOT a valid file to initialise SuPy!'.format(
                path=path_init_x))
            sys.exit()
        return df_state_init


def load_forcing_grid(path_runcontrol: str, grid: int)->pd.DataFrame:
    '''Load forcing data for a specific grid included in the index of `df_state_init </data-structure/supy-io.ipynb#df_state_init:-model-initial-states>`.

    Parameters
    ----------
    path_runcontrol : str
        Path to SUEWS :ref:`RunControl.nml <suews:RunControl.nml>`
    grid : int
        Grid number

    Returns
    -------
    df_forcing: pandas.DataFrame
        Forcing data. See `df_forcing_var` for details.

    Examples
    --------
    >>> path_runcontrol = "~/SUEWS_sims/RunControl.nml"  # a valid path to `RunControl.nml`
    >>> df_state_init = supy.init_supy(path_runcontrol) # get `df_state_init`
    >>> grid = df_state_init.index[0] # first grid number included in `df_state_init`
    >>> df_forcing = supy.load_forcing_grid(path_runcontrol, grid) # get df_forcing


    '''

    try:
        path_runcontrol = Path(path_runcontrol).expanduser().resolve()
    except FileNotFoundError:
        print('{path} does not exists!'.format(path=path_runcontrol))
    else:
        dict_mod_cfg = load_SUEWS_dict_ModConfig(path_runcontrol)
        df_state_init = init_supy(path_runcontrol)

        # load setting variables from dict_mod_cfg
        (
            filecode,
            kdownzen,
            tstep_met_in,
            tstep_ESTM_in,
            multiplemetfiles,
            multipleestmfiles,
            dir_input_cfg
        ) = (dict_mod_cfg[x] for x in
             [
            'filecode',
            'kdownzen',
            'resolutionfilesin',
            'resolutionfilesinestm',
            'multiplemetfiles',
            'multipleestmfiles',
            'fileinputpath'
        ]
        )
        tstep_mod, lat, lon, alt, timezone = df_state_init.loc[
            grid,
            [(x, '0') for x in ['tstep', 'lat', 'lng', 'alt', 'timezone']]
        ].values

        path_site = path_runcontrol.parent
        path_input = path_site / dict_mod_cfg['fileinputpath']

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

        # disable the AnOHM and ESTM components for now and for better performance
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # TS 28 Dec 2018
        # pack all records of `id` into `metforcingdata_grid` for AnOHM
        # df_grp = df_forcing_tstep.groupby('id')
        # dict_id_all = {xid: df_grp.get_group(xid)
        #                for xid in df_forcing_tstep['id'].unique()}
        # id_all = df_forcing_tstep['id'].apply(lambda xid: dict_id_all[xid])
        # df_forcing_tstep = df_forcing_tstep.merge(
        #     id_all.to_frame(name='metforcingdata_grid'),
        #     left_index=True,
        #     right_index=True)
        # # add Ts forcing for ESTM
        # if np.asscalar(df_state_init.iloc[0]['storageheatmethod'].values) == 4:
        #     # load ESTM forcing
        #     df_forcing_estm = load_SUEWS_Forcing_ESTM_df_raw(
        #         path_input, filecode, grid, tstep_ESTM_in, multipleestmfiles)
        #     # resample raw data from tstep_in to tstep_mod
        #     df_forcing_estm_tstep = resample_linear(
        #         df_forcing_estm, tstep_met_in, tstep_mod)
        #     df_forcing_tstep = df_forcing_tstep.merge(
        #         df_forcing_estm_tstep,
        #         left_on=['iy', 'id', 'it', 'imin'],
        #         right_on=['iy', 'id', 'it', 'imin'])
        #     # insert `ts5mindata_ir` into df_forcing_tstep
        #     ts_col = df_forcing_estm.columns[4:]
        #     df_forcing_tstep['ts5mindata_ir'] = (
        #         df_forcing_tstep.loc[:, ts_col].values.tolist())
        #     df_forcing_tstep['ts5mindata_ir'] = df_forcing_tstep[
        #         'ts5mindata_ir'].map(lambda x: np.array(x, order='F'))
        # else:
        #     # insert some placeholder values
        #     df_forcing_tstep['ts5mindata_ir'] = df_forcing_tstep['Tair']
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # disable the AnOHM and ESTM components for now and for better performance

        # coerced precision here to prevent numerical errors inside Fortran
        df_forcing = np.around(df_forcing_tstep, decimals=10)
        # new columns for later use in main calculation
        df_forcing[['iy', 'id', 'it', 'imin']] = df_forcing[[
            'iy', 'id', 'it', 'imin']].astype(np.int64)

    return df_forcing


# load sample data for quickly starting a demo run
def load_SampleData()->Tuple[pandas.DataFrame, pandas.DataFrame]:
    '''Load sample data for quickly starting a demo run.

    Returns
    -------
    df_state_init, df_forcing: Tuple[pandas.DataFrame, pandas.DataFrame]
        - df_state_init: `initial model states <df_state_var>`
        - df_forcing: `forcing data <df_forcing_var>`

    Examples
    --------

    >>> df_state_init, df_forcing = supy.load_SampleData()

    '''

    path_SampleData = Path(path_supy_module) / 'sample_run'
    path_runcontrol = path_SampleData / 'RunControl.nml'
    df_state_init = init_supy(path_runcontrol)
    # path_input = path_runcontrol.parent / ser_mod_cfg['fileinputpath']
    df_forcing = load_forcing_grid(
        path_runcontrol,
        df_state_init.index[0]
    )
    return df_state_init, df_forcing

# input processing code end here
##############################################################################


##############################################################################
# 2. compact wrapper for running a whole simulation
# # main calculation
# input as DataFrame

def run_supy(
        df_forcing: pandas.DataFrame,
        df_state_init: pandas.DataFrame,
        save_state=False,
)->Tuple[pandas.DataFrame, pandas.DataFrame]:
    '''Perform supy simulaiton.

    Parameters
    ----------
    df_forcing : pandas.DataFrame
        forcing data.
    df_state_init : pandas.DataFrame
        initial model states;
        or a collection of model states with multiple timestamps, whose last temporal record will be used as the initial model states.
    save_state : bool, optional
        flag for saving model states at each timestep, which can be useful in diagnosing model runtime performance or performing a restart run.
        (the default is False, which intructs supy not to save runtime model states).

    Returns
    -------
    df_output, df_state_final : Tuple[pandas.DataFrame, pandas.DataFrame]
        - df_output: `output results <df_output_var>`
        - df_state_final: `final model states <df_state_var>`

    Examples
    --------

    >>> df_output, df_state_final = supy.run_supy(df_forcing, df_state_init)


    '''

    # save df_init without changing its original data
    # df.copy() in pandas does work as a standard python deepcopy
    df_init = df_state_init.copy()

    # retrieve the last temporal record as `df_init`
    # if a `datetime` level existing in the index
    if df_init.index.nlevels > 1:
        idx_dt = df_init.index.get_level_values('datetime').unique()
        dt_last = idx_dt.max()
        df_init = df_init.loc[dt_last]
    # add placeholder variables for df_forcing
    # `metforcingdata_grid` and `ts5mindata_ir` are used by AnOHM and ESTM, respectively
    # they are now temporarily disabled in supy
    df_forcing = df_forcing\
        .assign(
            metforcingdata_grid=0,
            ts5mindata_ir=0,
        )\
        .rename(
            # remanae is a workaround to resolve naming inconsistency between
            # suews fortran code interface and input forcing file hearders
            columns={
                '%' + 'iy': 'iy',
                'id': 'id',
                'it': 'it',
                'imin': 'imin',
                'qn': 'qn1_obs',
                'qh': 'qh_obs',
                'qe': 'qe',
                'qs': 'qs_obs',
                'qf': 'qf_obs',
                'U': 'avu1',
                'RH': 'avrh',
                'Tair': 'temp_c',
                'pres': 'press_hpa',
                'rain': 'precip',
                'kdown': 'avkdn',
                'snow': 'snow_obs',
                'ldown': 'ldown_obs',
                'fcld': 'fcld_obs',
                'Wuh': 'wu_m3',
                'xsmd': 'xsmd',
                'lai': 'lai_obs',
                'kdiff': 'kdiff',
                'kdir': 'kdir',
                'wdir': 'wdir',
            }
        )
    # grid list determined by initial states
    grid_list = df_init.index

    # initialise dicts for holding results and model states
    dict_state = {}
    dict_output = {}

    # initial and final tsteps retrieved from forcing data
    tstep_init = df_forcing.index[0]
    tstep_final = df_forcing.index[-1]
    # tstep size retrieved from forcing data
    freq = df_forcing.index.freq

    # dict_state is used to save model states for later use
    dict_state = {
        # (t_start, grid): series_state_init.to_dict()
        (tstep_init, grid): pack_grid_dict(series_state_init)
        for grid, series_state_init
        in df_init.iterrows()
    }

    # remove 'problems.txt'
    if Path('problems.txt').exists():
        os.remove('problems.txt')

    if save_state:
        # use slower more functional single step wrapper

        # convert df to dict with `itertuples` for better performance
        dict_forcing = {row.Index: row._asdict()
                        for row in df_forcing.itertuples()}

        # # dict_state is used to save model states for later use
        # dict_state = {
        #     # (t_start, grid): series_state_init.to_dict()
        #     (tstep_init, grid): pack_grid_dict(series_state_init)
        #     for grid, series_state_init
        #     in df_init.iterrows()
        # }

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
                    dict_state_start, met_forcing_tstep)

                # update output & model state at tstep for the current grid
                dict_output.update({(tstep, grid): dict_output_tstep})
                dict_state.update({(tstep + 1*freq, grid): dict_state_end})

        # pack results as easier DataFrames
        df_output = pack_df_output(dict_output).swaplevel(0, 1)
        # drop unnecessary 'datetime' as it is already included in the index
        df_output = df_output.drop(columns=['datetime'], level=0)
        df_state_final = pack_df_state(dict_state).swaplevel(0, 1)

    else:
        # use higher level wrapper that calculate at a `block` level
        # for better performance

        # dict_state = {
        #     # grid: df_init.loc[grid]
        #     (tstep_init, grid): pack_grid_dict(series_state_init)
        #     for grid, series_state_init
        #     in df_init.iterrows()
        # }

        for grid in grid_list:
            dict_state_start_grid = dict_state[(tstep_init, grid)]
            dict_state_end, dict_output_array = suews_cal_tstep_multi(
                dict_state_start_grid, df_forcing)
            # update output & model state at tstep for the current grid
            dict_output.update({grid: dict_output_array})
            # model state for the next run
            dict_state.update({(tstep_final + freq, grid): dict_state_end})

        # save results as time-aware DataFrame
        df_output0 = pack_df_output_array(dict_output, df_forcing)
        df_output = df_output0.replace(-999., np.nan)
        df_state_final = pack_df_state(dict_state).swaplevel(0, 1)
        # df_state = pd.DataFrame(dict_state).T
        # df_state.index.set_names('grid')

    # drop ESTM for now as it is not supported yet
    # select only those supported output groups
    df_output = df_output.loc[:, ['SUEWS', 'snow', 'DailyState']]
    # trim multiindex based columns
    df_output.columns = df_output.columns.remove_unused_levels()

    # pack final model states into a proper dataframe
    df_state_final = pack_df_state_final(df_state_final, df_init)

    return df_output, df_state_final


##############################################################################
# 3. save results of a supy run
def save_supy(
        df_output: pandas.DataFrame,
        df_state_final: pandas.DataFrame,
        freq_s: int = 3600,
        site: str = '',
        path_dir_save: str = Path('.'),
        path_runcontrol: str = None,)->list:
    '''save supy run results to files

    Parameters
    ----------
    df_output : pandas.DataFrame
        DataFrame of output
    df_state_final : pandas.DataFrame
        DataFrame of final model states
    freq_s : int, optional
        Output frequency in seconds (the default is 3600, which indicates hourly output)
    site : str, optional
        Site identifier (the default is '', which indicates site identifier will be left empty)
    path_dir_save : str, optional
        Path to directory to saving the files (the default is Path('.'), which indicates the current working directory)
    path_runcontrol : str, optional
        Path to SUEWS :ref:`RunControl.nml <suews:RunControl.nml>` (the default is None, which is unset)

    Returns
    -------
    list
        a list of paths of saved files

    Examples
    --------
    1. save results of a supy run to the current working directory with default settings

    >>> list_path_save = supy.save_supy(df_output, df_state_final)


    2. save results according to settings in :ref:`RunControl.nml <suews:RunControl.nml>`

    >>> list_path_save = supy.save_supy(df_output, df_state_final, path_runcontrol='path/to/RunControl.nml')


    3. save results of a supy run at resampling frequency of 1800 s (i.e., half-hourly results) under the site code ``Test`` to a customised location 'path/to/some/dir'

    >>> list_path_save = supy.save_supy(df_output, df_state_final, freq_s=1800, site='Test', path_dir_save='path/to/some/dir')
    '''

    # get necessary information for saving procedure
    if path_runcontrol is not None:
        freq_s, dir_save, site = get_save_info(path_runcontrol)

    # save df_output to several files
    list_path_save = save_df_output(df_output, freq_s, site, path_dir_save)

    # save df_state
    path_state_save = save_df_state(df_state_final, site, path_dir_save)

    # update list_path_save
    list_path_save.append(path_state_save)
    return list_path_save
