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
# 28 Apr 2019: added support for parallel run
###########################################################################


import time
import dask.bag as db
# from multiprocessing import Pool, cpu_count, freeze_support

import os
import sys
# import functools
from pathlib import Path
from typing import Tuple

import pandas
import pathlib


import numpy as np
import pandas as pd

from ._env import path_supy_module
from ._load import (load_InitialCond_grid_df,
                    load_SUEWS_dict_ModConfig,
                    load_SUEWS_Forcing_ESTM_df_raw,
                    load_SUEWS_Forcing_met_df_raw,
                    load_df_state,
                    resample_forcing_met,
                    resample_linear)
from ._post import pack_df_output, pack_df_output_array, pack_df_state
from ._run import (pack_df_state_final,
                   pack_grid_dict,
                   suews_cal_tstep,
                   suews_cal_tstep_multi)
from ._save import get_save_info, save_df_output, save_df_state, save_initcond_nml


##############################################################################
# 1. compact wrapper for loading SUEWS settings
# @functools.lru_cache(maxsize=16)
def init_supy(path_init: str, force_reload=True) -> pd.DataFrame:
    '''Initialise supy by loading initial model states.

    Parameters
    ----------
    path_init : str
        Path to a file that can initialise SuPy, which can be either of the follows:
            * SUEWS :ref:`RunControl.nml<suews:RunControl.nml>`: a namelist file for SUEWS configurations
            * SuPy `df_state.csv`: a CSV file including model states produced by a SuPy run via :py:func:`supy.save_supy`

    force_reload: boolean, optional
        Flag to force reload all initialisation files by clearing all cached states, with default value `True` (i.e., force reload all files).
        Note: If the number of simulation grids is large (e.g., > 100), `force_reload=False` is strongly recommended for better performance.


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
            df_state_init = load_InitialCond_grid_df(
                path_init_x, force_reload=force_reload)
        elif path_init_x.suffix == '.csv':
            # SuPy `df_state.csv`:
            df_state_init = load_df_state(path_init_x)
        else:
            print('{path} is NOT a valid file to initialise SuPy!'.format(
                path=path_init_x))
            sys.exit()
        return df_state_init


# # TODO:
# def load_forcing(path_pattern: str, grid: int = 0) -> pd.DataFrame:
#     pass


# TODO:
# to be superseded by a more generic wrapper: load_forcing
def load_forcing_grid(path_runcontrol: str, grid: int) -> pd.DataFrame:
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
        df_forcing = df_forcing_tstep.round(10)

        # new columns for later use in main calculation
        df_forcing[['iy', 'id', 'it', 'imin']] = df_forcing[[
            'iy', 'id', 'it', 'imin']].astype(np.int64)

    return df_forcing


# load sample data for quickly starting a demo run
# TODO: to deprecate this by renaming for case consistency: load_SampleData-->load_sample_data
def load_SampleData() -> Tuple[pandas.DataFrame, pandas.DataFrame]:
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
        n_yr=10,
) -> Tuple[pandas.DataFrame, pandas.DataFrame]:
    '''Perform supy simulation.

    Parameters
    ----------
    df_forcing : pandas.DataFrame
        forcing data for all grids in `df_state_init`.
    df_state_init : pandas.DataFrame
        initial model states;
        or a collection of model states with multiple timestamps, whose last temporal record will be used as the initial model states.
    save_state : bool, optional
        flag for saving model states at each time step, which can be useful in diagnosing model runtime performance or performing a restart run.
        (the default is False, which instructs supy not to save runtime model states).
    n_yr : int, optional
        chunk size (`n_yr` years) to split simulation periods so memory usage can be reduced.
        (the default is 10, which implies 10-year forcing chunks used in simulations).

    Returns
    -------
    df_output, df_state_final : Tuple[pandas.DataFrame, pandas.DataFrame]
        - df_output: `output results <df_output_var>`
        - df_state_final: `final model states <df_state_var>`

    Examples
    --------

    >>> df_output, df_state_final = supy.run_supy(df_forcing, df_state_init)


    '''


    # set up a timer for simulation time
    start = time.time()

    # save df_init without changing its original data
    # df.copy() in pandas works as a standard python deepcopy
    df_init = df_state_init.copy()

    # print some diagnostic info
    print(f'====================')
    print(f'Simulation period:')
    print(f'  Start: {df_forcing.index[0]}')
    print(f'  End: {df_forcing.index[-1]}')
    print('')

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
            # rename is a workaround to resolve naming inconsistency between
            # suews fortran code interface and input forcing file headers
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
    # reorder columns of df_forcing to comply with SUEWS kernel convention in receiving the input
    # TODO: this re-ordering can be later put into the planned input checker
    list_var_forcing = [
        'iy',
        'id',
        'it',
        'imin',
        'qn1_obs',
        'qh_obs',
        'qe',
        'qs_obs',
        'qf_obs',
        'avu1',
        'avrh',
        'temp_c',
        'press_hpa',
        'precip',
        'avkdn',
        'snow_obs',
        'ldown_obs',
        'fcld_obs',
        'wu_m3',
        'xsmd',
        'lai_obs',
        'kdiff',
        'kdir',
        'wdir',
        'isec',
        'metforcingdata_grid',
        'ts5mindata_ir',
    ]
    df_forcing = df_forcing.loc[:, list_var_forcing]

    # grid list determined by initial states
    list_grid = df_init.index
    print(f'No. of grids: {list_grid.size}\n')

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

        for tstep in df_forcing.index:
            # temporal loop
            # initialise output of tstep:
            # load met_forcing if the same across all grids:
            met_forcing_tstep = dict_forcing[tstep]
            # spatial loop
            for grid in list_grid:
                dict_state_start = dict_state[(tstep, grid)]
                # calculation at one step:
                # series_state_end, series_output_tstep = suews_cal_tstep_df(
                #     series_state_start, met_forcing_tstep)
                try:
                    dict_state_end, dict_output_tstep = suews_cal_tstep(
                        dict_state_start, met_forcing_tstep)
                except:
                    raise RuntimeError('SUEWS kernel error')

                # update output & model state at tstep for the current grid
                dict_output.update({(tstep, grid): dict_output_tstep})
                dict_state.update({(tstep + 1*freq, grid): dict_state_end})

        # pack results as easier DataFrames
        df_output = pack_df_output(dict_output).swaplevel(0, 1)
        # drop unnecessary 'datetime' as it is already included in the index
        df_output = df_output.drop(columns=['datetime'], level=0)
        df_state_final = pack_df_state(dict_state).swaplevel(0, 1)

    else:
        # for multi-year run, reduce the whole df_forcing into {n_yr}-year chunks for less memory consumption
        grp_forcing_yr = df_forcing.groupby(df_forcing.index.year // n_yr)
        if len(grp_forcing_yr) > 1:
            df_state_init_yr = df_state_init.copy()
            list_df_output = []
            list_df_state = []
            for grp in grp_forcing_yr.groups:
                # get forcing of a specific year
                df_forcing_yr = grp_forcing_yr.get_group(grp)
                # run supy: actual execution done in the `else` clause below
                df_output_yr, df_state_final_yr = run_supy(
                    df_forcing_yr, df_state_init_yr)
                df_state_init_yr = df_state_final_yr.copy()
                # collect results
                list_df_output.append(df_output_yr)
                list_df_state.append(df_state_final_yr)

            # re-organise results of each year
            df_output = pd.concat(list_df_output).sort_index()
            df_state_final = pd.concat(
                list_df_state).sort_index().drop_duplicates()
            return df_output, df_state_final

        else:
            # for single-chunk run (1 chunk = {n_yr} years), directly put df_forcing into supy_driver for calculation
            # use higher level wrapper that calculate at a `block` level
            # for better performance

            # # construct input list for `Pool.starmap`
            # construct input list for `dask.bag`
            list_input = [
                # (dict_state[(tstep_init, grid)], df_forcing)
                dict_state[(tstep_init, grid)]
                for grid in list_grid
            ]

            # on windows `processes` has issues when importing
            # so set `threads` here
            method_parallel = 'threads' if os.name == 'nt' else 'processes'
            list_res = db.from_sequence(list_input)\
                .map(suews_cal_tstep_multi, df_forcing)\
                .compute(scheduler=method_parallel)
            try:
                list_state_end, list_output_array = zip(*list_res)
            except:
                raise RuntimeError('SUEWS kernel error')

            # collect output arrays
            dict_output = {
                grid: dict_output_array
                for grid, dict_output_array in zip(list_grid, list_output_array)
            }

            # collect final states
            dict_state_final_tstep = {
                (tstep_final + freq, grid): dict_state_end
                for grid, dict_state_end in zip(list_grid, list_state_end)
            }
            dict_state.update(dict_state_final_tstep)

            # save results as time-aware DataFrame
            df_output0 = pack_df_output_array(dict_output, df_forcing)
            df_output = df_output0.replace(-999., np.nan)
            df_state_final = pack_df_state(dict_state).swaplevel(0, 1)

    # drop ESTM for now as it is not supported yet
    # select only those supported output groups
    df_output = df_output.loc[:, ['SUEWS', 'snow', 'DailyState']]
    # trim multi-index based columns
    df_output.columns = df_output.columns.remove_unused_levels()

    # pack final model states into a proper dataframe
    df_state_final = pack_df_state_final(df_state_final, df_init)

    # show simulation time
    end = time.time()
    print(f'Execution time: {(end - start):.1f} s')
    print(f'====================\n')

    return df_output, df_state_final


##############################################################################
# 3. save results of a supy run
def save_supy(
        df_output: pandas.DataFrame,
        df_state_final: pandas.DataFrame,
        freq_s: int = 3600,
        site: str = '',
        path_dir_save: str = Path('.'),
        path_runcontrol: str = None,) -> list:
    '''Save SuPy run results to files

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
        Path to SUEWS :ref:`RunControl.nml <suews:RunControl.nml>`, which, if set, will be preferably used to derive `freq_s`, `site` and `path_dir_save`.
        (the default is None, which is unset)

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
        freq_s, path_dir_save, site = get_save_info(path_runcontrol)

    # save df_output to several files
    list_path_save = save_df_output(df_output, freq_s, site, path_dir_save)

    # save df_state
    if path_runcontrol is not None:
        # save as nml as SUEWS binary
        list_path_nml = save_initcond_nml(df_state_final, site, path_dir_save)
        list_path_save = list_path_save+list_path_nml
    else:
        # save as supy csv for later use
        path_state_save = save_df_state(df_state_final, site, path_dir_save)
        # update list_path_save
        list_path_save.append(path_state_save)

    return list_path_save
