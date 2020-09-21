from shutil import rmtree
import tempfile
import copy
import multiprocess
import os
import sys
import time

# import logging
import traceback
from ast import literal_eval
from pathlib import Path
from typing import Tuple
import pandas

import numpy as np
import pandas as pd
from supy_driver import suews_driver as sd

from ._load import (
    df_var_info,
    list_var_inout,
    list_var_inout_multitsteps,
    list_var_input,
    list_var_input_multitsteps,
    list_var_output,
    list_var_output_multitsteps,
)
from ._post import pack_df_output, pack_df_output_array, pack_df_state

from ._env import logger_supy


##############################################################################
# main calculation
# 1. calculation code for one time step
# 2. compact wrapper for running a whole simulation


# 1. calculation code for one time step


# test for performance
# dict_var_inout = {k: None for k in set_var_inout}


# high-level wrapper: suews_cal_tstep
def suews_cal_tstep(dict_state_start, dict_met_forcing_tstep):
    # save_state=False):
    # use single dict as input for suews_cal_main
    dict_input = copy.deepcopy(dict_state_start)
    dict_input.update(dict_met_forcing_tstep)
    dict_input = {k: dict_input[k] for k in list_var_input}

    # for var in dict_input:
    #     print(var)
    #     print(dict_input[var])
    #     print('\n')

    # main calculation:
    try:
        res_suews_tstep = sd.suews_cal_main(**dict_input)
    except Exception as ex:
        # show trace info
        logger_supy.exception(traceback.format_exc())
        # show SUEWS fatal error details produced by SUEWS kernel
        with open("problems.txt", "r") as f:
            logger_supy.critical(f.read())
        # clean slate
        # os.remove('problems.txt')
        # sys.exit()
        logger_supy.critical("SUEWS kernel error")
    else:
        # update state variables
        # if save_state:  # deep copy states results
        dict_state_end = copy.deepcopy(dict_state_start)
        dict_state_end.update(
            {var: copy.copy(dict_input[var]) for var in list_var_inout}
        )

        # update timestep info
        dict_state_end["tstep_prev"] = dict_state_end["tstep"]
        dict_state_end["dt_since_start"] += dict_state_end["tstep"]

        # pack output
        dict_output = {k: v for k, v in zip(list_var_output, res_suews_tstep)}

        return dict_state_end, dict_output


# high-level wrapper: suews_cal_tstep
# def suews_cal_tstep_multi(df_state_start_grid, df_met_forcing_block):
def suews_cal_tstep_multi(dict_state_start_grid, df_met_forcing_block):
    # use single dict as input for suews_cal_main
    dict_input = copy.deepcopy(dict_state_start_grid)
    dict_input.update(
        {
            "metforcingblock": np.array(
                df_met_forcing_block.drop(
                    columns=["metforcingdata_grid", "ts5mindata_ir", "isec",]
                ),
                order="F",
            ),
            "ts5mindata_ir": np.array(df_met_forcing_block["ts5mindata_ir"], order="F"),
            "len_sim": np.array(df_met_forcing_block.shape[0], dtype=int),
        }
    )
    dict_input = {k: dict_input[k] for k in list_var_input_multitsteps}

    # for var in dict_input:
    #     print(var)
    #     print(dict_input[var])
    #     print(dict_input[var].shape, '\n')

    # main calculation:
    try:
        res_suews_tstep_multi = sd.suews_cal_multitsteps(**dict_input)
    except Exception as ex:
        # show trace info
        # print(traceback.format_exc())
        # show SUEWS fatal error details produced by SUEWS kernel
        with open("problems.txt", "r") as f:
            logger_supy.critical(f.read())
        # clean slate
        # os.remove('problems.txt')
        # sys.exit()
        # raise RuntimeError("Something bad happened") from exs
        logger_supy.critical("SUEWS kernel error")
    else:
        # update state variables
        # dict_state_end = copy.copy(dict_input)
        dict_state_end = copy.deepcopy(dict_state_start_grid)
        dict_state_end.update(
            {var: dict_input[var] for var in list_var_inout_multitsteps}
        )

        # update timestep info
        dict_state_end["tstep_prev"] = dict_state_end["tstep"]
        idx_dt = df_met_forcing_block.index
        duration_s = int((idx_dt[-1] - idx_dt[0]).total_seconds())
        dict_state_end["dt_since_start"] += duration_s + dict_state_end["tstep"]

        # pack output
        dict_output_array = {
            k: v for k, v in zip(list_var_output[1:], res_suews_tstep_multi)
        }

        return dict_state_end, dict_output_array


# dataframe based wrapper
# serial mode:
def run_supy_ser(
    df_forcing: pandas.DataFrame,
    df_state_init: pandas.DataFrame,
    save_state=False,
    chunk_day=3660,
) -> Tuple[pandas.DataFrame, pandas.DataFrame]:
    """Perform supy simulation.

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
    chunk_day : int, optional
        chunk size (`chunk_day` days) to split simulation periods so memory usage can be reduced.
        (the default is 3660, which implies ~10-year forcing chunks used in simulations).

    Returns
    -------
    df_output, df_state_final : Tuple[pandas.DataFrame, pandas.DataFrame]
        - df_output: `output results <df_output_var>`
        - df_state_final: `final model states <df_state_var>`

    Examples
    --------

    >>> df_output, df_state_final = supy.run_supy(df_forcing, df_state_init)


    """

    # save df_init without changing its original data
    # df.copy() in pandas works as a standard python deepcopy
    df_init = df_state_init.copy()

    # retrieve the last temporal record as `df_init`
    # if a `datetime` level existing in the index
    if df_init.index.nlevels > 1:
        idx_dt = df_init.index.get_level_values("datetime").unique()
        dt_last = idx_dt.max()
        df_init = df_init.loc[dt_last]

    # add placeholder variables for df_forcing
    # `metforcingdata_grid` and `ts5mindata_ir` are used by AnOHM and ESTM, respectively
    # they are now temporarily disabled in supy
    df_forcing = df_forcing.assign(metforcingdata_grid=0, ts5mindata_ir=0,).rename(
        # rename is a workaround to resolve naming inconsistency between
        # suews fortran code interface and input forcing file headers
        columns={
            "%" + "iy": "iy",
            "id": "id",
            "it": "it",
            "imin": "imin",
            "qn": "qn1_obs",
            "qh": "qh_obs",
            "qe": "qe",
            "qs": "qs_obs",
            "qf": "qf_obs",
            "U": "avu1",
            "RH": "avrh",
            "Tair": "temp_c",
            "pres": "press_hpa",
            "rain": "precip",
            "kdown": "avkdn",
            "snow": "snowfrac_obs",
            "ldown": "ldown_obs",
            "fcld": "fcld_obs",
            "Wuh": "wu_m3",
            "xsmd": "xsmd",
            "lai": "lai_obs",
            "kdiff": "kdiff",
            "kdir": "kdir",
            "wdir": "wdir",
        }
    )
    # reorder columns of df_forcing to comply with SUEWS kernel convention in receiving the input
    # TODO: this re-ordering can be later put into the planned input checker
    list_var_forcing = [
        "iy",
        "id",
        "it",
        "imin",
        "qn1_obs",
        "qh_obs",
        "qe",
        "qs_obs",
        "qf_obs",
        "avu1",
        "avrh",
        "temp_c",
        "press_hpa",
        "precip",
        "avkdn",
        "snowfrac_obs",
        "ldown_obs",
        "fcld_obs",
        "wu_m3",
        "xsmd",
        "lai_obs",
        "kdiff",
        "kdir",
        "wdir",
        "isec",
        "metforcingdata_grid",
        "ts5mindata_ir",
    ]
    df_forcing = df_forcing.loc[:, list_var_forcing]

    # grid list determined by initial states
    list_grid = df_init.index

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
        for grid, series_state_init in df_init.iterrows()
    }

    # remove 'problems.txt'
    if Path("problems.txt").exists():
        os.remove("problems.txt")

    if save_state:
        # use slower more functional single step wrapper

        # convert df to dict with `itertuples` for better performance
        dict_forcing = {row.Index: row._asdict() for row in df_forcing.itertuples()}

        for tstep in df_forcing.index:
            # temporal loop
            # initialise output of tstep:
            # load met_forcing if the same across all grids:
            met_forcing_tstep = dict_forcing[tstep]
            # print(met_forcing_tstep.keys())
            # spatial loop
            for grid in list_grid:
                dict_state_start = dict_state[(tstep, grid)]
                # calculation at one step:
                try:
                    dict_state_end, dict_output_tstep = suews_cal_tstep(
                        dict_state_start, met_forcing_tstep
                    )
                except:
                    raise RuntimeError("SUEWS kernel error")

                # update output & model state at tstep for the current grid
                dict_output.update({(tstep, grid): dict_output_tstep})
                dict_state.update({(tstep + 1 * freq, grid): dict_state_end})

        # pack results as easier DataFrames
        df_output = pack_df_output(dict_output).swaplevel(0, 1)
        # drop unnecessary 'datetime' as it is already included in the index
        df_output = df_output.drop(columns=["datetime"], level=0)
        df_state_final = pack_df_state(dict_state).swaplevel(0, 1)

    else:
        # for multi-year run, reduce the whole df_forcing into {chunk_day}-day chunks for less memory consumption
        idx_start = df_forcing.index.min()
        idx_all = df_forcing.index
        grp_forcing_chunk = df_forcing.groupby(
            (idx_all - idx_start) // pd.Timedelta(chunk_day, "d")
        )
        if len(grp_forcing_chunk) > 1:
            df_state_init_chunk = df_state_init.copy()
            list_df_output = []
            list_df_state = []
            for grp in grp_forcing_chunk.groups:
                # get forcing of a specific year
                df_forcing_chunk = grp_forcing_chunk.get_group(grp)
                # run supy: actual execution done in the `else` clause below
                df_output_yr, df_state_final_chunk = run_supy_ser(
                    df_forcing_chunk, df_state_init_chunk, chunk_day=chunk_day
                )
                df_state_init_chunk = df_state_final_chunk.copy()
                # collect results
                list_df_output.append(df_output_yr)
                list_df_state.append(df_state_final_chunk)

            # re-organise results of each year
            df_output = pd.concat(list_df_output).sort_index()
            df_state_final = pd.concat(list_df_state).sort_index().drop_duplicates()
            return df_output, df_state_final

        else:
            # for single-chunk run (1 chunk = {chunk_day} years), directly put df_forcing into supy_driver for calculation
            # use higher level wrapper that calculate at a `block` level
            # for better performance

            # # construct input list for `Pool.starmap`
            # construct input list for `dask.bag`
            list_input = [
                # (dict_state[(tstep_init, grid)], df_forcing)
                dict_state[(tstep_init, grid)]
                for grid in list_grid
            ]

            try:
                list_res = [
                    suews_cal_tstep_multi(input_grid, df_forcing)
                    for input_grid in list_input
                ]
                list_state_end, list_output_array = zip(*list_res)

            except:
                raise RuntimeError("SUEWS kernel error")

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
            df_output = df_output0.replace(-999.0, np.nan)
            df_state_final = pack_df_state(dict_state).swaplevel(0, 1)

    # drop ESTM for now as it is not supported yet
    # select only those supported output groups
    list_group_use = [
        group for group in df_output.columns.levels[0] if group not in ["ESTM"]
    ]
    df_output = df_output.loc[:, list_group_use]
    # trim multi-index based columns
    df_output.columns = df_output.columns.remove_unused_levels()

    # pack final model states into a proper dataframe
    df_state_final = pack_df_state_final(df_state_final, df_init)

    # show simulation time
    # end = time.time()
    # print(f'Execution time: {(end - start):.1f} s')
    # print(f'====================\n')

    return df_output, df_state_final


def run_save_supy(
    df_forcing_tstep, df_state_init_m, ind, save_state, n_yr, path_dir_temp
):
    # run supy in serial mode
    df_output, df_state_final = run_supy_ser(
        df_forcing_tstep, df_state_init_m, save_state, n_yr
    )
    # save to path_dir_temp
    path_out = path_dir_temp / f"{ind}_out.h5"
    path_state = path_dir_temp / f"{ind}_state.h5"
    df_output.to_hdf(path_out, f"out_{ind}", mode="w")
    df_state_final.to_hdf(path_state, f"state_{ind}", mode="w")


# parallel mode: only used on Linux/macOS; Windows is not supported yet.
def run_supy_par(df_forcing_tstep, df_state_init_m, save_state, n_yr):
    n_grid = df_state_init_m.index.size
    list_forcing = [df_forcing_tstep for _ in range(n_grid)]
    list_state = [df_state_init_m.iloc[[i]] for i in np.arange(n_grid)]
    list_save_state = [save_state for _ in range(n_grid)]
    list_n_yr = [n_yr for _ in range(n_grid)]

    # create a temp directory for results
    with tempfile.TemporaryDirectory() as dir_temp:
        path_dir_temp = Path(dir_temp).resolve()
        # print(path_dir_temp)
        list_dir_temp = [path_dir_temp for _ in range(n_grid)]

        # parallel run
        with multiprocess.Pool() as pool:
            pool.starmap(
                run_save_supy,
                zip(
                    list_forcing,
                    list_state,
                    np.arange(n_grid),
                    list_save_state,
                    list_n_yr,
                    list_dir_temp,
                ),
            )

        # load dumped h5 files
        df_output = pd.concat(
            [
                pd.read_hdf(path_dir_temp / f"{n}_out.h5", f"out_{n}")
                for n in np.arange(n_grid)
            ]
        )
        df_state_final = pd.concat(
            [
                pd.read_hdf(path_dir_temp / f"{n}_state.h5", f"state_{n}")
                for n in np.arange(n_grid)
            ]
        )
        # print(list(path_dir_temp.glob('*')))

    return df_output, df_state_final


# main calculation end here
##############################################################################


# pack one Series of var into np.array
def pack_var(var_ser):
    dim = np.array(literal_eval(var_ser.index[-1])) + 1
    val = np.array(var_ser.values.reshape(dim), order="F")
    return val


# pack one Series of grid vars into dict of `np.array`s
def pack_grid_dict(ser_grid):
    ser_dtype = df_var_info.dtype
    list_var_int = df_var_info[(ser_dtype == "int") | (ser_dtype == "array('i')")].index
    list_var = ser_grid.index.levels[0].unique()
    # pack according to dimension info
    dict_var = {
        var: pack_var(ser_grid[var])  # .astype(np.float)
        for var in list_var
        if var not in ["file_init"]
    }
    # convert to int
    dict_var_int = {
        var: dict_var[var].astype(int) for var in list_var if var in list_var_int
    }
    dict_var.update(dict_var_int)
    return dict_var


# pack final state to the same format as initial state
def pack_df_state_final(df_state_end, df_state_start):
    ser_col_multi = df_state_start.columns.to_series()
    idx = df_state_end.index
    size_idx = idx.size

    dict_packed = {}
    for var in df_state_end.to_dict():
        # print(var)
        # print(df_state_end[var].values.shape)
        # reshape values to (number of columns, number of grids)
        val_flatten = np.concatenate(df_state_end[var].values).ravel()
        val = val_flatten.reshape((size_idx, -1)).T
        col_names = ser_col_multi[var].values
        dict_var = dict(zip(col_names, val))
        dict_packed.update(dict_var)

    df_state_end_packed = pd.DataFrame(dict_packed, index=idx)
    df_state_end_packed.columns.set_names(["var", "ind_dim"], inplace=True)

    # swap index levels to form: {datetime, grid}
    # so using loc to retrieve the last index can get a dataframe for a restart run
    df_state_end_packed = df_state_end_packed.swaplevel()
    # df_state_end_packed.index.set_names('grid', inplace=True)

    return df_state_end_packed
