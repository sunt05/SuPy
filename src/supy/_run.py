import copy
import sys
import os
import traceback
from ast import literal_eval

import numpy as np
import pandas as pd

from supy_driver import suews_driver as sd

from ._load import (df_var_info, list_var_inout,
                        list_var_inout_multitsteps, list_var_input,
                        list_var_input_multitsteps, list_var_output,
                        list_var_output_multitsteps)

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
        print(traceback.format_exc())
        # show SUEWS fatal error details produced by SUEWS kernel
        with open('problems.txt', 'r') as f:
            print(f.read())
        # clean slate
        # os.remove('problems.txt')
        # sys.exit()
        return 'SUEWS kernel error'
    else:
        # update state variables
        # if save_state:  # deep copy states results
        dict_state_end = copy.deepcopy(dict_state_start)
        dict_state_end.update(
            {
                var: copy.copy(dict_input[var])
                for var in list_var_inout
            }
        )

        # update timestep info
        dict_state_end['tstep_prev'] = dict_state_end['tstep']
        dict_state_end['dt_since_start'] += dict_state_end['tstep']

        # pack output
        dict_output = {k: v for k, v in zip(
            list_var_output, res_suews_tstep)}

        return dict_state_end, dict_output


# high-level wrapper: suews_cal_tstep
# def suews_cal_tstep_multi(df_state_start_grid, df_met_forcing_block):
def suews_cal_tstep_multi(dict_state_start_grid, df_met_forcing_block):
    # use single dict as input for suews_cal_main
    # dict_input = df_state_start_grid.copy().to_dict()
    dict_input = copy.deepcopy(dict_state_start_grid)
    dict_input.update({
        'metforcingblock': np.array(
            df_met_forcing_block.drop(
                columns=[
                    'metforcingdata_grid',
                    'ts5mindata_ir',
                    'isec',
                ]),
            order='F'
        ),
        'ts5mindata_ir': np.array(
            df_met_forcing_block['ts5mindata_ir'],
            order='F'
        ),
        'len_sim': np.array(df_met_forcing_block.shape[0], dtype=int)})
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
        with open('problems.txt','r') as f:
            print(f.read())
        # clean slate
        # os.remove('problems.txt')
        # sys.exit()
        # raise RuntimeError("Something bad happened") from exs
        return 'SUEWS kernel error'
    else:
        # update state variables
        # dict_state_end = copy.copy(dict_input)
        dict_state_end = copy.deepcopy(dict_state_start_grid)
        dict_state_end.update(
            {
                var: dict_input[var] for var in list_var_inout_multitsteps
            }
        )

        # update timestep info
        dict_state_end['tstep_prev'] = dict_state_end['tstep']
        idx_dt = df_met_forcing_block.index
        duration_s = int((idx_dt[-1] - idx_dt[0]).total_seconds())
        dict_state_end['dt_since_start'] += duration_s + dict_state_end['tstep']

        # pack output
        dict_output_array = {k: v for k, v in zip(
            list_var_output[1:], res_suews_tstep_multi)}

        return dict_state_end, dict_output_array

# main calculation end here
##############################################################################


# pack one Series of var into np.array
def pack_var(var_ser):
    dim = np.array(literal_eval(var_ser.index[-1])) + 1
    val = np.array(var_ser.values.reshape(dim), order='F')
    return val


# pack one Series of grid vars into dict of `np.array`s
def pack_grid_dict(ser_grid):
    ser_dtype = df_var_info.dtype
    list_var_int = df_var_info[
        (ser_dtype == 'int') | (ser_dtype == "array('i')")].index
    list_var = ser_grid.index.levels[0].unique()
    # pack according to dimension info
    dict_var = {
        var: pack_var(ser_grid[var])\
        # .astype(np.float)
        for var in list_var if var not in ['file_init']
    }
    # convert to int
    dict_var_int = {
        var: dict_var[var].astype(int)
        for var in list_var if var in list_var_int
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
    df_state_end_packed.columns.set_names(['var', 'ind_dim'], inplace=True)

    # swap index levels to form: {datetime, grid}
    # so using loc to retrieve the last index can get a dataframe for a restart run
    df_state_end_packed=df_state_end_packed.swaplevel()
    # df_state_end_packed.index.set_names('grid', inplace=True)

    return df_state_end_packed
