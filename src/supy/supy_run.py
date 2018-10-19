import copy
from ast import literal_eval

import numpy as np
from suews_driver import suews_driver as sd

from .supy_load import (gen_suews_arg_info_df, get_args_suews,
                        get_args_suews_multitsteps)

##############################################################################
# main calculation
# 1. calculation code for one time step
# 2. compact wrapper for running a whole simulation


# 1. calculation code for one time step
# store these lists for later use
list_var_input = list(get_args_suews()['var_input'])
list_var_inout = list(get_args_suews()['var_inout'])
list_var_output = list(get_args_suews()['var_output'])
set_var_input = set(list_var_input)
set_var_inout = set(list_var_inout)
set_var_ouput = set(list_var_output)

list_var_input_multitsteps = list(get_args_suews_multitsteps()['var_input'])
list_var_inout_multitsteps = list(get_args_suews_multitsteps()['var_inout'])
list_var_output_multitsteps = list(get_args_suews_multitsteps()['var_output'])
set_var_input_multitsteps = set(list_var_input_multitsteps)
set_var_inout_multitsteps = set(list_var_inout_multitsteps)
set_var_ouput_multitsteps = set(list_var_output_multitsteps)


df_info_suews_cal_main = gen_suews_arg_info_df(sd.suews_cal_main.__doc__)
df_info_suews_cal_multitsteps = gen_suews_arg_info_df(
    sd.suews_cal_multitsteps.__doc__)

df_var_info = df_info_suews_cal_multitsteps.merge(
    df_info_suews_cal_main, how='outer').set_index('name')

# test for performance
dict_var_inout = {k: None for k in set_var_inout}


# high-level wrapper: suews_cal_tstep
def suews_cal_tstep(dict_state_start, dict_met_forcing_tstep,
                    save_state=False):
    # use single dict as input for suews_cal_main
    dict_input = dict_state_start.copy()
    dict_input.update(dict_met_forcing_tstep)
    dict_input = {k: dict_input[k] for k in list_var_input}

    # for var in dict_input:
    #     print(var)
    #     print(dict_input[var])
    #     print('\n')

    # main calculation:
    res_suews_tstep = sd.suews_cal_main(**dict_input)

    # update state variables
    if save_state:  # deep copy states results
        dict_state_end = dict_state_start.copy()
        dict_state_end.update({var: copy.copy(dict_input[var])
                               for var in list_var_inout})
    else:  # only reference to dict_state_start
        dict_state_end = dict_state_start

    # update timestep info
    dict_state_end['tstep_prev'] = dict_state_end['tstep']
    dict_state_end['dt_since_start'] += dict_state_end['tstep']

    # pack output
    dict_output = {k: v for k, v in zip(
        list_var_output, res_suews_tstep)}

    return dict_state_end, dict_output


# # high-level wrapper: suews_cal_tstep
# def suews_cal_tstep_multi(df_state_start_grid, df_met_forcing_block):
#     # use single dict as input for suews_cal_main
#     dict_input = df_state_start_grid.copy().to_dict()
#     dict_input.update({
#         'metforcingblock': np.array(
#             df_met_forcing_block.drop(
#                 'metforcingdata_grid', axis=1), order='F'),
#         'ts5mindata_ir': np.array(
#             df_met_forcing_block['ts5mindata_ir'], order='F'),
#         'len_sim': df_met_forcing_block.shape[0]})
#     # dict_input = {k: dict_input[k] for k in list_var_input}
#
#     # main calculation:
#     res_suews_tstep_multi = sd.suews_cal_multitsteps(**dict_input)
#
#     # update state variables
#     dict_state_end = df_state_start_grid.copy().to_dict()
#     dict_state_end.update({var: dict_input[var] for var in list_var_inout})
#
#     # update timestep info
#     dict_state_end['tstep_prev'] = dict_state_end['tstep']
#     dict_state_end['dt_since_start'] += dict_state_end['tstep']
#
#     # pack output
#     dict_output_array = {k: v for k, v in zip(
#         list_var_output[1:], res_suews_tstep_multi)}
#
#     return dict_state_end, dict_output_array


# high-level wrapper: suews_cal_tstep
# def suews_cal_tstep_multi(df_state_start_grid, df_met_forcing_block):
def suews_cal_tstep_multi(dict_state_start_grid, df_met_forcing_block):
    # use single dict as input for suews_cal_main
    # dict_input = df_state_start_grid.copy().to_dict()
    dict_input = copy.copy(dict_state_start_grid)
    dict_input.update({
        'metforcingblock': np.array(
            df_met_forcing_block.drop(
                columns=[
                    'metforcingdata_grid',
                    # 'ts5mindata_ir',
                    # 'isec',
                ]),
            order='F'),
        'ts5mindata_ir': np.array(
            df_met_forcing_block['ts5mindata_ir'], order='F'),
        'len_sim': np.array(df_met_forcing_block.shape[0], dtype=int)})
    dict_input = {k: dict_input[k] for k in list_var_input_multitsteps}

    # for var in dict_input:
    #     print(var)
    #     print(dict_input[var])
    #     print(dict_input[var].shape, '\n')

    # main calculation:
    res_suews_tstep_multi = sd.suews_cal_multitsteps(**dict_input)

    # update state variables
    # dict_state_end = copy.copy(dict_input)
    dict_state_end = copy.copy(dict_input)
    dict_state_end.update({var: dict_input[var] for var in list_var_inout})

    # update timestep info
    dict_state_end['tstep_prev'] = dict_state_end['tstep']
    dict_state_end['dt_since_start'] += dict_state_end['tstep']

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
def pack_grid_dict(grid_ser):
    ser_dtype = df_var_info.dtype
    list_var_int = df_var_info[
        (ser_dtype == 'int') | (ser_dtype == "array('i')")].index
    list_var = grid_ser.index.levels[0].unique()
    # pack according to dimension info
    dict_var = {
        var: pack_var(grid_ser[var])\
        # .astype(np.float)
        for var in list_var
        if var not in ['file_init']
    }
    # convert to int
    dict_var_int = {
        var: dict_var[var].astype(int)
        for var in list_var
        if var in list_var_int
    }
    dict_var.update(dict_var_int)
    return dict_var


# # archived
# def run_suews_dict(df_forcing, df_init, save_state=False):
#     # initialise dicts for holding results and model states
#     dict_state = {}
#     dict_output = {}
#     # start tstep retrived from forcing data
#     t_start = df_forcing.index[0]
#     # convert df to dict with `itertuples` for better performance
#     dict_forcing = {row.Index: row._asdict()
#                     for row in df_forcing.itertuples()}
#     # dict_forcing = {tstep: met_forcing_tstep.to_dict()
#     #                 for tstep, met_forcing_tstep
#     #                 in df_forcing.iterrows()}
#     # grid list determined by initial states
#     grid_list = df_init.index
#
#     # dict_state is used to save model states for later use
#     dict_state = {(t_start, grid): series_state_init.to_dict()
#                   for grid, series_state_init
#                   in copy.deepcopy(df_init).iterrows()}
#
#     # temporal loop
#     for tstep in df_forcing.index:
#         # print 'tstep at', tstep
#         # initialise output of tstep:
#         # load met_forcing if the same across all grids:
#         met_forcing_tstep = dict_forcing[tstep]
#         # spatial loop
#         for grid in grid_list:
#             dict_state_start = dict_state[(tstep, grid)]
#             # calculation at one step:
#             # series_state_end, series_output_tstep = suews_cal_tstep_df(
#             #     series_state_start, met_forcing_tstep)
#             dict_state_end, dict_output_tstep = suews_cal_tstep(
#                 dict_state_start, met_forcing_tstep,
#                 save_state=save_state)
#             # update output & model state at tstep for the current grid
#             dict_output.update({(tstep, grid): dict_output_tstep})
#             dict_state.update({(tstep + 1, grid): dict_state_end})
#
#     # pack results as easier DataFrames
#     # df_output = pack_df_output(dict_output)
#     # df_state = pack_df_state(dict_state)
#     # df_output = dict_output
#     # df_state = dict_state
#
#     # return df_output, df_state
#     return dict_output, dict_state
