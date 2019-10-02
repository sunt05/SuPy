# functions to check validity of forcing and state DataFrames
#import yaml
#yaml.warnings({'YAMLLoadWarning': False})
import sys
import logging
from typing import List, Tuple
import xarray as xr
import numpy as np
import pandas as pd
import json
import glob
import os
from ._env import path_supy_module, logger_supy

# %%
# the check list file with ranges and logics
path_rules = path_supy_module/'checker_rules.json'


# opening the check list file
def load_rules(path_rules):

    with open(path_rules) as cf:
        dict_rules = json.load(cf)

    # making the keys lowercase to be consistent with supy
    dict_rules_lower = {}
    for key in dict_rules.keys():
        # for some reason pop() did not work here!!
        dict_rules_lower[key.lower()] = dict_rules[key]

    return dict_rules_lower


# checking the range of each parameter
def check_range(ser_to_check: pd.Series, rule_var: dict) -> Tuple:

    var = ser_to_check.name.lower()

    min_v = rule_var[var]['param']['min']
    max_v = rule_var[var]['param']['max']
    max_v = -np.inf if isinstance(min_v, str) else min_v
    max_v = np.inf if isinstance(max_v, str) else max_v
    description = f'{var} should be between {min_v} and {max_v}'
    is_accepted_flag = False

    for value in np.nditer(ser_to_check.values):
        if min_v <= value <= max_v:
            is_accepted_flag = True

    if(not is_accepted_flag):
        is_accepted = 'No'
        suggestion = 'change the parameter to fall into the acceptable range'
    else:
        is_accepted = 'Yes'
        suggestion = ''

    return var, is_accepted, description, suggestion


def check_zd_zh(var, values, cr):

    return 0


# # checks for suews parameters
# def check_var_suews(var, values, cr, df_sum):

#     logic = cr[var]['logic']

#     if logic == 'range':
#         out_list = check_range(var, values, cr)
#     elif logic == 'zd-zh':
#         out_list = check_zd_zh(var, values, cr)

#     df_sum.loc[len(df_sum)] = out_list

#     return df_sum


def check_forcing(df_forcing:pd.DataFrame):
    # collect issues
    list_issues = []
    flag_valid = True
    # check the following:
    # 1. correct columns
    list_col_valid = [
        'iy', 'id', 'it', 'imin',
        'qn', 'qh', 'qe', 'qs', 'qf',
        'U', 'RH', 'Tair', 'pres', 'rain',
        'kdown', 'snow', 'ldown', 'fcld',
        'Wuh', 'xsmd',
        'lai', 'kdiff', 'kdir', 'wdir',
        'isec',
    ]
    col_df = df_forcing.columns
    # 1.1 if all columns are present
    set_diff = set(list_col_valid).difference(col_df)
    if len(set_diff) > 0:
        str_issue = f'Missing columns found: {set_diff}'
        list_issues.append(str_issue)
        flag_valid = False
    # 1.2 if all columns are in right position
    for col_v, col in zip(list_col_valid, col_df):
        if col_v != col:
            str_issue = f'Column {col} is not in the valid position'
            list_issues.append(str_issue)
            flag_valid = False

    # 2. valid timestamps:
    ind_df = df_forcing.index
    # 2.1 must be a temporal index
    if not isinstance(ind_df, pd.DatetimeIndex):
        str_issue = f'Index must be {pd.DatetimeIndex}'
        list_issues.append(str_issue)
        flag_valid = False
    # 2.2 no duplicates
    if ind_df.has_duplicates:
        ind_dup = ind_df[ind_df.duplicated()]
        str_issue = f'Timestamps have duplicates: {ind_dup}'
        list_issues.append(str_issue)
        flag_valid = False

    # 2.3 monotonically increasing
    if not ind_df.is_monotonic_increasing:
        str_issue = f'Timestamps must be monotonically increasing'
        list_issues.append(str_issue)
        flag_valid = False

    # 2.4 must have a valid `freq` attribute
    if hasattr(ind_df,'freq'):
        if ind_df.freq is None:
            str_issue = f'Temporal index must have a valie `freq`'
            list_issues.append(str_issue)
            flag_valid = False
    else:
        str_issue = f'Temporal index must have `freq` attribute'
        list_issues.append(str_issue)
        flag_valid = False


    # 3. valid physical ranges
    dict_rules = load_rules(path_rules)
    for var in col_df:
        if var not in ['iy', 'id', 'it', 'imin', 'isec']:
            ser_var=df_forcing[var]
            res_check = check_range(ser_var, dict_rules)
            if not res_check[1]:
                str_issue = res_check[2]
                list_issues.append(str_issue)
                flag_valid = False

    if not flag_valid:
        logger_supy.error('Issues found:')
        str_issue = '\n'.join(list_issues)
        logger_supy.error(str_issue)
        return list_issues
    else:
        logger_supy.info('All checks passed!')
