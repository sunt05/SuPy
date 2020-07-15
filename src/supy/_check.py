# functions to check validity of forcing and state DataFrames
# import yaml
# yaml.warnings({'YAMLLoadWarning': False})
import json
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from ._env import logger_supy, path_supy_module
from ._load import dict_var_type_forcing

# the check list file with ranges and logics
path_rules_indiv = path_supy_module / "checker_rules_indiv.json"


# opening the check list file
def load_rules(path_rules) -> Dict:

    with open(path_rules) as cf:
        dict_rules = json.load(cf)

    # making the keys lowercase to be consistent with supy
    dict_rules_lower = {}
    for key in dict_rules.keys():
        # for some reason pop() did not work here!!
        dict_rules_lower[key.lower()] = dict_rules[key]

    return dict_rules_lower


# store rules as a dict
dict_rules_indiv = load_rules(path_rules_indiv)

# checking the range of each parameter
def check_range(ser_to_check: pd.Series, rule_var: dict) -> Tuple:

    var = ser_to_check.name.lower()

    min_v = rule_var[var]["param"]["min"]
    max_v = rule_var[var]["param"]["max"]
    min_v = -np.inf if isinstance(min_v, str) else min_v
    max_v = np.inf if isinstance(max_v, str) else max_v
    description = ""
    is_accepted_flag = False

    for ind, value in ser_to_check.items():
        if min_v <= value <= max_v:
            is_accepted_flag = True
        elif value == -999.0:
            # default `na` value as such option is unnecessary in SUEWS
            is_accepted_flag = True
        else:
            is_accepted_flag = False
            description = f"`{var}` should be between [{min_v}, {max_v}] but `{value}` is found at {ind}"
            break

    if not is_accepted_flag:
        is_accepted = is_accepted_flag
        suggestion = "change the parameter to fall into the acceptable range"
    else:
        is_accepted = is_accepted_flag
        suggestion = ""

    return var, is_accepted, description, suggestion


def check_zd_zh(var, values, cr):

    return 0


# check if a valid method is set
def check_method(ser_to_check: pd.Series, rule_var: dict) -> Tuple:
    var = ser_to_check.name.lower()

    list_val = rule_var[var]["param"]["allowed"]
    description = ""

    is_accepted_flag = False
    for value in np.nditer(ser_to_check.values):
        if value in list_val:
            is_accepted_flag = True
        else:
            description = f"`{var}` should be one of {list_val} but is set as `{value}`"

    if not is_accepted_flag:
        is_accepted = is_accepted_flag
        suggestion = "change the parameter to an allowed value"
    else:
        is_accepted = is_accepted_flag
        suggestion = ""

    return var, is_accepted, description, suggestion


# # checks for suews parameters
# def check_var_suews(var, values, cr, df_sum):

#     logic = cr[var]['logic']

#     if logic == 'range':
#         out_list = check_range(var, values, cr)
#     elif logic == 'zd-zh':
#         out_list = check_zd_zh(var, values, cr)

#     df_sum.loc[len(df_sum)] = out_list

#     return df_sum
list_col_forcing = list(dict_var_type_forcing.keys())


def check_forcing(df_forcing: pd.DataFrame):
    logger_supy.info("SuPy is validating `df_forcing`...")
    # collect issues
    list_issues = []
    flag_valid = True
    # check the following:
    # 1. correct columns
    col_df = df_forcing.columns
    # 1.1 if all columns are present
    set_diff = set(list_col_forcing).difference(col_df)
    if len(set_diff) > 0:
        str_issue = f"Missing columns found: {set_diff}"
        list_issues.append(str_issue)
        flag_valid = False
    # 1.2 if all columns are in right position
    for col_v, col in zip(list_col_forcing, col_df):
        if col_v != col:
            str_issue = f"Column {col} is not in the valid position"
            list_issues.append(str_issue)
            flag_valid = False

    # 2. valid timestamps:
    ind_df = df_forcing.index
    # 2.1 must be a temporal index
    if not isinstance(ind_df, pd.DatetimeIndex):
        str_issue = f"Index must be {pd.DatetimeIndex}"
        list_issues.append(str_issue)
        flag_valid = False
    # 2.2 no duplicates
    if ind_df.has_duplicates:
        ind_dup = ind_df[ind_df.duplicated()]
        str_issue = f"Timestamps have duplicates: {ind_dup}"
        list_issues.append(str_issue)
        flag_valid = False

    # 2.3 monotonically increasing
    if not ind_df.is_monotonic_increasing:
        str_issue = f"Timestamps must be monotonically increasing"
        list_issues.append(str_issue)
        flag_valid = False

    # 2.4 must have a valid `freq` attribute
    if hasattr(ind_df, "freq"):
        if ind_df.freq is None:
            str_issue = f"Temporal index must have a valid `freq`"
            list_issues.append(str_issue)
            flag_valid = False
    else:
        str_issue = f"Temporal index must have `freq` attribute"
        list_issues.append(str_issue)
        flag_valid = False

    # 3. valid physical ranges
    for var in col_df:
        if var not in ["iy", "id", "it", "imin", "isec"]:
            ser_var = df_forcing[var]
            res_check = check_range(ser_var, dict_rules_indiv)
            if not res_check[1]:
                str_issue = res_check[2]
                list_issues.append(str_issue)
                flag_valid = False

    if not flag_valid:
        str_issue = "\n".join(["Issues found in `df_forcing`:"] + list_issues)
        logger_supy.error(str_issue)
        return list_issues
    else:
        logger_supy.info("All checks for `df_forcing` passed!")


def check_state(df_state: pd.DataFrame) -> List:
    logger_supy.info("SuPy is validating `df_state`...")
    # collect issues
    list_issues = []
    flag_valid = True
    list_col_state = set(dict_rules_indiv.keys()).difference(
        [x.lower() for x in list_col_forcing]
    )

    # check the following:
    # 1. correct columns
    col_df = df_state.columns.get_level_values("var")
    # 1.1 if all columns are present
    set_diff = set(list_col_state).difference(col_df)
    if len(set_diff) > 0:
        str_issue = f"Mandatory columns missing from df_state: {set_diff}"
        list_issues.append(str_issue)
        flag_valid = False
    # 1.2 if all columns are included in the checking list
    set_diff = set(col_df).difference(list_col_state)
    if len(set_diff) > 0:
        str_issue = f"Columns not included in checking list: {set_diff}"
        list_issues.append(str_issue)
        flag_valid = False

    # 2. check based on logic types
    list_to_check = set(col_df).intersection(list_col_state)
    for var in list_to_check:
        # pack
        val = dict_rules_indiv[var]
        df_var = df_state[var]
        # 'NA' implies no checking required
        if val["logic"] != "NA":
            pass
        if val["logic"] == "range":
            for ind in df_var.index:
                ser_var = df_var.loc[ind].rename(var)
                res_check = check_range(ser_var, dict_rules_indiv)
                if not res_check[1]:
                    str_issue = res_check[2] + f" at index `{ind}`"
                    list_issues.append(str_issue)
                    flag_valid = False
        if val["logic"] == "method":
            for ind in df_var.index:
                ser_var = df_var.loc[ind].rename(var)
                res_check = check_method(ser_var, dict_rules_indiv)
                if not res_check[1]:
                    str_issue = res_check[2] + f" at index `{ind}`"
                    list_issues.append(str_issue)
                    flag_valid = False

    if not flag_valid:
        str_issue = "\n".join(["Issues found in `df_state`:"] + list_issues)
        logger_supy.error(str_issue)
        return list_issues
    else:
        logger_supy.info("All checks for `df_state` passed!")


# flatten columns from MultiIndex to Index with compound notation
def flatten_col(df_state: pd.DataFrame):
    # original MultiIndex columsn
    col_mi = df_state.columns
    # flattened columns
    col_flat = col_mi.map(
        lambda s: (
            "_".join(s)
            .replace("_0", "")
            .replace("(", "")
            .replace(", ", "_")
            .replace(",)", "")
            .replace(")", "")
        )
    )
    # replace columns with flattened ones
    df_state_flat = df_state.set_axis(col_flat)
    return df_state_flat
