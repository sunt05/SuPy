from __future__ import division, print_function

# import collections
import functools
# from scipy import interpolate
# import copy
import glob
# import inspect
import os
from datetime import timedelta

# import sys
import numpy as np
import pandas as pd
from suews_driver import suews_driver as sd

import f90nml

from .supy_env import path_supy_module
from .supy_misc import path_insensitive

# from pathlib import Path


########################################################################
# get_args_suews can get the interface informaiton
# of the f2py-converted Fortran interface
def get_args_suews(docstring=sd.suews_cal_main.__doc__):
    # split doc lines for processing
    docLines = np.array(
        docstring.splitlines(),
        dtype=str)

    # get the information of input variables for SUEWS_driver
    posInput = np.where(
        np.logical_or(
            docLines == 'Parameters', docLines == 'Returns'))
    varInputLines = docLines[posInput[0][0] + 2:posInput[0][1] - 1]
    varInputInfo = np.array([[xx.rstrip() for xx in x.split(':')]
                             for x in varInputLines])
    dict_InputInfo = {xx[0]: xx[1] for xx in varInputInfo}
    dict_InOutInfo = {xx[0]: xx[1] for xx in varInputInfo if 'in/out' in xx[1]}

    # get the information of output variables for SUEWS_driver
    posOutput = np.where(docLines == 'Returns')
    varOutputLines = docLines[posOutput[0][0] + 2:]
    varOutputInfo = np.array([[xx.rstrip() for xx in x.split(':')]
                              for x in varOutputLines])
    dict_OutputInfo = {xx[0]: xx[1] for xx in varOutputInfo}

    # pack in/out results:
    dict_inout_sd = {
        # 'input' and 'output' are dict's that store variable information:
        # 1. intent: e.g., input, in/output
        # 2. dimension: e.g., (366,7)
        'input': dict_InputInfo,
        'output': dict_OutputInfo,
        # 'var_input' and 'var_output' are tuples,
        # that keep the order of arguments as in the Fortran subroutine
        'var_input': tuple(varInputInfo[:, 0]),
        'var_inout': tuple(dict_InOutInfo.keys()),
        'var_output': tuple(varOutputInfo[:, 0])}

    return dict_inout_sd


# for `suews_cal_multitsteps`
def get_args_suews_multitsteps():
    # split doc lines for processing
    docLines = np.array(
        sd.suews_cal_multitsteps.__doc__.splitlines(),
        dtype=str)

    # get the information of input variables for SUEWS_driver
    posInput = np.where(
        np.logical_or(
            docLines == 'Parameters', docLines == 'Returns'))
    varInputLines = docLines[posInput[0][0] + 2:posInput[0][1] - 1]
    varInputInfo = np.array([[xx.rstrip() for xx in x.split(':')]
                             for x in varInputLines])
    dict_InputInfo = {xx[0]: xx[1] for xx in varInputInfo}
    dict_InOutInfo = {xx[0]: xx[1] for xx in varInputInfo if 'in/out' in xx[1]}

    # get the information of output variables for SUEWS_driver
    posOutput = np.where(docLines == 'Returns')
    varOutputLines = docLines[posOutput[0][0] + 2:]
    varOutputInfo = np.array([[xx.rstrip() for xx in x.split(':')]
                              for x in varOutputLines])
    dict_OutputInfo = {xx[0]: xx[1] for xx in varOutputInfo}

    # pack in/out results:
    dict_inout_sd = {
        # 'input' and 'output' are dict's that store variable information:
        # 1. intent: e.g., input, in/output
        # 2. dimension: e.g., (366,7)
        'input': dict_InputInfo,
        'output': dict_OutputInfo,
        # 'var_input' and 'var_output' are tuples,
        # that keep the order of arguments as in the Fortran subroutine
        'var_input': tuple(varInputInfo[:, 0]),
        'var_inout': tuple(dict_InOutInfo.keys()),
        'var_output': tuple(varOutputInfo[:, 0])}

    return dict_inout_sd

##############################################################################
# input processor
# 1. surface properties will be retrieved and packed together for later use
# 2. met forcing conditions will splitted into time steps and used to derive
# other information


# descriptive list/dicts for variables
# minimal required input files for configuration:
list_file_input = [
    'SUEWS_AnthropogenicHeat.txt',
    'SUEWS_BiogenCO2.txt',
    'SUEWS_Conductance.txt',
    'SUEWS_ESTMCoefficients.txt',
    'SUEWS_Irrigation.txt',
    'SUEWS_NonVeg.txt',
    'SUEWS_OHMCoefficients.txt',
    'SUEWS_Profiles.txt',
    'SUEWS_SiteSelect.txt',
    'SUEWS_Snow.txt',
    'SUEWS_Soil.txt',
    'SUEWS_Veg.txt',
    'SUEWS_Water.txt',
    'SUEWS_WithinGridWaterDist.txt'
]

# library of all properties
dict_libVar2File = {fileX.replace('.txt', '').replace(
    'SUEWS', 'lib'): fileX for fileX in list_file_input
    if fileX.endswith('.txt')}

# dictionaries:
# links between code in SiteSelect to properties in according tables
# this is described in SUEWS online manual:
# http://urban-climate.net/umep/SUEWS#SUEWS_SiteSelect.txt
path_code2file = os.path.join(path_supy_module, 'code2file.json')
dict_Code2File = pd.read_json(path_code2file, typ='series').to_dict()
# variable translation as done in Fortran-SUEWS
path_var2siteselect = os.path.join(path_supy_module, 'var2siteselect.json')
dict_var2SiteSelect = pd.read_json(path_var2siteselect, typ='series').to_dict()

# expand dict_Code2File for retrieving surface characteristics
dict_varSiteSelect2File = {
    x: 'SUEWS_SiteSelect.txt' for x in dict_var2SiteSelect.keys()}
dict_Code2File.update(dict_varSiteSelect2File)


# extract metainfo of var from dict_info
def extract_var_info(var_name, dict_info):
    dict_var_info = {'name': var_name}
    var_line = dict_info[var_name]
    list_var_line = var_line.split()
    if 'array' in var_line:
        dict_var_info.update({'intent': list_var_line[0]})
        # print(list_var_line)
        rank = int(list_var_line[1].split('-')[-1])
        dict_var_info.update({'rank': rank})
        dtype = list_var_line[2]
        dict_var_info.update({'dtype': dtype})
        if rank > 0:
            bounds = list_var_line[-1]
            dict_var_info.update({'bounds': bounds})
        else:
            dict_var_info.update({'bounds': 0})
    else:
        dict_var_info.update({'intent': list_var_line[0]})
        rank = 0
        dict_var_info.update({'rank': rank})
        dtype = list_var_line[1]
        dict_var_info.update({'dtype': dtype})
        dict_var_info.update({'bounds': 0})
    return dict_var_info


# generate DataFrame of docstring for suews_wrappers generated by f2py
def gen_suews_arg_info_df(docstring):
    dict_info = get_args_suews(docstring)['input']
    # dict_info.update(get_args_suews(docstring)['output'])
    dict_info = {var: extract_var_info(var, dict_info) for var in dict_info}
    df_info = pd.DataFrame(dict_info).T
    return df_info


# load model settings
# load configurations: mod_config
# process RunControl.nml
# this function can handle all SUEWS nml files
@functools.lru_cache(maxsize=128)
def load_SUEWS_nml(xfile):
    # remove case issues
    xfile = path_insensitive(xfile)
    df = pd.DataFrame(f90nml.read(xfile))
    return df


def load_SUEWS_RunControl(xfile):
    lib_RunControl = load_SUEWS_nml(xfile)
    # return DataFrame containing settings
    return lib_RunControl


# load all tables (xgrid.e., txt files)
def load_SUEWS_table(fileX):
    # remove case issues
    fileX = path_insensitive(fileX)
    rawdata = pd.read_table(fileX, delim_whitespace=True,
                            comment='!', error_bad_lines=False,
                            skiprows=1, index_col=0).dropna()
    return rawdata


# load all tables into variables staring with 'lib_' and filename
@functools.lru_cache(maxsize=16)
def load_SUEWS_Libs(path_input):
    dict_libs = {}
    for lib, lib_file in dict_libVar2File.items():
        # print(lib_file)
        lib_path = os.path.join(path_input, lib_file)
        dict_libs.update({lib: load_SUEWS_table(lib_path)})
    # return DataFrame containing settings
    return dict_libs


# look up properties according to code
@functools.lru_cache(maxsize=None)
def lookup_code_lib(libCode, codeKey, codeValue, path_input):
    str_lib = dict_Code2File[libCode].replace(
        '.txt', '').replace('SUEWS', 'lib')
    dict_libs = load_SUEWS_Libs(path_input)
    lib = dict_libs[str_lib]
    # print(str_lib,codeKey,codeValue)
    if codeKey == ':':
        res = lib.loc[int(np.unique(codeValue)), :].tolist()
    else:
        res = lib.loc[int(np.unique(codeValue)), codeKey]
        if isinstance(res, pd.Series):
            # drop duolicates, otherwise duplicate values
            # will be introduced with more complexity
            res = res.drop_duplicates()
            res = res.values[0] if res.size == 1 else res.tolist()
    return res


# a recursive function to retrieve value based on key sequences
# @hash_dict
# @functools.lru_cache()
def lookup_KeySeq_lib(indexKey, subKey, indexCode, path_input):
    # print('\ntest lookup_KeySeq_lib')
    # print(indexKey, subKey, indexCode, type(subKey))
    # print('test end\n')
    if isinstance(subKey, float):
        res = subKey
    # elif indexKey == 'const':
    #     res = subKey
    elif type(subKey) is str:
        res = lookup_code_lib(indexKey, subKey, indexCode, path_input)
    elif isinstance(subKey, dict):
        # elif isinstance(subKey, HDict):
        res = []
        for indexKeyX, subKeyX in subKey.items():
            # print(indexKeyX, subKeyX)
            if indexKeyX == 'const':
                resX = subKeyX
            else:
                # indexKeyX, subKeyX = list(subKey.items())[0]
                indexCodeX = lookup_code_lib(
                    indexKey, indexKeyX, indexCode, path_input)
                if isinstance(subKeyX, list):
                    resX = [
                        lookup_KeySeq_lib(
                            indexKeyX, x_subKeyX, indexCodeX, path_input)
                        for x_subKeyX in subKeyX]
                else:
                    resX = lookup_KeySeq_lib(
                        indexKeyX, subKeyX, indexCodeX, path_input)
            res.append(resX)
        res = res[0] if len(res) == 1 else res
    elif isinstance(subKey, list):
        # elif isinstance(subKey, HList):
        res = []
        for subKeyX in subKey:
            indexCodeX = indexCode
            resX = lookup_KeySeq_lib(indexKey, subKeyX, indexCodeX, path_input)
            res.append(resX)
        res = res[0] if len(res) == 1 else res
    # final result
    return res


# load surface characteristics
def load_SUEWS_SurfaceChar(path_input):
    # load all libraries
    dict_libs = load_SUEWS_Libs(path_input)
    list_grid = dict_libs['lib_SiteSelect'].index
    # construct a dictionary in the form: {grid:{var:value,...}}
    dict_gridSurfaceChar = {
        grid: {
            k: lookup_KeySeq_lib(k, v, grid, path_input)
            for k, v in dict_var2SiteSelect.items()
        }
        for grid in list_grid}
    # convert the above dict to DataFrame
    df_gridSurfaceChar = pd.DataFrame.from_dict(dict_gridSurfaceChar).T
    # empty dict to hold updated values
    # dict_x_grid = {}
    # modify some variables to be compliant with SUEWS requirement
    for xgrid in df_gridSurfaceChar.index:
        # transpoe snowprof:
        df_gridSurfaceChar.at[xgrid, 'snowprof_24hr'] = np.array(
            df_gridSurfaceChar.at[xgrid, 'snowprof_24hr'], order='F').T

        # transpoe laipower:
        df_gridSurfaceChar.at[xgrid, 'laipower'] = np.array(
            df_gridSurfaceChar.at[xgrid, 'laipower'], order='F').T

        # select non-zero values for waterdist of water surface:
        x = np.array(df_gridSurfaceChar.at[xgrid, 'waterdist'][-1])
        df_gridSurfaceChar.at[xgrid, 'waterdist'][-1] = (
            x[np.nonzero(x)])

        # surf order as F:
        df_gridSurfaceChar.at[xgrid, 'storedrainprm'] = np.array(
            df_gridSurfaceChar.at[xgrid, 'storedrainprm'], order='F')

        # convert to np.array
        df_gridSurfaceChar.at[xgrid, 'alb'] = np.array(
            df_gridSurfaceChar.at[xgrid, 'alb'])

        # convert unit of `surfacearea` from ha to m^2
        df_gridSurfaceChar.at[xgrid, 'surfacearea'] = np.array(
            df_gridSurfaceChar.at[xgrid, 'surfacearea'] * 10000.)

        # print type(df_gridSurfaceChar.loc[xgrid, 'alb'])

        # dict holding updated values that can be converted to DataFrame later
        # dict_x = df_gridSurfaceChar.loc[xgrid, :].to_dict()
        # print 'len(dict_x)',len(dict_x['laipower'])

        # profiles:
        list_varTstep = [
            'ahprof_24hr',
            'popprof_24hr',
            'traffprof_24hr',
            'humactivity_24hr',
            'wuprofm_24hr',
            'wuprofa_24hr',
        ]
        df_gridSurfaceChar.loc[xgrid, list_varTstep] = \
            df_gridSurfaceChar.loc[xgrid, list_varTstep].map(np.transpose)

    # convert to DataFrame
    # df_x_grid = pd.DataFrame.from_dict(dict_x_grid).T
    df_x_grid = df_gridSurfaceChar
    return df_x_grid


def func_parse_date(year, doy, hour, min):
    # dt = datetime.datetime.strptime(
    #     ' '.join([year, doy, hour, min]), '%Y %j %H %M')
    dt = pd.to_datetime(' '.join(
        [str(k) for k in [year, doy, hour, min]]),
        format='%Y %j %H %M')
    return dt


def func_parse_date_row(row):
    [year, doy, hour, tmin] = row.loc[['iy', 'id', 'it', 'imin']]
    # dt = datetime.datetime.strptime(
    #     ' '.join([year, doy, hour, min]), '%Y %j %H %M')
    dt = pd.to_datetime(' '.join(
        [str(k) for k in [year, doy, hour, tmin]]),
        format='%Y %j %H %M')
    return dt


# calculate decimal time
def dectime(timestamp):
    t = timestamp
    dectime = (t.dayofyear - 1) + (
        t.hour + (t.minute + t.second / 60.) / 60.
    ) / 24
    return dectime

# resample solar radiation by zenith correction and total amount distribution


def resample_kdn(data_raw_kdn, tstep_mod, timezone, lat, lon, alt):
    # adjust solar radiation
    datetime_mid_local = data_raw_kdn.index - timedelta(
        seconds=tstep_mod / 2)
    sol_elev = np.array([sd.suews_cal_sunposition(
        t.year, dectime(t), timezone, lat, lon, alt)[-1]
        for t in datetime_mid_local])
    sol_elev_reset = np.zeros_like(sol_elev)
    sol_elev_reset[sol_elev <= 90] = 1.
    data_tstep_kdn_adj = sol_elev_reset * data_raw_kdn.copy()

    # rescale daily amounts
    avg_raw = data_raw_kdn.resample('D').mean()
    avg_tstep = data_tstep_kdn_adj.resample('D').mean()
    ratio_SWdown = (avg_raw / avg_tstep).reindex(
        index=avg_tstep.index).resample(
        '{tstep}S'.format(tstep=tstep_mod)).mean().fillna(method='pad')
    data_tstep_kdn_adj = ratio_SWdown * \
        data_tstep_kdn_adj.fillna(method='pad')

    return data_tstep_kdn_adj


# correct precipitation by even redistribution over resampled periods
def resample_precip(data_raw_precip, tstep_mod, tstep_in):
    ratio_precip = 1. * tstep_mod / tstep_in
    data_tstep_precip_adj = ratio_precip * data_raw_precip.copy().shift(
        -tstep_in + tstep_mod, freq='S').resample(
        '{tstep}S'.format(tstep=tstep_mod)).mean().interpolate(
        method='polynomial', order=0)
    data_tstep_precip_adj = data_tstep_precip_adj.fillna(value=0.)
    return data_tstep_precip_adj


# resample input forcing by linear interpolation
def resample_linear(data_raw, tstep_in, tstep_mod):
    # reset index as timestamps
    # data_raw.index = data_raw.loc[:, ['iy', 'id', 'it', 'imin']].apply(
    #     func_parse_date_row, 1)
    # shift by half-tstep_in to generate a time series with instantaneous
    # values
    data_raw_shift = data_raw.copy().shift(-tstep_in / 2, freq='S')

    # downscale input data to desired time step
    data_raw_tstep = data_raw_shift.resample(
        '{tstep}S'.format(tstep=tstep_mod)).interpolate(
        method='polynomial', order=1).rolling(
        window=2, center=False).mean()

    # reindex data_tstep to valid range
    # ix = pd.date_range(
    #     data_raw.index[0] - timedelta(seconds=tstep_in - tstep_mod),
    #     data_raw.index[-1],
    #     freq='{tstep}S'.format(tstep=tstep_mod))
    # data_tstep = data_raw_tstep.copy().reindex(
    #     index=ix).bfill().ffill().dropna()
    data_tstep = data_raw_tstep.copy().bfill().ffill().dropna()

    # correct temporal information
    data_tstep['iy'] = data_tstep.index.year
    data_tstep['id'] = data_tstep.index.dayofyear
    data_tstep['it'] = data_tstep.index.hour
    data_tstep['imin'] = data_tstep.index.minute

    return data_tstep


# resample input met foring to tstep required by model
def resample_forcing_met(
        data_met_raw, tstep_in, tstep_mod, lat, lon, alt, timezone, kdownzen):
    # overall resample by linear interpolation
    data_met_tstep = resample_linear(data_met_raw, tstep_in, tstep_mod)

    # adjust solar radiation by zenith correction and total amount distribution
    if kdownzen == 1:
        data_met_tstep["avkdn"] = resample_kdn(
            data_met_tstep["avkdn"], tstep_mod, timezone, lat, lon, alt)

    # correct rainfall
    data_met_tstep['precip'] = resample_precip(
        data_met_raw['precip'], tstep_mod, tstep_in)

    # # reset index with numbers
    # data_met_tstep_out = data_met_tstep.copy().reset_index(drop=True)

    return data_met_tstep


# load raw data: met forcing
def load_SUEWS_Forcing_met_df_raw(
        path_input, filecode, grid, tstep_met_in, multiplemetfiles):
    # file name pattern for met files
    forcingfile_met_pattern = os.path.join(
        path_input,
        '{site}{grid}_*_{tstep}.txt'.format(
            site=filecode,
            grid=(grid if multiplemetfiles == 1 else ''),
            tstep=int(tstep_met_in / 60)))
    df_forcing_met = load_SUEWS_Forcing_met_df_pattern(forcingfile_met_pattern)
    return df_forcing_met


# caching loaded met df for better performance in initialisation
@functools.lru_cache(maxsize=None)
def load_SUEWS_Forcing_met_df_pattern(forcingfile_met_pattern):
    """Short summary.

    Parameters
    ----------
    forcingfile_met_pattern : type
        Description of parameter `forcingfile_met_pattern`.

    Returns
    -------
    type
        Description of returned object.

    """
    # list of met forcing files
    forcingfile_met_pattern = os.path.expanduser(forcingfile_met_pattern)
    forcingfile_met_pattern = os.path.abspath(forcingfile_met_pattern)
    list_file_MetForcing = sorted([
        f for f in glob.glob(forcingfile_met_pattern)
        if 'ESTM' not in f])

    # print(forcingfile_met_pattern)
    # print(list_file_MetForcing)
    # load raw data
    df_forcing_met = pd.concat(
        [pd.read_table(
            fileX,
            delim_whitespace=True,
            comment='!',
            error_bad_lines=True
        ).dropna() for fileX in list_file_MetForcing],
        ignore_index=True)
    col_suews_met_forcing = [
        'iy', 'id', 'it', 'imin',
        'qn1_obs', 'qh_obs', 'qe', 'qs_obs', 'qf_obs',
        'avu1', 'avrh', 'temp_c', 'press_hpa', 'precip', 'avkdn',
        'snow_obs', 'ldown_obs', 'fcld_obs',
        'wuh', 'xsmd', 'lai_obs', 'kdiff', 'kdir', 'wdir'
    ]
    # rename these columns to match variables via the driver interface
    df_forcing_met.columns = col_suews_met_forcing

    # convert unit from kPa to hPa
    df_forcing_met['press_hpa'] *= 10

    # add `isec` for WRF-SUEWS interface
    df_forcing_met['isec'] = 0

    # set correct data types
    df_forcing_met[['iy', 'id', 'it', 'imin', 'isec']] = df_forcing_met[[
        'iy', 'id', 'it', 'imin', 'isec']].astype(np.int64)

    # set timestamp as index
    idx_dt = pd.date_range(
        *df_forcing_met.iloc[[0, -1], :4].astype(int).astype(str).apply(
            lambda ser: ser.str.cat(sep=' '), axis=1).map(
            lambda dt: pd.Timestamp.strptime(dt, '%Y %j %H %M')),
        periods=df_forcing_met.shape[0])

    df_forcing_met = df_forcing_met.set_index(idx_dt)

    return df_forcing_met


# load raw data: met forcing
def load_SUEWS_Forcing_ESTM_df_raw(
        path_input, filecode, grid, tstep_ESTM_in, multipleestmfiles):
    # file name pattern for met files
    forcingfile_ESTM_pattern = os.path.join(
        path_input,
        '{site}{grid}*{tstep}*txt'.format(
            site=filecode,
            grid=(grid if multipleestmfiles == 1 else ''),
            tstep=tstep_ESTM_in / 60))

    # list of met forcing files
    list_file_MetForcing = [
        f for f in glob.glob(forcingfile_ESTM_pattern)
        if 'ESTM' in f]

    # load raw data
    df_forcing_estm = pd.concat(
        [pd.read_table(
            fileX,
            delim_whitespace=True,
            comment='!',
            error_bad_lines=True
            # parse_dates={'datetime': [0, 1, 2, 3]},
            # keep_date_col=True,
            # date_parser=func_parse_date
        ).dropna() for fileX in list_file_MetForcing],
        ignore_index=True).rename(
        columns={
            '%' + 'iy': 'iy',
            'id': 'id',
            'it': 'it',
            'imin': 'imin',
            'kdown': 'avkdn',
            'RH': 'avrh',
            'U': 'avu1',
            'fcld': 'fcld_obs',
            'lai': 'lai_obs',
            'ldown': 'ldown_obs',
            'rain': 'precip',
            'pres': 'press_hpa',
            'qh': 'qh_obs',
            'qn': 'qn1_obs',
            'snow': 'snow_obs',
            'Tair': 'temp_c',
            # 'all': 'metforcingdata_grid',
            'xsmd': 'xsmd'})

    # set correct data types
    df_forcing_estm[['iy', 'id', 'it', 'imin']] = df_forcing_estm[[
        'iy', 'id', 'it', 'imin']].astype(np.int64)

    return df_forcing_estm


# TODO: add support for loading multi-grid forcing datasets
# def load_SUEWS_Forcing_df(dir_site, ser_mod_cfg, df_state_init):
#     pass


###############################################################################
# new initialisation functions below
# for better performance


# expand dict into paired tuples
def exp_dict2tuple(rec):
    if isinstance(rec, str):
        return rec
    elif isinstance(rec, float):
        # expand this for consistency in value indexing
        return [(rec, 'base')]
    elif isinstance(rec, dict):
        res = [
            # expand profile values to all hour indices
            (k, [str(i) for i in range(24)])
            if v == ':' else (k, exp_dict2tuple(v)) for k, v in rec.items()]
        # print(res)
        return res
    elif isinstance(rec, list):
        return [exp_dict2tuple(v) for v in rec]


# test if a list exists in a tuple
def list_in_tuple_Q(rec):
    if isinstance(rec, tuple):
        return any(isinstance(x, list) for x in rec)
    else:
        return False


# expand list of nested tuples to simple chained tuples
# (a,[b1,b2])->[(a,b1),(a,b2)]
def exp_list2tuple(rec):
    if isinstance(rec, list):
        if len(rec) == 1:
            return exp_list2tuple(rec[0])
        else:
            return [exp_list2tuple(x) for x in rec]
    elif list_in_tuple_Q(rec):
        for i, v in enumerate(rec):
            if isinstance(v, list):
                res = [(*rec[:i], x) if isinstance(x, (str, int))
                       else (*rec[:i], *x)
                       for x in v]
                res = [exp_list2tuple(x) for x in res]
                return exp_list2tuple(res)
    else:
        return rec


# collect tracing entries under a reference `code`
def gather_code_set(code, dict_var2SiteSelect):
    set_res = set()
    # print('\n code', code)
    for k, v in dict_var2SiteSelect.items():
        # print(set_res, type(set_res), k, v)
        if k == code:
            if isinstance(v, str):
                set_res.update([v])
            elif isinstance(v, list):
                set_res.update(v)
            elif isinstance(v, dict):
                set_res.update(v.keys())
            else:
                print(k, v)
        elif isinstance(v, dict):
            # print(res,v)
            set_res.update(list(gather_code_set(code, v)))
        elif isinstance(v, list):
            # print(set_res, v)
            for v_x in v:
                if isinstance(v_x, dict):
                    set_res.update(list(gather_code_set(code, v_x)))
    return set_res


# generate DataFrame based on reference `code` with index from `df_base`
def build_code_df(code, path_input, df_base):
    # str_lib = dict_Code2File[libCode].replace(
    #     '.txt', '').replace('SUEWS', 'lib')
    dict_libs = load_SUEWS_Libs(path_input)

    keys_code = gather_code_set(code, dict_var2SiteSelect)
    if keys_code == {':'}:
        # for profile values, extract all
        list_keys = code
    else:
        list_keys = [(code, k) for k in list(keys_code)]
    lib_code = dict_Code2File[code]
    str_lib = lib_code.replace('.txt', '').replace('SUEWS', 'lib')
    # lib = dict_libs[str_lib]
    df_code0 = pd.concat([dict_libs[str_lib]], axis=1, keys=[code])
    # df_siteselect = dict_libs['lib_SiteSelect']
    if isinstance(df_base.columns, pd.core.index.MultiIndex):
        code_base = df_base.columns.levels[0][0]
        code_index = (code_base, code)
    else:
        code_index = code

    # print('list_keys:\n', list_keys)
    df_code = df_code0.loc[df_base[code_index].astype(int).values, list_keys]
    df_code.index = df_base.index
    # recover outmost level code if profile values are extracted
    if keys_code == {':'}:
        df_code = pd.concat([df_code], axis=1, keys=[code])
    return df_code


# if to expand
def to_exp_Q(code_str):
    return ('Code' in code_str) or ('Prof' in code_str)


# generate df with all code-references values
@functools.lru_cache(maxsize=16)
def gen_all_code_df(path_input):
    dict_libs = load_SUEWS_Libs(path_input)
    df_siteselect = dict_libs['lib_SiteSelect']
    list_code = [
        code for code in df_siteselect.columns
        if to_exp_Q(code)
    ]

    # dict with code:{vars}
    dict_code_var = {
        code: gather_code_set(code, dict_var2SiteSelect)
        for code in list_code
        if len(gather_code_set(code, dict_var2SiteSelect)) > 0
    }

    # stage one expansion
    df_all_code = pd.concat(
        [build_code_df(code, path_input, df_siteselect)
         for code in dict_code_var],
        axis=1)

    return df_all_code


# generate df for all const based columns
@functools.lru_cache(maxsize=16)
def gen_all_const_df(path_input):
    dict_libs = load_SUEWS_Libs(path_input)
    df_siteselect = dict_libs['lib_SiteSelect']

    df_cst_all = df_siteselect.copy()
    dict_var_tuple = exp_dict_full(dict_var2SiteSelect)

    set_const_col = {
        v
        for x in dict_var_tuple.values() if isinstance(x, list)
        for v in x if 'const' in v
    }

    list_const_col = [x[1] for x in set_const_col]
    for cst in set_const_col:
        df_cst_all[cst[1]] = cst[1]
    df_cst_all = df_cst_all.loc[:, list_const_col]
    df_cst_all = pd.concat([df_cst_all], axis=1, keys=['base'])
    df_cst_all = pd.concat([df_cst_all], axis=1, keys=['const'])
    df_cst_all = df_cst_all.swaplevel(i=1, j=-1, axis=1)

    return df_cst_all


# build code expanded ensemble libs
@functools.lru_cache(maxsize=16)
def build_code_exp_df(path_input, code_x):
    # df with all code-references values
    df_all_code = gen_all_code_df(path_input)

    # df with code_x for further exapansion
    df_base_x = df_all_code[code_x]

    # columns in df_base_x for further expansion
    list_code_x = [code for code in df_base_x.columns if to_exp_Q(code)]
    # print(list_code_x)
    df_all_code_x = pd.concat(
        [build_code_df(code, path_input, df_base_x)
         for code in list_code_x],
        axis=1)
    df_all_code_x = pd.concat([df_all_code_x], axis=1, keys=[code_x])

    # add innermost lables to conform dimensionality
    df_base_x_named = pd.concat([df_base_x], axis=1, keys=[code_x])
    df_base_x_named = pd.concat([df_base_x_named], axis=1, keys=['base'])
    df_base_x_named = df_base_x_named.swaplevel(
        i=0, j=1, axis=1).swaplevel(i=1, j=-1, axis=1)

    # print(df_base_x_named.columns)
    # print(df_all_code_x.columns)
    df_res_exp = pd.concat([df_base_x_named, df_all_code_x], axis=1)

    return df_res_exp


@functools.lru_cache(maxsize=16)
def gen_df_siteselect_exp(path_input):
    # df with all code-references values
    df_all_code = gen_all_code_df(path_input)

    # retrieve all `Code`-relaed names
    dict_libs = load_SUEWS_Libs(path_input)
    df_siteselect = dict_libs['lib_SiteSelect']
    # list_code = [code for code in df_siteselect.columns if 'Code' in code]

    # dict with code:{vars}
    dict_code_var = {
        code: gather_code_set(code, dict_var2SiteSelect)
        for code in df_all_code.columns.levels[0]
    }

    # code with nested references
    dict_code_var_nest = {
        code: dict_code_var[code]
        for code in dict_code_var.keys()
        if any(('Code' in v) or ('Prof' in v) for v in dict_code_var[code])
    }

    # code with simple values
    dict_code_var_simple = {
        code: dict_code_var[code]
        for code in dict_code_var.keys()
        if code not in dict_code_var_nest
    }

    # stage 2: expand all code to actual values
    # nested code: code -> code -> value
    df_code_exp_nest = pd.concat(
        [build_code_exp_df(path_input, code)
         for code in list(dict_code_var_nest)],
        axis=1
    )
    # simple code: code -> value
    df_code_exp_simple = pd.concat(
        [pd.concat(
            [build_code_df(code, path_input, df_siteselect)
             for code in list(dict_code_var_simple)],
            axis=1
        )], axis=1, keys=['base'])
    df_code_exp_simple = df_code_exp_simple.swaplevel(
        i=0, j=1, axis=1)
    df_code_exp_simple = df_code_exp_simple.swaplevel(i=1, j=-1, axis=1)

    # combine both nested and simple df's
    df_code_exp = pd.concat([df_code_exp_simple, df_code_exp_nest], axis=1)

    # generate const based columns
    df_cst_all = gen_all_const_df(path_input)

    # df_siteselect_pad: pad levels for pd.concat
    df_siteselect_pad = pd.concat([df_siteselect], axis=1, keys=['base'])
    df_siteselect_pad = df_siteselect_pad.swaplevel(i=0, j=1, axis=1)
    df_siteselect_pad = pd.concat([df_siteselect_pad], axis=1, keys=['base'])
    df_siteselect_pad = df_siteselect_pad.swaplevel(i=0, j=1, axis=1)

    df_siteselect_exp = pd.concat(
        [df_siteselect_pad, df_code_exp, df_cst_all],
        axis=1)

    return df_siteselect_exp


def flatten_list(l):
    res = []
    if isinstance(l, list):
        for x in l:
            if isinstance(x, list):
                res += flatten_list(x)
            else:
                res.append(x)
    else:
        res = l
    return res


# pad a single record to three-element tuple for easier indexing
def pad_rec_single(rec):
    base_pad = ['base' for i in range(3)]
    if isinstance(rec, tuple) and len(rec) < 3:
        # print(rec)
        if all(isinstance(v, (str, float, int)) for v in rec):
            v_pad = tuple(list(rec) + base_pad)[:3]
        else:
            v_pad = rec
    elif isinstance(rec, str):
        v_pad = tuple([rec] + base_pad)[:3]
    else:
        v_pad = rec
    return v_pad


# recursive version to generate three-element tuples
def pad_rec(rec):
    if isinstance(rec, list):
        rec_pad = [pad_rec(v) for v in rec]
    else:
        rec_pad = pad_rec_single(rec)
    return rec_pad


# expand a dict to fully referenced tuples
def exp_dict_full(dict_var2SiteSelect):
    # expand dict to list of tuples
    dict_list_tuple = {k: exp_dict2tuple(v)
                       for k, v in dict_var2SiteSelect.items()}
    # expand list of nested tuples to simple chained tuples
    dict_var_tuple = {k: exp_list2tuple(v)
                      for k, v in dict_list_tuple.items()}

    # pad records to three-element tuples
    dict_var_tuple = {k: pad_rec(v) for k, v in dict_var_tuple.items()}
    return dict_var_tuple


# save this as system variable
# dict_var_tuple = exp_dict_full(dict_var2SiteSelect)


# generate expanded surface characteristics df
# this df contains all the actual values retrived from all tables
@functools.lru_cache(maxsize=16)
def gen_df_gridSurfaceChar_exp(path_input):
    df_siteselect_exp = gen_df_siteselect_exp(path_input)
    dict_var_tuple = exp_dict_full(dict_var2SiteSelect)
    df_gridSurfaceChar_exp = pd.concat(
        {
            k: df_siteselect_exp.loc[:, flatten_list(dict_var_tuple[k])]
            for k in dict_var_tuple
        },
        axis=1)
    return df_gridSurfaceChar_exp


# quickly load surface characteristics as DataFrame
# the dimension info of each variable can be derived from level 1 in columns
@functools.lru_cache(maxsize=16)
def load_SUEWS_SurfaceChar_df(path_input):
    df_gridSurfaceChar_exp = gen_df_gridSurfaceChar_exp(path_input)
    dict_var_ndim = {
        'ahprof_24hr': (24, 2),
        'humactivity_24hr': (24, 2),
        'laipower': (4, 3),
        'ohm_coef': (8, 4, 3),
        'popprof_24hr': (24, 2),
        'snowprof_24hr': (24, 2),
        'storedrainprm': (6, 7),
        'traffprof_24hr': (24, 2),
        'waterdist': (8, 6),
        'wuprofa_24hr': (24, 2),
        'wuprofm_24hr': (24, 2)
    }
    # df_gridSurfaceChar_x=df_gridSurfaceChar_exp.copy()

    len_grid_orig = df_gridSurfaceChar_exp.shape[0]
    df_x = df_gridSurfaceChar_exp
    if len_grid_orig == 1:
        df_gridSurfaceChar_exp = pd.concat([df_x, df_x])

    # update len_grid for reshaping
    list_var = df_gridSurfaceChar_exp.columns.levels[0]
    len_grid = df_gridSurfaceChar_exp.shape[0]
    list_grid = df_gridSurfaceChar_exp.index

    # here we add append a dummy `pos` level to df_gridSurfaceChar_exp
    # to avoid subset difficulty in subseting duplicate entries in multiindex
    df_x = df_gridSurfaceChar_exp.T
    # `pos` has unique values for each entry in columns
    df_x['pos'] = np.arange(len(df_gridSurfaceChar_exp.columns))
    # trasnpose back to have an extra column level `pos`
    df_gridSurfaceChar_exp = df_x.set_index('pos', append=True).T

    # extract column values and construct corresponding variables
    dict_gridSurfaceChar = {
        var: np.squeeze(df_gridSurfaceChar_exp[var].values)
        for var in list_var
    }

    # correct the unit for 'surfacearea'
    dict_gridSurfaceChar['surfacearea'] *= 1E4

    for var, dim in dict_var_ndim.items():
        # print(var)
        val = df_gridSurfaceChar_exp.loc[:, var].values
        if '_24hr' in var:
            dim_x = dim[-1::-1]
            # val = df_gridSurfaceChar_exp.loc[:, var].values
            val_x = val.reshape((len_grid, *dim_x)).transpose((0, 2, 1))
            dict_gridSurfaceChar.update({var: val_x})
        elif var == 'laipower':
            dim_x = dim[-1::-1]
            # val = df_gridSurfaceChar_exp.loc[:, var].values
            val_x = val.reshape((len_grid, *dim_x)).transpose((0, 2, 1))
            dict_gridSurfaceChar.update({var: val_x})
        elif var == 'storedrainprm':
            dim_x = dim
            # val = df_gridSurfaceChar_exp.loc[:, var].values
            val_x = val.reshape((len_grid, *dim_x))
            dict_gridSurfaceChar.update({var: val_x})
        elif var == 'ohm_coef':
            dim_x = dim
            # val = df_gridSurfaceChar_exp.loc[:, var].values
            val_x = val.reshape((len_grid, *dim_x))
            dict_gridSurfaceChar.update({var: val_x})
        elif var == 'waterdist':
            dim_x = dict_var_ndim[var]  # [-1::-1]
            # val = df_gridSurfaceChar_exp.loc[:, var].values
            val_x0 = val.reshape((len_grid, 9, 6))
            # directly load the values for commen land covers
            val_x1 = val_x0[:, :7]
            # process the ToSoilStore and ToRunoff entries
            val_x2 = val_x0[:, 7:].reshape(
                len_grid, 6, 2).sum(axis=2).reshape(-1, 1, 6)
            # combine valuees of common land convers and special cases
            val_x = np.hstack((val_x1, val_x2))
            dict_gridSurfaceChar.update({var: val_x})

    # convert to DataFrame
    dict_var = {}
    for var in list_var:
        val = dict_gridSurfaceChar[var]
        # ind_list = list(np.ndindex(val.shape[1:]))
        ind_list = [str(x) for x in np.ndindex(val.shape[1:])]
        df_var = pd.DataFrame(
            val.reshape((val.shape[0], -1)),
            index=list_grid,
            columns=ind_list if len(ind_list) > 1 else ['0'])
        dict_var.update({var: df_var})

    df_sfc_char = pd.concat(dict_var, axis=1)
    if len_grid_orig == 1:
        df_sfc_char = df_sfc_char.drop_duplicates()

    return df_sfc_char


# load RunControl variables
@functools.lru_cache(maxsize=16)
def load_SUEWS_dict_ModConfig(path_runcontrol):
    lib_RunControl = load_SUEWS_nml(path_runcontrol)
    dict_RunControl = lib_RunControl.loc[:, 'runcontrol'].to_dict()

    # mod_config: static properties
    dict_ModConfig = {
        'aerodynamicresistancemethod': 2,
        'evapmethod': 2,
        'laicalcyes': 1,
        'veg_type': 1,
        'diagnose': 0,
        'diagqn': 0,
        'diagqs': 0
    }
    dict_ModConfig.update(dict_RunControl)
    return dict_ModConfig


# load initial conditions of all grids as a DataFrame
@functools.lru_cache(maxsize=16)
def load_SUEWS_InitialCond_df(path_runcontrol):
    # initialise InitialCond_df with default values
    nan = -999.
    dict_InitCond_default = {
        'dayssincerain':  0,
        # `temp_c0` defaults to daily mean air temperature of the first day
        # 'temp_c0':  df_forcing_met['temp_c'].iloc[:(24 * nsh) - 1].mean(),
        'temp_c0':  10,
        'leavesoutinitially':  int(nan),
        'gdd_1_0':  nan,
        'gdd_2_0':  nan,
        'laiinitialevetr':  nan,
        'laiinitialdectr':  nan,
        'laiinitialgrass':  nan,
        'albevetr0':  nan,
        'albdectr0':  nan,
        'albgrass0':  nan,
        'decidcap0':  nan,
        'porosity0':  nan,
        'pavedstate':  0,
        'bldgsstate':  0,
        'evetrstate':  0,
        'dectrstate':  0,
        'grassstate':  0,
        'bsoilstate':  0,
        # 'waterstate':  df_gridSurfaceChar.at[grid, 'waterdepth'],
        'waterstate':  10,
        'soilstorepavedstate':  nan,
        'soilstorebldgsstate':  nan,
        'soilstoreevetrstate':  nan,
        'soilstoredectrstate':  nan,
        'soilstoregrassstate':  nan,
        'soilstorebsoilstate':  nan,
        'soilstorewaterstate': 0,
        'snowinitially':  int(nan),
        'snowwaterpavedstate':  nan,
        'snowwaterbldgsstate':  nan,
        'snowwaterevetrstate':  nan,
        'snowwaterdectrstate':  nan,
        'snowwatergrassstate':  nan,
        'snowwaterbsoilstate':  nan,
        'snowwaterwaterstate':  nan,
        'snowpackpaved':  nan,
        'snowpackbldgs':  nan,
        'snowpackevetr':  nan,
        'snowpackdectr':  nan,
        'snowpackgrass':  nan,
        'snowpackbsoil':  nan,
        'snowpackwater':  nan,
        'snowfracpaved':  nan,
        'snowfracbldgs':  nan,
        'snowfracevetr':  nan,
        'snowfracdectr':  nan,
        'snowfracgrass':  nan,
        'snowfracbsoil':  nan,
        'snowfracwater':  nan,
        'snowdenspaved':  nan,
        'snowdensbldgs':  nan,
        'snowdensevetr':  nan,
        'snowdensdectr':  nan,
        'snowdensgrass':  nan,
        'snowdensbsoil':  nan,
        'snowdenswater':  nan,
        'snowalb0':  nan,
        'qn1_av': 0,
        'qn1_s_av': 0,
        'dqndt': 0,
        'dqnsdt': 0,
    }

    # load basic model settings
    dict_ModConfig = load_SUEWS_dict_ModConfig(path_runcontrol)
    # path for SUEWS input tables:
    path_input = path_runcontrol.parent / dict_ModConfig['fileinputpath']
    df_gridSurfaceChar = load_SUEWS_SurfaceChar_df(path_input)
    # only use the first year of each grid as base for initial conditions
    df_init = df_gridSurfaceChar.groupby('Grid').aggregate(np.min)

    # incorporate numeric entries from dict_ModConfig
    list_var_dim_from_dict_ModConfig = [
        (var, 0, val) for var, val in dict_ModConfig.items()
        if isinstance(val, (float, int))
    ]
    # set values according to `list_var_dim_from_dict_ModConfig`
    # and generate the secondary dimensions
    for var, dim, val in list_var_dim_from_dict_ModConfig:
        ind_dim = [str(i) for i in np.ndindex(int(dim))] if dim > 0 else ['0']
        val_x = df_init[val] if isinstance(val, str) else val
        for ind in ind_dim:
            df_init[(var, str(ind))] = val_x

    # create an empty df
    # list_grid = df_gridSurfaceChar_init.index
    # df_InitialCond_grid = pd.DataFrame(
    #     np.ones([len(list_grid), len(dict_InitCond_default)]),
    #     index=list_grid,
    #     columns=list(dict_InitCond_default)
    # )
    # df_InitialCond_grid.index.set_names('Grid', inplace=True)

    # initialise df_InitialCond_grid with default values
    for k in dict_InitCond_default:
        # df_InitialCond_grid[k] = dict_InitCond_default[k]
        df_init[(k, '0')] = dict_InitCond_default[k]

    # update `temp_c0` with met_forcing
    # # TODO: add temp_c0 from met_forcing

    # update `waterstate` with df_gridSurfaceChar
    # drop_duplicates in case multi-year run is set
    # df_InitialCond_grid['waterstate'] = df_gridSurfaceChar_init['waterdepth']
    df_init[('waterstate', '0')] = df_init['waterdepth']

    # generate proper Initial Condition file names
    base_str = 'initialconditions' + dict_ModConfig['filecode']
    grid_str = (df_init.index.astype(str)
                if dict_ModConfig['multipleinitfiles'] == 1
                else '')
    year_str = df_init[('year', '0')].astype(int).astype(str)
    # ser_path_init =base_str + grid_str + '_' + year_str + '.nml'
    df_init[('file_init', '0')] = base_str + grid_str + '_' + year_str + '.nml'
    # ser_path_init = ser_path_init.map(lambda fn: path_input / fn)
    df_init[('file_init', '0')] = df_init[('file_init', '0')]\
        .map(lambda fn: path_input / fn)

    return df_init


# load Initial Condition variables from namelist file
def add_file_init_df(df_init):
    # load all nml info from file names:('file_init', '0')
    df_init_file = pd.concat(
        df_init[('file_init', '0')].map(
            lambda fn: load_SUEWS_nml(fn))
        .to_dict()).unstack()['initialconditions']
    df_init_file = pd.concat([df_init_file], axis=1, keys=['0'])
    df_init_file = df_init_file.swaplevel(0, 1, axis=1)
    df_init_file.index.set_names(['Grid'], inplace=True)

    # merge only those appeard in base df
    df_init.update(df_init_file)

    # drop ('file_init', '0') as this may cause non-numeic errors later on
    df_init = df_init.drop(columns=[('file_init', '0')])
    return df_init


# add veg info into `df_init`
def add_veg_init_df(df_init):
    # get veg surface fractions
    df_init[('fr_veg_sum', '0')] = df_init['sfr'].iloc[:, 2:5].sum(axis=1)

    # determine veg related parameters based on `leavesoutinitially`
    ser_leaves_on_flag = df_init['leavesoutinitially'] == 1

    # get gdd_0
    veg_sfr = df_init['sfr'].iloc[:, 2:5].T
    veg_sfr = veg_sfr.reset_index(drop=True).T
    x_gddfull = df_init['gddfull'].T.reset_index(drop=True).T
    df_init[('gdd_0', '0')] = (veg_sfr * x_gddfull).sum(axis=1)
    df_init[('gdd_0', '0')] /= veg_sfr.sum(axis=1)
    df_init[('const', '0')] = 0

    # list of tuples for assignment of leave on/off parameters
    list_var_veg = [
        # aligned for:
        # (var, ,leavesout==1,leavesout==0)
        ('gdd_1_0', ('gdd_0', '0'), ('const', '0')),
        ('gdd_2_0', ('const', '0'), ('gdd_0', '0')),
        ('laiinitialevetr', ('laimax', '(0,)'), ('laimin', '(0,)')),
        ('laiinitialdectr', ('laimax', '(1,)'), ('laimin', '(1,)')),
        ('laiinitialgrass', ('laimax', '(2,)'), ('laimin', '(2,)')),
        ('albevetr0', 'albmax_evetr', 'albmin_evetr'),
        ('albdectr0', 'albmax_dectr', 'albmin_dectr'),
        ('albgrass0', 'albmax_grass', 'albmin_grass'),
        ('decidcap0', 'capmax_dec', 'capmin_dec'),
        ('porosity0', 'pormin_dec', 'pormax_dec'),
    ]

    # set veg related parameters
    for var_veg, var_leaves_on, var_leaves_off in list_var_veg:
        # print(var_veg, var_leaves_on, var_leaves_off)
        # get default values based on ser_leaves_on_flag
        val_dft = df_init[var_leaves_on]\
            .where(ser_leaves_on_flag, df_init[var_leaves_off])
        # replace with val_dft if set as -999 (i.e., nan)
        ser_nan_flag = ~(df_init[var_veg] == -999)
        df_init[var_veg] = df_init[var_veg].where(
            ser_nan_flag, val_dft, axis=1)

    return df_init


# add surface specific info into `df_init`
def add_sfc_init_df(df_init):
    # create snow flag
    ser_snow_init = df_init['snowinitially']
    ser_snow_flag = ser_snow_init.where(ser_snow_init == 0, 1)
    ser_snow_flag = (df_init['snowuse'] * ser_snow_flag)['0']

    # land cover based assignment
    list_sfc = ['paved', 'bldgs', 'evetr', 'dectr', 'grass', 'bsoil', 'water']
    # dict {major variable: [minor variables]}
    dict_var_sfc = {
        var: ['{var}{sfc}'.format(var=var, sfc=sfc) for sfc in list_sfc]
        for var in ['snowdens', 'snowfrac', 'snowpack']
    }
    dict_var_sfc.update({
        'meltwaterstore': ['snowwater{}state'.format(sfc) for sfc in list_sfc],
        'soilmoist_id': ['soilstore{}state'.format(sfc) for sfc in list_sfc],
        'state_id': ['{}state'.format(sfc) for sfc in list_sfc],
    })

    # set values according to `dict_var_sfc`
    # and generate the secondary dimensions
    for var in dict_var_sfc:
        for ind, var_sfc in zip(np.ndindex(len(list_sfc)), dict_var_sfc[var]):
            ind_str = str(ind)
            df_init[(var, ind_str)] = df_init[var_sfc]
            if var_sfc.startswith('snow'):
                df_init[(var, ind_str)] *= ser_snow_flag.values.flatten()

    return df_init


# add initial daily state into `df_init`
def add_state_init_df(df_init):
    # dimensional assignment
    # (var, dim, `default value` or `reference name`)
    list_var_dim = [
        ('icefrac', 7, 0.2),
        ('snowfallcum', 0, 0.),
        ('snowalb', 0, 'snowalb0'),

        ('porosity_id', 0, 'porosity0'),
        ('albdectr_id', 0, 'albdectr0'),
        ('albevetr_id', 0, 'albevetr0'),
        ('albgrass_id', 0, 'albgrass0'),
        ('decidcap_id', 0, 'decidcap0'),

        ('lai_id', 3, 0.),
        ('gdd_id', 5, 0.),
        ('hdd_id', 12, 0.),
        ('wuday_id', 9, 0.),

        ('dt_since_start', 0, 0.),
        ('tstep_prev', 0, 'tstep'),
        ('tair24hr', int(24 * 3600 / df_init['tstep'].values[0, 0]), 273.15),
    ]

    # set values according to `list_var_dim`
    # and generate the secondary dimensions
    for var, dim, val in list_var_dim:
        ind_dim = [str(i) for i in np.ndindex(int(dim))] if dim > 0 else ['0']
        val_x = df_init[val] if isinstance(val, str) else val
        for ind in ind_dim:
            df_init[(var, str(ind))] = val_x

    # modifications for special cases
    # `gdd_id` corrections:
    df_init[('gdd_id', '(2,)')] = 90
    df_init[('gdd_id', '(3,)')] = -90
    # `numcapita` calculation:
    df_init[('numcapita', '0')] = \
        df_init[['popdensdaytime', 'popdensnighttime']]\
        .mean(axis=1)
    # `gridiv` addition:
    df_init[('gridiv', '0')] = df_init.index

    return df_init


# add additional parameters to `df` produced by `load_SUEWS_InitialCond_df`
@functools.lru_cache(maxsize=16)
def load_InitialCond_grid_df(path_runcontrol):
    # load base df of InitialCond
    df_init = load_SUEWS_InitialCond_df(path_runcontrol)

    # add Initial Condition variables from namelist file
    df_init = add_file_init_df(df_init)

    # add veg info into `df_init`
    df_init = add_veg_init_df(df_init)

    # add surface specific info into `df_init`
    df_init = add_sfc_init_df(df_init)

    # add initial daily state into `df_init`
    df_init = add_state_init_df(df_init)

    return df_init


##############################################################################
# archived:
# # load initial conditions as dict's for SUEWS running:
# # 1. dict_mod_cfg: model settings
# # 2. dict_state_init: initial model states/conditions
# # 3. dict_met_forcing: met forcing conditions (NB:NOT USED)
# # return dict
# def init_SUEWS_dict(path_runcontrol):
#     path_runcontrol = Path(path_runcontrol)
#
#     # initialise dict_state_init
#     dict_state_init = {}
#
#     # load RunControl variables
#     lib_RunControl = load_SUEWS_nml(path_runcontrol)
#     dict_RunControl = lib_RunControl.loc[:, 'runcontrol'].to_dict()
#
#     # # path for SUEWS input tables:
#     path_input = path_runcontrol.parent / dict_RunControl['fileinputpath']
#
#     # mod_config: static properties
#     dict_ModConfig = {'aerodynamicresistancemethod': 2,
#                       'evapmethod': 2,
#                       'laicalcyes': 1,
#                       'veg_type': 1,
#                       'diagnose': 0,
#                       'diagqn': 0,
#                       'diagqs': 0}
#     dict_ModConfig.update(dict_RunControl)
#
#     # dict for temporally varying states
#     # load surface charasteristics
#     df_gridSurfaceChar = load_SUEWS_SurfaceChar(path_input)
#     for grid in df_gridSurfaceChar.index:
#         # load two dict's for `grid`:
#         # 1. dict_StateInit: initial states
#         # 2. dict_MetForcing: met forcing (NB:NOT USED)
#         dict_StateInit = init_SUEWS_dict_grid(
#             path_input, grid, dict_ModConfig, df_gridSurfaceChar)
#
#         # return dict_StateInit
#         # dict with all properties of one grid:
#         # 1. model settings: dict_ModConfig
#         # 2. surface properties
#         # combine them as dict_grid
#         dict_grid = df_gridSurfaceChar.loc[grid, :].to_dict()
#         dict_grid.update(dict_ModConfig)
#         dict_grid.update(dict_StateInit)
#         dict_grid.update(gridiv=grid)
#         # filter out unnecessary entries for main calculation
#         dict_state_init_grid = {
#             k: dict_grid[k] for k in list(
#                 set(dict_grid.keys()).intersection(
#                     set(get_args_suews()['var_input'])))}
#
#         # kepp other entries as model configuration
#         dict_mod_cfg = {
#             k: dict_grid[k] for k in list(
#                 set(dict_grid.keys()) - set(get_args_suews()['var_input']))}
#
#         # construct a dict with entries as:
#         # {grid: dict_state_init}
#         dict_state_init.update({grid: dict_state_init_grid})
#         # {grid: dict_MetForcing}
#         # dict_met_forcing.update({grid: dict_MetForcing})
#     # end grid loop
#
#     return dict_mod_cfg, dict_state_init
#
#
# # create initial states for one grid
# def init_SUEWS_dict_grid(path_input, grid,
#       dict_ModConfig, df_gridSurfaceChar):
#     # load tstep from dict_RunControl
#     tstep = dict_ModConfig['tstep']
#     # some constant values
#     nan = -999.
#     # ndays = 366
#     nsh = int(3600 / tstep)  # tstep from dict_RunControl
#
#     # load met forcing of `grid`:
#     filecode = dict_ModConfig['filecode']
#     tstep_in = dict_ModConfig['resolutionfilesin']
#     multiplemetfiles = dict_ModConfig['multiplemetfiles']
#     # load as DataFrame:
#     df_forcing_met = load_SUEWS_Forcing_met_df_raw(
#         path_input, filecode, grid, tstep_in, multiplemetfiles)
#
#     # define some met forcing determined variables:
#     # previous day index
#     # id_prev = int(df_forcing_met['id'].iloc[0] - 1)
#
#     # initialise dict_InitCond with default values
#     dict_InitCond = {
#         'dayssincerain':  0,
#         # `temp_c0` defaults to daily mean air temperature of the first day
#         'temp_c0':  df_forcing_met['temp_c'].iloc[:(24 * nsh) - 1].mean(),
#         'leavesoutinitially':  int(nan),
#         'gdd_1_0':  nan,
#         'gdd_2_0':  nan,
#         'laiinitialevetr':  nan,
#         'laiinitialdectr':  nan,
#         'laiinitialgrass':  nan,
#         'albevetr0':  nan,
#         'albdectr0':  nan,
#         'albgrass0':  nan,
#         'decidcap0':  nan,
#         'porosity0':  nan,
#         'pavedstate':  0,
#         'bldgsstate':  0,
#         'evetrstate':  0,
#         'dectrstate':  0,
#         'grassstate':  0,
#         'bsoilstate':  0,
#         'waterstate':  df_gridSurfaceChar.at[grid, 'waterdepth'],
#         'soilstorepavedstate':  nan,
#         'soilstorebldgsstate':  nan,
#         'soilstoreevetrstate':  nan,
#         'soilstoredectrstate':  nan,
#         'soilstoregrassstate':  nan,
#         'soilstorebsoilstate':  nan,
#         'snowinitially':  int(nan),
#         'snowwaterpavedstate':  nan,
#         'snowwaterbldgsstate':  nan,
#         'snowwaterevetrstate':  nan,
#         'snowwaterdectrstate':  nan,
#         'snowwatergrassstate':  nan,
#         'snowwaterbsoilstate':  nan,
#         'snowwaterwaterstate':  nan,
#         'snowpackpaved':  nan,
#         'snowpackbldgs':  nan,
#         'snowpackevetr':  nan,
#         'snowpackdectr':  nan,
#         'snowpackgrass':  nan,
#         'snowpackbsoil':  nan,
#         'snowpackwater':  nan,
#         'snowfracpaved':  nan,
#         'snowfracbldgs':  nan,
#         'snowfracevetr':  nan,
#         'snowfracdectr':  nan,
#         'snowfracgrass':  nan,
#         'snowfracbsoil':  nan,
#         'snowfracwater':  nan,
#         'snowdenspaved':  nan,
#         'snowdensbldgs':  nan,
#         'snowdensevetr':  nan,
#         'snowdensdectr':  nan,
#         'snowdensgrass':  nan,
#         'snowdensbsoil':  nan,
#         'snowdenswater':  nan,
#         'snowalb0':  nan
#     }
#
#     # load Initial Condition variables from namelist file
#     InitialCond_fn = 'initialconditions{site}{grid}_{year}.nml'.format(
#         site=dict_ModConfig['filecode'],
#         # grid info will be include in the nml file pattern
#         # if multiple init files are used
#         grid=(grid if dict_ModConfig['multipleinitfiles'] == 1 else ''),
#         year=int(np.min(df_gridSurfaceChar.loc[grid, 'year']))
#     )
#     lib_InitCond = load_SUEWS_nml(path_input / InitialCond_fn)
#     # update default InitialCond with values set in namelist
#     dict_InitCond.update(lib_InitCond['initialconditions'].to_dict())
#
#     # fr_veg_sum: total fraction of vegetation covers
#     fr_veg_sum = np.sum(df_gridSurfaceChar.loc[grid, 'sfr'][2:5])
#     # update vegetation-related variables according to LeavesOutInitially
#     if dict_InitCond['leavesoutinitially'] == 1:
#         dict_InitCond['gdd_1_0'] = (np.dot(
#             df_gridSurfaceChar.loc[grid, 'gddfull'],
#             df_gridSurfaceChar.loc[grid, 'sfr'][2:5]) / fr_veg_sum
#             if fr_veg_sum > 0
#             else 0)
#         dict_InitCond['gdd_2_0'] = 0
#         dict_InitCond['laiinitialevetr'] = df_gridSurfaceChar.loc[
#             grid, 'laimax'][0]
#         dict_InitCond['laiinitialdectr'] = df_gridSurfaceChar.loc[
#             grid, 'laimax'][1]
#         dict_InitCond['laiinitialgrass'] = df_gridSurfaceChar.loc[
#             grid, 'laimax'][2]
#         dict_InitCond['albevetr0'] = df_gridSurfaceChar.loc[
#             grid, 'albmax_evetr']
#         dict_InitCond['albdectr0'] = df_gridSurfaceChar.loc[
#             grid, 'albmax_dectr']
#         dict_InitCond['albgrass0'] = df_gridSurfaceChar.loc[
#             grid, 'albmax_grass']
#         dict_InitCond['decidcap0'] = df_gridSurfaceChar.loc[
#             grid, 'capmax_dec']
#         dict_InitCond['porosity0'] = df_gridSurfaceChar.loc[
#             grid, 'pormin_dec']
#
#     elif dict_InitCond['leavesoutinitially'] == 0:
#         dict_InitCond['gdd_1_0'] = 0
#         dict_InitCond['gdd_2_0'] = (np.dot(
#             df_gridSurfaceChar.loc[grid, 'sddfull'],
#             df_gridSurfaceChar.loc[grid, 'sfr'][2: 5]) / fr_veg_sum
#             if fr_veg_sum > 0
#             else 0)
#         dict_InitCond['laiinitialevetr'] = df_gridSurfaceChar.loc[
#             grid, 'laimin'][0]
#         dict_InitCond['laiinitialdectr'] = df_gridSurfaceChar.loc[
#             grid, 'laimin'][1]
#         dict_InitCond['laiinitialgrass'] = df_gridSurfaceChar.loc[
#             grid, 'laimin'][2]
#         dict_InitCond['albevetr0'] = df_gridSurfaceChar.loc[
#             grid, 'albmin_evetr']
#         dict_InitCond['albdectr0'] = df_gridSurfaceChar.loc[
#             grid, 'albmin_dectr']
#         dict_InitCond['albgrass0'] = df_gridSurfaceChar.loc[
#             grid, 'albmin_grass']
#         dict_InitCond['decidcap0'] = df_gridSurfaceChar.loc[
#             grid, 'capmin_dec']
#         dict_InitCond['porosity0'] = df_gridSurfaceChar.loc[
#             grid, 'pormax_dec']
#     # check when dict_InitCond['leavesoutinitially'] == nan
#     elif dict_InitCond['leavesoutinitially'] == nan:
#         list_var_veg = [
#             'gdd_1_0',
#             'gdd_2_0',
#             'laiinitialevetr',
#             'laiinitialdectr',
#             'laiinitialgrass',
#             'albevetr0',
#             'albdectr0',
#             'albgrass0',
#             'decidcap0',
#             'porosity0',
#         ]
#         for k in list_var_veg:
#             if dict_InitCond[k] == nan:
#                 print(
#                     '{var} is invalid for grid {grid}'
#                     .format(var=k, grid=grid)
#                 )
#   # END: update vegetation-related variables according to LeavesOutInitially
#
#     # snowflag = 1 if snow-related modules are enabled or 0 otherwise
#     snowflag = (
#         0 if (
#             dict_ModConfig['snowuse'] == 0
#             or dict_InitCond['snowinitially'] == 0
#         )
#         else 1)
#
#     # state_init: temporally-varying states from timestep to timestep
#     dict_StateInit = {
#         # TODO: water use patterns, which is currently not used
#         'wuday_id': 0. * np.ones(9, order='F'),
#         'numcapita': (
#             df_gridSurfaceChar.at[grid, 'popdensdaytime'] +\
#             df_gridSurfaceChar.at[grid, 'popdensnighttime']
#         ) * 0.5,
#         'qn1_av': 0 * np.ones(1),
#         'qn1_s_av': 0 * np.ones(1),
#         'dqndt': 0 * np.ones(1),
#         'dqnsdt': 0 * np.ones(1),
#         # 'qn1_s_store_grid': nan * np.ones(nsh),
#         # 'qn1_store_grid': nan * np.ones(nsh),
#
#         # snow:
#         'snowalb': dict_InitCond['snowalb0'],
#         'snowfallcum': 0.,
#         'icefrac': 0.2 * np.ones(7, order='F'),
#         # Initial liquid (melted) water for each surface
#         'meltwaterstore': snowflag * np.array(
#             [
#                 dict_InitCond[var] for var in
#                 [
#                     'snowwaterpavedstate',
#                     'snowwaterbldgsstate',
#                     'snowwaterevetrstate',
#                     'snowwaterdectrstate',
#                     'snowwatergrassstate',
#                     'snowwaterbsoilstate',
#                     'snowwaterwaterstate'
#                 ]
#             ],
#             order='F'),
#         'snowdens': snowflag * np.array(
#             [
#                 dict_InitCond[var] for var in
#                 [
#                     'snowdenspaved',
#                     'snowdensbldgs',
#                     'snowdensevetr',
#                     'snowdensdectr',
#                     'snowdensgrass',
#                     'snowdensbsoil',
#                     'snowdenswater'
#                 ]
#             ],
#             order='F'),
#         'snowfrac': snowflag * np.array(
#             [
#                 dict_InitCond[var] for var in
#                 [
#                     'snowfracpaved',
#                     'snowfracbldgs',
#                     'snowfracevetr',
#                     'snowfracdectr',
#                     'snowfracgrass',
#                     'snowfracbsoil',
#                     'snowfracwater'
#                 ]
#             ],
#             order='F'),
#         'snowpack': snowflag * np.array(
#             [
#                 dict_InitCond[var] for var in
#                 [
#                     'snowpackpaved',
#                     'snowpackbldgs',
#                     'snowpackevetr',
#                     'snowpackdectr',
#                     'snowpackgrass',
#                     'snowpackbsoil',
#                     'snowpackwater'
#                 ]
#             ],
#             order='F'),
#
#         # Initial soil stores for each surface (below ground)
#         'soilmoist_id': np.array(
#             [
#                 dict_InitCond[var] for var in [
#                     'soilstorepavedstate',
#                     'soilstorebldgsstate',
#                     'soilstoreevetrstate',
#                     'soilstoredectrstate',
#                     'soilstoregrassstate',
#                     'soilstorebsoilstate'
#                 ]
#             ]
#             + [0.], order='F'),
#
#         # Initial wetness status of each surface (above ground)
#         'state_id': np.array(
#             [
#                 dict_InitCond[var]
#                 for var in [
#                     'pavedstate',
#                     'bldgsstate',
#                     'evetrstate',
#                     'dectrstate',
#                     'grassstate',
#                     'bsoilstate',
#                     'waterstate'
#                 ]
#             ],
#             order='F'),
#
#         # mean Tair of past 24 hours, NOT used by supy as it is for ESTM
#         'tair24hr': 273.15 * np.ones(24 * nsh),
#
#
#         # vegetation related parameters:
#         'porosity_id': dict_InitCond['porosity0'] * np.ones(1, order='F'),
#         'albdectr_id': dict_InitCond['albdectr0'] * np.ones(1, order='F'),
#         'albevetr_id': dict_InitCond['albevetr0'] * np.ones(1, order='F'),
#         'albgrass_id': dict_InitCond['albgrass0'] * np.ones(1, order='F'),
#         'decidcap_id': dict_InitCond['decidcap0'] * np.ones(1, order='F'),
#         # leaf area index:
#         'lai_id': 1. * np.zeros(3),
#         # growing degree days:
#         'gdd_id': np.zeros(5),
#         # heating degree days:
#         'hdd_id': np.zeros(12)
#     }
#
#     # gdd-related parameters:
#     dict_StateInit['gdd_id'][2] = 90
#     dict_StateInit['gdd_id'][3] = -90
#
#     # update timestep info:
#     dict_StateInit.update(tstep_prev=tstep)
#     dict_StateInit.update(dt_since_start=0)
#
#     # pack final results with Fortran ordering
#     dict_StateInit = {k: np.array(v, order='F')
#                       for k, v in dict_StateInit.items()}
#
#     return dict_StateInit
#


# https://stackoverflow.com/a/44776960/920789
# class HDict(dict):
#     def __hash__(self):
#         # if not isinstance(self.items(), collections.Hashable):
#         #     print(self.items(), hash(ujson.dumps(self)))
#         return hash(ujson.dumps(self))
#
#
# class HList(list):
#     def __hash__(self):
#         list_proc = [HDict(x) if isinstance(x, dict)
#                      else x for x in self]
#         return hash(tuple(list_proc))


# def hash_dict(func):
#     """Transform mutable dictionnary
#     Into immutable
#     Useful to be compatible with cache
#     """
#
#     @functools.wraps(func)
#     def wrapped(*args, **kwargs):
#         # print('\n print args before')
#         # print(args)
#         args = tuple([HDict(arg) if isinstance(
#             arg, dict) else arg for arg in args])
#         args = tuple([HList(arg) if isinstance(
#             arg, list) else arg for arg in args])
#
#         kwargs = {k: HDict(v) if isinstance(v, dict)
#                   else v for k, v in kwargs.items()}
#         kwargs = {k: HList(v) if isinstance(v, list)
#                   else v for k, v in kwargs.items()}
#
#         # print('\n print args after')
#         # print(args)
#         # print(kwargs, '\n')
#         return func(*args, **kwargs)
#     return wrapped
