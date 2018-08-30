###########################################################################
# SUEWS for Python
# Authors:
# Ting Sun, ting.sun@reading.ac.uk
# History:
# 20 Jan 2018: first alpha release
# 01 Feb 2018: performance improvement
# 03 Feb 2018: improvement in output processing
# 08 Mar 2018: pypi packaging
###########################################################################

# load dependency modules
import os
import sys
import numpy as np
import pandas as pd
import f90nml
from pandas import DataFrame as df
import inspect
from scipy import interpolate
import copy
import glob
from datetime import timedelta

# import math
# import random
# import pkg_resources

# define local path for loading resources in this package
dir_path = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda: 0)))
# sys.path.append(dir_path)
# load f2py-based SUEWS calculation core
from .SUEWS_driver import suews_driver as sd


########################################################################
# get_args_suews can get the interface informaiton
# of the f2py-converted Fortran interface
def get_args_suews():
    # split doc lines for processing
    docLines = np.array(
        sd.suews_cal_main.__doc__.splitlines(),
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
list_file_input = ['SUEWS_AnthropogenicHeat.txt',
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
                   'SUEWS_WithinGridWaterDist.txt']

# library of all properties
dict_libVar2File = {fileX.replace('.txt', '').replace(
    'SUEWS', 'lib'): fileX for fileX in list_file_input
    if fileX.endswith('.txt')}

# dictionaries:
# links between code in SiteSelect to properties in according tables
# this is described in SUEWS online manual:
# http://urban-climate.net/umep/SUEWS#SUEWS_SiteSelect.txt
path_code2file = os.path.join(dir_path, 'code2file.json')
dict_Code2File = pd.read_json(path_code2file, typ='series').to_dict()
# variable translation as done in Fortran-SUEWS
path_var2siteselect = os.path.join(dir_path, 'var2siteselect.json')
dict_var2SiteSelect = pd.read_json(path_var2siteselect, typ='series').to_dict()

# expand dict_Code2File for retrieving surface characteristics
dict_varSiteSelect2File = {
    x: 'SUEWS_SiteSelect.txt' for x in dict_var2SiteSelect.keys()}
dict_Code2File.update(dict_varSiteSelect2File)


# load model settings
# load configurations: mod_config
# process RunControl.nml
# this function can handle all SUEWS nml files
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
                            comment='!', error_bad_lines=True,
                            skiprows=1, index_col=0).dropna()
    return rawdata


# load all tables into variables staring with 'lib_' and filename
def load_SUEWS_Libs(dir_input):
    dict_libs = {}
    for lib, lib_file in dict_libVar2File.items():
        lib_path = os.path.join(dir_input, lib_file)
        dict_libs.update({lib: load_SUEWS_table(lib_path)})
    # return DataFrame containing settings
    return dict_libs


# look up properties according to code
def lookup_code_lib(libName, codeKey, codeValue, dict_libs):
    str_lib = dict_Code2File[libName].replace(
        '.txt', '').replace('SUEWS', 'lib')
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
def lookup_KeySeq_lib(indexKey, subKey, indexCode, dict_libs):
    # print indexKey, subKey, indexCode
    if type(subKey) is float:
        res = subKey
    elif type(subKey) is str:
        res = lookup_code_lib(indexKey, subKey, indexCode, dict_libs)
    elif type(subKey) is dict:
        indexKeyX, subKeyX = list(subKey.items())[0]
        indexCodeX = lookup_code_lib(indexKey, indexKeyX, indexCode, dict_libs)
        res = lookup_KeySeq_lib(indexKeyX, subKeyX, indexCodeX, dict_libs)
    elif type(subKey) is list:
        res = []
        for subKeyX in subKey:
            indexCodeX = indexCode
            resX = lookup_KeySeq_lib(indexKey, subKeyX, indexCodeX, dict_libs)
            res.append(resX)
    # final result
    return res


# load surface characteristics
def load_SUEWS_SurfaceChar(dir_input):
    # load RunControl variables
    lib_RunControl = load_SUEWS_nml(os.path.join(dir_input, 'runcontrol.nml'))
    dict_RunControl = lib_RunControl.loc[:, 'runcontrol'].to_dict()
    tstep = dict_RunControl['tstep']
    dir_path = os.path.join(dir_input, dict_RunControl['fileinputpath'])
    # load all libraries
    dict_libs = load_SUEWS_Libs(dir_path)
    # construct a dictionary in the form: {grid:{var:value,...}}
    dict_gridSurfaceChar = {
        grid: {k: lookup_KeySeq_lib(k, v, grid, dict_libs)
               for k, v in dict_var2SiteSelect.items()}
        for grid in dict_libs['lib_SiteSelect'].index}
    # convert the above dict to DataFrame
    df_gridSurfaceChar = pd.DataFrame.from_dict(dict_gridSurfaceChar).T
    # empty dict to hold updated values
    dict_x_grid = {}
    # modify some variables to be compliant with SUEWS requirement
    for xgrid in df_gridSurfaceChar.index:
        # transpoe snowprof:
        df_gridSurfaceChar.loc[xgrid, 'snowprof_24hr'] = np.array(
            df_gridSurfaceChar.loc[xgrid, 'snowprof_24hr'], order='F').T

        # transpoe laipower:
        df_gridSurfaceChar.loc[xgrid, 'laipower'] = np.array(
            df_gridSurfaceChar.loc[xgrid, 'laipower'], order='F').T
        # print df_gridSurfaceChar.loc[xgrid, 'laipower'].shape

        # select non-zero values for waterdist of water surface:
        x = np.array(df_gridSurfaceChar.loc[xgrid, 'waterdist'][-1])
        df_gridSurfaceChar.loc[xgrid, 'waterdist'][-1] = (
            x[np.nonzero(x)])

        # surf order as F:
        df_gridSurfaceChar.loc[xgrid, 'storedrainprm'] = np.array(
            df_gridSurfaceChar.loc[xgrid, 'storedrainprm'], order='F')

        # convert to np.array
        df_gridSurfaceChar.loc[xgrid, 'alb'] = np.array(
            df_gridSurfaceChar.loc[xgrid, 'alb'])

        # convert unit of `surfacearea` from ha to m^2
        df_gridSurfaceChar.loc[xgrid, 'surfacearea'] = np.array(
            df_gridSurfaceChar.loc[xgrid, 'surfacearea'] * 10000.)

        # print type(df_gridSurfaceChar.loc[xgrid, 'alb'])

        # dict holding updated values that can be converted to DataFrame later
        dict_x = df_gridSurfaceChar.loc[xgrid, :].to_dict()
        # print 'len(dict_x)',len(dict_x['laipower'])

        # profiles:
        # t_tstep = np.linspace(0, 24, num=3600 / tstep * 24, endpoint=False)
        list_varTstep = ['ahprof_24hr',
                         'popprof_24hr',
                         'traffprof_24hr',
                         'humactivity_24hr',
                         'wuprofm_24hr',
                         'wuprofa_24hr']
        for var in list_varTstep:
            #     var_name = var.replace('_24hr', '')
            var_T = np.array(df_gridSurfaceChar.loc[xgrid, var]).T
        #     var0 = np.vstack((var0, var0))
        #     # interpolator:
        #     f = interpolate.interp1d(np.arange(0, 48), var0, axis=0)
            dict_x.update({var: var_T})
        #     # different normalisation processing:
        #     # For water use, SUM of the multipliers is equal to 1
        #     if var in ['wuprofm_24hr', 'wuprofa_24hr']:
        #         dict_x[var] = dict_x[var] / np.sum(dict_x[var], axis=0)
        #     # For traffic, AVERAGE of the multipliers is equal to 1
        # elif var in ['ahprof_24hr', 'traffprof_24hr']:
        #         dict_x[var] = (dict_x[var] / np.sum(dict_x[var], axis=0)
        #                        * len(dict_x[var]))

        # update dict to hold grids
        dict_x_grid.update({xgrid: dict_x})

    # convert to DataFrame
    df_x_grid = pd.DataFrame.from_dict(dict_x_grid).T
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
    dectime = (t.dayofyear - 1) + (t.hour +
                                   (t.minute + t.second / 60.) / 60.) / 24
    return dectime

# resample solar radiation by zenith correction and total amount distribution


def resample_kdn(data_raw_kdn, tstep_mod, timezone, lat, lon, alt):
    # adjust solar radiation
    datetime_mid_local = data_raw_kdn.index - timedelta(
        seconds=tstep_mod / 2)
    sol_elev = np.array([sd.cal_sunposition(
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
    data_raw.index = data_raw.loc[:, ['iy', 'id', 'it', 'imin']].apply(
        func_parse_date_row, 1)
    # shift by half-tstep_in to generate a time series with instantaneous
    # values
    data_raw_shift = data_raw.copy().shift(-tstep_in / 2, freq='S')

    # downscale input data to desired time step
    data_raw_tstep = data_raw_shift.resample(
        '{tstep}S'.format(tstep=tstep_mod)).interpolate(
        method='polynomial', order=1).rolling(
        window=2, center=False).mean()

    # reindex data_tstep to valid range
    ix = pd.date_range(
        data_raw.index[0] - timedelta(seconds=tstep_in - tstep_mod),
        data_raw.index[-1],
        freq='{tstep}S'.format(tstep=tstep_mod))
    data_tstep = data_raw_tstep.copy().reindex(
        index=ix).bfill().ffill().dropna()

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

    # reset index with numbers
    data_met_tstep_out = data_met_tstep.copy().reset_index(drop=True)

    return data_met_tstep_out


# load raw data: met forcing
def load_SUEWS_Forcing_met_df_raw(
        dir_input, filecode, grid, tstep_met_in, multiplemetfiles):
    # file name pattern for met files
    forcingfile_met_pattern = os.path.join(
        dir_input,
        '{site}{grid}*{tstep}*txt'.format(
            site=filecode,
            grid=(grid if multiplemetfiles == 1 else ''),
            tstep=int(tstep_met_in / 60)))

    # list of met forcing files
    list_file_MetForcing = sorted([
        f for f in glob.glob(forcingfile_met_pattern)
        if 'ESTM' not in f])

    # load raw data
    df_forcing_met = pd.concat(
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
        # rename these columns to match variables via the driver interface
        columns={
            '%' + 'iy': 'iy',
            'id': 'id',
            'it': 'it',
            'imin': 'imin',
            'Kdn': 'avkdn',
            'RH': 'avrh',
            'Wind': 'avu1',
            'fcld': 'fcld_obs',
            'lai_hr': 'lai_obs',
            'ldown': 'ldown_obs',
            'rain': 'precip',
            'press': 'press_hpa',
            'Qs': 'qs_obs',
            'Qf': 'qf_obs',
            'QH': 'qh_obs',
            'Q*': 'qn1_obs',
            'snow': 'snow_obs',
            'Td': 'temp_c',
            # 'all': 'metforcingdata_grid',
            'xsmd': 'xsmd'})

    # convert unit from kPa to hPa
    df_forcing_met['press_hpa'] *= 10

    # add `isec` for WRF-SUEWS interface
    df_forcing_met['isec'] = 0

    # set correct data types
    df_forcing_met[['iy', 'id', 'it', 'imin', 'isec']] = df_forcing_met[[
        'iy', 'id', 'it', 'imin', 'isec']].astype(np.int64)

    return df_forcing_met


# load raw data: met forcing
def load_SUEWS_Forcing_ESTM_df_raw(
        dir_input, filecode, grid, tstep_ESTM_in, multipleestmfiles):
    # file name pattern for met files
    forcingfile_ESTM_pattern = os.path.join(
        dir_input,
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
            'Kdn': 'avkdn',
            'RH': 'avrh',
            'Wind': 'avu1',
            'fcld': 'fcld_obs',
            'lai_hr': 'lai_obs',
            'ldown': 'ldown_obs',
            'rain': 'precip',
            'press': 'press_hpa',
            'QH': 'qh_obs',
            'Q*': 'qn1_obs',
            'snow': 'snow_obs',
            'Td': 'temp_c',
            # 'all': 'metforcingdata_grid',
            'xsmd': 'xsmd'})

    # set correct data types
    df_forcing_estm[['iy', 'id', 'it', 'imin']] = df_forcing_estm[[
        'iy', 'id', 'it', 'imin']].astype(np.int64)

    return df_forcing_estm

# TODO: add support for loading multi-grid forcing datasets


def load_SUEWS_Forcing_df(dir_site, ser_mod_cfg, df_state_init):
    pass


# load forcing datasets of `grid`
def load_SUEWS_Forcing_df_grid(dir_site, grid, ser_mod_cfg, df_state_init):
    # load setting variables from ser_mod_cfg
    (filecode, kdownzen,
     tstep_met_in, tstep_ESTM_in,
     multiplemetfiles, multipleestmfiles,
     dir_input_cfg) = ser_mod_cfg[
        ['filecode', 'kdownzen',
         'resolutionfilesin', 'resolutionfilesinestm',
         'multiplemetfiles', 'multipleestmfiles', 'fileinputpath']]
    tstep_mod, lat, lon, alt, timezone = df_state_init.loc[
        grid,
        ['tstep', 'lat', 'lng', 'alt', 'timezone']]
    dir_input = (os.path.join(dir_site, dir_input_cfg)
                 if not os.path.isabs(dir_input_cfg)
                 else dir_input_cfg)

    # load raw data
    # met forcing
    df_forcing_met = load_SUEWS_Forcing_met_df_raw(
        dir_input, filecode, grid, tstep_met_in, multiplemetfiles)

    # resample raw data from tstep_in to tstep_mod
    df_forcing_met_tstep = resample_forcing_met(
        df_forcing_met, tstep_met_in, tstep_mod, lat, lon, alt, timezone, kdownzen)

    # merge forcing datasets (met and ESTM)
    df_forcing_tstep = df_forcing_met_tstep.copy()
    # df_forcing_tstep = df_forcing_met_tstep.merge(
    #     df_forcing_estm_tstep,
    #     left_on=['iy', 'id', 'it', 'imin'],
    #     right_on=['iy', 'id', 'it', 'imin'])

    # pack all records of `id` into `metforcingdata_grid` for AnOHM and others
    df_grp = df_forcing_tstep.groupby('id')
    dict_id_all = {xid: df_grp.get_group(xid)
                   for xid in df_forcing_tstep['id'].unique()}
    id_all = df_forcing_tstep['id'].apply(lambda xid: dict_id_all[xid])
    df_forcing_tstep = df_forcing_tstep.merge(id_all.to_frame(name='metforcingdata_grid'),
                                              left_index=True,
                                              right_index=True)

    # add Ts forcing for ESTM
    if df_state_init.iloc[0]['storageheatmethod'] == 4:
        # load ESTM forcing
        df_forcing_estm = load_SUEWS_Forcing_ESTM_df_raw(
            dir_input, filecode, grid, tstep_ESTM_in, multipleestmfiles)
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
    df_forcing_tstep['dectime'] = ((df_forcing_tstep['id'] - 1) +
                                   (df_forcing_tstep['it']
                                    + df_forcing_tstep['imin'] / 60.) / 24.)
    df_forcing_tstep['id_prev_t'] = df_forcing_tstep['id'].shift(
        1).fillna(method='backfill')
    df_forcing_tstep['iy_prev_t'] = df_forcing_tstep['iy'].shift(
        1).fillna(method='backfill')

    return df_forcing_tstep


# load initial conditions as dict's for SUEWS running:
# 1. dict_mod_cfg: model settings
# 2. dict_state_init: initial model states/conditions
# 3. dict_met_forcing: met forcing conditions (NB:NOT USED)
def init_SUEWS_dict(dir_start):  # return dict
    # initialise dict_state_init
    dict_state_init = {}
    # dict_met_forcing = {}

    # load RunControl variables
    lib_RunControl = load_SUEWS_nml(os.path.join(dir_start, 'runcontrol.nml'))
    dict_RunControl = lib_RunControl.loc[:, 'runcontrol'].to_dict()

    # # path for SUEWS input tables:
    dir_input = os.path.join(dir_start,
                             dict_RunControl['fileinputpath'])

    # mod_config: static properties
    dict_ModConfig = {'aerodynamicresistancemethod': 2,
                      'evapmethod': 2,
                      'laicalcyes': 1,
                      'veg_type': 1,
                      'diagnose': 0,
                      'diagqn': 0,
                      'diagqs': 0}
    dict_ModConfig.update(dict_RunControl)

    # dict for temporally varying states
    # load surface charasteristics
    df_gridSurfaceChar = load_SUEWS_SurfaceChar(dir_start)
    for grid in df_gridSurfaceChar.index:
        # load two dict's for `grid`:
        # 1. dict_StateInit: initial states
        # 2. dict_MetForcing: met forcing (NB:NOT USED)
        dict_StateInit = init_SUEWS_dict_grid(
            dir_input, grid, dict_ModConfig, df_gridSurfaceChar)

        # return dict_StateInit
        # dict with all properties of one grid:
        # 1. model settings: dict_ModConfig
        # 2. surface properties
        # combine them as dict_grid
        dict_grid = df_gridSurfaceChar.loc[grid, :].to_dict()
        dict_grid.update(dict_ModConfig)
        dict_grid.update(dict_StateInit)
        dict_grid.update(gridiv=grid)
        # filter out unnecessary entries for main calculation
        dict_state_init_grid = {
            k: dict_grid[k] for k in list(
                set(dict_grid.keys()).intersection(
                    set(get_args_suews()['var_input'])))}

        # kepp other entries as model configuration
        dict_mod_cfg = {
            k: dict_grid[k] for k in list(
                set(dict_grid.keys()) - set(get_args_suews()['var_input']))}

        # construct a dict with entries as:
        # {grid: dict_state_init}
        dict_state_init.update({grid: dict_state_init_grid})
        # {grid: dict_MetForcing}
        # dict_met_forcing.update({grid: dict_MetForcing})
    # end grid loop

    return dict_mod_cfg, dict_state_init


# create initial states for one grid
def init_SUEWS_dict_grid(dir_input, grid,
                         dict_ModConfig, df_gridSurfaceChar):
    # load tstep from dict_RunControl
    tstep = dict_ModConfig['tstep']
    # some constant values
    nan = -999.
    ndays = 366
    nsh = int(3600 / tstep)  # tstep from dict_RunControl

    # load met forcing of `grid`:
    filecode = dict_ModConfig['filecode']
    tstep_in = dict_ModConfig['resolutionfilesin']
    multiplemetfiles = dict_ModConfig['multiplemetfiles']
    # load as DataFrame:
    df_forcing_met = load_SUEWS_Forcing_met_df_raw(
        dir_input, filecode, grid, tstep_in, multiplemetfiles)

    # define some met forcing determined variables:
    # previous day index
    id_prev = int(df_forcing_met['id'].iloc[0] - 1)

    # initialise dict_InitCond with default values
    dict_InitCond = {
        'dayssincerain':  0,
        # `temp_c0` defaults to daily mean air temperature of the first day
        'temp_c0':  df_forcing_met.loc[:(24 * nsh) - 1, 'temp_c'].mean(),
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
        'waterstate':  df_gridSurfaceChar.loc[grid, 'waterdepth'],
        'soilstorepavedstate':  nan,
        'soilstorebldgsstate':  nan,
        'soilstoreevetrstate':  nan,
        'soilstoredectrstate':  nan,
        'soilstoregrassstate':  nan,
        'soilstorebsoilstate':  nan,
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
        'snowalb0':  nan
    }

    # load Initial Condition variables from namelist file
    lib_InitCond = load_SUEWS_nml(os.path.join(
        dir_input, 'initialconditions{site}{grid}_{year}.nml'.format(
            site=dict_ModConfig['filecode'],
            # grid info will be include in the nml file pattern
            # if multiple init files are used
            grid=(grid if dict_ModConfig['multipleinitfiles'] == 1 else ''),
            year=int(np.min(df_gridSurfaceChar.loc[grid, 'year'])))))
    # update default InitialCond with values set in namelist
    dict_InitCond.update(
        lib_InitCond.loc[:, 'initialconditions'].to_dict())

    # fr_veg_sum: total fraction of vegetation covers
    fr_veg_sum = np.sum(df_gridSurfaceChar.loc[grid, 'sfr'][2:5])
    # update vegetation-related variables according to LeavesOutInitially
    if dict_InitCond['leavesoutinitially'] == 1:
        dict_InitCond['gdd_1_0'] = (np.dot(
            df_gridSurfaceChar.loc[grid, 'gddfull'],
            df_gridSurfaceChar.loc[grid, 'sfr'][2:5]) / fr_veg_sum
            if fr_veg_sum > 0
            else 0)
        dict_InitCond['gdd_2_0'] = 0
        dict_InitCond['laiinitialevetr'] = df_gridSurfaceChar.loc[
            grid, 'laimax'][0]
        dict_InitCond['laiinitialdectr'] = df_gridSurfaceChar.loc[
            grid, 'laimax'][1]
        dict_InitCond['laiinitialgrass'] = df_gridSurfaceChar.loc[
            grid, 'laimax'][2]
        dict_InitCond['albevetr0'] = df_gridSurfaceChar.loc[
            grid, 'albmax_evetr']
        dict_InitCond['albdectr0'] = df_gridSurfaceChar.loc[
            grid, 'albmax_dectr']
        dict_InitCond['albgrass0'] = df_gridSurfaceChar.loc[
            grid, 'albmax_grass']
        dict_InitCond['decidcap0'] = df_gridSurfaceChar.loc[
            grid, 'capmax_dec']
        dict_InitCond['porosity0'] = df_gridSurfaceChar.loc[
            grid, 'pormin_dec']

    elif dict_InitCond['leavesoutinitially'] == 0:
        dict_InitCond['gdd_1_0'] = 0
        dict_InitCond['gdd_2_0'] = (np.dot(
            df_gridSurfaceChar.loc[grid, 'sddfull'],
            df_gridSurfaceChar.loc[grid, 'sfr'][2: 5]) / fr_veg_sum
            if fr_veg_sum > 0
            else 0)
        dict_InitCond['laiinitialevetr'] = df_gridSurfaceChar.loc[
            grid, 'laimin'][0]
        dict_InitCond['laiinitialdectr'] = df_gridSurfaceChar.loc[
            grid, 'laimin'][1]
        dict_InitCond['laiinitialgrass'] = df_gridSurfaceChar.loc[
            grid, 'laimin'][2]
        dict_InitCond['albevetr0'] = df_gridSurfaceChar.loc[
            grid, 'albmin_evetr']
        dict_InitCond['albdectr0'] = df_gridSurfaceChar.loc[
            grid, 'albmin_dectr']
        dict_InitCond['albgrass0'] = df_gridSurfaceChar.loc[
            grid, 'albmin_grass']
        dict_InitCond['decidcap0'] = df_gridSurfaceChar.loc[
            grid, 'capmin_dec']
        dict_InitCond['porosity0'] = df_gridSurfaceChar.loc[
            grid, 'pormax_dec']
    # END: update vegetation-related variables according to LeavesOutInitially

    # snowflag = 1 if snow-related modules are enabled or 0 otherwise
    snowflag = (
        0 if (dict_ModConfig['snowuse'] == 0 or
              dict_InitCond['snowinitially'] == 0)
        else 1)

    # state_init: temporally-varying states from timestep to timestep
    dict_StateInit = {
        # TODO: water use patterns, which is currently not used
        'wuday_id': 0. * np.ones(9, order='F'),
        'numcapita': df_gridSurfaceChar.loc[
            grid,
            ['popdensdaytime', 'popdensnighttime']].mean(),
        'qn1_av': nan * np.ones(1),
        'qn1_s_av': nan * np.ones(1),
        'dqndt': nan * np.ones(1),
        'dqnsdt': nan * np.ones(1),
        # 'qn1_s_store_grid': nan * np.ones(nsh),
        # 'qn1_store_grid': nan * np.ones(nsh),

        # snow:
        'snowalb': dict_InitCond['snowalb0'],
        'snowfallcum': 0.,
        'icefrac': 0.2 * np.ones(7),
        # Initial liquid (melted) water for each surface
        'meltwaterstore': snowflag * np.array(
            [dict_InitCond[var] for var in ['snowwaterpavedstate',
                                            'snowwaterbldgsstate',
                                            'snowwaterevetrstate',
                                            'snowwaterdectrstate',
                                            'snowwatergrassstate',
                                            'snowwaterbsoilstate',
                                            'snowwaterwaterstate']],
            order='F'),
        'snowdens': snowflag * np.array(
            [dict_InitCond[var] for var in ['snowdenspaved',
                                            'snowdensbldgs',
                                            'snowdensevetr',
                                            'snowdensdectr',
                                            'snowdensgrass',
                                            'snowdensbsoil',
                                            'snowdenswater']],
            order='F'),
        'snowfrac': snowflag * np.array(
            [dict_InitCond[var] for var in ['snowfracpaved',
                                            'snowfracbldgs',
                                            'snowfracevetr',
                                            'snowfracdectr',
                                            'snowfracgrass',
                                            'snowfracbsoil',
                                            'snowfracwater']],
            order='F'),
        'snowpack': snowflag * np.array(
            [dict_InitCond[var] for var in ['snowpackpaved',
                                            'snowpackbldgs',
                                            'snowpackevetr',
                                            'snowpackdectr',
                                            'snowpackgrass',
                                            'snowpackbsoil',
                                            'snowpackwater']],
            order='F'),

        # Initial soil stores for each surface (below ground)
        'soilmoist_id': np.array(
            [dict_InitCond[var] for var in ['soilstorepavedstate',
                                            'soilstorebldgsstate',
                                            'soilstoreevetrstate',
                                            'soilstoredectrstate',
                                            'soilstoregrassstate',
                                            'soilstorebsoilstate']]
            + [0.], order='F'),

        # Initial wetness status of each surface (above ground)
        'state_id': np.array(
            [dict_InitCond[var] for var in ['pavedstate',
                                            'bldgsstate',
                                            'evetrstate',
                                            'dectrstate',
                                            'grassstate',
                                            'bsoilstate',
                                            'waterstate']],
            order='F'),

        # mean Tair of past 24 hours, NOT used by supy as it is for ESTM
        'tair24hr': 273.15 * np.ones(24 * nsh),


        # vegetation related parameters:
        'porosity_id': dict_InitCond['porosity0'] * np.ones(1, order='F'),
        'albdectr_id': dict_InitCond['albdectr0'] * np.ones(1, order='F'),
        'albevetr_id': dict_InitCond['albevetr0'] * np.ones(1, order='F'),
        'albgrass_id': dict_InitCond['albgrass0'] * np.ones(1, order='F'),
        'decidcap_id': dict_InitCond['decidcap0'] * np.ones(1, order='F'),
        # leaf area index:
        'lai_id': 1. * np.zeros(3),
        # growing degree days:
        'gdd_id': np.zeros(5),
        # heating degree days:
        'hdd_id': np.zeros(12)
    }

    # lai-related parameters:
    # dict_StateInit['lai'][4 + id_prev] = [
    #     dict_InitCond[var] for var in ['laiinitialevetr',
    #                                    'laiinitialdectr',
    #                                    'laiinitialgrass']]

    # hdd-related parameters:
    # print 'temp_c0', dict_InitCond['temp_c0']
    # dif_basethdd_temp_c0 = (
    #     df_gridSurfaceChar.loc[grid, 'basethdd']
    #     - dict_InitCond['temp_c0'])
    # gamma1 = (1 if dif_basethdd_temp_c0 >= 0 else 0)
    # gamma2 = (1 if dif_basethdd_temp_c0 <= 0 else 0)
    # hdd1 = gamma1 * dif_basethdd_temp_c0
    # hdd2 = gamma2 * (-1 * dif_basethdd_temp_c0)
    # update hdd:
    # dict_StateInit['hdd_id'][:6] = [
    #     hdd1, hdd2, dict_InitCond['temp_c0'],
    #     0, 0, dict_InitCond['dayssincerain']]
    # dict_StateInit['hdd_id'][4 + id_prev - 3:
    #                       4 + id_prev, 3 - 1] = dict_InitCond['temp_c0']
    # dict_StateInit['hdd_id'][4 + id_prev, 5] = dict_InitCond['dayssincerain']

    # gdd-related parameters:
    dict_StateInit['gdd_id'][2] = 90
    dict_StateInit['gdd_id'][3] = -90
    # dict_StateInit['gdd_id'][id_prev, 0] = dict_InitCond['gdd_1_0']
    # dict_StateInit['gdd_id'][id_prev, 1] = dict_InitCond['gdd_2_0']

    # update timestep info:
    dict_StateInit.update(tstep_prev=tstep)
    dict_StateInit.update(dt_since_start=0)

    # pack final results with Fortran ordering
    dict_StateInit = {k: np.array(v, order='F')
                      for k, v in dict_StateInit.items()}

    return dict_StateInit


# convert dict_InitCond to pandas Series and DataFrame
def init_SUEWS_pd(dir_input):  # return pd.DataFrame
    # dict_mod_cfg, dict_state_init = sp.init_SUEWS_dict(dir_start)
    dict_mod_cfg, dict_state_init = init_SUEWS_dict(dir_input)
    # ser_mod_cfg: all static model configuration info
    ser_mod_cfg = pd.Series(dict_mod_cfg)
    # df_state_init: initial conditions for SUEWS simulations
    df_state_init = df.from_dict(dict_state_init).T
    df_state_init.index.set_names('grid', inplace=True)

    return ser_mod_cfg, df_state_init

# input processing code end here
##############################################################################


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

# test for performance
dict_var_inout = {k: None for k in set_var_inout}


# high-level wrapper: suews_cal_tstep
def suews_cal_tstep(dict_state_start, dict_met_forcing_tstep):
    # use single dict as input for suews_cal_main
    dict_input = dict_state_start.copy()
    dict_input.update(dict_met_forcing_tstep)
    dict_input = {k: dict_input[k] for k in list_var_input}

    # main calculation:
    res_suews_tstep = sd.suews_cal_main(**dict_input)

    # update state variables
    dict_state_end = dict_state_start.copy()
    dict_state_end.update({var: copy.copy(dict_input[var])
                           for var in list_var_inout})

    # update timestep info
    dict_state_end['tstep_prev'] = dict_state_end['tstep']
    dict_state_end['dt_since_start'] += dict_state_end['tstep']

    # pack output
    dict_output = {k: v for k, v in zip(
        list_var_output, res_suews_tstep)}

    return dict_state_end, dict_output


# 2. compact wrapper for running a whole simulation
# # main calculation
# input as DataFrame
def run_suews_df(df_forcing, df_init):
    # initialise dicts for holding results and model states
    dict_state = {}
    dict_output = {}
    # start tstep retrived from forcing data
    t_start = df_forcing.index[0]
    # convert df to dict with `itertuples` for better performance
    dict_forcing = {row.Index: row._asdict()
                    for row in df_forcing.itertuples()}
    # dict_forcing = {tstep: met_forcing_tstep.to_dict()
    #                 for tstep, met_forcing_tstep
    #                 in df_forcing.iterrows()}
    # grid list determined by initial states
    grid_list = df_init.index
    # dict_state is used to save model states for later use
    dict_state = {(t_start, grid): series_state_init.to_dict()
                  for grid, series_state_init
                  in copy.deepcopy(df_init).iterrows()}

    # temporal loop
    for tstep in df_forcing.index:
        # print 'tstep at', tstep
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
            dict_state.update({(tstep + 1, grid): dict_state_end})

    # pack results as easier DataFrames
    df_output = pack_df_output(dict_output)
    df_state = pack_df_state(dict_state)
    # df_output = dict_output
    # df_state = dict_state

    return df_output, df_state


# main calculation end here
##############################################################################


##############################################################################
# post-processing part
# get variable information from Fortran
def get_output_info_df():
    size_var_list = sd.output_size()
    var_list_x = [np.array(sd.output_name_n(i))
                  for i in np.arange(size_var_list) + 1]
    df_var_list = pd.DataFrame(var_list_x, columns=['var', 'group', 'aggm'])
    df_var_list = pd.DataFrame(var_list_x, columns=['var', 'group', 'aggm'])
    df_var_list = df_var_list.applymap(lambda x: x.decode().strip())
    df_var_list_x = df_var_list.replace(r'^\s*$', np.nan, regex=True).dropna()
    var_dfm = df_var_list_x.set_index(['group', 'var'])
    return var_dfm


# get variable info as a DataFrame
# save `var_df` for later use
var_df = get_output_info_df()

# dict as var_df but keys in lowercase
var_df_lower = {group.lower(): group
                for group in var_df.index.levels[0].str.strip()}

# generate index for variables in different model groups


def gen_group_cols(group_x):
    # get correct group name by cleaning and swapping case
    group = group_x.replace('dataoutline', '').replace('line', '')
    # print group
    group = var_df_lower[group]
    header_group = np.apply_along_axis(
        list, 0, var_df.loc[['datetime', group]].index.values)[:, 1]

    # generate MultiIndex if not `datetimeline`
    if not group_x == 'datetimeline':
        index_group = pd.MultiIndex.from_product([[group], header_group],
                                                 names=['group', 'var'],
                                                 sortorder=None)
    else:
        index_group = header_group

    return index_group


# merge_grid: useful for both `dict_output` and `dict_state`
def pack_df_grid(dict_output):
    # pack all grid and times into index/columns
    df_xx = df.from_dict(dict_output, orient='index')
    # pack
    df_xx1 = df_xx.applymap(lambda s: pd.Series(s)).applymap(df.from_dict)
    df_xx2 = pd.concat({grid: pd.concat(
        df_xx1[grid].to_dict()).unstack().dropna(axis=1)
        for grid in df_xx1.columns})
    # drop redundant levels
    df_xx2.columns = df_xx2.columns.droplevel()
    # regroup by `grid`
    df_xx2.index.names = ['grid', 'time']
    gb_xx2 = df_xx2.groupby(level='grid')
    # merge results of each grid
    xx3 = gb_xx2.agg(lambda x: tuple(x.values)).applymap(np.array)

    return xx3


# generate MultiIndex for variable groups
def gen_index(varline_x):
    var_x = varline_x.replace('dataoutline', '').replace('line', '')
    group = var_df_lower[var_x]
    var = var_df.loc[group].index.tolist()
    mindex = pd.MultiIndex.from_product([[group], var], names=['group', 'var'])
    return mindex


# generate one MultiIndex from a whole dict
def gen_MultiIndex(dict_x):
    x_keys = dict_x.keys()
    mindex = pd.concat([gen_index(k).to_frame() for k in x_keys]).index
    return mindex


# generate one Series from a dict entry
def gen_Series(dict_x, varline_x):
    m_index = gen_index(varline_x)
    res_Series = pd.Series(dict_x[varline_x], index=m_index)
    return res_Series


# merge a whole dict into one Series
def comb_gen_Series(dict_x):
    x_keys = dict_x.keys()
    res_Series = pd.concat([gen_Series(dict_x, k) for k in x_keys])
    return res_Series


# pack up output of `run_suews`
def pack_df_output(dict_output):
    # TODO: add output levels as in the Fortran version
    df_output = pd.DataFrame(dict_output).T
    # df_output = pd.concat(dict_output).to_frame().unstack()
    # set index level names
    index = df_output.index.set_names(['tstep', 'grid'])
    # clean columns
    # df_output.columns = df_output.columns.droplevel()
    columns = gen_MultiIndex(df_output.iloc[0])
    values = np.apply_along_axis(np.hstack, 1, df_output.values)
    df_output = pd.DataFrame(values, index=index, columns=columns)
    return df_output


def pack_df_state(dict_state):
    df_state = pd.DataFrame(dict_state).T
    # df_state = pd.concat(dict_state).to_frame().unstack()
    # set index level names
    df_state.index = df_state.index.set_names(['tstep', 'grid'])
    # clean columns
    # df_state.columns = df_state.columns.droplevel()

    return df_state


##############################################################################
# auxiliary functions
# resolve path case issues
def path_insensitive(path):
    """
    Get a case-insensitive path for use on a case sensitive system.

    >>> path_insensitive('/Home')
    '/home'
    >>> path_insensitive('/Home/chris')
    '/home/chris'
    >>> path_insensitive('/HoME/CHris/')
    '/home/chris/'
    >>> path_insensitive('/home/CHRIS')
    '/home/chris'
    >>> path_insensitive('/Home/CHRIS/.gtk-bookmarks')
    '/home/chris/.gtk-bookmarks'
    >>> path_insensitive('/home/chris/.GTK-bookmarks')
    '/home/chris/.gtk-bookmarks'
    >>> path_insensitive('/HOME/Chris/.GTK-bookmarks')
    '/home/chris/.gtk-bookmarks'
    >>> path_insensitive("/HOME/Chris/I HOPE this doesn't exist")
    "/HOME/Chris/I HOPE this doesn't exist"
    """

    return _path_insensitive(path) or path


def _path_insensitive(path):
    """
    Recursive part of path_insensitive to do the work.
    """

    if path == '' or os.path.exists(path):
        return path

    base = os.path.basename(path)  # may be a directory or a file
    dirname = os.path.dirname(path)

    suffix = ''
    if not base:  # dir ends with a slash?
        if len(dirname) < len(path):
            suffix = path[:len(path) - len(dirname)]

        base = os.path.basename(dirname)
        dirname = os.path.dirname(dirname)

    if not os.path.exists(dirname):
        dirname = _path_insensitive(dirname)
        if not dirname:
            return

    # at this point, the directory exists but not the file

    try:  # we are expecting dirname to be a directory, but it could be a file
        files = os.listdir(dirname)
    except OSError:
        return

    baselow = base.lower()
    try:
        basefinal = next(fl for fl in files if fl.lower() == baselow)
    except StopIteration:
        return

    if basefinal:
        return os.path.join(dirname, basefinal) + suffix
    else:
        return


# convert pandas structures to Python native structures
def conv2PyData(df_x):
    try:  # if the data structures in dict_df_x have consistent dimensionality
        # convert to dict in a split way to expose data/values
        dict_x = df(df_x).to_dict('split')
        # convert data into native Python `list`
        dict_x['data'] = [np.array(var).tolist()
                          for var in dict_x['data']]
        # convert back to DataFrame and then to dict again, which is
        # Python-native
        dict_x_nat = df(**dict_x).to_dict()
    except:
        # if the data structures in dict_df_x have inconsistent dimensionality
        # convert to dict in a split way to expose data/values
        dict_x = df(df_x.T).to_dict('split')
        # convert data into native Python `list`
        dict_x['data'] = [np.array(var).tolist()
                          for var in dict_x['data']]
        # convert back to DataFrame and then to dict again, which is
        # Python-native
        dict_x_nat = df(**dict_x).to_dict()

    return dict_x_nat
##############################################################################
