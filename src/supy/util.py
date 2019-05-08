# supy utilities

import cdsapi
import time
from numpy import sqrt, cos, sin, deg2rad
from pathlib import Path
import atmosp
import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import interp1d


# locate the first position of period with in-between gaps
def loc_gap(ser_test, freq='1D', pattern='010'):
    rsmp = ser_test.resample(freq)
    ser_TF_10 = rsmp.apply(lambda ser: ser.isna().any()) * 1
    str_TF_10 = ser_TF_10.astype(str).str.cat()
    pos_gap = str_TF_10.find(pattern)
    loc_ser = ser_TF_10.iloc[pos_gap:pos_gap + len(pattern)].index
    return loc_ser


# fill gap with neighbouring days
def fill_gap_one(ser_test, freq='1D', pattern='010'):
    # resample into daily periods
    rsmp = ser_test.resample(freq)
    # locate the gaps according to gap pattern: 0 for NO gap, 1 for gapped
    loc_ser = loc_gap(ser_test, freq, pattern)

    # generator groups
    ser_find = (rsmp.get_group(x) for x in loc_ser)
    if len(loc_ser) == 0:
        return ser_test

    # assign series:
    # ser_prev: series prior to gapped period
    # ser_gap: series with gaps
    # ser_post: series after gapped period
    if pattern == '010':
        ser_prev, ser_gap, ser_post = ser_find
    elif pattern == '01':
        ser_prev, ser_gap = ser_find
        ser_post = pd.Series([])
    elif pattern == '10':
        ser_gap, ser_post = ser_find
        ser_prev = pd.Series([])

    # base series for gap filling
    ser_fill_base = pd.concat([ser_prev, ser_post])
    ser_fill = ser_fill_base.groupby([
        ser_fill_base.index.hour.rename('hr'),
        ser_fill_base.index.minute.rename('min')]).median().reset_index(
        drop=True)
    ser_fill.index = ser_gap.index

    # calculate rescaling factor
    scale_fill = (ser_fill / ser_gap).median()
    scale_fill = (1 if abs(scale_fill) > 10 else scale_fill)
    scale_fill = (1 if abs(scale_fill) < 0.1 else scale_fill)
    ser_fill_gap = ser_fill / scale_fill

    # fill in gaps with rescaled values of the
    ser_gap.loc[ser_gap.isna()] = ser_fill_gap.loc[ser_gap.isna()]
    ser_filled = pd.concat([ser_prev, ser_gap, ser_post])

    # fill the original gapped series
    ser_test_filled = ser_test.copy()
    ser_test_filled.loc[ser_filled.index] = ser_filled
    return ser_test_filled


# fill gaps iteratively
def fill_gap_all(ser_test, freq='1D'):
    ser_test_filled = ser_test.copy()
    ptn_list = ['010', '01', '10']
    while ser_test_filled.isna().any():
        # try to different gap patterns and fill gaps
        try:
            ptn_gap = next(ptn for ptn in ptn_list if len(
                loc_gap(ser_test_filled, freq, ptn)) > 0)
            ser_test_filled = fill_gap_one(ser_test_filled, freq, ptn_gap)
        except StopIteration:
            pass
    return ser_test_filled


#################################################################
# generate TMY dataframe from supy results
# weight class to determine constants for TMY generation
class _const:
    class ConstError(TypeError):
        pass

    class ConstCaseError(ConstError):
        pass

    def __setattr__(self, name, value):
        if name in self.__dict__:
            raise self.ConstError("can't change const %s" % name)
        if not name.isupper():
            raise self.ConstCaseError(
                'const name "%s" is not all uppercase' % name)
        self.__dict__[name] = value


const = _const()

const.T_MEAN = (2 / 24)
const.T_MAX = (1 / 24)
const.T_MIN = (1 / 24)
const.T_RANGE = 0  # todo
const.RH_MEAN = 2 / 24
const.RH_MAX = 1 / 24
const.RH_MIN = 1 / 24
const.RH_RANGE = 0  # todo
const.WIND_MEAN = 2 / 24
const.WIND_MAX = 2 / 24
const.WIND_MIN = 0
const.WIND_RANGE = 0  # todo
const.WIND_DIRECTION = 0  # todo
const.SOLAR_RADIATION_GLOBAL = 12 / 24
const.SOLAR_RADIATION_DIRECT = 0  # todo


def gen_score_list(length):
    list_score = (np.arange(length) + 0.5) / length
    return list_score


def gen_score_ser(ser_test):
    ser_score = ser_test.sort_values(ascending=True)
    length = ser_score.size
    list_score = (np.arange(length) + 0.5) / length
    ser_score.loc[:] = list_score
    # ser_score.loc[:] = gen_score_list(ser_score.size)
    return ser_score


def gen_FS_DF(df_output):
    """generate DataFrame of scores.

    Parameters
    ----------
    df_WS_data : type
        Description of parameter `df_WS_data`.

    Returns
    -------
    type
        Description of returned object.

    """
    df_day = pd.pivot_table(
        df_output,
        values=['T2', 'U10', 'Kdown', 'RH2'],
        index=['Year', 'Month', 'Day'],
        aggfunc=[min, max, np.mean, ])
    df_day_all_year = pd.pivot_table(
        df_output,
        values=['T2', 'U10', 'Kdown', 'RH2'],
        index=['Month', 'Day'],
        aggfunc=[min, max, np.mean, ])

    array_yr_mon = df_day.index.droplevel(
        'Day').to_frame().drop_duplicates().values

    df_fs = pd.DataFrame(
        {(yr, mon):
         (df_day.loc[(yr, mon)].apply(gen_score_ser) -
          df_day_all_year.loc[mon].apply(gen_score_ser)).abs().mean()
         for yr, mon in array_yr_mon})

    return df_fs


def gen_WS_DF(df_WS_data):
    """generate DataFrame of weighted sums.

    Parameters
    ----------
    df_WS_data : type
        Description of parameter `df_WS_data`.

    Returns
    -------
    type
        Description of returned object.

    """
    df_fs = gen_FS_DF(df_WS_data)

    list_index = [('mean', 'T2'), ('max', 'T2'), ('min', 'T2'),
                  ('mean', 'U10'), ('max', 'U10'), ('min', 'U10'),
                  ('mean', 'RH2'), ('max', 'RH2'), ('min', 'RH2'),
                  ('mean', 'Kdown')]

    list_const = [getattr(const, attr)
                  for attr in ['T_MEAN', 'T_MAX', 'T_MIN',
                               'WIND_MEAN', 'WIND_MAX', 'WIND_MIN',
                               'RH_MEAN', 'RH_MAX', 'RH_MIN',
                               'SOLAR_RADIATION_GLOBAL']]
    list_ws = [df_fs.loc[idx] * cst
               for idx, cst
               in zip(list_index, list_const)]
    df_ws = pd.concat(list_ws, axis=1).sum(axis=1).unstack().dropna()

    return df_ws


# def gen_WS_DF(df_WS_data):
#     """generate DataFrame of weighted sums.

#     Parameters
#     ----------
#     df_WS_data : type
#         Description of parameter `df_WS_data`.

#     Returns
#     -------
#     type
#         Description of returned object.

#     """
#     df_fs = gen_FS_DF(df_WS_data)

#     list_index = [('mean', 'T2'), ('max', 'T2'), ('min', 'T2'),
#                   ('mean', 'U10'), ('max', 'U10'), ('min', 'U10'),
#                   ('mean', 'RH2'), ('max', 'RH2'), ('min', 'RH2'),
#                   ('mean', 'Kdown')]

#     list_const = [getattr(const, attr)
#                   for attr in ['T_MEAN', 'T_MAX', 'T_MIN',
#                                'WIND_MEAN', 'WIND_MAX', 'WIND_MIN',
#                                'RH_MEAN', 'RH_MAX', 'RH_MIN',
#                                'SOLAR_RADIATION_GLOBAL']]
#     list_ws = [df_fs.loc[idx] * cst for idx,
#                cst in zip(list_index, list_const)]
#     df_ws = pd.concat(list_ws, axis=1).sum(axis=1).unstack().dropna()

#     return df_ws


def pick_year(df_ws, df_output, n=5):
    df_day = pd.pivot_table(
        df_output, values='Kdown',
        index=['Year', 'Month', 'Day'],
        aggfunc=[np.mean, ])

    df_day_all_year = pd.pivot_table(
        df_output, values='Kdown',
        index=['Month', 'Day'],
        aggfunc=[np.mean, ])

    array_yr_mon = df_day.index.droplevel(
        'Day').to_frame().drop_duplicates().values

    df_rmsd = pd.DataFrame(
        {(yr, mon):
         np.sqrt(
         np.square(
             df_day.loc[(yr, mon)] - df_day_all_year.loc[mon]).mean())
         for yr, mon in array_yr_mon}).stack().T.dropna()
    df_rmsd.columns = df_rmsd.columns.droplevel([0, 1])

    year_nsmallest = df_ws.apply(lambda ser: ser.nsmallest(n).index)

    year_sel = df_rmsd.apply(
        lambda ser: ser.loc[year_nsmallest[ser.name]]).idxmin()

    return year_sel


# headers of standard EPW files
header_EPW = '''
Year
Month
Day
Hour
Minute
Data Source and Uncertainty Flags
Dry Bulb Temperature
Dew Point Temperature
Relative Humidity
Atmospheric Station Pressure
Extraterrestrial Horizontal Radiation
Extraterrestrial Direct Normal Radiation
Horizontal Infrared Radiation Intensity
Global Horizontal Radiation
Direct Normal Radiation
Diffuse Horizontal Radiation
Global Horizontal Illuminance
Direct Normal Illuminance
Diffuse Horizontal Illuminance
Zenith Luminance
Wind Direction
Wind Speed
Total Sky Cover
Opaque Sky Cover
Visibility
Ceiling Height
Present Weather Observation
Present Weather Codes
Precipitable Water
Aerosol Optical Depth
Snow Depth
Days Since Last Snowfall
Albedo
Liquid Precipitation Depth
Liquid Precipitation Quantity
    '''

# list of variables in EPW
list_var_EPW = header_EPW.split('\n')[1:-1]

# dict: SuPy variables -> EPW standard names
dict_supy_epw = {
    'Kdown': 'Global Horizontal Radiation',
    'T2': 'Dry Bulb Temperature',
    'RH2': 'Relative Humidity',
    'U10': 'Wind Speed',

}
dict_epw_supy = {v: k for k, v in dict_supy_epw.items()}


def gen_TMY(df_output):
    '''generate TMY (typical meteorological year) from SuPy output.

    Parameters
    ----------
    df_output : pandas.DataFrame
        Output from `run_supy`: longterm (e.g., >10 years) simulation results, otherwise not very useful.

    '''

    # calculate weighted score
    ws = gen_WS_DF(df_output)

    # select year
    year_sel = pick_year(ws, df_output, n=5)

    # generate TMY data
    df_TMY = pd.concat(
        # shift(1) here is to conform the convention that
        # timestamps refer to the preceding period
        # [df_output.shift(1).groupby(['Month', 'Year']).get_group(grp)
        # shift(1) is not necessary
        [df_output.groupby(['Month', 'Year']).get_group(grp)
         for grp in year_sel.items()])

    # df_TMY = df_TMY.rename(columns=dict_supy_epw)

    return df_TMY


# function to read in EPW file
def read_epw(path_epw):
    df_tmy = pd.read_csv(path_epw, skiprows=8, sep=u',', header=None)
    df_tmy.columns = [x.strip() for x in header_EPW.split('\n')[1:-1]]
    df_tmy['DateTime'] = pd.to_datetime(
        pd.to_datetime(
            df_tmy['Year']*10000+df_tmy['Month']*100+df_tmy['Day'],
            format='%Y%m%d')+pd.to_timedelta(df_tmy['Hour'], unit='h'))
    df_tmy = df_tmy.set_index('DateTime')
    return df_tmy


# generate EPW file from `df_TMY`
def gen_epw(df_tmy, path_epw=Path('./uTMY.epw'), ratio_dif_dir=0.15):
    df_tmy = df_tmy.copy()
    # df_tmy = pd.concat([df_tmy.iloc[1:], df_tmy.iloc[[0]]])
    # adding necessary variables that can be derive from supy output
    df_tmy['Dew Point Temperature'] = atmosp.calculate(
        "Td",
        T=df_tmy['T2'].values + 273.15,
        qv=(df_tmy['Q2'].values), qv_unit='g/kg',
        RH=df_tmy['RH2'].values,
        rho=1.23)-273.15
    df_tmy['Atmospheric Station Pressure'] = atmosp.calculate(
        "p",
        T=df_tmy['T2'].values + 273.15,
        qv=(df_tmy['Q2'].values), qv_unit='g/kg',
        RH=df_tmy['RH2'].values,
        rho=1.23)-273.15
    # index = df_TMY.index
    # df_TMY = df_TMY.iloc[1:].append(df_TMY.ix[0])
    # df_TMY.index = index
    df_tmy['Year'] = df_tmy.index.year
    df_tmy['Month'] = df_tmy.index.month
    df_tmy['Day'] = df_tmy.index.day
    df_tmy['Hour'] = df_tmy.index.hour
    df_tmy['Minute'] = df_tmy.index.minute
    # convert air pressure to Pa
    df_tmy['Atmospheric Station Pressure'] *= 1000

    # df_TMY['Kdown'] = df_TMY['Kdown']*2.4
    # processing solar radiation components
    df_tmy.loc[df_tmy['Kdown'] < 0.001, 'Kdown'] = 0
    # direct beam
    frac_dir = (1/(1+ratio_dif_dir))
    df_tmy['Direct Normal Radiation'] = df_tmy['Kdown'] * frac_dir
    # diffuse radiation
    frac_dif = 1-(1/(1+ratio_dif_dir))
    df_tmy['Diffuse Horizontal Radiation'] = df_tmy['Kdown'] * frac_dif

    # conform column names to EPW standard
    df_TMY_x = df_tmy.rename(columns=dict_supy_epw)

    # initialise df_epw for EPW output
    df_epw = pd.DataFrame(columns=list_var_EPW, index=df_tmy.index)

    # dict of default values
    dict_var_dft = {
        'Data Source and Uncertainty Flags': -9992,
        'Extraterrestrial Horizontal Radiation': 9999,
        'Extraterrestrial Direct Normal Radiation': 9999,
        'Horizontal Infrared Radiation Intensity': 9999,
        'Direct Normal Radiation': 9999,
        'Global Horizontal Illuminance': 9999999,
        'Direct Normal Illuminance': 9999999,
        'Diffuse Horizontal Illuminance': 9999999,
        'Zenith Luminance': 9999,
        'Wind Direction': 999,
        'Total Sky Cover': 99,
        'Opaque Sky Cover': 99,
        'Visibility': 9999,
        'Ceiling Height': 99999,
        'Present Weather Observation': 9999,
        'Present Weather Codes': 9999,
        'Precipitable Water': 999,
        'Aerosol Optical Depth': 999,
        'Snow Depth': 999,
        'Days Since Last Snowfall': 99,
        'Albedo': 999,
        'Liquid Precipitation Depth': 999,
        'Liquid Precipitation Quantity': 999,
    }
    for var in list_var_EPW:
        try:
            df_epw[var] = df_TMY_x[var].values
        except:
            # print(f'{var} not existing! This variable will be filled with default value {dict_var_dft[var]}')
            try:
                df_epw[var] = np.ones(len(df_epw))*dict_var_dft[var]
            except:
                df_epw[var] = np.nan
    # fill 'Data Source and Uncertainty Flags'
    df_epw['Data Source and Uncertainty Flags'] = '?9?9?9?9E0?9?9?9*9*9?9*9*9?9*9*9?9?9*9*_*9*9*9*9*9'
    df_epw['Global Horizontal Radiation'] = np.ones(len(df_epw))*9999
    df_epw.index = df_TMY_x.index

    # convert `0h` to `24h` and take care of `day`
    loc_24h = df_epw.index == df_epw.index.normalize()
    ser_24h = df_epw.loc[loc_24h].index-pd.Timedelta('1h')
    df_epw.loc[loc_24h, 'Year'] = ser_24h.year
    df_epw.loc[loc_24h, 'Month'] = ser_24h.month
    df_epw.loc[loc_24h, 'Day'] = ser_24h.day
    df_epw.loc[loc_24h, 'Hour'] = 24

    df_epw = df_epw.sort_values(['Month', 'Day', 'Hour'], axis=0)

    # save pure data to a csv for formatting
    df_epw.to_csv(path_epw, index=None, header=None)
    text_data = path_epw.read_text().split('\n')
    # delete the csv file
    path_epw.unlink()

    text_meta = '''
LOCATION,Chongqing Shapingba,Chongqing,CHN,CSWD,575160,29.58,106.47,8,259.1
DESIGN CONDITIONS,1,Climate Design Data 2009 ASHRAE Handbook,,Heating,1,3.2,4.2,-0.2,3.8,6.5,1.3,4.3,6.2,4.9,7.6,4.3,7.5,1.4,0,Cooling,7,7.4,36.9,25.6,35.5,25.6,34.2,25.4,27.4,32.7,26.9,32.2,26.4,31.6,2.5,110,26.1,22.2,30.2,25.6,21.5,29.8,25.1,20.8,29.3,89.3,32.7,86.9,32.5,84.7,31.7,909,Extremes,5.1,4.3,3.6,35.4,1.1,38.8,1.3,1.6,0.1,40,-0.6,40.9,-1.4,41.8,-2.3,43
TYPICAL/EXTREME PERIODS,6,Summer - Week Nearest Max Temperature For Period,Extreme,7/27,8/ 2,Summer - Week Nearest Average Temperature For Period,Typical,7/ 6,7/12,Winter - Week Nearest Min Temperature For Period,Extreme,12/22,1/ 5,Winter - Week Nearest Average Temperature For Period,Typical,1/13,1/19,Autumn - Week Nearest Average Temperature For Period,Typical,10/13,10/19,Spring - Week Nearest Average Temperature For Period,Typical,4/12,4/18
GROUND TEMPERATURES,3,.5,,,,13.31,10.23,9.39,10.12,14.28,18.95,23.34,26.51,27.44,25.95,22.35,17.82,2,,,,16.09,13.20,11.82,11.77,13.97,17.16,20.59,23.54,25.06,24.77,22.74,19.63,4,,,,17.90,15.65,14.27,13.85,14.66,16.52,18.83,21.09,22.62,22.98,22.11,20.29
HOLIDAYS/DAYLIGHT SAVINGS,No,0,0,0
COMMENTS 1,  generated by SUEWS model
COMMENTS 2, todo
DATA PERIODS,1,1,Data,Sunday, 1/ 1,12/31
    '''
    text_meta = text_meta.split('\n')[1:-1]
    # lines = []
    text_epw = '\n'.join(text_meta+text_data)
    # with open(path_epw, 'r') as f:
    #     for line in f:
    #         lines.append(line)
    #     lines.insert(0, text_meta[1:])
    #     s = ''.join(lines)

    # write out the actual EPW file
    with open(path_epw, 'w') as fp:
        fp.write(text_epw)

    return df_epw, text_meta, path_epw


################################################
# more ERA-5 related functions
################################################

# utility functions
def roundPartial(value, resolution):
    return round(value / resolution) * resolution


"""Geopotential Functions on WGS84 Reference Ellipsoid

This module contains code for converting Geopotential to Geometric and vice-versa on the WGS84 reference ellipsoid

ERA-5 utility functions from Chris Roth
# https://pypi.org/project/eratools/

"""


Rmax_WGS84 = 6378137
Rmin_WGS84 = Rmax_WGS84 * (1 - 1/298.257223563)


def _geoid_radius(latitude: float) -> float:
    """Calculates the GEOID radius at a given latitude

    Parameters
    ----------
    latitude : float
        Latitude (degrees)

    Returns
    -------
    R : float
        GEOID Radius (meters)
    """
    lat = deg2rad(latitude)
    return sqrt(1/(cos(lat) ** 2 / Rmax_WGS84 ** 2 + sin(lat) ** 2 / Rmin_WGS84 ** 2))


def geometric2geopotential(z: float, latitude: float) -> float:
    """Converts geometric height to geopoential height

    Parameters
    ----------
    z : float
        Geometric height (meters)
    latitude : float
        Latitude (degrees)

    Returns
    -------
    h : float
        Geopotential Height (meters) above the reference ellipsoid
    """
    twolat = deg2rad(2 * latitude)
    g = 9.80616 * (1 - 0.002637*cos(twolat) + 0.0000059*cos(twolat)**2)
    re = _geoid_radius(latitude)
    return z * g * re / (re + z)


def geopotential2geometric(h: float, latitude: float) -> float:
    """Converts geopoential height to geometric height

    Parameters
    ----------
    h : float
        Geopotential height (meters)
    latitude : float
        Latitude (degrees)

    Returns
    -------
    z : float
        Geometric Height (meters) above the reference ellipsoid
    """
    twolat = deg2rad(2 * latitude)
    g = 9.80616 * (1 - 0.002637*cos(twolat) + 0.0000059*cos(twolat)**2)
    re = _geoid_radius(latitude)
    return h * re / (g * re - h)


# functions to interpolate the atmospheric variables to a specified height/altitude
def get_ser_val_alt(lat: float, lon: float,
                    da_alt_x: xr.DataArray,
                    da_alt: xr.DataArray, da_val: xr.DataArray)->pd.Series:
    '''interpolate atmospheric variable to a specified altitude

    Parameters
    ----------
    lat : float
        latitude of specified site
    lon : float
        longitude of specified site
    da_alt_x : xr.DataArray
        desired altitude to interpolate variable at
    da_alt : xr.DataArray
        altitude associated with `da_val`: variable array to interpolate
    da_val : xr.DataArray
        atmospheric varialble to interpolate

    Returns
    -------
    pd.Series
        interpolated values at the specified altitude of site positioned by [`lat`, `lon`]
    '''

    alt_t_1d = da_alt.sel(
        latitude=lat, longitude=lon, method='nearest')
    val_t_1d = da_val.sel(
        latitude=lat, longitude=lon, method='nearest')
    alt_x = da_alt_x.sel(
        latitude=lat, longitude=lon, method='nearest')[0]
    val_alt = np.array(
        [interp1d(alt_1d, val_1d)(alt_x)
         for alt_1d, val_1d
         in zip(alt_t_1d, val_t_1d)])
    ser_alt = pd.Series(
        val_alt,
        index=da_val.time.values,
        name=da_val.name,
    )
    return ser_alt


def get_df_val_alt(lat: float, lon: float, da_alt_meas: xr.DataArray, ds_val: xr.Dataset):
    '''interpolate atmospheric variables to a specified altitude

    Parameters
    ----------
    lat : float
        latitude of specified site
    lon : float
        longitude of specified site
    da_alt_x : xr.DataArray
        desired altitude to interpolate variable at
    da_alt : xr.DataArray
        altitude associated with `da_val`: variable array to interpolate
    da_val : xr.DataArray
        atmospheric varialble to interpolate

    Returns
    -------
    pd.DataFrame
        interpolated values at the specified altitude of site positioned by [`lat`, `lon`]
    '''
    da_alt = geopotential2geometric(ds_val.z, ds_val.latitude)
    # generate pressure series for grid x
    da_alt_x = da_alt.sel(
        latitude=lat, longitude=lon, method='nearest')
    alt_meas_x = da_alt_meas.sel(
        latitude=lat, longitude=lon, method='nearest')[0]

    val_pres = np.array([interp1d(alt, da_alt_x.level)(alt_meas_x)
                         for alt in da_alt_x])
    df_val_alt = pd.concat(
        [get_ser_val_alt(
            lat, lon, da_alt_meas, da_alt, ds_val[var])
         for var in ds_val.data_vars],
        axis=1
    )
    #     add pressure
    df_val_alt['p'] = val_pres
    df_val_alt.index = df_val_alt.index.set_names('time')
    df_val_alt.columns = df_val_alt.columns.set_names('var')

    return df_val_alt


# cds download related functions
def gen_dict_dt(dt_index):
    list_key = ['year', 'month', 'day', 'time']
    list_fmt = ['%Y', '%m', '%d', '%H:%M']
    dict_dt = {
        k: dt_index.dt.strftime(fmt).unique().tolist()
        for k, fmt in zip(list_key, list_fmt)
    }
    # dict_dt = {
    #     'year': dt_index.dt.strftime('%Y').unique().tolist(),
    #     'month': dt_index.dt.strftime('%m').unique().tolist(),
    #     'day': dt_index.dt.strftime('%d').unique().tolist(),
    #     'time': dt_index.dt.strftime('%H:%M').unique().tolist(),
    # }
    return dict_dt


def gen_dict_dt_sub(dt_index):
    # divide by [year, month] for surface level data
    ser_dict_sub = dt_index.groupby(
        dt_index.dt.strftime('%Y%m')
    ).apply(gen_dict_dt)
    dict_sub = ser_dict_sub.unstack().T.to_dict()
    return dict_sub


def gen_fn(dict_x):
    lat_x, lon_x = dict_x['area'][:2]
    yr_x = dict_x['year'][0]
    mon_x = dict_x['month'][0]
    type_x = 'sfc' if 'orography' in dict_x['variable'] else 'ml'
    fn_x = f'{lat_x}N{lon_x}E-{yr_x}{mon_x}-{type_x}.nc'
    return fn_x


# dict_x: a dict describing download elements
def gen_dict_proc(dict_x):
    type_x = 'sfc' if 'orography' in dict_x['variable'] else 'ml'
    dict_feed = {
        'sfc': 'reanalysis-era5-single-levels',
        'ml': 'reanalysis-era5-pressure-levels',
    }
    feed_x = dict_feed[type_x]
    dict_proc = dict(
        name=feed_x,
        request=dict_x,
        target=gen_fn(dict_x),
    )

    return dict_proc


list_var_sfc = [
    '10m_u_component_of_wind',
    '10m_v_component_of_wind',
    '2m_dewpoint_temperature',
    '2m_temperature',
    'orography',
    'surface_pressure',
    'surface_solar_radiation_downwards',
    'surface_thermal_radiation_downwards',
    'total_precipitation',
]

list_var_ml = [
    'geopotential',
    #     'relative_humidity',
    'specific_humidity',
    'temperature',
    'u_component_of_wind',
    'v_component_of_wind',
]

list_pres_level = [
    '1', '2', '3',
    '5', '7', '10',
    '20', '30', '50',
    '70', '100', '125',
    '150', '175', '200',
    '225', '250', '300',
    '350', '400', '450',
    '500', '550', '600',
    '650', '700', '750',
    '775', '800', '825',
    '850', '875', '900',
    '925', '950', '975',
    '1000',
]

# generate a dict of reqs kwargs for (lat_x,lon_x) spanning [start, end]


def gen_req_sfc(lat_x, lon_x, start, end, grid=[0.125, 0.125], scale=0):
    '''generate a dict of reqs kwargs for (lat_x,lon_x) spanning [start, end]

    Parameters
    ----------
    lat_x : [type]
        [description]
    lon_x : [type]
        [description]
    start : [type]
        [description]
    end : [type]
        [description]
    grid : list, optional
        [description] (the default is [0.125, 0.125], which [default_description])
    scale : int, optional
        [description] (the default is 0, which [default_description])

    Returns
    -------
    [type]
        [description]

    Examples
    --------
    >>> gen_req_sfc(28, 116, '2015-01', '2015-01-31 23', grid=[0.125, 0.125], scale=0)

    '''

    # scale is a factor to rescale grid size
    size = grid[0]*scale
    # generate pd.Series for timestamps
    ser_datetime = pd.date_range(start, end, freq='1h').to_series()
    # surface requests
    lat_c, lon_c = (roundPartial(x, grid[0]) for x in [lat_x, lon_x])
    area = [lat_c+size, lon_c-size, lat_c-size, lon_c+size]
    dict_req_sfc = {
        'variable': list_var_sfc,
        'product_type': 'reanalysis',
        'area': area,
        'grid': grid,
        'format': 'netcdf'
    }
    list_dict_req_sfc = [
        {**dict_req_sfc, **dict_dt}
        for dict_dt
        in list(gen_dict_dt_sub(ser_datetime).values())
    ]
    dict_req_sfc = {
        gen_fn(dict_req): gen_dict_proc(dict_req)
        for dict_req in list_dict_req_sfc
    }
    return dict_req_sfc


def sel_list_pres(ds_sfc_x):
    '''
    select proper levels for model level data download
    '''
    p_min, p_max = ds_sfc_x.sp.min().values, ds_sfc_x.sp.max().values
    list_pres_level = [
        '1', '2', '3',
        '5', '7', '10',
        '20', '30', '50',
        '70', '100', '125',
        '150', '175', '200',
        '225', '250', '300',
        '350', '400', '450',
        '500', '550', '600',
        '650', '700', '750',
        '775', '800', '825',
        '850', '875', '900',
        '925', '950', '975',
        '1000',
    ]
    ser_pres_level = pd.Series(list_pres_level).map(int)*100
    pos_lev_max, pos_lev_min = (
        ser_pres_level[ser_pres_level > p_max].idxmin(),
        ser_pres_level[ser_pres_level < p_min].idxmax()
    )
    list_pres_sel = ser_pres_level.loc[pos_lev_min:pos_lev_max]/100
    list_pres_sel = list_pres_sel.map(int).map(str).to_list()
    return list_pres_sel


# sel_list_pres(ds_sfc_x)


# for each sfc data file, determine the necessary vertical levels to model level data download
def gen_req_ml(fn_sfc, grid=[0.125, 0.125], scale=0):
    ds_sfc_x = xr.open_dataset(fn_sfc)
    list_pres_sel = sel_list_pres(ds_sfc_x)
    size = grid[0]*scale
    lat_x, lon_x = ds_sfc_x.latitude.values[0], ds_sfc_x.longitude.values[0]
    lat_c, lon_c = (roundPartial(x, grid[0]) for x in [lat_x, lon_x])
    area = [lat_c+size, lon_c-size, lat_c-size, lon_c+size]
    idx_time = ds_sfc_x.time.to_pandas()
    dict_dt = list(gen_dict_dt_sub(idx_time).values())[0]
    # model level requests
    dict_req_ml = {
        'variable': list_var_ml,
        'product_type': 'reanalysis',
        'area': area,
        'grid': grid,
        'format': 'netcdf'
    }
    dict_req_ml.update({'level': list_pres_sel})
    dict_req_ml.update(dict_dt)
    dict_req_ml = {
        gen_fn(dict_req_ml): gen_dict_proc(dict_req_ml)
    }
    return dict_req_ml


def download_era5(lat_x, lon_x, start, end, grid=[0.125, 0.125], scale=0):
    c = cdsapi.Client()
    dict_req_sfc = gen_req_sfc(lat_x, lon_x, start, end, grid=[
                               0.125, 0.125], scale=0)
    for fn_sfc, dict_req in dict_req_sfc.items():
        if not Path(fn_sfc).exists():
            print('To download:', fn_sfc)
            c.retrieve(**dict_req)
            time.sleep(.0100)

    dict_req_ml = {}
    for fn_sfc in dict_req_sfc.keys():
        if Path(fn_sfc).exists():
            print(f'{fn_sfc} exists!')
            dict_req = gen_req_ml(fn_sfc, grid, scale)
            dict_req_ml.update(dict_req)

    for fn_ml, dict_req in dict_req_ml.items():
        if Path(fn_ml).exists():
            print(f'{fn_ml} exists!')
            print('')
        else:
            print('To download:', fn_ml)
            c.retrieve(**dict_req)
            time.sleep(.0100)

    dict_req_all = {**dict_req_sfc, **dict_req_ml}
    return dict_req_all
