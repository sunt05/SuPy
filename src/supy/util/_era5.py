try:
    import cdsapi
    import xarray as xr
    import atmosp
except:
    pass

import time
from pathlib import Path

import numpy as np
import pandas as pd
from numpy import cos, deg2rad, sin, sqrt
from scipy.interpolate import interp1d

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
                    da_alt: xr.DataArray, da_val: xr.DataArray) -> pd.Series:
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


def download_era5(
        lat_x: float, lon_x: float,
        start: str, end: str,
        grid=[0.125, 0.125], scale=0)->dict:
    """Generate ERA-5 cdsapi-based requests and download data for area of interests.

    Parameters
    ----------
    lat_x : float
        Latitude of centre at the area of interest.
    lon_x : float
        Longitude of centre at the area of interest.
    start : str
        [description]
    end : str
        [description]
    grid : list, optional
        [description], by default [0.125, 0.125]
    scale : int, optional
        [description], by default 0

    Returns
    -------
    dict
        [description]
    """

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
