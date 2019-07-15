from scipy.optimize import least_squares
import numpy as np
import pandas as pd
from atmosp import calculate as ac

# atmospheric related utilities


def cal_des_dta(ta, pa, dta=1.0):
    """Calculate slope of es(Ta), i.e., saturation evaporation pressure `es` as function of air temperature `ta [K]`

    Parameters
    ----------
    ta : numeric
        Air temperature [K]
    pa : numeric
        Air pressure [Pa]
    dta : float, optional
        change in ta for calculating that in es, by default 1.0 K
    """

    des = ac('es', p=pa, T=ta + dta/2) - ac('es', p=pa, T=ta - dta/2)
    des_dta = des/dta
    try:
        # try to pack as Series
        des_dta = pd.Series(des_dta, index=ta.index)
    except AttributeError as ex:
        print(ex, 'cannot pack into pd.Series')
        pass
    return des_dta


def cal_rs_obs(qh, qe, ta, rh, pa):
    """Calculate surface resistance based on observations, notably turbulent fluxes.

    Parameters
    ----------
    qh : numeric
        sensible heat flux [W m-2]
    qe : numeric
        latent heat flux [W m-2]
    ta : numeric
        air temperature [K]
    rh : numeric
        relative humidity [%]
    pa : numeric
        air pressure [Pa]

    Returns
    -------
    numeric
        Surface resistance based on observations [mm s-1]
    """

    # surface resistance at water surface [s m-1]
    rav = 50

    # psychrometric constant [Pa K-1] as a function of air pressure
    ser_gamma = 0.665e-3 * pa

    # air density [kg m-3]
    val_rho = 1.27

    # heat capacity of air [J kg-1 K-1]
    val_cp = 1005

    # slope of es(Ta) curve at Ta
    ser_des_dTa = cal_des_dta(ta, pa, dta=1.0)
    #
    arr_e = ac('e', p=pa, T=ta, RH=rh)
    arr_es = ac('es', p=pa, T=ta)
    arr_vpd = arr_es-arr_e
    #
    ser_rs_1 = (ser_des_dTa / ser_gamma) * (qh / qe - 1) * rav
    ser_rs_2 = (val_rho * val_cp * arr_vpd / (ser_gamma * qe))
    ser_rs = ser_rs_1 + ser_rs_2

    try:
        # try to pack as Series
        ser_rs = pd.Series(ser_rs, index=ta.index)
    except AttributeError as ex:
        print(ex, 'cannot pack into pd.Series')
        pass

    return ser_rs


def cal_gs_obs(qh, qe, ta, rh, pa):
    """Calculate surface conductance based on observations, notably turbulent fluxes.

    Parameters
    ----------
    qh : numeric
        Sensible heat flux [W m-2]
    qe : numeric
        Latent heat flux [W m-2]
    ta : numeric
        Air temperature [K]
    rh : numeric
        Relative humidity [%]
    pa : numeric
        Air pressure [Pa]

    Returns
    -------
    numeric
        Surface conductance based on observations [mm s-1]
    """
    rs_obs = cal_rs_obs(qh, qe, ta, rh, pa)
    gs_obs = 1/rs_obs
    return gs_obs


def cal_g_lai(lai, g1, lai_max):
    """Calculate LAI-related correction coefficient for surface conductance.

    Parameters
    ----------
    lai : numeric
        Leaf area index [m2 m-2]
    g1 : numeric
        LAI-related correction parameter [-]
    lai_max : numeric
        Maximum LAI [m2 m-2]

    Returns
    -------
    numeric
        LAI-related correction coefficient [-]
    """
    g_lai = lai/lai_max*g1
    return g_lai


def cal_g_kd(kd, g2, kd_max=1200.):
    """Calculate solar radiation-related correction coefficient for surface conductance.

    Parameters
    ----------
    kd : numeric
        Incoming solar radiation [W m-2]
    g2 : numeric
        Solar radiation-related correction parameter [-]
    kd_max : numeric, optional
        Maximum incoming solar radiation [W m-2], by default 1200.

    Returns
    -------
    numeric
        Solar radiation-related correction coefficient [-]
    """
    g_kd_nom = kd/(g2+kd)
    g_kd_denom = kd_max/(g2+kd_max)
    g_kd = g_kd_nom/g_kd_denom
    return g_kd


def cal_g_dq(dq, g3, g4):
    """Calculate air humidity-related correction coefficient for surface conductance.

    Parameters
    ----------
    dq : numeric
        Specific humidity deficit [k kg-1]
    g3 : numeric
        Specific humidity-related correction parameter [-]
    g4 : numeric
        Specific humidity-related correction parameter [-]

    Returns
    -------
    numeric
        Air humidity-related correction coefficient
    """
    g_dq = g3+(1-g3)*g4**dq
    return g_dq


def cal_g_ta(ta_c, g5, tl=-10., th=55.):
    """Calculate air temperature-related correction coefficient for surface conductance.

    Parameters
    ----------
    ta_c : numeric
        Air temperature [degC]
    g5 : numeric
        Air temperature-related correction parameter
    tl : numeric, optional
        Low temperature limit [degC], by default -10.
    th : numeric, optional
        High temperature limit [degC], by default 55.

    Returns
    -------
    numeric
        Air temperature-related correction coefficient
    """


    tc = (th-g5)/(g5-tl)
    g_ta = ((ta_c-tl)*(th-ta_c)**tc)/((g5-tl)*(th-g5)**tc)
    return g_ta


def cal_g_smd(smd, g6, wp=120.):
    """Calculate soil moisture-related correction coefficient for surface conductance.

    Parameters
    ----------
    smd : numeric
        Soil moisture deficit [mm].
    g6 : numeric
        Soil moisture-related correction parameter.
    wp : numeric, optional
        Wilting point indicated as deficit [mm], by default 120.

    Returns
    -------
    numeric
        Soil moisture-related correction coefficient
    """
    g_smd_nom = 1-np.exp(g6*(smd-wp))
    g_smd_denom = 1-np.exp(g6*(0-wp))
    g_smd = g_smd_nom/g_smd_denom
    return g_smd


def cal_gs_mod(kd, ta_k, rh, pa, smd, lai, g_cst, g_max=30., lai_max=6.):
    """Model surface conductance/resistance using phenology and atmospheric forcing conditions.

    Parameters
    ----------
    kd : numeric
        Incoming solar radiation [W m-2]
    ta_k : numeric
        Air temperature [K]
    rh : numeric
        Relative humidity [%]
    pa : numeric
        Air pressure
    smd : numeric
        Soil moisture deficit [mm]
    lai : numeric
        Leaf area index [m2 m-2]
    g_cst : size-6 array
        Parameters to determine surface conductance/resistance:
        g1 (LAI related), g2 (solar radiation related),
        g3 (humidity related), g4 (humidity related),
        g5 (air temperature related),
        g6 (soil moisture related)
    g_max : numeric, optional
        Maximum surface conductance [mm s-1], by default 30
    lai_max : numeric, optional
        Maximum LAI [m2 m-2], by default 6

    Returns
    -------
    numeric
        Modelled surface conductance [mm s-1]
    """

    # broadcast g1 – g6
    g1, g2, g3, g4, g5, g6 = g_cst
    # print(g1, g2, g3, g4, g5, g6)
    # lai related
    g_lai = cal_g_lai(lai, g1, lai_max)
    # print('g_lai', g_lai)

    # kdown related
    g_kd = cal_g_kd(kd, g2)
    # print('g_kd', g_kd)
    # dq related
    # ta_k = ta_c+273.15
    dq = ac('qvs', T=ta_k, p=pa)-ac('qv', T=ta_k, p=pa, RH=rh)
    g_dq = cal_g_dq(dq, g3, g4)
    # print('g_dq', g_dq)
    # ta related
    ta_c = ta_k - 273.15
    g_ta = cal_g_ta(ta_c, g5)
    # print('g_ta', g_ta)
    # smd related
    g_smd = cal_g_smd(smd, g6)
    # print('g_smd', g_smd)
    # combine all corrections
    gs_c = g_lai*g_kd*g_dq*g_ta*g_smd
    gs = g_max*gs_c

    return gs


def calib_g(df_fc_suews, g_max=30., lai_max=6., debug=False):
    """Calibrate parameters for modelling surface conductance over vegetated surfaces.

    Parameters
    ----------
    df_fc_suews : pandas.DataFrame
        DataFrame in `SuPy forcing <https://supy.readthedocs.io/en/latest/data-structure/df_forcing.html>`_ format
    g_max : numeric, optional
        Maximum surface conductance [mm s-1], by default 30
    lai_max : numeric, optional
        Maximum LAI [m2 m-2], by default 6
    debug : bool, optional
        Option to output final calibrated `scipy.least_squares`, by default False

    Returns
    -------
    `numpy.array`, or `scipy.least_squares` if `debug==True`
        Calibrated parameters for surface conductance:
        g1 (LAI related), g2 (solar radiation related),
        g3 (humidity related), g4 (humidity related),
        g5 (air temperature related),
        g6 (soil moisture related)
    """

    df_obs = df_fc_suews.copy()
    df_obs.pres *= 100
    df_obs.Tair += 273.15

    # function for least_square optimiser
    def cal_prm_g(prm_g):
        gs_obs = cal_gs_obs(df_obs.qh, df_obs.qe, df_obs.Tair,
                            df_obs.RH, df_obs.pres)
        gs_mod = cal_gs_mod(df_obs.kdown, df_obs.Tair,
                            df_obs.RH, df_obs.pres, df_obs.xsmd,
                            df_obs.lai, prm_g, g_max, lai_max)
        resid_gs = gs_obs - gs_mod
        resid_gs = resid_gs.dropna()
        return resid_gs

    # initial guess
    prm_g_0 = [3.5, 200.0, 0.13, 0.7, 30.0, 0.05]
    # calibrated model
    res_ls = least_squares(cal_prm_g,
                           prm_g_0,
                           bounds=([0, 0, 0, 0, 0, 0], np.inf))

    res = res_ls if debug else res_ls.x

    return res
