

# from scipy.optimize import least_squares
import numpy as np
import pandas as pd



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
    from atmosp import calculate as ac
    des = ac("es", p=pa, T=ta + dta / 2) - ac("es", p=pa, T=ta - dta / 2)
    des_dta = des / dta
    try:
        # try to pack as Series
        des_dta = pd.Series(des_dta, index=ta.index)
    except AttributeError as ex:
        print(ex, "cannot pack into pd.Series")
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
        Surface resistance based on observations [s m-1]
    """
    from atmosp import calculate as ac
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
    arr_e = ac("e", p=pa, T=ta, RH=rh)
    arr_es = ac("es", p=pa, T=ta)
    arr_vpd = arr_es - arr_e
    #
    ser_rs_1 = (ser_des_dTa / ser_gamma) * (qh / qe - 1) * rav
    ser_rs_2 = val_rho * val_cp * arr_vpd / (ser_gamma * qe)
    ser_rs = ser_rs_1 + ser_rs_2

    try:
        # try to pack as Series
        ser_rs = pd.Series(ser_rs, index=ta.index)
    except AttributeError as ex:
        print(ex, "cannot pack into pd.Series")
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
    gs_obs = 1e3 / rs_obs
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
    g_lai = lai / lai_max * g1
    return g_lai


def cal_g_kd(kd, g2, kd_max=1200.0):
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
    g_kd_nom = kd / (g2 + kd)
    g_kd_denom = kd_max / (g2 + kd_max)
    g_kd = g_kd_nom / g_kd_denom
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
    g_dq = g3 + (1 - g3) * g4 ** dq
    return g_dq


def cal_g_ta(ta_c, g5, tl=-10.0, th=55.0):
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

    tc = (th - g5) / (g5 - tl)
    # set a threshold for avoiding numeric difficulty
    tc = np.min([tc, 20])
    # g_ta = ((ta_c-tl)*(th-ta_c)**tc)/((g5-tl)*(th-g5)**tc)
    g_ta_nom = (ta_c - tl) * np.power((th - ta_c), tc)
    g_ta_denom = (g5 - tl) * np.power((th - g5), tc)
    g_ta = g_ta_nom / g_ta_denom

    return g_ta


def cal_g_smd(smd, g6, s1=5.56):
    """Calculate soil moisture-related correction coefficient for surface conductance.

    Parameters
    ----------
    smd : numeric
        Soil moisture deficit [mm].
    g6 : numeric
        Soil moisture-related correction parameter.
    s1 : numeric, optional
        Wilting point (WP=s1/g6, indicated as deficit [mm]) related parameter, by default 5.56

    Returns
    -------
    numeric
        Soil moisture-related correction coefficient
    """
    # Wilting point calculated following SUEWS
    wp = s1 / g6

    g_smd_nom = 1 - np.exp(g6 * (smd - wp))
    g_smd_denom = 1 - np.exp(g6 * (0 - wp))
    g_smd = g_smd_nom / g_smd_denom
    return g_smd


def cal_gs_mod(kd, ta_k, rh, pa, smd, lai, g_cst, g_max=30.0, lai_max=6.0, s1=5.56):
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
    s1 : numeric, optional
        Wilting point (WP=s1/g6, indicated as deficit [mm]) related parameter, by default 5.56


    Returns
    -------
    numeric
        Modelled surface conductance [mm s-1]
    """
    from atmosp import calculate as ac
    # broadcast g1 – g6
    # print('g_cst', g_cst)
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
    dq = ac("qvs", T=ta_k, p=pa) - ac("qv", T=ta_k, p=pa, RH=rh)
    g_dq = cal_g_dq(dq, g3, g4)
    # print('g_dq', g_dq)
    # ta related
    ta_c = ta_k - 273.15
    g_ta = cal_g_ta(ta_c, g5)
    # print('g_ta', g_ta)
    # smd related
    g_smd = cal_g_smd(smd, g6, s1)
    # print('g_smd', g_smd)
    # combine all corrections
    gs_c = g_lai * g_kd * g_dq * g_ta * g_smd
    gs = g_max * gs_c

    return gs


def calib_g(
    df_fc_suews,
    g_max=33.1,
    lai_max=5.9,
    s1=5.56,
    method="cobyla",
    prms_init=None,
    debug=False,
):
    """Calibrate parameters for modelling surface conductance over vegetated surfaces using `LMFIT <https://lmfit.github.io/lmfit-py/model.html>`.

    Parameters
    ----------
    df_fc_suews : pandas.DataFrame
        DataFrame in `SuPy forcing <https://supy.readthedocs.io/en/latest/data-structure/df_forcing.html>`_ format
    g_max : numeric, optional
        Maximum surface conductance [mm s-1], by default 30
    lai_max : numeric, optional
        Maximum LAI [m2 m-2], by default 6
    s1 : numeric, optional
        Wilting point (WP=s1/g6, indicated as deficit [mm]) related parameter, by default 5.56
    method: str, optional
        Method used in minimisation by `lmfit.minimize`: details refer to its `method<lmfit:minimize>`.
    prms_init: lmfit.Parameters
        Initial parameters for calibration
    debug : bool, optional
        Option to output final calibrated `ModelResult <lmfit:ModelResult>`, by default False

    Returns
    -------
    dict, or `ModelResult <lmfit:ModelResult>` if `debug==True`
        1. dict: {parameter_name -> best_fit_value}
        2. `ModelResult`

        Note:
            Parameters for surface conductance:
            g1 (LAI related), g2 (solar radiation related),
            g3 (humidity related), g4 (humidity related),
            g5 (air temperature related),
            g6 (soil moisture related)

    Note
    ----
    For calibration validity, turbulent fluxes, QH and QE, in `df_fc_suews` should ONLY be observations, i.e., interpolated values should be avoided.
    To do so, please place `np.nan` as missing values for QH and QE.

    """
    from lmfit import Model, Parameters, Parameter
    list_var_sel = ["qh", "qe", "Tair", "RH", "pres", "kdown", "xsmd", "lai"]
    df_obs = df_fc_suews[list_var_sel].copy().dropna()
    df_obs.pres *= 100
    df_obs.Tair += 273.15

    gs_obs = cal_gs_obs(df_obs.qh, df_obs.qe, df_obs.Tair, df_obs.RH, df_obs.pres)

    def func_fit_g(kd, ta, rh, pa, smd, lai, g1, g2, g3, g4, g5, g6):
        return cal_gs_mod(
            kd, ta, rh, pa, smd, lai, [g1, g2, g3, g4, g5, g6], g_max, lai_max, s1
        )

    gmodel = Model(
        func_fit_g,
        independent_vars=["lai", "kd", "ta", "rh", "pa", "smd"],
        param_names=["g1", "g2", "g3", "g4", "g5", "g6"],
    )
    if prms_init is None:
        print("Preset parameters will be loaded!")
        print("Please use with caution.")
        prms = Parameters()
        prm_g_0 = [3.5, 200.0, 0.13, 0.7, 30.0, 0.05]
        list_g = (
            Parameter(f"g{i+1}", prm_g_0[i], True, 0, None, None, None)
            for i in range(6)
        )
        prms.add_many(*list_g)
        # set specific bounds:
        # g3, g4: specific humidity related
        prms["g3"].set(min=0, max=1)
        prms["g4"].set(min=0, max=1)
        # g5: within reasonable temperature ranges
        prms["g5"].set(min=-10, max=55)
        # g6: within sensitive ranges of SMD
        prms["g6"].set(min=0.02, max=0.1)
    else:
        print("User provided parameters are loaded!")
        prms = prms_init

    # pack into a DataFrame for filtering out nan
    df_fit = pd.concat([gs_obs.rename("gs_obs"), df_obs], axis=1).dropna()

    res_fit = gmodel.fit(
        df_fit.gs_obs,
        kd=df_fit.kdown,
        ta=df_fit.Tair,
        rh=df_fit.RH,
        pa=df_fit.pres,
        smd=df_fit.xsmd,
        lai=df_fit.lai,
        params=prms,
        # useful ones: ['nelder', 'powell', 'cg', 'cobyla', 'bfgs', 'trust-tnc']
        method=method,
        #     nan_policy='omit',
        verbose=True,
    )

    # provide full fitted model if debug == True otherwise only a dict with best fit parameters
    res = res_fit if debug else res_fit.best_values

    return res


# calculate specific humidity using relative humidity
def cal_qa(rh_pct, theta_K, pres_hPa):
    from atmosp import calculate as ac
    qa = ac("qv", RH=rh_pct, p=pres_hPa * 100, theta=theta_K)
    return qa


# calculate relative humidity using specific humidity
def cal_rh(qa_kgkg, theta_K, pres_hPa):
    from atmosp import calculate as ac
    RH = ac("RH", av=qa_kgkg, p=pres_hPa * 100, theta=theta_K)
    return RH


# calculate latent heat of vaporisation
def cal_lat_vap(qa_kgkg, theta_K, pres_hPa):
    from atmosp import calculate as ac
    # wel-bulb temperature
    tw = ac(
        "Tw", qv=qa_kgkg, p=pres_hPa, theta=theta_K, remove_assumptions=("constant Lv")
    )
    # latent heat [J kg-1]
    Lv = 2.501e6 - 2370.0 * (tw - 273.15)
    return Lv


# calculate heat capacity of air
def cal_cp(qa_kgkg, theta_K, pres_hPa):
    from atmosp import calculate as ac

    temp_C = ac("T", theta=theta_K, p=pres_hPa * 100) - 273.15

    rh_pct = ac("RH", qv=qa_kgkg, theta=theta_K, p=pres_hPa * 100)

    # Garratt equation a20(1992)
    cpd = 1005.0 + ((temp_C + 23.16) ** 2) / 3364.0

    # Beer(1990) for water vapour
    cpm = (
        1859
        + 0.13 * rh_pct
        + (19.3 + 0.569 * rh_pct) * (temp_C / 100.0)
        + (10.0 + 0.5 * rh_pct) * (temp_C / 100.0) ** 2
    )

    # air density
    rho = ac("rho", qv=qa_kgkg, theta=theta_K, p=pres_hPa * 100)

    # water vapour mixing ratio
    rv = ac("rv", qv=qa_kgkg, theta=theta_K, p=pres_hPa * 100)

    # dry air density
    rho_d = rv / (1 + rv) * rho

    # water vapour density
    rho_v = rho - rho_d

    # heat capacity of air
    cp = cpd * (rho_d / (rho_d + rho_v)) + cpm * (rho_v / (rho_d + rho_v))

    return cp


# stability correction for momentum
def cal_psi_mom(zoL):
    # limit for neutral condition
    lim_ntrl = 1e-5

    zoL = np.where(np.abs(zoL) > 5, 5 * np.sign(zoL), zoL)

    # stable, zoL>0
    zoL_stab = np.where(zoL > lim_ntrl, zoL, 0)
    psim_stab = (-6) * np.log(1 + zoL_stab)

    # unstable, zoL<0
    zoL_unstab = np.where(zoL < -lim_ntrl, zoL, 0)
    psim_unstab = 0.6 * (2) * np.log((1 + (1 - 16 * zoL_unstab) ** 0.5) / 2)

    # populate values with respect to stability
    psim = np.where(zoL > lim_ntrl, psim_stab, psim_unstab)
    psim = np.where(np.abs(zoL) <= lim_ntrl, 0, psim)

    return psim


# stability correction for heat
def cal_psi_heat(zoL):
    # limit for neutral condition
    lim_ntrl = 1e-5

    zoL = np.where(np.abs(zoL) > 5, 5 * np.sign(zoL), zoL)

    # stable, zoL>0
    zoL_stab = np.where(zoL > lim_ntrl, zoL, 0)
    psih_stab = -4.5 * zoL_stab

    # unstable, zoL<0
    zoL_unstab = np.where(zoL < -lim_ntrl, zoL, 0)
    psih_unstab = (2) * np.log((1 + (1 - 16 * zoL_unstab) ** 0.5) / 2)

    # populate values with respect to stability
    psih = np.where(zoL > lim_ntrl, psih_stab, psih_unstab)
    psih = np.where(np.abs(zoL) <= lim_ntrl, 0, psih)

    return psih
