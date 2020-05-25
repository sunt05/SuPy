import numpy as np

from .._env import logger_supy
from ._atm import cal_cp

# saturation vapour pressure [hPa]
def cal_vap_sat(Temp_C, Press_hPa):
    # temp_c= 0.001 if np.abs(Temp_C)<0.001 else Temp_C
    Press_kPa = Press_hPa / 10

    if 0.001000 <= Temp_C < 50:
        e_mb = 6.1121 * np.exp(((18.678 - Temp_C / 234.5) * Temp_C) / (Temp_C + 257.14))
        f = 1.00072 + Press_kPa * (3.2e-6 + 5.9e-10 * Temp_C ** 2)
        es_hPa = e_mb * f

    if -40 < Temp_C <= -0.001000:
        e_mb = 6.1115 * np.exp(((23.036 - Temp_C / 333.7) * Temp_C) / (Temp_C + 279.82))
        f = 1.00022 + Press_kPa * (3.83e-6 + 6.4e-10 * Temp_C ** 2)
        es_hPa = e_mb * f

    if -0.001 < Temp_C < 0.001:
        es_hPa = cal_vap_sat(0.001, Press_hPa)

    return es_hPa


# density of dry air [kg m-3]
def cal_dens_dry(RH_pct, Temp_C, Press_hPa):
    gas_ct_dry = 8.31451 / 0.028965  # dry_gas/molar
    es_hPa = cal_vap_sat(Temp_C, Press_hPa)
    Ea_hPa = RH_pct / 100 * es_hPa
    dens_dry = ((Press_hPa - Ea_hPa) * 100) / (gas_ct_dry * (273.16 + Temp_C))
    return dens_dry


# density of vapour [kg m-3]
def cal_dens_vap(RH_pct, Temp_C, Press_hPa):
    gas_ct_wv = 8.31451 / 0.0180153  # dry_gas/molar_wat_vap
    es_hPa = cal_vap_sat(Temp_C, Press_hPa)
    Ea_hPa = RH_pct / 100 * es_hPa
    vap_dens = Ea_hPa * 100 / ((Temp_C + 273.16) * gas_ct_wv)
    return vap_dens


#
# # specific heat capacity of air mass [J kg-1 K-1]
# def cal_cpa(Temp_C, RH_pct, Press_hPa):
#     # heat capacity of dry air depending on air temperature
#     cpd = 1005.0 + ((Temp_C + 23.16) ** 2) / 3364.0
#     # heat capacity of vapour
#     cpm = (
#         1859
#         + 0.13 * RH_pct
#         + (19.3 + 0.569 * RH_pct) * (Temp_C / 100.0)
#         + (10.0 + 0.5 * RH_pct) * (Temp_C / 100.0) ** 2
#     )
#
#     # density of dry air
#     rho_d = cal_dens_dry(RH_pct, Temp_C, Press_hPa)
#
#     # density of vapour
#     rho_v = cal_dens_vap(RH_pct, Temp_C, Press_hPa)
#
#     # specific heat
#     cpa = cpd * (rho_d / (rho_d + rho_v)) + cpm * (rho_v / (rho_d + rho_v))
#     return cpa


# air density [kg m-3]
def cal_dens_air(Press_hPa, Temp_C):
    # dry_gas/molar
    gas_ct_dry = 8.31451 / 0.028965

    # air density [kg m-3]
    dens_air = (Press_hPa * 100) / (gas_ct_dry * (Temp_C + 273.16))
    return dens_air


# Obukhov length
def cal_Lob(QH, UStar, Temp_C, RH_pct, Press_hPa, g=9.8, k=0.4):
    # gravity constant/(Temperature*Von Karman Constant)
    G_T_K = (g / (Temp_C + 273.16)) * k

    # air density [kg m-3]
    rho = cal_dens_air(Press_hPa, Temp_C)

    # specific heat capacity of air mass [J kg-1 K-1]
    cpa = cal_cp(Temp_C, RH_pct, Press_hPa)

    # Kinematic sensible heat flux [K m s-1]
    H = QH / (rho * cpa)

    # temperature scale
    uStar = np.max([0.01, UStar])
    TStar = -H / uStar

    # Obukhov length
    Lob = (uStar ** 2) / (G_T_K * TStar)

    return Lob


def cal_neutral(df_val, z_meas, h_sfc):
    """ Calculates the rows associated with neutral condition (threshold=0.01)


    Parameters
    ----------
    df_val: pd.DataFrame
        Index should be time with columns: 'H', 'USTAR', 'TA', 'RH', 'PA', 'WS'
    z_meas
        measurement height in m
    h_sfc
        vegetation height in m

    Returns
    -------
    ser_ws: pd.series
        observation time series of WS (Neutral conditions)
    ser_ustar: pd.series
        observation time series of u* (Neutral conditions)
    """

    # calculate Obukhov length
    ser_Lob = df_val.apply(
        lambda ser: cal_Lob(ser.H, ser.USTAR, ser.TA, ser.RH, ser.PA * 10), axis=1
    )

    # zero-plane displacement: estimated using rule f thumb `d=0.7*h_sfc`

    z_d = 0.7 * h_sfc

    if z_d >= z_meas:
        logger_supy.exception(
            "vegetation height is greater than measuring height. Please fix this before continuing . . ."
        )

    # calculate stability scale
    ser_zL = (z_meas - z_d) / ser_Lob

    # determine periods under quasi-neutral conditions
    limit_neutral = 0.01
    ind_neutral = ser_zL.between(-limit_neutral, limit_neutral)

    ind_neutral = ind_neutral[ind_neutral]

    df_sel = df_val.loc[ind_neutral.index, ["WS", "USTAR"]].dropna()
    ser_ustar = df_sel.USTAR
    ser_ws = df_sel.WS

    return ser_ws, ser_ustar


# Optimization for calculating z0 and d
def optimize_MO(df_val, z_meas, h_sfc):
    """Calculates surface roughness and zero plane displacement height.
    Refer to https://suews-parameters-docs.readthedocs.io/en/latest/steps/roughness-SuPy.html for example

    Parameters
    ----------
    df_val: pd.DataFrame
        Index should be time with columns: 'H', 'USTAR', 'TA', 'RH', 'PA', 'WS'
    z_meas
        measurement height in m
    h_sfc
        vegetation height in m

    Returns
    -------
    z0
        surface roughness
    d
        zero displacement height
    ser_ws: pd.series
        observation time series of WS (Neutral conditions)
    ser_ustar: pd.series
        observation time series of u* (Neutral conditions)
    """

    from platypus.core import Problem
    from platypus.types import Real, random
    from platypus.algorithms import NSGAIII

    # Calculates rows related to neutral conditions
    ser_ws, ser_ustar = cal_neutral(df_val, z_meas, h_sfc)

    # function to optimize
    def func_uz(params):
        z0 = params[0]
        d = params[1]
        z = z_meas
        k = 0.4
        uz = (ser_ustar / k) * np.log((z - d) / z0)  # logarithmic law

        o1 = abs(1 - np.std(uz) / np.std(ser_ws))  # objective 1: normalized STD
        # objective 2: normalized MAE
        o2 = np.mean(abs(uz - ser_ws)) / (np.mean(ser_ws))

        return [o1, o2], [uz.min(), d - z0]

    problem = Problem(2, 2, 2)
    problem.types[0] = Real(0, 10)  # bounds for first parameter (z0)
    problem.types[1] = Real(0, h_sfc)  # bounds for second parameter (zd)

    problem.constraints[0] = ">=0"  # constrain for first parameter
    problem.constraints[1] = ">=0"  # constrain for second parameter

    problem.function = func_uz
    random.seed(12345)
    algorithm = NSGAIII(problem, divisions_outer=50)
    algorithm.run(30000)

    z0s = []
    ds = []
    os1 = []
    os2 = []
    # getting the solution vaiables
    for s in algorithm.result:
        z0s.append(s.variables[0])
        ds.append(s.variables[1])
        os1.append(s.objectives[0])
        os2.append(s.objectives[1])
    # getting the solution associated with minimum obj2 (can be changed)
    idx = os2.index(min(os2, key=lambda x: abs(x - np.mean(os2))))
    z0 = z0s[idx]
    d = ds[idx]

    return z0, d, ser_ws, ser_ustar
