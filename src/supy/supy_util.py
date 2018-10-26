# supy utilities
import os

import numpy as np
import pandas as pd


# locate the first position of period with in-between gaps
def loc_gap(ser_test, freq='1D', pattern='010'):
    rsmp = ser_test.resample(freq)
    ser_TF_10 = rsmp.apply(lambda ser: ser.isna().any()) * 1
    str_TF_10 = ser_TF_10.astype(str).str.cat()
    pos_gap = str_TF_10.find(pattern)
    loc_ser = ser_TF_10.iloc[pos_gap:pos_gap + len(pattern)].index
    return loc_ser


# fill gap with neighboring days
def fill_gap_one(ser_test, freq='1D', pattern='010'):
    # resmaple into daily periods
    rsmp = ser_test.resample(freq)
    # locate the gaps according to gap pattern: 0 for NO gap, 1 for gapped
    loc_ser = loc_gap(ser_test, freq, pattern)

    # generator groups
    ser_find = (rsmp.get_group(x) for x in loc_ser)
    if len(loc_ser) == 0:
        return ser_test
    # print len(loc_ser), loc_ser

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
