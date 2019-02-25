# %%
from typing import Tuple
from .supy_post import resample_output
from .supy_load import load_SUEWS_dict_ModConfig
import os
import pandas as pd
from pathlib import Path

# %%


def gen_df_save(df_grid_group: pd.DataFrame)->pd.DataFrame:
    '''generate a dataframe for saving

    Parameters
    ----------
    df_output_grid_group : pd.DataFrame
        an output dataframe of a single group and grid

    Returns
    -------
    pd.DataFrame
        a dataframe with date time info prepended for saving
    '''
    # generate df_datetime for prepending
    idx_dt = df_grid_group.index
    ser_year = pd.Series(idx_dt.year, index=idx_dt, name='Year')
    ser_DOY = pd.Series(idx_dt.dayofyear, index=idx_dt, name='DOY')
    ser_hour = pd.Series(idx_dt.hour, index=idx_dt, name='Hour')
    ser_min = pd.Series(idx_dt.minute, index=idx_dt, name='Min')
    df_datetime = pd.concat([
        ser_year,
        ser_DOY,
        ser_hour,
        ser_min,
    ], axis=1)
    df_datetime['Dectime'] = ser_DOY-1+idx_dt.to_perioddelta(
        'd').total_seconds()/(24*60*60)
    df_save = pd.concat([df_datetime, df_grid_group], axis=1)
    return df_save



def format_df_save(df_save):
    # format datetime columns
    for var in df_save.columns[:4]:
        width_var_name = max([3, len(var)])
        df_save[var] = df_save[var].map(
            lambda s: '{s:{c}>{n}}'.format(s=s, n=width_var_name, c=' '))
    # df_save.Year = df_save.Year.map(
    #     lambda s: '{s:{c}>{n}}'.format(s=s, n=4, c=' ')
    # )
    # df_save.DOY = df_save.DOY.map(
    #     lambda s: '{s:{c}>{n}}'.format(s=s, n=3, c=' ')
    # )
    # df_save.Hour = df_save.DOY.map(
    #     lambda s: '{s:{c}>{n}}'.format(s=s, n=4, c=' ')
    # )
    # df_save.Min = df_save.Min.map(
    #     lambda s: '{s:{c}>{n}}'.format(s=s, n=6, c=' ')
    # )
    df_save.Dectime = df_save.Dectime.map(
        lambda s: '{s:{c}>{n}.4f}'.format(s=s, n=8, c=' ')
    )
    # fill nan values
    df_save = df_save.fillna(-999.)
    # format value columns

    for var in df_save.columns[5:]:
        width_var_name = max([8, len(var)])
        df_save[var] = df_save[var].map(
            lambda s: '{s:{c}>{n}.2f}'.format(s=s, n=width_var_name, c=' '))

    # format column names
    col_fmt = df_save.columns.to_series()
    col_fmt[4:] = col_fmt[4:].map(
        lambda s: '{s:{c}>{n}}'.format(s=s, n=max([8, len(s)]), c=' ')
    )
    df_save.columns = col_fmt

    return df_save


# df_save = format_df_save(df_save)
# df_save.iloc[0]
# df_save.iloc[0].str.len()


def save_df_grid_group(df_grid_group, grid, group, site='test', dir_save='.'):
    # processing path
    path_dir = Path(dir_save)

    # output frequency in min
    freq_min = int(df_grid_group.index.freq.delta.total_seconds()/60)
    # staring year
    year = df_grid_group.index[0].year
    # sample file name: 'Kc98_2012_SUEWS_60.txt'
    file_out = '{site}{grid}_{year}_{group}_{freq_min}.txt'.format(
        site=site, grid=grid, year=year, group=group, freq_min=freq_min)
    # 'DailyState_1440' will be trimmed
    file_out = file_out.replace('DailyState_1440', 'DailyState')
    path_out = path_dir/file_out
    print('writing out: {path_out}'.format(path_out=path_out))
    import time
    t_start = time.time()
    # generate df_save with datetime info prepended to each row
    df_save = gen_df_save(df_grid_group)
    t_end = time.time()
    # print(t_end-t_start)

    t_start = time.time()
    # format df_save with right-justified view
    df_save = format_df_save(df_save)
    t_end = time.time()
    # print(t_end-t_start)

    t_start = time.time()
    # save to txt file
    df_save.to_csv(
        path_out,
        index=False,
        sep='\t',
    )
    t_end = time.time()
    # print(t_end-t_start)
    return path_out


# save output files
def save_df_output(
        df_output: pd.DataFrame,
        freq_s: int = 3600,
        site: str = '',
        path_dir_save: Path = Path('.'),)->list:
    '''save supy output dataframe to txt files

    Parameters
    ----------
    df_output : pd.DataFrame
        output dataframe of supy simulation
    freq_s : int, optional
        output frequency in second (the default is 3600, which indicates the a txt with hourly values)
    path_dir_save : Path, optional
        directory to save txt files (the default is '.', which the current working directory)
    site : str, optional
        site code used for filename (the default is '', which indicates no site name prepended to the filename)
    path_runcontrol : str or anything that can be parsed as `Path`, optional
        path to SUEWS 'RunControl.nml' file (the default is None, which indicates necessary saving options should be specified via other parameters)

    Returns
    -------
    list
        a list of `Path` objects for saved txt files
    '''

    list_path_save = []
    list_group = df_output.columns.get_level_values('group').unique()
    list_grid = df_output.index.get_level_values('grid').unique()
    for grid in list_grid:
        for group in list_group:
            df_output_grid_group = df_output\
                .loc[grid, group]\
                .dropna(how='all', axis=0)
            # save output at the runtime frequency (usually 5 min)
            # 'DailyState' group will be save a daily frequency
            path_save = save_df_grid_group(
                df_output_grid_group, grid, group,
                site=site, dir_save=path_dir_save)
            list_path_save.append(path_save)

    # resample output if freq_s is different from runtime freq (usually 5 min)
    freq_save = pd.Timedelta(freq_s, 'second')
    # 'DailyState' group will be dropped in `resample_output` as resampling is not needed
    df_rsmp = resample_output(df_output, freq_save)

    list_group = df_rsmp.columns.get_level_values('group').unique()
    list_grid = df_rsmp.index.get_level_values('grid').unique()
    # save output at the resampling frequency
    for grid in list_grid:
        for group in list_group:
            df_output_grid_group = df_rsmp.loc[grid, group]
            path_save = save_df_grid_group(
                df_output_grid_group, grid, group,
                site=site, dir_save=path_dir_save)
            list_path_save.append(path_save)

    return list_path_save


# %%
# save model state for restart runs
def save_df_state(
        df_state: pd.DataFrame,
        site: str = '',
        path_dir_save: Path = Path('.'),)->Path:
    '''save `df_state` to a csv file

    Parameters
    ----------
    df_state : pd.DataFrame
        a dataframe of model states produced by a supy run
    site : str, optional
        site identifier (the default is '', which indicates an empty site code)
    path_dir_save : Path, optional
        path to directory to save results (the default is Path('.'), which the current working directory)

    Returns
    -------
    Path
        path to the saved csv file
    '''

    file_state_save = 'df_state_{site}.csv'.format(site=site)
    # trim filename if site == ''
    file_state_save = file_state_save.replace('_.csv', '.csv')
    path_state_save = path_dir_save/file_state_save
    df_state_final.to_csv(path_state_save)
    return path_state_save


# %%
# get information for saving results
def get_save_info(path_runcontrol: str)->Tuple[int, Path, str]:
    '''get necessary information for saving supy results, which are (freq_s, dir_save, site)

    Parameters
    ----------
    path_runcontrol : Path
        Path to SUEWS :ref:`RunControl.nml <suews:RunControl.nml>`

    Returns
    -------
    tuple
        A tuple including (freq_s, dir_save, site):
        freq_s: output frequency in seconds
        dir_save: directory name to save results
        site: site identifier
    '''

    try:
        path_runcontrol = Path(path_runcontrol).expanduser().resolve()
    except FileNotFoundError:
        print('{path} does not exists!'.format(path=path_runcontrol))
    else:
        dict_mod_cfg = load_SUEWS_dict_ModConfig(path_runcontrol)
        freq_s, dir_save, site = [
            dict_mod_cfg[x] for x in
            [
                'resolutionfilesout',
                'fileoutputpath',
                'filecode',
            ]
        ]
        dir_save = path_runcontrol.parent/dir_save
        if not dir_save.exists():
            dir_save.mkdir()
        return freq_s, dir_save, site
