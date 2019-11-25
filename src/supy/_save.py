import logging
import f90nml
from typing import Tuple
from ._post import resample_output
from ._load import load_SUEWS_dict_ModConfig, dict_InitCond_out
import os
import pandas as pd
from pathlib import Path
from ._env import logger_supy


def gen_df_save(df_grid_group: pd.DataFrame) -> pd.DataFrame:
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
    df_datetime['Dectime'] = ser_DOY - 1 + idx_dt.to_perioddelta(
        'd').total_seconds() / (24 * 60 * 60)
    df_save = pd.concat([df_datetime, df_grid_group], axis=1)
    return df_save


def format_df_save(df_save):
    # format datetime columns
    for var in df_save.columns[:4]:
        width_var_name = max([3, len(var)])
        df_save[var] = df_save[var].map(
            lambda s: '{s:{c}>{n}}'.format(s=s, n=width_var_name, c=' '))

    df_save.Dectime = df_save.Dectime.map(
        lambda s: '{s:{c}>{n}.4f}'.format(s=s, n=8, c=' '))
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
        lambda s: '{s:{c}>{n}}'.format(s=s, n=max([8, len(s)]), c=' '))
    df_save.columns = col_fmt

    return df_save


def save_df_grid_group(df_grid_group, grid, group, site='test', dir_save='.'):
    # processing path
    path_dir = Path(dir_save)

    # pandas bug here: monotonic datetime index would lose `freq` once `pd.concat`ed
    if df_grid_group.shape[0] > 0 and df_grid_group.index.size > 2:
        ind = df_grid_group.index
        freq_cal = ind[1] - ind[0]
        df_grid_group = df_grid_group.asfreq(freq_cal)
    else:
        df_grid_group = df_grid_group.asfreq('5T')
    # output frequency in min
    freq_min = int(df_grid_group.index.freq.delta.total_seconds() / 60)
    # staring year
    year = df_grid_group.index[0].year
    # sample file name: 'Kc98_2012_SUEWS_60.txt'
    file_out = '{site}{grid}_{year}_{group}_{freq_min}.txt'.format(
        site=site, grid=grid, year=year, group=group, freq_min=freq_min)
    # 'DailyState_1440' will be trimmed
    file_out = file_out.replace('DailyState_1440', 'DailyState')
    path_out = path_dir / file_out
    # print('writing out: {path_out}'.format(path_out=path_out))
    import time
    t_start = time.time()
    # generate df_save with datetime info prepended to each row
    df_save = gen_df_save(df_grid_group)
    t_end = time.time()
    logger_supy.debug(f'df_save generated in {t_end-t_start:.2f} s')

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
    logger_supy.debug(f'{path_out} saved in {t_end-t_start:.2f} s')
    return path_out


# save output files
def save_df_output(
        df_output: pd.DataFrame,
        freq_s: int = 3600,
        site: str = '',
        path_dir_save: Path = Path('.'),
        save_tstep=False,
) -> list:
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
    save_tstep : bool, optional
        whether to save results in temporal resolution as in simulation (which may result very large files and slow progress), by default False.

    Returns
    -------
    list
        a list of `Path` objects for saved txt files
    '''

    list_path_save = []
    # list_group = df_output.columns.get_level_values('group').unique()
    # list_grid = df_output.index.get_level_values('grid').unique()
    # for grid in list_grid:
    #     for group in list_group:
    #         df_output_grid_group = df_output\
    #             .loc[grid, group]\
    #             .dropna(how='all', axis=0)
    #         # save output at the runtime frequency (usually 5 min)
    #         # 'DailyState' group will be save a daily frequency
    #         path_save = save_df_grid_group(df_output_grid_group,
    #                                        grid,
    #                                        group,
    #                                        site=site,
    #                                        dir_save=path_dir_save)
    #         list_path_save.append(path_save)

    # resample output if freq_s is different from runtime freq (usually 5 min)
    freq_save = pd.Timedelta(freq_s, 's')
    # resample `df_output` at `freq_save`
    df_rsmp = resample_output(df_output, freq_save)
    # 'DailyState' group will be dropped in `resample_output` as resampling is not needed
    df_rsmp = df_rsmp.drop(columns='DailyState')

    list_df_save = [df_output, df_rsmp] if save_tstep else [df_rsmp, df_output.loc[:,['DailyState']]]
    # save output at the resampling frequency
    for df_save in list_df_save:
        list_grid = df_save.index.get_level_values('grid').unique()
        for grid in list_grid:
            list_group = df_save.columns.get_level_values('group').unique()
            for group in list_group:
                df_output_grid_group = df_save\
                    .loc[grid, group]\
                    .dropna(how='all', axis=0)
                path_save = save_df_grid_group(
                    df_output_grid_group,
                    grid,
                    group,
                    site=site,
                    dir_save=path_dir_save,
                )

                # remove freq info from `DailyState` file
                if 'DailyState' in path_save.name:
                    str_fn_dd = str(path_save).replace('DailyState_5',
                                                       'DailyState')
                    path_save.rename(Path(str_fn_dd))
                    path_save = Path(str_fn_dd)

                list_path_save.append(path_save)

    return list_path_save


# save model state for restart runs
def save_df_state(
        df_state: pd.DataFrame,
        site: str = '',
        path_dir_save: Path = Path('.'),
) -> Path:
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
    path_state_save = Path(path_dir_save) / file_state_save
    # print('writing out: {path_out}'.format(path_out=path_state_save))
    df_state.to_csv(path_state_save)
    return path_state_save


# get information for saving results
def get_save_info(path_runcontrol: str) -> Tuple[int, Path, str]:
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
        logger_supy.exception(f'{path_runcontrol} does not exists!')
    else:
        dict_mod_cfg = load_SUEWS_dict_ModConfig(path_runcontrol)
        freq_s, dir_save, site, save_tstep = [
            dict_mod_cfg[x] for x in [
                'resolutionfilesout',
                'fileoutputpath',
                'filecode',
                'keeptstepfilesout',
            ]
        ]
        dir_save = path_runcontrol.parent / dir_save
        if not dir_save.exists():
            dir_save.mkdir()
        return freq_s, dir_save, site, save_tstep


# dict for {nml_save:(df_state_var,index)}
dict_init_nml = {
    'dayssincerain': ('hdd_id', '(11,)'),
    'temp_c0': ('hdd_id', '(8,)'),
    'gdd_1_0': ('gdd_id', '(0,)'),
    'gdd_2_0': ('gdd_id', '(1,)'),
    'laiinitialevetr': ('lai_id', '(0,)'),
    'laiinitialdectr': ('lai_id', '(1,)'),
    'laiinitialgrass': ('lai_id', '(2,)'),
    'albevetr0': ('albevetr_id', '0'),
    'albdectr0': ('albdectr_id', '0'),
    'albgrass0': ('albgrass_id', '0'),
    'decidcap0': ('decidcap_id', '0'),
    'porosity0': ('porosity_id', '0'),
    'soilstorepavedstate': ('soilstore_id', '(0,)'),
    'soilstorebldgsstate': ('soilstore_id', '(1,)'),
    'soilstoreevetrstate': ('soilstore_id', '(2,)'),
    'soilstoredectrstate': ('soilstore_id', '(3,)'),
    'soilstoregrassstate': ('soilstore_id', '(4,)'),
    'soilstorebsoilstate': ('soilstore_id', '(5,)'),
    'pavedstate': ('state_id', '(0,)'),
    'bldgsstate': ('state_id', '(1,)'),
    'evetrstate': ('state_id', '(2,)'),
    'dectrstate': ('state_id', '(3,)'),
    'grassstate': ('state_id', '(4,)'),
    'bsoilstate': ('state_id', '(5,)'),
    'waterstate': ('state_id', '(6,)'),
    'snowwaterpavedstate': ('snowwater', '(0,)'),
    'snowwaterbldgsstate': ('snowwater', '(1,)'),
    'snowwaterevetrstate': ('snowwater', '(2,)'),
    'snowwaterdectrstate': ('snowwater', '(3,)'),
    'snowwatergrassstate': ('snowwater', '(4,)'),
    'snowwaterbsoilstate': ('snowwater', '(5,)'),
    'snowwaterwaterstate': ('snowwater', '(6,)'),
    'snowpackpaved': ('snowpack', '(0,)'),
    'snowpackbldgs': ('snowpack', '(1,)'),
    'snowpackevetr': ('snowpack', '(2,)'),
    'snowpackdectr': ('snowpack', '(3,)'),
    'snowpackgrass': ('snowpack', '(4,)'),
    'snowpackbsoil': ('snowpack', '(5,)'),
    'snowpackwater': ('snowpack', '(6,)'),
    'snowfracpaved': ('snowfrac', '(0,)'),
    'snowfracbldgs': ('snowfrac', '(1,)'),
    'snowfracevetr': ('snowfrac', '(2,)'),
    'snowfracdectr': ('snowfrac', '(3,)'),
    'snowfracgrass': ('snowfrac', '(4,)'),
    'snowfracbsoil': ('snowfrac', '(5,)'),
    'snowfracwater': ('snowfrac', '(6,)'),
    'snowdenspaved': ('snowdens', '(0,)'),
    'snowdensbldgs': ('snowdens', '(1,)'),
    'snowdensevetr': ('snowdens', '(2,)'),
    'snowdensdectr': ('snowdens', '(3,)'),
    'snowdensgrass': ('snowdens', '(4,)'),
    'snowdensbsoil': ('snowdens', '(5,)'),
    'snowdenswater': ('snowdens', '(6,)'),
    'snowalb0': ('snowalb', '0'),
}


# save initcond namelist as SUEWS binary
def save_initcond_nml(
        df_state: pd.DataFrame,
        site: str = '',
        path_dir_save: Path = Path('.'),
) -> Path:
    # get last time step
    try:
        tstep_last = df_state.index.levels[0].max()
    except AttributeError:
        logger_supy.exception(
            ('incorrect structure detected;' +
             ' check if `df_state` is the final model state.'))
        return

    # get year for filename formatting
    year_last = tstep_last.year
    # generate a df with records of the last tstep
    df_state_last_tstep = df_state.loc[tstep_last]
    # get grid list
    list_grid = df_state_last_tstep.index

    # list holder for paths written out in nml
    list_path_nml = []
    for grid in list_grid:
        # generate nml filename
        filename_out_grid = f'InitialConditions{site}{grid}_{year_last}_EndofRun.nml'
        # derive a save path
        path_nml = path_dir_save / filename_out_grid
        # retrieve initcond values from `df_state_last_tstep`
        nml = {
            'InitialConditions': {
                key: df_state_last_tstep.loc[grid, var]
                for key, var in dict_init_nml.items()
            }
        }
        # save nml
        f90nml.write(nml, path_nml, force=True)
        # f90nml.write(nml, nml_file,force=True)
        list_path_nml.append(path_nml)
    return list_path_nml
