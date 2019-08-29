# %% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataScience.changeDirOnImportExport setting
from ast import literal_eval
from pathlib import Path
import os
try:
    os.chdir(os.path.join(os.getcwd(), 'docs/proc_var_info'))
    print(os.getcwd())
except:
    pass

# %%
from nml_rst_proc import gen_url_option, set_input, set_site, set_initcond, set_runcontrol, set_input_initcond, set_input_runcontrol, parse_option_rst
from urlpath import URL
import urllib.request
import numpy as np
import pandas as pd
import supy as sp
import os
os.getcwd()

# %% [markdown]
# ## generate site characteristics related dataframe


# %%
print('loading in', 'gen_df_state', '...')
# list of useful URLs
url_repo_base = ('https://github.com/'
                 + 'Urban-Meteorology-Reading/'
                 + 'SUEWS/raw/master/docs/source')
url_repo_input = URL(url_repo_base)/'input_files'
url_repo_input_site = url_repo_input/'SUEWS_SiteInfo/csv-table'

# list of `SUEWS_**.csv `tables
list_table = [file.replace('.txt', '.csv')
              for file in sp._load.list_file_input
              if 'Profile' not in file]

# list of namelist based variables
list_var_nml = list(set_initcond)+list(set_runcontrol)
list_url_rst_x = list(gen_url_option(var, source='github')
                      for var in sorted(list_var_nml))
# list of URLs for rst files of namelist variables
list_url_rst = list(set(URL('/'.join(url.parts)) for url in list_url_rst_x))


# preload sample data for later use
df_init_sample, df_forcing_sample = sp.load_SampleData()

print('loading in', 'gen_df_state', 'done')
# %%
# generate df_site: site characteristics related dataframe


def extract_var_suews(dict_var_full: dict, var_supy: str)->list:
    '''extract related SUEWS variables for a supy variable `var_supy`

    Parameters
    ----------
    dict_var_full : dict
        dict_var_full = sp._load.exp_dict_full(sp._load.dict_var2SiteSelect)
    var_supy : str
        supy variable name

    Returns
    -------
    list
        related SUEWS variables for `var_supy`
    '''

    x = sp._load.flatten_list(dict_var_full[var_supy])
    x = np.unique(x)
    x = [
        xx for xx in x
        if xx not in ['base', 'const', '0.0'] + [str(x) for x in range(24)]
    ]
    x = [xx for xx in x if 'Code' not in xx]
    return x


# retrieve variable dimensionality
def gen_df_dim(df_init_sample)->pd.DataFrame:
    df_var_dim_x = df_init_sample\
        .columns\
        .to_frame()\
        .reset_index(drop=True)\
        .groupby('var')\
        .last()

    df_var_dim = df_var_dim_x\
        .applymap(lambda x: literal_eval(x))\
        .applymap(
            lambda x:
            tuple(np.array(x) + 1)
            if isinstance(x, tuple) else x)\
        .rename(columns={'ind_dim': 'Dimensionality'})
    return df_var_dim


# gen_df_dim(df_init_sample)

# %%
# generate site characteristics dataframe
def gen_df_site(
        list_csv_in=list_table,
        url_base=url_repo_input_site)->pd.DataFrame:
    '''Generate description info of supy output results as a dataframe

    Parameters
    ----------
    path_csv_out : str, optional
        path to the output csv file (the default is 'df_output.csv')
    list_csv_in : list, optional
        list of file names for csv files with meta info (the default is url_repo_input_site, which is defined at the top of this file)
    url_base : URL, optional
        URL to the input dir of repo base (the default is url_repo_input, which is defined at the top of this file)

    Returns
    -------
    pd.DataFrame
        full path to the output csv file
    '''

    # list of URLs
    list_url_table = [
        url_base/table for table in list_csv_in
    ]
    try:
        df_var_info = pd.concat([pd.read_csv(f) for f in list_url_table])

        # df_var_info = pd.concat(
        #     [pd.read_csv(f) for f in list_url_table],
        #     sort=False)
    except:
        for url in list_url_table:
            if not url.get().ok:
                print(f'{url} not existing!')
    else:
        # clean meta info
        df_var_info_x = df_var_info\
            .drop(['No.', 'Use'], axis=1)\
            .set_index('Column Name')
        df_var_info_x.index = df_var_info_x.index.map(
            lambda x: x.replace('`', ''))

        # retrieve SUEWS-related variables
        dict_var_full = sp._load.exp_dict_full(
            sp._load.dict_var2SiteSelect)
        dict_var_ref_suews = {
            k: extract_var_suews(dict_var_full, k)
            for k in dict_var_full
        }
        df_var_ref_suews = pd.DataFrame(
            {k: ', '.join(dict_var_ref_suews[k])
             for k in dict_var_ref_suews},
            index=[0]).T.rename({
                0: 'SUEWS-related variables'
            }, axis=1)

        # retrive supy variable description
        dict_var_desc = {
            k: '\n'.join(df_var_info_x.loc[v].values.flatten())
            for k, v in dict_var_ref_suews.items()
        }
        df_var_desc = pd.DataFrame(dict_var_desc, index=[0]).T\
            .rename(columns={0: 'Description'})

        # retrieve variable dimensionality
        df_var_dim = gen_df_dim(df_init_sample)

        df_var_site_raw = pd.concat(
            [df_var_dim, df_var_desc, df_var_ref_suews],
            axis=1, sort=False)

        df_var_site = df_var_site_raw.filter(items=set_input, axis=0).dropna()

        return df_var_site


# generate dataframe with info of all namelist related variables
def gen_df_nml(set_initcond, set_runcontrol):
    list_var_nml = list(set_initcond)+list(set_runcontrol)
    list_url_rst_x = list(gen_url_option(var, source='github')
                          for var in sorted(list_var_nml))
    list_url_rst = list(set(URL('/'.join(url.parts))
                            for url in list_url_rst_x))

    #
    try:
        df_var_nml = pd.concat(parse_option_rst(url) for url in list_url_rst)
    except:
        for url in list_url_rst:
            if not url.get().ok:
                print(f'{url} not existing!')
    else:
        df_var_nml_x = df_var_nml\
            .reset_index()\
            .loc[:, ['Description', 'index']]\
            .rename(columns={'index': 'SUEWS-related variables'})
        index_nml = df_var_nml_x.loc[:, 'SUEWS-related variables'].str.lower()
        df_var_nml_x.index = index_nml.rename('var')
        return df_var_nml_x


# generate dataframe with info of all runcontrol related variables
def gen_df_runcontrol(set_initcond, set_runcontrol, set_input_runcontrol):
    df_var_nml_x = gen_df_nml(set_initcond, set_runcontrol)
    df_var_runcontrol_x = df_var_nml_x.filter(
        items=set_input_runcontrol, axis=0)

    # retrieve variable dimensionality
    df_var_dim = gen_df_dim(df_init_sample)
    df_var_runcontrol_dim = df_var_dim.filter(
        items=set_input_runcontrol, axis=0)

    df_var_runcontrol = pd.concat(
        [df_var_runcontrol_dim, df_var_runcontrol_x],
        axis=1,
        sort=False)
    return df_var_runcontrol


# generate dataframe with info of all initial condition related variables
def gen_df_initcond(set_initcond, set_runcontrol):
    df_var_nml_x = gen_df_nml(set_initcond, set_runcontrol)
    dict_initcond_related = {
        'albdectr_id': 'albdectr0',
        'albevetr_id': 'albevetr0',
        'albgrass_id': 'albgrass0',
        'decidcap_id': 'decidcap0',
        # 'dqndt': 'dqndt',
        # 'dqnsdt': 'dqnsdt',
        # 'gdd_id': ['gdd_1_0', 'gdd_2_0', ],
        # 'hdd_id':,
        # 'icefrac':,
        'lai_id': [
            'laiinitialdectr',
            'laiinitialevetr',
            'laiinitialgrass',
        ],
        # 'numcapita':,
        'porosity_id': 'porosity0',
        # 'qn1_av': 'qn1_av',
        # 'qn1_s_av': 'qn1_s_av',
        'snowalb': 'snowalb0',
        'snowwater': [
            'snowwaterbldgsstate',
            'snowwaterpavedstate',
            'snowwaterdectrstate',
            'snowwaterevetrstate',
            'snowwatergrassstate',
            'snowwaterbsoilstate',
            'snowwaterwaterstate',
        ],
        'snowdens': [
            'snowdensbldgs',
            'snowdenspaved',
            'snowdensdectr',
            'snowdensevetr',
            'snowdensgrass',
            'snowdensbsoil',
            'snowdenswater',
        ],
        # 'snowfallcum':,
        'snowfrac': [
            'snowfracbldgs',
            'snowfracpaved',
            'snowfracdectr',
            'snowfracevetr',
            'snowfracgrass',
            'snowfracbsoil',
            'snowfracwater',
        ],
        'snowpack': [
            'snowpackbldgs',
            'snowpackpaved',
            'snowpackdectr',
            'snowpackevetr',
            'snowpackgrass',
            'snowpackbsoil',
            'snowpackwater',
        ],
        'soilstore_id': [
            'soilstorebldgsstate',
            'soilstorepavedstate',
            'soilstoredectrstate',
            'soilstoreevetrstate',
            'soilstoregrassstate',
            'soilstorebsoilstate',
            # 'soilstorewaterstate',
        ],
        'state_id': [
            'bldgsstate',
            'pavedstate',
            'dectrstate',
            'evetrstate',
            'grassstate',
            'bsoilstate',
            'waterstate',
        ],
        # 'tair24hr':,
        # 'tstep_prev':,
        # 'wuday_id':,
    }
    dict_var_initcond_desc = {}

    for var in dict_initcond_related:
        var_related = dict_initcond_related[var]
        desc = df_var_nml_x.loc[var_related, 'Description']
        desc = desc if isinstance(desc, str) else '\n'.join(desc.tolist())
        # print(var_related)
        var_related = df_var_nml_x.loc[var_related, 'SUEWS-related variables']
        var_related = var_related if isinstance(
            var_related, str) else ', '.join(var_related.tolist())
        # print(var_related)
        dict_var = {
            var: pd.Series({
                'Description': desc,
                'SUEWS-related variables': var_related,
            }),
        }

        dict_var_initcond_desc.update(dict_var)
    df_initcond_desc = pd.concat(dict_var_initcond_desc).unstack()

    # retrieve dimension info
    df_var_dim = gen_df_dim(df_init_sample)
    df_initcond_dim = df_var_dim.loc[df_initcond_desc.index]

    # generate init info
    df_var_initcond = pd.concat([df_initcond_dim, df_initcond_desc], axis=1)
    return df_var_initcond


# %%
# split multiple opts into a list of them
# then generate the links
# def gen_rst_url_opt(opt_str):
#     url_str = str(gen_url_option(opt_str.lower()))
#     rst_str = f'`{opt_str} <{url_str}>`_'
#     return rst_str

def gen_rst_url_split_opts(opts_str):
    """generate option list for RST docs

    Parameters
    ----------
    opts_str : str
        a string including all SUEWS related options/variables.
        e.g. 'SUEWS_a, SUEWS_b'


    Returns
    -------
    list
        a list of parsed RST `:ref:` roles.
        e.g. [':option:`SUEWS_a <suews:SUEWS_a>`']
    """
    if opts_str is not 'None':
        list_opts = opts_str.split(',')
        # list_rst = [gen_rst_url_opt(opt.strip()) for opt in list_opts]
        list_rst = [opt.strip() for opt in list_opts]
        # list_rst = [f'`{opt}`' for opt in list_rst]
        # more properly handle SUEWS options by explicitly adding prefix `suews`:
        list_rst = [f':option:`{opt} <suews:{opt}>`' for opt in list_rst]
        list_url_rst = ', '.join(list_rst)
    else:
        list_url_rst = 'None'
    return list_url_rst


def proc_df_state(
        df_var_site: pd.DataFrame,
        df_var_runcontrol: pd.DataFrame,
        df_var_initcond: pd.DataFrame)->pd.DataFrame:

    # dataframe of descriptions
    df_var_desc = pd.concat([df_var_site, df_var_runcontrol, df_var_initcond])
    # split variables of different dimensions
    df_1Dplus = df_var_desc.loc[df_var_desc.Dimensionality != 0]
    pos_1D = df_1Dplus.Dimensionality.map(len) == 1
    df_1D = df_1Dplus.loc[pos_1D]
    pos_23D = df_1Dplus.Dimensionality.map(len) > 1
    df_23D = df_1Dplus.loc[pos_23D]

    # pre-defined dim annotations
    indent = r'    '  # 4 spaces
    sep = ('\n\n'+indent*2)
    dict_dim_anno = {
        0:
        'Scalar',
        (2, ):
        '2: {Weekday, Weekend}',
        (3, ):
        '3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}',
        (7, ):
        '7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}',
        (8, ):
        '8: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`, one extra land cover type (currently NOT used)} ',
        (24, 2):
        sep.join([
            '24: hours of a day',
            '2: {Weekday, Weekend}',
        ]),

        (8, 4, 3):
        sep.join([
            '8: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`, one extra land cover type (currently NOT used)}',
            '4: {SummerWet, SummerDry, WinterWet, WinterDry}',
            '3: {a1, a2, a3}',
        ]),
        (4, 3):
        sep.join([
            '4: {`LeafGrowthPower1`, `LeafGrowthPower2`, `LeafOffPower1`, `LeafOffPower2`}',
            '3: { :term:`EveTr`, :term:`DecTr`, :term:`Grass`}',
        ]),
        (8, 6):
        sep.join([
            '8: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`, one extra land cover type (currently NOT used)}',
            '6: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`}',
        ]),
        (6, 7):
        sep.join([
            '6: { `StorageMin`, `DrainageEq`, `DrainageCoef1`, `DrainageCoef2`, `StorageMax`, current storage}',
            '7: { :term:`Paved`, :term:`Bldgs`, :term:`EveTr`, :term:`DecTr`, :term:`Grass`, :term:`BSoil`, :term:`Water`}',
        ]),
    }

    df_var_dim_anno = df_var_desc['Dimensionality'].map(
        dict_dim_anno).rename('Dimensionality Remarks').to_frame()

    # correct annotation for seven-day dimension
    anno_week = '7: {Sunday, Monday, Tuesday, Wednesday, Thursday, Friday, Saturday}'
    pos_daywat = df_var_dim_anno.index.str.contains('daywat')
    df_var_dim_anno.loc[pos_daywat, 'Dimensionality Remarks'] = anno_week

    # combine description and dimension annotations
    df_var_state_x = pd.concat([df_var_desc, df_var_dim_anno], axis=1)
    df_var_state_x = df_var_state_x.copy().fillna('None')
    df_var_state_x['SUEWS-related variables'] = df_var_state_x['SUEWS-related variables']\
        .map(gen_rst_url_split_opts)
    df_var_state_x.loc[:, 'Description'] = df_var_state_x.loc[:, 'Description'].str.split('\n').map(
        lambda x: pd.Series(x).drop_duplicates().values).str.join(';;')

    # set default values
    df_var_state_x = df_var_state_x.replace(
        {'Description': {'None': 'Internal use. Please DO NOT modify'}})
    df_var_state_x.index = df_var_state_x.index.rename('variable')

    # filter manual-work-demanding entries
    pos_entry_manual = df_var_state_x.loc[:, 'Description'].str.contains(';;')
    df_var_supy_manual = df_var_state_x.copy()[pos_entry_manual]

    # dict for manual edits of several entries
    dict_var_desc_manual = {
        # surface fractions
        'sfr': 'Surface cover fractions.',
        # ohm coefficients:
        'ohm_coef': 'Coefficients for OHM calculation.',
        # profiles:
        'ahprof_24hr': 'Hourly profile values used in energy use calculation.',
        'traffprof_24hr': 'Hourly profile values used in traffic activity calculation.',
        'popprof_24hr': 'Hourly profile values used in dynamic population estimation.',
        'humactivity_24hr': 'Hourly profile values used in human activity calculation.',
        'wuprofa_24hr': 'Hourly profile values used in automatic irrigation.',
        'wuprofm_24hr': 'Hourly profile values used in manual irrigation.',
        'snowprof_24hr': 'Hourly profile values used in snow clearing.',
        # QF related:
        'ah_min': 'Minimum QF values.',
        'ah_slope_heating': 'Heating slope of QF calculation.',
        'ah_slope_cooling': 'Cooling slope of QF calculation.',
        'qf_a': 'Base value for QF calculation.',
        'qf_b': 'Parameter related to heating degree days.',
        'qf_c': 'Parameter related to heating degree days.',
        't_critic_heating': 'Critical heating temperature.',
        't_critic_cooling': 'Critical cooling temperature.',
        'trafficrate': 'Traffic rate used for CO2 flux calculation.',
        # irrigation related:
        'daywat': 'Irrigation flag: 1 for on and 0 for off.',
        'daywatper': 'Fraction of properties using irrigation for each day of a week.',
        'ie_a': 'Coefficient for automatic irrigation model.',
        'ie_m': 'Coefficient for manual irrigation model.',
        # water redistribution related:
        'storedrainprm': 'Coefficients used in drainage calculation.',
        'waterdist': 'Fraction of water redistribution',
        # snow related:
        'snowdens': 'Initial snow density of each land cover.',
        'snowfrac': 'Initial plan area fraction of snow on each land cover`',
        'snowpack': 'Initial snow water equivalent on each land cover',
        'snowwater': 'Initial amount of liquid water in the snow on each land cover',
        # LAI related:
        'laipower': 'parameters required by LAI calculation.',
        'lai_id': 'Initial LAI values.',
        # wetness related:
        'soilstore_id': 'Initial water stored in soil beneath each land cover',
        'state_id': 'Initial wetness condition on each land cover',
    }
    ser_desc_manual = pd.Series(dict_var_desc_manual).rename('Description')

    # update manual entries
    df_var_state = df_var_state_x.copy()
    df_var_state['Description'].update(ser_desc_manual)

    return df_var_state


# %%
# generate dataframe of all state variables used by supy
def gen_df_state(
        list_table: list,
        set_initcond: set,
        set_runcontrol: set,
        set_input_runcontrol: set)->pd.DataFrame:
    '''generate dataframe of all state variables used by supy

    Parameters
    ----------
    list_table : list
        csv files for site info: `SUEWS_xx.csv` on github SUEWS-docs repo
    set_initcond : set
        initial condition related variables
    set_runcontrol : set
        runcontrol related variables
    set_input_runcontrol : set
        runcontrol related variables used as supy input

    Returns
    -------
    pd.DataFrame
        Description of all state variables used by supy
    '''

    # generate a base df for site characteristics related variables
    df_var_site = gen_df_site(list_table)
    # generate a base df for runcontrol related variables
    df_var_runcontrol = gen_df_runcontrol(
        set_initcond, set_runcontrol, set_input_runcontrol)
    # generate a base df for initial condition related variables
    df_var_initcond = gen_df_initcond(set_initcond, set_runcontrol)
    # further processing by modifying several entries
    df_var_state = proc_df_state(
        df_var_site, df_var_runcontrol, df_var_initcond)

    # reorganising the result:
    df_var_state = df_var_state.sort_index()
    # delete duplicates while considering the variable name (stored as index)
    df_var_state = df_var_state.reset_index()
    df_var_state = df_var_state.drop_duplicates()
    # convert index back
    df_var_state = df_var_state.set_index('variable')
    return df_var_state


df_var_state = gen_df_state(
    list_table, set_initcond, set_runcontrol, set_input_runcontrol)
df_var_state.to_csv('df_state.csv')
