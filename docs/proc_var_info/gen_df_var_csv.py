
#%%
from nml_rst_proc import gen_url_option, set_input, set_site, set_initcond, set_runcontrol, set_input_initcond, set_input_runcontrol, parse_option_rst
from urlpath import URL
import urllib.request
import numpy as np
import pandas as pd
import supy as sp
import os
os.getcwd()

#%% [markdown]
# ## generate site characteristics related dataframe
#%% [markdown]
# ### load `SUEWS_***.txt` related tables

#%%
list_table = [file.replace('.txt', '.csv')
              for file in sp.supy_load.list_file_input
              if 'Profile' not in file]  # exclude `Profile` txt: a special case

url_repo_base = 'https://github.com/Urban-Meteorology-Reading/SUEWS-Docs/raw/master/docs/source'
url_table_base = (
    url_repo_base
    + '/input_files/SUEWS_SiteInfo/csv-table/'
)

# url_table_base = (
#     '/Users/sunt05/Dropbox/8-Research/98.ReadingWork/SUWES-Docs/'
#     + 'docs/source/input_files/SUEWS_SiteInfo/csv-table/')

list_url_table = [url_table_base + table for table in list_table]

#%% [markdown]
# #### retrieve SUEWS variable descriptions

#%%
df_var_info = pd.concat([pd.read_csv(f) for f in list_url_table])
df_var_info_x = df_var_info.drop(['No.', 'Use'], axis=1)
df_var_info_x = df_var_info_x.set_index('Column Name')
df_var_info_x.index = df_var_info_x.index.map(lambda x: x.replace('`', ''))
df_var_info_x.head()

#%% [markdown]
# #### retrieve SUEWS-related variable names

#%%
# retrieve SUEWS-related variables
dict_var2SiteSelect = sp.supy_load.dict_var2SiteSelect

dict_var_full = sp.supy_load.exp_dict_full(dict_var2SiteSelect)


def extract_var_suews(dict_var_full, k):
    x = sp.supy_load.flatten_list(dict_var_full[k])
    x = np.unique(x)
    x = [
        xx for xx in x
        if xx not in ['base', 'const', '0.0'] + [str(x) for x in range(24)]
    ]
    x = [xx for xx in x if 'Code' not in xx]
    return x


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

# df_var_ref_suews.head(10)
# set_input_site_exp0=df_var_ref_suews.loc[:,'SUEWS-related variables'].str.lower().str.split(',')

# set_input_site_exp= set(x.lower() for x in np.concatenate(set_input_site_exp0.values))


#%%
# retrive supy variable description
dict_var_desc = {k: '\n'.join(df_var_info_x.loc[v].values.flatten())
                 for k, v in dict_var_ref_suews.items()}

df_var_desc = pd.DataFrame(dict_var_desc, index=[0]).T.rename({
    0: 'Description'}, axis=1)

df_var_desc.head()

#%% [markdown]
# #### list supy variables that are combinations of multiple SUEWS variables

#%%
for k, v in dict_var_ref_suews.items():
    if len(v) > 1:
        print('')
        print(k, v)
    # df_var_info_x.loc[v].values.flatten()

#%% [markdown]
# #### retrieve variable dimensionality

#%%
df_init_sample, df_forcing_sample = sp.load_SampleData()
df_var_dim_x = df_init_sample.columns.to_frame()
df_var_dim_x.index = df_var_dim_x.index.droplevel(-1)

df_var_dim_x = df_var_dim_x.groupby('var').last()

df_var_dim = df_var_dim_x.applymap(eval).applymap(
    lambda x: tuple(np.array(x) + 1) if isinstance(x, tuple) else x)

df_var_dim = df_var_dim.rename({'ind_dim': 'Dimensionality'}, axis=1)
df_var_dim = pd.DataFrame(df_var_dim)
df_var_dim.filter(like='method', axis=0)
df_var_dim.head()

#%% [markdown]
# ### select only those useful for input to `supy`

#%%

df_var_site_raw = pd.concat(
    [df_var_dim, df_var_desc, df_var_ref_suews], axis=1)
df_var_site = df_var_site_raw.filter(items=set_input, axis=0).dropna()
# df_var_site.to_csv('df_var_site.csv')


#%%
# input options other than those from site characteristics
set_input_site = set(df_var_site.index)
set_input_nonsite = set_input.difference(set_input_site)
set_input_nonsite

#%% [markdown]
# ## process `runcontrol` and `initialcondition` related variables
#%% [markdown]
# ### these variables are based on `nml` files so the processing logic is different

#%%
list_var_nml = list(set_initcond)+list(set_runcontrol)
list_url_rst_x = list(gen_url_option(var, source='github')
                      for var in sorted(list_var_nml))

# gen_url_option('tstep')
sorted(list_var_nml)


#%%
# URL('/'.join(list_var_url[0].parts)).get().ok
list_url_rst = list(set(URL('/'.join(url.parts)) for url in list_url_rst_x))


#%%

df_var_nml = pd.concat(parse_option_rst(url) for url in list_url_rst)
df_var_nml_x = df_var_nml.reset_index().loc[:, ['Description', 'index']].rename(
    columns={'index': 'SUEWS-related variables'})
index_nml = df_var_nml_x.loc[:, 'SUEWS-related variables'].str.lower()
df_var_nml_x.index = index_nml.rename('var')
df_var_nml_x.sort_index()

#%% [markdown]
# ### generate `df_var_runcontrol` for supy

#%%
df_var_runcontrol_x = df_var_nml_x.filter(items=set_input_runcontrol, axis=0)
df_var_runcontrol_dim = df_var_dim.filter(items=set_input_runcontrol, axis=0)
df_var_runcontrol = pd.concat(
    [df_var_runcontrol_dim, df_var_runcontrol_x], axis=1)
df_var_runcontrol

#%% [markdown]
# ### generate `df_var_initcond` for supy
#%% [markdown]
# #### define related variables for `initcond` variables

#%%
set_input_nonsite - set_input_runcontrol


#%%
set_initcond


#%%
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

#%% [markdown]
# #### examine related variables

#%%
for var in dict_initcond_related:
    print(dict_initcond_related[var])
    df_var_nml_x.loc[dict_initcond_related[var], 'SUEWS-related variables']


#%%
df_var_nml_x.loc['albdectr0']
df_var_nml_x.index

#%% [markdown]
# #### generate `df_var_initcond`

#%%
dict_var_initcond_desc = {}
for var in dict_initcond_related:
    var_related = dict_initcond_related[var]
    desc = df_var_nml_x.loc[var_related, 'Description']
    desc = desc if isinstance(desc, str) else '\n'.join(desc.tolist())
    print(var_related)
    var_related = df_var_nml_x.loc[var_related, 'SUEWS-related variables']
    var_related = var_related if isinstance(
        var_related, str) else ', '.join(var_related.tolist())
    print(var_related)
    dict_var = {
        var: pd.Series({
            'Description': desc,
            'SUEWS-related variables': var_related,
        }),
    }

    dict_var_initcond_desc.update(dict_var)
df_initcond_desc = pd.concat(dict_var_initcond_desc).unstack()
df_initcond_dim = df_var_dim.loc[df_initcond_desc.index]
df_var_initcond = pd.concat([df_initcond_dim, df_initcond_desc], axis=1)
df_var_initcond


#%%
df_initcond_desc
# dict_initcond_related['albdectr_id']
df_var_nml_x.sort_index()

#%% [markdown]
# ## summarise the `DataFrame`s
#%% [markdown]
# ### combine all `df_var_xx`s

#%%
df_var_supy_docs = pd.concat([df_var_site, df_var_runcontrol, df_var_initcond])
df_var_supy_docs.head()

#%% [markdown]
# ### these variables are for internal use and thus not released for users

#%%
set_input - set(df_var_supy_docs.index)

#%% [markdown]
# ### deal with dim info

#%%
df_var_supy_docs.loc[df_var_supy_docs.Dimensionality != 0]

#%% [markdown]
# #### 1D+ data

#%%
df_1Dplus = df_var_supy_docs.loc[df_var_supy_docs.Dimensionality != 0]
pos_1D = df_1Dplus.Dimensionality.map(len) == 1
df_1D = df_1Dplus.loc[pos_1D]
pos_23D = df_1Dplus.Dimensionality.map(len) > 1
df_23D = df_1Dplus.loc[pos_23D]


#%%
df_1D


#%%
df_23D


#%%
df_1D.Dimensionality.unique()


#%%
df_23D.Dimensionality.unique()


#%%
dict_dim_anno = {
    0:
    'Scalar',
    (2, ):
    '2: Weekday and Weekend',
    (3, ):
    '3: See variable description for specifics',
    (7, ):
    '7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html',
    (8, ):
    '8: Seven SUEWS land cover types and one extra land cover type (currently NOT used)',
    (24, 2):
    '24: hours of a day; 2: Weekday and Weekend',
    (8, 4, 3):
    '; '.join([
        '8: Seven SUEWS land cover types and one extra land cover type (currently NOT used)',
        '4: SummerWet, SummerDry, WinterWet, WinterDry',
        '3: a1, a2, a3'
    ]),
    (4, 3):
    '4: See variable description for specifics; 3: Three vegetated land cover types (`EveTr`, `DecTr`, `Grass`)',
    (8, 6):
    '8: Seven SUEWS land cover types and Runoff/SoilStore as water receiver; 6: SUEWS land cover types other than water as water contributors',
    (6, 7):
    '6: See variable description for specifics; 7: Seven SUEWS land cover types in the order of [paved; buildings; evergreen trees/shrubs; deciduous trees/shrubs; grass; bare soil and water]: https://suews-docs.readthedocs.io/en/latest/introduction.html',
}


#%%
df_var_dim_anno = df_var_supy_docs['Dimensionality'].map(
    dict_dim_anno).rename('Dimensionality Remarks').to_frame()
# correct the `(7,)` for week days
anno_week = '7: Seven days of a week: from Sunday to Saturday'
df_var_dim_anno.loc[
    df_var_dim_anno.index.str.contains('daywat'), 
    'Dimensionality Remarks'] = anno_week
df_var_dim_anno.sort_index().iloc[:36]

#%% [markdown]
# #### attach dim info to df_docs

#%%
df_var_supy_docs_dim = pd.concat([df_var_supy_docs, df_var_dim_anno], axis=1)
df_var_supy_docs_dim

#%% [markdown]
# ### add docs URL to SUEWS related variables

#%%
opt_str = 'SnowAlb0'
url_str = str(gen_url_option(opt_str.lower()))
f'`{opt_str} <{url_str}>`_'


#%%
def gen_rst_url_opt(opt_str):
    url_str = str(gen_url_option(opt_str.lower()))
    rst_str = f'`{opt_str} <{url_str}>`_'
    return rst_str


gen_rst_url_opt('QF0_BEU_WD')


#%%
'QF0_BEU_WD, QF0_BEU_WE'.split(',')


#%%
# split multiple opts into a list of them
# then generate the links
def gen_rst_url_split_opts(opts_str):
    if opts_str is not 'None':
        list_opts = opts_str.split(',')
        list_rst = [gen_rst_url_opt(opt.strip()) for opt in list_opts]
        list_url_rst = ', '.join(list_rst)
    else:
        list_url_rst = 'None'
    return list_url_rst


gen_rst_url_split_opts('None')


#%%
xx = df_var_supy_docs_dim.loc[:, 'SUEWS-related variables'].sort_index()
xx.loc[xx.isna()]


#%%
df_var_supy_docs_dim_url = df_var_supy_docs_dim.copy().fillna('None')
df_var_supy_docs_dim_url['SUEWS-related variables'] = df_var_supy_docs_dim_url['SUEWS-related variables'].map(
    gen_rst_url_split_opts)
# df_var_supy_docs_dim_url=df_var_supy_docs_dim.transform({'SUEWS-related variables':gen_rst_url_split_opts})
df_var_supy_docs_dim_url

#%% [markdown]
# ### clean description

#%%
# drop duplicates
xx = df_var_supy_docs_dim_url.loc[:, 'Description']    .str.split('\n')    .map(lambda x: pd.Series(x).drop_duplicates().values)    .str.join(';;')
df_var_supy_docs_dim_url.loc[:, 'Description'] = xx
df_var_supy_docs_dim_url.head()


#%%
# set default values
df_var_supy_docs_dim_url = df_var_supy_docs_dim_url.replace(
    {'Description': {'None': 'Internal use. Please DO NOT modify'}})
df_var_supy_docs_dim_url.index = df_var_supy_docs_dim_url.index.rename('variable')

#%% [markdown]
# ### save to csv files for rst generation

#%%
# filter manual-work-demanding entries
pos_entry_manual = df_var_supy_docs_dim_url.loc[:, 'Description'].str.contains(';;')
df_var_supy_manual = df_var_supy_docs_dim_url.copy()[pos_entry_manual]
df_var_supy_manual.to_csv('df_var_supy_manual.csv')

# save auto-generation-safe entries
df_var_supy_auto = df_var_supy_docs_dim_url.copy()[~pos_entry_manual]
df_var_supy_auto.to_csv('df_var_supy_auto.csv')


