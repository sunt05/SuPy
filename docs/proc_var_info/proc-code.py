# this helper code is used to retrieve basic info from SUEWS-docs
# for variable meaning of `df_init`

import numpy as np
import pandas as pd

import supy as sp

list_table = [file.replace('.txt', '.csv')
              for file in sp.supy_load.list_file_input
              if 'Profile' not in file]

# url_table_base = (
#     'https://github.com/Urban-Meteorology-Reading/SUEWS-Docs' +
#     '/raw/master/docs/source/input_files/SUEWS_SiteInfo/csv-table/')

url_table_base = (
    '/Users/sunt05/Dropbox/8-Research/98.ReadingWork/SUWES-Docs/' +
    'docs/source/input_files/SUEWS_SiteInfo/csv-table/')

list_url_table = [url_table_base + table for table in list_table]

# generate basic var info
df_var_info = pd.concat([pd.read_csv(f) for f in list_url_table])
df_var_info_x = df_var_info.drop(['No.', 'Use'], axis=1)
df_var_info_x = df_var_info_x.set_index('Column Name')
# df_var_info_x=df_var_info_x.drop_duplicates()
df_var_info_x.index = df_var_info_x.index.map(lambda x: x.replace('`', ''))
df_var_info_x.head()


# retrieve SUEWS-related variables
dict_var2SiteSelect = sp.supy_load.dict_var2SiteSelect

dict_var_full = sp.supy_load.exp_dict_full(dict_var2SiteSelect)


def extract_var_suews(dict_var_full, k):
    x = sp.supy_load.flatten_list(dict_var_full[k])
    x = np.unique(x)
    x = [xx for xx in x if xx not in [
        'base', 'const', '0.0'] + [str(x) for x in range(24)]]
    x = [xx for xx in x if 'Code' not in xx]
    return x


dict_var_ref_suews = {
    k: extract_var_suews(dict_var_full, k)
    for k in dict_var_full}

df_var_ref_suews = pd.DataFrame(
    {k: ', '.join(dict_var_ref_suews[k])
     for k in dict_var_ref_suews},
    index=[0]).T.rename(
    {0: 'SUEWS-related variables'}, axis=1)

df_var_ref_suews.head()


# retrive variable description
dict_var_desc = {k: '\n'.join(df_var_info_x.loc[v].values.flatten())
                 for k, v in dict_var_ref_suews.items()}

df_var_desc = pd.DataFrame(dict_var_desc, index=[0]).T.rename({
    0: 'Description'}, axis=1)

df_var_desc.head()


# retrieve variable dimensionality
_, res_sample, _ = sp.load_SampleData()
df_var_dim_x = res_sample.columns.to_frame()
df_var_dim_x.index = df_var_dim_x.index.droplevel(-1)

df_var_dim_x = df_var_dim_x.groupby('var').last()

df_var_dim = df_var_dim_x.applymap(eval).applymap(
    lambda x: tuple(np.array(x) + 1) if isinstance(x, tuple) else x)

df_var_dim = df_var_dim.rename({'ind_dim': 'Dimensionality'}, axis=1)


# save raw info as csv
# select only those useful
set_input = sp.supy_run.set_var_input
set_input.update(sp.supy_run.set_var_input_multitsteps)
df_var_csv_raw = pd.concat([df_var_dim, df_var_desc, df_var_ref_suews], axis=1)
df_var_csv_raw = df_var_csv_raw.filter(items=set_input, axis=0)
df_var_csv_raw.to_csv('df_var_csv_raw.csv')
