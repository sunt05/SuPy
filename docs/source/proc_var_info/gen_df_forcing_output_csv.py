# %% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataSciece.changeDirOnImportExport setting

import os
try:
    os.chdir(os.path.join(os.getcwd(), 'docs/proc_var_info'))
    print(os.getcwd())
except:
    pass

# %%
from urlpath import URL
from pathlib import Path
import numpy as np
import pandas as pd
import supy as sp
import os
os.getcwd()
# %% sample run
print('loading in', 'gen_df_forcing', '...')
df_state_init_sample, df_forcing_sample = sp.load_SampleData()
df_output_sample, df_state_end_sample = sp.run_supy(
    df_forcing_sample.iloc[:10], df_state_init_sample)
print('loading in', 'gen_df_forcing', '...')

# %% [markdown]
# ## generate forcing related dataframe
# %% [markdown]
# ### load `SUEWS_***.txt` related tables

# %%
url_repo_base = ('https://github.com/'
                 + 'Urban-Meteorology-Reading/'
                 + 'SUEWS-Docs/raw/master/docs/source')
url_repo_input = URL(url_repo_base)/'input_files'
url_repo_output = URL(url_repo_base)/'output_files'


def gen_df_forcing(
        path_csv_in='SSss_YYYY_data_tt.csv',
        url_base=url_repo_input,)->pd.DataFrame:
    '''Generate description info of supy forcing data into a dataframe

    Parameters
    ----------
    path_csv_in : str, optional
        path to the input csv file relative to url_base (the default is '/input_files/SSss_YYYY_data_tt.csv'])
    url_base : urlpath.URL, optional
        URL to the input files of repo base (the default is url_repo_input, which is defined at the top of this file)

    Returns
    -------
    pd.DataFrame
        Description info of supy forcing data
    '''

    try:
        # load info from SUEWS docs repo
        # this is regarded as the official source
        urlpath_table = url_base/path_csv_in
        df_var_info = pd.read_csv(urlpath_table)
    except:
        print(f'{urlpath_table} not existing!')
    else:
        # clean info dataframe
        df_var_forcing = df_var_info.drop(['No.', 'Use'], axis=1)

        # set index with `Column name`
        df_var_forcing = df_var_forcing.set_index('Column Name')
        df_var_forcing.index = df_var_forcing.index\
            .map(lambda x: x.replace('`', ''))\
            .rename('variable')

        # add `Second` info
        df_var_forcing.loc['isec'] = 'Second [S]'

        return df_var_forcing


# %% [markdown]
# ## generate output related dataframe


# %%
def gen_df_output(
        list_csv_in=[
            'SSss_YYYY_SUEWS_TT.csv',
            'SSss_DailyState.csv',
            'SSss_YYYY_snow_TT.csv',
        ],
        url_base=url_repo_output)->Path:
    '''Generate description info of supy output results into dataframe

    Parameters
    ----------
    list_csv_in : list, optional
        list of file names for csv files with meta info (the default is ['SSss_YYYY_SUEWS_TT.csv','SSss_DailyState.csv','SSss_YYYY_snow_TT.csv',], which [default_description])
    url_base : [type], optional
        URL to the output dir of repo base (the default is url_repo_output, which is defined at the top of this file)

    Returns
    -------
    pd.DataFrame
         Description info of supy output results
    '''

    # list of URLs
    list_url_table = [
        url_base/table for table in list_csv_in
    ]
    try:
        df_var_info = pd.concat(
            [pd.read_csv(f) for f in list_url_table],
            sort=False)
    except:
        for url in list_url_table:
            if not url.get().ok:
                print(f'{url} not existing!')
    else:
        # clean meta info
        df_var_info_x = df_var_info\
            .set_index('Name')\
            .loc[:, ['Description']]\
            .drop_duplicates()

        df_var_output = df_var_info_x\
            .copy()\
            .assign(lower=df_var_info_x.index.str.lower())\
            .reset_index()\
            .set_index('lower')

        df_var_group = df_output_sample.columns.to_frame()
        df_var_group.index = df_var_group.index.droplevel(0).rename('Name')

        # wrap into a dataframe
        df_var_output = df_var_group\
            .merge(
                df_var_output.set_index('Name'),
                left_on='Name',
                right_on='Name')\
            .rename(columns={
                'var': 'variable',
                'group': 'Group',
            })\
            .set_index('variable')\
            .drop_duplicates()

        return df_var_output


# %% [markdown]
# ## generate csv files for meta info
# %%
# df_forcing=gen_df_forcing('SSss_YYYY_data_tt.csv')

# df_output=gen_df_output(
#     [
#         'SSss_YYYY_SUEWS_TT.csv',
#         'SSss_DailyState.csv',
#         'SSss_YYYY_snow_TT.csv',
#     ],
# )


# # %%
# df_forcing.head()


# #%%
# df_output.head()


# #%%
