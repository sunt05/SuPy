# %% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataSciece.changeDirOnImportExport setting
import os
try:
    os.chdir(os.path.join(os.getcwd(), 'docs/proc_var_info'))
    print(os.getcwd())
except:
    pass

# %%
from pathlib import Path
import pandas as pd
import supy as sp
import os
try:
    os.chdir(os.path.join(os.getcwd(), 'docs/proc_var_info'))
    print(os.getcwd())
except:
    pass
# get_ipython().run_line_magic('load_ext', 'autoreload')
# get_ipython().run_line_magic('autoreload', '2')


# %%
from gen_df_state_csv import (gen_df_state, list_table, set_initcond,
                              set_runcontrol, set_input_runcontrol, gen_df_dim)
from gen_df_forcing_output_csv import gen_df_forcing, gen_df_output

# %% [markdown]
# # generate option rst files
# %% [markdown]
# ## generate dataframes for variable groups

# %%
print('generating df_state.csv ...')
df_state = gen_df_state(
    list_table, set_initcond, set_runcontrol, set_input_runcontrol)
df_state.to_csv('df_state.csv')
print('df_state.csv done!')


# #%%
# get_ipython().run_line_magic('load_ext', 'snakeviz')
# get_ipython().run_line_magic('snakeviz', 'gen_df_state(list_table, set_initcond, set_runcontrol, set_input_runcontrol)')


# %%
print('generating df_forcing.csv ...')
df_forcing = gen_df_forcing('SSss_YYYY_data_tt.csv')
df_forcing.to_csv('df_forcing.csv')
print('df_forcing.csv done!')


# %%
print('generating df_output.csv ...')
df_output = gen_df_output(
    [
        'SSss_YYYY_SUEWS_TT.csv',
        'SSss_DailyState.csv',
        'SSss_YYYY_snow_TT.csv',
    ],
)
df_output.to_csv('df_output.csv')
print('df_output.csv done!')

# %% [markdown]
# ## generate option string for rst option file

# %%


def gen_opt_str(ser_rec: pd.Series)->str:
    '''generate rst option string

    Parameters
    ----------
    ser_rec : pd.Series
        record for specifications

    Returns
    -------
    str
        rst string
    '''

    name = ser_rec.name
    indent = r'    '
    str_opt = f'.. option:: {name}'+'\n\n'
    for spec in ser_rec.sort_index().index:
        str_opt += indent+f':{spec}:'+'\n'
        spec_content = ser_rec[spec]
        str_opt += indent+indent+f'{spec_content}'+'\n'
    return str_opt


# xx=df_var_info.set_index('variable').iloc[10]
# print(gen_opt_str(xx))


# %%
def gen_rst(path_rst, path_df_csv, rst_title):
    df_var_info = pd.read_csv(path_df_csv).set_index('variable')
    df_var_info['rst'] = df_var_info.copy().apply(gen_opt_str, axis=1)
    df_var_info = df_var_info.sort_index().reset_index(drop=True)
    rst_txt_x = '\n\n'.join(df_var_info.rst)
    rst_txt = '\n'.join([rst_title, rst_txt_x])
    with open(path_rst, 'w') as f:
        print(rst_txt, file=f)

    return path_rst


# gen_rst(
#     '../source/data-structure/test.rst',
#     df_state,
#     'xx\n')


# %%
def gen_group_dict(
    group,
    path_rst_base=Path('../data-structure/')
)->dict:

    rst_title = f'''
.. _df_{group}_var:

``df_{group}`` variables
============================

'''

    dict_group = {
        'path_rst': path_rst_base/('df_'+group+'.rst'),
        'path_df_csv': 'df_'+group+'.csv',
        'rst_title': rst_title,
    }

    return dict_group


# print(gen_group_dict('state'))


# %%

dict_rst_out = {group: gen_group_dict(group)
                for group in ['state', 'forcing', 'output']}
# dict_rst_out


# %%
for group in dict_rst_out:
    print('working on group:', group)
    print('file generated:', gen_rst(**dict_rst_out[group]), '\n')
