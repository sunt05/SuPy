# %%
import webbrowser
import os
import supy as sp
import numpy as np
from urlpath import URL
import pandas as pd

os.getcwd()
# %% [markdown]
# ### filter input variables
# %%
set_input = sp.supy_load.set_var_input.copy()
set_input.update(sp.supy_load.set_var_input_multitsteps)
df_init_sample, df_forcing_sample = sp.load_SampleData()
set_input.difference_update(set(df_forcing_sample.columns))
# set_input, len(set_input)

# %% [markdown]
# #### retrieve SUEWS-related variables

# %%
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

ser_input_site_exp = df_var_ref_suews.filter(
    items=set_input, axis=0).loc[:, 'SUEWS-related variables'].str.lower().str.split(',')

set_site = set(
    x.lower().strip()
    for x in np.concatenate(ser_input_site_exp.values))

# set_site, len(set_site)
# %% [markdown]
# ### filter `runcontrol` related variables
# %%
# runcontrol variables for supy input
path_runcontrol = sp.supy_env.path_supy_module / 'sample_run'/'Runcontrol.nml'
dict_runcontrol = sp.supy_load.load_SUEWS_dict_ModConfig(
    path_runcontrol).copy()
set_runcontrol = set(dict_runcontrol.keys())
set_input_runcontrol = set_runcontrol.intersection(set_input)
set_input_runcontrol, len(set_input_runcontrol)

# %% [markdown]
# ### filter `initialcondition` related variables
# %%
# initcond variables for supy input
dict_initcond = sp.supy_load.dict_InitCond_default.copy()
set_initcond = set(dict_initcond.keys())
set_input_initcond = set_initcond.intersection(set_input)
# set_input_initcond, len(set_input_initcond)


# %% [markdown]
# ### functions to process `nml` related variables
# %%
def form_option(str_opt):
    '''generate option name based suffix for URL

    :param str_opt: opt name
    :type str_opt: str
    :return: URL suffix for the specified option
    :rtype: str
    '''

    str_base = '#cmdoption-arg-'
    str_opt_x = str_base+str_opt.lower()\
        .replace('_', '-')\
        .replace('(', '-')\
        .replace(')', '')
    return str_opt_x


# form_option('snowinitially')


# %%
def choose_page(str_opt, set_site=set_site, set_runcontrol=set_runcontrol, set_initcond=set_initcond, source='docs'):
    # print('str_opt', str_opt)
    suffix_page = ('html' if source is 'docs' else 'rst')
    # runcontrol variable:
    if str_opt in set_runcontrol:
        str_base = 'RunControl'
        if str_opt.startswith('tstep'):
            name_page = 'Time_related_options'
        else:
            name_page = 'scheme_options'

    # initcondition variable:
    elif str_opt in set_initcond:
        str_base = 'Initial_Conditions'
        # the following test sequence is IMPORTANT!
        if str_opt.startswith('soilstore'):
            name_page = 'Soil_moisture_states'
        elif str_opt.startswith('snow'):
            name_page = 'Snow_related_parameters'
        elif str_opt.endswith('state'):
            name_page = 'Above_ground_state'
        elif str_opt in ('dayssincerain', 'temp_c0'):
            name_page = 'Recent_meteorology'
        else:
            name_page = 'Vegetation_parameters'

    # site characteristics variable:
    elif str_opt in set_site:
        str_base = 'SUEWS_SiteInfo'
        name_page = 'Input_Options'

    # defaults to empty strings
    else:
        str_base = ''
        name_page = ''

    str_page = '.'.join([name_page, suffix_page])
    str_page_full = str_base+'/'+str_page
    return str_page_full

# for source in ['docs','github']:
#     print(source)
#     for x in sorted(list(set_site)+list(set_runcontrol)+list(set_initcond)):
#         print()
#         print(choose_page(x, source=source))
# choose_page('tstep', set_site, set_runcontrol, set_initcond)
# choose_page('snowinitially', set_site, set_runcontrol, set_initcond)


# %%
def gen_url_option(
    str_opt, 
    set_site=set_site, 
    set_runcontrol=set_runcontrol, 
    set_initcond=set_initcond, 
    source='docs'):
    '''construct a URL for option based on source 

    :param str_opt: option name, defaults to ''
    :param str_opt: str, optional
    :param source: URL source: 'docs' for readthedocs.org; 'github' for github repo, defaults to 'docs'
    :param source: str, optional
    :return: a valid URL pointing to the option related resources
    :rtype: urlpath.URL
    '''
    dict_base = {
        'docs': URL('https://suews-docs.readthedocs.io/en/latest/input_files/'),
        'github': URL('https://github.com/Urban-Meteorology-Reading/SUEWS-Docs/raw/master/docs/source/input_files/'),
    }
    url_base = dict_base[source]

    url_page = choose_page(
        str_opt, set_site, set_runcontrol, set_initcond, source=source)
    # print('str_opt', str_opt, url_base, url_page)
    str_opt_x = form_option(str_opt)
    url_opt = url_base/(url_page+str_opt_x)
    return url_opt


# for source in [
#     # 'docs',
#     'github',
#     ]:
#     print(source)
#     for x in sorted(list(set_site)+list(set_runcontrol)+list(set_initcond)):
#         print()
#         print(gen_url_option(x, source=source))

#         # webbrowser.open(str(gen_url_option(x, source=source)))

# gen_url_option('sss', source='github')
# gen_url_option('sss', source='github').get().ok
# %%
# test connectivity of all generated option URLs
# for opt in list(set_initcond)+list(set_runcontrol):
#     for source in ['github', 'docs']:
#         url = gen_url_option(opt, source=source)
#         if not url.get().ok:
#             print()
#             print(opt)
#             print(url)

# %%
def parse_block(block):
    xx = block.reset_index(drop=True)
    name_block = xx.loc[0].replace('.. option::', '').strip()
    ind_field = xx.index[xx.str.startswith('\t:')]
    list_field = [xx.iloc[slice(*x)].str.strip().reset_index(drop=True)
                  for x in zip(ind_field, list(ind_field[1:])+[None])]
    name_field = [field.loc[0].replace(':', '') for field in list_field]
    content_field = [field.loc[1:].str.join('') for field in list_field]
    ser_field = pd.Series({field.loc[0].replace(':', ''): ' '.join(
        field.loc[1:]) for field in list_field}).rename(name_block)
    return ser_field


def parse_option_rst(path_rst):
    ser_opts = pd.read_csv(
        path_rst,
        sep='\n', skipinitialspace=True)
    ser_opts = ser_opts.iloc[:, 0]
    ind_opt = ser_opts.index[ser_opts.str.contains('.. option::')]
    ser_opt_name = ser_opts[ind_opt].str.replace('.. option::', '').str.strip()
    list_block_opt = [ser_opts.loc[slice(*x)]
                      for x in zip(ind_opt, list(ind_opt[1:])+[None])]
    df_opt = pd.concat([parse_block(block)
                        for block in list_block_opt], axis=1).T
    return df_opt


# url_test = gen_url_option('pavedstate', source='github')
# parse_option_rst(url_test)


# %%
