#%% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataSciece.changeDirOnImportExport setting
import os
try:
	os.chdir(os.path.join(os.getcwd(), 'docs/proc_var_info'))
	print(os.getcwd())
except:
	pass
#%% [markdown]
# # generate option rst files

#%%
import pandas as pd

#%% [markdown]
# ## load processed csv files of df_state info
# %%
df_var_info = pd.read_csv('./df_state.csv')

#%% [markdown]
# ## generate option string for rst option file
#%%
# generate option string for rst option file
def gen_opt_str(ser_rec):

    name = ser_rec['variable']
    desc = ser_rec['Description']
    vars = ser_rec['SUEWS-related variables']
    dim = ser_rec['Dimensionality']
    remark = ser_rec['Dimensionality Remarks']
    str_opt = '''
.. option:: {name}

    :Description:
        {desc}
    :SUEWS-related variables:
        {vars}
    :Dimensionality:
        {dim}
    :Dimensionality remarks:
        {remark}'''.format(
        name=name,
        desc=desc,
        vars=vars,
        dim=dim,
        remark=remark,
    )
    return str_opt


gen_opt_str(df_var_info.iloc[10])


#%%

df_var_info['rst'] = df_var_info.apply(gen_opt_str, axis=1)
df_var_info = df_var_info.sort_values('variable').reset_index(drop=True)

rst_txt_x = '\n\n'.join(df_var_info.rst)

rst_title='''
.. _df_state_options:

Model state variables
============================

'''

rst_txt = '\n'.join([rst_title, rst_txt_x])

# print()
with open('../source/data-structure/df_state.rst', 'w') as f:
    print(rst_txt, file=f)

df_var_info.tail()




#%%
