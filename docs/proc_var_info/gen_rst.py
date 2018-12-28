# generate option rst files
# %%
import pandas as pd

df_var_info = pd.concat([
    pd.read_csv(file)
    for file in [
        'df_var_supy_auto.csv',
        'df_var_supy_manual_mod.csv',
    ]
])
df_var_info


# %%
# xx = df_var_info.iloc[100]
#
# xx
# ser_rec = xx
# name = ser_rec.variable
# dim = ser_rec.Dimensionality
# vars = ser_rec['SUEWS-related variables']
# desc = ser_rec.Description

# %%


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

# %%

df_var_info['rst'] = df_var_info.apply(gen_opt_str, axis=1)
df_var_info = df_var_info.sort_values('var').reset_index(drop=True)
rst_txt = '\n\n'.join(df_var_info.rst)
# print()
with open('df_state_raw.rst', 'w') as f:
    print(rst_txt, file=f)

df_var_info.tail()


#%%
