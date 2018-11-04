# generate option rst files

import pandas as pd

df_var_info = pd.read_csv('df_var_csv.csv').dropna(how='all')


# xx = df_var_info.iloc[100]
#
# xx
# ser_rec = xx
# name = ser_rec.variable
# dim = ser_rec.Dimensionality
# vars = ser_rec['SUEWS-related variables']
# desc = ser_rec.Description


def gen_opt_str(ser_rec):
    name = ser_rec.variable
    dim = ser_rec.Dimensionality
    vars = ser_rec['SUEWS-related variables']
    desc = ser_rec.Description
    str_opt = '''
.. option:: {name}

	:Dimensionality:
		{dim}
	:Description:
		{desc}
	:SUEWS-related variables:
		{vars}
    '''.format(
        name=name,
        dim=dim,
        desc=desc,
        vars=vars,
    )
    return str_opt


df_var_info['rst'] = df_var_info.apply(gen_opt_str, axis=1)
df_var_info=df_var_info.sort_values('variable').reset_index(drop=True)
rst_txt = '\n\n'.join(df_var_info.rst)
# print()
with open('df_state_raw.rst', 'w') as f:
    print(rst_txt, file=f)

df_var_info.tail()
