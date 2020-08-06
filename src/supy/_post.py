import numpy as np
import pandas as pd
from supy_driver import suews_driver as sd


##############################################################################
# post-processing part
# get variable information from Fortran
def get_output_info_df():
    size_var_list = sd.output_size()
    var_list_x = [np.array(sd.output_name_n(i)) for i in np.arange(size_var_list) + 1]

    df_var_list = pd.DataFrame(var_list_x, columns=["var", "group", "aggm", "outlevel"])
    df_var_list = df_var_list.applymap(lambda x: x.decode().strip())
    df_var_list_x = df_var_list.replace(r"^\s*$", np.nan, regex=True).dropna()
    var_dfm = df_var_list_x.set_index(["group", "var"])
    return var_dfm


# get variable info as a DataFrame
# save `var_df` for later use
var_df = get_output_info_df()

# dict as var_df but keys in lowercase
var_df_lower = {group.lower(): group for group in var_df.index.levels[0].str.strip()}

#  generate dict of functions to apply for each variable
dict_func_aggm = {
    "T": "last",
    "A": "mean",
    "S": "sum",
    "L": "last",
}
var_df["func"] = var_df.aggm.apply(lambda x: dict_func_aggm[x])

# dict of resampling ruls:
#  {group: {var: agg_method}}
dict_var_aggm = {
    group: var_df.loc[group, "func"].to_dict() for group in var_df.index.levels[0]
}


# generate index for variables in different model groups
def gen_group_cols(group_x):
    # get correct group name by cleaning and swapping case
    group = group_x.replace("dataoutline", "").replace("line", "")
    # print group
    group = var_df_lower[group]
    header_group = np.apply_along_axis(
        list, 0, var_df.loc[["datetime", group]].index.values
    )[:, 1]

    # generate MultiIndex if not `datetimeline`
    if not group_x == "datetimeline":
        index_group = pd.MultiIndex.from_product(
            [[group], header_group], names=["group", "var"], sortorder=None
        )
    else:
        index_group = header_group

    return index_group


# merge_grid: useful for both `dict_output` and `dict_state`
def pack_df_grid(dict_output):
    # pack all grid and times into index/columns
    df_xx = pd.DataFrame.from_dict(dict_output, orient="index")
    # pack
    df_xx0 = df_xx.applymap(pd.Series)
    df_xx1 = df_xx0.applymap(pd.DataFrame.from_dict)
    df_xx2 = pd.concat(
        {
            grid: pd.concat(df_xx1[grid].to_dict()).unstack().dropna(axis=1)
            for grid in df_xx1.columns
        }
    )
    # drop redundant levels
    df_xx2.columns = df_xx2.columns.droplevel(0)
    # regroup by `grid`
    df_xx2.index.names = ["grid", "time"]
    gb_xx2 = df_xx2.groupby(level="grid")
    # merge results of each grid
    xx3 = gb_xx2.agg(lambda x: tuple(x.values)).applymap(np.array)

    return xx3


# generate MultiIndex for variable groups
def gen_index(varline_x):
    var_x = varline_x.replace("dataoutline", "").replace("line", "")
    group = var_df_lower[var_x]
    var = var_df.loc[group].index.tolist()
    mindex = pd.MultiIndex.from_product([[group], var], names=["group", "var"])
    return mindex


# generate one MultiIndex from a whole dict
def gen_MultiIndex(dict_x):
    x_keys = dict_x.keys()
    mindex = pd.concat([gen_index(k).to_frame() for k in x_keys]).index
    return mindex


# generate one Series from a dict entry
def gen_Series(dict_x, varline_x):
    m_index = gen_index(varline_x)
    res_Series = pd.Series(dict_x[varline_x], index=m_index)
    return res_Series


# merge a whole dict into one Series
def comb_gen_Series(dict_x):
    x_keys = dict_x.keys()
    res_Series = pd.concat([gen_Series(dict_x, k) for k in x_keys])
    return res_Series


# pack up output of `run_suews`
def pack_df_output(dict_output):
    # TODO: add output levels as in the Fortran version
    df_output = pd.DataFrame(dict_output).T
    # df_output = pd.concat(dict_output).to_frame().unstack()
    # set index level names
    index = df_output.index.set_names(["datetime", "grid"])
    # clean columns
    columns = gen_MultiIndex(df_output.iloc[0])
    values = np.apply_along_axis(np.hstack, 1, df_output.values)
    df_output = pd.DataFrame(values, index=index, columns=columns)
    return df_output


def pack_df_state(dict_state):
    df_state = pd.DataFrame(dict_state).T
    # df_state = pd.concat(dict_state).to_frame().unstack()
    # set index level names
    df_state.index = df_state.index.set_names(["datetime", "grid"])

    return df_state


def pack_df_output_array(dict_output_array, df_forcing):
    grid_list = list(dict_output_array.keys())
    grid_start = grid_list[0]
    col_df = gen_MultiIndex(dict_output_array[grid_start])
    dict_df = {}
    for grid in grid_list:
        array_grid = np.hstack([v[:, 5:] for v in dict_output_array[grid].values()])
        df_grid = pd.DataFrame(array_grid, columns=col_df, index=df_forcing.index)

        dict_df.update({grid: df_grid})

    # join results of all grids
    df_grid_res = pd.concat(dict_df, keys=dict_df.keys())

    # set index level names
    df_grid_res.index.set_names(["grid", "datetime"], inplace=True)

    return df_grid_res


# resample supy output
def resample_output(df_output, freq="60T", dict_aggm=dict_var_aggm):

    # get grid and group names
    list_grid = df_output.index.get_level_values("grid").unique()
    list_group = df_output.columns.get_level_values("group").unique()

    # resampling output according to different rules defined in dict_aggm
    # note the setting in .resample: (closed='right',label='right')
    # which is to conform to SUEWS convention
    # that timestamp refer to the ending of previous period
    df_rsmp = pd.concat(
        {
            grid: pd.concat(
                {
                    group: df_output.loc[grid, group]
                    .resample(freq, closed="right", label="right",)
                    .agg(dict_aggm[group])
                    for group in list_group
                },
                axis=1,
                names=["group", "var"],
            )
            for grid in list_grid
        },
        names=["grid"],
    )

    # clean results
    df_rsmp = df_rsmp.dropna(how="all", axis=0)

    return df_rsmp


def proc_df_rsl(df_output, debug=False):
    try:
        # if we work on the whole output with multi-index columns
        df_rsl_raw = df_output["RSL"].copy()
    except:
        # if we directly work on the RSL output
        df_rsl_raw = df_output.copy()

    try:
        # drop unnecessary columns if existing
        df_rsl_data = df_rsl_raw.drop(["Year", "DOY", "Hour", "Min", "Dectime"], axis=1)
    except:
        df_rsl_data = df_rsl_raw

    # retrieve data for plotting
    df_rsl = df_rsl_data.iloc[:, : 30 * 4]
    df_rsl.columns = (
        df_rsl.columns.str.split("_")
        .map(lambda l: tuple([l[0], int(l[1])]))
        .rename(["var", "level"])
    )
    df_rsl_proc = df_rsl.stack()
    if debug:
        # retrieve debug variables
        df_rsl_debug = df_rsl_data.iloc[:, 120:]
        return df_rsl_proc, df_rsl_debug
    else:
        return df_rsl_proc
