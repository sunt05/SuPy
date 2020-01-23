# IQR filling plot:


import numpy as np
import pandas as pd




def plot_day_clm(df_var, fig=None, ax=None):
    """Produce a ensemble diurnal climatologies with uncertainties shown in inter-quartile ranges.

    Parameters
    ----------
    df_var : pd.DataFrame
        DataFrame containing variables to plot with datetime as index

    Returns
    -------
    MPL.figure
        figure showing median lines and IQR in shadings

    """
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    if fig is None and ax is None:
        fig, ax = plt.subplots()
    elif fig is None:
        fig = ax.get_figure()
    elif ax is None:
        ax = fig.gca()

    # group by hour and minute
    grp_sdf_var = df_var.groupby(
        [df_var.index.hour.rename("hr"), df_var.index.minute.rename("min")]
    )
    # get index
    year = df_var.index.year.min()
    month = df_var.index.month.min()
    day = df_var.index.day.min()
    idx = [
        pd.datetime(year, month, day, h, m)
        for h, m in sorted(grp_sdf_var.groups.keys())
    ]
    idx = pd.date_range(idx[0], idx[-1], periods=len(idx))
    idx = mdates.date2num(idx)

    # calculate quartiles
    quar_sel_pos_clm = grp_sdf_var.apply(lambda grp: grp.quantile([0.75, 0.5, 0.25]))
    # rearrangement
    quar_sel_pos_clm = quar_sel_pos_clm.unstack()
    # indexing with proper datetime
    quar_sel_pos_clm = quar_sel_pos_clm.set_index(idx)
    # quar_sel_pos_clm = grp_sdf_var.quantile(
    #     [.75, .5, .25]).unstack().set_index(idx)
    # fig, ax = plt.subplots(1)

    for var in quar_sel_pos_clm.columns.levels[0]:
        df_x = quar_sel_pos_clm.loc[:, var]
        y0 = df_x[0.5]
        y1, y2 = df_x[0.75], df_x[0.25]
        y0.plot(ax=ax, label=var).fill_between(
            quar_sel_pos_clm.index, y1, y2, alpha=0.3
        )
    # add legend
    ax.legend(title="variable")
    # adjust xticks formar
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=np.arange(0, 23, 3)))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

    return fig, ax


# comparison plot with 1:1 line added:
def plot_comp(df_var, fig=None, ax=None):
    """Produce a scatter plot with linear regression line to compare simulation results and observations.

    Parameters
    ----------
    df_var : pd.DataFrame
        DataFrame containing variables to plot with datetime as index.
        Two columns, 'Obs' and 'Sim' for observations and simulation results, respectively, must exist.

    Returns
    -------
    MPL.figure
        figure showing 1:1 line plot

    """
    # import when used for better performance in loading supy
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy import stats

    if fig is None and ax is None:
        fig, ax = plt.subplots()
    elif fig is None:
        fig = ax.get_figure()
    elif ax is None:
        ax = fig.gca()

    # add regression expression
    df_var_fit = df_var.dropna(how="any")
    val_x = df_var_fit["Obs"]
    val_y = df_var_fit["Sim"]
    slope, intercept, r_value, p_value, std_err = stats.linregress(val_x, val_y)
    mae = (val_y - val_x).abs().mean()

    sns.regplot(
        x="Obs",
        y="Sim",
        data=df_var,
        ax=ax,
        fit_reg=True,
        line_kws={
            "label": "y={0:.2f}x{1}{2:.2f}".format(
                slope, "+" if intercept > 0 else "", intercept
            )
            + "\n"
            + "$R^2$={0:.4f}".format(r_value)
            + "\n"
            + "MAE={0:.2f}".format(mae)
            + "\n"
            + "n={}".format(df_var.shape[0])
        },
    )

    ax.legend()

    # set equal plotting range
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    lim_low, lim_high = np.min([x0, y0]), np.max([x1, y1])
    ax.set_xlim(lim_low, lim_high)
    ax.set_ylim(lim_low, lim_high)

    # set 1:1 aspect ratio
    ax.set_aspect("equal")

    # add 1:1 line
    ax.plot(
        [lim_low, lim_high], [lim_low, lim_high], color="red", linewidth=1, zorder=0
    )

    return fig, ax
