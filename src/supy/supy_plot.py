# IQR filling plot:
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


def plot_day_clm(df_var, fig=None, ax=None, **kwargs):
    """Short summary.

    Parameters
    ----------
    df_var : pd.DataFrame
        DataFrame containing variables to plot with datetime as index

    Returns
    -------
    MPL.figure
        figure showing median lines and IQR in shadings

    """
    if fig is None and ax is None:
        fig, ax = plt.subplots()
    elif fig is None:
        fig = ax.get_figure()
    elif ax is None:
        ax = fig.gca()
    # plt.clf()
    # group by hour and minute
    grp_sdf_var = df_var.groupby(
        [df_var.index.hour.rename('hr'),
         df_var.index.minute.rename('min')])
    # get index
    idx = [pd.datetime(2014, 1, 1, h, m)
           for h, m in sorted(grp_sdf_var.groups.keys())]
    idx = pd.date_range(idx[0], idx[-1], periods=len(idx))
    idx = mdates.date2num(idx)

    # calculate quartiles
    quar_sel_pos_clm = grp_sdf_var.quantile(
        [.75, .5, .25]).unstack().set_index(idx)
    # fig, ax = plt.subplots(1)

    for var in quar_sel_pos_clm.columns.levels[0]:
        df_x = quar_sel_pos_clm.loc[:, var]
        y0 = df_x[0.5]
        y1, y2 = df_x[0.75], df_x[0.25]
        y0.plot(ax=ax, label=var).fill_between(
            quar_sel_pos_clm.index, y1, y2, alpha=0.3)
    # add legend
    ax.legend(title='variable')
    # adjust xticks formar
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=np.arange(0, 23, 3)))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    return fig, ax

# comparison plot with 1:1 line added:


def plot_comp(df_var, fig=None, ax=None, **kwargs):
    """Short summary.

    Parameters
    ----------
    df_var : pd.DataFrame
        DataFrame containing variables to plot with datetime as index

    Returns
    -------
    MPL.figure
        figure showing 1:1 line plot

    """
    if fig is None and ax is None:
        fig, ax = plt.subplots()
    elif fig is None:
        fig = ax.get_figure()
    elif ax is None:
        ax = fig.gca()
    # plt.clf()
    # plt.cla()
    # ax = sns.regplot(
    #     x='Obs', y='Sim',
    #     data=df_var,
    #     fit_reg=True)

    # add regression expression
    df_var_fit = df_var.dropna(how='any')
    # regr = linear_model.LinearRegression()
    # val_x = df_var_fit['Obs'].values.reshape(-1, 1)
    # val_y = df_var_fit['Sim'].values.reshape(-1, 1)
    # regr.fit(val_x, val_y)
    val_x = df_var_fit['Obs']
    val_y = df_var_fit['Sim']
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        val_x, val_y)
    mae = (val_y - val_x).abs().mean()

    sns.regplot(
        x='Obs', y='Sim',
        data=df_var,
        ax=ax,
        fit_reg=True,
        line_kws={
            'label': "y={0:.2f}x+{1:.2f}".format(slope, intercept) +
            '\n' + '$R^2$={0:.4f}'.format(r_value) +
            '\n' + 'MAE={0:.2f}'.format(mae) +
            '\n' + 'n={}'.format(df_var.shape[0])
        },
        **kwargs
    )
    # ax.plot(val_x, y_pred, color='red', linewidth=2,
    #         label='r2= ' + str("%.3f" % r2) + '\n' +
    #         'y=' + str("%.3f" % a[0][0]) + 'x+' + str("%.2f" % b[0]))

    # ax.legend(fontsize=15)
    ax.legend()
    # ax.set_title(var + '_' + title)

    # set equal plotting range
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    lim_low, lim_high = np.min([x0, y0]), np.max([x1, y1])
    ax.set_xlim(lim_low, lim_high)
    ax.set_ylim(lim_low, lim_high)

    # set 1:1 aspect ratio
    ax.set_aspect('equal')

    # add 1:1 line
    ax.plot([lim_low, lim_high], [lim_low, lim_high],
            color='red', linewidth=1, zorder=0)

    # fig = ax.figure

    return fig, ax
