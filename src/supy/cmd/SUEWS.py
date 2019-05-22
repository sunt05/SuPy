# command line tools
import click
import sys

from .._supy_module import (
    init_supy, run_supy, save_supy, load_forcing_grid,
    pd, Path)

# run the whole supy workflow mimicking SUEWS binary
@click.command()
@click.argument('path_runcontrol', default='./RunControl.nml')
def SUEWS(path_runcontrol):
    try:
        path_runcontrol=Path(path_runcontrol).resolve()
    except:
        click.echo(f'{str(path_runcontrol)} not existing!')
        sys.exit()

    # init supy
    df_state_init = init_supy(path_runcontrol)

    # load forcing
    list_grid = df_state_init.index
    df_forcing = pd.concat(
        {
            grid: load_forcing_grid(path_runcontrol, grid)
            for grid in list_grid
        },
        names=['grid', 'datetime']
    )

    # run supy
    df_output, df_state_final = run_supy(df_forcing, df_state_init)

    # save result
    list_out_files = save_supy(df_output, df_state_final)
    # show output files
    for file in list_out_files:
        click.echo(file)

    # return
    click.echo('SUEWS run successfully done!')
