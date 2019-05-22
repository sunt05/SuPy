# command line tools
import click
import sys

from ..util._converter import convert_table, list_ver_from, list_ver_to

# run the whole supy workflow mimicking SUEWS binary
@click.command(short_help='convert SUEWS input tables between different versions.')
@click.option('-f', '--from', 'fromVer',
              help='Version to convert from',
              type=click.Choice(list_ver_from),
              required=True)
@click.option('-t', '--to', 'toVer',
              help='Version to convert to',
              type=click.Choice(list_ver_to),
              required=True)
@click.option('-i', '--input', 'fromDir',
              help='Original directory to convert',
              type=click.Path(),
              required=True)
@click.option('-o', '--output', 'toDir',
              help='New directory to create for converted tables',
              type=click.Path(),
              required=True)
def convert_table_cmd(fromDir, toDir, fromVer, toVer):
    '''convert SUEWS input tables between different versions.'''
    convert_table(fromDir, toDir, fromVer, toVer)
