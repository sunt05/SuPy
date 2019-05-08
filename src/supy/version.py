# version info for supy

from supy_driver import __version__ as sd_ver
ver_milestone = 2019
ver_major = 5
ver_minor = 8
ver_remark = 'dev'
__version__ = '{ver_milestone}.{ver_major}.{ver_minor}{ver_remark}'.format(
    ver_milestone=ver_milestone,
    ver_major=ver_major,
    ver_minor=ver_minor,
    ver_remark=ver_remark,
)

__version_driver__ = sd_ver

def show_version():
    print(f'supy: {__version__}')
    print(f'supy_driver: {__version_driver__}')
