# version info for supy

from supy_driver import __version__ as sd_ver
from ._env import path_supy_module

import pandas as pd

ser_ver = pd.read_json(
    path_supy_module / "supy_version.json", typ="series", convert_dates=False
)
__version__ = f"{ser_ver.ver_milestone}.{ser_ver.ver_major}.{ser_ver.ver_minor}{ser_ver.ver_remark}"
__version_driver__ = sd_ver


def show_version():
    """print `supy` and `supy_driver` version information.
    """
    print("SuPy versions")
    print("-------------")
    print(f"supy: {__version__}")
    print(f"supy_driver: {__version_driver__}")
    print("\n=================")
    print("SYSTEM DEPENDENCY")
    pd.show_versions()