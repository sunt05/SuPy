import inspect
from pathlib import Path

# define local path for loading resources in this package
path_supy_module = Path(inspect.getsourcefile(lambda: 0)).resolve().parent
