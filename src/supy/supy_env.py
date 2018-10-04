import inspect
import os

# define local path for loading resources in this package
path_supy_module = os.path.dirname(
    os.path.abspath(inspect.getsourcefile(lambda: 0)))
