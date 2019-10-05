from logging.handlers import TimedRotatingFileHandler
import sys
import logging
import inspect
from pathlib import Path

########################################################################
# this file provides variable and functions useful for the whole module.
########################################################################


# define local path for loading resources in this package
path_supy_module = Path(inspect.getsourcefile(lambda: 0)).resolve().parent


# set up logger
FORMATTER = logging.Formatter(
    "%(asctime)s — %(name)s — %(levelname)s — %(message)s")
LOG_FILE = "SuPy.log"


def get_console_handler():
   console_handler = logging.StreamHandler(sys.stdout)
   console_handler.setFormatter(FORMATTER)
   return console_handler


def get_file_handler():
   file_handler = TimedRotatingFileHandler(LOG_FILE, when='midnight')
   file_handler.setFormatter(FORMATTER)
   return file_handler


def get_logger(logger_name, level=logging.DEBUG):
   logger = logging.getLogger(logger_name)
   # better to have too much log than not enough
   logger.setLevel(level)
   logger.addHandler(get_console_handler())
   logger.addHandler(get_file_handler())
   # with this pattern, it's rarely necessary to propagate the error up to parent
   logger.propagate = False
   return logger


logger_supy = get_logger("SuPy",logging.INFO)
logger_supy.debug("a debug message from SuPy")
