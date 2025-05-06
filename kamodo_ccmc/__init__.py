import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import warnings

import kamodo # get installed kamodo

# insert Kamodo class into kamodo namespace
from kamodo.kamodo import Kamodo

import readers
import tools
import flythrough
import kamodo-filedriver

__version__ = '23.3.2'

