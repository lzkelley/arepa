"""
"""

import os
import sys
import socket 

_CDIR = os.path.split(os.path.abspath(__file__))[0]

_PATH_OUTPUT = "analysis_output"
PATH_OUTPUT = os.path.join(_CDIR, _PATH_OUTPUT)

_check_paths = [PATH_OUTPUT]

# Determine Hostname (i.e. host system) and set host-dependent parameters
# -------------------------------------------------------------------------------

_hostname = socket.gethostname().lower()
ILLUSTRIS_PYTHON_PATH = None
if 'daedalus' in _hostname:
    HOST = 'daedalus'
elif 'harvard' in _hostname:
    HOST = 'odyssey'
    ILLUSTRIS_PYTHON_PATHS = "/n/home00/lkelley/illustris/"
elif 'ufhpc' in _hostname:
    HOST = 'hipergator'
    ILLUSTRIS_PYTHON_PATH = "/home/lkelley/arepo_sims"

    import matplotlib
    matplotlib.use('Agg')

else:
    print("WARNING: unrecognized host!")

# Add Path to `illustris_python` module
'''
ILLUSTRIS_PYTHON_PATH = "/n/home00/lkelley/illustris/"
if ILLUSTRIS_PYTHON_PATH not in sys.path:
    sys.path.append(ILLUSTRIS_PYTHON_PATH)
    print("Added '{}' `to sys.path`".format(ILLUSTRIS_PYTHON_PATH))
'''
if (ILLUSTRIS_PYTHON_PATH is not None) and (ILLUSTRIS_PYTHON_PATH not in sys.path):
    if os.path.exists(ILLUSTRIS_PYTHON_PATH):
        sys.path.append(ILLUSTRIS_PYTHON_PATH)
        print("Added '{}' `to sys.path`".format(ILLUSTRIS_PYTHON_PATH))
    else:
        print("WARNING: illustris path '{}' does not exist!".format(ILLUSTRIS_PYTHON_PATH))

# Make sure required paths exist
for path in _check_paths:
    if not os.path.exists(path):
        os.makedirs(path)
        print("Created '{}'".format(path))
