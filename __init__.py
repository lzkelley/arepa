"""
"""

import os
import sys

_CDIR = os.path.split(os.path.abspath(__file__))[0]

_PATH_OUTPUT = "analysis_output"
PATH_OUTPUT = os.path.join(_CDIR, _PATH_OUTPUT)

_check_paths = [PATH_OUTPUT]

ILLUSTRIS_PYTHON_PATH = "/n/home00/lkelley/illustris/"
if ILLUSTRIS_PYTHON_PATH not in sys.path:
    sys.path.append(ILLUSTRIS_PYTHON_PATH)
    print("Added '{}' `to sys.path`".format(ILLUSTRIS_PYTHON_PATH))

# Make sure required paths exist
for path in _check_paths:
    if not os.path.exists(path):
        os.makedirs(path)
        print("Created '{}'".format(path))
