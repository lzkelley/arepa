#!/usr/bin/env python
"""Routines to determine and manipulate the file name formats of various
output files."""

__author__ = "Paul Torrey and contributing authors"
__copyright__ = "Copyright 2014, The Authors"
__credits__ = ["Paul Torrey and contributing authors"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Paul Torrey"
__email__ = "ptorrey@mit.harvard.edu"
__status__ = "Beta -- forever."

import os
import util

def python_path():
    p_path = os.path.dirname(util.__file__)
    return p_path[:p_path.index('/util')]

def repository_path():
    p_path = os.path.dirname(util.__file__)
    print p_path
    r_path = p_path[:p_path.index('/Python/util')]
    print r_path
    return r_path

def c_routines_dir():
    r_path = repository_path()
    c_path = r_path+'/C/c_libraries'
    print c_path
    return c_path

def sps_dir():
    r_path = repository_path()
    c_path = r_path+'/Python/physicalmodels/stellarproperties/data/'
    print c_path
    return c_path

def atten_dir():
    r_path = repository_path()
    c_path = r_path+'/C/c_libraries/Attenuation/'
    print c_path
    return c_path