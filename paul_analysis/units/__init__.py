#!/usr/bin/env python
""" init file for GADGET/AREPO reading modules 

Init procedure for simread module for reading GADGET and AREPO simulation output.
Imports all modules, but does no other init functions.

"""


__author__ = "Paul Torrey"
__copyright__ = "Copyright 2014, The Authors"
__credits__ = ["Paul Torrey"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Paul Torrey"
__email__ = "ptorrey@mit.harvard.edu"
__status__ = "Beta -- forever."


#from units import constants as constants
# access these vals with:
# 
# import units
# print units.constants["time_code_to_gyr"]
# 
##

constants =     {
    
                        "time_code_to_gyr":     0.9778,
                        "gamma":                5.0/3.0,
                        "g_minus_1":                2.0/3.0,
                        "PROTONMASS":           1.6726e-24,
                        "ElectronMass":		9.109e-28,		# cgs
                        "BoltzMann_ergs":       1.3806e-16,		# cgs
                        "c":			2.997e10,		# cgs
                        "UnitMass_in_g":        1.989e43,
                        "UnitEnergy_in_cgs":    1.989e53,
                        "UnitDist_in_cm":	3.086e21,
                        "dummy":                1234,
                        "sigma_thomson":        6.6524e-25,		# cm^2	
                }



