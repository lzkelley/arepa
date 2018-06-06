##
## function to return the AGN spectrum
##   given nu (in Hz), and log(L_bol/L_sun), returns 
##   log(nuL_nu/L_sun)
##
##	keywords to shortcut for various bands:
##	  BB = B-band, IR = mid-IR (15 microns), 
##    SX = soft X-rays (0.5-2 keV), 
##    HX = hard X-rays (2-10 keV)
##
##  keywords for different spectra:
##    MARCONI = Marconi et al. 2004-ish spectrum -- it's not *really* 
##      the same (updated w. more recent estimates of the X-ray 
##      input parameters, i.e. Gamma & reflection angle, from Tozzi et al. 2006 
##		as well as the more recent Steffen et al. 2006 calibration of the 
##      alpha_ox relation (specifically their l_uv-l_2keV bisector measurement, 
##      which is a subtle difference but is very important for the bolometric 
##      corrections of the most luminous X-ray sources
##    HRH = Hopkins, Richards, & Hernquist 2006 template spectrum
##	    compiled from a number of observations therein. It's a little closer to 
##      a typical observed spectrum, with a more detailed modeling of a number of 
##      parts of the spectrum (i.e. more continuum shape). the major difference is 
##      that it includes a hot dust component, which is important to account for 
##      since much of the energy comes out there. however, if you want to model 
##      that generation self-consistently, you should use the "input" marconi et al. 
##      spectrum. still, the galaxy-scale obscuration doesn't necessarily produce this 
##      and the alpha-ox calibrations are for objects with this feature, so it is 
##      generally more representative
##    RICHARDS = Richards et al. 2006 all-quasar mean SED: if you want a 
##      non-luminosity-dependent spectrum for whatever reason, this is an improved 
##      version of Elvis et al. 1994 (which is too X-ray bright)
##
##    SDSS can be added as a keyword with any of the model spectra, and it will 
##      overlay the vanden Berk et al. 2001 median SDSS SED over the relevant 
##      portion of the spectrum :: basically does so with a sliding continuum determination
##      that then overlays this spectrum, such that the integrated luminosity in the 
##      entire bandpass and the continuum spectrum are conserved (i.e. you can put this 
##      over an arbitrary continuum
##
##	unless you're using the Richards et al. 2006 spectra, if you want a reference for this, 
##    it should be Hopkins, Richards, & Hernquist 2006. Even the "MARCONI" key spectrum 
##    is substantially modified as described. If you want the additional relevant 
##    observational compilations on which it's all based, the list is : 
##    Richards et al. 2006, Steffen et al. 2006, Hatziminaoglou et al. 2005, 
##    Tozzi et al. 2006, Strateva et al. 2005, Telfer et al. 2002, Vanden Berk et al. 2001, 
##    George et al. 1998, Elvis et al. 1994, Marconi et al. 2004 (not really an observational 
##    paper but the methodology does follow them), with the X-ray reflection component 
##    following Ueda et al. 2003 in the PEXRAV code, Magdziarz & Zdziarski 1995
##
##
##
def agn_spectrum( nu_in_Hz, log_l_bol, \
	BB=0, IR=0, SX=0, HX=0, \
	HRH=0, MARCONI=0, RICHARDS=0, \
	SDSS=0 ):

    import utilities as util
    import numpy as np
    import math
    import ctypes

    ## location of shared library
    exec_call=util.return_python_routines_homedir()+'/agn_spectrum/agn_spectrum_py.so'
    lib=ctypes.cdll[exec_call];

    if (1 == BB) : nu_in_Hz = -1.0
    if (1 == IR) : nu_in_Hz = -2.0
    if (1 == SX) : nu_in_Hz = -3.0
    if (1 == HX) : nu_in_Hz = -4.0	
    spectrum_key = 0
    if (1 == HRH) 		: spectrum_key = 0
    if (1 == MARCONI) 	: spectrum_key = 1
    if (1 == RICHARDS) 	: spectrum_key = 2	
    sloan_key    = 0
    if (1 == SDSS)		: sloan_key = 1

    nu_in_Hz=np.array(nu_in_Hz,ndmin=1,dtype='d');
    log_l_bol=np.array(log_l_bol,ndmin=1,dtype='d');
    
    N_nu = nu_in_Hz.shape[0]
    N_lum = log_l_bol.shape[0]
    l_band_all = np.zeros((N_nu,N_lum),dtype='d')

    out_cast = ctypes.c_double*N_nu
    for i in range(N_lum):
        #l_band_vec = np.zeros(N_nu,dtype='d')
        l_band_vec = out_cast()
        log_l_bol_pass = log_l_bol[i]

        lib.main( ctypes.c_int(N_nu),\
            nu_in_Hz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
            ctypes.c_double(log_l_bol_pass),\
            ctypes.c_int(spectrum_key),\
            ctypes.c_int(sloan_key),\
            ctypes.byref(l_band_vec) );
            #l_band_vec.ctypes.data_as(ctypes.POINTER(ctypes.c_double)));
        ## now put the output arrays into a useful format 
        l_band_vec = np.copy(np.ctypeslib.as_array(l_band_vec));
        l_band_all[:,i] = l_band_vec

    return l_band_all

